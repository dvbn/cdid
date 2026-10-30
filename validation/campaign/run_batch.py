#!/usr/bin/env python3
"""Portable reliability campaign: single-batch supervisor and fake-worker self-test suite.
Uses only the Python 3 standard library. Spawns one fresh Rscript process per; operation via subprocess without shell=True, enforces per-operation timeouts and
batch wall budgets with a 30-second cleanup reserve, manages host/PID-verified; locks without age-based stealing, reconciles proven-dead interrupted attempts
with linked replacement attempts, and rebuilds the attempt index from disk records.
"""
import argparse, csv, datetime, errno, hashlib, json, os, signal, socket, subprocess, sys, tempfile, time, uuid, platform, fcntl; from pathlib import Path
def now_utc() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()
def write_json_atomic(data: dict, path: Path, overwrite: bool = True):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite: raise FileExistsError(f"Refusing to overwrite immutable record: {path}")
    tmp_path = path.parent / f".{path.name}.tmp.{os.getpid()}.{uuid.uuid4().hex[:8]}"; tmp_path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8"); 
    if overwrite: os.replace(tmp_path, path)
    else:
        try: os.link(tmp_path, path)
        finally: tmp_path.unlink(missing_ok=True)
def check_pid_liveness(pid) -> str:
    """Check process liveness on the current host using os.kill(pid, 0) without /proc.
    Returns:
      'dead'    if ProcessLookupError (ESRCH) proves the process does not exist; 'alive'   if os.kill(pid, 0) succeeds; 'unknown' if PermissionError (EPERM) or invalid PID prevents proving death
    """
    if pid is None: return "dead"
    try:
        pid_int = int(pid)
    except (TypeError, ValueError):
        return "unknown"
    if pid_int <= 0: return "unknown"
    try:
        os.kill(pid_int, 0); return "alive"
    except ProcessLookupError:
        return "dead"
    except PermissionError:
        return "unknown"
    except OSError as e:
        if e.errno == errno.ESRCH: return "dead"
        return "unknown"
def check_pgid_liveness(pgid) -> str:
    if pgid is None: return "dead"
    try:
        pgid_int = int(pgid)
    except (TypeError, ValueError):
        return "unknown"
    if pgid_int <= 0: return "unknown"
    try:
        os.killpg(pgid_int, 0); return "alive"
    except ProcessLookupError:
        return "dead"
    except PermissionError:
        return "unknown"
    except OSError as e:
        if e.errno == errno.ESRCH: return "dead"
        return "unknown"
def kill_and_reap_worker(proc: subprocess.Popen, grace_seconds: float = 2.0):
    # start_new_session fixes pgid=pid even if the leader exits before its descendants.
    pgid = proc.pid
    try: os.killpg(pgid, signal.SIGTERM)
    except ProcessLookupError: pass
    deadline = time.monotonic() + grace_seconds
    while time.monotonic() < deadline:
        if check_pgid_liveness(pgid) == "dead": break
        proc.poll()
        time.sleep(0.02)
    try: os.killpg(pgid, signal.SIGKILL)
    except ProcessLookupError: pass
    proc.wait(timeout=2.0)
def rng_stream_index(stream_index: int, rep: int, purpose: int) -> int:
    return ((int(stream_index) - 1) * 100 + (int(rep) - 1)) * 5 + int(purpose)
def expand_manifest_tasks(manifest_path: Path, max_reps_override: int = 0):
    rows = list(csv.DictReader(manifest_path.open(encoding="utf-8")))
    if not rows: return []
    # Check whether rows use replicate_scope="1..N" or explicit replicate numbers
    scope0 = rows[0].get("replicate_scope", "1..100").strip()
    if ".." in scope0:
        lo_str, hi_str = scope0.split("..", 1); lo_rep, hi_rep = int(lo_str), int(hi_str)
        if max_reps_override > 0:
            hi_rep = min(hi_rep, max_reps_override)
        expanded = []
        for rep in range(lo_rep, hi_rep + 1):
            for r in rows:
                t_id = f"r{rep:03d}__{r['task_id']}"; p_id = f"r{rep:03d}__{r['parent_task_id']}" if r.get("parent_task_id") else ""; item = dict(r)
                item["template_task_id"] = r["task_id"]; item["task_id"] = t_id; item["parent_task_id"] = p_id; item["replicate"] = rep; expanded.append(item)
        return expanded
    else:
        expanded = []
        for r in rows:
            rep = int(r.get("replicate", scope0))
            if max_reps_override > 0 and rep > max_reps_override: continue
            item = dict(r); item["template_task_id"] = r["task_id"]; item["replicate"] = rep; expanded.append(item)
        return expanded
def scan_existing_attempts(records_dir: Path):
    """Scan all immutable .result.json and .open.json files on disk."""
    by_task = {}; open_records = []
    if not records_dir.exists(): return by_task, open_records
    for p in sorted(records_dir.glob("*.result.json")):
        rec = json.loads(p.read_text(encoding="utf-8")); tid = rec["task_id"]; by_task.setdefault(tid, []).append(rec)
    for tid, recs in by_task.items():
        recs.sort(key=lambda x: int(x.get("attempt_number", 1)))
    for p in sorted(records_dir.glob("*.open.json")):
        open_records.append((p, json.loads(p.read_text(encoding="utf-8"))))
    return by_task, open_records
TERMINAL_DISPOSITIONS = frozenset({"success", "partial-support", "expected-refusal", "unexpected-success", "unexpected-error", "timeout", "parent-unavailable", "unavailable-comparator", "not-applicable" })
def task_terminal_record(attempts_for_task: list):
    for rec in reversed(attempts_for_task):
        if rec.get("disposition") in TERMINAL_DISPOSITIONS: return rec
    return None
def reconcile_interrupted_attempts(records_dir: Path, owner_id: str, reason: str):
    """Mark any open attempt files as interrupted atomically and remove the .open.json marker."""
    reconciled = []
    if not records_dir.exists(): return reconciled
    for open_path in sorted(records_dir.glob("*.open.json")):
        open_rec = json.loads(open_path.read_text(encoding="utf-8")); att_id = open_rec["attempt_id"]; res_path = records_dir / f"{att_id}.result.json"
        if not res_path.exists():
            interrupted_rec = {**open_rec, "disposition": "interrupted", "settled_utc": now_utc(), "reconciled_by_owner_id": owner_id, "reconciliation_reason": reason, "censored_timeout": False}
            write_json_atomic(interrupted_rec, res_path, overwrite=False); reconciled.append(att_id)
        open_path.unlink(missing_ok=True)
    return reconciled
def _acquire_or_reconcile_lock(out_dir: Path, owner_id: str, reconcile_cross_host: str = ""):
    lock_path = out_dir / "lock.json"; records_dir = out_dir / "records"; cur_host = socket.gethostname(); cur_pid = os.getpid(); cur_pgid = os.getpgrp(); reconciled_attempts = []
    if lock_path.exists():
        existing = json.loads(lock_path.read_text(encoding="utf-8")); lock_host = existing.get("hostname", ""); lock_owner = existing.get("owner_id", "")
        if lock_host != cur_host:
            if reconcile_cross_host and reconcile_cross_host == lock_owner:
                reconciled_attempts = reconcile_interrupted_attempts(records_dir, owner_id, f"explicit_cross_host_reconciliation_from_{lock_host}")
            else:
                raise RuntimeError(f"Refusing lock acquisition: lock at {lock_path} is held by another host ({lock_host}, owner_id={lock_owner}). Cross-host process liveness is unknown; verify that host is stopped and pass --reconcile-cross-host {lock_owner}.")
        else:
            sup_state = check_pid_liveness(existing.get("supervisor_pid")) if existing.get("supervisor_pid") else "unknown"; wrk_state = check_pid_liveness(existing.get("worker_pid")); pg_state = check_pgid_liveness(existing.get("worker_pgid"))
            if sup_state != "dead" or wrk_state != "dead" or pg_state != "dead": raise RuntimeError(f"Refusing lock acquisition: existing same-host lock (owner_id={lock_owner}) is not proven dead "
                f"(supervisor={sup_state}, worker_pid={wrk_state}, worker_pgid={pg_state}). " f"Age is diagnostic only; never stealing an active or unknown lock.")
            if existing.get("launch_pending") and reconcile_cross_host != lock_owner:
                raise RuntimeError("Interrupted launch has unknown worker identity; reconcile the recorded owner explicitly")
            reconciled_attempts = reconcile_interrupted_attempts(records_dir, owner_id, "same_host_proven_dead_reclamation")
    else:
        # A missing lock never proves that an orphan worker is dead.
        for op in records_dir.glob("*.open.json"):
            rec = json.loads(op.read_text())
            if reconcile_cross_host and reconcile_cross_host == rec.get("owner_id"): continue
            if not rec.get("supervisor_pid") or rec.get("hostname") != cur_host or rec.get("launch_pending") or any(
                check_pid_liveness(rec.get(k)) != "dead" for k in ("supervisor_pid", "worker_pid")) or check_pgid_liveness(rec.get("worker_pgid")) != "dead":
                raise RuntimeError("Unowned open attempt has alive/unknown liveness; explicit owner reconciliation required")
        reconciled_attempts = reconcile_interrupted_attempts(records_dir, owner_id, "verified_unlocked_open_attempt_reconciliation")
    lock_payload = {"owner_id": owner_id, "hostname": cur_host, "supervisor_pid": cur_pid, "supervisor_pgid": cur_pgid, "worker_pid": None, "worker_pgid": None, "active_task_id": None, "active_attempt_id": None,
        "created_utc": now_utc(), "updated_utc": now_utc(), "reconciled_interrupted_attempts": reconciled_attempts}
    write_json_atomic(lock_payload, lock_path, overwrite=True); return lock_path, reconciled_attempts
def acquire_or_reconcile_lock(out_dir: Path, owner_id: str, reconcile_cross_host: str = ""):
    # Stable-inode advisory lock serializes check/reconcile/create, including two simultaneous starters.
    with (out_dir / ".lock-transition").open("a+") as guard:
        try: fcntl.flock(guard, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError("Another lock transition is active")
        return _acquire_or_reconcile_lock(out_dir, owner_id, reconcile_cross_host)

def update_lock_worker(lock_path: Path, owner_id: str, worker_pid, worker_pgid, task_id, attempt_id):
    if not lock_path.exists():
        raise RuntimeError("Lock disappeared while worker was active")
    cur = json.loads(lock_path.read_text(encoding="utf-8"))
    if cur.get("owner_id") != owner_id:
        raise RuntimeError(f"Lock ownership lost while updating worker state (expected {owner_id})")
    cur["launch_pending"] = worker_pid == "launching"; cur["worker_pid"] = None if worker_pid == "launching" else worker_pid; cur["worker_pgid"] = worker_pgid; cur["active_task_id"] = task_id; cur["active_attempt_id"] = attempt_id; cur["updated_utc"] = now_utc(); write_json_atomic(cur, lock_path, overwrite=True)
def release_lock(lock_path: Path, owner_id: str):
    if not lock_path.exists():
        return
    try:
        cur = json.loads(lock_path.read_text(encoding="utf-8"))
        if cur.get("owner_id") == owner_id:
            lock_path.unlink(missing_ok=True)
    except Exception:
        pass
def verify_provenance_on_resume(out_dir: Path, campaign_dir: Path, manifest_path: Path, allow_unfrozen: bool = False, repo_dir: str = ""):
    from freeze import checkout_commit, source_snapshot, compute_harness_hashes
    prov_path = out_dir / "provenance.json"
    if not prov_path.exists():
        if allow_unfrozen: return None
        raise RuntimeError(f"Missing {prov_path}. Run freeze.py --freeze with the supplied source commit SHA before running batches.")
    prov = json.loads(prov_path.read_text(encoding="utf-8"))
    if not allow_unfrozen and (not prov.get("source_sha") or prov.get("source_sha") == "FREEZE_SOURCE_SHA"):
        raise RuntimeError("Refusing to run: provenance.json has placeholder FREEZE_SOURCE_SHA.")
    recorded_hashes = prov.get("harness_hashes", {})
    if not recorded_hashes: raise RuntimeError("Missing frozen harness hashes")
    for rel, exp_hash in recorded_hashes.items():
        cand = campaign_dir / rel
        if not cand.is_file() or sha256_file(cand) != exp_hash:
            raise RuntimeError(f"Refusing resume: frozen file hash mismatch or missing file: {rel}")
    if not allow_unfrozen:
        if set(recorded_hashes) != set(compute_harness_hashes(campaign_dir)): raise RuntimeError("Incomplete frozen harness inventory")
        for name in ("manifest.csv", "targets.json"):
            if not (out_dir / name).is_file() or sha256_file(out_dir / name) != recorded_hashes[name]:
                raise RuntimeError(f"Frozen output hash mismatch: {name}")
        if sha256_file(manifest_path) != recorded_hashes["manifest.csv"]: raise RuntimeError("Manifest hash mismatch")
        if sha256_file(out_dir / "seeds_29500.rds") != prov["seeds_rds_sha256"]: raise RuntimeError("Seed hash mismatch")
        if sha256_file(out_dir / "r_runtime.json") != prov["runtime_json_sha256"]: raise RuntimeError("Runtime receipt hash mismatch")
        if sys.version != prov["python_version"] or platform.platform() != prov["python_platform"]: raise RuntimeError("Frozen Python/runtime profile mismatch")
        if not repo_dir or checkout_commit(repo_dir) != prov["source_sha"] or source_snapshot(repo_dir) != prov["source_hashes"]:
            raise RuntimeError("Frozen checkout source identity mismatch")
    return prov
def rebuild_attempts_index(out_dir: Path):
    records_dir = out_dir / "records"; by_task, _ = scan_existing_attempts(records_dir); all_attempts = []; counts_by_disposition = {}
    for tid, recs in sorted(by_task.items()):
        for r in recs:
            disp = r.get("disposition", "unknown"); counts_by_disposition[disp] = counts_by_disposition.get(disp, 0) + 1
            all_attempts.append({"attempt_id": r["attempt_id"], "task_id": r["task_id"], "attempt_number": r["attempt_number"], "linked_previous_attempt_id": r.get("linked_previous_attempt_id", ""), "disposition": disp,
                "call_type": r.get("call_type", ""), "recipe_id": r.get("recipe_id", ""), "cell_id": r.get("cell_id", ""), "replicate": r.get("replicate", 0),
                "operation_seconds": r.get("operation_seconds"), "total_attempt_seconds": r.get("total_attempt_seconds") })
    index_payload = {"updated_utc": now_utc(), "n_distinct_tasks_touched": len(by_task), "n_total_attempts": len(all_attempts), "counts_by_disposition": counts_by_disposition, "attempts": all_attempts}
    write_json_atomic(index_payload, out_dir / "attempts_index.json", overwrite=True); return index_payload
def run_batch(repo_dir: str, manifest_path: Path, out_dir: Path, max_calls: int = 25, wall_seconds: float = 1200.0, op_timeout: float = 300.0, cleanup_reserve: float = 30.0, rscript: str = "Rscript", lib_paths: list = None,
    reconcile_cross_host: str = "", allow_unfrozen: bool = False, max_reps_override: int = 0, fake_worker_cmd: list = None, ):
    if not fake_worker_cmd:
        if allow_unfrozen: raise RuntimeError("Unfrozen execution is restricted to fake workers")
        if op_timeout != 300 or cleanup_reserve < 30: raise ValueError("Campaign requires a 300-second operation cap and at least 30 seconds cleanup")
    if max_calls < 0 or wall_seconds <= 0 or not 0 <= max_reps_override <= 100: raise ValueError("Invalid batch limits")
    t_batch0 = time.monotonic(); deadline = t_batch0 + float(wall_seconds); owner_id = f"batch_{datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')}_{os.getpid()}_{uuid.uuid4().hex[:6]}"
    campaign_dir = Path(__file__).resolve().parent; out_dir = out_dir.resolve(); out_dir.mkdir(parents=True, exist_ok=True); records_dir = out_dir / "records"
    records_dir.mkdir(parents=True, exist_ok=True); data_cache_dir = out_dir / "data_cache"
    data_cache_dir.mkdir(parents=True, exist_ok=True); parents_dir = out_dir / "parents"; parents_dir.mkdir(parents=True, exist_ok=True); work_tmp_dir = out_dir / "tmp"
    work_tmp_dir.mkdir(parents=True, exist_ok=True); targets_path = out_dir / "targets.json"
    if not targets_path.exists():
        targets_path = campaign_dir / "targets.json"
    targets_payload = json.loads(targets_path.read_text(encoding="utf-8")) if targets_path.exists() else {"tasks": {}}; tasks_meta = targets_payload.get("tasks", {})
    prov = verify_provenance_on_resume(out_dir, campaign_dir, manifest_path, allow_unfrozen=allow_unfrozen, repo_dir=repo_dir); frozen_runtime = prov.get("r_runtime") if prov else None
    eff_libs = list(lib_paths) if lib_paths else (prov.get("lib_paths", []) if prov else []); seeds_rds_path = str(out_dir / "seeds_29500.rds")
    lock_path, reconciled_interrupted = acquire_or_reconcile_lock(out_dir, owner_id, reconcile_cross_host=reconcile_cross_host); launched_calls = 0; resolved_no_worker = 0; stopped_reason = "all_pending_completed"
    proc = None
    old_term = signal.signal(signal.SIGTERM, lambda *_: (_ for _ in ()).throw(KeyboardInterrupt()))
    try:
        expanded_tasks = expand_manifest_tasks(manifest_path, max_reps_override=max_reps_override); by_task, _ = scan_existing_attempts(records_dir)
        for item in expanded_tasks:
            if deadline-time.monotonic() < cleanup_reserve:
                stopped_reason = "cleanup_reserve_reached"; break
            tid = item["task_id"]; tpl_id = item["template_task_id"]; rep = int(item["replicate"]); s_idx = int(item["stream_index"]); rec_id = item["recipe_id"]; cid = item["cell_id"]
            call_type = item["call_type"]; parent_tid = item.get("parent_task_id", ""); existing_attempts = by_task.get(tid, [])
            if task_terminal_record(existing_attempts) is not None: continue
            att_num = len(existing_attempts) + 1; att_id = f"{tid}__att{att_num:02d}"; prev_att_id = existing_attempts[-1]["attempt_id"] if existing_attempts else ""
            prev_disp = existing_attempts[-1].get("disposition", "") if existing_attempts else ""
            if existing_attempts and prev_disp != "interrupted":
                # Immutable non-interrupted attempts are never retried
                continue
            tmeta = tasks_meta.get(tpl_id, {}); p_coord_keys = [k for k in item.get("primary_coordinate_set", "").split(";") if k]; parent_rds_path = parents_dir / f"{parent_tid}.rds" if parent_tid else None
            # Check if this is an AGG child whose parent FIT failed or refused
            if call_type == "AGG" and parent_tid:
                parent_attempts = by_task.get(parent_tid, []); parent_term = task_terminal_record(parent_attempts)
                if parent_term is None:
                    # Parent has not settled yet
                    continue
                p_disp = parent_term.get("disposition", "")
                if p_disp not in ("success", "partial-support", "unexpected-success") or not (parent_rds_path and parent_rds_path.exists()):
                    no_worker_rec = {"attempt_id": att_id, "task_id": tid, "template_task_id": tpl_id, "attempt_number": att_num, "linked_previous_attempt_id": prev_att_id, "owner_id": owner_id,
                        "replicate": rep, "stream_index": s_idx, "stream_key": item["stream_key"], "recipe_id": rec_id, "family": item["family"], "cell_id": cid, "call_type": call_type, "package": item["package"],
                        "parent_task_id": parent_tid, "parent_disposition": p_disp, "agg_kind": item.get("agg_kind", ""), "literal_call": item["literal_call"],
                        "primary_target_id": item["primary_target_id"], "primary_truth": tmeta.get("primary_truth_map", {}),
                        "secondary_reference_id": item.get("secondary_reference_id", ""), "secondary_reference": tmeta.get("secondary_reference_map", {}),
                        "secondary_status": item.get("secondary_status", ""), "target_match": item.get("target_match", ""), "causal_validity": item.get("causal_validity", ""),
                        "inference_eligible": item.get("inference_eligible", "") == "TRUE", "moment_regime": item.get("moment_regime", ""), "linked_cell_ids": item.get("linked_cell_ids", cid).split(";"),
                        "disposition": "parent-unavailable", "worker_launched": False, "censored_timeout": False, "supervisor_overhead_seconds": 0.0, "startup_seconds": 0.0, "data_seconds": 0.0,
                        "operation_seconds": 0.0, "total_attempt_seconds": 0.0, "settled_utc": now_utc(),}
                    write_json_atomic(no_worker_rec, records_dir / f"{att_id}.result.json", overwrite=False); by_task.setdefault(tid, []).append(no_worker_rec); resolved_no_worker += 1; continue
            # Check batch call cap and wall-time admission budget before launching a worker
            if launched_calls >= max_calls:
                stopped_reason = "max_calls_reached"; break
            remaining = deadline - time.monotonic()
            if remaining < (op_timeout + cleanup_reserve):
                stopped_reason = "admission_reserve_reached"; break
            t_attempt0 = time.monotonic(); pop_rng_idx = rng_stream_index(s_idx, rep, 1); obs_rng_idx = rng_stream_index(s_idx, rep, 2)
            op_rng_idx = rng_stream_index(s_idx, rep, int(item["rng_purpose"])); data_rds_path = data_cache_dir / f"r{rep:03d}__{rec_id}.rds"
            data_json_path = data_cache_dir / f"r{rep:03d}__{rec_id}.json"; save_parent = bool(tmeta.get("save_parent_rds", False))
            own_parent_rds = str(parents_dir / f"{tid}.rds") if save_parent else ""; task_req_path = work_tmp_dir / f"{att_id}.task.json"
            worker_out_path = work_tmp_dir / f"{att_id}.worker_out.json"; open_marker_path = records_dir / f"{att_id}.open.json"
            task_req = {"attempt_id": att_id, "task_id": tid, "template_task_id": tpl_id, "attempt_number": att_num, "linked_previous_attempt_id": prev_att_id, "replicate": rep, "stream_index": s_idx, "recipe_id": rec_id,
                "cell_id": cid, "call_type": call_type, "package": item["package"], "literal_call": item["literal_call"], "alp": 0.05, "compute_inffunc": item.get("inference_eligible", "TRUE") == "TRUE",
                "point_only_rcs": bool(tmeta.get("point_only_rcs", False)) or (call_type == "DID_RCS"), "expected_disposition": item.get("expected_disposition", "success"),
                "expected_condition": item.get("expected_condition", ""), "primary_coordinates": p_coord_keys, "pop_rng_index": pop_rng_idx, "obs_rng_index": obs_rng_idx, "op_rng_index": op_rng_idx,
                "seeds_rds_path": seeds_rds_path, "data_rds_path": str(data_rds_path), "data_json_path": str(data_json_path), "save_parent_rds": save_parent,
                "parent_rds_path": str(parent_rds_path) if (call_type == "AGG" and parent_rds_path) else own_parent_rds, "worker_out_json": str(worker_out_path), "lib_paths": eff_libs, "frozen_runtime": frozen_runtime, "source_sha": prov["source_sha"] if prov else "FAKE", "repo_path": str(Path(repo_dir).resolve()), "formals_path": str(campaign_dir / "reference/formals.json"),}
            write_json_atomic(task_req, task_req_path, overwrite=True)
            open_record = {"attempt_id": att_id, "task_id": tid, "template_task_id": tpl_id, "attempt_number": att_num, "linked_previous_attempt_id": prev_att_id, "owner_id": owner_id, "hostname": socket.gethostname(), "supervisor_pid": os.getpid(), "worker_pid": None, "worker_pgid": None, "launch_pending": True, "replicate": rep,
                "stream_index": s_idx, "stream_key": item["stream_key"], "recipe_id": rec_id, "family": item["family"], "cell_id": cid, "call_type": call_type, "package": item["package"], "parent_task_id": parent_tid,
                "agg_kind": item.get("agg_kind", ""), "literal_call": item["literal_call"], "pop_rng_index": pop_rng_idx, "obs_rng_index": obs_rng_idx, "op_rng_index": op_rng_idx,
                "data_rds_path": str(data_rds_path), "parent_rds_path": task_req["parent_rds_path"], "primary_target_id": item["primary_target_id"], "primary_truth": tmeta.get("primary_truth_map", {}),
                "auxiliary_truth": tmeta.get("auxiliary_truth_map", {}), "secondary_reference_id": item.get("secondary_reference_id", ""), "secondary_reference": tmeta.get("secondary_reference_map", {}),
                "secondary_status": item.get("secondary_status", ""), "target_match": item.get("target_match", ""), "causal_validity": item.get("causal_validity", ""), "inference_eligible": item.get("inference_eligible", "") == "TRUE",
                "moment_regime": item.get("moment_regime", ""), "equality_ref": item.get("equality_ref", ""), "linked_cell_ids": item.get("linked_cell_ids", cid).split(";"), "started_utc": now_utc(),}
            write_json_atomic(open_record, open_marker_path, overwrite=True)
            cmd = (list(fake_worker_cmd) + ["--task", str(task_req_path)] if fake_worker_cmd else [rscript, "--vanilla", str(campaign_dir / "worker.R"), "--task", str(task_req_path)])
            update_lock_worker(lock_path, owner_id, "launching", None, tid, att_id)
            t_spawn0 = time.monotonic()
            try: proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, start_new_session=True)
            except OSError as exc:
                failed = {**open_record, "launch_pending": False, "worker_launched": False, "disposition": "unexpected-error",
                          "total_attempt_seconds": time.monotonic()-t_attempt0, "error": {"class": [type(exc).__name__], "message": str(exc)}, "settled_utc": now_utc()}
                write_json_atomic(failed, records_dir / f"{att_id}.result.json", overwrite=False)
                open_marker_path.unlink(missing_ok=True)
                update_lock_worker(lock_path, owner_id, None, None, None, None)
                raise RuntimeError("Worker launch failed; immutable failure retained") from exc
            launched_calls += 1
            try:
                w_pgid = os.getpgid(proc.pid)
            except Exception:
                w_pgid = proc.pid
            update_lock_worker(lock_path, owner_id, proc.pid, w_pgid, tid, att_id)
            open_record.update(worker_pid=proc.pid, worker_pgid=w_pgid, launch_pending=False, worker_launched=True)
            write_json_atomic(open_record, open_marker_path)
            timed_out = False
            try:
                stdout_str, stderr_str = proc.communicate(timeout=op_timeout)
            except subprocess.TimeoutExpired:
                timed_out = True; kill_and_reap_worker(proc)
                try:
                    stdout_str, stderr_str = proc.communicate(timeout=1.0)
                except Exception:
                    stdout_str, stderr_str = "", ""
            t_done = time.monotonic(); total_attempt_secs = t_done - t_attempt0; update_lock_worker(lock_path, owner_id, None, None, None, None)
            (work_tmp_dir / f"{att_id}.stdout.log").write_text(stdout_str)
            (work_tmp_dir / f"{att_id}.stderr.log").write_text(stderr_str)
            identity_failure = False
            if timed_out:
                final_rec = {**open_record, "disposition": "timeout", "worker_launched": True, "censored_timeout": True, "timeout_limit_seconds": float(op_timeout), "supervisor_overhead_seconds": round(max(0.0,
                    t_spawn0 - t_attempt0), 6), "startup_seconds": None, "data_seconds": None, "operation_seconds": float(op_timeout), "total_attempt_seconds": round(total_attempt_secs, 6),
                    "exit_code": proc.returncode, "error": {"class": ["timeout_error"], "message": f"Operation exceeded {op_timeout}s wall timeout"}, "warnings": [], "n_warnings": 0, "settled_utc": now_utc()}
            elif worker_out_path.exists():
                try:
                    w_out = json.loads(worker_out_path.read_text(encoding="utf-8"))
                    if not isinstance(w_out, dict): raise ValueError("Worker receipt is not an object")
                except (ValueError, TypeError) as exc:
                    w_out = {"identity_ok": False, "error": {"class": ["invalid_worker_receipt"], "message": str(exc)}}
                identity_failure = (not w_out.get("identity_ok", False) or w_out.get("attempt_id") != att_id or w_out.get("task_id") != tid or proc.returncode != 0)
                if identity_failure:
                    w_out["disposition"] = "unexpected-error"
                    w_out["error"] = {"class": ["worker_identity_or_exit_error"], "message": str(w_out.get("error") or w_out.get("checks"))}
                elif w_out.get("disposition") not in TERMINAL_DISPOSITIONS:
                    w_out["disposition"] = "unexpected-error"
                    w_out["error"] = {"class": ["invalid_worker_disposition"], "message": "Worker disposition is missing or invalid"}
                s_secs = float(w_out.get("startup_seconds", 0.0) or 0.0); d_secs = float(w_out.get("data_seconds", 0.0) or 0.0)
                o_secs = float(w_out.get("operation_seconds", 0.0) or 0.0); sup_secs = max(0.0, total_attempt_secs - (s_secs + d_secs + o_secs))
                final_rec = {**open_record, "disposition": w_out["disposition"], "worker_launched": True, "censored_timeout": False, "supervisor_overhead_seconds": round(sup_secs, 6),
                    "startup_seconds": round(s_secs, 6), "data_seconds": round(d_secs, 6), "data_generated_now": bool(w_out.get("data_generated_now", False)),
                    "operation_seconds": round(o_secs, 6), "operation_user_seconds": w_out.get("operation_user_seconds"),
                    "operation_sys_seconds": w_out.get("operation_sys_seconds"), "total_attempt_seconds": round(total_attempt_secs, 6), "rng_consumed": w_out.get("rng_consumed"), "input_unchanged": w_out.get("input_unchanged"),
                    "data_panel_md5": w_out.get("data_panel_md5"), "data_frame_md5": w_out.get("data_frame_md5"), "parent_saved": bool(w_out.get("parent_saved", False)), "exit_code": proc.returncode,
                    "error": w_out.get("error"), "warnings": w_out.get("warnings", []), "n_warnings": int(w_out.get("n_warnings", 0) or 0), "fit_result": w_out.get("fit_result"),
                    "agg_result": w_out.get("agg_result"), "pretest_result": w_out.get("pretest_result"), "sample": w_out.get("sample"), "if_audit": w_out.get("if_audit"), "settled_utc": now_utc()}
            else:
                final_rec = {**open_record, "disposition": "unexpected-error", "worker_launched": True, "censored_timeout": False, "supervisor_overhead_seconds": round(max(0.0,
                    t_spawn0 - t_attempt0), 6), "startup_seconds": 0.0, "data_seconds": 0.0, "operation_seconds": round(max(0.0, t_done - t_spawn0), 6),
                    "total_attempt_seconds": round(total_attempt_secs, 6), "exit_code": proc.returncode, "error": { "class": ["worker_crash_error"],
                    "message": f"Worker exited with code {proc.returncode} without writing worker_out_json: {stderr_str.strip()}", }, "warnings": [], "n_warnings": 0, "settled_utc": now_utc()}
            write_json_atomic(final_rec, records_dir / f"{att_id}.result.json", overwrite=False); open_marker_path.unlink(missing_ok=True)
            by_task.setdefault(tid, []).append(final_rec)
            if identity_failure: raise RuntimeError(f"Worker runtime/package identity check failed on {tid}; failure retained")
        idx_summary = rebuild_attempts_index(out_dir)
    finally:
        if proc is not None and proc.poll() is None: kill_and_reap_worker(proc)
        rebuild_attempts_index(out_dir)
        if proc is None or proc.poll() is not None: release_lock(lock_path, owner_id)
        signal.signal(signal.SIGTERM, old_term)
    batch_elapsed = time.monotonic() - t_batch0
    receipt = {"owner_id": owner_id, "launched_calls": launched_calls, "resolved_without_worker": resolved_no_worker, "reconciled_interrupted_attempts": reconciled_interrupted, "stopped_reason": stopped_reason,
        "batch_elapsed_seconds": round(batch_elapsed, 4), "index_summary": idx_summary}
    write_json_atomic(receipt, out_dir / "batches" / f"{owner_id}.json", overwrite=False)
    return receipt
# ---- Fake worker and self-test suite (no R, no RNG, no fits) ----
def run_fake_worker_cli(task_json_path: str):
    task = json.loads(Path(task_json_path).read_text(encoding="utf-8")); tid = task["task_id"]; call_type = task["call_type"]; out_path = Path(task["worker_out_json"])
    # Simulate behaviour based on cell_id / task_id tags in the test manifest
    cid = task["cell_id"]
    if "SLEEP_TIMEOUT" in cid:
        time.sleep(5.0); return
    if call_type == "AGG":
        p_rds = Path(task["parent_rds_path"])
        if not p_rds.exists() or p_rds.read_text(encoding="utf-8") == "CORRUPT_RDS":
            rec = {"attempt_id": task["attempt_id"], "task_id": tid, "call_type": call_type, "identity_ok": True, "disposition": "unexpected-error", "startup_seconds": 0.01, "data_seconds": 0.01,
                "data_generated_now": False, "operation_seconds": 0.0, "error": {"class": ["serialization_error"], "message": "Failed to read parent RDS"}, "warnings": [], "n_warnings": 0}
            write_json_atomic(rec, out_path, overwrite=True); return
        if "I_NO_IF" in cid:
            rec = {"attempt_id": task["attempt_id"], "task_id": tid, "call_type": call_type, "identity_ok": True,
                "disposition": "expected-refusal", "startup_seconds": 0.01, "data_seconds": 0.01, "data_generated_now": False, "operation_seconds": 0.02,
                "error": {"class": ["cdid_chained_request_no_inffunc_error", "cdid_error"], "message": "No inffunc"}, "warnings": [], "n_warnings": 0}
            write_json_atomic(rec, out_path, overwrite=True); return
        rec = {"attempt_id": task["attempt_id"], "task_id": tid, "call_type": call_type, "identity_ok": True, "disposition": "success", "startup_seconds": 0.01, "data_seconds": 0.01, "data_generated_now": False,
            "operation_seconds": 0.02, "error": None, "warnings": [], "n_warnings": 0, "agg_result": { "type": "simple", "overall_att": 1.25, "overall_se": 0.05, "overall_ci_lower": 1.152,
            "overall_ci_upper": 1.348, "egt": [], "att_egt": [], "se_egt": [], }}
        write_json_atomic(rec, out_path, overwrite=True); return
    # Simulate dataset caching once per (recipe, rep)
    d_rds = Path(task["data_rds_path"]); d_json = Path(task["data_json_path"]); gen_now = False
    if not d_rds.exists():
        gen_now = True; d_rds.write_text(f"FAKE_DATA_{task['recipe_id']}_{task['replicate']}", encoding="utf-8"); d_json.write_text(json.dumps({"panel_md5": "md5_p", "frame_md5": "md5_f"}), encoding="utf-8")
    if "PARENT_FAIL" in cid:
        rec = {"attempt_id": task["attempt_id"], "task_id": tid, "call_type": call_type, "identity_ok": True, "disposition": "unexpected-error", "startup_seconds": 0.01, "data_seconds": 0.01,
            "data_generated_now": gen_now, "operation_seconds": 0.02, "error": {"class": ["simpleError"], "message": "Simulated parent FIT error"}, "warnings": [], "n_warnings": 0}
        write_json_atomic(rec, out_path, overwrite=True); return
    if task.get("save_parent_rds") and task.get("parent_rds_path"):
        p_rds = Path(task["parent_rds_path"])
        if "CORRUPT_PARENT_RDS" in cid:
            p_rds.write_text("CORRUPT_RDS", encoding="utf-8")
        else:
            p_rds.write_text("VALID_PARENT_RDS", encoding="utf-8")
    rec = {"attempt_id": task["attempt_id"], "task_id": tid, "call_type": call_type, "identity_ok": True, "disposition": "success", "startup_seconds": 0.01, "data_seconds": 0.02 if gen_now else 0.005,
        "data_generated_now": gen_now, "operation_seconds": 0.03, "data_panel_md5": "md5_p", "data_frame_md5": "md5_f", "parent_saved": bool(task.get("save_parent_rds")), "error": None, "warnings": [], "n_warnings": 0,
        "fit_result": {"cells": { "group": [3.0], "t": [3.0], "att": [1.2], "se": [0.05], "ci_lower": [1.102], "ci_upper": [1.298], "band_lower": [1.08], "band_upper": [1.32], }, "c_val": 1.96, "bstrap": False, "cband": False},}
    write_json_atomic(rec, out_path, overwrite=True)
def run_selftest():
    """Exercise atomic records, timeouts, admission reserve, lock liveness, linked recovery, and parent/child rules."""
    fake_cmd = [sys.executable, str(Path(__file__).resolve()), "--fake-worker"]
    with tempfile.TemporaryDirectory(prefix="cdid_campaign_selftest_") as tmp:
        tmp_path = Path(tmp); manifest_path = tmp_path / "test_manifest.csv"; out_dir = tmp_path / "out"; out_dir.mkdir(parents=True, exist_ok=True)
        fieldnames = ["task_id", "replicate_scope", "stream_index", "stream_key", "recipe_id", "family", "cell_id", "logical_source", "role", "call_type", "package", "parent_task_id", "agg_kind", "literal_call",
            "rng_purpose", "primary_target_id", "primary_truth", "primary_population_masses", "primary_coordinate_set", "target_match", "causal_validity", "inference_eligible",
            "expected_disposition", "expected_condition", "secondary_reference_id", "secondary_reference", "secondary_status", "moment_regime", "equality_ref", "linked_cell_ids"]
        def make_row(tid, cid, call_type, parent_tid="", exp_disp="success", exp_cond=""):
            return {"task_id": tid, "replicate_scope": "1..1", "stream_index": 1, "stream_key": "COMP_adjacent", "recipe_id": "COMP_adjacent", "family": "COMPOSITION", "cell_id": cid, "logical_source": "targeted",
                "role": "PRIMARY-valid", "call_type": call_type, "package": "cdid", "parent_task_id": parent_tid, "agg_kind": "simple" if call_type == "AGG" else "",
                "literal_call": 'cdid::aggte(MP=FIT, type="simple")' if call_type == "AGG" else 'cdid::att_gt(yname="y")', "rng_purpose": 4 if call_type == "AGG" else 3, "primary_target_id": "ORIGINAL_COHORT_ATT",
                "primary_truth": '{"gt:3:3":1.2}', "primary_population_masses": '{"0":0.3,"3":0.25,"4":0.25,"5":0.2}', "primary_coordinate_set": "overall" if call_type == "AGG" else "gt:3:3", "target_match": "matched",
                "causal_validity": "VALID", "inference_eligible": "TRUE", "expected_disposition": exp_disp, "expected_condition": exp_cond, "secondary_reference_id": "", "secondary_reference": "",
                "secondary_status": "not_applicable", "moment_regime": "finite_fourth_moment", "equality_ref": "", "linked_cell_ids": cid,}
        rows = [make_row("CELL_OK__FIT", "CELL_OK", "FIT"), make_row("CELL_OK__AGG_S", "CELL_OK", "AGG", parent_tid="CELL_OK__FIT"), make_row("CELL_PARENT_FAIL__FIT", "CELL_PARENT_FAIL", "FIT"),
            make_row("CELL_PARENT_FAIL__AGG_S", "CELL_PARENT_FAIL", "AGG", parent_tid="CELL_PARENT_FAIL__FIT"), make_row("CELL_CORRUPT_PARENT_RDS__FIT", "CELL_CORRUPT_PARENT_RDS", "FIT"),
            make_row("CELL_CORRUPT_PARENT_RDS__AGG_S", "CELL_CORRUPT_PARENT_RDS", "AGG", parent_tid="CELL_CORRUPT_PARENT_RDS__FIT"), make_row("CELL_I_NO_IF__FIT", "CELL_I_NO_IF", "FIT"),
            make_row("CELL_I_NO_IF__AGG_S", "CELL_I_NO_IF", "AGG", parent_tid="CELL_I_NO_IF__FIT", exp_disp="expected-refusal", exp_cond="cdid_chained_request_no_inffunc_error"),
            make_row("CELL_SLEEP_TIMEOUT__FIT", "CELL_SLEEP_TIMEOUT", "FIT"),]
        with manifest_path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames); w.writeheader(); w.writerows(rows)
        targets_tasks = {}
        for r in rows:
            targets_tasks[r["task_id"]] = {**r, "primary_truth_map": {"overall": 1.25} if r["call_type"] == "AGG" else {"gt:3:3": 1.2},
                "auxiliary_truth_map": {}, "secondary_reference_map": {}, "population_masses_map": {"0": 0.3, "3": 0.25, "4": 0.25, "5": 0.2}, "save_parent_rds": r["call_type"] == "FIT" and r["cell_id"] != "CELL_SLEEP_TIMEOUT",}
        (out_dir / "targets.json").write_text(json.dumps({"tasks": targets_tasks}, indent=2), encoding="utf-8")
        # 1. Test alive same-host lock refusal
        write_json_atomic({"owner_id": "alive_owner", "hostname": socket.gethostname(), "supervisor_pid": os.getpid(), "worker_pid": None, "worker_pgid": None }, out_dir / "lock.json")
        try:
            run_batch("", manifest_path, out_dir, allow_unfrozen=True, fake_worker_cmd=fake_cmd); raise AssertionError("Expected alive same-host lock to refuse acquisition")
        except RuntimeError as e:
            assert "not proven dead" in str(e)
        # 2. Test cross-host lock refusal and explicit operator reconciliation
        write_json_atomic({"owner_id": "remote_owner_99", "hostname": "remote-cluster-node-99.internal", "supervisor_pid": 999999, "worker_pid": None, "worker_pgid": None }, out_dir / "lock.json")
        try:
            run_batch("", manifest_path, out_dir, allow_unfrozen=True, fake_worker_cmd=fake_cmd); raise AssertionError("Expected cross-host lock to refuse acquisition without --reconcile-cross-host")
        except RuntimeError as e:
            assert "held by another host" in str(e)
        # 3. Simulate dead same-host lock + open interrupted attempt on CELL_OK__FIT
        records_dir = out_dir / "records"; records_dir.mkdir(parents=True, exist_ok=True)
        write_json_atomic({"attempt_id": "r001__CELL_OK__FIT__att01", "task_id": "r001__CELL_OK__FIT", "template_task_id": "CELL_OK__FIT", "attempt_number": 1, "linked_previous_attempt_id": "",
            "owner_id": "dead_owner_1", "replicate": 1, "stream_index": 1, "recipe_id": "COMP_adjacent", "cell_id": "CELL_OK", "call_type": "FIT", "pop_rng_index": 1, "obs_rng_index": 2, "op_rng_index": 3
            }, records_dir / "r001__CELL_OK__FIT__att01.open.json")
        write_json_atomic({
            "owner_id": "dead_owner_1", "hostname": socket.gethostname(), "supervisor_pid": 99999999,
            # non-existent PID "worker_pid": 99999998, "worker_pgid": 99999998
        }, out_dir / "lock.json")
        # 4. Test admission reserve refusal when wall_seconds < op_timeout + cleanup_reserve
        res_admit = run_batch("", manifest_path, out_dir, max_calls=25, wall_seconds=1.0, op_timeout=0.8, cleanup_reserve=0.5, allow_unfrozen=True, fake_worker_cmd=fake_cmd)
        assert res_admit["launched_calls"] == 0; assert res_admit["stopped_reason"] == "admission_reserve_reached"; assert "r001__CELL_OK__FIT__att01" in res_admit["reconciled_interrupted_attempts"]
        # Verify att01 was marked interrupted
        att01_rec = json.loads((records_dir / "r001__CELL_OK__FIT__att01.result.json").read_text(encoding="utf-8")); assert att01_rec["disposition"] == "interrupted"
        # 5. Run batch with sufficient budget (op_timeout=0.3s so SLEEP_TIMEOUT times out quickly)
        res_run = run_batch("", manifest_path, out_dir, max_calls=25, wall_seconds=30.0, op_timeout=0.3, cleanup_reserve=0.2, allow_unfrozen=True, fake_worker_cmd=fake_cmd); assert res_run["stopped_reason"] == "all_pending_completed"
        assert res_run["resolved_without_worker"] == 1  # CELL_PARENT_FAIL__AGG_S
        # Verify linked replacement attempt r001__CELL_OK__FIT__att02 succeeded and linked to att01
        att02_rec = json.loads((records_dir / "r001__CELL_OK__FIT__att02.result.json").read_text(encoding="utf-8")); assert att02_rec["disposition"] == "success"; assert att02_rec["attempt_number"] == 2
        assert att02_rec["linked_previous_attempt_id"] == "r001__CELL_OK__FIT__att01"; assert att02_rec["pop_rng_index"] == att01_rec["pop_rng_index"]
        # Verify parent failure -> child parent-unavailable without worker
        pfail_agg = json.loads((records_dir / "r001__CELL_PARENT_FAIL__AGG_S__att01.result.json").read_text(encoding="utf-8")); assert pfail_agg["disposition"] == "parent-unavailable"; assert pfail_agg["worker_launched"] is False
        # Verify corrupt parent RDS -> child unexpected-error without refitting parent
        corrupt_agg = json.loads((records_dir / "r001__CELL_CORRUPT_PARENT_RDS__AGG_S__att01.result.json").read_text(encoding="utf-8"))
        assert corrupt_agg["disposition"] == "unexpected-error"; assert "serialization_error" in corrupt_agg["error"]["class"]; assert not (records_dir / "r001__CELL_CORRUPT_PARENT_RDS__FIT__att02.result.json").exists()
        # Verify I_NO_IF child expected-refusal
        noif_agg = json.loads((records_dir / "r001__CELL_I_NO_IF__AGG_S__att01.result.json").read_text(encoding="utf-8")); assert noif_agg["disposition"] == "expected-refusal"
        # Verify timeout handling and no retry on second batch invocation
        to_rec = json.loads((records_dir / "r001__CELL_SLEEP_TIMEOUT__FIT__att01.result.json").read_text(encoding="utf-8")); assert to_rec["disposition"] == "timeout"; assert to_rec["censored_timeout"] is True
        res_rerun = run_batch("", manifest_path, out_dir, max_calls=25, wall_seconds=30.0, op_timeout=0.3, cleanup_reserve=0.2, allow_unfrozen=True, fake_worker_cmd=fake_cmd)
        assert res_rerun["launched_calls"] == 0, "Terminal timeout/error/success tasks must never be retried"
        # 6. Verify provenance hash mismatch refusal on resume
        write_json_atomic({"source_sha": "abc123456789", "harness_hashes": {"manifest.csv": "0000000000000000000000000000000000000000000000000000000000000000"}, }, out_dir / "provenance.json")
        (out_dir / "manifest.csv").write_bytes(manifest_path.read_bytes())
        try:
            run_batch("", manifest_path, out_dir, allow_unfrozen=False, fake_worker_cmd=fake_cmd); raise AssertionError("Expected provenance hash mismatch to refuse resume")
        except RuntimeError as e:
            assert "hash mismatch" in str(e)
    return {"status": "SELFTEST_OK", "verified": [ "alive_same_host_lock_refusal", "cross_host_lock_refusal", "proven_dead_lock_reclamation", "interrupted_attempt_marking_and_linked_replacement",
        "admission_budget_cleanup_reserve", "parent_failure_to_parent_unavailable_without_worker", "child_serialization_failure_without_parent_refit", "expected_refusal_classification_I_NO_IF",
        "per_operation_timeout_process_group_kill_and_no_retry", "provenance_mismatch_resume_refusal", ]}
def main():
    parser = argparse.ArgumentParser(description="Run one bounded reliability campaign batch and exit.")
    parser.add_argument("--repo", default="", help="Path to source checkout."); parser.add_argument("--manifest", default="", help="Path to manifest.csv.")
    parser.add_argument("--out", default="", help="Path to campaign results directory."); parser.add_argument("--max-calls", type=int, default=25, help="Max launched operations in this batch (default 25).")
    parser.add_argument("--wall-seconds", type=float, default=1200.0, help="Total batch wall-clock cap in seconds (default 1200).")
    parser.add_argument("--op-timeout", type=float, default=300.0, help="Per-operation wall timeout in seconds (default 300).")
    parser.add_argument("--cleanup-reserve", type=float, default=30.0, help="Cleanup reserve in seconds (default 30).")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable."); parser.add_argument("--lib", action="append", default=[], help="R library path(s).")
    parser.add_argument("--reconcile-cross-host", default="", help="Explicit owner_id to reconcile after operator confirms remote host is stopped.")
    parser.add_argument("--max-reps", type=int, default=0, help="Optional replicate cap for bounded smoke runs."); parser.add_argument("--selftest", action="store_true", help="Run fake-worker unit test suite (no R, no RNG, no fits).")
    parser.add_argument("--fake-worker", action="store_true", help=argparse.SUPPRESS); parser.add_argument("--task", default="", help=argparse.SUPPRESS); args = parser.parse_args()
    if args.fake_worker:
        run_fake_worker_cli(args.task); return
    if args.selftest:
        res = run_selftest(); print(json.dumps(res, indent=2)); return
    if not args.manifest or not args.out:
        raise SystemExit("Usage: run_batch.py --repo CHECKOUT --manifest manifest.csv --out RESULTS [--max-calls 25] [--wall-seconds 1200]")
    summary = run_batch(repo_dir=args.repo, manifest_path=Path(args.manifest).resolve(), out_dir=Path(args.out).resolve(), max_calls=args.max_calls, wall_seconds=args.wall_seconds, op_timeout=args.op_timeout,
        cleanup_reserve=args.cleanup_reserve, rscript=args.rscript, lib_paths=args.lib, reconcile_cross_host=args.reconcile_cross_host, max_reps_override=args.max_reps)
    print(json.dumps(summary, indent=2))
if __name__ == "__main__":
    main()
