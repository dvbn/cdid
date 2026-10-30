#!/usr/bin/env python3
"""Portable reliability campaign: tidy CSV and Markdown summary generator.
Reads immutable disk attempt records (.result.json), manifest.csv, and targets.json.
Never modifies campaign state or executes R. Reports primary causal errors, bias,
MCSE (where moments permit), empirical SD, RMSE, mean SE, SE/SD ratio, interval; width, conditional coverage and all-attempt containment with exact 95% Clopper-Pearson
intervals, full-family simultaneous containment over all frozen coordinates,
separate secondary reference diagnostics, paired method differences on common; completed draws, and mean/median runtimes separated into successful latency vs; all-attempt cost vs setup/data cost vs right-censored 300s timeouts.
"""
import argparse, csv, datetime, functools, json, math, statistics, tempfile; from pathlib import Path
@functools.lru_cache(maxsize=None)
def cp_interval(k: int, n: int, alpha: float = 0.05):
    """Exact 95% Clopper-Pearson (stats::binom.test) interval."""
    if n <= 0: return (None, None)
    if not 0 <= k <= n or not 0 < alpha < 1: raise ValueError("Invalid binomial parameters")
    cdf = lambda p, kk: sum(math.comb(n, i) * (p ** i) * ((1.0 - p) ** (n - i)) for i in range(kk + 1))
    def solve(fn):
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if fn(mid) < 0.0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)
    lower = 0.0 if k <= 0 else solve(lambda p: (1.0 - cdf(p, k - 1)) - 0.5 * alpha); upper = 1.0 if k >= n else solve(lambda p: 0.5 * alpha - cdf(p, k)); return (lower, upper)
def fmt_val(x, digits: int = 6) -> str:
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))): return ""
    if isinstance(x, bool): return "TRUE" if x else "FALSE"
    if isinstance(x, int): return str(x)
    if isinstance(x, float): return f"{x:.{digits}f}"
    return str(x)
def fmt_num_key(v: float) -> str:
    if float(v).is_integer(): return str(int(v))
    return f"{v:.6g}"
def mean_or_none(vals: list):
    return statistics.fmean(vals) if vals else None
def median_or_none(vals: list):
    return statistics.median(vals) if vals else None
def sd_or_none(vals: list):
    return statistics.stdev(vals) if len(vals) >= 2 else None
def rmse_or_none(errs: list):
    return math.sqrt(statistics.fmean([e * e for e in errs])) if errs else None
def measured_setup_cost(rec: dict):
    """Return measured setup/data cost; incomplete crash receipts remain unknown."""
    if not rec.get("worker_launched", True): return None
    parts = [rec.get(k) for k in ("supervisor_overhead_seconds", "startup_seconds", "data_seconds")]
    if any(v is None for v in parts): return None
    return sum(float(v) for v in parts)
def moment_qualification_note(moment_regime: str) -> str:
    if moment_regime == "finite_variance_only_t3":
        return "t3_finite_variance_only: bias_mcse_valid; no_fourth_moment_rmse_or_sd_uncertainty_claim"
    if moment_regime == "infinite_variance_t1_5":
        return "t1.5_infinite_variance: cautious_finite_sample_descriptive_only; moment_mcse_undefined"
    if moment_regime == "zero_variance_deterministic": return "zero_noise_deterministic_scale"
    return "finite_fourth_moment"
def extract_fit_coord_map(rec: dict):
    """Return mapping coord_key -> {est, se, se_analytical, ci_lo, ci_hi, band_lo, band_hi}."""
    fr = rec.get("fit_result")
    if not fr or not isinstance(fr, dict): return {}
    cells = fr.get("cells") or {}; groups = cells.get("group") or []; ts = cells.get("t") or []; atts = cells.get("att") or []; ses = cells.get("se") or []; se_an = cells.get("se_analytical") or []
    ci_los = cells.get("ci_lower") or []; ci_his = cells.get("ci_upper") or []; b_los = cells.get("band_lower") or []; b_his = cells.get("band_upper") or []; use_band = bool(fr.get("cband")); out = {}
    for i in range(min(len(groups), len(ts), len(atts))):
        g_v, t_v = groups[i], ts[i]
        if g_v is None or t_v is None: continue
        ck = f"gt:{fmt_num_key(float(g_v))}:{fmt_num_key(float(t_v))}"; est = atts[i] if i < len(atts) else None; se = ses[i] if i < len(ses) else None; se_a = se_an[i] if i < len(se_an) else None
        ci_l = ci_los[i] if i < len(ci_los) else None; ci_u = ci_his[i] if i < len(ci_his) else None; bd_l = b_los[i] if i < len(b_los) else ci_l; bd_u = b_his[i] if i < len(b_his) else ci_u
        out[ck] = {"est": float(est) if (est is not None and math.isfinite(float(est))) else None, "se": float(se) if (se is not None and math.isfinite(float(se))) else None,
            "se_analytical": float(se_a) if (se_a is not None and math.isfinite(float(se_a))) else None, "ci_lo": float(ci_l) if (ci_l is not None and math.isfinite(float(ci_l))) else None,
            "ci_hi": float(ci_u) if (ci_u is not None and math.isfinite(float(ci_u))) else None, "fam_lo": float(bd_l) if (use_band and bd_l is not None and math.isfinite(float(bd_l))) else (
            float(ci_l) if (ci_l is not None and math.isfinite(float(ci_l))) else None), "fam_hi": float(bd_u) if (use_band and bd_u is not None and math.isfinite(float(bd_u))) else (
            float(ci_u) if (ci_u is not None and math.isfinite(float(ci_u))) else None),}
    return out
def extract_agg_coord_map(rec: dict):
    ar = rec.get("agg_result")
    if not ar or not isinstance(ar, dict): return {}, False
    out = {}; ov_est = ar.get("overall_att"); ov_se = ar.get("overall_se"); ov_lo = ar.get("overall_ci_lower"); ov_hi = ar.get("overall_ci_upper")
    out["overall"] = {"est": float(ov_est) if (ov_est is not None and math.isfinite(float(ov_est))) else None, "se": float(ov_se) if (ov_se is not None and math.isfinite(float(ov_se))) else None, "se_analytical": None,
        "ci_lo": float(ov_lo) if (ov_lo is not None and math.isfinite(float(ov_lo))) else None, "ci_hi": float(ov_hi) if (ov_hi is not None and math.isfinite(float(ov_hi))) else None,
        "fam_lo": float(ov_lo) if (ov_lo is not None and math.isfinite(float(ov_lo))) else None, "fam_hi": float(ov_hi) if (ov_hi is not None and math.isfinite(float(ov_hi))) else None,}
    egts = ar.get("egt") or []; att_egts = ar.get("att_egt") or []; se_egts = ar.get("se_egt") or []; ci_los = ar.get("ci_lower_egt") or []; ci_his = ar.get("ci_upper_egt") or []
    b_los = ar.get("band_lower_egt") or []; b_his = ar.get("band_upper_egt") or []; use_band = bool(ar.get("cband"))
    for i in range(min(len(egts), len(att_egts))):
        ev = egts[i]
        if ev is None: continue
        ck = f"egt:{fmt_num_key(float(ev))}"; est = att_egts[i]; se = se_egts[i] if i < len(se_egts) else None; ci_l = ci_los[i] if i < len(ci_los) else None; ci_u = ci_his[i] if i < len(ci_his) else None
        bd_l = b_los[i] if i < len(b_los) else ci_l; bd_u = b_his[i] if i < len(b_his) else ci_u
        out[ck] = {"est": float(est) if (est is not None and math.isfinite(float(est))) else None, "se": float(se) if (se is not None and math.isfinite(float(se))) else None, "se_analytical": None,
            "ci_lo": float(ci_l) if (ci_l is not None and math.isfinite(float(ci_l))) else None, "ci_hi": float(ci_u) if (ci_u is not None and math.isfinite(float(ci_u))) else None,
            "fam_lo": float(bd_l) if (use_band and bd_l is not None and math.isfinite(float(bd_l))) else (float(ci_l) if (ci_l is not None and math.isfinite(float(ci_l))) else None),
            "fam_hi": float(bd_u) if (use_band and bd_u is not None and math.isfinite(float(bd_u))) else (float(ci_u) if (ci_u is not None and math.isfinite(float(ci_u))) else None),}
    excluded = ar.get("excluded_cells") or []; support_changed = bool(excluded) or bool(ar.get("support_changed")); return out, support_changed
def summarize_campaign(manifest_path: Path, targets_path: Path, out_dir: Path, n_scheduled_reps: int = 100):
    manifest_rows = list(csv.DictReader(manifest_path.open(encoding="utf-8"))); targets_payload = json.loads(targets_path.read_text(encoding="utf-8"))
    tasks_meta = targets_payload.get("tasks", {}); records_dir = out_dir / "records"; raw_attempts = []
    if records_dir.exists():
        for p in sorted(records_dir.glob("*.result.json")):
            raw_attempts.append(json.loads(p.read_text(encoding="utf-8")))
    if not 1 <= n_scheduled_reps <= 100: raise ValueError("Scheduled denominator must be between 1 and 100")
    if any(int(r.get("replicate", 1)) > n_scheduled_reps for r in raw_attempts):
        raise ValueError("Summary denominator excludes recorded draws")
    # Group attempts by (template_task_id, replicate)
    attempts_by_tpl_rep = {}; attempts_by_tpl = {}
    for rec in raw_attempts:
        tpl_id = rec.get("template_task_id") or rec["task_id"].split("__", 1)[-1]; rep = int(rec.get("replicate", 1)); attempts_by_tpl_rep.setdefault((tpl_id, rep), []).append(rec); attempts_by_tpl.setdefault(tpl_id, []).append(rec)
    for key in attempts_by_tpl_rep:
        attempts_by_tpl_rep[key].sort(key=lambda r: int(r.get("attempt_number", 1)))
    # Terminal record per (template_task_id, replicate)
    terminal_by_tpl_rep = {}
    for (tpl_id, rep), recs in attempts_by_tpl_rep.items():
        for r in reversed(recs):
            if r.get("disposition") != "interrupted":
                terminal_by_tpl_rep[(tpl_id, rep)] = r; break
    # 1. Build attempts.csv rows
    attempts_csv_rows = []
    for rec in sorted(raw_attempts, key=lambda r: (int(r.get("replicate", 1)), r.get("task_id", ""), int(r.get("attempt_number", 1)))):
        err = rec.get("error") or {}; err_cls = ";".join(err.get("class", [])) if isinstance(err.get("class"), list) else str(err.get("class", ""))
        attempts_csv_rows.append({"attempt_id": rec.get("attempt_id", ""), "task_id": rec.get("task_id", ""), "template_task_id": rec.get("template_task_id", ""), "attempt_number": rec.get("attempt_number", 1),
            "linked_previous_attempt_id": rec.get("linked_previous_attempt_id", ""), "replicate": rec.get("replicate", 1), "stream_index": rec.get("stream_index", ""),
            "recipe_id": rec.get("recipe_id", ""), "family": rec.get("family", ""), "cell_id": rec.get("cell_id", ""), "call_type": rec.get("call_type", ""), "package": rec.get("package", ""),
            "disposition": rec.get("disposition", ""), "worker_launched": fmt_val(bool(rec.get("worker_launched", True))), "censored_timeout": fmt_val(bool(rec.get("censored_timeout", False))),
            "supervisor_overhead_seconds": fmt_val(rec.get("supervisor_overhead_seconds")), "startup_seconds": fmt_val(rec.get("startup_seconds")), "data_seconds": fmt_val(rec.get("data_seconds")),
            "operation_seconds": fmt_val(rec.get("operation_seconds")), "total_attempt_seconds": fmt_val(rec.get("total_attempt_seconds")), "n_warnings": rec.get("n_warnings", 0),
            "error_class": err_cls, "error_message": err.get("message", "") if isinstance(err, dict) else "", })
    # 2. Build cells.csv, aggregates.csv, pretests.csv
    cells_csv_rows = []; aggregates_csv_rows = []; pretests_csv_rows = []
    for mrow in manifest_rows:
        tpl_id = mrow["task_id"]; call_type = mrow["call_type"]; rec_id = mrow["recipe_id"]; cid = mrow["cell_id"]; tmeta = tasks_meta.get(tpl_id, {}); m_regime = mrow.get("moment_regime", "finite_fourth_moment")
        m_note = moment_qualification_note(m_regime); inf_elig = mrow.get("inference_eligible", "TRUE") == "TRUE"; infinite_var = (m_regime == "infinite_variance_t1_5"); all_att_recs = attempts_by_tpl.get(tpl_id, [])
        term_recs = [terminal_by_tpl_rep[(tpl_id, rep)] for rep in range(1, n_scheduled_reps + 1) if (tpl_id, rep) in terminal_by_tpl_rep]
        # Disposition counts across attempts and terminal records
        n_completed = len(term_recs); n_success = sum(1 for r in term_recs if r.get("disposition") == "success"); n_partial = sum(1 for r in term_recs if r.get("disposition") == "partial-support")
        n_refusals = sum(1 for r in term_recs if r.get("disposition") == "expected-refusal"); n_unexp_err = sum(1 for r in term_recs if r.get("disposition") == "unexpected-error")
        n_unexp_succ = sum(1 for r in term_recs if r.get("disposition") == "unexpected-success"); n_timeouts = sum(1 for r in term_recs if r.get("disposition") == "timeout")
        n_parent_unavail = sum(1 for r in term_recs if r.get("disposition") == "parent-unavailable"); n_unavail_comp = sum(1 for r in term_recs if r.get("disposition") == "unavailable-comparator")
        n_not_app = sum(1 for r in term_recs if r.get("disposition") == "not-applicable"); n_interrupted = sum(1 for r in all_att_recs if r.get("disposition") == "interrupted")
        n_warn_draws = sum(1 for r in term_recs if int(r.get("n_warnings", 0) or 0) > 0)
        # Runtime separation: successful operation latency vs all-attempt cost vs setup/data cost vs timeouts
        succ_op_secs = [float(r["operation_seconds"]) for r in term_recs if r.get("disposition") == "success" and r.get("operation_seconds") is not None]
        all_op_secs = [float(r["operation_seconds"]) for r in all_att_recs if r.get("operation_seconds") is not None]
        all_tot_secs = [float(r["total_attempt_seconds"]) for r in all_att_recs if r.get("total_attempt_seconds") is not None]
        setup_data_secs = [v for r in all_att_recs if (v := measured_setup_cost(r)) is not None]
        rt_fields = {"mean_success_op_seconds": fmt_val(mean_or_none(succ_op_secs)), "median_success_op_seconds": fmt_val(median_or_none(succ_op_secs)),
            "mean_all_attempt_op_seconds": fmt_val(mean_or_none(all_op_secs)), "median_all_attempt_op_seconds": fmt_val(median_or_none(all_op_secs)), "total_all_attempt_seconds": fmt_val(sum(all_tot_secs) if all_tot_secs else 0.0),
            "mean_all_attempt_total_seconds": fmt_val(mean_or_none(all_tot_secs)), "median_all_attempt_total_seconds": fmt_val(median_or_none(all_tot_secs)),
            "median_setup_and_data_seconds": fmt_val(median_or_none(setup_data_secs)), "n_unknown_attempt_costs": sum(r.get("total_attempt_seconds") is None for r in all_att_recs),
            "mean_setup_and_data_seconds": fmt_val(mean_or_none(setup_data_secs)), "n_censored_timeouts_300s": n_timeouts,}
        if call_type in ("FIT", "DID_CC", "DID_UNBALANCED") and inf_elig:
            # Package Wald fields are separate from the conditional CvM workflow.
            vals = [(r.get("fit_result") or {}).get("pretest") or {} for r in term_recs]
            pvals = [v["Wpval"] for v in vals if v.get("Wpval") is not None and math.isfinite(v["Wpval"])]
            stats = [v["W"] for v in vals if v.get("W") is not None and math.isfinite(v["W"])]
            n_avail = len(pvals); n_rej = sum(v < .05 for v in pvals)
            ci = cp_interval(n_rej,n_avail); aa = cp_interval(n_rej,n_scheduled_reps)
            aux = tmeta.get("auxiliary_truth_map", {})
            pretests_csv_rows.append({"test_type":"package_Wald", "task_id":tpl_id, "recipe_id":rec_id,
                "family":mrow["family"], "cell_id":cid, "role":mrow["role"],
                "conditional_pt_null_true":"", "unconditional_mean_zero":fmt_val(all(v==0 for v in aux.values())) if aux else "",
                "n_scheduled":n_scheduled_reps,"n_completed":n_completed,"n_pending":n_scheduled_reps-n_completed,
                "complete_100_draws":n_completed==100,"n_available":n_avail,"n_test_unavailable":n_completed-n_avail,
                "n_not_applicable":sum(v.get("status")=="not_applicable" for v in vals),
                "n_refusals":n_refusals,"n_partial_support":n_partial,"n_errors":n_unexp_err,
                "n_timeouts":n_timeouts,"n_interrupted":n_interrupted,"n_parent_unavailable":n_parent_unavail,
                "n_unavailable_comparator":n_unavail_comp,"n_warning_draws":n_warn_draws,
                "n_reject_005":n_rej,"rejection_rate_005":fmt_val(n_rej/n_avail if n_avail else None),
                "rej_ci_low":fmt_val(ci[0]),"rej_ci_high":fmt_val(ci[1]),
                "all_attempt_rejection_rate":fmt_val(n_rej/n_scheduled_reps),"all_rej_ci_low":fmt_val(aa[0]),"all_rej_ci_high":fmt_val(aa[1]),
                "mean_W":fmt_val(mean_or_none(stats)),"median_pval":fmt_val(median_or_none(pvals)),
                "causal_validity":mrow["causal_validity"], "target_match":mrow["target_match"], **rt_fields})
        if call_type == "PRETEST":
            p_truth_map = tmeta.get("primary_truth_map", {}); pvals = []; cvm_stats = []
            for r in term_recs:
                pr = r.get("pretest_result") or {}; pv = pr.get("CvMpval"); st = pr.get("CvM")
                if pv is not None and math.isfinite(float(pv)):
                    pvals.append(float(pv))
                if st is not None and math.isfinite(float(st)):
                    cvm_stats.append(float(st))
            n_avail = len(pvals); n_rej = sum(1 for pv in pvals if pv < 0.05); rej_rate = (n_rej / n_avail) if n_avail > 0 else None; rej_lo, rej_hi = cp_interval(n_rej, n_avail) if n_avail > 0 else (None, None)
            aa_lo, aa_hi = cp_interval(n_rej, n_scheduled_reps)
            pretests_csv_rows.append({"test_type": "conditional_CvM", "task_id": tpl_id, "recipe_id": rec_id, "family": mrow["family"], "cell_id": cid, "role": mrow["role"], "conditional_pt_null_true": fmt_val(bool(p_truth_map.get("conditional_pt_null_true", True))),
                "unconditional_mean_zero": fmt_val(bool(p_truth_map.get("unconditional_mean_zero", True))), "n_scheduled": n_scheduled_reps, "n_completed": n_completed, "n_available": n_avail,
                "n_pending": n_scheduled_reps - n_completed, "complete_100_draws": n_completed == 100, "n_interrupted": n_interrupted, "n_refusals": n_refusals, "n_partial_support": n_partial, "n_unavailable_comparator": n_unavail_comp, "n_parent_unavailable": n_parent_unavail, "n_warning_draws": n_warn_draws,
                "all_attempt_rejection_rate": fmt_val(n_rej / n_scheduled_reps), "all_rej_ci_low": fmt_val(aa_lo), "all_rej_ci_high": fmt_val(aa_hi), "n_not_applicable": n_not_app, "n_errors": n_unexp_err, "n_timeouts": n_timeouts, "n_reject_005": n_rej, "rejection_rate_005": fmt_val(rej_rate), "rej_ci_low": fmt_val(rej_lo),
                "rej_ci_high": fmt_val(rej_hi), "mean_cvm_stat": fmt_val(mean_or_none(cvm_stats)), "median_cvm_pval": fmt_val(median_or_none(pvals)), **rt_fields })
            continue
        # Extract per-rep coordinate maps for FIT, DID_*, or AGG
        rep_coord_maps = {}; n_support_changed = 0
        for r in term_recs:
            rep = int(r.get("replicate", 1))
            if call_type == "AGG":
                cmap, sup_chg = extract_agg_coord_map(r)
                if sup_chg:
                    n_support_changed += 1
            else:
                cmap = extract_fit_coord_map(r)
            if cmap:
                rep_coord_maps[rep] = cmap
        p_truth_map = tmeta.get("primary_truth_map", {}); aux_truth_map = tmeta.get("auxiliary_truth_map", {}); sec_ref_map = tmeta.get("secondary_reference_map", {})
        p_coord_keys = [k for k in mrow.get("primary_coordinate_set", "").split(";") if k]
        def build_coord_summary_row(coord_key: str, coord_role: str, truth_val: float, sec_val):
            ests, errs, ses, ses_an, widths = [], [], [], [], []; n_cov_cond = 0; n_avail_int = 0; sec_errs = []; n_sec_cov = 0
            for rep, cmap in rep_coord_maps.items():
                cinfo = cmap.get(coord_key)
                if not cinfo or cinfo["est"] is None: continue
                est = cinfo["est"]; err = est - truth_val; ests.append(est); errs.append(err)
                if sec_val is not None:
                    sec_errs.append(est - sec_val)
                if cinfo["se"] is not None:
                    ses.append(cinfo["se"])
                if cinfo["se_analytical"] is not None:
                    ses_an.append(cinfo["se_analytical"])
                if inf_elig and cinfo["ci_lo"] is not None and cinfo["ci_hi"] is not None:
                    lo, hi = cinfo["ci_lo"], cinfo["ci_hi"]; widths.append(hi - lo); n_avail_int += 1
                    if lo <= truth_val <= hi:
                        n_cov_cond += 1
                    if sec_val is not None and (lo <= sec_val <= hi):
                        n_sec_cov += 1
            n_pts = len(ests); emp_sd = sd_or_none(errs); bias_val = mean_or_none(errs); mcse_val = (emp_sd / math.sqrt(n_pts)) if (emp_sd is not None and n_pts >= 2 and not infinite_var) else None
            rmse_val = rmse_or_none(errs); mean_se_val = mean_or_none(ses) if inf_elig else None; se_sd_val = (mean_se_val / emp_sd) if (inf_elig and mean_se_val is not None and emp_sd and emp_sd > 0) else None
            mean_se_an_val = mean_or_none(ses_an) if inf_elig else None; se_sd_an_val = (mean_se_an_val / emp_sd) if (inf_elig and mean_se_an_val is not None and emp_sd and emp_sd > 0) else None
            mean_w_val = mean_or_none(widths) if inf_elig else None
            if inf_elig and n_avail_int > 0:
                cond_cov = n_cov_cond / n_avail_int; cc_lo, cc_hi = cp_interval(n_cov_cond, n_avail_int)
            else:
                cond_cov, cc_lo, cc_hi = (None, None, None)
            if inf_elig:
                all_att_cov = n_cov_cond / n_scheduled_reps; aa_lo, aa_hi = cp_interval(n_cov_cond, n_scheduled_reps)
            else:
                all_att_cov, aa_lo, aa_hi = (None, None, None)
            sec_bias = mean_or_none(sec_errs) if sec_errs else None; sec_rmse = rmse_or_none(sec_errs) if sec_errs else None; sec_cov = (n_sec_cov / n_avail_int) if (inf_elig and sec_val is not None and n_avail_int > 0) else None
            sec_ci = cp_interval(n_sec_cov, n_avail_int) if sec_cov is not None else (None,None)
            base_row = {"task_id": tpl_id, "recipe_id": rec_id, "family": mrow["family"], "cell_id": cid, "linked_cell_ids": mrow.get("linked_cell_ids", cid), "call_type": call_type, "package": mrow["package"],
                "agg_kind": mrow.get("agg_kind", ""), "coordinate": coord_key, "coordinate_role": coord_role, "primary_target_id": mrow["primary_target_id"], "primary_truth": fmt_val(truth_val), "target_match": mrow["target_match"],
                "causal_validity": mrow["causal_validity"], "inference_eligible": "TRUE" if inf_elig else "FALSE", "moment_regime": m_regime,
                "moment_note": m_note, "n_scheduled": n_scheduled_reps, "n_pending": n_scheduled_reps-n_completed, "complete_100_draws": n_completed == 100, "n_attempts_including_interrupted": len(all_att_recs), "n_completed_attempts": n_completed, "n_available_points": n_pts,
                "n_available_intervals": n_avail_int if inf_elig else 0, "n_success": n_success, "n_partial_support": n_partial,
                "n_refusals": n_refusals, "n_unexpected_errors": n_unexp_err, "n_unexpected_success": n_unexp_succ, "n_timeouts": n_timeouts,
                "n_parent_unavailable": n_parent_unavail, "n_unavailable_comparator": n_unavail_comp, "n_interrupted": n_interrupted, "n_warning_draws": n_warn_draws, "n_support_changed": n_support_changed if call_type == "AGG" else 0,
                "mean_estimate": fmt_val(mean_or_none(ests)), "bias": fmt_val(bias_val), "mcse_bias": fmt_val(mcse_val), "empirical_sd": fmt_val(emp_sd), "rmse": fmt_val(rmse_val), "mean_se": fmt_val(mean_se_val),
                "se_sd_ratio": fmt_val(se_sd_val), "mean_se_analytical": fmt_val(mean_se_an_val), "se_sd_ratio_analytical": fmt_val(se_sd_an_val), "mean_width": fmt_val(mean_w_val), "conditional_coverage": fmt_val(cond_cov),
                "cond_cov_ci_low": fmt_val(cc_lo), "cond_cov_ci_high": fmt_val(cc_hi), "all_attempt_containment": fmt_val(all_att_cov), "all_att_ci_low": fmt_val(aa_lo), "all_att_ci_high": fmt_val(aa_hi),
                "secondary_reference_id": mrow.get("secondary_reference_id", ""), "secondary_status": mrow.get("secondary_status", ""),
                "secondary_reference_value": fmt_val(sec_val), "secondary_bias": fmt_val(sec_bias), "secondary_rmse": fmt_val(sec_rmse), "secondary_conditional_coverage": fmt_val(sec_cov), "secondary_cond_ci_low": fmt_val(sec_ci[0]), "secondary_cond_ci_high": fmt_val(sec_ci[1]), **rt_fields,}
            return base_row
        target_list = aggregates_csv_rows if call_type == "AGG" else cells_csv_rows
        for ck in p_coord_keys:
            if ck not in p_truth_map: continue
            t_val = float(p_truth_map[ck]); s_val = float(sec_ref_map[ck]) if ck in sec_ref_map and sec_ref_map[ck] is not None else None; target_list.append(build_coord_summary_row(ck, "primary_causal", t_val, s_val))
        # Full-family simultaneous containment row over all frozen primary coordinates
        if p_coord_keys:
            n_full_elig = 0; n_full_cov = 0
            for rep, cmap in rep_coord_maps.items():
                if not inf_elig: continue
                all_present = True; all_contained = True
                for ck in p_coord_keys:
                    cinfo = cmap.get(ck); t_val = float(p_truth_map[ck])
                    if not cinfo or cinfo["fam_lo"] is None or cinfo["fam_hi"] is None:
                        all_present = False; all_contained = False; break
                    if not (cinfo["fam_lo"] <= t_val <= cinfo["fam_hi"]):
                        all_contained = False
                if all_present:
                    n_full_elig += 1
                    if all_contained:
                        n_full_cov += 1
            ff_cond = (n_full_cov / n_full_elig) if (inf_elig and n_full_elig > 0) else None; ff_cc_lo, ff_cc_hi = cp_interval(n_full_cov, n_full_elig) if (inf_elig and n_full_elig > 0) else (None, None)
            ff_all = (n_full_cov / n_scheduled_reps) if inf_elig else None; ff_aa_lo, ff_aa_hi = cp_interval(n_full_cov, n_scheduled_reps) if inf_elig else (None, None)
            target_list.append({"task_id": tpl_id, "recipe_id": rec_id, "family": mrow["family"], "cell_id": cid, "linked_cell_ids": mrow.get("linked_cell_ids", cid), "call_type": call_type, "package": mrow["package"],
                "agg_kind": mrow.get("agg_kind", ""), "coordinate": "FULL_FAMILY", "coordinate_role": "full_family_simultaneous", "primary_target_id": mrow["primary_target_id"], "primary_truth": "ALL_FROZEN_COORDINATES",
                "target_match": mrow["target_match"], "causal_validity": mrow["causal_validity"], "inference_eligible": "TRUE" if inf_elig else "FALSE", "moment_regime": m_regime, "moment_note": m_note,
                "n_scheduled": n_scheduled_reps, "n_pending": n_scheduled_reps-n_completed, "complete_100_draws": n_completed == 100, "n_attempts_including_interrupted": len(all_att_recs), "n_completed_attempts": n_completed, "n_available_points": n_full_elig, "n_available_intervals": n_full_elig if inf_elig else 0, "n_success": n_success, "n_partial_support": n_partial,
                "n_refusals": n_refusals, "n_unexpected_errors": n_unexp_err, "n_unexpected_success": n_unexp_succ, "n_timeouts": n_timeouts,
                "n_parent_unavailable": n_parent_unavail, "n_unavailable_comparator": n_unavail_comp, "n_interrupted": n_interrupted,
                "n_warning_draws": n_warn_draws, "n_support_changed": n_support_changed if call_type == "AGG" else 0, "mean_estimate": "",
                "bias": "", "mcse_bias": "", "empirical_sd": "", "rmse": "", "mean_se": "", "se_sd_ratio": "", "mean_se_analytical": "", "se_sd_ratio_analytical": "", "mean_width": "", "conditional_coverage": fmt_val(ff_cond),
                "cond_cov_ci_low": fmt_val(ff_cc_lo), "cond_cov_ci_high": fmt_val(ff_cc_hi), "all_attempt_containment": fmt_val(ff_all), "all_att_ci_low": fmt_val(ff_aa_lo), "all_att_ci_high": fmt_val(ff_aa_hi),
                "secondary_reference_id": mrow.get("secondary_reference_id", ""), "secondary_status": mrow.get("secondary_status", ""),
                "secondary_reference_value": "", "secondary_bias": "", "secondary_rmse": "", "secondary_conditional_coverage": "", **rt_fields, })
        # Auxiliary pre-period contrasts (reported separately as law-implied auxiliary contrasts, not primary causal ATT)
        for aux_k, aux_val in sorted(aux_truth_map.items()):
            target_list.append(build_coord_summary_row(aux_k, "pre_auxiliary_contrast", float(aux_val), None))
    # 3. Build paired.csv: paired method differences on common completed draws
    paired_csv_rows = []; manifest_by_tpl = {m["task_id"]: m for m in manifest_rows}
    def add_paired_comparison(pair_type: str, lhs_tpl: str, rhs_tpl: str, context_cell_id: str):
        lhs_meta = tasks_meta.get(lhs_tpl, {}); p_truth_map = lhs_meta.get("primary_truth_map", {}); p_coords = [k for k in manifest_by_tpl[lhs_tpl].get("primary_coordinate_set", "").split(";") if k in p_truth_map]
        m_regime = manifest_by_tpl[lhs_tpl].get("moment_regime", "finite_fourth_moment"); inf_var = (m_regime == "infinite_variance_t1_5")
        rhs_truth = tasks_meta.get(rhs_tpl, {}).get("primary_truth_map", {})
        if manifest_by_tpl[lhs_tpl]["recipe_id"] != manifest_by_tpl[rhs_tpl]["recipe_id"]:
            # Explicit aliases share one preserved stream and law by contract.
            if manifest_by_tpl[lhs_tpl]["stream_key"] != manifest_by_tpl[rhs_tpl]["stream_key"]: raise ValueError("Unpaired RNG streams")
        for ck in p_coords:
            if ck not in rhs_truth or not math.isclose(p_truth_map[ck], rhs_truth[ck], abs_tol=1e-12, rel_tol=0):
                raise ValueError(f"Paired primary truth mismatch: {lhs_tpl}, {rhs_tpl}, {ck}")
            t_val = float(p_truth_map[ck]); diffs = []; se_diffs = []; abs_err_diffs = []; lhs_errs = []; rhs_errs = []
            for rep in range(1, n_scheduled_reps + 1):
                r_l = terminal_by_tpl_rep.get((lhs_tpl, rep)); r_r = terminal_by_tpl_rep.get((rhs_tpl, rep))
                if not r_l or not r_r: continue
                map_l = extract_fit_coord_map(r_l); map_r = extract_fit_coord_map(r_r); cl = map_l.get(ck); cr = map_r.get(ck)
                if not cl or not cr or cl["est"] is None or cr["est"] is None: continue
                el, er = cl["est"], cr["est"]; diffs.append(el - er)
                if cl["se"] is not None and cr["se"] is not None:
                    se_diffs.append(cl["se"] - cr["se"])
                err_l = el - t_val; err_r = er - t_val; lhs_errs.append(err_l); rhs_errs.append(err_r); abs_err_diffs.append(abs(err_l) - abs(err_r))
            n_comm = len(diffs); sd_d = sd_or_none(diffs); mcse_d = (sd_d / math.sqrt(n_comm)) if (sd_d is not None and n_comm >= 2 and not inf_var) else None; rmse_l = rmse_or_none(lhs_errs); rmse_r = rmse_or_none(rhs_errs)
            paired_csv_rows.append({"pair_type": pair_type, "context_cell_id": context_cell_id, "recipe_id": manifest_by_tpl[lhs_tpl]["recipe_id"], "lhs_task_id": lhs_tpl, "rhs_task_id": rhs_tpl,
                "stream_key": manifest_by_tpl[lhs_tpl]["stream_key"], "pairing_key": "stream_key;recipe/law;replicate;coordinate;primary_truth", "moment_note": moment_qualification_note(m_regime),
                "coordinate": ck, "primary_target_id": manifest_by_tpl[lhs_tpl]["primary_target_id"], "primary_truth": fmt_val(t_val), "n_scheduled": n_scheduled_reps, "n_common_completed": n_comm,
                "mean_diff_est": fmt_val(mean_or_none(diffs)), "sd_diff_est": fmt_val(sd_d), "mcse_diff_est": fmt_val(mcse_d), "max_abs_diff_est": fmt_val(max(abs(x) for x in diffs) if diffs else None),
                "max_abs_diff_se": fmt_val(max(abs(x) for x in se_diffs) if se_diffs else None), "mean_abs_err_diff": fmt_val(mean_or_none(abs_err_diffs)), "lhs_rmse": fmt_val(rmse_l), "rhs_rmse": fmt_val(rmse_r),
                "rmse_diff": fmt_val((rmse_l - rmse_r) if (rmse_l is not None and rmse_r is not None) else None), })
    # Link comparator tasks back to each logical cell in linked_cell_ids
    for mrow in manifest_rows:
        if mrow["call_type"] in ("DID_CC", "DID_UNBALANCED", "DID_RCS"):
            comp_tpl = mrow["task_id"]
            for linked_cid in mrow.get("linked_cell_ids", "").split(";"):
                linked_cid = linked_cid.strip(); fit_tpl = f"{linked_cid}__FIT"
                if fit_tpl in manifest_by_tpl:
                    add_paired_comparison(f"CDID_vs_{mrow['call_type']}", fit_tpl, comp_tpl, linked_cid)
        # Also check equality_ref on FIT rows
        if mrow["call_type"] == "FIT" and mrow.get("equality_ref", "").strip():
            eq_raw = mrow["equality_ref"].strip(); cid = mrow["cell_id"]; lhs_tpl = mrow["task_id"]
            if eq_raw.startswith("did explicit-RCS"):
                # Match the recipe's DID_RCS_01 comparator
                rhs_tpl = f"{mrow['recipe_id']}__DID_RCS__01"
            else:
                ref_cid = eq_raw.split("(")[0].strip(); rhs_tpl = f"{ref_cid}__FIT"
            if rhs_tpl in manifest_by_tpl:
                add_paired_comparison(f"EQUALITY_REF ({eq_raw})", lhs_tpl, rhs_tpl, cid)
    # Write CSV files
    def write_csv(path: Path, rows: list, default_fields: list):
        fieldnames = list(dict.fromkeys(k for row in rows for k in row)) if rows else default_fields
        with path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames); w.writeheader()
            if rows:
                w.writerows(rows)
    write_csv(out_dir / "attempts.csv", attempts_csv_rows, ["attempt_id", "task_id", "disposition"]); write_csv(out_dir / "cells.csv", cells_csv_rows, ["task_id", "cell_id", "coordinate", "bias", "rmse"])
    write_csv(out_dir / "aggregates.csv", aggregates_csv_rows, ["task_id", "cell_id", "agg_kind", "coordinate", "bias", "rmse"]); write_csv(out_dir / "pretests.csv", pretests_csv_rows, ["task_id", "cell_id", "rejection_rate_005"])
    write_csv(out_dir / "paired.csv", paired_csv_rows, ["pair_type", "context_cell_id", "lhs_task_id", "rhs_task_id", "mean_diff_est"])
    # 4. Build human-readable SUMMARY.md
    disp_counts = {}
    for r in raw_attempts:
        d = r.get("disposition", "unknown"); disp_counts[d] = disp_counts.get(d, 0) + 1
    succ_ops_all = [float(r["operation_seconds"]) for r in raw_attempts if r.get("disposition") == "success" and r.get("operation_seconds") is not None]
    all_ops_all = [float(r["operation_seconds"]) for r in raw_attempts if r.get("operation_seconds") is not None]
    setup_all = [v for r in raw_attempts if (v := measured_setup_cost(r)) is not None]
    n_timeouts_all = disp_counts.get("timeout", 0)
    totals_all = [float(r["total_attempt_seconds"]) for r in raw_attempts if r.get("total_attempt_seconds") is not None]
    md_lines = ["# Portable Reliability Campaign Summary", "", f"- **Generated UTC**: `{datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}`",
        f"- **Manifest Tasks per Replicate**: `{len(manifest_rows)}` (`{len(manifest_rows) * n_scheduled_reps}` across `{n_scheduled_reps}` scheduled replicates)",
        f"- **Total Recorded Attempts on Disk**: `{len(raw_attempts)}`", f"- **Distinct (Task, Replicate) Slots Settled**: `{len(terminal_by_tpl_rep)}`",
        f"- **Pending slots**: `{len(manifest_rows)*n_scheduled_reps-len(terminal_by_tpl_rep)}`. Incomplete runs retain the scheduled denominator; containment is not a completed 100-draw result.",
        "- Primary causal ATT and auxiliary precontrasts occupy separate labelled rows. Unknown secondary limits remain blank and never suppress primary comparisons.",
        "- t3: finite variance, infinite fourth moment; bias MCSE only. t1.5: infinite variance; SD/RMSE are finite-sample descriptions and moment MCSE is undefined. Binomial intervals remain defined.",
        "- Pairing uses preserved stream, law, replicate, coordinate and matching primary truth; common-completed counts are explicit. Interrupted costs without a measured receipt remain unknown and are never treated as zero.",
        "- **Scope**: Descriptive 100-draw reliability campaign against common prespecified full-study causal ATT truths; no nominal PASS/qualification gate.",
        "", "## 1. Attempt Dispositions", "", "| Disposition | Attempt Count |", "| :--- | ---: |"]
    for d_key in ("success", "partial-support", "expected-refusal", "unexpected-success", "unexpected-error", "timeout", "interrupted", "parent-unavailable", "unavailable-comparator", "not-applicable", ):
        md_lines.append(f"| `{d_key}` | {disp_counts.get(d_key, 0)} |")
    md_lines.extend(["", "## 2. Runtime Summary (Separated by Role)", "", "| Metric | Value (seconds) | Count |", "| :--- | ---: | ---: |",
        f"| Successful operation latency (mean) | {fmt_val(mean_or_none(succ_ops_all)) or 'NA'} | {len(succ_ops_all)} |",
        f"| Successful operation latency (median) | {fmt_val(median_or_none(succ_ops_all)) or 'NA'} | {len(succ_ops_all)} |", f"| All-attempt operation time (mean) | {fmt_val(mean_or_none(all_ops_all)) or 'NA'} | {len(all_ops_all)} |",
        f"| All-attempt operation time (median) | {fmt_val(median_or_none(all_ops_all)) or 'NA'} | {len(all_ops_all)} |",
        f"| Measured all-attempt total cost (mean) | {fmt_val(mean_or_none(totals_all)) or 'NA'} | {len(totals_all)} |",
        f"| Measured all-attempt total cost (median) | {fmt_val(median_or_none(totals_all)) or 'NA'} | {len(totals_all)} |",
        f"| Unknown interrupted/crash costs | NA | {sum(r.get('total_attempt_seconds') is None for r in raw_attempts)} |",
        f"| Supervisor + R startup + data load/gen (mean) | {fmt_val(mean_or_none(setup_all)) or 'NA'} | {len(setup_all)} |",
        f"| Supervisor + R startup + data load/gen (median) | {fmt_val(median_or_none(setup_all)) or 'NA'} | {len(setup_all)} |",
        f"| Right-censored timeouts (at 300s cap) | 300.000000 | {n_timeouts_all} |", "", "## 3. Output Artifacts", "",
        "- `attempts.csv`: Immutable per-attempt ledger including interrupted attempts, timings, warnings, and condition classes.",
        "- `cells.csv`: Per-coordinate and `FULL_FAMILY` primary causal ATT accuracy, SE/SD, conditional coverage, all-attempt containment, exact 95% Clopper-Pearson intervals, `V_analytical` clustered SE diagnostics, and separate secondary comparator reference metrics.",
        "- `aggregates.csv`: Summary of all `aggte` tasks (`AGG_S`, `AGG_6` 6 variants, `AGG_NA` 2 variants) against frozen population cohort masses.",
        "- `pretests.csv`: Conditional pretrend test (`conditional_did_pretest`) rejection rates and Clopper-Pearson intervals under null and conditional alternatives.",
        "- `paired.csv`: Paired differences on common completed draws (`cdid` vs `DID_CC`, `DID_UNBALANCED`, `DID_RCS`, plus `equality_ref` verification)." ])
    summary_md_path = out_dir / "SUMMARY.md"; summary_md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    return {"status": "SUMMARIZE_OK", "n_attempts": len(attempts_csv_rows), "n_cells_rows": len(cells_csv_rows), "n_aggregates_rows": len(aggregates_csv_rows), "n_pretests_rows": len(pretests_csv_rows),
        "n_paired_rows": len(paired_csv_rows), "summary_md": str(summary_md_path)}
def run_selftest():
    """Verify summarize.py on synthetic disk records including t3/t1.5 moment regimes, V_analytical, and paired diffs."""
    # Verify exact Clopper-Pearson interval against known stats::binom.test values
    for k, n, expected in ((0,1,(0,0.975)), (1,1,(0.025,1)),
                            (1,2,(1-math.sqrt(0.975), math.sqrt(0.975))),
                            (5,10,(0.1870860284473985,0.8129139715526015)),
                            (100,100,(0.025**0.01,1))):
        got = cp_interval(k,n)
        assert all(abs(a-b)<1e-12 for a,b in zip(got,expected)), (k,n,got)
    lo_95, hi_95 = cp_interval(95, 100, 0.05); assert abs(lo_95 - 0.8871651) < 1e-6 and abs(hi_95 - 0.9835681) < 1e-6; lo_0, hi_0 = cp_interval(0, 100, 0.05); assert lo_0 == 0.0 and abs(hi_0 - 0.03621669) < 1e-6
    with tempfile.TemporaryDirectory(prefix="cdid_summarize_selftest_") as tmp:
        tmp_path = Path(tmp); manifest_path = tmp_path / "manifest.csv"; targets_path = tmp_path / "targets.json"; out_dir = tmp_path / "out"; records_dir = out_dir / "records"; records_dir.mkdir(parents=True, exist_ok=True)
        m_rows = [{"task_id": "I_MOMENT_t1.5__FIT", "replicate_scope": "1..2", "stream_index": 46, "stream_key": "MOMENT_t1.5", "recipe_id": "MOMENT_t1.5", "family": "ERROR_MOMENTS",
            "cell_id": "I_MOMENT_t1.5", "logical_source": "targeted", "role": "moment-failure", "call_type": "FIT", "package": "cdid", "parent_task_id": "", "agg_kind": "", "literal_call": 'cdid::att_gt(yname="y")', "rng_purpose": 3,
            "primary_target_id": "ORIGINAL_COHORT_ATT", "primary_truth": '{"gt:3:3":1.2}', "primary_population_masses": '{"0":0.3, "3":0.25, "4":0.25, "5":0.2}',
            "primary_coordinate_set": "gt:3:3", "target_match": "matched", "causal_validity": "VALID", "inference_eligible": "TRUE", "expected_disposition": "success", "expected_condition": "",
            "secondary_reference_id": "RCS_COMP_ARBITRARY", "secondary_reference": '{"gt:3:3":1.05}',
            "secondary_status": "known_biased_limit", "moment_regime": "infinite_variance_t1_5", "equality_ref": "", "linked_cell_ids": "I_MOMENT_t1.5"}, {
            "task_id": "MOMENT_t1.5__DID_UNBALANCED__01", "replicate_scope": "1..2", "stream_index": 46, "stream_key": "MOMENT_t1.5", "recipe_id": "MOMENT_t1.5", "family": "ERROR_MOMENTS",
            "cell_id": "I_MOMENT_t1.5", "logical_source": "comparator", "role": "comparator",
            "call_type": "DID_UNBALANCED", "package": "did", "parent_task_id": "", "agg_kind": "", "literal_call": 'did::att_gt(yname="y")', "rng_purpose": 3,
            "primary_target_id": "ORIGINAL_COHORT_ATT", "primary_truth": '{"gt:3:3":1.2}', "primary_population_masses": '{"0":0.3, "3":0.25, "4":0.25, "5":0.2}',
            "primary_coordinate_set": "gt:3:3", "target_match": "mismatched", "causal_validity": "NO_GUARANTEE", "inference_eligible": "TRUE", "expected_disposition": "success", "expected_condition": "",
            "secondary_reference_id": "RCS_COMP_ARBITRARY", "secondary_reference": '{"gt:3:3":1.05}',
            "secondary_status": "known_biased_limit", "moment_regime": "infinite_variance_t1_5", "equality_ref": "", "linked_cell_ids": "I_MOMENT_t1.5"},]
        with manifest_path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(m_rows[0].keys())); w.writeheader(); w.writerows(m_rows)
        targets_tasks = {r["task_id"]: {**r, "primary_truth_map": {"gt:3:3": 1.2}, "auxiliary_truth_map": {"gt:3:2": 0.0},
            "secondary_reference_map": {"gt:3:3": 1.05}, "population_masses_map": {"0": 0.3, "3": 0.25, "4": 0.25, "5": 0.2}} for r in m_rows}
        targets_path.write_text(json.dumps({"tasks": targets_tasks}), encoding="utf-8")
        for rep, est_cdid, est_did in ((1, 1.22, 1.06), (2, 1.18, 1.04)):
            rec_c = {"attempt_id": f"r{rep:03d}__I_MOMENT_t1.5__FIT__att01", "task_id": f"r{rep:03d}__I_MOMENT_t1.5__FIT", "template_task_id": "I_MOMENT_t1.5__FIT",
                "attempt_number": 1, "replicate": rep, "stream_index": 46, "recipe_id": "MOMENT_t1.5", "family": "ERROR_MOMENTS", "cell_id": "I_MOMENT_t1.5", "call_type": "FIT", "package": "cdid",
                "disposition": "success", "worker_launched": True, "censored_timeout": False, "supervisor_overhead_seconds": 0.01, "startup_seconds": 0.02, "data_seconds": 0.01, "operation_seconds": 0.10, "total_attempt_seconds": 0.14,
                "fit_result": { "cells": { "group": [3.0, 3.0], "t": [2.0, 3.0], "att": [0.01, est_cdid], "se": [0.05, 0.05], "se_analytical": [0.048, 0.049], "ci_lower": [-0.09, est_cdid - 0.10],
                "ci_upper": [0.11, est_cdid + 0.10], }, "cband": False, }}
            (records_dir / f"{rec_c['attempt_id']}.result.json").write_text(json.dumps(rec_c), encoding="utf-8")
            rec_d = {**rec_c, "attempt_id": f"r{rep:03d}__MOMENT_t1.5__DID_UNBALANCED__01__att01", "task_id": f"r{rep:03d}__MOMENT_t1.5__DID_UNBALANCED__01",
                "template_task_id": "MOMENT_t1.5__DID_UNBALANCED__01", "call_type": "DID_UNBALANCED", "package": "did", "fit_result": { "cells": { "group": [3.0, 3.0], "t": [2.0, 3.0], "att": [-0.02,
                est_did], "se": [0.05, 0.05], "ci_lower": [-0.12, est_did - 0.10], "ci_upper": [0.08, est_did + 0.10], }, "cband": False, }}
            (records_dir / f"{rec_d['attempt_id']}.result.json").write_text(json.dumps(rec_d), encoding="utf-8")
        # An interrupted predecessor has no measured setup receipt and must not lower the mean.
        first_path = records_dir / "r001__I_MOMENT_t1.5__FIT__att01.result.json"
        completed = json.loads(first_path.read_text()); completed["attempt_number"] = 2
        completed["attempt_id"] = "r001__I_MOMENT_t1.5__FIT__att02"
        interrupted = {**completed, "attempt_number": 1, "attempt_id": "r001__I_MOMENT_t1.5__FIT__att01",
            "disposition": "interrupted", "total_attempt_seconds": None, "operation_seconds": None,
            "startup_seconds": None, "data_seconds": None, "supervisor_overhead_seconds": None}
        first_path.write_text(json.dumps(interrupted))
        (records_dir / "r001__I_MOMENT_t1.5__FIT__att02.result.json").write_text(json.dumps(completed))
        assert measured_setup_cost(interrupted) is None
        assert measured_setup_cost({"supervisor_overhead_seconds": 0, "startup_seconds": 0, "data_seconds": 0}) == 0
        res = summarize_campaign(manifest_path, targets_path, out_dir, n_scheduled_reps=2); cells_out = list(csv.DictReader((out_dir / "cells.csv").open(encoding="utf-8")))
        t15_row = [r for r in cells_out if r["task_id"] == "I_MOMENT_t1.5__FIT" and r["coordinate"] == "gt:3:3"][0]; assert t15_row["mcse_bias"] == "", "MOMENT_t1.5 must leave moment-based MCSE undefined"
        assert float(t15_row["mean_setup_and_data_seconds"]) == .04
        assert int(t15_row["n_unknown_attempt_costs"]) == 1
        assert float(t15_row["conditional_coverage"]) == 1.0; assert float(t15_row["mean_se_analytical"]) > 0.0; paired_out = list(csv.DictReader((out_dir / "paired.csv").open(encoding="utf-8"))); assert len(paired_out) == 1
        assert int(paired_out[0]["n_common_completed"]) == 2; assert res["status"] == "SUMMARIZE_OK"
    return {"status": "SUMMARIZE_SELFTEST_OK"}
def main():
    parser = argparse.ArgumentParser(description="Generate tidy CSV and Markdown summaries from immutable attempt records."); parser.add_argument("--manifest", default="", help="Path to manifest.csv.")
    parser.add_argument("--targets", default="", help="Path to targets.json (defaults to alongside manifest.csv)."); parser.add_argument("--out", default="", help="Path to campaign results directory.")
    parser.add_argument("--reps", type=int, choices=(100,), default=100, help="Scheduled denominator is fixed at 100; incomplete runs retain it."); parser.add_argument("--selftest", action="store_true", help="Run summarize.py self-test suite."); args = parser.parse_args()
    if args.selftest:
        print(json.dumps(run_selftest(), indent=2)); return
    if not args.out: raise SystemExit("Usage: summarize.py --manifest manifest.csv --out RESULTS")
    out_dir = Path(args.out).resolve(); manifest_path = Path(args.manifest).resolve() if args.manifest else (out_dir / "manifest.csv")
    if not manifest_path.exists():
        manifest_path = Path(__file__).resolve().parent / "manifest.csv"
    targets_path = Path(args.targets).resolve() if args.targets else (out_dir / "targets.json")
    if not targets_path.exists():
        targets_path = Path(__file__).resolve().parent / "targets.json"
    res = summarize_campaign(manifest_path, targets_path, out_dir, n_scheduled_reps=args.reps); print(json.dumps(res, indent=2))
if __name__ == "__main__":
    main()
