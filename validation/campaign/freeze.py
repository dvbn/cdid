#!/usr/bin/env python3
"""Portable reliability campaign: manifest builder, exact truth calculator, and freeze helper.
Uses only the Python 3 standard library. Does not execute DGP draws, estimator fits,
or R RNG during manifest generation or pre-freeze verification (--check).
"""
import argparse, csv, datetime, fractions, hashlib, itertools, json, math, os, platform, re, subprocess, sys; from pathlib import Path; SEED_MASTER = 61012026; N_STREAMS = 59; N_REPS = 100; N_PURPOSES = 5
# Fixed public formal names; verified against bundled data and installed functions at freeze.
CDID_ATT_GT_FORMALS = ("yname", "tname", "idname", "gname", "xformla", "data", "panel", "allow_unbalanced_panel", "control_group",
    "anticipation", "weightsname", "fix_weights", "alp", "bstrap", "cband", "biters", "clustervars", "est_method", "base_period", "faster_mode", "print_details", "pl", "cores", "compute_inffunc", "mode", "frame", "selection",
    "selection_model", "selection_pooling", "selection_xformla", "selection_universe", "selection_rho",
    "selection_justification", "selection_est_method", "imputation_pooling", "gmm_weighting", "gmm_covariance", "pt_assumption", "pt_start", "control")
DID_ATT_GT_FORMALS = ("yname", "tname", "idname", "gname", "xformla", "data", "panel", "allow_unbalanced_panel", "control_group",
    "anticipation", "weightsname", "fix_weights", "alp", "bstrap", "cband", "biters", "clustervars", "est_method", "base_period", "faster_mode", "print_details", "pl", "cores", "compute_inffunc")
CDID_AGGTE_FORMALS = ("MP", "type", "balance_e", "min_e", "max_e", "na.rm", "bstrap", "biters", "cband", "alp", "clustervars")
CDID_PRETEST_FORMALS = ("yname", "tname", "idname", "gname", "xformla", "data", "panel", "allow_unbalanced_panel", "control_group", "weightsname", "alp", "bstrap", "cband", "biters", "clustervars", "est_method", "print_details", "pl",
    "cores", "anticipation", "base_period", "fix_weights", "mode", "frame", "selection", "selection_model", "selection_pooling", "selection_xformla", "selection_universe", "selection_rho", "selection_justification",
    "selection_est_method", "imputation_pooling", "control")
CHAINED_ONLY_ARGS = ("frame", "selection", "selection_model", "selection_pooling", "selection_xformla", "selection_universe", "selection_rho", "selection_justification", "selection_est_method", "imputation_pooling", "gmm_weighting",
    "gmm_covariance", "pt_assumption", "pt_start", "control")
PRETEST_EXCLUDED_ARGS = ("pt_assumption", "pt_start", "gmm_weighting", "gmm_covariance", "compute_inffunc", "allow_unbalanced_panel")
def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()
def expit(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))
def split_top_level_commas(s: str):
    parts = []; cur = []; depth = 0; in_quote = False
    for ch in s:
        if ch == '"':
            in_quote = not in_quote; cur.append(ch)
        elif not in_quote and ch in "([":
            depth += 1; cur.append(ch)
        elif not in_quote and ch in ")]":
            depth -= 1; cur.append(ch)
        elif not in_quote and depth == 0 and ch == ",":
            tok = "".join(cur).strip()
            if tok:
                parts.append(tok)
            cur = []
        else:
            cur.append(ch)
    tok = "".join(cur).strip()
    if tok:
        parts.append(tok)
    return parts
def parse_r_args_string(arg_str: str):
    """Parse a comma-separated R argument string from RECIPES.csv or CELLS.csv."""
    flags = set(); kv = {}
    if not arg_str or not arg_str.strip(): return flags, kv
    for part in split_top_level_commas(arg_str.strip()):
        if "=" not in part:
            flags.add(part.strip())
        else:
            k, v = part.split("=", 1); kv[k.strip()] = v.strip()
    return flags, kv
def default_fit_args(xformla_default: str):
    """Base FIT argument dictionary from CAMPAIGN_DESIGN.md Section 1."""
    return {"yname": '"y"', "tname": '"t"', "idname": '"id"', "gname": '"G"', "xformla": xformla_default.strip(), "data": "PANEL", "panel": "TRUE", "allow_unbalanced_panel": "TRUE", "control_group": '"nevertreated"',
        "anticipation": "0", "alp": "0.05", "bstrap": "FALSE", "cband": "FALSE", "est_method": '"dr"', "base_period": '"varying"', "compute_inffunc": "TRUE", "pl": "FALSE", "cores": "1", "mode": '"chained"',
        "frame": "F", "selection": '"none"', "gmm_weighting": '"optimal"', "pt_assumption": '"post"'}
UPSTREAM_DEFAULTS = {"anticipation": "0", "weightsname": "NULL", "fix_weights": "NULL", "alp": "0.05",
    "bstrap": "TRUE", "cband": "TRUE", "biters": "1000", "clustervars": "NULL", "est_method": '"dr"',
    "base_period": '"varying"', "faster_mode": "TRUE", "print_details": "FALSE", "pl": "FALSE", "cores": "1", "compute_inffunc": "TRUE"}
def effective_upstream_args(fit_args, rcs=False):
    args = {**UPSTREAM_DEFAULTS, **{k:v for k,v in fit_args.items() if k in DID_ATT_GT_FORMALS}}
    if rcs:
        args.pop("idname", None)
        args.update(panel="FALSE", allow_unbalanced_panel="FALSE", bstrap="FALSE", cband="FALSE", compute_inffunc="FALSE")
        args["clustervars"] = "NULL"
    else:
        args.update(panel="TRUE", allow_unbalanced_panel="TRUE")
    if args["compute_inffunc"] == "FALSE":
        args.update(bstrap="FALSE", cband="FALSE")
    return args
def upstream_key(fit_args, rcs=False):
    return tuple(sorted(effective_upstream_args(fit_args, rcs).items()))
def build_effective_fit_args(recipe_row: dict, cell_arg_str: str):
    args = default_fit_args(recipe_row["xformla_default"]); _, sweep_kv = parse_r_args_string(recipe_row.get("sweep_fixed_args", ""))
    for k, v in sweep_kv.items():
        args[k] = v
    flags, cell_kv = parse_r_args_string(cell_arg_str)
    if "UPSTREAM_ONLY" in flags:
        for drop_k in ("allow_unbalanced_panel",) + CHAINED_ONLY_ARGS:
            args.pop(drop_k, None)
    for k, v in cell_kv.items():
        args[k] = v
    return flags, args
def format_r_call(fn_name: str, arg_dict: dict, formal_order: tuple) -> str:
    ordered_keys = [k for k in formal_order if k in arg_dict]; extra_keys = [k for k in arg_dict if k not in formal_order]
    if extra_keys: raise ValueError(f"Unknown formal argument(s) for {fn_name}: {extra_keys}")
    pieces = [f"{k}={arg_dict[k]}" for k in ordered_keys]; return f"{fn_name}({', '.join(pieces)})"
def build_pretest_args(fit_args: dict):
    pre_args = {}
    for k, v in fit_args.items():
        if k in PRETEST_EXCLUDED_ARGS: continue
        if k in CDID_PRETEST_FORMALS:
            pre_args[k] = v
    pre_args["biters"] = "999"; return pre_args
def expand_children_specs(children_code: str):
    code = (children_code or "").strip()
    if not code or code == "none": return []
    if code == "AGG_S": return [("AGG_S", "simple", {"MP": "FIT", "type": '"simple"'})]
    if code == "AGG_6":
        return [("AGG_6_simple", "simple", {"MP": "FIT", "type": '"simple"'}), ("AGG_6_group", "group", {"MP": "FIT", "type": '"group"'}), ("AGG_6_dynamic", "dynamic", {"MP": "FIT", "type": '"dynamic"'}),
            ("AGG_6_dynamic_bal1", "dynamic_bal1", {"MP": "FIT", "type": '"dynamic"', "balance_e": "1"}), ("AGG_6_dynamic_0_2", "dynamic_0_2", {"MP": "FIT", "type": '"dynamic"', "min_e": "0", "max_e": "2"}),
            ("AGG_6_calendar", "calendar", {"MP": "FIT", "type": '"calendar"'})]
    if code == "AGG_NA": return [("AGG_NA_false", "simple_narm_false", {"MP": "FIT", "type": '"simple"', "na.rm": "FALSE"}), ("AGG_NA_true", "simple_narm_true", {"MP": "FIT", "type": '"simple"', "na.rm": "TRUE"}),]
    raise ValueError(f"Unknown children specification: {code}")
def parse_calendar(recipe_row: dict):
    T = int(recipe_row["n_dates"]); cal_str = recipe_row["calendar"].strip()
    if cal_str == "1:T": return [float(k) for k in range(1, T + 1)]
    m = re.match(r"^c\(([\d,\s]+)\)$", cal_str)
    if m:
        vals = [float(x.strip()) for x in m.group(1).split(",") if x.strip()]; assert len(vals) == T; return vals
    raise ValueError(f"Unsupported calendar string: {cal_str}")
def continuous_cohort_mass(quadratic=False):
    # Composite Simpson integration of a smooth fixed density (no draws or estimates).
    n = 4096
    return sum((1 if i in (0, n) else 4 if i % 2 else 2) *
               expit(0.2 + 0.6 * ((-1 + 2*i/n)**2 if quadratic else (-1 + 2*i/n)))
               for i in range(n + 1)) / (3*n)
def recipe_population_info(recipe_id: str, recipe_row: dict, weightsname: str, fix_weights: str):
    cal = parse_calendar(recipe_row); weighted = bool(weightsname and weightsname != "NULL")
    default_0345_recipes = {"COMP_adjacent", "COMP_arbitrary", "CC_BASE", "CC_COMPOSITION", "BAL_BASE", "C_MCAR_DATE", "C_MCAR_PAIR", "COMP_LADDER_01", "COMP_LADDER_02", "COMP_LADDER_03", "COMP_LADDER_04", "COMP_LADDER_05",
        "COMP_LADDER_06", "COMP_LADDER_07", "COMP_LADDER_08", "CC_LADDER_01", "CC_LADDER_02", "CC_LADDER_03", "CC_LADDER_04", "DATE_BASE", "DATE_HET", "DATE_COHORT_SLOPES", "REPRESENTATIVE_ROSTERS",
        "DATE_CONSTANT_EFFECT", "DATE_INTERACTION_Q", "DATE_DISTINCT_0", "DATE_DISTINCT_1", "DATE_WRONG_Q_J", "INNOVATION_SELECTION", "PT_POST_INVALID", "ANTICIPATION_ONE", "DATED_WEIGHTS", "W_STATIC",
        "CLUSTER_N2000_C100", "CLUSTER_N8000_C400", "CLUSTER_N2000_C20", "CLUSTER_N2000_C5", "MOMENT_t3", "MOMENT_t1.5", "GMM_SPARSE", "GMM_ZERO_SCALE", "BAL_ALIAS", "F_NO_INTERCEPT", "DATE_NO_INTERCEPT",
        "MCAR_HALF", "E_OR_NEAR_COLLINEAR", "DATED_MISSING", "E_CHANGING_WEIGHT_NULL", "E_SELECTION_SEPARATED", "SUPPORT_DISCONNECTED", "AR1_MCAR_DATE", "ATTRITION_MTAR"}
    if recipe_id in default_0345_recipes or recipe_id in ("COMP_LATE", "GAP_LABELS"):
        if recipe_id == "COMP_LATE":
            codes_pos = [0, 3, 7, 11]; labels = [0.0, 3.0, 7.0, 11.0]
        elif recipe_id == "GAP_LABELS":
            codes_pos = [0, 3, 4, 5]; labels = [0.0, cal[2], cal[3], cal[4]]
        else:
            codes_pos = [0, 3, 4, 5]; labels = [0.0, 3.0, 4.0, 5.0]
        pi_base = [0.30, 0.25, 0.25, 0.20]; mu_base = [1.0 / 3.0, 0.4, 0.6, 0.75]; q0 = [0.4, 0.3, 0.2, 0.1]; q1 = [0.2, 0.2, 0.3, 0.3]
        if weighted and recipe_id == "W_STATIC":
            pi_raw = [p * (1.0 + 0.5 * m) for p, m in zip(pi_base, mu_base)]
        elif weighted and recipe_id == "DATED_WEIGHTS":
            pi_raw = [p * (1.0 + 0.1 * m) for p, m in zip(pi_base, mu_base)]
        else:
            pi_raw = list(pi_base)
        return {"calendar": cal, "codes_pos": codes_pos, "labels": labels, "pi_base": pi_base, "pi_raw": pi_raw, "mu_base": mu_base, "q0": q0, "q1": q1}
    if recipe_id in ("CONT_PS_only", "CONT_OR_only", "CONT_both_wrong", "CONT_both_correct", "TEST_null", "TEST_conditional_alternative", "E_PREDICTION_SUPPORT"):
        p3 = continuous_cohort_mass(recipe_id in ("CONT_OR_only", "CONT_both_wrong")) if recipe_id.startswith("CONT_") else 0.5
        return {"calendar": cal, "codes_pos": [0, 3], "labels": [0.0, 3.0], "pi_base": [1-p3, p3], "pi_raw": [1-p3, p3], "mu_base": [0.5, 0.5], "q0": [1-p3, p3], "q1": [1-p3, p3]}
    if recipe_id == "E_OVERLAP":
        return {"calendar": cal, "codes_pos": [0, 3], "labels": [0.0, 3.0], "pi_base": [0.5, 0.5], "pi_raw": [0.5, 0.5], "mu_base": [0.0, 1.0], "q0": [1.0, 0.0], "q1": [0.0, 1.0]}
    if recipe_id == "EARLY_DRIFT":
        return {"calendar": cal, "codes_pos": [0, 4], "labels": [0.0, 4.0], "pi_base": [0.5, 0.5], "pi_raw": [0.5, 0.5], "mu_base": [0.5, 0.5], "q0": [0.5, 0.5], "q1": [0.5, 0.5]}
    if recipe_id == "E_TWO_DATES_A":
        return {"calendar": cal, "codes_pos": [0, 2], "labels": [0.0, 2.0], "pi_base": [0.5, 0.5], "pi_raw": [0.5, 0.5], "mu_base": [0.5, 0.5], "q0": [0.5, 0.5], "q1": [0.5, 0.5]}
    raise ValueError(f"Unhandled recipe_id in recipe_population_info: {recipe_id}")
def fmt_num(v: float) -> str:
    if float(v).is_integer(): return str(int(v))
    return f"{v:.6g}"
def primary_target_identifier(recipe_id: str, weightsname: str, fix_weights: str, call_type: str = "FIT") -> str:
    if call_type == "PRETEST": return "PRETEST_CONDITIONAL_PT_NULL"
    weighted = bool(weightsname and weightsname != "NULL")
    if not weighted: return "ORIGINAL_COHORT_ATT"
    if recipe_id == "W_STATIC": return f"{recipe_id}__WEIGHTED_STATIC_ATT"
    if recipe_id == "DATED_WEIGHTS":
        if fix_weights == '"base_period"': return f"{recipe_id}__WEIGHTED_BASE_PERIOD_ATT"
        return f"{recipe_id}__WEIGHTED_FIRST_PERIOD_ATT"
    return "ORIGINAL_COHORT_ATT"
def compute_cell_primary_and_aux_truth(recipe_id: str, recipe_row: dict, base_period: str, anticipation: int, weightsname: str, fix_weights: str, ):
    info = recipe_population_info(recipe_id, recipe_row, weightsname, fix_weights); cal = info["calendar"]; T = len(cal); weighted = bool(weightsname and weightsname != "NULL"); base_period_clean = base_period.strip('"')
    # Map sorted union of date and positive cohort labels to consecutive positions q (1-based)
    pos_labels = [lbl for lbl in info["labels"] if lbl > 0]; u_labels = sorted(set(cal) | set(pos_labels))
    q_map = {lbl: idx + 1 for idx, lbl in enumerate(u_labels)}; cells = []; primary_coords = {}; aux_coords = {}; unweighted_diag_coords = {}
    for rank_idx in range(1, len(info["codes_pos"])):
        g_pos = info["codes_pos"][rank_idx]; g_lbl = info["labels"][rank_idx]; mu_g = info["mu_base"][rank_idx]; theta_g = 1.0 + 0.1 * (g_pos - 3.0)
        # b_g = max{d_j: d_j + a < g_lbl}
        eligible_b = [(j + 1, d_j) for j, d_j in enumerate(cal) if d_j + anticipation < g_lbl]; assert eligible_b, f"No base date before cohort {g_lbl} with anticipation={anticipation}"; b_g_pos, b_g_lbl = eligible_b[-1]
        for k_idx, t_lbl in enumerate(cal):
            t_pos = k_idx + 1
            if base_period_clean == "varying" and t_pos == 1:
                continue  # varying base omits first calendar date
            if base_period_clean == "universal":
                act_b_pos, act_b_lbl = b_g_pos, b_g_lbl
            else:
                if t_lbl < g_lbl:
                    act_b_pos, act_b_lbl = t_pos - 1, cal[t_pos - 2]
                else:
                    act_b_pos, act_b_lbl = b_g_pos, b_g_lbl
            coord_key = f"gt:{fmt_num(g_lbl)}:{fmt_num(t_lbl)}"; is_causal_post = (t_pos >= g_pos) or (recipe_id == "ANTICIPATION_ONE" and anticipation == 1 and t_pos == g_pos - 1)
            # Unweighted causal ATT at (g, t)
            if t_pos < g_pos:
                unw_att = 0.3 if (recipe_id == "ANTICIPATION_ONE" and t_pos == g_pos - 1) else 0.0
            elif recipe_id == "DATE_CONSTANT_EFFECT":
                unw_att = theta_g + 0.25 * (t_pos - g_pos)
            elif recipe_id == "GMM_ZERO_SCALE":
                unw_att = theta_g
            elif recipe_id == "DATE_NO_INTERCEPT":
                unw_att = 0.5 * mu_g
            elif recipe_id == "E_TWO_DATES_A":
                unw_att = 1.0 + 0.5 * mu_g
            else:
                unw_att = theta_g + 0.25 * (t_pos - g_pos) + 0.5 * mu_g
            # Weighted causal ATT at (g, t) when weightsname is active
            if weighted and recipe_id == "W_STATIC":
                mu_eff = 1.5 * mu_g / (1.0 + 0.5 * mu_g)
            elif weighted and recipe_id == "DATED_WEIGHTS":
                anchor = float(b_g_lbl) if fix_weights == '"base_period"' else 1.0; mu_eff = mu_g * (1.0 + 0.1 * anchor) / (1.0 + 0.1 * anchor * mu_g)
            else:
                mu_eff = mu_g
            if t_pos < g_pos:
                wt_att = 0.3 if (recipe_id == "ANTICIPATION_ONE" and anticipation == 1 and t_pos == g_pos - 1) else 0.0
            elif recipe_id == "DATE_CONSTANT_EFFECT":
                wt_att = theta_g + 0.25 * (t_pos - g_pos)
            elif recipe_id == "GMM_ZERO_SCALE":
                wt_att = theta_g
            elif recipe_id == "DATE_NO_INTERCEPT":
                wt_att = 0.5 * mu_eff
            elif recipe_id == "E_TWO_DATES_A":
                wt_att = 1.0 + 0.5 * mu_eff
            else:
                wt_att = theta_g + 0.25 * (t_pos - g_pos) + 0.5 * mu_eff
            if is_causal_post:
                role = "primary_post"; p_truth = wt_att; aux_truth = None; primary_coords[coord_key] = round(p_truth, 12)
                if weighted:
                    unweighted_diag_coords[coord_key] = round(unw_att, 12)
            elif base_period_clean == "universal" and t_pos == b_g_pos:
                role = "reference_zero"; p_truth = None; aux_truth = 0.0; aux_coords[coord_key] = 0.0
            else:
                role = "pre_auxiliary_contrast"; p_truth = None
                # Law-implied auxiliary pre-period contrast E[Y_t - Y_b | G=g] - E[Y_t - Y_b | C]
                if recipe_id == "EARLY_DRIFT" and g_pos == 4:
                    drift_t = 0.4 * max(3 - t_pos, 0); drift_b = 0.4 * max(3 - act_b_pos, 0); aux_truth = drift_t - drift_b
                elif recipe_id == "ANTICIPATION_ONE" and anticipation == 0:
                    eff_t = 0.3 if t_pos == (g_pos - 1) else 0.0; eff_b = 0.3 if act_b_pos == (g_pos - 1) else 0.0; aux_truth = eff_t - eff_b
                else:
                    aux_truth = 0.0
                aux_coords[coord_key] = round(aux_truth, 12)
            cells.append({"coord_key": coord_key, "group": g_lbl, "t": t_lbl, "g_pos": g_pos, "t_pos": t_pos, "b_pos": act_b_pos, "b_lbl": act_b_lbl, "q_g": q_map[g_lbl], "q_t": q_map[t_lbl], "role": role,
                "primary_truth": p_truth, "auxiliary_truth": aux_truth, "unweighted_truth": unw_att if is_causal_post else None, })
    masses_map = {fmt_num(lbl): round(w, 12) for lbl, w in zip(info["labels"], info["pi_raw"])}
    return {"info": info, "cells": cells, "primary_coords": primary_coords, "aux_coords": aux_coords, "unweighted_diag_coords": unweighted_diag_coords, "masses_map": masses_map}
def compute_aggte_truth_py(truth_bundle: dict, agg_kind: str):
    info = truth_bundle["info"]; cal = info["calendar"]; cells = truth_bundle["cells"]; w_map = {lbl: w for lbl, w in zip(info["labels"], info["pi_raw"]) if lbl > 0}
    post = [c for c in cells if c["q_t"] >= c["q_g"] and c["primary_truth"] is not None]
    if agg_kind in ("simple", "simple_narm_false", "simple_narm_true"):
        num = sum(w_map[c["group"]] * c["primary_truth"] for c in post); den = sum(w_map[c["group"]] for c in post); overall = round(num / den, 12)
        return {"overall": overall, "coordinates": {"overall": overall}, "coordinate_keys": ["overall"],}
    if agg_kind == "group":
        groups = sorted({c["group"] for c in post}); coords = {}; att_g = {}
        for g in groups:
            sub = [c["primary_truth"] for c in post if c["group"] == g]; val = sum(sub) / len(sub); att_g[g] = val; coords[f"egt:{fmt_num(g)}"] = round(val, 12)
        overall = round(sum(w_map[g] * att_g[g] for g in groups) / sum(w_map[g] for g in groups), 12); coords = {"overall": overall, **coords}; return {"overall": overall, "coordinates": coords, "coordinate_keys": list(coords.keys()),}
    if agg_kind in ("dynamic", "dynamic_bal1", "dynamic_0_2"):
        d1, dT = min(cal), max(cal); balance_e = 1.0 if agg_kind == "dynamic_bal1" else None; min_e = 0.0 if agg_kind == "dynamic_0_2" else -math.inf; max_e = 2.0 if agg_kind == "dynamic_0_2" else math.inf; elig = []
        for c in cells:
            val = c["primary_truth"] if c["primary_truth"] is not None else c["auxiliary_truth"]
            if val is None: continue
            e = c["t"] - c["group"]
            if balance_e is not None:
                if (dT - c["group"]) < balance_e: continue
                if e < (balance_e - dT + d1) or e > balance_e: continue
            if e < min_e or e > max_e: continue
            elig.append((e, c["group"], val))
        e_vals = sorted({item[0] for item in elig}); att_e = {}; coords = {}
        for ev in e_vals:
            sub = [(g, v) for (e_item, g, v) in elig if e_item == ev]; val = sum(w_map[g] * v for g, v in sub) / sum(w_map[g] for g, _ in sub); att_e[ev] = val; coords[f"egt:{fmt_num(ev)}"] = round(val, 12)
        nonneg = [att_e[ev] for ev in e_vals if ev >= 0]; overall = round(sum(nonneg) / len(nonneg), 12); coords = {"overall": overall, **coords}
        aux = {key: value for key, value in coords.items() if key.startswith("egt:") and float(key[4:]) < 0}
        primary = {key: value for key, value in coords.items() if key not in aux}
        return {"overall": overall, "coordinates": primary, "coordinate_keys": list(primary), "auxiliary_coordinates": aux}
    if agg_kind == "calendar":
        t_vals = sorted({c["t"] for c in post}); att_t = {}; coords = {}
        for tv in t_vals:
            sub = [c for c in post if c["t"] == tv]; val = sum(w_map[c["group"]] * c["primary_truth"] for c in sub) / sum(w_map[c["group"]] for c in sub); att_t[tv] = val; coords[f"egt:{fmt_num(tv)}"] = round(val, 12)
        overall = round(sum(att_t.values()) / len(att_t), 12); coords = {"overall": overall, **coords}; return {"overall": overall, "coordinates": coords, "coordinate_keys": list(coords.keys()),}
    raise ValueError(f"Unknown agg_kind: {agg_kind}")
def compute_rcs_comp_coords(truth_bundle: dict, control_group: str, anticipation: int, edges_kind: str):
    """Exact RCS_COMP probability limit from CAMPAIGN_DESIGN.md Section 4.3."""
    info = truth_bundle["info"]; cal = info["calendar"]; T = len(cal); q0 = info["q0"]; q1 = info["q1"]; lam = [0.0, 0.6, 0.4, 0.2]; cg_clean = control_group.strip('"'); H = {}
    for d in range(1, T + 1):
        if edges_kind == "adjacent":
            if d == 1:
                H[d] = 1.0
            elif d == T:
                H[d] = -1.0
            else:
                H[d] = (T + 1.0 - 2.0 * d) / (T - 2.0)
        else:
            H[d] = ((T + 1.0) / 2.0 - d) / (T - 1.0)
    out = {}
    for c in truth_bundle["cells"]:
        if c["role"] != "primary_post": continue
        g_lbl = c["group"]; t_lbl = c["t"]; b_lbl = c["b_lbl"]; g_rank = info["labels"].index(g_lbl); mu_g = info["mu_base"][g_rank]; lam_g = lam[g_rank]
        if cg_clean == "nevertreated":
            pool_ranks = [0]
        else:
            thresh = max(t_lbl, b_lbl) + anticipation; pool_ranks = [idx for idx, lbl in enumerate(info["labels"]) if idx == 0 or (lbl != g_lbl and lbl > thresh)]
        clambda = 0.0
        for x_val, px, qx in ((0, 1.0 - mu_g, q0), (1, mu_g, q1)):
            denom = sum(qx[r] for r in pool_ranks); num = sum(qx[r] * lam[r] for r in pool_ranks); clambda += px * (num / denom)
        lim = c["primary_truth"] + (lam_g - clambda) * (H[c["t_pos"]] - H[c["b_pos"]]); out[c["coord_key"]] = round(lim, 12)
    return out
def compute_mu_shifted_att(truth_bundle: dict, mu_override_by_rank: dict, extra_shift_fn=None):
    out = {}
    for c in truth_bundle["cells"]:
        if c["role"] != "primary_post": continue
        g_pos = c["g_pos"]; t_pos = c["t_pos"]; g_rank = truth_bundle["info"]["codes_pos"].index(g_pos); mu_new = mu_override_by_rank[g_rank]; theta_g = 1.0 + 0.1 * (g_pos - 3.0); val = theta_g + 0.25 * (t_pos - g_pos) + 0.5 * mu_new
        if extra_shift_fn is not None:
            val += extra_shift_fn(c, g_rank, mu_new)
        out[c["coord_key"]] = round(val, 12)
    return out
def compute_secondary_reference(recipe_id: str, recipe_row: dict, cell_id: str, call_type: str, control_group: str, anticipation: int, weightsname: str, fix_weights: str, truth_bundle: dict, ):
    """Compute secondary comparator reference or cell secondary_target diagnostic."""
    info = truth_bundle["info"]; T = len(info["calendar"]); weighted = bool(weightsname and weightsname != "NULL")
    # 1. Weighted profiles: unweighted reference is secondary diagnostic only
    if weighted:
        if recipe_id == "DATED_WEIGHTS" and call_type in ("DID_UNBALANCED", "DID_RCS"):
            return ("UNKNOWN_SECONDARY_LIMIT", {}, "unknown_limit")
        if recipe_id == "DATED_WEIGHTS" and call_type == "DID_CC":
            return ("UNAVAILABLE_NO_COMPLETE_UNITS", {}, "unavailable_no_complete_units")
        return ("UNWEIGHTED_DIAGNOSTIC_ATT", dict(truth_bundle["unweighted_diag_coords"]), "unweighted_diagnostic_for_weighted_study")
    # 2. Cell-specific secondary targets on FIT calls
    if call_type == "FIT":
        if cell_id == "X_NONE_SELECTIVE":
            r0, r1 = 0.35, 0.65; mu_map = {r: (m * (r1 ** 2)) / ((1.0 - m) * (r0 ** 2) + m * (r1 ** 2)) for r, m in enumerate(info["mu_base"])}
            return ("PAIR_SELECTED_ATT", compute_mu_shifted_att(truth_bundle, mu_map), "known_biased_limit")
        if cell_id == "X_TRANSPORT_A":
            r0, r1 = 0.35, 0.65; uB0, uB1 = 1.0 - (1.0 - r0) ** T, 1.0 - (1.0 - r1) ** T; uA0 = uB0 - T * r0 * ((1.0 - r0) ** (T - 1)); uA1 = uB1 - T * r1 * ((1.0 - r1) ** (T - 1))
            mu_map = {r: (m * uA1) / ((1.0 - m) * uA0 + m * uA1) for r, m in enumerate(info["mu_base"])}; return ("A_COMPOSITION_ATT", compute_mu_shifted_att(truth_bundle, mu_map), "known_biased_limit")
        if cell_id == "X_TRANSPORT_B":
            r0, r1 = 0.35, 0.65; uB0, uB1 = 1.0 - (1.0 - r0) ** T, 1.0 - (1.0 - r1) ** T; mu_map = {r: (m * uB1) / ((1.0 - m) * uB0 + m * uB1) for r, m in enumerate(info["mu_base"])}
            return ("B_COMPOSITION_ATT", compute_mu_shifted_att(truth_bundle, mu_map), "known_biased_limit")
        if cell_id == "PT_POST_INVALID":
            coords = {c["coord_key"]: round(c["primary_truth"] + 0.3 * (c["t_pos"] - c["g_pos"] + 1), 12) for c in truth_bundle["cells"] if c["role"] == "primary_post"}
            return ("PT_POST_DRIFT_PSEUDO_TARGET", coords, "known_biased_limit")
        if cell_id == "ANTICIPATION_0":
            coords = {c["coord_key"]: round(c["primary_truth"] - 0.3, 12) for c in truth_bundle["cells"] if c["role"] == "primary_post"}; return ("OMITTED_ANTICIPATION_PSEUDO_TARGET", coords, "known_biased_limit")
        if cell_id == "W_DISTINCT_imputation_1":
            return ("POPULATION_ATT_X5_CANCELS", dict(truth_bundle["primary_coords"]), "same_as_primary")
        return ("", {}, "not_applicable")
    if call_type in ("AGG", "PRETEST"): return ("", {}, "not_applicable")
    # 3. Comparator calls: DID_CC, DID_UNBALANCED, DID_RCS
    comp_recipes = {"COMP_adjacent", "COMP_arbitrary", "COMP_LADDER_01", "COMP_LADDER_02", "COMP_LADDER_03", "COMP_LADDER_04", "COMP_LADDER_05", "COMP_LADDER_06", "COMP_LADDER_07", "COMP_LADDER_08", "COMP_LATE", "GAP_LABELS",
        "CLUSTER_N2000_C100", "CLUSTER_N8000_C400", "CLUSTER_N2000_C20", "CLUSTER_N2000_C5", "MOMENT_t3", "MOMENT_t1.5", "GMM_SPARSE"}
    if recipe_id in comp_recipes:
        if call_type == "DID_CC": return ("UNAVAILABLE_NO_COMPLETE_UNITS", {}, "unavailable_no_complete_units")
        edges_kind = "adjacent" if recipe_id == "COMP_adjacent" else "arbitrary"; return (f"RCS_COMP_{edges_kind.upper()}", compute_rcs_comp_coords(truth_bundle, control_group, anticipation, edges_kind), "known_biased_limit")
    if recipe_id in ("CC_BASE", "CC_LADDER_01", "CC_LADDER_02", "CC_LADDER_03", "CC_LADDER_04"):
        c0, c1 = 0.2, 0.6
        if call_type == "DID_CC":
            mu_map = {r: (m * c1) / ((1.0 - m) * c0 + m * c1) for r, m in enumerate(info["mu_base"])}; return ("ATT_CC", compute_mu_shifted_att(truth_bundle, mu_map), "known_biased_limit")
        r0 = c0 + (1.0 - c0) * (2.0 / T); r1 = c1 + (1.0 - c1) * (2.0 / T); mu_map = {r: (m * r1) / ((1.0 - m) * r0 + m * r1) for r, m in enumerate(info["mu_base"])}
        return ("ATT_R", compute_mu_shifted_att(truth_bundle, mu_map), "known_biased_limit")
    if recipe_id == "CC_COMPOSITION":
        c0, c1 = 0.2, 0.6
        if call_type == "DID_CC":
            mu_map = {r: (m * c1) / ((1.0 - m) * c0 + m * c1) for r, m in enumerate(info["mu_base"])}; return ("ATT_CC", compute_mu_shifted_att(truth_bundle, mu_map), "known_biased_limit")
        r0 = c0 + (1.0 - c0) * (2.0 / T); r1 = c1 + (1.0 - c1) * (2.0 / T); mu_map = {r: (m * r1) / ((1.0 - m) * r0 + m * r1) for r, m in enumerate(info["mu_base"])}
        att_r = compute_mu_shifted_att(truth_bundle, mu_map); lam = [0.0, 0.6, 0.4, 0.2]; H = {d: ((T + 1.0) / 2.0 - d) / (T - 1.0) for d in range(1, T + 1)}; cg_clean = control_group.strip('"'); out = {}
        for c in truth_bundle["cells"]:
            if c["role"] != "primary_post": continue
            g_lbl = c["group"]; t_lbl = c["t"]; b_lbl = c["b_lbl"]; g_rank = info["labels"].index(g_lbl); mu_r = mu_map[g_rank]
            pool_ranks = [0] if cg_clean == "nevertreated" else [idx for idx, lbl in enumerate(info["labels"]) if idx == 0 or (lbl != g_lbl and lbl > max(t_lbl, b_lbl) + anticipation)]; shift = 0.0
            for x_val, px_r, cx, rx, qx in ((0, 1.0 - mu_r, c0, r0, info["q0"]), (1, mu_r, c1, r1, info["q1"]), ):
                scale_x = (1.0 - cx) * (2.0 / T) / rx; denom = sum(qx[r] for r in pool_ranks); lam_pool = sum(qx[r] * lam[r] for r in pool_ranks) / denom
                shift += px_r * scale_x * (lam[g_rank] - lam_pool) * (H[c["t_pos"]] - H[c["b_pos"]])
            out[c["coord_key"]] = round(att_r[c["coord_key"]] + shift, 12)
        return ("ATT_R_CC_COMPOSITION", out, "known_biased_limit")
    if recipe_id in ("BAL_BASE", "C_MCAR_DATE", "REPRESENTATIVE_ROSTERS", "BAL_ALIAS", "F_NO_INTERCEPT", "MCAR_HALF", "AR1_MCAR_DATE", "GMM_ZERO_SCALE",
        "CONT_PS_only", "CONT_OR_only", "CONT_both_correct", "EARLY_DRIFT", "TEST_null"):
        return ("POPULATION_ATT", dict(truth_bundle["primary_coords"]), "same_as_primary")
    if recipe_id == "C_MCAR_PAIR":
        if call_type == "DID_CC": return ("UNAVAILABLE_NO_COMPLETE_UNITS", {}, "unavailable_no_complete_units")
        return ("POPULATION_ATT", dict(truth_bundle["primary_coords"]), "same_as_primary")
    if recipe_id == "DATE_BASE":
        r0, r1 = 0.35, 0.65; p_exp = T if call_type == "DID_CC" else 1; mu_map = {r: (m * (r1 ** p_exp)) / ((1.0 - m) * (r0 ** p_exp) + m * (r1 ** p_exp)) for r, m in enumerate(info["mu_base"])}
        tag = "ATT_CC_DATE_BASE" if call_type == "DID_CC" else "ATT_R_DATE_BASE"; return (tag, compute_mu_shifted_att(truth_bundle, mu_map), "known_biased_limit")
    if recipe_id == "DATE_NO_INTERCEPT":
        r0, r1 = 0.35, 0.65; p_exp = T if call_type == "DID_CC" else 1; mu_map = {r: (m * (r1 ** p_exp)) / ((1.0 - m) * (r0 ** p_exp) + m * (r1 ** p_exp)) for r, m in enumerate(info["mu_base"])}; out = {}
        for c in truth_bundle["cells"]:
            if c["role"] != "primary_post": continue
            g_rank = info["codes_pos"].index(c["g_pos"]); out[c["coord_key"]] = round(0.5 * mu_map[g_rank], 12)
        tag = "ATT_CC_DATE_NO_INTERCEPT" if call_type == "DID_CC" else "ATT_R_DATE_NO_INTERCEPT"; return (tag, out, "known_biased_limit")
    if recipe_id == "DATE_CONSTANT_EFFECT":
        return ("POPULATION_ATT_CONSTANT_EFFECT", dict(truth_bundle["primary_coords"]), "same_as_primary")
    if recipe_id == "DATE_HET":
        if call_type == "DID_CC":
            a_t = [-1.0, -0.6, -0.2, 0.2, 0.6]; w0 = math.prod(expit(a) for a in a_t); w1 = math.prod(expit(a + 0.6) for a in a_t); mu_map = {r: (m * w1) / ((1.0 - m) * w0 + m * w1) for r, m in enumerate(info["mu_base"])}
            return ("ATT_CC_DATE_HET", compute_mu_shifted_att(truth_bundle, mu_map), "known_biased_limit")
        return ("UNKNOWN_SECONDARY_LIMIT", {}, "unknown_limit")
    if recipe_id == "DATE_COHORT_SLOPES":
        b_g = [0.2, 0.6, 1.0, 0.8]; p_exp = T if call_type == "DID_CC" else 1; mu_map = {}
        for r, m in enumerate(info["mu_base"]):
            rg0 = expit(-0.5) ** p_exp; rg1 = expit(-0.5 + b_g[r]) ** p_exp; mu_map[r] = (m * rg1) / ((1.0 - m) * rg0 + m * rg1)
        tag = "ATT_CC_COHORT_SLOPES" if call_type == "DID_CC" else "ATT_R_COHORT_SLOPES"; return (tag, compute_mu_shifted_att(truth_bundle, mu_map), "known_biased_limit")
    if recipe_id in ("DATE_INTERACTION_Q", "DATE_DISTINCT_0", "DATE_DISTINCT_1", "DATE_WRONG_Q_J"):
        p_exp = T if call_type == "DID_CC" else 1; has_inter_q = recipe_id in ("DATE_INTERACTION_Q", "DATE_WRONG_Q_J"); has_x1x5_eff = (recipe_id == "DATE_WRONG_Q_J"); out = {}
        for c in truth_bundle["cells"]:
            if c["role"] != "primary_post": continue
            g_pos = c["g_pos"]; t_pos = c["t_pos"]; g_rank = info["codes_pos"].index(g_pos); mu_g = info["mu_base"][g_rank]; theta_g = 1.0 + 0.1 * (g_pos - 3.0); num, den = 0.0, 0.0
            for x_val, px in ((0.0, 1.0 - mu_g), (1.0, mu_g)):
                for z_val in (-0.5, 0.5):
                    p_xz = 0.5 * px; eta = -0.5 + 0.6 * x_val + 0.8 * z_val + (1.2 * x_val * z_val if has_inter_q else 0.0); w_xz = expit(eta) ** p_exp
                    tau_xz = theta_g + 0.25 * (t_pos - g_pos) + 0.5 * x_val + (0.8 * x_val * z_val if has_x1x5_eff else 0.0); num += p_xz * w_xz * tau_xz; den += p_xz * w_xz
            out[c["coord_key"]] = round(num / den, 12)
        tag = f"ATT_CC_{recipe_id}" if call_type == "DID_CC" else f"ATT_R_{recipe_id}"; return (tag, out, "known_biased_limit")
    if recipe_id == "INNOVATION_SELECTION":
        if call_type == "DID_CC": return ("UNAVAILABLE_NO_COMPLETE_UNITS", {}, "unavailable_no_complete_units")
        n_shift = {(3, 2): -0.075 / math.sqrt(math.pi), (3, 3): +0.075 / math.sqrt(math.pi)}; out = {}
        for c in truth_bundle["cells"]:
            if c["role"] != "primary_post": continue
            g_pos, t_pos, b_pos = c["g_pos"], c["t_pos"], c["b_pos"]; s_t = n_shift.get((g_pos, t_pos), 0.0); s_b = n_shift.get((g_pos, b_pos), 0.0); out[c["coord_key"]] = round(c["primary_truth"] + (s_t - s_b), 12)
        return ("RCS_INNOVATION_SELECTION", out, "known_biased_limit")
    if recipe_id == "PT_POST_INVALID":
        coords = {c["coord_key"]: round(c["primary_truth"] + 0.3 * (c["t_pos"] - c["g_pos"] + 1), 12) for c in truth_bundle["cells"] if c["role"] == "primary_post"}; return ("PT_POST_DRIFT_PSEUDO_TARGET", coords, "known_biased_limit")
    if recipe_id == "ANTICIPATION_ONE":
        if anticipation == 0:
            coords = {c["coord_key"]: round(c["primary_truth"] - 0.3, 12) for c in truth_bundle["cells"] if c["role"] == "primary_post"}; return ("OMITTED_ANTICIPATION_PSEUDO_TARGET", coords, "known_biased_limit")
        return ("POPULATION_ATT_ANTICIPATION_1", dict(truth_bundle["primary_coords"]), "same_as_primary")
    if recipe_id == "TEST_conditional_alternative":
        return ("POPULATION_ATT_MEAN_ZERO_COND_ALT", dict(truth_bundle["primary_coords"]), "same_as_primary")
    if recipe_id == "SUPPORT_DISCONNECTED": return ("UNAVAILABLE_DISCONNECTED_SUPPORT", {}, "unavailable_no_complete_units")
    if recipe_id == "ATTRITION_MTAR":
        if call_type == "DID_CC":
            gam = [0.0, 0.6, 0.4, 0.2]; mu_map = {}
            for r, m in enumerate(info["mu_base"]):
                rho0 = 0.5 * (expit(gam[r]) + expit(-gam[r])); rho1 = 0.5 * (expit(0.6 + gam[r]) + expit(0.6 - gam[r])); mu_map[r] = (m * rho1) / ((1.0 - m) * rho0 + m * rho1)
            return ("ATT_CC_ATTRITION_MTAR", compute_mu_shifted_att(truth_bundle, mu_map), "known_biased_limit")
        return ("UNKNOWN_SECONDARY_LIMIT", {}, "unknown_limit")
    return ("UNKNOWN_SECONDARY_LIMIT", {}, "unknown_limit")
def moment_regime_for_recipe(recipe_id: str) -> str:
    if recipe_id == "MOMENT_t3": return "finite_variance_only_t3"
    if recipe_id == "MOMENT_t1.5": return "infinite_variance_t1_5"
    if recipe_id == "GMM_ZERO_SCALE": return "zero_variance_deterministic"
    return "finite_fourth_moment"
def verify_reference_files(root_dir: Path):
    """Verify exact rational agreement with reference/CHECKS.json and POPULATION_TARGETS.csv."""
    ref_dir = root_dir / "reference"; checks_path = ref_dir / "CHECKS.json"; targets_path = ref_dir / "POPULATION_TARGETS.csv"
    assert checks_path.exists() and targets_path.exists(), "Missing reference/ verification files"; F = fractions.Fraction; T = 5; Q = {0: (F(4, 10), F(3, 10), F(2, 10), F(1, 10)), 1: (F(2, 10), F(2, 10), F(3, 10), F(3, 10))}
    groups = [0, 3, 4, 5]; lam = {0: F(0), 3: F(6, 10), 4: F(4, 10), 5: F(2, 10)}
    def tau(g, t, x):
        return F(0) if g == 0 or t < g else F(1) + F(g - 3, 10) + F(t - g, 4) + F(x, 2)
    def y_fn(g, t, x, a):
        return F(a) + F(t, 10) + F(3, 10) * t * x + tau(g, t, x)
    computed_checks = []; computed_rows = []
    for family, edges in (("adjacent", [(s, s + 1) for s in range(1, T)]), ("all_pairs", list(itertools.combinations(range(1, T + 1), 2))), ):
        K = len(edges); h = {e: F(T + 1 - s - t, T - 2) for e in edges for s, t in [e]}; states = []
        for x in (0, 1):
            for g, pg in zip(groups, Q[x]):
                for a in (-1, 1):
                    for e in edges:
                        p = (1 + lam[g] * a * h[e]) / K; states.append({"x": x, "g": g, "a": a, "e": e, "p": F(1, 2) * pg * F(1, 2) * p})
        def avg(pred, val):
            z = [s for s in states if pred(s)]; d = sum(s["p"] for s in z); return sum(s["p"] * val(s) for s in z) / d
        eq_count = 0
        for g in groups:
            for x in (0, 1):
                for e in edges:
                    s, t = e; sel = avg(lambda z: z["g"] == g and z["x"] == x and z["e"] == e, lambda z: y_fn(g, t, x, z["a"]) - y_fn(g, s, x, z["a"]))
                    pop = avg(lambda z: z["g"] == g and z["x"] == x, lambda z: y_fn(g, t, x, z["a"]) - y_fn(g, s, x, z["a"]))
                    assert sel == pop, f"Conditional pair-change equality failed for {family}, g={g}, x={x}, e={e}"; eq_count += 1
        computed_checks.append({"family": family, "finite_states": len(states), "conditional_pair_change_equalities": eq_count })
    ref_checks = json.loads(checks_path.read_text(encoding="utf-8")); assert computed_checks[0]["conditional_pair_change_equalities"] == 32; assert computed_checks[1]["conditional_pair_change_equalities"] == 80
    for c_got, c_exp in zip(computed_checks, ref_checks["checks"]):
        assert c_got["family"] == c_exp["family"]; assert c_got["finite_states"] == c_exp["finite_states"]; assert c_got["conditional_pair_change_equalities"] == c_exp["conditional_pair_change_equalities"]
    # Also verify our compute_cell_primary_and_aux_truth and compute_rcs_comp_coords against POPULATION_TARGETS.csv
    ref_rows = list(csv.DictReader(targets_path.open(encoding="utf-8"))); assert len(ref_rows) == 24; dummy_recipe_row = {"n_dates": "5", "calendar": "1:T"}
    for r in ref_rows:
        rec_id = "COMP_adjacent" if r["family"] == "adjacent" else "COMP_arbitrary"; edges_kind = "adjacent" if r["family"] == "adjacent" else "arbitrary"
        tb = compute_cell_primary_and_aux_truth(rec_id, dummy_recipe_row, '"varying"', 0, "", ""); rcs_map = compute_rcs_comp_coords(tb, f'"{r["control_group"]}"', 0, edges_kind)
        ck = f"gt:{r['group']}:{r['time']}"; got_att = tb["primary_coords"][ck]; got_rcs = rcs_map[ck]; exp_att = float(F(r["exact_ATT"])); exp_rcs = float(F(r["exact_rcs_limit"]))
        assert abs(got_att - exp_att) < 1e-11, f"Primary ATT mismatch on {r}: {got_att} vs {exp_att}"; assert abs(got_rcs - exp_rcs) < 1e-11, f"RCS_COMP mismatch on {r}: {got_rcs} vs {exp_rcs}"
    return {"adjacent_pair_equalities": 32, "all_pairs_pair_equalities": 80, "rcs_comp_targets_verified": 24,}
def verify_generate_r_dispatch(campaign_dir: Path, recipes_by_id: dict):
    gen_path = campaign_dir / "generate.R"; text = gen_path.read_text(encoding="utf-8"); found_specs = set(re.findall(r"^\s*([A-Za-z0-9_.]+)\s*=\s*spec\(", text, flags=re.MULTILINE))
    expected_specs = set(recipes_by_id.keys()); missing = sorted(expected_specs - found_specs); extra = sorted(found_specs - expected_specs)
    if missing or extra: raise AssertionError(f"generate.R recipe dispatch mismatch: missing={missing}, extra={extra}")
    if "refusing fallback to any base law" not in text:
        raise AssertionError("generate.R is missing explicit unknown-recipe refusal guard")
    return {"n_recipes_in_contract": len(expected_specs), "n_recipes_in_generate_r": len(found_specs), "missing": missing, "extra": extra}
def build_campaign_manifest(root_dir: Path, campaign_dir: Path):
    contract_dir = campaign_dir / "contract"; recipes_rows = list(csv.DictReader((contract_dir / "RECIPES.csv").open(encoding="utf-8"))); cells_rows = list(csv.DictReader((contract_dir / "CELLS.csv").open(encoding="utf-8")))
    subsumed_rows = list(csv.DictReader((contract_dir / "SUBSUMED_MAP.csv").open(encoding="utf-8")))
    assert len(recipes_rows) == 64, f"Expected 64 recipes, got {len(recipes_rows)}"; assert len(cells_rows) == 130, f"Expected 130 targeted cells, got {len(cells_rows)}"
    assert len(subsumed_rows) == 49, f"Expected 49 subsumed rows, got {len(subsumed_rows)}"; recipes_by_id = {r["recipe_id"]: r for r in recipes_rows}; recipe_order = {r["recipe_id"]: idx + 1 for idx, r in enumerate(recipes_rows)}
    streams = {int(r["stream_index"]): r["stream_key"] for r in recipes_rows}; assert len(streams) == 59 and set(streams.keys()) == set(range(1, 60)), "Expected 59 RNG streams 1..59"
    dispatch_check = verify_generate_r_dispatch(campaign_dir, recipes_by_id); ref_check = verify_reference_files(campaign_dir)
    fixed_formals = json.loads((campaign_dir / "reference/formals.json").read_text())
    for name, expected in fixed_formals.items():
        assert list(globals()[name]) == expected, f"Fixed formal mismatch: {name}"
    # 1. Expand 9 core48 recipes x 48 cells = 432 logical cells (CAMPAIGN_DESIGN.md Section 6.1)
    logical_cells = []; core48_ids = set()
    for r in recipes_rows:
        if r["core_sweep"].strip() != "core48": continue
        rec_id = r["recipe_id"]
        for cg_code, cg_val in (("NT", '"nevertreated"'), ("NYT", '"notyettreated"')):
            for pt_code, pt_val, pt_start in (("post", '"post"', None), ("all", '"all"', "1")):
                for bp_code, bp_val in (("var", '"varying"'), ("uni", '"universal"')):
                    for w_code, w_val in (("id", '"identity"'), ("opt", '"optimal"')):
                        for est in ("ipw", "reg", "dr"):
                            cid = f"{rec_id}__{cg_code}_{pt_code}_{bp_code}_{w_code}_{est}"; core48_ids.add(cid); arg_pieces = [f"control_group={cg_val}", f"pt_assumption={pt_val}",]
                            if pt_start is not None:
                                arg_pieces.append(f"pt_start={pt_start}")
                            arg_pieces.extend([f"base_period={bp_val}", f"gmm_weighting={w_val}", f'est_method="{est}"' ])
                            logical_cells.append({"cell_id": cid, "recipe_id": rec_id, "family": r["family"], "logical_source": "core48", "role": "PRIMARY-valid", "calls": "FIT", "args": ", ".join(arg_pieces), "children": "AGG_S", "expected_validity": "VALID", "expected_disposition": "success", "expected_condition": "",
                                "primary_target_override": "", "secondary_target": "", "comparators": "keyed", "equality_ref": "" })
    assert len(core48_ids) == 432, f"Expected 432 core48 cells, got {len(core48_ids)}"
    # Verify SUBSUMED_MAP.csv against core48_ids and CELLS.csv
    targeted_ids = {r["cell_id"] for r in cells_rows}
    for srow in subsumed_rows:
        assert srow["core_sweep_cell_id"] in core48_ids, f"Missing core48 target {srow['core_sweep_cell_id']}"; assert srow["input_cell_id"] not in targeted_ids, f"Subsumed cell {srow['input_cell_id']} still in CELLS.csv"
    # Append 130 targeted cells
    for crow in cells_rows:
        logical_cells.append({"cell_id": crow["cell_id"], "recipe_id": crow["recipe_id"], "family": crow["family"], "logical_source": "targeted", "role": crow["role"], "calls": crow["calls"], "args": crow["args"],
            "children": crow["children"], "expected_validity": crow["expected_validity"], "expected_disposition": crow["expected_disposition"], "expected_condition": crow["expected_condition"],
            "primary_target_override": crow["primary_target_override"], "secondary_target": crow["secondary_target"], "comparators": crow["comparators"], "equality_ref": crow["equality_ref"] })
    assert len({c["cell_id"] for c in logical_cells}) == 562, "Duplicate logical ID"
    assert len(logical_cells) == 562, f"Expected 562 total logical cells, got {len(logical_cells)}"
    # Sort logical cells by (stream_index, recipe_order, logical_cell_index)
    for idx, lc in enumerate(logical_cells):
        r = recipes_by_id[lc["recipe_id"]]; lc["_stream_index"] = int(r["stream_index"]); lc["_recipe_order"] = recipe_order[lc["recipe_id"]]; lc["_cell_order"] = idx + 1
    logical_cells_sorted = sorted(logical_cells, key=lambda c: (c["_stream_index"], c["_recipe_order"], c["_cell_order"]),); manifest_rows = []
    targets_payload = {"schema": "cdid_portable_campaign_targets_v1", "seed_master": SEED_MASTER, "n_streams": N_STREAMS, "n_replicates": N_REPS, "n_purposes": N_PURPOSES,
        "subsumed_map": {s["input_cell_id"]: s["core_sweep_cell_id"] for s in subsumed_rows}, "tasks": {},}
    # Group logical cells by recipe to append deduplicated recipe-level comparators right after the recipe's cells
    cells_by_recipe = {}
    for lc in logical_cells_sorted:
        cells_by_recipe.setdefault(lc["recipe_id"], []).append(lc)
    recipes_sorted = sorted(recipes_rows, key=lambda r: (int(r["stream_index"]), recipe_order[r["recipe_id"]]),); n_fit = 0; n_agg = 0; n_pretest = 0; n_did_cc = 0; n_did_unbal = 0; n_did_rcs = 0
    for r in recipes_sorted:
        rec_id = r["recipe_id"]; s_idx = int(r["stream_index"]); s_key = r["stream_key"]; m_regime = moment_regime_for_recipe(rec_id); rec_cells = cells_by_recipe.get(rec_id, [])
        # Track comparator deduplication per recipe
        comp_cc_unbal_groups = {}; comp_rcs_groups = {}
        for lc in rec_cells:
            cid = lc["cell_id"]; flags, fit_args = build_effective_fit_args(r, lc["args"]); lit_fit = format_r_call("cdid::att_gt", fit_args, CDID_ATT_GT_FORMALS)
            bp = fit_args.get("base_period", '"varying"'); ant = int(fit_args.get("anticipation", "0")); wname = fit_args.get("weightsname", ""); fixw = fit_args.get("fix_weights", "")
            cg = fit_args.get("control_group", '"nevertreated"'); comp_inffunc = fit_args.get("compute_inffunc", "TRUE") == "TRUE"
            panel_true = fit_args.get("panel", "TRUE") == "TRUE"; tb = compute_cell_primary_and_aux_truth(rec_id, r, bp, ant, wname, fixw); p_target_id = primary_target_identifier(rec_id, wname, fixw, "FIT")
            sec_id, sec_coords, sec_status = compute_secondary_reference(rec_id, r, cid, "FIT", cg, ant, wname, fixw, tb); is_refusal = (lc["role"] == "expected-refusal") or (lc["expected_disposition"] == "expected-refusal")
            exp_disp = "expected-refusal" if is_refusal else "success"; exp_cond = lc["expected_condition"].strip() or ("cdid_error" if is_refusal else ""); inf_elig = comp_inffunc and panel_true
            t_match = ("refusal" if is_refusal else ("matched" if lc["expected_validity"].startswith("VALID") or lc["role"] in (
                "PRIMARY-valid", "valid-design", "valid-branch", "calibration", "cluster-sensitivity", "few-cluster-diagnostic", "valid-heavy-tail-sensitivity",
                "moment-failure", "valid-null", "catalog-support-negative" ) else "mismatched"))
            fit_task_id = f"{cid}__FIT"; children_specs = expand_children_specs(lc["children"]); p_coord_keys = list(tb["primary_coords"].keys())
            fit_row = {"task_id": fit_task_id, "replicate_scope": f"1..{N_REPS}", "stream_index": s_idx, "stream_key": s_key, "recipe_id": rec_id, "family": lc["family"], "cell_id": cid,
                "logical_source": lc["logical_source"], "role": lc["role"], "call_type": "FIT", "package": "cdid", "parent_task_id": "", "agg_kind": "", "literal_call": lit_fit, "rng_purpose": 3,
                "primary_target_id": p_target_id, "primary_truth": json.dumps(tb["primary_coords"], separators=(",", ":")), "primary_population_masses": json.dumps(tb["masses_map"], separators=(",", ":")),
                "primary_coordinate_set": ";".join(p_coord_keys), "target_match": t_match, "causal_validity": lc["expected_validity"] or "VALID", "inference_eligible": "TRUE" if inf_elig else "FALSE", "expected_disposition": exp_disp,
                "expected_condition": exp_cond, "secondary_reference_id": sec_id, "secondary_reference": json.dumps(sec_coords, separators=(",", ":")) if sec_coords else "",
                "secondary_status": sec_status, "moment_regime": m_regime, "equality_ref": lc["equality_ref"], "linked_cell_ids": cid,}
            manifest_rows.append(fit_row)
            targets_payload["tasks"][fit_task_id] = {**fit_row, "primary_truth_map": tb["primary_coords"], "auxiliary_truth_map": tb["aux_coords"], "secondary_reference_map": sec_coords, "population_masses_map": tb["masses_map"],
                "save_parent_rds": bool(children_specs), "point_only_rcs": not panel_true, "wald_pretest_null_true": all(v == 0 for v in tb["aux_coords"].values())}
            n_fit += 1
            # Optional PRETEST call on the same logical cell
            call_tokens = [t.strip() for t in lc["calls"].split(";") if t.strip()]
            if "PRETEST" in call_tokens:
                pre_args = build_pretest_args(fit_args); lit_pre = format_r_call("cdid::conditional_did_pretest", pre_args, CDID_PRETEST_FORMALS); pre_task_id = f"{cid}__PRETEST"
                pre_exp_disp = "not-applicable" if cid == "F_CALLBACK_EXISTING" else "success"; pre_truth = {"conditional_pt_null_true": (rec_id == "TEST_null"), "unconditional_mean_zero": True,}
                pre_row = {"task_id": pre_task_id, "replicate_scope": f"1..{N_REPS}", "stream_index": s_idx, "stream_key": s_key, "recipe_id": rec_id, "family": lc["family"], "cell_id": cid,
                    "logical_source": lc["logical_source"], "role": lc["role"], "call_type": "PRETEST", "package": "cdid", "parent_task_id": "", "agg_kind": "", "literal_call": lit_pre, "rng_purpose": 5,
                    "primary_target_id": "PRETEST_CONDITIONAL_PT_NULL", "primary_truth": json.dumps(pre_truth, separators=(",", ":")), "primary_population_masses": json.dumps(tb["masses_map"], separators=(",", ":")),
                    "primary_coordinate_set": "CvMpval", "target_match": "pretest", "causal_validity": lc["expected_validity"] or "VALID", "inference_eligible": "FALSE" if cid == "F_CALLBACK_EXISTING" else "TRUE",
                    "expected_disposition": pre_exp_disp, "expected_condition": "", "secondary_reference_id": "",
                    "secondary_reference": "", "secondary_status": "not_applicable", "moment_regime": m_regime, "equality_ref": "", "linked_cell_ids": cid,}
                manifest_rows.append(pre_row)
                targets_payload["tasks"][pre_task_id] = {**pre_row, "primary_truth_map": pre_truth, "auxiliary_truth_map": {}, "secondary_reference_map": {}, "population_masses_map": tb["masses_map"], "save_parent_rds": False}
                n_pretest += 1
            # Child AGG calls
            for agg_suffix, agg_kind, agg_arg_dict in children_specs:
                lit_agg = format_r_call("cdid::aggte", agg_arg_dict, CDID_AGGTE_FORMALS); agg_task_id = f"{cid}__{agg_suffix}"; agg_tb = compute_aggte_truth_py(tb, agg_kind)
                if cid == "I_NO_IF":
                    agg_exp_disp = "expected-refusal"; agg_exp_cond = "cdid_chained_request_no_inffunc_error"; agg_inf_elig = False
                else:
                    agg_exp_disp = "success"; agg_exp_cond = ""; agg_inf_elig = inf_elig
                agg_row = {"task_id": agg_task_id, "replicate_scope": f"1..{N_REPS}", "stream_index": s_idx, "stream_key": s_key, "recipe_id": rec_id, "family": lc["family"], "cell_id": cid,
                    "logical_source": lc["logical_source"], "role": lc["role"], "call_type": "AGG", "package": "cdid", "parent_task_id": fit_task_id, "agg_kind": agg_kind, "literal_call": lit_agg, "rng_purpose": 4,
                    "primary_target_id": p_target_id, "primary_truth": json.dumps(agg_tb["coordinates"], separators=(",", ":")), "primary_population_masses": json.dumps(tb["masses_map"], separators=(",", ":")),
                    "primary_coordinate_set": ";".join(agg_tb["coordinate_keys"]), "target_match": t_match, "causal_validity": lc["expected_validity"] or "VALID",
                    "inference_eligible": "TRUE" if agg_inf_elig else "FALSE", "expected_disposition": agg_exp_disp, "expected_condition": agg_exp_cond, "secondary_reference_id": "", "secondary_reference": "",
                    "secondary_status": "not_applicable", "moment_regime": m_regime, "equality_ref": lc["equality_ref"], "linked_cell_ids": cid,}
                manifest_rows.append(agg_row)
                targets_payload["tasks"][agg_task_id] = {**agg_row, "primary_truth_map": agg_tb["coordinates"], "auxiliary_truth_map": agg_tb.get("auxiliary_coordinates", {}),
                    "secondary_reference_map": {}, "population_masses_map": tb["masses_map"], "save_parent_rds": False}
                n_agg += 1
            # Register comparator keys for keyed non-refusal, non-callback cells
            if lc["comparators"].strip() == "keyed" and not is_refusal:
                up_key = upstream_key(fit_args)
                if up_key not in comp_cc_unbal_groups:
                    comp_cc_unbal_groups[up_key] = {"fit_args": dict(fit_args), "tb": tb, "p_target_id": p_target_id, "p_coord_keys": p_coord_keys, "cell_ids": []}
                comp_cc_unbal_groups[up_key]["cell_ids"].append(cid)
                if fit_args.get("est_method") == '"dr"' and fit_args.get("base_period") == '"varying"':
                    rcs_key = upstream_key(fit_args, rcs=True)
                    if rcs_key not in comp_rcs_groups:
                        comp_rcs_groups[rcs_key] = {"fit_args": dict(fit_args), "tb": tb, "p_target_id": p_target_id, "p_coord_keys": p_coord_keys, "cell_ids": []}
                    comp_rcs_groups[rcs_key]["cell_ids"].append(cid)
        # Emit deduplicated DID_CC and DID_UNBALANCED tasks for this recipe
        has_complete_units = not r["did_cc_secondary_reference"].strip().startswith("unavailable")
        for idx_k, (up_key, grp) in enumerate(comp_cc_unbal_groups.items(), start=1):
            f_args = grp["fit_args"]; tb = grp["tb"]; p_target_id = grp["p_target_id"]; p_coord_keys = grp["p_coord_keys"]; linked_cids = ";".join(grp["cell_ids"])
            bp = f_args.get("base_period", '"varying"'); ant = int(f_args.get("anticipation", "0")); wname = f_args.get("weightsname", ""); fixw = f_args.get("fix_weights", ""); cg = f_args.get("control_group", '"nevertreated"')
            comp_inffunc = f_args.get("compute_inffunc", "TRUE") == "TRUE"; variants = [("DID_UNBALANCED", "TRUE")]
            if has_complete_units:
                variants.insert(0, ("DID_CC", "FALSE"))
            for call_type, allow_unbal_str in variants:
                did_args = effective_upstream_args(f_args)
                did_args["panel"] = "TRUE"; did_args["allow_unbalanced_panel"] = allow_unbal_str; lit_did = format_r_call("did::att_gt", did_args, DID_ATT_GT_FORMALS)
                sec_id, sec_coords, sec_status = compute_secondary_reference(rec_id, r, grp["cell_ids"][0], call_type, cg, ant, wname, fixw, tb); t_match = "matched" if sec_status == "same_as_primary" else "mismatched"
                c_valid = "VALID" if sec_status == "same_as_primary" else "NO_GUARANTEE"
                if rec_id == "CONT_PS_only" and f_args.get("est_method") == '"reg"':
                    sec_id, sec_coords, sec_status = "UNKNOWN_SECONDARY_LIMIT", {}, "unknown_limit"
                    t_match = "unknown"
                if rec_id == "CONT_PS_only" and f_args.get("est_method") == '"reg"': c_valid = "NO_GUARANTEE"
                if rec_id == "CONT_OR_only" and f_args.get("est_method") == '"ipw"': c_valid = "NO_GUARANTEE"
                did_task_id = f"{rec_id}__{call_type}__{idx_k:02d}"
                did_row = {"task_id": did_task_id, "replicate_scope": f"1..{N_REPS}", "stream_index": s_idx, "stream_key": s_key, "recipe_id": rec_id, "family": r["family"], "cell_id": grp["cell_ids"][0],
                    "logical_source": "comparator", "role": "comparator", "call_type": call_type, "package": "did", "parent_task_id": "", "agg_kind": "", "literal_call": lit_did, "rng_purpose": 3, "primary_target_id": p_target_id,
                    "primary_truth": json.dumps(tb["primary_coords"], separators=(",", ":")), "primary_population_masses": json.dumps(tb["masses_map"], separators=(",", ":")),
                    "primary_coordinate_set": ";".join(p_coord_keys), "target_match": t_match, "causal_validity": c_valid, "inference_eligible": "TRUE" if comp_inffunc else "FALSE",
                    "expected_disposition": "success", "expected_condition": "", "secondary_reference_id": sec_id, "secondary_reference": json.dumps(sec_coords, separators=(",", ":")) if sec_coords else "",
                    "secondary_status": sec_status, "moment_regime": m_regime, "equality_ref": "", "linked_cell_ids": linked_cids,}
                manifest_rows.append(did_row)
                targets_payload["tasks"][did_task_id] = {**did_row, "primary_truth_map": tb["primary_coords"], "auxiliary_truth_map": tb["aux_coords"],
                    "secondary_reference_map": sec_coords, "population_masses_map": tb["masses_map"], "save_parent_rds": False, "point_only_rcs": False}
                if call_type == "DID_CC":
                    n_did_cc += 1
                else:
                    n_did_unbal += 1
        # Emit deduplicated explicit point-only DID_RCS tasks for this recipe
        for idx_r, (rcs_key, grp) in enumerate(comp_rcs_groups.items(), start=1):
            f_args = grp["fit_args"]; tb = grp["tb"]; p_target_id = grp["p_target_id"]; p_coord_keys = grp["p_coord_keys"]; linked_cids = ";".join(grp["cell_ids"])
            ant = int(f_args.get("anticipation", "0")); wname = f_args.get("weightsname", ""); fixw = f_args.get("fix_weights", ""); cg = f_args.get("control_group", '"nevertreated"')
            rcs_args = effective_upstream_args(f_args, rcs=True)
            lit_rcs = format_r_call("did::att_gt", rcs_args, DID_ATT_GT_FORMALS); sec_id, sec_coords, sec_status = compute_secondary_reference(rec_id, r, grp["cell_ids"][0], "DID_RCS", cg, ant, wname, fixw, tb)
            t_match = "matched" if sec_status == "same_as_primary" else "mismatched"; c_valid = "VALID" if sec_status == "same_as_primary" else "NO_GUARANTEE"; rcs_task_id = f"{rec_id}__DID_RCS__{idx_r:02d}"
            rcs_row = {"task_id": rcs_task_id, "replicate_scope": f"1..{N_REPS}", "stream_index": s_idx, "stream_key": s_key, "recipe_id": rec_id, "family": r["family"], "cell_id": grp["cell_ids"][0],
                "logical_source": "comparator", "role": "comparator-point-only", "call_type": "DID_RCS", "package": "did", "parent_task_id": "", "agg_kind": "", "literal_call": lit_rcs, "rng_purpose": 3,
                "primary_target_id": p_target_id, "primary_truth": json.dumps(tb["primary_coords"], separators=(",", ":")), "primary_population_masses": json.dumps(tb["masses_map"], separators=(",", ":")),
                "primary_coordinate_set": ";".join(p_coord_keys), "target_match": t_match, "causal_validity": c_valid,
                "inference_eligible": "FALSE", "expected_disposition": "success", "expected_condition": "", "secondary_reference_id": sec_id, "secondary_reference": json.dumps(sec_coords, separators=(",", ":")) if sec_coords else "",
                "secondary_status": sec_status, "moment_regime": m_regime, "equality_ref": "", "linked_cell_ids": linked_cids,}
            manifest_rows.append(rcs_row)
            targets_payload["tasks"][rcs_task_id] = {**rcs_row, "primary_truth_map": tb["primary_coords"], "auxiliary_truth_map": tb["aux_coords"],
                "secondary_reference_map": sec_coords, "population_masses_map": tb["masses_map"], "save_parent_rds": False, "point_only_rcs": True}
            n_did_rcs += 1
    counts = {"n_recipes": len(recipes_rows), "n_streams": len(streams), "n_replicates": N_REPS, "n_core48_cells": 432, "n_targeted_cells": 130, "n_subsumed_cells": 49, "n_logical_cells": len(logical_cells),
        "per_replicate": { "FIT": n_fit, "PRETEST": n_pretest, "AGG": n_agg, "cdid_total": n_fit + n_pretest + n_agg, "DID_CC": n_did_cc, "DID_UNBALANCED": n_did_unbal,
        "DID_RCS": n_did_rcs, "did_total": n_did_cc + n_did_unbal + n_did_rcs, "all_tasks_per_rep": len(manifest_rows), }, "campaign_total_100_reps": { "FIT": n_fit * N_REPS,
        "PRETEST": n_pretest * N_REPS, "AGG": n_agg * N_REPS, "cdid_total": (n_fit + n_pretest + n_agg) * N_REPS, "DID_CC": n_did_cc * N_REPS, "DID_UNBALANCED": n_did_unbal * N_REPS, "DID_RCS": n_did_rcs * N_REPS,
        "did_total": (n_did_cc + n_did_unbal + n_did_rcs) * N_REPS, "all_tasks_total": len(manifest_rows) * N_REPS, }, "dispatch_verification": dispatch_check, "reference_verification": ref_check}
    targets_payload["counts"] = counts; return manifest_rows, targets_payload
def write_manifest_and_targets(manifest_rows: list, targets_payload: dict, campaign_dir: Path):
    manifest_path = campaign_dir / "manifest.csv"; targets_path = campaign_dir / "targets.json"; fieldnames = list(manifest_rows[0].keys())
    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames); writer.writeheader(); writer.writerows(manifest_rows)
    targets_path.write_text(json.dumps(targets_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"); return manifest_path, targets_path
def compute_harness_hashes(campaign_dir: Path) -> dict:
    files = ["contract/CAMPAIGN_DESIGN.md", "contract/RECIPES.csv", "contract/CELLS.csv", "contract/SUBSUMED_MAP.csv", "generate.R", "worker.R", "freeze.py", "run_batch.py", "run.py", "summarize.py", "check.py", "manifest.csv", "targets.json", "README.md", "reference/CHECKS.json", "reference/POPULATION_TARGETS.csv", "reference/formals.json"]
    out = {}
    for rel in files:
        p = campaign_dir / rel
        if not p.is_file(): raise RuntimeError(f"Missing required toolkit file: {rel}")
        out[rel] = sha256_file(p)
    return out
def checkout_commit(repo):
    repo = Path(repo).resolve(); gd = repo / ".git"
    if not gd.is_dir(): raise RuntimeError("A clean checkout with local .git metadata is required")
    head = (gd / "HEAD").read_text().strip()
    if head.startswith("ref: "):
        ref = head[5:]
        if not ref.startswith("refs/") or ".." in Path(ref).parts: raise RuntimeError("Invalid HEAD reference")
        rp = gd / ref
        if rp.exists(): head = rp.read_text().strip()
        else:
            packed = (gd / "packed-refs").read_text().splitlines() if (gd / "packed-refs").exists() else []
            head = next((line.split()[0] for line in packed if line.endswith(" " + ref)), "")
    if not re.fullmatch(r"[0-9a-f]{40}", head): raise RuntimeError("Unresolved checkout commit")
    return head

def source_snapshot(repo):
    repo = Path(repo).resolve(); files = [repo / n for n in ("DESCRIPTION", "NAMESPACE")]
    for name in ("R", "src", "inst", "data", "man"):
        d = repo / name
        if d.exists(): files.extend(p for p in sorted(d.rglob("*")) if p.is_file())
    if not (repo / "R").is_dir(): raise RuntimeError("Missing package R source")
    return {p.relative_to(repo).as_posix(): sha256_file(p) for p in files}

def run_actual_freeze(args, root_dir: Path, campaign_dir: Path, manifest_rows: list, targets_payload: dict):
    if sys.version_info < (3, 9): raise RuntimeError(f"Python >= 3.9 is required (found {sys.version})")
    source_sha = (args.source_sha or "").strip()
    if not source_sha or source_sha == "FREEZE_SOURCE_SHA": raise RuntimeError(
        "Refusing freeze: --source-sha is missing or still set to placeholder FREEZE_SOURCE_SHA. " "Supply the reviewed private repository commit SHA at actual freeze time.")
    out_dir = Path(args.out).resolve(); out_dir.mkdir(parents=True, exist_ok=True); provenance_path = out_dir / "provenance.json"
    if any(out_dir.iterdir()):
        raise RuntimeError("Freeze requires an empty output directory; preserved runs cannot be overwritten")
    repo = Path(args.repo).resolve()
    if checkout_commit(repo) != source_sha: raise RuntimeError("Supplied source SHA differs from checkout HEAD")
    source_hashes = source_snapshot(repo)
    write_manifest_and_targets(manifest_rows, targets_payload, campaign_dir)
    # Also copy frozen manifest.csv and targets.json into out_dir
    (out_dir / "manifest.csv").write_bytes((campaign_dir / "manifest.csv").read_bytes()); (out_dir / "targets.json").write_bytes((campaign_dir / "targets.json").read_bytes())
    seeds_rds_path = out_dir / "seeds_29500.rds"; runtime_json_path = out_dir / "r_runtime.json"; freeze_cfg_path = out_dir / "freeze_request.json"; lib_paths = [str(Path(p).resolve()) for p in args.lib] if args.lib else []
    freeze_cfg = {"source_sha": source_sha, "repo_path": str(repo), "formals_path": str(campaign_dir / "reference/formals.json"), "lib_paths": lib_paths, "seeds_rds_path": str(seeds_rds_path), "runtime_out_json": str(runtime_json_path)}
    freeze_cfg_path.write_text(json.dumps(freeze_cfg, indent=2) + "\n", encoding="utf-8"); cmd = [args.rscript, "--vanilla", str(campaign_dir / "worker.R"), "--freeze", str(freeze_cfg_path)]
    before_harness = compute_harness_hashes(campaign_dir)
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=300)
    if proc.returncode != 0: raise RuntimeError(f"worker.R --freeze failed (exit {proc.returncode}):\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    r_runtime = json.loads(runtime_json_path.read_text(encoding="utf-8"))
    if not r_runtime.get("ok") or r_runtime.get("source_sha") != source_sha or r_runtime.get("n_states") != 29500:
        raise RuntimeError("Invalid R freeze receipt")
    if r_runtime.get("seeds_rds_hash") != "md5:" + hashlib.md5(seeds_rds_path.read_bytes()).hexdigest():
        raise RuntimeError("R/Python seed receipt hash mismatch")
    harness_hashes = compute_harness_hashes(campaign_dir)
    if harness_hashes != before_harness or source_snapshot(repo) != source_hashes or checkout_commit(repo) != source_sha:
        raise RuntimeError("Source/harness changed during freeze")
    provenance = {"schema": "cdid_portable_campaign_provenance_v1", "frozen_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source_sha": source_sha, "repo_path": str(Path(args.repo).resolve()) if args.repo else "", "lib_paths": lib_paths, "python_version": sys.version, "python_platform": platform.platform(),
        "seed_master": SEED_MASTER, "n_streams": N_STREAMS, "n_replicates": N_REPS, "n_purposes": N_PURPOSES, "seeds_rds_path": str(seeds_rds_path), "seeds_rds_sha256": sha256_file(seeds_rds_path),
        "source_hashes": source_hashes, "runtime_json_sha256": sha256_file(runtime_json_path), "harness_hashes": harness_hashes, "r_runtime": r_runtime, "counts": targets_payload["counts"],}
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"); return provenance
def main():
    parser = argparse.ArgumentParser(description="Build and verify campaign manifest, targets, and freeze state.")
    parser.add_argument("--check", action="store_true", help="Run static/rational truth and inventory checks (no R execution).")
    parser.add_argument("--write-manifest", action="store_true", help="Write deterministic manifest.csv and targets.json.")
    parser.add_argument("--freeze", action="store_true", help="Run full campaign freeze (requires real --source-sha and R packages)."); parser.add_argument("--repo", default="", help="Path to source checkout.")
    parser.add_argument("--lib", action="append", default=[], help="R library path(s) to prepend in worker.R."); parser.add_argument("--out", default="", help="Campaign output directory for freeze artifacts.")
    parser.add_argument("--source-sha", default="FREEZE_SOURCE_SHA", help="Committed source SHA (placeholder FREEZE_SOURCE_SHA is refused on --freeze).")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable for --freeze."); parser.add_argument("--force-new", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args(); campaign_dir = Path(__file__).resolve().parent; root_dir = campaign_dir.parent.parent
    manifest_rows, targets_payload = build_campaign_manifest(root_dir, campaign_dir)
    if args.write_manifest or args.check or not args.freeze:
        m_path, t_path = write_manifest_and_targets(manifest_rows, targets_payload, campaign_dir)
    if args.freeze:
        if not args.out or not args.repo or not args.lib: raise SystemExit("--out, --repo and --lib are required for freeze")
        prov = run_actual_freeze(args, root_dir, campaign_dir, manifest_rows, targets_payload); print(json.dumps({"status": "FROZEN", "source_sha": prov["source_sha"], "counts": prov["counts"]}, indent=2)); return
    print(json.dumps({"status": "CHECK_OK", "manifest_csv": str(campaign_dir / "manifest.csv"), "targets_json": str(campaign_dir / "targets.json"), "counts": targets_payload["counts"] }, indent=2))
if __name__ == "__main__":
    main()
