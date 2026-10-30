#!/usr/bin/env python3
"""
Mechanical stdlib-only JSON-to-CSV exporter and immutable CSV copy tool for
CP07 Chained DR independent reference fixtures (validation/reference/dr_finite_v1).

- Uses ONLY Python 3 standard library (csv, hashlib, json, os, shutil, sys).
- Verifies SHA-256 hashes of all frozen source files in validation/reference/dr_finite_v1.
- Never calls any production or reference estimator; never parses or reformats
  floating-point numbers (uses json.load(..., parse_float=str, parse_int=str)
  so every expected number string is preserved byte-for-byte).
- Copies atoms.csv, unit_if_rows.csv, nuisance_if_rows.csv, and
  population_examples.csv byte-for-byte and verifies identical SHA-256 hashes.
"""

import csv
import hashlib
import json
import os
import shutil
import sys

sys.dont_write_bytecode = True

EXPECTED_INPUT_HASHES = {
    "FIXED_CHECK_PLAN.md": "34a59ec3260b12d1d25ea6fac16d566372cc1620d6593b51d46225fc8173ba64",
    "REFERENCE_REPORT.md": "f37cb313f087eb07393111d41d0b7651a32fafc19f8fb528a83ccd4262657014",
    "atoms.csv": "f6869b76895e3c0c91b8d7152bf38f6dbff65ab5400f033c07730c378a03d8d1",
    "backend_fixture.json": "e36b4575d97b52ac5a691de07f3b6439c403836d0dcff917ad1f15ae3bbaf928",
    "check_results.json": "f0e4be45bf195823132ca09a5c48dd01fc284a018cbfe3ba0dfa610bf33d8a08",
    "generate_dr_reference.py": "7cb9f66c9ed6d250c734796508c18eedee6d8aca328614ed6751e20190c98282",
    "nuisance_if_rows.csv": "a607bf85c4603ef3ee0d276e4cb2e50c41bb0d5edc788b72de450ff44a68ad9b",
    "population_examples.csv": "0378ee9ff9c5dc85794a4680eb272d3705c862779eb31ba9092c91728b27f770",
    "unit_if_rows.csv": "aa5618592cd9bfa345c32cbd1b8af73a286c98e14bcd783778f592ead2356048",
}

COPIED_CSV_FILES = (
    "atoms.csv",
    "unit_if_rows.csv",
    "nuisance_if_rows.csv",
    "population_examples.csv",
)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(path, fieldnames, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    out_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(out_dir, "..", "..", "..", ".."))
    src_dir = os.path.join(repo_root, "validation", "reference", "dr_finite_v1")

    # 1. Verify all frozen input file SHA-256 pins
    verified_inputs = {}
    for fname, expected_sha in EXPECTED_INPUT_HASHES.items():
        fpath = os.path.join(src_dir, fname)
        if not os.path.isfile(fpath):
            raise RuntimeError(f"Missing frozen reference file: {fpath}")
        actual_sha = sha256_file(fpath)
        if actual_sha != expected_sha:
            raise RuntimeError(
                f"SHA-256 mismatch for {fname}: expected {expected_sha}, got {actual_sha}"
            )
        verified_inputs[fname] = actual_sha

    # 2. Byte-for-byte copy of independent CSV files
    for fname in COPIED_CSV_FILES:
        src_path = os.path.join(src_dir, fname)
        dst_path = os.path.join(out_dir, fname)
        shutil.copyfile(src_path, dst_path)
        if sha256_file(dst_path) != EXPECTED_INPUT_HASHES[fname]:
            raise RuntimeError(f"Copied CSV hash mismatch for {fname}")

    # 3. Parse JSON files preserving literal number strings byte-for-byte
    with open(os.path.join(src_dir, "backend_fixture.json"), "r", encoding="utf-8") as f:
        bf = json.load(f, parse_float=str, parse_int=str)
    with open(os.path.join(src_dir, "check_results.json"), "r", encoding="utf-8") as f:
        cr = json.load(f, parse_float=str, parse_int=str)

    # 4. Export calendar.csv
    calendar_rows = [
        {"calendar_index": str(i + 1), "time": str(d)}
        for i, d in enumerate(bf["calendar_dates"])
    ]
    write_csv(
        os.path.join(out_dir, "calendar.csv"),
        ["calendar_index", "time"],
        calendar_rows,
    )

    # 5. Export cohorts.csv
    cohort_rows = [
        {"cohort_index": str(i + 1), "cohort": str(c)}
        for i, c in enumerate(bf["cohorts"])
    ]
    write_csv(
        os.path.join(out_dir, "cohorts.csv"),
        ["cohort_index", "cohort"],
        cohort_rows,
    )

    # 6. Export catalog_links.csv
    link_rows = []
    for lk in bf["catalog_links"]:
        link_rows.append(
            {
                "link_id": lk["link_id"],
                "g": lk["g"],
                "K": "|".join(lk["K"]),
                "pair": lk["pair"],
                "s": lk["s"],
                "t": lk["t"],
                "L_e": lk["L_e"],
                "raw_gap": lk["raw_gap"],
                "or_block": lk["or_block"],
                "functional_key": f"cmpcat:{lk['link_id']}|g={lk['g']}|K={','.join(lk['K'])}|pair={lk['s']}->{lk['t']}",
            }
        )
    write_csv(
        os.path.join(out_dir, "catalog_links.csv"),
        [
            "link_id",
            "g",
            "K",
            "pair",
            "s",
            "t",
            "L_e",
            "raw_gap",
            "or_block",
            "functional_key",
        ],
        link_rows,
    )

    # 7. Export or_blocks.csv
    block_rows = []
    for blk in bf["or_blocks"]:
        block_rows.append(
            {
                "block_id": blk["block_id"],
                "pair": blk["pair"],
                "s": blk["s"],
                "t": blk["t"],
                "L_e": blk["L_e"],
                "raw_gap": blk["raw_gap"],
                "K": "|".join(blk["K"]),
            }
        )
    write_csv(
        os.path.join(out_dir, "or_blocks.csv"),
        ["block_id", "pair", "s", "t", "L_e", "raw_gap", "K"],
        block_rows,
    )

    # 8. Export fixture_specs.csv, points.csv, covariances.csv, link_summaries.csv
    spec_rows = []
    point_rows = []
    cov_rows = []
    link_summary_rows = []
    link_ids = [lk["link_id"] for lk in bf["catalog_links"]]

    for fid, fobj in bf["fixtures"].items():
        sp = fobj["spec"]
        cs = cr["fixture_summaries"][fid]
        spec_rows.append(
            {
                "fixture_id": sp["fixture_id"],
                "route": sp["route"],
                "selection_model": sp["selection_model"],
                "U": sp["U"],
                "selection_pooling": sp["selection_pooling"],
                "sel_intercept": "1" if sp["sel_intercept"] else "0",
                "sel_slopes": "|".join(sp["sel_slopes"]),
                "or_branch": sp["or_branch"],
                "or_intercept": "0" if sp["or_branch"] == "no_intercept" else "1",
                "or_slopes": "|".join(sp["or_slopes"]),
                "trt_cols": "|".join(sp["trt_cols"]),
                "beta_fd_err_1e6": cs["beta_fd_err_1e6"],
                "beta_omit_sel_forcing_err": cs["beta_omit_sel_forcing_err"],
                "core_variance_identity_err": cs["core_variance_identity_err"],
                "cluster_decomposition_identity_err": cs["cluster_decomposition_identity_err"],
                "cluster_vs_unit_omega_diff": cs["cluster_vs_unit_omega_diff"],
                "gmm_target_preservation_err": cs["gmm_target_preservation_err"],
                "gmm_fixed_W_fd_err": cs["gmm_fixed_W_fd_err"],
            }
        )

        # beta_hat
        for idx, bval in enumerate(fobj["beta_hat"]):
            col_name = sp["or_slopes"][idx]
            point_rows.append(
                {
                    "fixture_id": fid,
                    "group": "beta_hat",
                    "name": f"beta_{col_name}",
                    "index": str(idx + 1),
                    "value": bval,
                }
            )

        # alpha_hat
        for idx, (blk_id, aval) in enumerate(fobj["alpha_hat"].items()):
            point_rows.append(
                {
                    "fixture_id": fid,
                    "group": "alpha_hat",
                    "name": blk_id,
                    "index": str(idx + 1),
                    "value": aval,
                }
            )

        # mu_dr, mu_or, mu_ipw
        for grp in ("mu_dr", "mu_or", "mu_ipw"):
            for idx, lid in enumerate(link_ids):
                point_rows.append(
                    {
                        "fixture_id": fid,
                        "group": grp,
                        "name": lid,
                        "index": str(idx + 1),
                        "value": fobj[grp][lid],
                    }
                )

        # theta_gmm_fixed_W
        for idx, tval in enumerate(fobj["theta_gmm_fixed_W"]):
            point_rows.append(
                {
                    "fixture_id": fid,
                    "group": "theta_gmm_fixed_W",
                    "name": f"theta_{idx + 1}",
                    "index": str(idx + 1),
                    "value": tval,
                }
            )

        # V_unit, V_cluster, Omega_cluster
        for mtype in ("V_unit", "V_cluster", "Omega_cluster"):
            mat = fobj[mtype]
            for r_i, row_vals in enumerate(mat):
                for c_i, val in enumerate(row_vals):
                    cov_rows.append(
                        {
                            "fixture_id": fid,
                            "matrix_type": mtype,
                            "row_index": str(r_i + 1),
                            "col_index": str(c_i + 1),
                            "row_link": link_ids[r_i],
                            "col_link": link_ids[c_i],
                            "value": val,
                        }
                    )

        # Per-link summary metrics from check_results.json
        for lid in link_ids:
            linfo = cs["links"][lid]
            row_dict = {"fixture_id": fid, "link_id": lid}
            for k, v in linfo.items():
                row_dict[k] = v
            link_summary_rows.append(row_dict)

    write_csv(
        os.path.join(out_dir, "fixture_specs.csv"),
        [
            "fixture_id",
            "route",
            "selection_model",
            "U",
            "selection_pooling",
            "sel_intercept",
            "sel_slopes",
            "or_branch",
            "or_intercept",
            "or_slopes",
            "trt_cols",
            "beta_fd_err_1e6",
            "beta_omit_sel_forcing_err",
            "core_variance_identity_err",
            "cluster_decomposition_identity_err",
            "cluster_vs_unit_omega_diff",
            "gmm_target_preservation_err",
            "gmm_fixed_W_fd_err",
        ],
        spec_rows,
    )
    write_csv(
        os.path.join(out_dir, "points.csv"),
        ["fixture_id", "group", "name", "index", "value"],
        point_rows,
    )
    write_csv(
        os.path.join(out_dir, "covariances.csv"),
        [
            "fixture_id",
            "matrix_type",
            "row_index",
            "col_index",
            "row_link",
            "col_link",
            "value",
        ],
        cov_rows,
    )

    link_summary_fields = ["fixture_id", "link_id"] + [
        k for k in cr["fixture_summaries"]["F1_baseline_none_intercept_slope"]["links"]["l1"].keys()
    ]
    write_csv(
        os.path.join(out_dir, "link_summaries.csv"),
        link_summary_fields,
        link_summary_rows,
    )

    # 9. Export qualifications.csv (global summary metrics, 4-term witnesses, retained failed predictions)
    qual_rows = []
    for k, v in cr["summary_metrics"].items():
        qual_rows.append(
            {
                "metric_group": "summary_metrics",
                "fixture_id": "ALL",
                "metric_name": k,
                "value": str(v),
            }
        )
    for fid, wdict in cr["four_term_nonzero_witnesses"].items():
        for k, v in wdict.items():
            qual_rows.append(
                {
                    "metric_group": "four_term_nonzero_witnesses",
                    "fixture_id": fid,
                    "metric_name": k,
                    "value": str(v),
                }
            )
    for pred in cr["retained_failed_predictions"]:
        pid = pred["prediction_id"]
        qual_rows.append(
            {
                "metric_group": "retained_failed_predictions",
                "fixture_id": "ALL",
                "metric_name": f"{pid}:original_predicate_status",
                "value": str(pred["original_predicate_status"]),
            }
        )
        qual_rows.append(
            {
                "metric_group": "retained_failed_predictions",
                "fixture_id": "ALL",
                "metric_name": f"{pid}:corrected_identity_status",
                "value": str(pred["corrected_identity_status"]),
            }
        )
        if "actual_max_abs_psi_or" in pred:
            for fid, val in pred["actual_max_abs_psi_or"].items():
                qual_rows.append(
                    {
                        "metric_group": "retained_failed_predictions",
                        "fixture_id": fid,
                        "metric_name": f"{pid}:actual_max_abs_psi_or",
                        "value": str(val),
                    }
                )
        if "actual_max_d_bar_diff" in pred:
            for fid, val in pred["actual_max_d_bar_diff"].items():
                qual_rows.append(
                    {
                        "metric_group": "retained_failed_predictions",
                        "fixture_id": fid,
                        "metric_name": f"{pid}:actual_max_d_bar_diff",
                        "value": str(val),
                    }
                )
        if "actual_or_bias_equal_laws" in pred:
            qual_rows.append(
                {
                    "metric_group": "retained_failed_predictions",
                    "fixture_id": "P3B_PS_ONLY_EQUAL_LAWS_INTERCEPT_CONSISTENCY",
                    "metric_name": f"{pid}:actual_or_bias_equal_laws",
                    "value": str(pred["actual_or_bias_equal_laws"]),
                }
            )
            qual_rows.append(
                {
                    "metric_group": "retained_failed_predictions",
                    "fixture_id": "P3A_PS_ONLY_DISTINCT_LAWS",
                    "metric_name": f"{pid}:actual_or_bias_distinct_laws",
                    "value": str(pred["actual_or_bias_distinct_laws"]),
                }
            )

    write_csv(
        os.path.join(out_dir, "qualifications.csv"),
        ["metric_group", "fixture_id", "metric_name", "value"],
        qual_rows,
    )

    # 10. Write export_receipt.json with all input and output hashes
    exported_files = sorted(
        [
            "atoms.csv",
            "unit_if_rows.csv",
            "nuisance_if_rows.csv",
            "population_examples.csv",
            "calendar.csv",
            "cohorts.csv",
            "catalog_links.csv",
            "or_blocks.csv",
            "fixture_specs.csv",
            "points.csv",
            "covariances.csv",
            "link_summaries.csv",
            "qualifications.csv",
        ]
    )
    output_hashes = {
        fname: sha256_file(os.path.join(out_dir, fname)) for fname in exported_files
    }
    receipt = {
        "status": "EXPORTED_MECHANICALLY_FROM_FROZEN_INDEPENDENT_REFERENCE",
        "source_directory": "validation/reference/dr_finite_v1",
        "output_directory": "tests/testthat/fixtures/dr_finite_v1",
        "expected_values_from_production_backend": False,
        "byte_identical_number_strings": True,
        "n_fixtures": int(bf["fixtures"].__len__()),
        "n_units": int(bf["n_units"]),
        "n_A": int(bf["n_A"]),
        "input_hashes": verified_inputs,
        "output_hashes": output_hashes,
    }
    with open(os.path.join(out_dir, "export_receipt.json"), "w", encoding="utf-8") as f:
        json.dump(receipt, f, indent=2)
        f.write("\n")


if __name__ == "__main__":
    main()
