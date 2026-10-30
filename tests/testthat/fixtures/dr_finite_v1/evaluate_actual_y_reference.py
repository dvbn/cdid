"""Eight fixed evaluations of the accepted functional, not its main/FD campaign."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import time

import numpy as np

base = Path(__file__).resolve().parent
source = base.parents[3] / "validation/reference/dr_finite_v1/generate_dr_reference.py"
atoms_path = base.parents[3] / "validation/reference/dr_finite_v1/atoms.csv"
expected = {
    source: "7cb9f66c9ed6d250c734796508c18eedee6d8aca328614ed6751e20190c98282",
    atoms_path: "f6869b76895e3c0c91b8d7152bf38f6dbff65ab5400f033c07730c378a03d8d1",
}
for path, value in expected.items():
    assert hashlib.sha256(path.read_bytes()).hexdigest() == value
assert np.finfo(np.longdouble).nmant == 112
spec = importlib.util.spec_from_file_location("accepted_dr_reference", source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
integer_fields = {"unit_id", "cluster_id", "cohort", "F", "B", "A"}
atoms = []
for row in csv.DictReader(atoms_path.open()):
    atoms.append({k: int(v) if k in integer_fields or k.startswith(("R_", "S_")) else float(v)
                  for k, v in row.items()})
changes = []
for row in atoms:
    for pair, s, t in [("e12", 2000, 2002), ("e23", 2002, 2005), ("e13", 2000, 2005)]:
        if row["S_" + pair] == 1:
            # Compute the same binary64 subtraction as actual masked-Y owner input.
            actual_delta = row[f"Y_{t}"] - row[f"Y_{s}"]
            changes.append(abs(actual_delta - row["Delta_" + pair]))
            row["Delta_" + pair] = actual_delta
        else:
            row["Delta_" + pair] = 0.0  # Structural indicator, never an observed value.
assert len(atoms) == 36 and sum(x["A"] for x in atoms) == 24
assert [x["unit_id"] for x in atoms] == list(range(1, 37))
out = base / "actual_y_reference"
out.mkdir(exist_ok=False)
records = []
def add(fid, group, key1, key2, value):
    records.append([fid, group, str(key1), str(key2), repr(float(value))])
started = time.perf_counter()
for fixture in module.FIXTURE_SPECS:
    result = module.run_dr_pipeline(atoms, fixture, dtype=np.longdouble)
    fid = fixture["fixture_id"]
    for name, link in result["links"].items():
        for group in ["mu_dr", "mu_or", "mu_ipw"]:
            add(fid, group, name, "", link[group])
        for group in ["psi_dr", "psi_dir", "psi_or", "psi_ps", "psi_sel", "psi_ipw"]:
            for i, value in enumerate(link[group]):
                add(fid, group, atoms[i]["unit_id"], name, value)
    for j, value in enumerate(result["beta_hat"]):
        add(fid, "beta_hat", "beta_" + fixture["or_slopes"][j], "", value)
    for key, value in result["alpha_hat"].items():
        add(fid, "alpha_hat", key, "", value)
    for group, matrix, columns in [
        ("beta_if", result["xi_beta"], ["beta_" + x for x in fixture["or_slopes"]]),
        ("tau_if", np.column_stack([result["xi_tau"][x] for x in module.COHORTS]), ["tau_" + str(x) for x in module.COHORTS]),
        ("eta_if", np.column_stack([result["xi_eta"][x] for x in module.COHORTS]), ["eta_" + str(x) for x in module.COHORTS]),
    ]:
        for i, row in enumerate(matrix):
            for j, value in enumerate(row):
                add(fid, group, atoms[i]["unit_id"], columns[j], value)
    names = [link["link_id"] for link in module.CATALOG_LINKS]
    phi = np.column_stack([result["links"][name]["psi_dr"] for name in names])
    cluster = np.stack([phi[[i for i, row in enumerate(atoms) if row["cluster_id"] == key]].sum(axis=0)
                        for key in sorted({row["cluster_id"] for row in atoms})])
    for group, matrix in [("V_unit", phi.T @ phi / 36**2),
                          ("V_cluster", cluster.T @ cluster / 36**2),
                          ("Omega_cluster", cluster.T @ cluster / 36)]:
        for i, row in enumerate(matrix):
            for j, value in enumerate(row):
                add(fid, group, names[i], names[j], value)
elapsed = time.perf_counter() - started
with (out / "values.csv").open("w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["fixture_id", "group", "key1", "key2", "value"])
    writer.writerows(records)
receipt = {
    "scope": "eight prespecified accepted fixtures evaluated at actual binary64 observed-level differences",
    "source_generator_sha256": expected[source], "source_atoms_sha256": expected[atoms_path],
    "accepted_inputs_mutated": False, "original_Delta_reference_replaced": False,
    "main_generator_calls": 0, "finite_difference_calls": 0, "RNG_draws": 0,
    "production_imports_or_fits": 0, "functional_evaluations": len(module.FIXTURE_SPECS),
    "compute_seconds": elapsed, "numpy": np.__version__,
    "dtype": "binary128 with binary64 atom inputs and binary64 observed-level subtraction",
    "longdouble_epsilon": str(np.finfo(np.longdouble).eps),
    "maximum_observed_input_difference": max(changes), "rows": len(records),
    "values_sha256": hashlib.sha256((out / "values.csv").read_bytes()).hexdigest(),
}
(out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps(receipt, indent=2))
