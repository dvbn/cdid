# Conditional pretest fixtures (CP08 O14/O15)

Verbatim copies of accepted, independently prepared inputs (saved values only;
nothing here is recomputed by the tests). Short names keep archive paths short.

| file | source | sha256 |
|------|--------|--------|
| fA.csv | fixture_A.csv | dd577ccf66ed092b33ee3f2a0c15696242c9a3bd0f2473c23906f6287e125467 |
| fB.csv | fixture_B.csv | 8deb1e41df5a03c9f71ede59160857b097b200929b550c599899b3a4c4a07a57 |
| fBd.csv | fixture_B_dated.csv | 5a242b7e1db0333607a2edd4d41e0d76910db65ce8c7d87e3de259c1e37524bd |
| full.rds | full_registry_golden_reference.rds (A1/B4/B8, full production registry) | 4fa2bd7efe3145054e39b3a67e8a3e49bf59cb0aa5d6da186d78ab814ec6db68 |
| imp.rds | full_registry_imputation_golden_reference.rds (B2/B3) | e14ee3f05f107847c28ff346e981c6b53516cd14d1efd0840323f163b5cdcac6 |
| priv.rds | private_golden_reference.rds (fixed signs, clusters, callback references) | 5f83a7fcf379bd55d3444b490912295198de668638aa77d7ae967154c62ce85f |

The private reference has smaller nuisance catalogs: it is used only for the
fixed sign matrices and cluster labels (private CvM arithmetic) and for callback
rows whose inputs match the production fit (baseline rates per link; separate
pair logits). Public RNG reproducibility is tested with the caller's RNG.
