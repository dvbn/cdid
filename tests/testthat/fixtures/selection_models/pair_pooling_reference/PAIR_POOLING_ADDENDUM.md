# Independent reference for explicit pooling of direct-pair logits

The flexible **separate** specification is the default for both observation
models. Date logits retain independent date/cohort
coefficients; pair logits retain independent actual-pair/cohort coefficients.
`within_cohort` and `across_cohorts` sharing remain explicit alternatives, with
the appropriate date/cohort or pair/cohort intercepts. There is no
automatic pooling fallback. This reference supplies fixed-input expectations
for the packaged tests.

## Same block algebra with a different response

Use the accepted U in {F,B,A}, fixed law `P(f)=H(wf)/H(w)`, actual pair
`e=(s,t)`, response `S_e=R_s R_t`, and unique cells `(e,k)`. Fix the union of
the scientific fitting catalogs before reporting selection. Each cell enters
once, regardless of how many links consume it. Group cells into parameter
blocks m: one cell per block for separate fits; one cohort's pairs per block
for within-cohort sharing; or all selected pair/cohort cells for across-cohort
sharing. With intercept-containing formulas, write

`p_e,k^U(X)=expit(alpha_e,k+z(X)' beta_m)`.

For each original F unit, the stacked block score and negative derivative are

`h_i,m=sum_(e,k in m) w_i I_ik U_i (S_ie-p_e,k^U(X_i)) Z_e,k,i`,

`J_m=H[sum_(e,k in m) w I_k U p_e,k^U(1-p_e,k^U) Z_e,k Z_e,k']`.

Thus `xi_theta_m,i=solve(J_m,h_i,m)`. The date-model block elimination applies
verbatim with cell `(date,cohort)` replaced by `(actual pair,cohort)`, response
R replaced by S, and cell weights multiplied by U. Its information-weighted
means give the same profiled slope/intercept formulas. The fully separate
cell has exactly the accepted CP04 score, information and IF. Ordinary formula
intercept semantics and homogeneous predictor-span reductions carry over;
pooling does not insert a removed intercept.

Retain `a_k^U=H(w I_k U A)/H(w I_k U)` and
`xi_a_k^U,i=w_i I_ik U_i(A_i-a_k^U)/H(w I_k U)`. For `U=A`, the rate is one
and its row is zero. The known population-normalized transport rho is fixed,
with zero derivative. On an actual pair, the existing observation weight is
`v_e,k=a_k^U/(rho_k^U p_e,k^U)`, and its complete log-weight row is

`xi_a_k^U/a_k^U - (1-p_e,k^U(X)) Z_e,k(X)' xi_theta_m`.

Substitute this loading into the accepted IPW/shared-OR stack; no downstream
rederivation is needed. Sum all reused pair scores **by original unit** before
outer products or cluster sums. Shared slopes transmit other pairs' scores
into a pair prediction; across-cohort slopes also transmit other cohorts'
scores. Even separate coefficient fits can have correlated original-unit
influence rows because pairs overlap. Preserve all F identities in normalization
and inference, including off-A nuisance rows when U permits them.

No date-independence assumption is added to the direct-pair model. Correctness
of the selected pooled conditional probability on its actual U is required;
correct pooled F logits do not generally imply correct pooled A/B logits.
The accepted actual-pair outcome-mean restriction, transports and estimator-
specific nuisance/identification premises remain. Pooling is a user-specified
restriction on coefficients, not a consequence of missingness terminology.

## Frozen finite evidence

`reference.py` uses the existing frozen NumPy dense logit solver in memory with
S as its response. It enumerates 32 records: two cohorts, two X values and all
eight three-date histories. Fixed positive weights define an explicit weighted
law P under which direct-pair logits have common slope 0.3 and cohort intercepts
-2.2/-1.9. For each (k,X), let p be that pair probability, assign mass p/2 to
history 111 and each of 110,101,011, singleton masses `(0.2,0.15,0.15)*p`, and
mass `1-2.5p` to 000. Pair marginals are p, A mass is 2p and B mass is 2.5p.
Consequently conditional pair probabilities are 0.4 on B and 0.5 on A; all
nine selected mode/U models are correct for this constructed law. Known rho
is normalized once under P and frozen during all case perturbations.

H is uniform over the atoms, so the fixed weights encode this law and depend
on history. This is deterministic functional evidence, not sampled DGP,
survey-design or coverage evidence. Conditional date independence is not
assumed for this direct-pair table. Its point weights agree across the nine
cases by construction; this is not general finite-function/IF equivalence.

All coefficients, empirical a, pair predictions and combined v rows agree with
original-unit refits at steps 1e-4 and 5e-5, maximum **9.21869e-8 < 3e-7**.
Separate fits also match independent binary-cell log odds and the direct-pair block
score. Dropping other-pair score effects fails both shared modes; its effect is
correctly zero in separate fits. Cross-cohort prediction effects are nonzero
only with across-cohort sharing. Independent-pair-copy covariance fails on F/B.
In this symmetric A law the S_e have pairwise covariance zero, so that omission
is truthfully zero there. The initial contrary hand-prediction and its failure
are retained; no fixed atom, formula or tolerance was changed to obtain it.

`separate_date_reference.py` adds only the missing separate-date F branch,
using the same 32 original atoms and frozen solver with explicitly distinct
date/cohort slopes. Full coefficient/a/product/weight derivatives pass at both
steps, maximum **9.39284e-8 < 3e-7**. Unit covariance retains all date scores;
the independent-date-copy covariance mutant fails. A nonendpoint date has zero
predictor loading when coefficients are fully separate, while reused endpoints
still couple pair predictions. The original three-case date reference is
unchanged. Centering/root tolerance remains 2e-11.

`results.json`, `separate_date_results.json`, atom tables, IF CSVs and
`receipt.json` retain the fixed-input evidence and historical input hashes.
No random draws, Monte Carlo, package execution, installation or compilation
were used to construct this reference.
