# Independent date-logit finite reference

This reference was derived and frozen independently of production date-model
outputs and code. It uses NumPy only, no package code, RNG, sampled
DGPs, Monte Carlo, build, compiler or downstream outcome fits. It supplies independent fixed-input
expectations for the packaged tests.

## Fixed records and weighted law

There are 32 distinct original records: two cohorts, X in {-1,1}, and all eight
three-date observation histories. All records, including no-observation units
and singletons, belong to F. A indicates at least two observed dates.

H is uniform on these 32 atoms. Each specification fixes positive weights
`w_i=32 q(G_i,X_i) product_t p_it^R_it (1-p_it)^(1-R_it)` with joint (G,X)
masses `(0.21,0.29;0.31,0.19)`. Thus the *declared weighted law*
`P(f)=H(wf)/H(w)` has conditionally mutually independent dates by construction.
The weights depend on observation history under H. This is an intentional
finite weighted-functional witness, not an unqualified sampling or coverage
claim. Independence under P does not remove date-score covariance under H when
the complete score contains w. Under case perturbations w is frozen, all roots
are refitted, and the perturbed law need not remain in the model family.

The three small cases have date/cohort intercepts and: (1) two cohort slopes
shared over dates, (2) one slope shared over dates and cohorts, or (3) no X slope.
Their design dimensions are 8, 7 and 6. These are reference labels, not new API
spellings or a default pooling choice. Known coefficients generate the weighted
atom law, while a separate dense Newton solve starting at zero recovers each
coefficient vector. No production fitter or pooled-WLS implementation is used.

## Original-unit equations and derivatives

Write `H_N f=sum_F f_i/N`, `p_it=expit(z_it' theta)`. The zero equation is

`H_N h_theta=0`, with `h_theta,i=w_i sum_t z_it (R_it-p_it)`.

Its negative derivative is
`J=H_N[w sum_t p_it(1-p_it) z_it z_it']`, so the original-unit row is
`xi_theta,i=solve(J,h_theta,i)`. All dates of an original unit move together.
The information contains one w, while the score outer product contains w^2.
The normalized case functional is `H_N^c f=sum_i c_i f_i/sum_i c_i`.
The checked identity is `N d theta(c)/d c_i=xi_theta,i`.

For a fixed (G,X) prediction point and pair (s,t), the F probability is
`pi_st=p_s p_t`. Its gradient is
`pi_st[(1-p_s)z_s+(1-p_t)z_t]`, giving
`xi_pi_st,i=gradient' xi_theta,i`. The same fitted theta and original-unit row
are reused for all three pairs. An additional date can affect a pair through
shared slopes even when it is not an endpoint; that channel correctly vanishes
in the date-only case. Pair product validity requires the appropriate conditional
pair independence under P; this reference explicitly constructs the stronger
mutual independence and does not infer it in other data.

Keep the accepted empirical repeat rate, rather than substituting a model-based
probability: `a_k=H_N(w I_k A)/H_N(w I_k)` and
`xi_a_k,i=w_i I_ki(A_i-a_k)/H_N(w I_k)`. For the F correction
`b_st(k,x)=a_k/pi_st(k,x)`, its row is
`xi_b=xi_a/pi_st-a_k xi_pi_st/pi_st^2`. All these coordinates are jointly checked.
Both empirical repeat and date-score channels can be nonzero on A=0 records.
Only F probabilities are multiplied. No U-, B- or A-conditional probabilities
are multiplied and then presented as population pair probabilities.

## Checks and limits

`python3 reference.py` (from this fixture directory) passes all
32 original-unit perturbations at each of 1e-4 and 5e-5 in all three cases.
The complete coefficient/rate/pair/multiplier maximum error is **5.77299e-8**,
below unchanged 3e-7; root and centering checks use 2e-11. `results.json` retains
the actual errors, coefficients, information matrices, covariance and witnesses;
`atoms.csv` and the three IF CSVs retain the exact fixed data and complete rows.

Every case rejects a dropped endpoint derivative, a one-date copy score used as
the original-unit score, and covariance formed by summing separate date outer
products. The latter loses within-unit cross terms; its coefficient-covariance
error is at least 0.02936. Shared-pair covariance is nonzero. Removing a
nonendpoint date score fails the slope-sharing cases and correctly gives zero
in the date-only case. Original A=0 coefficient, rate and multiplier rows are
explicitly nonzero and retained.

With three independent Bernoulli(1/2) dates on F, conditioning on A leaves four
histories. Each date marginal is 3/4, the joint pair probability is 1/2, and the
product of those conditional marginals is 9/16. This exact witness separates
population date independence from independence after selection on repeat status.

This does not verify outcome ignorability, causal identification, nuisance
correctness for a future dataset, asymptotics, inference coverage, arbitrary
fitting universes, estimated weights/transports or the eventual public backend.
`receipt.json` retains the historical input hashes and pins the standalone
outputs. The fixed-weight, full-F normalization is defined above.
