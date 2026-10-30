# Deterministic two-date illustration, not simulation or coverage evidence.
# Each cohort has equal population X cells; pair probabilities are .25/.75.
# Pair trends and the heterogeneous treatment effect X are functions of X,
# so pair mean-MAR and untreated parallel trends hold by construction.
cp04_demo_fixture <- function() {
  frame <- expand.grid(replicate=1:24, x=c(0,1), G=c(2,0))
  frame$id <- seq_len(nrow(frame))
  pair_n <- ifelse(frame$x==0,6L,18L)
  frame$n_obs <- ifelse(frame$replicate<=pair_n,2L,
    ifelse(frame$replicate<=pair_n+(24L-pair_n)/2,1L,0L))
  panel <- data.frame(id=rep(frame$id,each=2L),period=rep(1:2,nrow(frame)),
    G=rep(frame$G,each=2L),x=rep(frame$x,each=2L))
  panel$y <- panel$id/100 + (panel$period-1)*(panel$x/3 + (panel$G==2)*panel$x)
  count <- frame$n_obs[match(panel$id,frame$id)]
  panel$y[count==0L | (count==1L & panel$period==2L)] <- NA_real_
  # Linked zero-outcome units are supplied through frame, rather than fake
  # outcome rows. Missing rows and observed outcome counts are distinct.
  panel <- panel[count>0L,]
  list(panel=panel,frame=frame[,c("id","G","x","n_obs")],population_att=.5)
}
