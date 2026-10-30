# Private packaged regressions. Every closure comes from actual cdid namespace.
multiplier_ns <- base::asNamespace("cdid")
private_multiplier <- base::get("chained_cluster_multiplier",multiplier_ns,inherits=FALSE)
stopifnot(identical(environment(private_multiplier),multiplier_ns))

test_that("private multiplier resource and shape guards precede data scans", {
  b <- list(max_units=4,max_clusters=2,max_targets=1,max_draws=4,
    max_if_batch_rows=4,max_if_batch_cells=4,max_sign_batch_rows=4,
    max_sign_batch_cells=8,max_workspace_cells=4096,max_workspace_bytes=32768)
  calls <- 0L
  args <- list(target_ids="t",expected_unit_ids=c("u2","u1"),
    unit_cluster_map=list(unit_ids=c("u1","u2"),cluster_ids=c("a","b")),cluster_ids=c("b","a"),
    n_denom=2,center=FALSE,if_batch_fn=function(i){calls<<-calls+1L;list(unit_ids=c("u1","u2"),target_ids="t",if_matrix=matrix(c(1,-1),2,1))},
    n_if_batches=1,sign_batch_fn=function(i)list(draw_ids=letters[1:4],cluster_ids=c("b","a"),sign_matrix=matrix(c(-1,-1,1,1,-1,1,-1,1),4,2)),
    n_sign_batches=1,expected_draw_ids=letters[1:4],n_draws=4,alpha=.25,return_raw_draws=TRUE,budget=b)
  run <- function(change=list()) {a<-args;a[names(change)]<-change;do.call(private_multiplier,a)}
  expect_error(run(list(expected_unit_ids=rep("duplicate",1000))),class="cdid_multiplier_budget_error")
  tight <- b;tight$max_workspace_cells<-1
  expect_error(run(list(expected_unit_ids=c("u1","u1"),budget=tight)),class="cdid_multiplier_budget_error")
  expect_identical(calls,0L)
  expect_error(run(list(n_denom=matrix(2,1,1))),class="cdid_multiplier_input_error")
  expect_error(run(list(center=structure(FALSE,class="opaque"))),class="cdid_multiplier_input_error")
  expect_error(run(list(budget=structure(b,class="opaque"))),class="cdid_multiplier_budget_error")
  expect_error(run(list(unit_cluster_map=structure(args$unit_cluster_map,class="opaque"))),class="cdid_multiplier_unit_alignment_error")
  expect_error(run(list(if_batch_fn=function(i) c(args$if_batch_fn(i),list(if_matrix=matrix(0,2,1))))),class="cdid_multiplier_if_dimension_error")
  expect_error(run(list(sign_batch_fn=function(i)c(args$sign_batch_fn(i),list(draw_ids=letters[1:4])))),class="cdid_multiplier_sign_dimension_error")
  expect_error(run(list(if_batch_fn=function(i)list(unit_ids=letters[1:3],target_ids="t",if_matrix=matrix(NA_real_,3,1)))),class="cdid_multiplier_unit_alignment_error")
  expect_error(run(list(sign_batch_fn=function(i)list(draw_ids=letters[1:5],cluster_ids=c("b","a"),sign_matrix=matrix(0,5,2)))),class="cdid_multiplier_budget_error")
})

test_that("loop/vector lengths and payload integer limits are distinct from matrix dimensions", {
  scalar <- get("validate_positive_int_scalar",multiplier_ns,inherits=FALSE)
  estimate <- get("compute_multiplier_workspace_estimate",multiplier_ns,inherits=FALSE)
  expect_identical(scalar(2147483648,"n_if_batches"),2147483648)
  expect_error(scalar(2147483648,"n_draws",matrix_dimension=TRUE),class="cdid_multiplier_budget_error")
  expect_error(estimate(1e15,1,1,1,1,1),class="cdid_multiplier_budget_error")
  expected <- 2*4*2+6*2^2+2*16*2+6*6+4*4+8*16+32*2+2*6*2+4*16*4
  expect_identical(expected,628)
  expect_identical(estimate(6,4,2,16,64,128)$total_persistent_cells,expected)
})
