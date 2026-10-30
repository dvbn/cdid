test_that("owner cohort storage and preflight respect G exceeding calendar T", {
  ns <- asNamespace("cdid")
  ctrl <- function(ret=1e6,mat=1e6,work=1e6) ns$ipw_owner_control(4,ret,mat,work,1e6)
  make <- function(n) {
    d <- data.frame(id=rep(seq_len(n),each=2),time=rep(c(0,1),n),
                    y=rep(c(1,2),n),g=rep(2012+seq_len(n)/10,each=2))
    p <- ns$construct_ipw_panel_owner(d,"id","time","y","g",ctrl())
    f <- ns$construct_ipw_frame_owner(p,seq_len(n),matrix(1,n,1),control=ctrl())
    v <- ns$construct_ipw_link_view(f,0,1,g=d$g[1],K=d$g[3],control=ctrl())
    list(panel=p,frame=f,view=v)
  }
  # Preserve the initially mistaken 512@1m success prediction as the correct
  # work-bound refusal, with the exact failed script/output under audit/.
  expect_error(make(512L),class="cdid_ipw_owner_bound_error")
  a <- make(128L); b <- make(64L)
  expect_identical(length(a$panel$a_cohorts),128L)
  expect_identical(a$panel$n_calendar_dates,2L)
  expect_identical(length(a$view$cohorts$owner_cohort_index),2L)
  expect_identical(length(b$view$cohorts$owner_cohort_index),2L)
  expect_null(a$view$cohorts$a_cohort_to_view)
  expect_identical(a$view$preflight$total_retained_entries,b$view$preflight$total_retained_entries)
  expect_lte(a$frame$accounting$retained_vector_entries+a$frame$accounting$retained_matrix_entries,
             a$frame$preflight$total_retained_entries)
  # Exact original failing cap: preflight must now refuse before poisoned Z.
  expect_error(ns$construct_ipw_frame_owner(a$panel,1:128,matrix(NaN,128,1),control=ctrl(422)),
               class="cdid_ipw_owner_bound_error")
  expect_error(ns$construct_ipw_link_view(a$frame,0,1,g=2012.1,K=2012.2,control=ctrl(49)),
               class="cdid_ipw_owner_bound_error")
  expect_error(ns$query_ipw_panel_cell_status(a$panel,1,0,control=ctrl(mat=1)),
               class="cdid_ipw_owner_bound_error")
  ch <- ns$query_ipw_view_chunk(a$view,1,4,ctrl(),read_outcomes=TRUE)
  expect_identical(ch$in_risk,c(TRUE,TRUE,FALSE,FALSE))
  expect_identical(ch$observed_risk$a_index,1:2)
  # Copying an equal array creates a foreign outer object; do not recursively
  # compare its potentially huge content to accept it as the trusted snapshot.
  copy <- a$panel; copy$unit_id_key <- as.double(a$panel$unit_id_key)+0
  expect_false(ns$gmm_sexp_aliases(copy,a$panel))
  expect_error(ns$validate_ipw_panel_owner(copy),class="cdid_ipw_owner_integrity_error")
  alias <- a$panel
  expect_silent(ns$validate_ipw_panel_owner(alias))
  # Refuse fractional/Inf dense cohort indices rather than truncating to valid
  # integers or converting Inf to the same NA used outside risk.
  ci <- rep(NA_real_,128);ci[1:2]<-c(1.5,2)
  expect_error(ns$verify_ipw_view_dense_masks(a$view,ctrl(),cohort_index=ci),
               class="cdid_ipw_owner_dense_mask_error")
  ci[1]<-1;ci[3]<-Inf
  expect_error(ns$verify_ipw_view_dense_masks(a$view,ctrl(),cohort_index=ci),
               class="cdid_ipw_owner_dense_mask_error")
})

test_that("long character comparison temporaries are preflighted in bytes", {
  ns <- asNamespace("cdid")
  keys <- c("a",paste0("m",strrep("x",16384)),"z","zz")
  large <- ns$ipw_owner_control(4,8e6,8e6,1e8,8e6)
  tiny <- ns$ipw_owner_control(1,1e6,1000,10000,1e6)
  d <- data.frame(id=rep(keys,each=2),time=rep(0:1,4),y=as.double(1:8),g=rep(c(2,2,Inf,Inf),each=2))
  p <- ns$construct_ipw_panel_owner(d,"id","time","y","g",large)
  expect_identical(p$id_byte_metadata$max,16385)
  expect_identical(p$id_byte_metadata$total,32778)
  expect_gt(p$preflight$id_byte_preflight$additional_scratch_entries,0)
  expect_error(ns$query_ipw_panel_cell_status(p,"a",0,tiny),class="cdid_ipw_owner_bound_error")
  scratch_tiny <- ns$ipw_owner_control(1,8e6,1000,1e8,8e6)
  expect_error(ns$query_ipw_panel_cell_status(p,"a",0,scratch_tiny),class="cdid_ipw_owner_bound_error")
  work_tiny <- ns$ipw_owner_control(1,8e6,8e6,10000,8e6)
  expect_error(ns$query_ipw_panel_cell_status(p,"a",0,work_tiny),class="cdid_ipw_owner_bound_error")
  expect_identical(ns$query_ipw_panel_cell_status(p,"a",0,large)$status,"observed")
  # A foreign query ID also needs its own byte-length term, not merely the
  # stored-key bound; output remains missing-row under the explicit large cap.
  unknown <- paste0("y",strrep("x",32768))
  expect_error(ns$query_ipw_panel_cell_status(p,unknown,0,tiny),class="cdid_ipw_owner_bound_error")
  expect_identical(ns$query_ipw_panel_cell_status(p,unknown,0,large)$status,"missing_row")
  f <- ns$construct_ipw_frame_owner(p,keys,matrix(1,4,1),control=large)
  expect_gt(f$preflight$id_byte_preflight$additional_scratch_entries,0)
  v <- ns$construct_ipw_link_view(f,0,1,g=2,K=Inf,control=large)
  tr <- new.env(parent=emptyenv());tr$calls <- 7L
  expected_tiny <- ns$ipw_owner_control(4,1e6,10000,1e7,1e6)
  expect_error(ns$query_ipw_view_chunk(v,1,4,expected_tiny,expected_a_id=keys,trace_env=tr),class="cdid_ipw_owner_bound_error")
  expect_identical(tr$calls,7L) # byte preflight precedes trace writes too
  ch <- ns$query_ipw_view_chunk(v,1,4,large,expected_a_id=keys,read_outcomes=TRUE)
  expect_identical(ch$id,keys)
  expect_gt(ch$preflight$id_byte_preflight$additional_scratch_entries,0)
  # Ordering code itself is unchanged and normalizes compared scalars to UTF8.
  latin <- rawToChar(as.raw(233));Encoding(latin) <- "latin1"
  utf <- enc2utf8(latin)
  expect_false(ns$ipw_owner_chr_lt_radix(latin,utf))
  expect_false(ns$ipw_owner_chr_lt_radix(utf,latin))
})

test_that("sorted lookups preserve mixed encodings and long ALTREP positions", {
  ns <- asNamespace("cdid")
  x <- seq_len(2^31+5)
  before <- capture.output(.Internal(inspect(x,0L,1L)))[1L]
  expect_true(ns$gmm_sexp_aliases(x,x))
  expect_identical(ns$ipw_owner_match_sorted_double(2^31+3,x),2^31+3)
  after <- capture.output(.Internal(inspect(x,0L,1L)))[1L]
  expect_match(before,"compact")
  expect_match(after,"compact")
  # Byte comparison must normalize each compared scalar, as radix character
  # ordering treats equal Latin1/UTF8 identifiers as the same character value.
  latin <- rawToChar(as.raw(233));Encoding(latin)<-"latin1"
  utf <- enc2utf8(latin)
  keys <- sort(c("z",utf),method="radix")
  expect_identical(ns$ipw_owner_match_sorted_chr(latin,keys),match(latin,keys))
  expect_identical(ns$ipw_owner_match_sorted_double(c(1,2,4),c(0,1,2,3,4),2,4),c(2L,3L,NA_integer_))
})

test_that("diagnostic traces retain only bounded latest-call positions", {
  ns <- asNamespace("cdid")
  ctrl <- ns$ipw_owner_control(4,1e6,1e6,1e6,1e6)
  d <- data.frame(id=rep(1:4,each=2),time=rep(0:1,4),y=1:8,g=rep(c(2,2,Inf,Inf),each=2))
  p <- ns$construct_ipw_panel_owner(d,"id","time","y","g",ctrl)
  f <- ns$construct_ipw_frame_owner(p,1:4,matrix(1,4,1),control=ctrl)
  v <- ns$construct_ipw_link_view(f,0,1,g=2,K=Inf,control=ctrl)
  tr <- new.env(parent=emptyenv())
  ns$query_ipw_view_chunk(v,1,4,ctrl,read_outcomes=TRUE,trace_env=tr)
  expect_identical(tr$calls,1L)
  expect_identical(tr$pos_s_index,c(1L,3L,5L,7L))
  ns$query_ipw_view_chunk(v,3,3,ctrl,read_outcomes=TRUE,trace_env=tr)
  expect_identical(tr$calls,1L)
  expect_identical(tr$pos_s_index,5L)
  expect_null(tr$pos_s_history)
  tiny <- ns$ipw_owner_control(4,1e6,1,1e6,1e6)
  expect_error(ns$query_ipw_view_chunk(v,1,4,tiny,read_outcomes=TRUE,trace_env=tr),
               class="cdid_ipw_owner_bound_error")
  expect_identical(tr$pos_s_index,5L) # preflight refusal leaves trace unchanged
  ns$query_ipw_view_chunk(v,1,1,ctrl,read_outcomes=FALSE,trace_env=tr)
  expect_identical(tr$calls,0L)
  expect_identical(tr$pos_s_index,integer(0L))
  expect_identical(tr$pos_t_index,integer(0L))
  nu <- ns$construct_ipw_nuisance_view(f,2,Inf,0,1,ctrl)
  unchecked <- ns$construct_ipw_link_view(f,0,1,nuisance_view=nu,control=ctrl,verify_nuisance_mask=FALSE)
  expect_false(unchecked$nuisance_mask_verified)
  expect_error(ns$validate_ipw_link_view(unchecked,nuisance_view=nu),
               class="cdid_ipw_owner_mask_unverified_error")
  expect_true(ns$verify_ipw_view_mask_compatibility(nu,unchecked,ctrl)$compatible)
  checked <- ns$construct_ipw_link_view(f,0,1,nuisance_view=nu,control=ctrl)
  expect_true(checked$nuisance_mask_verified)
  expect_silent(ns$validate_ipw_link_view(checked,nuisance_view=nu))
})

test_that("optional candidate normalization and duplicate matrix fields are budgeted", {
  ns <- asNamespace("cdid")
  ctrl <- ns$ipw_owner_control(4,1e6,1e6,1e6,1e6)
  pre <- ns$preflight_ipw_link_view(2,64,64,ctrl,has_candidates=TRUE)
  expect_identical(pre$retained_matrix_entries,3*64^2) # A_beta_r, L, H_chol
  expect_gte(pre$constructor_scratch_entries,4*pre$retained_matrix_entries)
  tiny <- ns$ipw_owner_control(4,1e6,50000,1e6,1e6)
  expect_error(ns$preflight_ipw_link_view(2,64,64,tiny,has_candidates=TRUE),class="cdid_ipw_owner_bound_error")
  nu <- ns$preflight_ipw_nuisance_view(2,64,2,ctrl,has_candidates=TRUE)
  expect_gte(nu$constructor_scratch_entries,4*64^2)
  empty <- ns$construct_ipw_panel_owner(data.frame(id=1,time=0,y=NA_real_,g=2),"id","time","y","g",ctrl)
  tiny_empty <- ns$ipw_owner_control(1,1,1,1,1)
  expect_error(ns$query_ipw_panel_pair_positions(empty,0,1,control=tiny_empty),class="cdid_ipw_owner_bound_error")
})
