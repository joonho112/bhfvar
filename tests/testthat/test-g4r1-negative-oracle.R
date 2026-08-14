g4r1_direct_summary <- function(a_probability,b_probability,shares) {
  summarize <- function(p) {
    center <- sum(shares*p)
    sum(shares*(p-center)^2)
  }
  a <- summarize(a_probability)
  b <- summarize(b_probability)
  c(a_between=a,b_between=b,delta=a-b,r=abs(a-b)/max(a,1e-12))
}

make_g4r1_selector_ledger <- function() {
  rows <- list()
  position <- 0L
  for (profile in c("low","high")) {
    for (candidate in seq_len(20L)) {
      position <- position+1L
      regime <- if (candidate<=4L) "grey" else if (candidate<=10L) {
        "detectable"
      } else "equivalence"
      direction <- if (candidate%%2L) "negative" else "positive"
      seed <- bhf_g4r_candidate_seeds(profile,candidate)
      rows[[position]] <- data.frame(
        profile=profile,candidate_index=candidate,
        effect_seed=unname(seed[["effect_seed"]]),
        delta_truth=if(direction=="negative") -.2 else .2,
        r_truth=if(regime=="detectable") .2 else if(
          regime=="equivalence"
        ) .005 else .05,
        regime=regime,direction=direction,stringsAsFactors=FALSE
      )
    }
  }
  do.call(rbind,rows)
}

test_that("negative summary oracle is detected exactly", {
  estimands <- list(
    A=list(summary=c(between=.01)),
    B=list(summary=c(between=.02))
  )
  observed <- bhf_g4r_gap_metrics(estimands)
  expect_equal(observed$delta_truth,-.01,tolerance=0)
  expect_equal(observed$r_truth,1,tolerance=0)
  expect_identical(observed$regime,"detectable")
  expect_identical(observed$direction,"negative")
})

test_that("independent probability microcase has negative detectable gap", {
  shares <- c(A=.5,B=.5)
  oracle <- g4r1_direct_summary(c(A=.4,B=.6),c(A=.2,B=.8),shares)
  expect_equal(oracle[["a_between"]],.01,tolerance=1e-14)
  expect_equal(oracle[["b_between"]],.09,tolerance=1e-14)
  expect_equal(oracle[["delta"]],-.08,tolerance=1e-14)
  expect_equal(oracle[["r"]],8,tolerance=1e-12)
  estimands <- list(
    A=list(summary=c(between=oracle[["a_between"]])),
    B=list(summary=c(between=oracle[["b_between"]]))
  )
  observed <- bhf_g4r_gap_metrics(estimands)
  expect_identical(observed$direction,"negative")
  expect_identical(observed$regime,"detectable")
})

test_that("negative microcase is name/order permutation invariant", {
  shares <- c(A=.5,B=.5)
  first <- g4r1_direct_summary(c(A=.4,B=.6),c(A=.2,B=.8),shares)
  second <- g4r1_direct_summary(
    c(B=.6,A=.4),c(B=.8,A=.2),shares[c("B","A")]
  )
  expect_equal(second,first,tolerance=0)
})

test_that("posterior sign adjudicator handles negative truth and .80 boundary", {
  boundary <- bhf_g4r1_regime_evidence(
    delta_draw=c(-.2,-.1,.1,-.3,-.4),
    a_between_draw=rep(.1,5),delta_truth=-.08,regime="detectable"
  )
  expect_identical(boundary$metric,"correct_sign_probability")
  expect_equal(boundary$probability,.8,tolerance=0)
  expect_true(boundary$pass)
  failed <- bhf_g4r1_regime_evidence(
    delta_draw=c(-.2,-.1,.1,.2,-.4),
    a_between_draw=rep(.1,5),delta_truth=-.08,regime="detectable"
  )
  expect_equal(failed$probability,.6,tolerance=0)
  expect_false(failed$pass)
})

test_that("equivalence evidence uses frozen five-percent ROPE", {
  observed <- bhf_g4r1_regime_evidence(
    delta_draw=c(.001,.002,.006,-.001,.004),
    a_between_draw=rep(.1,5),delta_truth=0,regime="equivalence"
  )
  expect_identical(observed$metric,"practical_equivalence_probability")
  expect_equal(observed$probability,.8,tolerance=0)
  expect_true(observed$pass)
})

test_that("G4-R1 selector takes first six detectable regardless of sign", {
  ledger <- make_g4r1_selector_ledger()
  selected <- bhf_g4r1_select_candidates(ledger[rev(seq_len(nrow(ledger))),])
  expect_equal(nrow(selected),24L)
  counts <- table(selected$profile,selected$selection_group)
  expect_identical(unname(counts[,"detectable"]),c(6L,6L))
  expect_identical(unname(counts[,"equivalence"]),c(6L,6L))
  low <- subset(selected,profile=="low" & selection_group=="detectable")
  expect_identical(low$candidate_index,5:10)
  expect_true(all(c("outcome_seed","mcmc_seed")%in%names(selected)))
})

test_that("G4-R1 selector and evidence fail closed", {
  ledger <- make_g4r1_selector_ledger()
  ledger$posterior <- 0
  expect_error(bhf_g4r1_select_candidates(ledger),
               class="bhf_g4r_preflight_error")
  expect_error(
    bhf_g4r1_regime_evidence(c(-1,1),c(.1),-.1,"detectable"),
    class="bhf_g4r_preflight_error"
  )
  expect_error(
    bhf_g4r1_regime_evidence(c(-1,1),c(.1,.1),0,"detectable"),
    class="bhf_g4r_preflight_error"
  )
})
