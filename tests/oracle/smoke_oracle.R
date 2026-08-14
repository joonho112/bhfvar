before_namespaces <- loadedNamespaces()
sys.source("tests/oracle/estimand_oracle.R", envir = environment())

L <- log(3)
design <- list(
  row_id = paste0("r", 1:8),
  domain_labels = c("A", "B"),
  stratum_labels = c("H1", "H2"),
  psu_labels = c("J1", "J2", "J3", "J4"),
  state_id = c(rep(1L, 4), rep(2L, 4)),
  stratum_id = rep(c(1L, 1L, 2L, 2L), 2),
  psu_flat_id = rep(1:4, 2),
  w_lik = c(2.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 2.5),
  population_share = c(A = 0.5, B = 0.5)
)
draws <- list(
  draw_id = "d1",
  alpha = 0,
  u_state = matrix(c(-L, L), 1, dimnames = list("d1", c("A", "B"))),
  u_stratum = matrix(c(L / 2, -L / 2), 1,
                     dimnames = list("d1", c("H1", "H2"))),
  u_psu = matrix(c(L / 2, -L / 2, L / 2, -L / 2), 1,
                 dimnames = list("d1", c("J1", "J2", "J3", "J4"))),
  sigma_state = L,
  sigma_stratum = L / 2,
  sigma_psu = L / 2
)
result <- bhf_reference_oracle(
  design,
  draws,
  vhat = c(B = 3 / 32, A = 1 / 32)
)

stopifnot(
  isTRUE(all.equal(unname(result$A$p_state[1, ]), c(1 / 4, 3 / 4),
                   tolerance = 1e-12)),
  isTRUE(all.equal(unname(result$B$p_state[1, ]), c(31 / 80, 49 / 80),
                   tolerance = 1e-12)),
  isTRUE(all.equal(unname(result$B$within_bernoulli), 343 / 1600,
                   tolerance = 1e-12)),
  isTRUE(all.equal(unname(result$B$within_mixture), 147 / 6400,
                   tolerance = 1e-12)),
  isTRUE(all.equal(result$A_star$summary[1, "between"], 0,
                   tolerance = 1e-12)),
  identical(sort(before_namespaces), sort(loadedNamespaces())),
  !"bhfvar" %in% loadedNamespaces()
)

writeLines("clean-base-oracle-gold: PASS")
