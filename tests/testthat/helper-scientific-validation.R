.bhf_science_seed <- 7101005L
.bhf_science_tolerance <- c(absolute = 1e-10, relative = 1e-8)
.bhf_science_frozen_hashes <- c(
  "estimand_oracle.R" =
    "03ff2fb87e1ca1124d671607774394b31da09a9254b79da8fd4cedf035897d72",
  "smoke_oracle.R" =
    "80c81b621ec0e780fa9eb7d4a17bee2eac71d8f5d7b4dc564e07e5a70d652be8",
  "manifest.txt" =
    "885a44a0e2fa92ca028b15301c4c17202031afbcd55504c206bd9058157dbb8a",
  "test-estimand-oracle.R" =
    "ab8b54f6c20c819053aa81a8b815938051c93b0115b36a1ba376a002454d3dd9"
)

.bhf_science_oracle <- local({
  oracle <- new.env(parent = baseenv())
  sys.source(
    testthat::test_path("..", "oracle", "estimand_oracle.R"),
    envir = oracle
  )
  oracle$bhf_reference_oracle
})

.bhf_science_seeded <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.bhf_science_close <- function(actual, expected,
                               tolerance = .bhf_science_tolerance) {
  actual <- as.numeric(actual)
  expected <- as.numeric(expected)
  if (length(actual) != length(expected)) return(FALSE)
  finite_pair <- is.finite(actual) & is.finite(expected)
  same_missing <- is.na(actual) & is.na(expected)
  bound <- tolerance[["absolute"]] +
    tolerance[["relative"]] * pmax(abs(actual), abs(expected))
  all(same_missing | (finite_pair & abs(actual - expected) <= bound))
}

.bhf_science_expect_close <- function(actual, expected, info,
                                      tolerance = .bhf_science_tolerance) {
  expect_true(
    .bhf_science_close(actual, expected, tolerance),
    info = paste0(
      info, "; tolerance abs=", tolerance[["absolute"]],
      ", rel=", tolerance[["relative"]]
    )
  )
}

.bhf_science_random_fixture <- function(seed = .bhf_science_seed,
                                        draws = 32L) {
  .bhf_science_seeded(seed, {
    S <- 4L
    H <- 3L
    J_h <- c(2L, 3L, 1L)
    J <- sum(J_h)
    psu_stratum <- rep(seq_len(H), J_h)
    grid <- expand.grid(
      state = seq_len(S),
      psu = seq_len(J),
      KEEP.OUT.ATTRS = FALSE
    )
    state_id <- as.integer(grid$state)
    psu_flat_id <- as.integer(grid$psu)
    stratum_id <- as.integer(psu_stratum[psu_flat_id])
    N <- nrow(grid)

    state_labels <- paste0("S", seq_len(S))
    stratum_labels <- paste0("H", seq_len(H))
    psu_labels <- paste0("J", seq_len(J))
    draw_id <- paste0("d", seq_len(draws))

    share <- stats::setNames(runif(S, 0.2, 1.5), state_labels)
    share <- share / sum(share)
    w_lik <- runif(N, 0.2, 2.5)
    w_lik <- w_lik * N / sum(w_lik)

    z_state <- matrix(rnorm(draws * S), draws, S)
    u_state <- sweep(
      z_state,
      1L,
      as.vector(z_state %*% unname(share)),
      "-"
    )
    z_stratum <- matrix(rnorm(draws * H), draws, H)
    u_stratum <- sweep(z_stratum, 1L, rowMeans(z_stratum), "-")
    z_psu <- matrix(rnorm(draws * J), draws, J)
    u_psu <- z_psu
    start <- c(1L, cumsum(J_h)[-H] + 1L)
    for (h in seq_len(H)) {
      index <- start[h]:(start[h] + J_h[h] - 1L)
      u_psu[, index] <- sweep(
        z_psu[, index, drop = FALSE],
        1L,
        rowMeans(z_psu[, index, drop = FALSE]),
        "-"
      )
    }

    dimnames(u_state) <- list(draw_id, state_labels)
    dimnames(u_stratum) <- list(draw_id, stratum_labels)
    dimnames(u_psu) <- list(draw_id, psu_labels)
    design <- list(
      row_id = paste0("r", seq_len(N)),
      domain_labels = state_labels,
      stratum_labels = stratum_labels,
      psu_labels = psu_labels,
      state_id = state_id,
      stratum_id = stratum_id,
      psu_flat_id = psu_flat_id,
      w_lik = w_lik,
      population_share = share
    )
    posterior <- list(
      draw_id = draw_id,
      alpha = rnorm(draws, -0.3, 0.8),
      u_state = u_state,
      u_stratum = u_stratum,
      u_psu = u_psu,
      sigma_state = runif(draws, 0, 1.2),
      sigma_stratum = runif(draws, 0, 0.9),
      sigma_psu = runif(draws, 0, 0.8)
    )
    vhat <- stats::setNames(runif(S, 0, 0.002), state_labels)
    list(
      design = design,
      draws = posterior,
      vhat = vhat,
      J_h = J_h,
      psu_start = start
    )
  })
}

.bhf_science_exact_fixture <- function(informative = FALSE) {
  L <- log(3)
  design <- list(
    row_id = paste0("r", 1:8),
    domain_labels = c("A", "B"),
    stratum_labels = c("H1", "H2"),
    psu_labels = c("J1", "J2", "J3", "J4"),
    state_id = c(rep(1L, 4L), rep(2L, 4L)),
    stratum_id = rep(c(1L, 1L, 2L, 2L), 2L),
    psu_flat_id = rep(1:4, 2L),
    w_lik = if (informative) {
      c(2.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 2.5)
    } else {
      rep(1, 8)
    },
    population_share = c(A = 0.5, B = 0.5)
  )
  draws <- list(
    draw_id = "d1",
    alpha = 0,
    u_state = matrix(
      c(-L, L), 1L,
      dimnames = list("d1", c("A", "B"))
    ),
    u_stratum = matrix(
      if (informative) c(L / 2, -L / 2) else c(0, 0),
      1L,
      dimnames = list("d1", c("H1", "H2"))
    ),
    u_psu = matrix(
      if (informative) c(L / 2, -L / 2, L / 2, -L / 2) else rep(0, 4),
      1L,
      dimnames = list("d1", c("J1", "J2", "J3", "J4"))
    ),
    sigma_state = L,
    sigma_stratum = if (informative) L / 2 else 0,
    sigma_psu = if (informative) L / 2 else 0
  )
  list(design = design, draws = draws)
}

.bhf_science_centering_residuals <- function(fixture) {
  shares <- unname(fixture$design$population_share)
  state <- as.vector(fixture$draws$u_state %*% shares)
  stratum <- rowMeans(fixture$draws$u_stratum)
  D <- nrow(fixture$draws$u_psu)
  H <- length(fixture$J_h)
  psu <- matrix(
    0,
    D,
    H,
    dimnames = list(fixture$draws$draw_id,
                    fixture$design$stratum_labels)
  )
  for (h in seq_len(H)) {
    index <- fixture$psu_start[h]:(
      fixture$psu_start[h] + fixture$J_h[h] - 1L
    )
    psu[, h] <- rowMeans(fixture$draws$u_psu[, index, drop = FALSE])
  }
  list(state = state, stratum = stratum, psu = psu)
}

.bhf_science_assert_centered <- function(fixture, tolerance = 1e-12) {
  residual <- .bhf_science_centering_residuals(fixture)
  bad <- which(abs(residual$state) > tolerance)[1L]
  if (!is.na(bad)) {
    stop(
      "state centering failed at draw ", fixture$draws$draw_id[bad],
      "; residual=", residual$state[bad],
      call. = FALSE
    )
  }
  bad <- which(abs(residual$stratum) > tolerance)[1L]
  if (!is.na(bad)) {
    stop(
      "stratum centering failed at draw ", fixture$draws$draw_id[bad],
      "; residual=", residual$stratum[bad],
      call. = FALSE
    )
  }
  bad <- which(abs(residual$psu) > tolerance, arr.ind = TRUE)
  if (nrow(bad)) {
    stop(
      "PSU centering failed at draw ",
      rownames(residual$psu)[bad[1L, "row"]],
      ", stratum ", colnames(residual$psu)[bad[1L, "col"]],
      "; residual=", residual$psu[bad[1L, "row"], bad[1L, "col"]],
      call. = FALSE
    )
  }
  invisible(residual)
}
