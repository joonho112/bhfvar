# Test tiers

## Tier 1 — default local/package checks

Runs all deterministic unit, invariant, property, schema, documentation,
migration, and static Stan tests through `testthat`. Native Stan sampling is
disabled. This tier must have no failures, errors, warnings, or unexplained
skips.

Command:

```sh
Rscript -e 'devtools::test()'
```

## Tier 2 — native Stan parity smoke

Sets `BHFVAR_RUN_STAN_PARITY=true` and runs the tiny frozen oracle parity test.
It compiles and samples one small fixture. This is an implementation parity
check, not recovery validation.

```sh
BHFVAR_RUN_STAN_PARITY=true Rscript -e 'devtools::test(filter="stan-oracle-parity")'
```

## Tier 3 — functional clean-install smoke

Builds and installs the source archive into an empty temporary library, loads
the package, prepares data, compiles the bundled source, runs a tiny fit, and
extracts results. It checks path/DSO/runtime operation only and cannot be cited
as scientific evidence.

## Tier 4 — simulation studies

Recovery and sensitivity jobs are not run by default, by package checks, or by
CI. They are long-running simulation studies driven from `dev/`, and their
results are reported in `SCIENTIFIC-LIMITATIONS.md` rather than asserted by the
test suite.
