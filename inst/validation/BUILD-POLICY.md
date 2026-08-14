# Source and generated-content policy

The installable source archive includes `R/`, `man/`, `data/`, `inst/stan/`
Stan source, `inst/validation/`, tests, vignettes, package metadata, and license.

The following are local development/evidence surfaces and are excluded by
`.Rbuildignore`: `.git`, `.github`, IDE state, `dev/`, `log/`,
`output/`, `data-raw/`, rendered `docs/`, local Quarto/build/check directories,
archives, and the platform-specific `inst/stan/bhf_hybrid.rds` cache.

The compiled RDS remains a local pre-existing cache but is not part of the
source archive. Installation and runtime resolve and compile the bundled
`inst/stan/bhf_hybrid.stan`; no stale DSO or absolute path may be required.

Generated Rd and NAMESPACE files are source-controlled build inputs generated
from roxygen. The local pkgdown site is reproducible output, not part of the
source archive and must not be deployed without separate authorization.

`renv` adoption is deferred. Reproducibility is instead recorded with source
hashes, DESCRIPTION constraints, R/session information, exact commands, and a
source-archive manifest. This decision can be revisited before an external
release.
