# DID-INT in Stata

Estimate a policy effect when covariate relationships can differ across schools and years.

DID-INT fits school-year outcome models and evaluates comparisons at the treated school’s covariate distribution. It then aggregates the supported school-time effects. Standardizing fitted outcomes is important: differences between regression intercepts alone would depend on how covariates were centered.

## Install

Requires **Stata 16 or later**. No R, Python, or additional Stata estimation
packages are required to use this command.

Install the command directly from GitHub:

```stata
net install suite_didint, from("https://raw.githubusercontent.com/mattdwebb/suitedidintstata/main/") replace
discard
which suite_didint
```

For a source-only trial, `adopath ++ "."` makes this checkout available for the
current session. Use one known version and record it for reproducible work.

## A school-policy example

The included fictional data contain 960 pupil observations: ten pupils in each of 12 schools and eight years (2010–2017). Pupils are repeated cross sections, so their row numbers are not longitudinal identifiers.
Four schools adopt in 2013, four in 2015, and four never adopt. The outcome is
`test_score`, measured in score points. `baseline_score` is a predetermined
covariate. Schools are the independent inference clusters.

The data generator adds a two-point policy effect, alongside school, year,
covariate, and varying disturbance components. This is a software teaching
example, not evidence about an actual policy. The [data guide](examples/DATA_GUIDE.md)
defines every variable and links the reproducible generator.

Run from the extracted repository folder:

```stata
clear
set type double
insheet using "examples/school_policy.csv", comma names clear
suite_didint test_score baseline_score, time(year) cluster(school) gvar(first_treat) inference(jackknife)
display e(estimate)
display e(jackknife_se)
```

The executed example gives:

| Estimate (score points) | Native clustered SE | Jackknife SE |
|---:|---:|---:|
| 2.394744 | Not reported | 0.519843 |

The estimated effect is about **2.39 score points**. The jackknife
standard error is **0.52 points**, calculated by leaving each of the twelve
schools out in turn. It describes uncertainty around this estimate; it does not
change the point estimate or establish the identifying assumptions of a policy
study. Standalone jackknife p-values use Student t with 11 degrees of freedom here.

`do examples/school_policy.do` runs this example, prints a compact table, and
saves `examples/school_results.csv`. Its [recorded output](examples/expected_output.txt)
is checked against the release example. Small last-digit differences may depend
on numerical libraries; they should not change the substantive result.

## Adapt it to your data

Use zero for never-treated schools and repeat each school’s first treatment year on all its rows. Treatment must be absorbing and common to observations in a school-year.

Covariates follow the outcome in the variable list. Omit `baseline_score` for an unadjusted specification. Use numeric covariates whose timing makes sense for your study;
the command does not make post-treatment controls causally appropriate.

Request `inference(jackknife)` for the delete-school calculation. The native analytic option is unavailable for this aggregate.

When a deletion removes a cohort, the jackknife recomputes the estimator’s
supported aggregate and its weights using the surviving estimable effects.
All cluster deletions must yield a scalar estimate. A genuinely unidentified
deletion is reported as a failure; it is never silently omitted.

Coefficient randomization inference and the supported studentized variants use
alternative cohort assignments while preserving treatment counts. These require
more input than the introductory example; see [inference options](INFERENCE.md).
Changing the treatment assignment must also change any treatment-dependent fits.

## Supported scope

This implementation supports the fully interacted covariate specification (CCC int), never-treated controls, and finite required row values. It aggregates identified surviving state-year comparisons; an absent cell does not automatically invalidate the entire aggregate. Other CCC modes and observation weights are unsupported. There is no separate analytic clustered standard error for this aggregate.

Different DiD estimators can target different weighted averages when effects
vary. Comparing their coefficients does not by itself show that one is correct.
Unsupported options return explicit errors rather than substituting another model.

## Help, methods, and validation

- [Estimator definition and references](METHOD_SPECIFICATION.md)
- [Inference definitions](INFERENCE.md)
- [Technical development reference](docs/DEVELOPMENT_REFERENCE.md)
- [Runnable example](examples/school_policy.do) and [data definitions](examples/DATA_GUIDE.md)
- [License](LICENSE)

Use the repository’s Issues page for a reproducible problem report. Include your
software version, a small artificial dataset, and the exact command. Do not post
confidential research data. For developers, `tests/run_release_tests.ps1` runs the
native Stata checks and records their status with a licensed Stata installation.
