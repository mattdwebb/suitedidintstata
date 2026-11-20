# Native Stata DIDINT suite adapter

Install using `net install suite_didint, from("<local repository path>")`.
No R, Python, Julia, Java, reghdfe, or ftools runtime dependency.
Independent native target-standardized fitted DID-INT implementation, not the official Julia wrapper.



Use `treatment(D)` instead of or consistent with gvar. Never treated = 0.
Use `ivar(id) datamode(panel)` for panels, `datamode(rcs)` for repeated cross
sections. Covariates enter as supplied; no arbitrary baseline freezing.
DID-INT supports CCC int only. Other CCC modes error explicitly.

Run `do tests/suite_tests.do` from the repository root. Validation is a
supported-subset check, not full official feature parity or certification.
See METHOD_SPECIFICATION.md, INFERENCE.md and VALIDATION.md.
