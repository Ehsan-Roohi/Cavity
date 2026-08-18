# Open research status

Last audit: 2026-08-18.

| PR | Work | Current status | Gate before merge |
| ---: | --- | --- | --- |
| [#1](https://github.com/Ehsan-Roohi/Cavity/pull/1) | native full-cavity Shakhov DVM for the JFM thermal case | draft; syntax, parser, and wall-flux formula checks passed, but no local executable Fortran smoke run was available | complete the Unity production run and require metadata `converged: true`; retain the distinction between the BGK/Shakhov model equation and a full hard-sphere Boltzmann/DSMC reference |

Do not merge this draft solely to reduce the PR count.  Preserve any
non-converged archive as diagnostic evidence.
