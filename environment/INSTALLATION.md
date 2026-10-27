# Software needed for the supplied code

## Plotting

Figure 1 uses Python 3 and its standard library. It writes SVG from the supplied JSON layout. Open the SVG in a vector-graphics application if a different format is needed.

R plotting was tested with R 4.5.3. Required packages are `data.table`, `ggplot2`, `patchwork`, `gtable`, `ragg`, `systemfonts`, `scales`, `ggrepel`, and `R.utils`. Supplementary figures S2 and S3 additionally use `ggrastr` for point-layer rasterisation in PDF. Install them in a normal R library, for example:

```r
install.packages(c("data.table", "ggplot2", "patchwork", "gtable", "ragg",
                   "systemfonts", "scales", "ggrepel", "R.utils", "ggrastr"))
```

This installs current packages, not an exact historical environment. Recorded session files describe the packages actually used for analysis and the separate plotting portability check. System requirements for compiling these packages differ by operating system; an installation failure should be resolved before starting a figure. The independent plotting scripts do not silently install or upgrade packages. Some retained historical analysis/preparation scripts contain conditional installation attempts; install the documented versions first and inspect the selected step before running it.

Arial is required to preserve typography. Install a legally obtained copy in the operating system and check it with `systemfonts::match_fonts("Arial")`. Scripts check the resolved font family instead of assuming a Windows font-file path. Cairo support is needed for R PDF output; check `capabilities("cairo")`.

R discovers packages through its standard library configuration (`.libPaths()`, `R_LIBS_USER`). No private project library directory or internal Skill package is required. `--check` reports missing packages, input files or Arial before plotting.

## Statistical analyses

See the recorded analysis session files. The raw-count sensitivity analysis used R 4.5.3 and edgeR 4.8.2. The individual primary analysis steps include additional dependencies for single-cell processing, normalization and statistical testing. They retain the original staged project layout and intermediate-file requirements; the plotting check is not proof of an end-to-end primary-analysis rerun.

Do not use the plotting environment as evidence that every historical analysis dependency has been restored. Analysis packages, data, parameters and expected stage outputs must be reviewed separately before a new computational run.
