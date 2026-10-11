# Bundled GAPIT source

`gapit_functions.txt` is the GAPIT source supplied by the EFP code author for this repository. It is included and has not been altered during the upload review. EFP_analysis.R loads this file directly; an installed GAPIT package and the earlier EFP script are not required.

- Internal version reported by `GAPIT.0000()`: **2026.5.19, GAPIT 4.1**.
- SHA-256: `bf65c581e3cce23310e58e12c0f8e9f34670c91a9049c2ca4fabec782a611889`.
- Upstream project and authors: [GAPIT](https://github.com/jiabowang/GAPIT), Jiabo Wang and Zhiwu Zhang, with method contributors credited in the source.
- Upstream source distribution: https://zzlab.net/GAPIT/gapit_functions.txt (a changing URL, not a version pin).
- The precise upstream commit and original download date of the supplied file were not provided. The checksum identifies this exact copy. The internal version alone does not establish which source was used in earlier manuscript analyses.

## Loading and dependencies

The default location is `vendor/gapit_functions.txt` relative to the project root. An existing `EFP_GAPIT_SOURCE` environment variable overrides it. To ensure this bundled copy is used, run `Sys.unsetenv("EFP_GAPIT_SOURCE")` before the installer or analysis. The installer preserves the existing file and sets the variable to this location for that R session.

The bundled file contains package-installation helper functions and checks packages when sourced. Install dependencies first with `install_packages.R`. The installer explicitly includes the core source dependencies plus MASS, RcppArmadillo, snowfall, plotly, htmltools and magrittr for GAPIT routines. It does not edit this source or automatically fetch a newer GAPIT version. The optional `download_gapit=TRUE` is for a deliberately new analysis when this file is absent.

Each use through the EFP loader records its actual path, MD5 and reported version in `outputs/Reproducibility/GAPIT_source_used.csv`. SHA256SUMS.txt also covers the bundled source.

## Attribution and licence notices

The [upstream package DESCRIPTION](https://github.com/jiabowang/GAPIT/blob/master/DESCRIPTION), checked on 9 October 2026, declares **GPL (>= 2)**. A copy of the GPL version 2 text is supplied as [GAPIT_GPL-2.0.txt](GAPIT_GPL-2.0.txt). The source also contains an Apache-licence attribution to MVP developers; that attribution is retained, with the Apache 2.0 text supplied as [MVP_Apache-2.0.txt](MVP_Apache-2.0.txt). These notices concern third-party code and do not select a licence for the EFP scripts or data. The supplied monolithic file does not itself contain a complete licence header or commit identifier; preserve any additional notices received with the original source distribution.

Follow the [upstream citation instructions](https://github.com/jiabowang/GAPIT#citation) for the version actually used. GAPIT 4: Wang and Zhang (2026), *GAPIT Version 4: Integration of GWAS into Genomic Prediction*, [doi:10.1093/molbev/msag107](https://doi.org/10.1093/molbev/msag107). A manuscript statement of GAPIT 3.1.0 does not describe the bundled 4.1 source.
