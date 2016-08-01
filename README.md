cytofkit: an integrated mass cytometry data analysis pipeline
============

**NOTE**: <u>This is the development version of cytofkit package</u>

### cytofkit

This package is designed to facilitate the analysis workflow of mass cytometry data with automatic subset identification and mapping of cellular progression. Both command line and a GUI client are provided for executing the workflow easily.

### Installation

The offical and stable version, please refer to 

- [Bioconductor](https://www.bioconductor.org/packages/cytofkit/)
- [github](https://github.com/JinmiaoChenLab/cytofkit)

Install the stable version from Bioconductor, use:

``` r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("cytofkit")
```

Install this development version, use:

``` r
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github("haoeric/cytofkit_devel")
```

### Usage

- [cytofkit: Analysis Pipeline](https://www.bioconductor.org/packages/release/bioc/vignettes/cytofkit/inst/doc/cytofkit_workflow.html)    
- [cytofkit: Quick Start](https://www.bioconductor.org/packages/release/bioc/vignettes/cytofkit/inst/doc/cytofkit_example.html)   
- [cytofkit: ShinyAPP Tutorial](https://www.bioconductor.org/packages/release/bioc/vignettes/cytofkit/inst/doc/cytofkit_shinyAPP.html)    
- [shiny APP](https://chenhao.shinyapps.io/cytofkitShinyAPP/) link:  https://chenhao.shinyapps.io/cytofkitShinyAPP/