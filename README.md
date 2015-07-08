# coveRage
Tools to explore sequencing coverage in high-throughput sequencing projects

[![Travis-CI Build Status](https://travis-ci.org/knausb/coveRage.png?branch=master)](https://travis-ci.org/knausb/coveRage)

This package was called `covR` for a brief time.
Due to a name conflict with a package already at CRAN ([covr](https://github.com/search?utf8=%E2%9C%93&q=covR)) I've renamed it to coveRage.

While this project is in development it can be installed through github:

    devtools::install_github(repo="knausb/coveRage")
    library(coveRage)


If you would like the vignettes use:

  devtools::install_github(repo="knausb/vcfR", build_vignettes=TRUE)

You can browse the vignettes with:

  browseVignettes(package="coveRage")


If you'd like to try the development branch (which may not be stable) use:

  devtools::install_github(repo="knausb/coveRage@devel")

