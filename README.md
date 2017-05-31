# Streamlined workflow for the quality control, normalization and bias-free analysis of Illumina methylation array data - The Leiden approach #


![alt text](http://www.molepi.nl/images/logo.png)
	  
[![Travis-CI Build Status](https://travis-ci.org/molepi/DNAmArray.svg?branch=master)](https://travis-ci.org/molepi/DNAmArray)
	     
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.158908.svg)](https://doi.org/10.5281/zenodo.158908)


# Workflow #

Check out our compleet workflow: [workflow](https://molepi.github.io/DNAmArray_workflow/)

# Installation #

The **DNAmArray**-package can be installed in several ways which are
described below. The package has been installed successfully for >=
R-3.2.0 on different linux-builds. Currently, two branches are
available:
1. master for >= R-3.4.0 compatible with >= 1.22.1 minfi
2. R-3.3.0  for older R and minfi version's

The package depends on many other packages from
[BioConductor](https://www.bioconductor.org) or
[cran](https://cran.r-project.org/). Usually these are installed
automatically, otherwise, we refer to
[BioConductor](https://www.bioconductor.org/install/) or
[cran](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages)
documentation for the installation of these packages.

## Install using the **devtools**-package ##

First install [**devtools**](https://github.com/hadley/devtools). Next
use:

```{r devtools, eval=FALSE}
library(devtools)
install_github("molepi/DNAmArray") ##for master
install_github("molepi/DNAmArray", ref="R-3.3.0") ##for other branches
```

Sometimes `install_github` fails with CA cert error. Try running `httr::set_config(httr::config( ssl_verifypeer = 0L))` before running `install_github`!

## Install from source using `git/R` ##

Using [git](https://git-scm.com/), e.g, use `git clone` and then build
and install the package from source:

```{r git, engine='bash', eval=FALSE}
git clone git@git.lumc.nl:molepi/DNAmArray.git
R CMD build DNAmArray
R CMD INSTALL DNAmArray_x.y.z.tar.gz
```
Change `_x.y.z.` to the proper version you downloaded!
    
