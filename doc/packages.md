# Package installation

## R packages
This analysis requires `tools`, `parallel`, `data.table`, `devtools`, `phyloRNA` and `beter` packages.
Packages `tools` and `parallel` are part of standard install of R do not need to be installed manually.
Packages `data.table` and `devtools` are available from CRAN. First, start `R` console by typing `R` into terminal:

```{bash}
R
```

and now type:

```{r}
install.packages("data.table")
install.packages("devtools")
```

Packages might have additional system requirements that will not be automatically installed by R. You would need to find them in the error log and install them manually. On `ubuntu`, start terminal and run:
```{bash}
sudo apt install name-of-library
```
Alternatively, you could install these packages with the `conda`, start terminal and type:
```{bash}
conda install r-devtools r-data.table
```
This should automatically install all required system libraries.

### phyloRNA and beter
These packages are not available on CRAN. You can however easily install them directly from github using the `devtools` package. Start the `R` console and type:

```{r}
devtools::install_github("biods/phyloRNA")
devtools::install_github("biods/beter")
```

this will automatically install these packages and required dependencies.

## Python packages
To install the `pysam` package, simply type into the terminal:
```{bash}
pip install pysam
```
