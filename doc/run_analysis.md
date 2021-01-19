# Running the analysis
Make sure you have installed required [software](required_software.md), [packages](packages.md) and downloaded necessary [files](required_files.md).

If you want to replicate the analysis as published, navigate to the `phyloRNAanalysis` folder and type:
```{bash}
Rscript run.r
```

This will remap BAM files to reference genome, detect expression levels and SNV, prepare FASTAs and perform phylogenetic reconstruction.

## Running on your own data
If you want to use this for your own data, put your data into the `data` directory and delete the old data.
Then modify the `chemistry`, `densities`, `hdi` and `selection` to suit your needs.

* `chemistry` -- While the automatic detection of chemistry is preferred, `cellranger` might fail to detect chemistry for low-quality data, for this reason, the chemistry was fixed in the analysis. Either change it to `auto` or to chemistry of your data.
* `densities` -- Filter data to a particular data density, set to a single value if you do not care about comparing phylogenies from different data densities.
* `hdi` -- Highest Density Interval method for discretization of expression values. **This method will only work for homogeneous population of cells**. If your population has heterogeneous expression levels, this method of discretization will not work.
* `selection` -- Named vector to select the best performing cells from each sample for the alternative filtering method used in the study.
