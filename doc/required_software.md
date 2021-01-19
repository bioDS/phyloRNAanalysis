# Software installation
The analysis requires [R](https://cran.r-project.org/), [python3](https://www.python.org/), [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger), [bamtofastq](https://github.com/10XGenomics/bamtofastq), [GATK](https://gatk.broadinstitute.org/hc/en-us), [VCFtools](https://vcftools.github.io/), [IQtree](http://www.iqtree.org/) and [BEAST2](https://www.beast2.org/)

Required software can be installed independently, but we suggest the [conda](https://docs.conda.io/en/latest/) package manager. `Cellranger` and `bamtofastq` are not available through the `conda` package manager and need to be installed separately.

## Conda
To install the `conda` package manager, either follow the instructions bellow or, for more information, instructions on the conda [website](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html).

To download `miniconda` (minimal installation of `conda`), open your terminal and type the following:

```{bash}
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

This will download the `miniconda` into your current directory. Now type:
```{bash}
bash Miniconda3-latest-Linux-x86_64.sh
```

and follow the installation instructions in your terminal. This will install `conda` into your computer. To verify that `conda` is installed, open new terminal window (or type `exec bash`) and type:

```{bash}
conda --help
```

You can now install required software by typing:

```{bash}
conda install r-base gatk4 vcftools iqtree beast2
```

There is no need to install `python3` as it is a part of the `conda` installation.

## Cellranger
Cellranger is a proprietary software and its not available through the `conda` manager.

Visit the `Cellranger` [homepage](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) and fill in required information. You will be provided with a download link and [detailed installation instructions](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation).

Download provided link with `curl` or `wget`. You will download a compressed archive `cellranger-5.0.1.tar.gz` (or a newer version). To unpack this archive, type:
```{bash}
tar -xzvf cellranger-5.0.1.tar.gz
```

This will unpack the archive in your current directory. Personally, You can now remove `cellranger-5.0.1.tar.gz` by typing:

```{bash}
rm cellranger-5.0.1.tar.gz
```

Now open `.bashrc` or `.profile` in your favourite editor and add `export PATH=~/cellranger-5.0.1:$PATH` or type:

```
echo 'export PATH=~/cellranger-5.0.1:$PATH' >> .bashrc
```

Restart your terminal by typing `exec bash` and type:

```
cellranger --help
```
to verify the installation.

## bamtofastq
To install `bamtofastq` visit the project [github](https://github.com/10XGenomics/bamtofastq) page or follow these instructions:

To download `bamtofastq`, type:

```{bash}
wget -O bamtofastq https://github.com/10XGenomics/bamtofastq/releases/download/v1.3.5/bamtofastq_linux
```

Now you need to add execute (x) permission to `bamtofastq`:

```{bash}
chmod +x bamtofastq
```

Finally, add `bamtofastq` to your PATH like you did with Cellranger:

```{bash}
echo 'export PATH=~/cellranger-5.0.1:$PATH' >> .bashrc
```

Restart your terminal (`exec bash`) and type:

```{bash}
bamtofastq --help
```
to verify the installation.
