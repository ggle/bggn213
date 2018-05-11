Class 12
================
Gary Le
May 11, 2018

Hands-on
========

-   Loading in Bio3D and HIV-1 Protease

``` r
library(bio3d)
```

    ## Warning: package 'bio3d' was built under R version 3.4.4

``` r
file.name <- get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

-   Reading PDB file into R

``` r
hiv <- read.pdb(file.name)
```

-   Trimming large PDB object into separate components for easier access

``` r
#Trimming files by protein/ligand
prot <- trim.pdb(hiv, "protein")
lig <- trim.pdb(hiv, "ligand")

#Writing out trimmed files
write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(lig, "1hsg_ligand.pdb")
```

-   Preparing PDB files in AutoDocTools.
    -   AutoDocTools used to add Hydrogen and charge data back to the file since the visualization technique for the protein didn't have enough resolution to resolve hydrogens.
    -   But, we already know where hydrogens *should* be based off of biochemical knowledge.
-   AutoDocTools Grid &gt; GridBox settings
    -   Spacing (angstrom) = 1.000
    -   (x,y,z) Box size = (30, 30, 30)
    -   (x,y,z) center = (16, 25, 4)
-   Scored docking with 1hsg\_ligand using AutoDoc Vina
    -   config = config.txt
    -   log = log.txt
    -   Results stored in all.pdbqt
-   Using Bio3D to prepare a PDB file for VMD to read

``` r
library(bio3d) 
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```

-   Computing Root Mean Square Deviation for each Result with ligand

``` r
res <- read.pdb("all.pdbqt", multi=TRUE) 
ori <- read.pdb("1hsg_ligand.pdbqt")
rmsd(ori, res)
```

    ##  [1]  0.590 11.163 10.531  4.364 11.040  3.682  5.741  3.864  5.442 10.920
    ## [11]  4.318  6.249 11.084  8.929

-   RMSD shows that 1st result deviates by less than 0.6 angstroms.

Q6. RMSD Based on heavy atoms only?

``` r
#Grabbing indices of all non-hydrogen atoms
inds <- atom.select(ori, "noh")

#use $xyz to compare the x,y,z coordinates of the heavy-atoms in results vs the original ligand (in angstroms)
rmsd(ori$xyz[, inds$xyz], res$xyz[, inds$xyz])
```

    ##  [1]  0.458 11.021 10.374  4.301 10.891  3.717  5.764  3.791  5.498 10.759
    ## [11]  4.224  6.308 10.889  8.776

``` r
version <- sessionInfo()

print(version)
```

    ## R version 3.4.3 (2017-11-30)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Sierra 10.12.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] bio3d_2.3-4
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_3.4.3  backports_1.1.2 magrittr_1.5    rprojroot_1.3-2
    ##  [5] parallel_3.4.3  tools_3.4.3     htmltools_0.3.6 yaml_2.1.19    
    ##  [9] Rcpp_0.12.16    stringi_1.2.2   rmarkdown_1.9   grid_3.4.3     
    ## [13] knitr_1.20      stringr_1.3.1   digest_0.6.15   evaluate_0.10.1
