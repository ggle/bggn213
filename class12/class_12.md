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
