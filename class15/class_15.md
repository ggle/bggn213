Class 15
================
Gary Le
May 23, 2018

-   With gene enrichment analyses, matching ID's from one database to another can be time consuming.
    -   Some applications like unitprot and excel have tools to help, but we should invest in a more robust method for conversion.

Hands-on
--------

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, cbind, colMeans,
    ##     colnames, colSums, do.call, duplicated, eval, evalq, Filter,
    ##     Find, get, grep, grepl, intersect, is.unsorted, lapply,
    ##     lengths, Map, mapply, match, mget, order, paste, pmax,
    ##     pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
    ##     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
    ##     tapply, union, unique, unsplit, which, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following object is masked from 'package:base':
    ## 
    ##     apply

-   Importing data for DE seq

``` r
countLink <- "https://bioboot.github.io/bimm143_W18/class-material/GSE37704_featurecounts.csv"
metaLink <- "https://bioboot.github.io/bimm143_W18/class-material/GSE37704_metadata.csv"

#Importing metadata
colData <-  read.csv(metaLink, row.names = 1)
#Importing counts data
countData <- read.csv(countLink, row.names = 1)
```

-   Formatting count data to remove unnecessary "length" column and removing non-zero entries.

``` r
#Remove 1st column
countData <- as.matrix(countData[, -1])
#Remove rows with a sum of zero
countData <- countData[ rowSums(countData) > 1 ,]
```

-   Setting up DESeq

``` r
dds <-  DESeqDataSetFromMatrix(countData=countData,
                             colData=colData,
                             design=~condition)
dds <-  DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

-   Checking results

``` r
res <- results(dds)
res
```

    ## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
    ## Wald test p-value: condition hoxa1 kd vs control sirna 
    ## DataFrame with 15280 rows and 6 columns
    ##                   baseMean log2FoldChange      lfcSE        stat
    ##                  <numeric>      <numeric>  <numeric>   <numeric>
    ## ENSG00000279457   29.91358     0.17927483 0.32459294   0.5523066
    ## ENSG00000187634  183.22965     0.42644724 0.14017817   3.0421802
    ## ENSG00000188976 1651.18808    -0.69272061 0.05484412 -12.6307172
    ## ENSG00000187961  209.63794     0.72975918 0.13178350   5.5375609
    ## ENSG00000187583   47.25512     0.04055411 0.27169055   0.1492658
    ## ...                    ...            ...        ...         ...
    ## ENSG00000273748  35.302652      0.6743994  0.3034582   2.2223801
    ## ENSG00000278817   2.423024     -0.3889516  1.1295943  -0.3443286
    ## ENSG00000278384   1.101796      0.3328870  1.6590966   0.2006435
    ## ENSG00000276345  73.644956     -0.3561673  0.2075751  -1.7158482
    ## ENSG00000271254 181.595903     -0.6096640  0.1412340  -4.3166951
    ##                       pvalue         padj
    ##                    <numeric>    <numeric>
    ## ENSG00000279457 5.807383e-01 6.846746e-01
    ## ENSG00000187634 2.348712e-03 5.109223e-03
    ## ENSG00000188976 1.429691e-36 1.745816e-35
    ## ENSG00000187961 3.067131e-08 1.109758e-07
    ## ENSG00000187583 8.813439e-01 9.191354e-01
    ## ...                      ...          ...
    ## ENSG00000273748 2.625763e-02 4.756160e-02
    ## ENSG00000278817 7.305992e-01 8.086868e-01
    ## ENSG00000278384 8.409773e-01 8.927559e-01
    ## ENSG00000276345 8.618983e-02 1.389975e-01
    ## ENSG00000271254 1.583827e-05 4.470014e-05

-   Reformating results

``` r
mcols(res, use.names = TRUE) 
```

    ## DataFrame with 6 rows and 2 columns
    ##                        type
    ##                 <character>
    ## baseMean       intermediate
    ## log2FoldChange      results
    ## lfcSE               results
    ## stat                results
    ## pvalue              results
    ## padj                results
    ##                                                                description
    ##                                                                <character>
    ## baseMean                         mean of normalized counts for all samples
    ## log2FoldChange log2 fold change (MLE): condition hoxa1 kd vs control sirna
    ## lfcSE                  standard error: condition hoxa1 kd vs control sirna
    ## stat                   Wald statistic: condition hoxa1 kd vs control sirna
    ## pvalue              Wald test p-value: condition hoxa1 kd vs control sirna
    ## padj                                                  BH adjusted p-values

-   Summarize results

``` r
summary(res)
```

    ## 
    ## out of 15280 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)     : 4352, 28% 
    ## LFC < 0 (down)   : 4400, 29% 
    ## outliers [1]     : 0, 0% 
    ## low counts [2]   : 590, 3.9% 
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

-   Too many results, raising thresholds

``` r
#Defining significant results as alpha < 0.05% and log2 threshold > 2
resSig <- results(dds, alpha = 0.05, lfcThreshold = 2)
summary(resSig)
```

    ## 
    ## out of 15280 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)     : 99, 0.65% 
    ## LFC < 0 (down)   : 134, 0.88% 
    ## outliers [1]     : 0, 0% 
    ## low counts [2]   : 1482, 9.7% 
    ## (mean count < 2)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

-   Can see that only 99 results are sig. rather than 4352

-   Alternative way to screen the results for the same alpha and log2 params.

``` r
#Results that have alpha < 0.05%
table(resSig$padj < 0.05)
```

    ## 
    ## FALSE  TRUE 
    ## 13565   233

``` r
#Results that have log2 threshold > 2
table(resSig$log2FoldChange > 2)
```

    ## 
    ## FALSE  TRUE 
    ## 14723   557

``` r
#Crosstabulation of results
table(resSig$padj < 0.05, resSig$log2FoldChange > 2)
```

    ##        
    ##         FALSE  TRUE
    ##   FALSE 13292   273
    ##   TRUE    134    99

-   Redefine resSig to help order results for further analysis

``` r
# Make a new results object 'resSig' with only significant genes
resSig <- subset(res, padj < 0.05)

# Order resSig by log2FoldChange
# Print the first 10 strongest DOWN genes
ord.down <- order(resSig$log2FoldChange)
head(resSig[ ord.down, ], 10)
```

    ## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
    ## Wald test p-value: condition hoxa1 kd vs control sirna 
    ## DataFrame with 10 rows and 6 columns
    ##                  baseMean log2FoldChange     lfcSE      stat       pvalue
    ##                 <numeric>      <numeric> <numeric> <numeric>    <numeric>
    ## ENSG00000171587  2.561658      -4.902903  1.514295 -3.237746 1.204782e-03
    ## ENSG00000155011  7.971242      -4.601385  1.147782 -4.008936 6.099288e-05
    ## ENSG00000101306  1.888372      -4.425978  1.644059 -2.692105 7.100260e-03
    ## ENSG00000179855  3.760403      -4.424858  1.400644 -3.159161 1.582241e-03
    ## ENSG00000188581  1.795033      -4.405196  1.637739 -2.689804 7.149406e-03
    ## ENSG00000164342  3.517655      -4.323085  1.415825 -3.053403 2.262620e-03
    ## ENSG00000117598  1.680397      -4.261884  1.686659 -2.526821 1.151003e-02
    ## ENSG00000180332  9.377022      -4.258232  0.966178 -4.407296 1.046691e-05
    ## ENSG00000136943  3.367218      -4.250652  1.454036 -2.923347 3.462904e-03
    ## ENSG00000204183  1.542573      -4.168844  1.679468 -2.482241 1.305591e-02
    ##                         padj
    ##                    <numeric>
    ## ENSG00000171587 0.0027291049
    ## ENSG00000155011 0.0001615844
    ## ENSG00000101306 0.0142919727
    ## ENSG00000179855 0.0035179536
    ## ENSG00000188581 0.0143825517
    ## ENSG00000164342 0.0049307055
    ## ENSG00000117598 0.0224039088
    ## ENSG00000180332 0.0000301340
    ## ENSG00000136943 0.0073542377
    ## ENSG00000204183 0.0251331755

``` r
# Print the first 10 strongest UP genes
ord.up <- order(resSig$log2FoldChange, decreasing = TRUE)
head(resSig[ ord.up, ], 10)
```

    ## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
    ## Wald test p-value: condition hoxa1 kd vs control sirna 
    ## DataFrame with 10 rows and 6 columns
    ##                    baseMean log2FoldChange     lfcSE      stat
    ##                   <numeric>      <numeric> <numeric> <numeric>
    ## ENSG00000128052  158.249550       8.822086 1.0347214  8.526049
    ## ENSG00000141668   43.941562       8.820638 1.2040490  7.325813
    ## ENSG00000004799   97.260234       8.117454 1.0373919  7.824868
    ## ENSG00000109321 3907.190542       7.479166 0.2880675 25.963240
    ## ENSG00000162892   25.278929       7.052106 1.2210031  5.775666
    ## ENSG00000073756 1357.627066       6.337747 0.3399843 18.641296
    ## ENSG00000239590    7.399064       6.250557 1.2943354  4.829164
    ## ENSG00000151322    7.309569       6.237074 1.2988709  4.801920
    ## ENSG00000135363    6.035727       5.955831 1.3350911  4.460992
    ## ENSG00000175868    5.981592       5.945427 1.3207432  4.501577
    ##                        pvalue          padj
    ##                     <numeric>     <numeric>
    ## ENSG00000128052  1.514306e-17  9.338857e-17
    ## ENSG00000141668  2.374553e-13  1.195414e-12
    ## ENSG00000004799  5.081914e-15  2.783494e-14
    ## ENSG00000109321 1.288826e-148 1.147446e-146
    ## ENSG00000162892  7.664930e-09  2.907251e-08
    ## ENSG00000073756  1.485963e-77  4.409858e-76
    ## ENSG00000239590  1.371077e-06  4.320274e-06
    ## ENSG00000151322  1.571514e-06  4.919145e-06
    ## ENSG00000135363  8.158097e-06  2.375470e-05
    ## ENSG00000175868  6.745119e-06  1.982349e-05

-   Organizing from ENSEMBLE to KEGG IDs

``` r
# N.B. these are BioConductor pkgs.
library("AnnotationDbi")
library("org.Hs.eg.db")
```

    ## 

``` r
columns(org.Hs.eg.db)
```

    ##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
    ##  [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
    ##  [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
    ## [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
    ## [17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
    ## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    ## [25] "UNIGENE"      "UNIPROT"

-   This is the organism annotation package (“org”) for Homo sapiens (“Hs”), organized as an AnnotationDbi database package (“db”), using Entrez Gene IDs (“eg”) as primary key.

-   Adding SYMBOL, ENTREZ, and GENENAME IDs to results

``` r
res$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
res$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
res$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
# Displaying top 10 again
head(res, 10)
```

    ## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
    ## Wald test p-value: condition hoxa1 kd vs control sirna 
    ## DataFrame with 10 rows and 9 columns
    ##                   baseMean log2FoldChange      lfcSE        stat
    ##                  <numeric>      <numeric>  <numeric>   <numeric>
    ## ENSG00000279457   29.91358     0.17927483 0.32459294   0.5523066
    ## ENSG00000187634  183.22965     0.42644724 0.14017817   3.0421802
    ## ENSG00000188976 1651.18808    -0.69272061 0.05484412 -12.6307172
    ## ENSG00000187961  209.63794     0.72975918 0.13178350   5.5375609
    ## ENSG00000187583   47.25512     0.04055411 0.27169055   0.1492658
    ## ENSG00000187642   11.97975     0.54275443 0.52117592   1.0414035
    ## ENSG00000188290  108.92213     2.05704993 0.19679801  10.4525953
    ## ENSG00000187608  350.71687     0.25737707 0.10267368   2.5067482
    ## ENSG00000188157 9128.43942     0.38990881 0.04673847   8.3423530
    ## ENSG00000131591  156.47908     0.19659758 0.14551980   1.3510022
    ##                       pvalue         padj       symbol      entrez
    ##                    <numeric>    <numeric>  <character> <character>
    ## ENSG00000279457 5.807383e-01 6.846746e-01 LOC102723897   102723897
    ## ENSG00000187634 2.348712e-03 5.109223e-03       SAMD11      148398
    ## ENSG00000188976 1.429691e-36 1.745816e-35        NOC2L       26155
    ## ENSG00000187961 3.067131e-08 1.109758e-07       KLHL17      339451
    ## ENSG00000187583 8.813439e-01 9.191354e-01      PLEKHN1       84069
    ## ENSG00000187642 2.976883e-01 4.016657e-01        PERM1       84808
    ## ENSG00000188290 1.425695e-25 1.226198e-24         HES4       57801
    ## ENSG00000187608 1.218475e-02 2.358287e-02        ISG15        9636
    ## ENSG00000188157 7.282667e-17 4.347110e-16         AGRN      375790
    ## ENSG00000131591 1.766947e-01 2.604501e-01     C1orf159       54991
    ##                                                                     name
    ##                                                              <character>
    ## ENSG00000279457                        WAS protein family homolog 2-like
    ## ENSG00000187634                 sterile alpha motif domain containing 11
    ## ENSG00000188976 NOC2 like nucleolar associated transcriptional repressor
    ## ENSG00000187961                              kelch like family member 17
    ## ENSG00000187583                 pleckstrin homology domain containing N1
    ## ENSG00000187642             PPARGC1 and ESRR induced regulator, muscle 1
    ## ENSG00000188290                   hes family bHLH transcription factor 4
    ## ENSG00000187608                            ISG15 ubiquitin-like modifier
    ## ENSG00000188157                                                    agrin
    ## ENSG00000131591                      chromosome 1 open reading frame 159

### Section 2

-   First time install

``` r
#source("http://bioconductor.org/biocLite.R")
#biocLite( c("pathview", "gage", "gageData") )

library(pathview)
```

    ## ##############################################################################
    ## Pathview is an open source software package distributed under GNU General
    ## Public License version 3 (GPLv3). Details of GPLv3 is available at
    ## http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
    ## formally cite the original Pathview paper (not just mention it) in publications
    ## or products. For details, do citation("pathview") within R.
    ## 
    ## The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
    ## license agreement (details at http://www.kegg.jp/kegg/legal.html).
    ## ##############################################################################

``` r
library(gage)
library(gageData)

data(kegg.sets.hs)
data(sigmet.idx.hs)

# Signaling and Metabolic Pathway subset of KEGG pathways by ENTREZ ID
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)
```

    ## $`hsa00232 Caffeine metabolism`
    ## [1] "10"   "1544" "1548" "1549" "1553" "7498" "9"   
    ## 
    ## $`hsa00983 Drug metabolism - other enzymes`
    ##  [1] "10"     "1066"   "10720"  "10941"  "151531" "1548"   "1549"  
    ##  [8] "1551"   "1553"   "1576"   "1577"   "1806"   "1807"   "1890"  
    ## [15] "221223" "2990"   "3251"   "3614"   "3615"   "3704"   "51733" 
    ## [22] "54490"  "54575"  "54576"  "54577"  "54578"  "54579"  "54600" 
    ## [29] "54657"  "54658"  "54659"  "54963"  "574537" "64816"  "7083"  
    ## [36] "7084"   "7172"   "7363"   "7364"   "7365"   "7366"   "7367"  
    ## [43] "7371"   "7372"   "7378"   "7498"   "79799"  "83549"  "8824"  
    ## [50] "8833"   "9"      "978"   
    ## 
    ## $`hsa00230 Purine metabolism`
    ##   [1] "100"    "10201"  "10606"  "10621"  "10622"  "10623"  "107"   
    ##   [8] "10714"  "108"    "10846"  "109"    "111"    "11128"  "11164" 
    ##  [15] "112"    "113"    "114"    "115"    "122481" "122622" "124583"
    ##  [22] "132"    "158"    "159"    "1633"   "171568" "1716"   "196883"
    ##  [29] "203"    "204"    "205"    "221823" "2272"   "22978"  "23649" 
    ##  [36] "246721" "25885"  "2618"   "26289"  "270"    "271"    "27115" 
    ##  [43] "272"    "2766"   "2977"   "2982"   "2983"   "2984"   "2986"  
    ##  [50] "2987"   "29922"  "3000"   "30833"  "30834"  "318"    "3251"  
    ##  [57] "353"    "3614"   "3615"   "3704"   "377841" "471"    "4830"  
    ##  [64] "4831"   "4832"   "4833"   "4860"   "4881"   "4882"   "4907"  
    ##  [71] "50484"  "50940"  "51082"  "51251"  "51292"  "5136"   "5137"  
    ##  [78] "5138"   "5139"   "5140"   "5141"   "5142"   "5143"   "5144"  
    ##  [85] "5145"   "5146"   "5147"   "5148"   "5149"   "5150"   "5151"  
    ##  [92] "5152"   "5153"   "5158"   "5167"   "5169"   "51728"  "5198"  
    ##  [99] "5236"   "5313"   "5315"   "53343"  "54107"  "5422"   "5424"  
    ## [106] "5425"   "5426"   "5427"   "5430"   "5431"   "5432"   "5433"  
    ## [113] "5434"   "5435"   "5436"   "5437"   "5438"   "5439"   "5440"  
    ## [120] "5441"   "5471"   "548644" "55276"  "5557"   "5558"   "55703" 
    ## [127] "55811"  "55821"  "5631"   "5634"   "56655"  "56953"  "56985" 
    ## [134] "57804"  "58497"  "6240"   "6241"   "64425"  "646625" "654364"
    ## [141] "661"    "7498"   "8382"   "84172"  "84265"  "84284"  "84618" 
    ## [148] "8622"   "8654"   "87178"  "8833"   "9060"   "9061"   "93034" 
    ## [155] "953"    "9533"   "954"    "955"    "956"    "957"    "9583"  
    ## [162] "9615"

-   GAGE fxn requires named vector of fold changes

``` r
#Vector of fold changes from results
foldchanges <-  res$log2FoldChange
#Naming foldchanges by entrez ids from results
names(foldchanges) <-  res$entrez

head(foldchanges)
```

    ##   102723897      148398       26155      339451       84069       84808 
    ##  0.17927483  0.42644724 -0.69272061  0.72975918  0.04055411  0.54275443

-   GAGE results on the KEGG subset of genes

``` r
# Get the results
keggres <-  gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
```

-   Looking at head of each attribute in keggres

``` r
lapply(keggres, head)
```

    ## $greater
    ##                                         p.geomean stat.mean       p.val
    ## hsa04640 Hematopoietic cell lineage   0.002709366  2.857393 0.002709366
    ## hsa04630 Jak-STAT signaling pathway   0.005655916  2.557207 0.005655916
    ## hsa04142 Lysosome                     0.008948808  2.384783 0.008948808
    ## hsa00140 Steroid hormone biosynthesis 0.009619717  2.432105 0.009619717
    ## hsa04740 Olfactory transduction       0.010846448  2.357877 0.010846448
    ## hsa04916 Melanogenesis                0.017613266  2.123151 0.017613266
    ##                                           q.val set.size        exp1
    ## hsa04640 Hematopoietic cell lineage   0.3470863       49 0.002709366
    ## hsa04630 Jak-STAT signaling pathway   0.3470863      103 0.005655916
    ## hsa04142 Lysosome                     0.3470863      117 0.008948808
    ## hsa00140 Steroid hormone biosynthesis 0.3470863       26 0.009619717
    ## hsa04740 Olfactory transduction       0.3470863       39 0.010846448
    ## hsa04916 Melanogenesis                0.4696871       85 0.017613266
    ## 
    ## $less
    ##                                      p.geomean stat.mean        p.val
    ## hsa04110 Cell cycle               1.004024e-05 -4.353447 1.004024e-05
    ## hsa03030 DNA replication          8.909718e-05 -3.968605 8.909718e-05
    ## hsa03013 RNA transport            1.333746e-03 -3.038773 1.333746e-03
    ## hsa04114 Oocyte meiosis           2.496109e-03 -2.840329 2.496109e-03
    ## hsa03440 Homologous recombination 2.942017e-03 -2.868137 2.942017e-03
    ## hsa00240 Pyrimidine metabolism    5.800212e-03 -2.549616 5.800212e-03
    ##                                         q.val set.size         exp1
    ## hsa04110 Cell cycle               0.001606438      120 1.004024e-05
    ## hsa03030 DNA replication          0.007127774       36 8.909718e-05
    ## hsa03013 RNA transport            0.071133118      143 1.333746e-03
    ## hsa04114 Oocyte meiosis           0.094144560       98 2.496109e-03
    ## hsa03440 Homologous recombination 0.094144560       28 2.942017e-03
    ## hsa00240 Pyrimidine metabolism    0.138500584       95 5.800212e-03
    ## 
    ## $stats
    ##                                       stat.mean     exp1
    ## hsa04640 Hematopoietic cell lineage    2.857393 2.857393
    ## hsa04630 Jak-STAT signaling pathway    2.557207 2.557207
    ## hsa04142 Lysosome                      2.384783 2.384783
    ## hsa00140 Steroid hormone biosynthesis  2.432105 2.432105
    ## hsa04740 Olfactory transduction        2.357877 2.357877
    ## hsa04916 Melanogenesis                 2.123151 2.123151

-   Sanity check

``` r
## Sanity check displaying all pathways data
pathways = data.frame(id=rownames(keggres$greater), keggres$greater)
head(pathways)
```

    ##                                                                          id
    ## hsa04640 Hematopoietic cell lineage     hsa04640 Hematopoietic cell lineage
    ## hsa04630 Jak-STAT signaling pathway     hsa04630 Jak-STAT signaling pathway
    ## hsa04142 Lysosome                                         hsa04142 Lysosome
    ## hsa00140 Steroid hormone biosynthesis hsa00140 Steroid hormone biosynthesis
    ## hsa04740 Olfactory transduction             hsa04740 Olfactory transduction
    ## hsa04916 Melanogenesis                               hsa04916 Melanogenesis
    ##                                         p.geomean stat.mean       p.val
    ## hsa04640 Hematopoietic cell lineage   0.002709366  2.857393 0.002709366
    ## hsa04630 Jak-STAT signaling pathway   0.005655916  2.557207 0.005655916
    ## hsa04142 Lysosome                     0.008948808  2.384783 0.008948808
    ## hsa00140 Steroid hormone biosynthesis 0.009619717  2.432105 0.009619717
    ## hsa04740 Olfactory transduction       0.010846448  2.357877 0.010846448
    ## hsa04916 Melanogenesis                0.017613266  2.123151 0.017613266
    ##                                           q.val set.size        exp1
    ## hsa04640 Hematopoietic cell lineage   0.3470863       49 0.002709366
    ## hsa04630 Jak-STAT signaling pathway   0.3470863      103 0.005655916
    ## hsa04142 Lysosome                     0.3470863      117 0.008948808
    ## hsa00140 Steroid hormone biosynthesis 0.3470863       26 0.009619717
    ## hsa04740 Olfactory transduction       0.3470863       39 0.010846448
    ## hsa04916 Melanogenesis                0.4696871       85 0.017613266

-   Using pathways() to visualize our pathways

``` r
pathview(gene.data=foldchanges, pathway.id="hsa04110")
```

    ## 'select()' returned 1:1 mapping between keys and columns

    ## Info: Working in directory C:/Users/gary_/Desktop/R_Bootcamp/bggn213_github/class15

    ## Info: Writing image file hsa04110.pathview.png

-   Looking only at the 5 genes that had the greatest change

``` r
## Focus on top 5 upregulated pathways here for demo purposes only
keggrespathways <- rownames(keggres$greater)[1:5]

# Extract the IDs part of each string
keggresids <-  substr(keggrespathways, start=1, stop=8)
keggresids
```

    ## [1] "hsa04640" "hsa04630" "hsa04142" "hsa00140" "hsa04740"

``` r
#pathview(gene.data=foldchanges, pathway.id=keggresids, species="hsa")
```

![](hsa04640.pathview.png)
