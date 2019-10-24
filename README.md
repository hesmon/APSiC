# APSiC
The whole data analayzed in this paper are published as a part of project DRIVE (https://data.mendeley.com/datasets/y3ds55n88r/5) and at the CCLE (https://portals.broadinstitute.org/ccle). The breast cancer data is available as an R object in this repository (`R/BreastData.RData`). 

To find non-genetic tumor suppressors using APSiC, use the `identifyDependencies` function with argument `dependencyType = "non-genetic-tsg`:
```
source("apsic_common_functions.r")
load("BreastData.RData")

# The APSiC for detecting non-genetic tumor suppressors
identifyDependencies(breastData, dependencyType = "non-genetic-tsg")
```
The `dependencyType` argument takes values from `non-genetic-tsg`, `non-genetic-oncogene`, 
`mutation-tsg`, `mutation-oncogene`, `amplification-oncogene` to find the cancer dependencies defined in the manuscript.


To plot the rank profiles, the `waterfallForGene` function can be used:
```
gene = "TP53"
waterfallForGene(breastData, gene=gene, title=paste("Breast cancer: rank profile of", gene), rank=TRUE)
```

# Manuscript
The manuscript is available at:
[https://www.biorxiv.org/content/10.1101/807248v1](https://www.biorxiv.org/content/10.1101/807248v1)

# P-value files
APSiC p-values for 26 cancers as well as pan-cancer data for identification of genetic and non-genetic drivers are available [here.](hits/)


# Shiny app
A web portal using the Shiny framework in R has been developed to visualize rank profiles of the DRIVE shRNA screen and corresponding gene expression data from TCGA at https://apsic.shinyapps.io/APSiC/. 

### Contributions
- [Hesam Montazeri](http://lcbb.ut.ac.ir/)
- [Salvatore Piscuoglio](http://oncogenomicslab.org/lab-members/)
- [Charlotte K Y Ng](http://oncogenomicslab.org/lab-members/)

### Contact
```
hesam.montazeri (at) ut.ac.ir
charlotte.ng (at) dbmr.unibe.ch
```
