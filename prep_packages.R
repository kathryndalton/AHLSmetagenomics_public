
## Packages Needed for AHLS Data Analysis ##

#########################################################################
#### INSTALL ####

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
BiocManager::install("phyloseq")
devtools::install_github("jbisanz/qiime2R")
BiocManager::install("decontam")
BiocManager::install("DESeq2")
install.packages("vegan")
devtools::install_github('twbattaglia/btools')
install.packages("plyr")
install.packages("tidyverse")
devtools::install_github("tidymodels/tidymodels")
install.packages("reshape2")
devtools::install_github("wilkelab/cowplot")
install.packages("gridExtra")
install.packages("cluster")
install.packages("qwraps2")
install.packages("knitr")
install.packages("stargazer")
install.packages("compareGroups")
devtools::install_github("rapporter/pander")
devtools::install_github("dcomtois/summarytools", build_vignettes = TRUE)
install.packages("gtsummary")
install.packages("GGally")
install.packages("car")
install.packages("Hmisc")
install.packages("foreign") # for STATA and SPSS files
install.packages("haven") # for SAS files
    


#########################################################################
#### LOAD ####

library(phyloseq)
library("qiime2R")
library(decontam)
library("DESeq2")
library("vegan")
library(btools)
library(plyr) ## need to upload before dpylr 
library("tidyverse")
# library(tidyr) ## part of tidyverse
# library("dplyr") ## part of tidyverse
# library(readr) ## part of tidyverse
# library("ggplot2") ## part of tidyverse
# library(tibble) ## part of tidyverse
# library(readr)  ## part of tidyverse 
library(tidymodels) # says error but I think it works
# library(broom)  ## part of tidymodels - also re-installs tidyverse packages
library("reshape2")
library("cowplot")
library(gridExtra)
library("cluster")
library(qwraps2)
library(knitr)
library(stargazer)
library(compareGroups)
library("pander")
library(summarytools)
library(gtsummary)
library(GGally)
library(car)
library(Hmisc)
library(foreign)
library(haven)
##occasionally dplyr acts up with plyr attached - need to detach(package:plyr) 

