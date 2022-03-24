## Step 1 import datasets, clean and explore variables ##

### Import ## 
## use same code from previous Phase 1 - Step 1 - Import
wgs_data <- readRDS("smb://wine/metagenome/Data/AGHL_taxa_abn_filter_1.rds")
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6528 taxa and 781 samples ]
#sample_data() Sample Data:       [ 781 samples by 49 sample variables ]
#tax_table()   Taxonomy Table:    [ 6528 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 6528 tips and 6527 internal nodes ]
ahls_data <- read_sas("https://sftp.s-3.net/#/LHS27001_BodyComp/DATA/lhs27001_master_22feb07.sas7bdat", 
                      NULL)

### Explore Variables ###

wgs_var<-as.data.frame(colnames(sample_data(wgs_data)))
## what is "orig_sample_id.y"? - guess is it's something that Ziyue created
ahls_var<-as.data.frame(colnames(ahls_data))
## many variables - need to use full data dictionary to determine important ones
## to include only MB participants - use "In_Final_Analysis_n_879" and  "Excluded_41"
table(ahls_data$In_Final_Analysis_n_879)
table(ahls_data$Excluded_41) ## 879-41 = 838 **plus more excluded during bioinformatics QC

## Edit some Variable Labels ##
sd_wgs_data<-as.matrix(sample_data(wgs_data)) %>%
  as.data.frame() %>%
  rownames_to_column(., var="otu") %>%
  mutate(State = factor(State, levels = c("NO", "YES"), labels = c("IA", "NC")),
         Farmer = factor(Farmer, levels = c("NO", "YES")),
         Male = factor(Male, levels = c("NO", "YES")),
         Age = as.numeric(Age),
         Season = factor(Season, levels = c("Spring", "Summer", "Fall", "Winter")),
         DogsOrCats = factor(DogsOrCats, levels = c("NO", "YES")),
         HomeCondition = factor(HomeCondition, levels=c("NO", "YES", ""), labels = c("OK", "Poor", "NA")),
         Carpeting = factor(Carpeting, levels=c("N", "Y"), labels = c("NO", "YES")),
         LiveFarm = factor(LiveFarm, levels=c("NO", "YES")),
         CropAnimalFarming = factor(CropAnimalFarming, levels=c("NoCrop.NoAnimals","NoCrop.YesAnimals",
                                                                "YesCrop.NoAnimals","YesCrop.YesAnimals")))
## In 16S dataset - have Asthma variable (yes/no) --- pull in 
data_16s_asthma<-read.delim(
  "smb://wine/metagenome/Data/ALHS_Factors_Asthma_Metagenome_MKL.txt",
  colClasses = "character") %>%
  select(., c("SampleID", "Asthma"))
## match by SampleID (asthma) and SampleID.y
sd_wgs_data<-merge(sd_wgs_data, data_16s_asthma, by="SampleID", all.x=T) %>%
  column_to_rownames(., var="otu")

sd_wgs_data_ps<- phyloseq(sample_data(sd_wgs_data))
#sd_wgs_data<-phyloseq(sample_data(sd_wgs_data))
sample_data(wgs_data)<-sd_wgs_data_ps

###################################################################

### Simple EDA / Table 1 for Key Variables ###

key_var<-sd_wgs_data[, c("Male", "Age", "State", "Farmer", "Season", "HomeCondition", "Carpeting", "LiveFarm", "DogsOrCats",
                         "Dogs", "Cats", "CropAnimalFarming", "BeefCattle", "DairyCattle", "Hogs", "Poultry", "Asthma")]
dfSummary(key_var)

## Skewed groups - Carpeting only 6% Yes, DairyCattle 6% Yes 

