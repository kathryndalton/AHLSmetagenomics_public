
## Step 2 do alpha diversity descriptive tests and models ##

## do phase2_step1_import first - end with wgs_data (with updated variables) & sd_wgs_data (data.frame not phyloseq)

###############################################################################################

# Get Alpha Metrics - Shannon, Richness/Observed and Faiths PD

plot_richness(wgs_data, x="State", measures=c("Observed", "Shannon"))
plot_richness(wgs_data, x="Season", color="CropAnimalFarming", measures=c("Observed", "Shannon")) 
# no obvious trends at least with these

#prep - need to convert sd_wgs_data rownames to column
sd_wgs_data_analysis<-rownames_to_column(sd_wgs_data, var="SampleID_seq")

alpha<-estimate_richness(wgs_data, measures=c("Observed", "Shannon"))%>%
  rownames_to_column(., var="SampleID")
## get Faiths
fpd<-estimate_pd(wgs_data) %>%
  rownames_to_column(., var="SampleID") %>%
  mutate(Faiths = PD)
alpha<-merge(alpha, fpd, by="SampleID") %>%
  select(., SampleID, Observed, Shannon, Faiths) %>%
  merge(., sd_wgs_data_analysis, by.x="SampleID", by.y="SampleID_seq")


###############################################################################################

## Basic EDA of Shannon and Richness ##

# main covariates = Male, Age, State, Farmer, Season, HomeCondition, Carpeting, LiveFarm, DogsOrCats, Dogs, Cats, 
# CropAnimalFarming, BeefCattle, DairyCattle, Hogs, Poultry 

hist(alpha$Observed)
hist(alpha$Shannon) # very left skewed - Ziyue transformed it
hist(alpha$Faiths)

ggplot(alpha, aes(y=Observed, x=Male)) + geom_boxplot() + geom_point(aes(color=Age), position = "jitter")
ggplot(alpha, aes(y=Observed, x=State, color=Farmer)) + geom_boxplot() # maybe diff within farmers by state 
ggplot(alpha, aes(y=Observed, x=Season)) + geom_boxplot()
ggplot(alpha, aes(y=Observed, x=HomeCondition, color=Carpeting)) + geom_boxplot() ## may be diff in homecondition (esp without carpeting)
ggplot(alpha, aes(y=Observed, x=LiveFarm, color=DogsOrCats)) + geom_boxplot() # live on farm higher
ggplot(alpha, aes(y=Shannon, x=Male)) + geom_boxplot() + geom_point(aes(color=Age), position = "jitter")
ggplot(alpha, aes(y=Shannon, x=State, color=Farmer)) + geom_boxplot()
ggplot(alpha, aes(y=Shannon, x=Season)) + geom_boxplot()
ggplot(alpha, aes(y=Shannon, x=HomeCondition, color=Carpeting)) + geom_boxplot()
ggplot(alpha, aes(y=Shannon, x=LiveFarm, color=DogsOrCats)) + geom_boxplot() # live on farm higher
ggplot(alpha, aes(y=Faiths, x=Male)) + geom_boxplot() + geom_point(aes(color=Age), position = "jitter")
ggplot(alpha, aes(y=Faiths, x=State, color=Farmer)) + geom_boxplot() # very similar to Observed
ggplot(alpha, aes(y=Faiths, x=Season)) + geom_boxplot()
ggplot(alpha, aes(y=Faiths, x=HomeCondition, color=Carpeting)) + geom_boxplot() # very similar to Observed
ggplot(alpha, aes(y=Faiths, x=LiveFarm, color=DogsOrCats)) + geom_boxplot() # very similar to Observed
## harder to see differences in Shannon since so skewed ##

## descriptive KW test

compareGroups(Male ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable() 
#methods=2 means non-normal distribution, present median [IQR] and p-value from Kruskall-Wallis test
# other option method=1 normal distribution, would show mean (sd) and p-value for t-test
compareGroups(State ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable() ##Shannon sig p=0.04
compareGroups(Farmer ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable() 
compareGroups(Season ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable()
compareGroups(HomeCondition ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable()  ##Shannon sig p=0.02
compareGroups(Carpeting ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable() 
compareGroups(LiveFarm ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable()  ##Shannon sig p<0.001
compareGroups(DogsOrCats ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable()  ##Shannon sig p<0.001
compareGroups(Dogs ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable()
compareGroups(Cats ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable() ##Shannon sig p=0.02
compareGroups(CropAnimalFarming ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable() ##Shannon sig p<0.001, Obs & FPD also sig p=0.03 !!
compareGroups(BeefCattle ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable()  ##Shannon sig p<0.001 (Obs/FPD p=0.07)
compareGroups(DairyCattle ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable()
compareGroups(Hogs ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable()  ##Shannon sig p<0.001
compareGroups(Poultry ~ Shannon + Observed + Faiths, alpha, method=2) %>% createTable()  ##Shannon sig p=0.001
# need to do Age as correlation
cor.test(alpha$Age, alpha$Observed, method="spearman")
cor.test(alpha$Age, alpha$Shannon, method="spearman") # p=0.01
cor.test(alpha$Age, alpha$Faiths, method="spearman")


###############################################################################################

## Univariate Modeling of Richness/Observed, Shannon and Faiths ##

# main covariates = Male, Age, State, Farmer, Season, HomeCondition, Carpeting, LiveFarm, DogsOrCats, Dogs, Cats, 
# CropAnimalFarming, BeefCattle, DairyCattle, Hogs, Poultry 

tidy_ulm_obs<-rbind(
  tidy(glm(Observed ~ Male, data=alpha)),
  tidy(glm(Observed ~ Age, data=alpha)),
  tidy(glm(Observed ~ State, data=alpha)),
  tidy(glm(Observed ~ Farmer, data=alpha)),
  tidy(glm(Observed ~ Season, data=alpha)),
  tidy(glm(Observed ~ HomeCondition, data=alpha)),
  tidy(glm(Observed ~ Carpeting, data=alpha)),
  tidy(glm(Observed ~ LiveFarm, data=alpha)),
  tidy(glm(Observed ~ DogsOrCats, data=alpha)),
  tidy(glm(Observed ~ Dogs, data=alpha)),
  tidy(glm(Observed ~ Cats, data=alpha)),
  tidy(glm(Observed ~ CropAnimalFarming, data=alpha)),
  tidy(glm(Observed ~ BeefCattle, data=alpha)),
  tidy(glm(Observed ~ DairyCattle, data=alpha)),
  tidy(glm(Observed ~ Hogs, data=alpha)),
  tidy(glm(Observed ~ Poultry, data=alpha))) %>%
  subset(term!="(Intercept)")
## sig vars (p<0.05) = HomeCondition (Poor vs OK, p=0.03), LiveFarm (Yes vs No, p=0.01), 
#                      BeefCattle (Yes vs No, p=0.05), Hogs (Yes vs No, p=0.005)

tidy_ulm_shan<-rbind(
  tidy(glm(Shannon ~ Male, data=alpha)),
  tidy(glm(Shannon ~ Age, data=alpha)),
  tidy(glm(Shannon ~ State, data=alpha)),
  tidy(glm(Shannon ~ Farmer, data=alpha)),
  tidy(glm(Shannon ~ Season, data=alpha)),
  tidy(glm(Shannon ~ HomeCondition, data=alpha)),
  tidy(glm(Shannon ~ Carpeting, data=alpha)),
  tidy(glm(Shannon ~ LiveFarm, data=alpha)),
  tidy(glm(Shannon ~ DogsOrCats, data=alpha)),
  tidy(glm(Shannon ~ Dogs, data=alpha)),
  tidy(glm(Shannon ~ Cats, data=alpha)),
  tidy(glm(Shannon ~ CropAnimalFarming, data=alpha)),
  tidy(glm(Shannon ~ BeefCattle, data=alpha)),
  tidy(glm(Shannon ~ DairyCattle, data=alpha)),
  tidy(glm(Shannon ~ Hogs, data=alpha)),
  tidy(glm(Shannon ~ Poultry, data=alpha))) %>%
  subset(term!="(Intercept)")
## sig vars (p<0.05) = HomeCondition (Poor vs OK, p=0.05), LiveFarm (Yes vs No, p<0.001), DogsOrCats (Yes vs No, p=0.007),
#                      Cats (Yes vs No, p=0.01), CropAnimalFarming (NoCrop.YesAnimals vs No/No p=0.04, YesCrop.YesAnimals vs No/No p<0.001),
#                      BeefCattle (Yes vs No, p<0.001), Hogs (Yes vs No, p<0.001), Poultry (Yes vs No, p=0.001)

tidy_ulm_fpd<-rbind(
  tidy(glm(Faiths ~ Male, data=alpha)),
  tidy(glm(Faiths ~ Age, data=alpha)),
  tidy(glm(Faiths ~ State, data=alpha)),
  tidy(glm(Faiths ~ Farmer, data=alpha)),
  tidy(glm(Faiths ~ Season, data=alpha)),
  tidy(glm(Faiths ~ HomeCondition, data=alpha)),
  tidy(glm(Faiths ~ Carpeting, data=alpha)),
  tidy(glm(Faiths ~ LiveFarm, data=alpha)),
  tidy(glm(Faiths ~ DogsOrCats, data=alpha)),
  tidy(glm(Faiths ~ Dogs, data=alpha)),
  tidy(glm(Faiths ~ Cats, data=alpha)),
  tidy(glm(Faiths ~ CropAnimalFarming, data=alpha)),
  tidy(glm(Faiths ~ BeefCattle, data=alpha)),
  tidy(glm(Faiths ~ DairyCattle, data=alpha)),
  tidy(glm(Faiths ~ Hogs, data=alpha)),
  tidy(glm(Faiths ~ Poultry, data=alpha))) %>%
  subset(term!="(Intercept)")
## sig vars (p<0.05) = LiveFarm (Yes vs No, p=0.01), Hogs (Yes vs No, p=0.007)


###############################################################################################

## Multivariate Modeling of Richness/Observed, Shannon and Faiths ##


## testing model development

# m1 - check basic vars
m1<-glm(Observed ~ State + Farmer + Male + Age + Season + DogsOrCats + HomeCondition + Carpeting + LiveFarm + CropAnimalFarming, data=alpha)
car::vif(m1)  ##farmer & Male VIF >10
# m1A - simple, w/o Male
m1A<-glm(Observed ~ State + Farmer + Age + Season + DogsOrCats + HomeCondition + Carpeting + LiveFarm + CropAnimalFarming, data=alpha) ##farmer & Male VIF >10
# m1B - simple, w/o Farmer
m1B<-glm(Observed ~ State + Male + Age + Season + DogsOrCats + HomeCondition + Carpeting + LiveFarm + CropAnimalFarming, data=alpha) ##farmer & Male VIF >10
car::vif(m1A) # now all VIFs below 1
car::vif(m1B) 
## check AIC in stargazer
stargazer(m1, m1A, m1B, type="text")
# no difference in AIC = no difference in predictive power 
## In MK's paper - did not include Farmer as covariate (just Male) --- KEEP AS THAT MODEL == m1B
# Interesting when removing either Farmer or Male (m1A or m1B) - LiveFarm comes up as significant

## decide between cropanimalfarming vs cropfarming, animalfarming
m1C<-glm(Observed ~ State + Farmer + Male + Age + Season + DogsOrCats + HomeCondition + Carpeting + LiveFarm + CropFarming + AnimalFarming, data=alpha)
car::vif(m1C)
stargazer(m1, m1C, type="text") -# little difference in AIC - keep cropanimalfarming

m1C## testing variables assumed to have collinearity (dogsorcats / dogs / cat and cropanimalfarming / diff animal production)
# m2 - changing DogsOrCats to Dogs+Cats
m2<-glm(Observed ~ State + Male + Age + Season + Dogs + Cats + HomeCondition + Carpeting + LiveFarm + CropAnimalFarming, data=alpha)
# m3 - changing cropanimal to types of animals
m3<-glm(Observed ~ State + Male + Age + Season + DogsOrCats + HomeCondition + Carpeting + LiveFarm + BeefCattle + DairyCattle + Hogs + Poultry, data=alpha)
# m4 - changing BOTH DogsOrCats to Dogs+Cats & cropanimal to types of animals
m4<-glm(Observed ~ State + Male + Age + Season + Dogs + Cats + HomeCondition + Carpeting + LiveFarm + BeefCattle + DairyCattle + Hogs + Poultry, data=alpha)
# m5 - including both DogsOrCats AND Dogs+Cats
m5<-glm(Observed ~ State + Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition + Carpeting + LiveFarm + CropAnimalFarming, data=alpha) 
# m6 - including both cropanimal AND types of animals
m6<-glm(Observed ~ State + Male + Age + Season + DogsOrCats + HomeCondition + Carpeting + LiveFarm + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry, data=alpha)
# m7 - including ALL variables DogsOrCats AND Dogs+Cats + cropanimal AND types of animals
m7<-glm(Observed ~ State + Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition + Carpeting + LiveFarm + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry, data=alpha)
car::vif(m2) ## all vifs 1
car::vif(m3) ## all vifs 1
car::vif(m4) ## all vifs 1
car::vif(m5) ## dogsorcats = 6, dogs = 4, cats=2
car::vif(m6) ## cropanimal = 4. beef = 2, rest animals = 1
car::vif(m7) ## same as above 
stargazer(m1B, m2, m3, m4, m5, m6, m7, type="text", column.labels = c("m1", "m2", "m3", "m4", "m5", "m6", "m7"), model.numbers=F)

## double checking cropfarm & animalfarm vs cropanimalfarm
m8<-glm(Observed ~ State + Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition + Carpeting + LiveFarm + CropFarming + AnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry, data=alpha)
car::vif(m8) #Animal farming now 3, beef = 2
stargazer(m3, m7, m8, type="text", column.labels = c("m3", "m7", "m8"), model.numbers=F) ## AIC for m7 and m8 ~equal, so keep as cropanimal

## based on these results --
# LiveFarm sig in all models
# Hog strong sig in models that include it (not m2, m5)
# HomeCondition is only sig in models without individ animal production adjusted (JUST m1, m2, m5) !!
# YesCrop/YesAnimal is only sig in m1 - goes away when separating out dogs and cats, and including ind animal production

## AIC values = m3 < m1 < m4 < m2 < m5 < m6 < m7, range 13,281-13,289 --- not that much of a change

## USE MODEL 7 (all vars)
##################################

##### final models ####

mlm_obs<- glm(Observed ~ State + Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition +
                          Carpeting + LiveFarm + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry, data=alpha)
tidy_mlm_obs<- tidy(mlm_obs) %>% subset(term!="(Intercept)")
mlm_shan<- glm(Shannon ~ State + Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition +
                           Carpeting + LiveFarm + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry, data=alpha)
tidy_mlm_shan<- tidy(mlm_shan) %>% subset(term!="(Intercept)")
mlm_fpd<- glm(Faiths ~ State + Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition +
                          Carpeting + LiveFarm + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry, data=alpha)
tidy_mlm_fpd<- tidy(mlm_fpd) %>% subset(term!="(Intercept)")
stargazer(mlm_obs, mlm_shan, mlm_fpd, type="text", column.labels = c("Observed", "Shannon", "Faiths"), model.numbers=F)

### LiveFarm sig all three metrics (positive increase diversity with yes farm live)
### Male sig neg assoc with Shannon
### YesCropYesAnimal sig pos assoc with Shannon
### Hogs sig pos assoc with Obs & Faiths


###############################################################################################
###############################################################################################

## Overview of Results:
## Key consistent variables associated with alpha diversity seem to be: 
# Living on a Farm most consistent -- Yes to living on a farm ~~ higher alpha diversity across different metrics
# Crop and Animal production - yes to both compared to neither, associated with higher diversity (and potentially Animal > Crop)
# Hog seems to be importantly associated with higher diversity
# other potentially important variables = home condition, dogs/cats (cats > dogs)


## left to do:

##### log transformed Shannon ####
## add in asthma
## interaction of state with all variables or only potentially important / relevant ones? --- did she end up doing an interaction model or stratified models?
## compare results of Ziyue's models with mine

# Questions for ZW
# transformation "exponential" vs log - ln or log10?
# interaction vs stratified
# whats your reference for season?


###############################################################################################
###############################################################################################


## Adding in Asthma to Models ##

# Rerun final multivariate models with asthma ##
mlm_as_obs<- glm(Observed ~ State + Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition +
                Carpeting + LiveFarm + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry + Asthma, data=alpha)
tidy_mlm_as_obs<- tidy(mlm_as_obs) %>% subset(term!="(Intercept)")
mlm_as_shan<- glm(Shannon ~ State + Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition +
                 Carpeting + LiveFarm + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry + Asthma, data=alpha)
tidy_mlm_as_shan<- tidy(mlm_as_shan) %>% subset(term!="(Intercept)")
mlm_as_fpd<- glm(Faiths ~ State + Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition +
                Carpeting + LiveFarm + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry + Asthma, data=alpha)
tidy_mlm_as_fpd<- tidy(mlm_as_fpd) %>% subset(term!="(Intercept)")
stargazer(mlm_as_obs, mlm_as_shan, mlm_as_fpd, type="text", column.labels = c("Observed", "Shannon", "Faiths"), model.numbers=F)

## Asthma itself not sig, controlling for it decrease sig for Hogs in Shannon and increased Shann assoc cropanimalfarm, male, and livefarm, rest not changed


###############################################################################################

## State Stratified Models ##

alpha_nc <- alpha %>% filter(., State=="NC")
alpha_ia <- alpha %>% filter(., State=="IA")

# Rerun final multivariate models within state stratified datasets ##
mlm_nc_obs<- glm(Observed ~ Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition +
                   Carpeting + LiveFarm + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry + Asthma, data=alpha_nc)
tidy_mlm_nc_obs<- tidy(mlm_nc_obs) %>% subset(term!="(Intercept)")
mlm_nc_shan<- glm(Shannon ~ Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition +
                    Carpeting + LiveFarm + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry + Asthma, data=alpha_nc)
tidy_mlm_nc_shan<- tidy(mlm_nc_shan) %>% subset(term!="(Intercept)")
mlm_nc_fpd<- glm(Faiths ~ Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition +
                   Carpeting + LiveFarm + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry + Asthma, data=alpha_nc)
tidy_mlm_nc_fpd<- tidy(mlm_nc_fpd) %>% subset(term!="(Intercept)")

mlm_ia_obs<- glm(Observed ~ Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition +
                      Carpeting + LiveFarm + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry + Asthma, data=alpha_ia)
tidy_mlm_ia_obs<- tidy(mlm_ia_obs) %>% subset(term!="(Intercept)")
mlm_ia_shan<- glm(Shannon ~ Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition +
                       Carpeting + LiveFarm + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry + Asthma, data=alpha_ia)
tidy_mlm_ia_shan<- tidy(mlm_ia_shan) %>% subset(term!="(Intercept)")
mlm_ia_fpd<- glm(Faiths ~ Male + Age + Season + DogsOrCats + Dogs + Cats + HomeCondition +
                      Carpeting + LiveFarm + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry + Asthma, data=alpha_ia)
tidy_mlm_ia_fpd<- tidy(mlm_ia_fpd) %>% subset(term!="(Intercept)")

stargazer(mlm_as_obs, mlm_as_shan, mlm_as_fpd,
          mlm_nc_obs, mlm_nc_shan, mlm_nc_fpd,
          mlm_ia_obs, mlm_ia_shan, mlm_ia_fpd, type="text",
          column.labels = c("Observed", "Shannon", "Faiths",
                            "Observed NC", "Shannon NC", "Faiths NC", 
                            "Observed IA", "Shannon IA", "Faiths IA"), model.numbers=F)

## Biggest difference - LiveFarm only sig in IA, not in NC (with stronger effect)
## HomeCondition sig in NC not IA (Obs/FPD), 
## also effect male & yescrop.yesanimal only in IA

## Plan - test for interaction with LiveFarm and HomeCondition

mlm_int_obs<- glm(Observed ~ State * HomeCondition * LiveFarm + Male + Age + Season + DogsOrCats + Dogs + Cats +
                   Carpeting + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry + Asthma, data=alpha)
tidy_mlm_int_obs<- tidy(mlm_int_obs) %>% subset(term!="(Intercept)")
mlm_int_shan<- glm(Shannon ~ State * HomeCondition * LiveFarm + Male + Age + Season + DogsOrCats + Dogs + Cats +
                     Carpeting + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry + Asthma, data=alpha)
tidy_mlm_int_shan<- tidy(mlm_int_shan) %>% subset(term!="(Intercept)")
mlm_int_fpd<- glm(Faiths ~ State * HomeCondition * LiveFarm + Male + Age + Season + DogsOrCats + Dogs + Cats +
                    Carpeting + CropAnimalFarming + BeefCattle + DairyCattle + Hogs + Poultry + Asthma, data=alpha)
tidy_mlm_int_fpd<- tidy(mlm_int_fpd) %>% subset(term!="(Intercept)")
stargazer(mlm_as_obs, mlm_as_shan, mlm_as_fpd,
          mlm_int_obs, mlm_int_shan, mlm_int_fpd, type="text", column.labels = c("Observed", "Shannon", "Faiths",
                                                                                 "Observed Int", "Shannon Int", "Faiths Int"), model.numbers=F)
## none of the interaction terms sig ##
## LiveFarm not sig for any metrics in interaction models! 
## But when removing interaction term for *Home Condition, still sig (State*LiveFarm not sig) - would use this more so


