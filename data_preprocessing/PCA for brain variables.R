#########=============extracting statistically significant brain variables=================

mor<-read.csv(file='./mor_merge.csv')
CP<-read.csv(file='./pvalue_corrected_cov_CP.csv')
EA<-read.csv(file='./pvalue_corrected_cov_EA.csv')
subject_list<-read.csv(file='./GPS_phenotype.csv')
mor$subjectkey<-gsub("_","",mor$subjectkey)
subject_list<-subject_list[,1:2]
mor<-merge(subject_list,mor,by='subjectkey',all.x=TRUE)
mor<-mor[,-2]

CP_extract<-CP[CP$BONF<0.05,]
CP_list<-CP_extract$glm.result.brain
CP_mor<-mor[,CP_list]
subjectkey<-mor[,1]
CP_mor<-cbind(subjectkey,CP_mor)

EA_extract<-EA[EA$BONF<0.05,]
EA_list<-EA_extract$glm.result.brain
EA_mor<-mor[,EA_list]
subjectkey<-mor[,1]
EA_mor<-cbind(subjectkey,EA_mor)

#save
write.csv(CP_mor,file='./CP_mor_extracted_bonf.csv',row.names=FALSE,quote=F)
write.csv(EA_mor,file='./EA_mor_extracted_bonf.csv',row.names=FALSE,quote=F)


######By performing PCA, reduce statistically significant brain variables. And assign Principal Component(reduced brain variable) to subjects 
CP_mor_filled<-CP_mor
subjectkey<-CP_mor[,1]
for (i in 2:ncol(CP_mor)){
  CP_mor_filled[,i][is.na(CP_mor_filled[,i])]<-mean(CP_mor_filled[,i],na.rm=T)
} ##fill na with mean of each variables
CP_mor_filled_pca<-CP_mor_filled[,-1]  ##remove subjectkey column for PCA
CP_pca<-prcomp(CP_mor_filled_pca,center=T,scale. = T)
CP_pca_var<-as.matrix(CP_mor_filled_pca)%*%CP_pca$rotation[,1:16]#check cumulative Proportion to select how many PCs you need
CP_PC_transformed<-cbind(subjectkey,as.data.frame(CP_pca_var))

EA_mor_filled<-EA_mor
subjectkey<-EA_mor[,1]
for (i in 2:ncol(EA_mor)){
  EA_mor_filled[,i][is.na(EA_mor_filled[,i])]<-mean(EA_mor_filled[,i],na.rm=T)
} ##fill na with mean of each variables
EA_mor_filled_pca<-EA_mor_filled[,-1]  ##remove subjectkey column for PCA
EA_pca<-prcomp(EA_mor_filled_pca,center=T,scale. = T)
EA_pca_var<-as.matrix(EA_mor_filled_pca)%*%EA_pca$rotation[,1:4]#check cumulative Proportion to select how many PCs you need
EA_PC_transformed<-cbind(subjectkey,as.data.frame(EA_pca_var))

##save
write.csv(CP_PC_transformed,file='./CP_PCA_transformed.csv',row.names=FALSE,quote=F)
write.csv(EA_PC_transformed,file='./EA_PCA_transformed.csv',row.names=FALSE,quote=F)



##======PCA plotting 
library(ggplot2)
library(ggfortify)
cov_data<-read.csv(file='./GPS_phenotype.csv')
list<-c('subjectkey','age', 'sex', 'race.ethnicity', 'high.educ', 'income', 'abcd_site', 'height', 'weight', 'BMI')
cov_list<-c('age', 'sex', 'race.ethnicity', 'high.educ', 'income', 'abcd_site', 'height', 'weight', 'BMI')
cov_data<-cov_data[,list]

#for CP
autoplot(CP_pca,data=cov_data,colour='age',main='age')
autoplot(CP_pca,data=cov_data,colour='sex',main='sex')
autoplot(CP_pca,data=cov_data,colour='race.ethnicity',main='race.ethnicity')
autoplot(CP_pca,data=cov_data,colour='high.educ',main='high.educ')
autoplot(CP_pca,data=cov_data,colour='income',main='income')
autoplot(CP_pca,data=cov_data,colour='abcd_site',main='abcd_site')
autoplot(CP_pca,data=cov_data,colour='height',main='height')
autoplot(CP_pca,data=cov_data,colour='weight',main='weight')
autoplot(CP_pca,data=cov_data,colour='BMI',main='BMI')

dev.off()

#for EA
autoplot(EA_pca,data=cov_data,colour='age',main='age')
autoplot(EA_pca,data=cov_data,colour='sex',main='sex')
autoplot(EA_pca,data=cov_data,colour='race.ethnicity',main='race.ethnicity')
autoplot(EA_pca,data=cov_data,colour='high.educ',main='high.educ')
autoplot(EA_pca,data=cov_data,colour='income',main='income')
autoplot(EA_pca,data=cov_data,colour='abcd_site',main='abcd_site')
autoplot(EA_pca,data=cov_data,colour='height',main='height')
autoplot(EA_pca,data=cov_data,colour='weight',main='weight')
autoplot(EA_pca,data=cov_data,colour='BMI',main='BMI')

#########=============extracting statistically significant brain variables with sex classification=================
#=========CP============
####about CP for male
library(dplyr)
# mor<-read.csv(file='./mor_merge.csv')
# CP<-read.csv(file='./pvalue_corrected_cov_CP.csv')
# EA<-read.csv(file='./pvalue_corrected_cov_EA.csv')
subject_list<-read.csv(file='./GPS_phenotype.csv')
mor$subjectkey<-gsub("_","",mor$subjectkey)
subject_list<-subject_list[,c('subjectkey','sex')]
mor<-merge(subject_list,mor,by='subjectkey',all.x=TRUE)

male_mor<-mor[mor$sex==1,]
male_mor<-male_mor[!is.na(male_mor$sex),]
male_mor<-male_mor[,-2]

CP_extract<-CP[CP$BONF<0.05,]
CP_list<-CP_extract$glm.result.brain
CP_mor_male<-subset(male_mor,select=c('subjectkey',CP_list))

CP_mor_filled_male<-CP_mor_male
subjectkey<-CP_mor_male[,1]
for (i in 2:ncol(CP_mor_male)){
  CP_mor_filled_male[,i][is.na(CP_mor_filled_male[,i])]<-mean(CP_mor_filled_male[,i],na.rm=T)
} ##fill na with mean of each variables
CP_mor_filled_pca_male<-CP_mor_filled_male[,-1]  ##remove subjectkey column for PCA
CP_pca_male<-prcomp(CP_mor_filled_pca_male,center=T,scale. = T)
CP_pca_var_male<-as.matrix(CP_mor_filled_pca_male)%*%CP_pca_male$rotation[,1:16]#check cumulative Proportion to select how many PCs you need
CP_PC_transformed_male<-cbind(subjectkey,as.data.frame(CP_pca_var_male))

#plotting for male
#library(ggplot2)
#library(ggfortify)
cov_data<-read.csv(file='./GPS_phenotype.csv')
list<-c('subjectkey','age', 'sex', 'race.ethnicity', 'high.educ', 'income', 'abcd_site', 'height', 'weight', 'BMI')
cov_list<-c('age', 'sex', 'race.ethnicity', 'high.educ', 'income', 'abcd_site', 'height', 'weight', 'BMI')
cov_data<-cov_data[,list]
PC_data<-merge(cov_data,CP_PC_transformed_male,by='subjectkey')

autoplot(CP_pca_male,data=PC_data,colour='age',main='age')
autoplot(CP_pca_male,data=PC_data,colour='race.ethnicity',main='race.ethnicity')
autoplot(CP_pca_male,data=PC_data,colour='high.educ',main='high.educ')
autoplot(CP_pca_male,data=PC_data,colour='income',main='income')
autoplot(CP_pca_male,data=PC_data,colour='abcd_site',main='abcd_site')
autoplot(CP_pca_male,data=PC_data,colour='height',main='height')
autoplot(CP_pca_male,data=PC_data,colour='weight',main='weight')
autoplot(CP_pca_male,data=PC_data,colour='BMI',main='BMI')



####About CP for female
#library(dplyr)
mor<-read.csv(file='./mor_merge.csv')
CP<-read.csv(file='./pvalue_corrected_cov_CP.csv')
EA<-read.csv(file='./pvalue_corrected_cov_EA.csv')
subject_list<-read.csv(file='./GPS_phenotype.csv')
mor$subjectkey<-gsub("_","",mor$subjectkey)
subject_list<-subject_list[,c('subjectkey','sex')]
mor<-merge(subject_list,mor,by='subjectkey',all.x=TRUE)
female_mor<-mor[mor$sex==2,]
female_mor<-female_mor[!is.na(female_mor$sex),]
female_mor<-female_mor[,-2]



CP_extract<-CP[CP$BONF<0.05,]
CP_list<-CP_extract$glm.result.brain
CP_mor_female<-subset(female_mor,select=c('subjectkey',CP_list))

CP_mor_filled_female<-CP_mor_female
subjectkey<-CP_mor_female[,1]
for (i in 2:ncol(CP_mor_female)){
  CP_mor_filled_female[,i][is.na(CP_mor_filled_female[,i])]<-mean(CP_mor_filled_female[,i],na.rm=T)
} ##fill na with mean of each variables
CP_mor_filled_pca_female<-CP_mor_filled_female[,-1]  ##remove subjectkey column for PCA
CP_pca_female<-prcomp(CP_mor_filled_pca_female,center=T,scale. = T)
CP_pca_var_female<-as.matrix(CP_mor_filled_pca_female)%*%CP_pca_female$rotation[,1:16]#check cumulative Proportion to select how many PCs you need
CP_PC_transformed_female<-cbind(subjectkey,as.data.frame(CP_pca_var_female))

#plotting for female
library(ggplot2)
library(ggfortify)
cov_data<-read.csv(file='./GPS_phenotype.csv')
list<-c('subjectkey','age', 'sex', 'race.ethnicity', 'high.educ', 'income', 'abcd_site', 'height', 'weight', 'BMI')
cov_list<-c('age', 'sex', 'race.ethnicity', 'high.educ', 'income', 'abcd_site', 'height', 'weight', 'BMI')
cov_data<-cov_data[,list]
PC_data<-merge(cov_data,CP_PC_transformed_female,by='subjectkey')

autoplot(CP_pca_female,data=PC_data,colour='age',main='age',loadings=T)
autoplot(CP_pca_female,data=PC_data,colour='race.ethnicity',main='race.ethnicity')
autoplot(CP_pca_female,data=PC_data,colour='high.educ',main='high.educ')
autoplot(CP_pca_female,data=PC_data,colour='income',main='income')
autoplot(CP_pca_female,data=PC_data,colour='abcd_site',main='abcd_site')
autoplot(CP_pca_female,data=PC_data,colour='height',main='height')
autoplot(CP_pca_female,data=PC_data,colour='weight',main='weight')
autoplot(CP_pca_female,data=PC_data,colour='BMI',main='BMI')


##save
write.csv(CP_PC_transformed_male,file='./CP_PCA_transformed_male.csv',row.names=FALSE,quote=F)
write.csv(CP_PC_transformed_female,file='./CP_PCA_transformed_female.csv',row.names=FALSE,quote=F)



#=========EA============
####about EA for male
#library(dplyr)
mor<-read.csv(file='./mor_merge.csv')
CP<-read.csv(file='./pvalue_corrected_cov_CP.csv')
EA<-read.csv(file='./pvalue_corrected_cov_EA.csv')
subject_list<-read.csv(file='./GPS_phenotype.csv')
mor$subjectkey<-gsub("_","",mor$subjectkey)
subject_list<-subject_list[,c('subjectkey','sex')]
mor<-merge(subject_list,mor,by='subjectkey',all.x=TRUE)

male_mor<-mor[mor$sex==1,]
male_mor<-male_mor[!is.na(male_mor$sex),]
male_mor<-male_mor[,-2]

EA_extract<-EA[EA$BONF<0.05,]
EA_list<-EA_extract$glm.result.brain
EA_mor_male<-subset(male_mor,select=c('subjectkey',EA_list))

EA_mor_filled_male<-EA_mor_male
subjectkey<-EA_mor_male[,1]
for (i in 2:ncol(EA_mor)){
  EA_mor_filled_male[,i][is.na(EA_mor_filled_male[,i])]<-mean(EA_mor_filled_male[,i],na.rm=T)
} ##fill na with mean of each variables
EA_mor_filled_pca_male<-EA_mor_filled_male[,-1]  ##remove subjectkey column for PCA
EA_pca_male<-prcomp(EA_mor_filled_pca_male,center=T,scale. = T)
EA_pca_var_male<-as.matrix(EA_mor_filled_pca_male)%*%EA_pca_male$rotation[,1:16]#check cumulative Proportion to select how many PCs you need
EA_PC_transformed_male<-cbind(subjectkey,as.data.frame(EA_pca_var_male))

#plotting for male
#library(ggplot2)
#library(ggfortify)
cov_data<-read.csv(file='./GPS_phenotype.csv')
list<-c('subjectkey','age', 'sex', 'race.ethnicity', 'high.educ', 'income', 'abcd_site', 'height', 'weight', 'BMI')
cov_list<-c('age', 'sex', 'race.ethnicity', 'high.educ', 'income', 'abcd_site', 'height', 'weight', 'BMI')
cov_data<-cov_data[,list]
PC_data<-merge(cov_data,EA_PC_transformed_male,by='subjectkey')

autoplot(EA_pca_male,data=PC_data,colour='age',main='age')
autoplot(EA_pca_male,data=PC_data,colour='race.ethnicity',main='race.ethnicity')
autoplot(EA_pca_male,data=PC_data,colour='high.educ',main='high.educ')
autoplot(EA_pca_male,data=PC_data,colour='income',main='income')
autoplot(EA_pca_male,data=PC_data,colour='abcd_site',main='abcd_site')
autoplot(EA_pca_male,data=PC_data,colour='height',main='height')
autoplot(EA_pca_male,data=PC_data,colour='weight',main='weight')
autoplot(EA_pca_male,data=PC_data,colour='BMI',main='BMI')



####About EA for female
#library(dplyr)
mor<-read.csv(file='./mor_merge.csv')
CP<-read.csv(file='./pvalue_corrected_cov_CP.csv')
EA<-read.csv(file='./pvalue_corrected_cov_EA.csv')
subject_list<-read.csv(file='./GPS_phenotype.csv')
mor$subjectkey<-gsub("_","",mor$subjectkey)
subject_list<-subject_list[,c('subjectkey','sex')]
mor<-merge(subject_list,mor,by='subjectkey',all.x=TRUE)
female_mor<-mor[mor$sex==2,]
female_mor<-female_mor[!is.na(female_mor$sex),]
female_mor<-female_mor[,-2]



EA_extract<-EA[EA$BONF<0.05,]
EA_list<-EA_extract$glm.result.brain
EA_mor_female<-subset(female_mor,select=c('subjectkey',EA_list))

EA_mor_filled_female<-EA_mor_female
subjectkey<-EA_mor_female[,1]
for (i in 2:ncol(EA_mor_female)){
  EA_mor_filled_female[,i][is.na(EA_mor_filled_female[,i])]<-mean(EA_mor_filled_female[,i],na.rm=T)
} ##fill na with mean of each variables
EA_mor_filled_pca_female<-EA_mor_filled_female[,-1]  ##remove subjectkey column for PCA
EA_pca_female<-prcomp(EA_mor_filled_pca_female,center=T,scale. = T)
EA_pca_var_female<-as.matrix(EA_mor_filled_pca_female)%*%EA_pca_female$rotation[,1:16]#check cumulative Proportion to select how many PCs you need
EA_PC_transformed_female<-cbind(subjectkey,as.data.frame(EA_pca_var_female))

#plotting for female
library(ggplot2)
library(ggfortify)
cov_data<-read.csv(file='./GPS_phenotype.csv')
list<-c('subjectkey','age', 'sex', 'race.ethnicity', 'high.educ', 'income', 'abcd_site', 'height', 'weight', 'BMI')
cov_list<-c('age', 'sex', 'race.ethnicity', 'high.educ', 'income', 'abcd_site', 'height', 'weight', 'BMI')
cov_data<-cov_data[,list]
PC_data<-merge(cov_data,EA_PC_transformed_female,by='subjectkey')

autoplot(EA_pca_female,data=PC_data,colour='age',main='age')
autoplot(EA_pca_female,data=PC_data,colour='race.ethnicity',main='race.ethnicity')
autoplot(EA_pca_female,data=PC_data,colour='high.educ',main='high.educ')
autoplot(EA_pca_female,data=PC_data,colour='income',main='income')
autoplot(EA_pca_female,data=PC_data,colour='abcd_site',main='abcd_site')
autoplot(EA_pca_female,data=PC_data,colour='height',main='height')
autoplot(EA_pca_female,data=PC_data,colour='weight',main='weight')
autoplot(EA_pca_female,data=PC_data,colour='BMI',main='BMI')


##save
write.csv(EA_PC_transformed_male,file='./EA_PCA_transformed_male.csv',row.names=FALSE,quote=F)
write.csv(EA_PC_transformed_female,file='./EA_PCA_transformed_female.csv',row.names=FALSE,quote=F)

GPS_phe <- read.csv(file='./GPS_phenotype.csv')
list <- c('subjectkey', 'EA', 'ELS_total', 'nihtbx_totalcomp_uncorrected', 'age')
inter <- GPS_phe[, list]
all_male <- merge(EA_PC_transformed_male, inter, by='subjectkey')
all_male$ELS_total[is.na(all_male$ELS_total)] <- mean(all_male$ELS_total, na.rm=T)
  

all_female <- merge(EA_PC_transformed_female, inter, by='subjectkey')
all_female$ELS_total[is.na(all_female$ELS_total)] <- mean(all_female$ELS_total, na.rm=T)


nih <- read.table('abcd_tbss01.txt', header=T)
nih <- nih[-1,]
nih$nihtbx_cryst_agecorrected
nih_list <- c('subjectkey', 'nihtbx_totalcomp_uncorrected','nihtbx_totalcomp_agecorrected', 'nihtbx_fluidcomp_uncorrected', 'nihtbx_fluidcomp_agecorrected', 'nihtbx_cryst_uncorrected', 'nihtbx_cryst_agecorrected')
nih <- nih[, nih_list]
nih$subjectkey<-gsub("_","",nih$subjectkey)

for (i in 2:ncol(nih)){
  nih[,i]<-as.numeric(nih[,i])
  nih[,i][is.na(nih[,i])]<-mean(nih[,i],na.rm=T)
  nih[,i]<-scale(nih[,i])
}
all_male<- merge(all_male, nih, by='subjectkey')
all_female<-merge(all_female, nih, by='subjectkey')
write.csv(all_male,file='./EA_PCA_GPS_phenotype_merged_male.csv',row.names=FALSE,quote=F)
write.csv(all_female,file='./EA_PCA_GPS_phenotype_merged_female.csv',row.names=FALSE,quote=F)

# Code for new model: ELS_total<->NIH_reading, Ab_PA<->NIH_total
