load(file = 'Phenotype_selected.Rdata')
load(file = 'riskscore_validation.Rdata')
load(file = 'Survival_selected.Rdata')
Phenotype_validation=Phenotype_selected[names(riskscore),]
Cox_factor=data.frame(Phenotype_validation[,c('demographic.ethnicity','demographic.gender','demographic.race','diagnoses.age_at_diagnosis')])
Cox_factor[1:5,]
Cox_factor[,5]=riskscore_multiple
Cox_factor[,5]=riskscore
colnames(Cox_factor)[5]='riskscore'
Survival_validation=Survival_selected[rownames(Cox_factor),]
Cox_factor[,c(6,7)]=Survival_validation[,c(2,4)]
library(survival)
library(survminer)
Cox_factor[1:5]
Cox_factor[1:5,]
Cox_factor=Cox_factor[,c(2:7)]
Cox_factor[which(Cox_factor[,1]=='female'),1]=0
Cox_factor[which(Cox_factor[,1]=='male'),1]=1
Cox_factor[which(Cox_factor[,2]=='white'),2]=0
Cox_factor[which(Cox_factor[,2]=='black or african american'),2]=1
Cox_factor[which(Cox_factor[,2]=='asian'),2]=2
Cox_factor[which(Cox_factor[,2]=='other'),2]=3
Cox_factor[which(Cox_factor[,2]=='no reported'),2]=NA
Cox_factor[which(Cox_factor[,2]=='not allowed to collect'),2]=NA
Cox_factor[which(Cox_factor[,2]=='not reported'),2]=NA
cox=coxph(Surv(OS.time,OS)~.,data = Cox_factor)
Cox_factor[,1]=as.numeric(Cox_factor[,1])
Cox_factor[,2]=as.numeric(Cox_factor[,2])
cox=coxph(Surv(OS.time,OS)~.,data = Cox_factor)
summary(cox)
Cox_factor1=as.data.frame(cbind(Cox_factor[,c(1:3)],Phenotype_validation[,7],Cox_factor[,c(4:6)]))
Cox_factor1[which(Cox_factor1[,4]=='not reported'),4]=NA
Cox_factor1[which(Cox_factor1[,4]=='not hispanic or latino'),4]=1
Cox_factor1[which(Cox_factor1[,4]=='hispanic or latino'),4]=2
Cox_factor1[,1]=Cox_factor1[,1]+1
Cox_factor1[,2]=Cox_factor1[,2]+1
colnames(Cox_factor1)[4]='ethnicity'
Cox_factor2=Cox_factor1
riskscore_normalize=(riskscore-min(riskscore))/(max(riskscore)-min(riskscore))
Cox_factor2[,5]=riskscore_normalize
age=Cox_factor2[,3]/365
Cox_factor2[,3]=age
cox2=coxph(Surv(OS.time,OS)~.,data = Cox_factor2)
summary(cox2)