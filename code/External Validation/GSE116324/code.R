library(ggplot2)
load(file = 'PSImatrix_GSE116324.Rdata')
load(file = 'Min.Rdata')
load(file = 'D.Rdata')
load(file = 'Median_norm.Rdata')
load(file = 'M_0.3.Rdata')
###normalize the target
PSI_norm=apply(PSImatrix,2,as.numeric)
rownames(PSI_norm)=c(1:49076)
for (i in c(1:49076)){
  p=PSI_norm[i,]
  p1=p[which(p>=0)]
  p1_norm=((p1-Min[i])/D[i])
  PSI_norm[i,which(p>=0)]=p1_norm
  print(i)
  
}
###Get target score 
PSIscore=PSI_norm
for (i in c(1:49076)){
  s=PSIscore[i,]
  s1=s[!is.na(s)]
  m=Median_norm[i]
  PSIscore[i,!is.na(PSIscore[i,])]=(s1-m)
  print(i)
}
###Calculate RBP activity
RBP_Activity=matrix(0,150,44)
rownames(RBP_Activity)=rownames(M_0.3)
colnames(RBP_Activity)=colnames(PSIscore)
for (i in c(1:150)){
  for (j in c(1:44)){
    p=intersect(which(M_0.3[i,]==1),which(PSIscore[,j]>-2))
    n=intersect(which(M_0.3[i,]==(-1)),which(PSIscore[,j]>-2))
    RBP_Activity[i,j]=(sum(PSIscore[p,j])-sum(PSIscore[n,j]))/(length(p)+length(n))
  }
  print(i)
} 
###Based on the minimum meanErr to select the best hyper-parameter: ntree=50,mtry=15,nodesize=5
###Train with the whole 762 MMRF data; Validate with 44 MM samples in GSE116324 and do boxplot.
load(file='RBP_train_Activity.Rdata')
load(file='Survival_GSE116324.Rdata')
load(file = 'Survival_selected.Rdata')
library('randomForestSRC')
data_train=data.frame(t(RBP_train_Activity),Survival_selected[,c(2,4)])
data_validation=data.frame(t(RBP_Activity))
model=rfsrc(Surv(OS.time, OS) ~ .,data_train,na.action = 'na.impute',
            ntree=50,mtry=10,nodesize=5)
predict=predict.rfsrc(model,newdata=data_validation,na.action = 'na.impute')
riskscore=predict$predicted
names(riskscore)=rownames(data_validation)
rownames(clinic)=clinic[,1]
Stage_selected_validation=as.character(clinic[names(riskscore),4])
names(Stage_selected_validation)=names(riskscore)
###Boxplot
boxplot(riskscore)
riskscore_selected=riskscore[-c(1,13,22,36)]
stage=Stage_selected_validation[-c(1,13,22,36)]
riskscore_I=riskscore_selected[which(stage=='I')]
riskscore_II=riskscore_selected[which(stage=='II')]
riskscore_III=riskscore_selected[which(stage=='III')]
t.test(c(riskscore_I,riskscore_II),riskscore_III,alternative = 'less')
stage[which(stage=='I')]='Stage I'
stage[which(stage=='II')]='Stage II'
stage[which(stage=='III')]='Stage III'
GSE116324=data.frame(stage,riskscore_selected)
p1=ggplot(GSE116324,aes(x=stage,y=riskscore_selected,fill=stage))+geom_boxplot(outlier.shape=NA)
p1=p1+ggtitle('Plot of Different Stage Riskscore in GSE116324')
p1=p1+theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(limits=c('Stage I','Stage II','Stage III'))
GSE116324[-which(stage=='Stage III'),1]='Before Stage III'
p2=ggplot(GSE116324,aes(x=stage,y=riskscore_selected,fill=stage))+geom_boxplot(outlier.shape=NA)
p2=p2+ggtitle('Plot of Different Stage Riskscore in GSE116324\nt-test P value =0.046')
p2=p2+theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(limits=c('Before Stage III','Stage III'))
