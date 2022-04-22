library(ggplot2)
load(file = 'PSImatrix_GSE110486.Rdata')
load(file = 'Min.Rdata')
load(file = 'D.Rdata')
load(file = 'Median_norm.Rdata')
load(file = 'M_0.3.Rdata')
load(file = 'RBP_GSE110486.Rdata')
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
RBP_Activity=matrix(0,150,17)
rownames(RBP_Activity)=rownames(M_0.3)
colnames(RBP_Activity)=colnames(PSIscore)
for (i in c(1:150)){
  for (j in c(1:17)){
    p=intersect(which(M_0.3[i,]==1),which(PSIscore[,j]>-2))
    n=intersect(which(M_0.3[i,]==(-1)),which(PSIscore[,j]>-2))
    RBP_Activity[i,j]=(sum(PSIscore[p,j])-sum(PSIscore[n,j]))/(length(p)+length(n))
  }
  print(i)
} 
###See RBP Activity&Gene correlation
Correlation=rep(0,1,150)
Pvalue=rep(0,1,150)
for (i in c(1:150)){
  if(sum(RBP_Activity[i,],na.rm=T)!=0){
    Correlation[i]=cor.test(RBP_Activity[i,],RBP_GSE110486[i,])$estimate
    Pvalue[i]=cor.test(RBP_Activity[i,],RBP_GSE110486[i,])$p.value
  }
  else{Correlation[i]=NA
  Pvalue[i]=NA}
}
###Based on minimum meanErr to selec best hyper-parameter: ntree=50,mtry=10,nodesize=5
###Train with the whole 762 training samples; Validate with 17 MM samples in GSE110486 and do boxplot.
load(file='RBP_train_Activity.Rdata')
load(file='Survival_GSE110486.Rdata')
load(file = 'Median_norm.Rdata')
load(file = 'Min.Rdata')
load(file = 'D.Rdata')
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
Stage_selected_validation=as.character(clinic[names(riskscore),2])
names(Stage_selected_validation)=names(riskscore)
stage=Stage_selected_validation
stage[6:17]='MM'
GSE110486=data.frame(stage,riskscore)
p=ggplot(GSE110486,aes(x=stage,y=riskscore,fill=stage))+geom_boxplot()
p=p+ggtitle('Plot of Different Stage Riskscore in GSE110486')
p=p+theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(limits=c('Normal','MM'))+guides(fill = guide_legend(reverse=TRUE))
