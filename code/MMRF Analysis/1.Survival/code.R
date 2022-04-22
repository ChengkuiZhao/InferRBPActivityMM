###Split Data
load(file = 'PSImatrix_selected.Rdata')
load(file = 'RBP_selected.Rdata')
train=sample(c(1:762),600)
validation=c(1:762)[-train]
RBP_train=as.matrix(RBP_selected[,train])
PSI_train=as.matrix(PSImatrix_selected[,train])
RBP_validation=as.matrix(RBP_selected[,validation])
PSI_validation=as.matrix(PSImatrix_selected[,validation])
###Find CLIP Targets
ENCODE=read.csv('CorrelationMatrix')
rownames(ENCODE)=ENCODE[,1]
ENCODE=t(as.matrix(ENCODE[,-1]))
load(file = 'Transcriptname.Rdata')
M_CLIP=matrix(0,150,49076)
rownames(M_CLIP)=rownames(ENCODE)
colnames(M_CLIP)=ENSGname
for(i in c(1:150)){
  encode=colnames(ENCODE)[which(ENCODE[i,]!=0)]
  for(j in c(1:49076)){
    if(length(intersect(ENSGname[j],encode))==1){
      M_CLIP[i,j]=1
    }
  }
  print(i)
}
Targetnumber=apply(M_CLIP,1,sum)
boxplot(Targetnumber,col = c('lightblue'),outline = F,xlab = 'RBPs',ylab='Targetsnumber',main='Target number of RBP in CLIP')
text(c(0.9),c(10000),labels = round(median(Targetnumber),0))
###Finall Targets Selected By Spearman Correlation
M=matrix(0,150,49076)
rownames(M)=rownames(M_CLIP)
colnames(M)=colnames(M_CLIP)
for (i in c(1:150)){
  index=which(M_CLIP[i,]==1)
  correlation=rep(0,length(index))
  for (j in 1:(length(index))){
    nonnoise=which(PSI_train[index[j],]>=0)
    correlation[j]=cor.test(RBP_train[i,nonnoise],PSI_train[index[j],nonnoise],method  = 'spearman')$estimate
  }
  M[i,index[which(correlation>=0.3)]]=1
  M[i,index[which(correlation<=(-0.3))]]=-1
  print(i)
}
s=apply(abs(M),1,sum)
boxplot(s,col = c('lightblue'),outline = F,xlab = 'RBPs',ylab='Targetsnumber',main='Target number of RBP Finally')
text(c(0.9),c(700),labels = round(median(s),0))
###Filter&Normalize&Calculate
###list all the range difference
D=rep(0,1,49076)
IQR=rep(0,1,49076)
for (i in c(1:49076)){
  l=as.numeric(PSI_train[i,])
  l=l[which(l>=0)]
  d=diff(range(l))
  iqr=IQR(l)
  D[i]=d
  IQR[i]=iqr
}
###Get target matrix
M_IQR=M
M_IQR[,which(IQR<summary(IQR)[2])]=0
###get min of PSI_train& normalize targets
Min=rep(0,1,49076)
PSI_train_norm=apply(PSI_train,2,as.numeric)
rownames(PSI_train_norm)=c(1:49076)
for (i in c(1:49076)){
  p=PSI_train_norm[i,]
  p1=p[which(p>=0)]
  p1_norm=((p1-min(p1,na.rm = T))/D[i])
  Min[i]=min(p1,na.rm=T)
  PSI_train_norm[i,which(p>=0)]=p1_norm
}
###Get target score
PSIscore_train=PSI_train_norm
Median=rep(0,49076)
for (i in c(1:49076)){
  s=PSIscore_train[i,]
  s1=s[!is.na(s)]
  m=median(s1)
  Median[i]=m
  PSIscore_train[i,which(!is.na(s))]=(s1-m)
}
###Calculate RBP activity
RBP_train_Activity=matrix(0,150,600)
rownames(RBP_train_Activity)=rownames(M_IQR)
colnames(RBP_train_Activity)=colnames(PSIscore_train)
for (i in c(1:150)){
  for (j in c(1:600)){
    p=intersect(which(M_IQR[i,]==1),which(PSIscore_train[,j]>-2))
    n=intersect(which(M_IQR[i,]==(-1)),which(PSIscore_train[,j]>-2))
    RBP_train_Activity[i,j]=(sum(PSIscore_train[p,j])-sum(PSIscore_train[n,j]))/(length(p)+length(n))
  }
  print(i)
} 
###Normalize Targets
PSI_validation_norm=apply(PSI_validation,2,as.numeric)
rownames(PSI_validation_norm)=c(1:49076)
for (i in c(1:49076)){
  q=PSI_validation_norm[i,]
  q1=q[which(q>=0)]
  q1_norm=((q1-Min[i])/D[i])
  PSI_validation_norm[i,which(q>=0)]=q1_norm
  print(i)
}
###Get validation target score 
PSIscore_validation=PSI_validation_norm
for (i in c(1:49076)){
  s=PSIscore_validation[i,]
  s1=s[!is.na(s)]
  m=Median[i]
  PSIscore_validation[i,!is.na(PSIscore_validation[i,])]=(s1-m)
  print(i)
}
###Calculate validation RBP activity
RBP_validation_Activity=matrix(0,150,162)
rownames(RBP_validation_Activity)=rownames(M_IQR)
colnames(RBP_validation_Activity)=colnames(PSIscore_validation)
for (i in c(1:150)){
  for (j in c(1:162)){
    p=intersect(which(M_IQR[i,]==1),which(PSIscore_validation[,j]>-2))
    n=intersect(which(M_IQR[i,]==(-1)),which(PSIscore_validation[,j]>-2))
    RBP_validation_Activity[i,j]=(sum(PSIscore_validation[p,j])-sum(PSIscore_validation[n,j]))/(length(p)+length(n))
  }
  print(i)
} 
######Conduct 5-fold cross-validation for training data
library('randomForestSRC')
load(file = 'Survival_selected.Rdata')
Survival_selected_train=Survival_selected[colnames(RBP_train_Activity),]
data_train=data.frame(t(RBP_train_Activity),Survival_selected_train[,c(2,4)])
a=sample(c(1:600))
a1=a[1:120]
a2=a[121:240]
a3=a[241:360]
a4=a[361:480]
a5=a[481:600]
X=rbind(a1,a2,a3,a4,a5)
ntree=c(50,100,200,500,1000)
mtry=c(10,13,15,20,30)
nodesize=c(5,10,20,50,100)
Err=array(0,dim = c(5,5,5,5))
for (i in c(1:5)){
  for (j in c(1:5)) {
    for (k in c(1:5)) {
      for (n in c(1:5)) {
        model=rfsrc(Surv(OS.time, OS) ~ .,data_train[-X[n,],],na.action = 'na.impute',
                    ntree=ntree[i],mtry=mtry[j],nodesize=nodesize[k])
        predict=predict.rfsrc(model,newdata=data_train[X[n,],],na.action = 'na.impute')
        Err[i,j,k,n]=predict$err.rate[ntree[i]]
        print(c(i,j,k,n))
      }
    }
  }
}
meanErr=array(0,dim = c(5,5,5))
for (i in c(1:5)){
  for (j in c(1:5)) {
    for (k in c(1:5)){
      meanErr[i,j,k]=mean(Err[i,j,k,])
    }
  }
}
###Based on the minimum meanErr to select the best hyper-parameter: ntree=50,mtry=15,nodesize=5
###Train with the whole 600 training samples; Vlidating with the 162 validation samples  Plot survival curve for the two predicted groups (High&Low risk group)
Survival_selected_validation=Survival_selected[colnames(RBP_validation_Activity),]
data_validation=data.frame(t(RBP_validation_Activity),Survival_selected_validation[,c(2,4)])
model=rfsrc(Surv(OS.time, OS) ~ .,data_train,na.action = 'na.impute',
            ntree=50,mtry=15,nodesize=5)
predict=predict.rfsrc(model,newdata=data_validation,na.action = 'na.impute')
###Plot survival by group
library("survival")
library("survminer")
riskscore=predict$predicted
names(riskscore)=rownames(Survival_selected_validation)
status=riskscore
status[which(riskscore>median(riskscore))]='HighRisk'
status[which(riskscore<=median(riskscore))]='LowRisk'
data_suv=cbind(Survival_selected_validation,riskscore,status)
surv_object=Surv(time=data_suv$OS.time,event = data_suv$OS)
fit=survfit(surv_object~status,data=data_suv)
ggsurvplot(fit,data=data_suv,pval=TRUE)
###Check ISS data
load(file='Phenotype_selected.Rdata')
ISSstage_selected_validation=Phenotype_selected[rownames(Survival_selected_validation),24]
names(ISSstage_selected_validation)=rownames(Survival_selected_validation)
stage1=riskscore[which(ISSstage_selected_validation=='I')]
stage2=riskscore[which(ISSstage_selected_validation=='II')]
stage3=riskscore[which(ISSstage_selected_validation=='III')]
stage12=c(stage1,stage2)
boxplot(stage12,stage3,names = c('Stage I&II','Stage III'),col = c('lightblue','green'),ylab='Riskscore',main='Riskscore in different ISS stage')
ttest=round(t.test(stage12,stage3,'less')$p.value,4)
text(c(0.8,1.8,1.5),c(27,40,45),labels = c(round(median(stage12),2),round(median(stage3),2),'p=0.003'))
###ggplot
riskscore_selected=riskscore[which(riskscore<40)]
stage=ISSstage_selected_validation[which(riskscore<40)]
stage[which(stage=='I')]='Stage I'
stage[which(stage=='II')]='Stage II'
stage[which(stage=='III')]='Stage III'
VALIDATION=data.frame(stage,riskscore_selected)
p1=ggplot(VALIDATION,aes(x=stage,y=riskscore_selected,fill=stage))+geom_boxplot(outlier.shape=NA)
p1=p1+ggtitle('Plot of Different Stage Riskscore in VALIDATION')
p1=p1+theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(limits=c('Stage I','Stage II','Stage III'))
riskscore_I=riskscore_selected[which(stage=='Stage I')]
riskscore_II=riskscore_selected[which(stage=='Stage II')]
riskscore_III=riskscore_selected[which(stage=='Stage III')]
VALIDATION[-which(stage=='Stage III'),1]='Before Stage III'
t.test(c(riskscore_I,riskscore_II),riskscore_III,alternative = 'less')
p2=ggplot(VALIDATION,aes(x=stage,y=riskscore_selected,fill=stage))+geom_boxplot(outlier.shape=NA)
p2=p2+ggtitle('Plot of Different Stage Riskscore in VALIDATION\nt-test P value =0.033')
p2=p2+theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(limits=c('Before Stage III','Stage III'))
p1
p2