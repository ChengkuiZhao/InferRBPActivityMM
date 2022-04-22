library("survival")
library("survminer")
load(file = 'TPM_selected.Rdata')
load(file = 'TPM.Rdata')
load(file='status_validation.Rdata')
names_validation=names(status)
TPM_validation=TPM_selected[,names_validation]
TPM_validation=as.matrix(TPM_validation)
rownames(TPM_validation)=rownames(TPM)
load(file = 'hgnc_symbol.Rdata')
TPM_validation_hgncsymbol=TPM_validation[hgnc_symbol[,1],]
rownames(TPM_validation_hgncsymbol)=hgnc_symbol[,2]
status=data.frame(c(1:162),status)
library('ArrayTools')
output.cls(status,'status',filename = 'phenotype')
write.table(TPM_validation_hgncsymbol,'TPM_validation_hgncsymbol.txt',sep='\t',row.names = T,col.names = T)
