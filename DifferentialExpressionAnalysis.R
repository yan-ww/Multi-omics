# 三种方法：limma,DESeq2,EdgeR
#limma
#limma算法原理：核心在于lmfit()和eBays()，参考https://www.jianshu.com/p/716cb81d1f61
library(limma)
library(dplyr)
count<-read.csv("count-FvN.csv",row.names=1)
# count的行名是基因名，列名是样本名
condition = c(rep("F",3),rep("N",3)) %>% factor(.,levels = c("F","N"),ordered = F) 
design = model.matrix(~0+condition)
colnames(design) = c("F","N")
df.fit = lmFit(count,design) #数据与list进行匹配
df.matrix = makeContrasts(F-N,levels=design)
fit = contrasts.fit(df.fit,df.matrix)
fit = eBayes(fit)
res = topTable(fit,number=Inf)
res = res[order(res$P.Value),]
