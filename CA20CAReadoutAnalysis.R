## This is proper calculation of CA20Score ##
CA20_MedianCentred <- t(CA20_xai)
colnames(CA20_MedianCentred) <- CA20

library(matrixStats)
CA20_MedianCentred <- as.data.frame(scale(as.matrix(CA20_MedianCentred), center = colMedians(as.matrix(CA20_MedianCentred)), scale = F))
CA20Sums <- rowSums(as.matrix(CA20_MedianCentred))
CA20Mean <- mean(CA20Sums)
CA20Score_nci <- ifelse(CA20Sums < CA20Mean, 0, 1)
#########################################


CACellLines <- xlsx::read.xlsx2("C:\\LispHome\\OneDrive - National University of Ireland, Galway\\BreastCancer\\NCI60CA.xlsx",1,header = TRUE, as.data.frame = TRUE, colIndex = c(1, 3))
names <- CACellLines[,1]
CACellLines <- CACellLines[,2]
names(CACellLines) <- names
CACellLines <- CACellLines[names(CA20Sums)]
## Find corrlation between CA20Score and CA from biological readout
cor.test(CACellLines, CA20Sums)


CACellLine_bi <- ifelse(CACellLines > 13,1,0)
cor.test(CA20Score_nci , CACellLine_bi, method = "kendall")

CACellLine_bi <- ifelse(CACellLines > 13,1,0)
cor.test(CA20Score_nci , CACellLine_bi, method = "spearman")

#Correlation between expression of each individual gene in 
#CA20 and percentage of CA (from biological readout)

cr <- crossCors(CACellLines, t(CA20_MedianCentred))
corrplot(cr$cor, cl.lim = c(-1,1), p.mat = cr$pval, sig.level = 0.1, insig =  "blank",  tl.col = "black", tl.srt = 45, type = "full",  title = "Sum of CA20 expression values vs. CJ30 genes")


## Logistic Regression
data <- as.data.frame(cbind(CA20_MedianCentred, CACellLine_bi))
model <- glm(CACellLine_bi ~.,family=binomial(link='logit'),data=data)
anova(model, test = "LRT")

a <- CACellLine_bi[which(CA20Score_nci == 1)]
b <- CA20Score_nci[which(CA20Score_nci == 1)]
