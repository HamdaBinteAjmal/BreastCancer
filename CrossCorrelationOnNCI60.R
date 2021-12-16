if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("rcellminer")
BiocManager::install("rcellminerData")

library(rcellminer)
library(rcellminerData)
library(corrplot)
drugAct <- exprs(getAct(rcellminerData::drugData))
molData <- getMolDataMatrices()

#### Starting OVer ## Measure Cross Correlation 
CA20 <- c("AURKA", "CCNA2","CCND1", "CCNE2", "CDK1",  "CEP63", "CEP152", "E2F1", "E2F2", "LMO4" 
               , "MDM2", "MYCN", "NDRG1", "NEK2", "PIN1", "PLK1","PLK4","SASS6", "STIL","TUBG1")
CJ30  = c("CLDN1","CLDN3", "CLDN4", "CLDN7", "OCLN", "TJP1", "F11R", "CGN", "CXADR", "MARVELD3","CRB3", "PARD3", "SCRIB", "SFN", 
        "CDH1","CTNNB1", "CTNNA1", "JUP", "DSG2", "DSC3", "GJA1","PECAM1", "EPCAM","NCAM1", "ICAM1","VCAM1", "SELE", "ITGB3", "ITGB1", "ITGAV" )

CA20_xai <- molData$xai[paste0("xai", CA20),-8]

CJ30_xai <- molData$xai[paste0("xai", CJ30), -8 ]
CA20_MedianCentred <- t(CA20_xai)
colnames(CA20_MedianCentred) <- CA20

library(matrixStats)
CA20_MedianCentred <- as.data.frame(scale(as.matrix(CA20_MedianCentred), center = colMedians(as.matrix(CA20_MedianCentred)), scale = F))
CA20Sums <- rowSums(as.matrix(CA20_MedianCentred))
CA20Mean <- mean(CA20Sums)
CA20Score <- ifelse(CA20Sums < CA20Mean, "Lower", "Higher")


## Perhaps also median centre the CJ30 data (should not really matter)
CJ30_MedianCentred <- t(CJ30_xai)
colnames(CJ30_MedianCentred) <- CJ30
CJ30_MedianCentred <- as.data.frame(scale(as.matrix(CJ30_MedianCentred), center = colMedians(as.matrix(CJ30_MedianCentred)), scale = F))


cr <- crossCors(t(CA20_MedianCentred) , t(CJ30_MedianCentred))
corrplot(cr$cor, cl.lim = c(-1,1), p.mat = cr$pval, sig.level = 0.1, insig =  "blank",  tl.col = "black", tl.srt = 45, type = "full",  title = "Sum of CA20 expression values vs. CJ30 genes")
corrplot(cr$cor, cl.lim = c(-1,1), p.mat = cr$pval, sig.level = 0.1, insig =  "blank",  tl.col = "black", tl.srt = 90, type = "full", mar = c(2, 2, 1, 0))
mtext(text = "CJ30 Panel", side = 3, at=15, line=2, cex=1.25, las =1 )
mtext(text = "CA20 Panel", side = 2, line = 1.5, at = 12 , cex = 1.25)

## Separate positive and negative 
positives <- which(cr$cor > 0)
negatives <- which(cr$cor <=0)

crPos <- cr
crPos$cor[negatives] <- 0
crPos$pval[negatives] <- 1
corrplot(crPos$cor, cl.lim = c(0,1), p.mat = crPos$pval, sig.level = 0.1, insig =  "blank",  tl.col = "black", tl.srt = 45, type = "full",  title = "Sum of CA20 expression values vs. CJ30 genes")

crNeg <- cr
crNeg$cor[positives] <-0
crNeg$pval[positives] <-1
corrplot(crNeg$cor, cl.lim = c(-1,0), p.mat = crNeg$pval, sig.level = 0.1, insig =  "blank",  tl.col = "black", tl.srt = 45, type = "full",  title = "Sum of CA20 expression values vs. CJ30 genes")


## CA20 Sum and CJ30 indv. Pos and Negative

cr <- crossCors(CA20Sums , t(CJ30_MedianCentred))
corrplot(cr$cor, cl.lim = c(-1,1), p.mat = cr$pval, sig.level = 0.1, insig =  "blank",  tl.col = "black", tl.srt = 45, type = "full",  title = "Sum of CA20 expression values vs. CJ30 genes")
dimnames(cr$cor)[[1]] = "CA20 Score"
dimnames(cr$pval)[[1]] = "CA20 Score"
rownames(cr$cor) = c("CA20 Score")
rownames(cr$pval) = c("CA20S core")

corrplot(cr$cor, cl.lim = c(-1,1), p.mat = cr$pval, sig.level = 0.05, insig =  "blank",  tl.col = "black", tl.srt = 90, type = "full", cl.pos = "n",mar = c(2, 1, 1, 1), xaxs="i", yaxs="i")
mtext(text = "CJ30 Panel", side = 3, at=15, line=-5, cex=1.25, las =1 )


CJ30_Score <- CJ30_xai

CA20_exp <- molData$exp[paste0("exp", CA20),-8]
CJ30_exp <- molData$exp[paste0("exp", CJ30),-8]

CA20_exp_sums <- colSums(CA20_exp)
CA20Score <- CA20_exp_sums
cr <- crossCors(CA20_exp_sums , CJ30_exp)
corrplot(cr$cor, cl.lim = c(-1,1), p.mat = cr$pval, sig.level = 0.1, insig =  "blank",  tl.col = "black", tl.srt = 45, type = "full",  title = "Sum of CA20 expression values vs. CJ30 genes")
CJ30_Score <- CJ30_MedianCentred

## CA20 and CJ30 ind.. scores
# Correlation between CA20 sums and CJ genes expression
means <- colMeans(CJ30_MedianCentred, na.rm = TRUE)
for (i in 1:ncol(CJ30_Score))
{
  
  mean <- means[i]
  print(mean)
  for (j in 1:nrow(CJ30_Score))
  {
    if (!is.na(CJ30_Score[j,i]))
    {
      if(CJ30_Score[j,i] > mean)
      {
        CJ30_Score[j,i] <- 1
      }
      else
      {
        CJ30_Score[j,i] <- 0
      }
    }
    else
    {
      CJ30_Score[j,i] <- NA
    }
  }
}

highs <- which(CA20Score > mean(CA20Score))
lows <- which(CA20Score < mean(CA20Score))
CA20Score[highs] = 1
CA20Score[lows] = 0
cr1<- apply(CJ30_Score, 1, function(x) chisq.test(x,CA20Score))



## Instead of ChiSq check spearmans monotonic 

#cr1<- apply(CJ30_Score, 1, function(x) cor.test(x = x, y = CA20Score, method ="spearman"))
CA20Score <- ifelse(CA20Score == "Higher", 1 , 0)
cr1 <- crossCors(Y = t(CJ30_Score), X = CA20Score, method = "spearman")
corrplot(cr1$cor, cl.lim = c(-1,1), p.mat = cr1$pval, sig.level = 0.1, insig =  "blank",  tl.col = "black", tl.srt = 45, type = "full",  title = "Sum of CA20 expression values vs. CJ30 genes")

corrplot(cr1$cor, cl.lim = c(-1,1), p.mat = cr$pval, sig.level = 0.05, insig =  "blank",  tl.col = "black", tl.srt = 90, type = "full", cl.pos = "n",mar = c(2, 1, 1, 1), xaxs="i", yaxs="i")
mtext(text = "CJ30 Panel", side = 3, at=15, line=-5, cex=1.25, las =1 )

## Correlation between CA20 score and CJ30  Score
CJ30Sums <- rowSums(as.matrix(CJ30_MedianCentred))
CJ30Mean <- mean(CJ30Sums)
CJ30Score <- ifelse(CJ30Sums < CJ30Mean, 0, 1)
cor.test(CJ30Score, CA20Score, method = "kendall")

## Linear Regression for all CJ30Genes and CA20Sums
data <- as.data.frame(cbind(CJ30_MedianCentred, CA20Sums))
model <- lm(CA20Sums ~., data = data)

##logistic regreesion for all CJ30 genes and CA20Score
data <- as.data.frame(cbind(CJ30_MedianCentred, CA20Score))
model <- glm(CA20Score ~.,family=binomial(link='logit'),data=data)
anova(model, test = "LRT")


