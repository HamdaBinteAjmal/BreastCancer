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
CJ  = c("CLDN1","CLDN3", "CLDN4", "CLDN7", "OCLN", "TJP1", "F11R", "CGN", "CXADR", "MARVELD3","CRB3", "PARD3", "SCRIB", "SFN", 
        "CDH1","CTNNB1", "CTNNA1", "JUP", "DSG2", "DSG3", "GJA1","PECAM1", "EPCAM","NCAM1", "ICAM1","VCAM1", "SELE", "ITGB3", "ITGB1", "ITGAV" )

CA20_xai <- molData$xai[paste0("xai", CA20), -8]
CJ30_xai <- molData$xai[paste0("xai", CJ), -8]

rownames(CA20_xai) <- CA20
rownames(CJ30_xai) <- CJ

CA20_exp <- molData$exp[paste0("exp", CA20),]
CJ30_exp <- molData$exp[paste0("exp", CJ),]

rownames(CA20_exp) <- CA20
rownames(CJ30_exp) <- CJ


cr <- crossCors(CA20_xai , CJ30_xai)
corrplot(abs(cr$cor), cl.lim = c(0,1), p.mat = cr$pval, sig.level = 0.1, insig =  "blank",  tl.col = "black", tl.srt = 45, type = "full",  title = "Sum of CA20 expression values vs. CJ30 genes")



cr <- crossCors(CA20_exp , CJ30_exp)
corrplot(abs(cr$cor), cl.lim = c(0,1), p.mat = cr$pval, sig.level = 0.1, insig =  "blank",  tl.col = "black", tl.srt = 45, type = "full",  title = "Sum of CA20 expression values vs. CJ30 genes")
corrplot(cr$cor, cl.lim = c(-1,1), p.mat = cr$pval, sig.level = 0.1, insig =  "blank",  tl.col = "black", tl.srt = 45, type = "full",  title = "Sum of CA20 expression values vs. CJ30 genes")


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

CA20_xai_sums <- colSums(CA20_xai)
CA20Score <- CA20_xai_sums
cr <- crossCors(CA20_xai_sums , CJ30_xai)
corrplot(cr$cor, cl.lim = c(-1,1), p.mat = cr$pval, sig.level = 0.1, insig =  "blank",  tl.col = "black", tl.srt = 45, type = "full",  title = "Sum of CA20 expression values vs. CJ30 genes")
CJ30_Score <- CJ30_xai


CA20_exp_sums <- colSums(CA20_exp)
CA20Score <- CA20_exp_sums
cr <- crossCors(CA20_exp_sums , CJ30_exp)
corrplot(cr$cor, cl.lim = c(-1,1), p.mat = cr$pval, sig.level = 0.1, insig =  "blank",  tl.col = "black", tl.srt = 45, type = "full",  title = "Sum of CA20 expression values vs. CJ30 genes")
CJ30_Score <- CJ30_exp

## CA20 and CJ30 ind.. scores
# Correlation between CA20 sums and CJ genes expression
means <- rowMeans(CJ30_exp, na.rm = TRUE)
for (i in 1:nrow(CJ30_Score))
{
  
  mean <- means[i]
  print(mean)
  for (j in 1:ncol(CJ30_Score))
  {
    if (!is.na(CJ30_Score[i,j]))
    {
      if(CJ30_Score[i,j] > mean)
      {
        CJ30_Score[i,j] <- 1
      }
      else
      {
        CJ30_Score[i,j] <- 0
      }
    }
    else
    {
      CJ30_Score[i,j] <- NA
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
cr1 <- crossCors(Y = CJ30_Score, X = CA20Score, method = "spearman")
corrplot(cr1$cor, cl.lim = c(-1,1), p.mat = cr1$pval, sig.level = 0.1, insig =  "blank",  tl.col = "black", tl.srt = 45, type = "full",  title = "Sum of CA20 expression values vs. CJ30 genes")


## Correlation between CA20 score and CJ30  Score

CJ30Score_full <- colMeans(CJ30_exp, na.rm = TRUE)
CJ30Score_full <- ifelse(CJ30Score_full > mean(CJ30Score_full), 1 , 0)
cor.test(CJ30Score_full, CA20Score, method = "kendall")