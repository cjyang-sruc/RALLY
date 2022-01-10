library(AlphaSimR)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(GWASpoly)
library(ggrepel)

### set working directory.
setwd()

### LL function for normal distribution pdf.
LL.norm <- function(p, x){
  a <- 1/(p[2]*sqrt(2*pi))
  L1 <- exp(-(x-p[1])^2/(2*p[2]^2))
  L2 <- exp(-(-x-p[1])^2/(2*p[2]^2))
  LL <- -( sum(log(L1/2 + L2/2)) + length(x)*log(a) )
  return(LL)
}


### function to center-jitter the values on x-axis.
center.jitter <- function(x, y, width=0.5){
  y2 <- unique(y)
  z <- rep(0, length(y))
  for(i in y2){
    temp <- which(y == i)
    z[temp] <- seq(-(length(temp)-1)/2, (length(temp)-1)/2, 1)
  }
  z <- z/max(z)*width
  return(x+z)
}


### we will simulate 32 founders.
### the first 10 generations will be derived from F2 + 4 SSD from random pairs of founders without selection.
### the QTLs will be set up such that the QTLs are found in 16 to 1 founder(s) - 32 QTLs total.
### run one with selection, another without selection.
### run RALLY and GWAS at the end.
### loop this for 100 times.

### create 10 founder varieties.
### set up the genetic map positions (markers evenly spaced at 0.1 cM).
glen <- seq(1.0, 2.8, 0.2)
gmap <- lapply(glen, FUN=function(x) seq(0, x-0.001, 0.001))

outT <- list()
for(h in 1:100){
  
  ### set up the founder haplotypes.
  founder <- runMacs(nInd=32,
                     nChr=10,
                     segSites=glen*1000,
                     inbred=TRUE,
                     species="GENERIC",
                     ploidy=2L,
                     manualGenLen=glen-0.001)
  SP <- SimParam$new(founder)
  founder <- newPop(founder, simParam=SP)
  fgeno <- pullSegSiteGeno(founder, simParam=SP)
  fac <- colSums(fgeno==2)
  
  ### identify the QTL positions.
  ### 32 QTLs total, 2 QTLs in 16, 15, ..., 1 founders.
  qtl <- vector()
  for(i in 1:16){
    temp <- which(fac==i)
    qtl <- c(qtl, sample(temp, 2, replace=FALSE))
  }
  
  ### set up the qtl effects.
  ### scale the qtl effects such that Vg=1 when the allele frequencies are 0.5.
  meff <- rep(0, ncol(fgeno))
  meff[qtl] <- c(rep(4, 2), rep(3, 4), rep(2, 8), rep(1, 18))
  meff <- meff/sd(c(fgeno%*%meff))
  
  ### create the first 10 generations of unselected populations.
  pop <- list()
  for(i in 1:10){
    pop[[i]] <- randCross(pop=founder,
                          nCrosses=32,
                          nProgeny=1,
                          simParam=SP)
    for(j in 1:5) pop[[i]] <- self(pop=pop[[i]], nProgeny=1, simParam=SP) # make F2 and 4 SSD.
  }
  
  ### run the population over 50 generations (with and without phenotypic selection).
  ### in each generation, we randomly select 16 pairs to make F2 with 4 SSD.
  ### these 16 pairs (32 parent varieties) are chosen from 6-10 years ago.
  popS <- pop # with selection.
  popU <- pop # without selection.
  
  ### create the genetic and phenotypic values for the populations in the first 10 years.
  GV <- PV <- matrix(0, nrow=65, ncol=32)
  for(i in 1:10) GV[i,] <- c(pullSegSiteGeno(popS[[i]], simParam=SP)%*%meff)
  temp <- rnorm(n=32*10, mean=0, sd=1)
  temp <- temp/sd(temp)
  PV[1:10,] <- GV[1:10,] + temp
  
  ### run the population from generation 11 to 65.
  for(i in 11:65){
    
    ### column g is the randomly chosen year (6 - 10 years ago).
    ### column v is the randomly chosen individual within year.
    xinfo <- cbind(g1=sample((i-10):(i-6), 16, replace=TRUE),
                   g2=sample((i-10):(i-6), 16, replace=TRUE),
                   v1=sample(1:32, 16, replace=TRUE),
                   v2=sample(1:32, 16, replace=TRUE))
    
    ### prevent selfing in the chosen parent pairs.
    temp <- which(xinfo[,1] == xinfo[,2] & xinfo[,3] == xinfo[,4])
    if(length(temp) > 0){
      for(j in temp){
        if(xinfo[j, 4] == 32) xinfo[j, 4] <- xinfo[j, 4] - 1 else xinfo[j, 4] <- xinfo[j, 4] + 1
      }
    }
    
    ### with selection.
    temp.pop <- list()
    for(j in 1:16){
      temp <- makeCross2(females=popS[[ xinfo[j,1] ]],
                         males=popS[[ xinfo[j,2] ]],
                         crossPlan=xinfo[j, 3:4, drop=FALSE],
                         nProgeny=1,
                         simParam=SP)
      temp <- self(pop=temp, nProgeny=100, simParam=SP) # 100 F2 individuals.
      for(k in 1:4) temp <- self(pop=temp, nProgeny=1, simParam=SP) # self each for 4 generations of SSD.
      
      temp.GV <- c(pullSegSiteGeno(temp, simParam=SP)%*%meff)
      temp.PV <- temp.GV + rnorm(n=100, mean=0, sd=1)
      temp.sel <- order(-temp.PV)[1:2] # select top 2 individuals.
      temp.pop[[j]] <- selectInd(pop=temp,
                                 nInd=2,
                                 use="rand",
                                 candidates=temp.sel,
                                 simParam=SP)
      temp.sel <- sapply(1:2, FUN=function(x) which(temp@iid == temp.pop[[j]]@iid[x]))
      GV[i, ((j-1)*2 + 1):(j*2)] <- temp.GV[temp.sel] # re-adjust because selectInd randomize the individuals.
      PV[i, ((j-1)*2 + 1):(j*2)] <- temp.PV[temp.sel] # re-adjust because selectInd randomize the individuals.
    }
    popS[[i]] <- mergePops(temp.pop)
    
    ### without selection.
    temp.pop <- list()
    for(j in 1:16){
      temp <- makeCross2(females=popU[[ xinfo[j,1] ]],
                         males=popU[[ xinfo[j,2] ]],
                         crossPlan=xinfo[j, 3:4, drop=FALSE],
                         nProgeny=1,
                         simParam=SP)
      temp <- self(pop=temp, nProgeny=2, simParam=SP) # 2 F2 individuals.
      for(k in 1:4) temp <- self(pop=temp, nProgeny=1, simParam=SP) # self each for 4 generations of SSD.
      temp.pop[[j]] <- temp
    }
    popU[[i]] <- mergePops(temp.pop)
    
  }
  

  ### create the marker name.
  marker <- data.frame(chr=unlist(lapply(1:10, FUN=function(x) rep(x, glen[x]*1000))),
                       pos=unlist(lapply(1:10, FUN=function(x) 0:(glen[x]*1000-1))))
  marker <- cbind(marker=paste("M", marker$chr, "_", formatC(x=marker$pos, width=4, flag="0"), sep=""), marker)
  marker$pos <- marker$pos/10 # convert it to cM.
  
  ### extract all 32 lines from each generation (16 to 65) - with selection.
  genoS <- vector()
  phenoS <- vector()
  for(i in 16:65){
    genoS <- rbind(genoS, pullSegSiteGeno(pop=popS[[i]], simParam=SP))
    phenoS <- rbind(phenoS, data.frame(year=i-15, GV=GV[i, ], PV=PV[i, ]))
  }
  phenoS <- cbind(ID=paste("S", 1:1600, sep=""), phenoS)
  phenoS$PV2 <- phenoS$GV + rnorm(n=1600, mean=0, sd=1) # reset Ve to 1 so it resembles a new trial with all varieties.
  rownames(genoS) <- phenoS$ID
  colnames(genoS) <- marker$marker

  ### extract all 32 lines from each generation (16 to 65) - without selection.
  genoU <- vector()
  phenoU <- vector()
  for(i in 16:65){
    genoU <- rbind(genoU, pullSegSiteGeno(pop=popU[[i]], simParam=SP))
    phenoU <- rbind(phenoU, data.frame(year=rep(i-15, 32)))
  }
  phenoU$GV <- c(genoU%*%meff)
  phenoU$PV <- phenoU$GV + rnorm(n=1600, mean=0, sd=1)
  phenoU <- cbind(ID=paste("U", 1:1600, sep=""), phenoU)
  rownames(genoU) <- phenoU$ID
  colnames(genoU) <- marker$marker
  
  ### set all the heterozygous markers to NA, and marker 2 to 1.
  genoS2 <- genoS
  genoS2[genoS2==1] <- NA
  genoS2 <- genoS2/2
  genoU2 <- genoU
  genoU2[genoU2==1] <- NA
  genoU2 <- genoU2/2
  
  ### subset 8 random lines from each generation.
  temp <- c(sapply(1:50, FUN=function(x) 32*(x-1) + sample(1:32, 8, replace=FALSE) ))
  genoS2 <- genoS2[temp,]
  phenoS <- phenoS[temp,]
  genoU2 <- genoU2[temp,]
  phenoU <- phenoU[temp,]
  
  # recode the marker so that 0 and 1 are equally balanced.
  temp <- which(sample(c(0,1), nrow(marker), replace=TRUE)==1)
  for(i in temp) genoS2[,i] <- 1 - genoS2[,i]
  temp <- which(sample(c(0,1), nrow(marker), replace=TRUE)==1)
  for(i in temp) genoU2[,i] <- 1 - genoU2[,i]

  ### identify the minor allele frequency (maf).
  mafS <- colSums(genoS2==1, na.rm=TRUE)/colSums(!is.na(genoS2))
  mafS[mafS > 0.5] <- 1- mafS[mafS > 0.5]
  mafU <- colSums(genoU2==1, na.rm=TRUE)/colSums(!is.na(genoU2))
  mafU[mafU > 0.5] <- 1- mafU[mafU > 0.5]
  
  # keep markers with r2 <= 0.99 with QTLs and maf >= 0.01.
  idxS <- vector()
  for(i in 1:ncol(genoS2)){
    temp <- sapply(qtl, FUN=function(x) suppressWarnings(cor(genoS2[,x], genoS2[,i], use="complete.obs")^2))
    if(all(temp <= 0.99, na.rm=TRUE) & mafS[i] >= 0.01) idxS <- c(idxS, i)
  }
  
  idxU <- vector()
  for(i in 1:ncol(genoU2)){
    temp <- sapply(qtl, FUN=function(x) suppressWarnings(cor(genoU2[,x], genoU2[,i], use="complete.obs")^2))
    if(all(temp <= 0.99, na.rm=TRUE) & mafU[i] >= 0.01) idxU <- c(idxU, i)
  }

  ### run RALLY for the population with selection.
  rcS <- vector()
  afcS <- vector()
  for(i in idxS){
    
    temp <- cbind(phenoS, allele=genoS2[, i])
    
    lm1 <- glm(allele~year, data=temp, family=binomial(link="logit"))
    if(lm1$converged){
      rcS <- rbind(rcS, summary(lm1)$coefficients[2,])
      afcS <- c(afcS, abs(diff(unname(predict(object=lm1, newdata=data.frame(year=c(1,50)), type="response")))))
    } else {
      rcS <- rbind(rcS, rep(NA,4))
      afcS <- c(afcS, NA)
    }

  }
  rcS <- data.frame(marker[idxS,], rcS, afc=afcS)
  colnames(rcS)[4:7] <- c("B", "SE", "Z", "P")

  ### run RALLY for the population without selection.
  rcU <- vector()
  afcU <- vector()
  for(i in idxU){
    
    temp <- cbind(phenoU, allele=genoU2[, i])
    
      lm1 <- glm(allele~year, data=temp, family=binomial(link="logit"))
      if(lm1$converged){
        rcU <- rbind(rcU, summary(lm1)$coefficients[2,])
        afcU <- c(afcU, abs(diff(unname(predict(object=lm1, newdata=data.frame(year=c(1,50)), type="response")))))
      } else {
        rcU <- rbind(rcU, rep(NA,4))
        afcU <- c(afcU, NA)
      }
    
  }
  rcU <- data.frame(marker[idxU,], rcU, afc=afcU)
  colnames(rcU)[4:7] <- c("B", "SE", "Z", "P")
  
  ### apply genomic control (GC) and delta control (DC) using different thresholds for null.
  adj <- data.frame(t=seq(0.05, 0.5, 0.01), deltaS=NA, sigmaS=NA, deltaU=NA, sigmaU=NA)
  pS <- vector()
  pU <- vector()
  
  for(i in 1:length(adj$t)){
    temp <- rcS$Z[!is.na(rcS$Z) & rcS$afc < adj$t[i]]
    temp <- nlm(f=LL.norm, p=c(mean(temp), sd(temp)), x=temp, print.level=0, gradtol=1e-36)
    pS <- cbind(pS, 2*pnorm(q=abs((rcS$Z - temp$estimate[1])/temp$estimate[2]), mean=0, sd=1, lower.tail=FALSE))
    adj$deltaS[i] <- temp$estimate[1]
    adj$sigmaS[i] <- temp$estimate[2]
    rm(temp)
    
    temp <- rcU$Z[!is.na(rcU$Z) & rcU$afc < adj$t[i]]
    temp <- nlm(f=LL.norm, p=c(mean(temp), sd(temp)), x=temp, print.level=0, gradtol=1e-36)
    pU <- cbind(pU, 2*pnorm(q=abs((rcU$Z - temp$estimate[1])/temp$estimate[2]), mean=0, sd=1, lower.tail=FALSE))
    adj$deltaU[i] <- temp$estimate[1]
    adj$sigmaU[i] <- temp$estimate[2]
    rm(temp)
    
  }
  
  adj <- cbind(adj,
               S=colSums(pS < 0.05/nrow(rcS), na.rm=TRUE)/nrow(rcS),
               U=colSums(pU < 0.05/nrow(rcU), na.rm=TRUE)/nrow(rcU))
  
  ### include the un-adjusted results.
  adj <- rbind(data.frame(t=0,
                          deltaS=0,
                          sigmaS=1,
                          deltaU=0,
                          sigmaU=1,
                          S=sum(rcS$P < 0.05/nrow(rcS), na.rm=TRUE)/nrow(rcS),
                          U=sum(rcU$P < 0.05/nrow(rcU), na.rm=TRUE)/nrow(rcU)),
               adj)
  pS <- cbind(rcS$P, pS)
  pU <- cbind(rcU$P, pU)
  
  ### calculate the distance and r2 between a marker and its closest QTL (within chromosome).
  MS <- data.frame(marker, qtl=NA, dist=NA, r2=NA, GWAS=NA)
  for(i in idxS){
    temp <- qtl[ marker$chr[qtl] == marker$chr[i] ]
    if(length(temp) > 0){
      temp <- temp[which.min(abs(temp - i))]
      MS$qtl[i] <- temp
      MS$dist[i] <- abs(temp - i)/10
      MS$r2[i] <- suppressWarnings(cor(genoS2[,temp], genoS2[,i], use="complete.obs")^2)
    }
  }
  MS <- MS[idxS,] # retain only the markers used in RALLY/GWAS.

  MU <- data.frame(marker, qtl=NA, dist=NA, r2=NA, GWAS=NA)
  for(i in idxU){
    temp <- qtl[ marker$chr[qtl] == marker$chr[i] ]
    if(length(temp) > 0){
      temp <- temp[which.min(abs(temp - i))]
      MU$qtl[i] <- temp
      MU$dist[i] <- abs(temp - i)/10
      MU$r2[i] <- suppressWarnings(cor(genoU2[,temp], genoU2[,i], use="complete.obs")^2)
    }
  }
  MU <- MU[idxU,] # retain only the markers used in RALLY/GWAS.
  
  ### add the p-values from RALLY at different null threshold to MS and MU.
  colnames(pS) <- colnames(pU) <- paste("t", formatC(x=c(0,seq(5,50,1)), width=3, flag="0"), sep="")
  MS <- cbind(MS, pS)
  MU <- cbind(MU, pU)
  
  ### impute the marker data because GWASPOLY will not test marker with any missing data.
  ### compute the GRM in S and U population.
  temp <- rrBLUP::A.mat(2*genoS2[, idxS] - 1, return.imputed=TRUE)
  AS <- temp$A
  genoS3 <- temp$imputed + 1
  
  temp <- rrBLUP::A.mat(2*genoU2[, idxU] - 1, return.imputed=TRUE)
  AU <- temp$A
  genoU3 <- temp$imputed + 1
  
  ### GWAS on the population with selection.
  temp <- rcS[,1:3]
  colnames(temp) <- c("Marker", "Chrom", "Position")
  temp$Chrom <- as.factor(temp$Chrom) # Because GWASPOLY requires this.
  temp$Ref <- 0
  temp$Alt <- 1
  dat <- new("GWASpoly.K",
             K=list(all=AS),
             map=temp,
             pheno=phenoS[, c(1,5)],
             geno=genoS3,
             fixed=phenoS[, 2, drop=FALSE],
             ploidy=2)
  z <- set.params(fixed=c("year"), fixed.type=c("numeric"))
  gwas <- GWASpoly(dat, models="additive", traits="PV2", params=z, n.core=8)
  MS$GWAS <- 10^-gwas@scores$PV2$additive
  
  ### envGWAS - not included.
  #dat <- new("GWASpoly.K",
  #           K=list(all=AS),
  #           map=temp,
  #           pheno=phenoS[, c(1,2)],
  #           geno=genoS3,
  #           ploidy=2)
  #gwas <- GWASpoly(dat, models="additive", traits="year", n.core=8)

  ### GWAS on the population without selection.
  temp <- rcU[, 1:3]
  colnames(temp) <- c("Marker", "Chrom", "Position")
  temp$Chrom <- as.factor(temp$Chrom) # Because GWASPOLY requires this.
  temp$Ref <- 0
  temp$Alt <- 1
  dat <- new("GWASpoly.K",
             K=list(all=AU),
             map=temp,
             pheno=phenoU[, c(1,4)],
             geno=genoU3,
             fixed=phenoU[, 2, drop=FALSE],
             ploidy=2)
  z <- set.params(fixed=c("year"), fixed.type=c("numeric"))
  gwas <- GWASpoly(dat, models="additive", traits="PV", params=z, n.core=8)
  MU$GWAS <- 10^-gwas@scores$PV$additive

  outT <- c(outT, list(adj))
  
  save(fgeno, founder, genoS, genoS2, genoS3, genoU, genoU2, genoU3, GV, MS, MU, phenoS, phenoU,
       pop, popS, popU, PV, rcS, rcU, fac, idxS, idxU, mafS, mafU, meff, qtl, SP,
       file=paste("output/sim/sim", formatC(x=h, width=3, flag="0"), ".RData", sep=""))
  
  message(h)
}

saveRDS(outT, "output/sim/outT.RDS")


### check the effect of different thresholds on S and U.
df <- vector()
for(i in 1:100) df <- rbind(df, cbind(outT[[i]], sim=i))

df <- df[abs(df$deltaS) < 0.1 & abs(df$deltaU) < 0.1, ] # remove those that did not converge correctly.

### Figure 2A - proportion of significant markers - selection vs drift.
#####
df.plot <- data.frame(t=unique(df$t),
                      x=tapply(df$U, df$t, mean),
                      y=tapply(df$S-df$U, df$t, mean),
                      adj=c("no",rep("yes",46)))
ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=df.plot[df.plot$adj=="yes", ], aes(x=x, y=y), size=1, color="#999999") +
  geom_point(data=df.plot[df.plot$adj=="no", ], aes(x=x, y=y), size=1, color="#FF0000") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.title=element_text(size=8), axis.text=element_text(size=8)) +
  scale_y_continuous(breaks=c(0, 0.02, 0.04)) +
  xlab("drift") +
  ylab("selection") +
  coord_fixed()

ggsave(filename="output/2A_sel_vs_drift.png",
       width=1.75,
       height=1.5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/2A_sel_vs_drift.svg",
       width=1.75,
       height=1.5,
       units="in",
       scale=4/3)
#####


### identify the simulations for threshold=0.20 with converged MLE.
idx <- df$sim[df$t==df$t[17]] 

### calculate the proportion of significant markers for all simulations.
df1 <- df2 <- df3 <- vector()
for(i in idx){
  load(paste("output/sim/sim", formatC(x=i, width=3, flag="0"), ".RData", sep=""))
  tempS <- MS[, c(1:7, 24)]
  tempU <- MU[, c(1:7, 24)]
  
  df1 <- rbind(df1,
               data.frame(S1=sum(tempS$t020 <= 0.05/nrow(tempS) & tempS$GWAS <= 0.05/nrow(tempS) & tempS$dist <= 5, na.rm=TRUE),
                          S2=sum(tempS$t020 <= 0.05/nrow(tempS) & tempS$GWAS <= 0.05/nrow(tempS) & tempS$dist > 5, na.rm=TRUE),
                          S3=sum(tempS$t020 <= 0.05/nrow(tempS) & tempS$GWAS <= 0.05/nrow(tempS), na.rm=TRUE),
                          U1=sum(tempU$t020 <= 0.05/nrow(tempU) & tempU$GWAS <= 0.05/nrow(tempU) & tempU$dist <= 5, na.rm=TRUE),
                          U2=sum(tempU$t020 <= 0.05/nrow(tempU) & tempU$GWAS <= 0.05/nrow(tempU) & tempU$dist > 5, na.rm=TRUE),
                          U3=sum(tempU$t020 <= 0.05/nrow(tempU) & tempU$GWAS <= 0.05/nrow(tempU), na.rm=TRUE)))
  
  df2 <- rbind(df2, 
               data.frame(S1=sum(tempS$t020 <= 0.05/nrow(tempS) & tempS$GWAS > 0.05/nrow(tempS) & tempS$dist <= 5, na.rm=TRUE),
                          S2=sum(tempS$t020 <= 0.05/nrow(tempS) & tempS$GWAS > 0.05/nrow(tempS) & tempS$dist > 5, na.rm=TRUE),
                          S3=sum(tempS$t020 <= 0.05/nrow(tempS) & tempS$GWAS > 0.05/nrow(tempS), na.rm=TRUE),
                          U1=sum(tempU$t020 <= 0.05/nrow(tempU) & tempU$GWAS > 0.05/nrow(tempU) & tempU$dist <= 5, na.rm=TRUE),
                          U2=sum(tempU$t020 <= 0.05/nrow(tempU) & tempU$GWAS > 0.05/nrow(tempU) & tempU$dist > 5, na.rm=TRUE),
                          U3=sum(tempU$t020 <= 0.05/nrow(tempU) & tempU$GWAS > 0.05/nrow(tempU), na.rm=TRUE)))
  
  df3 <- rbind(df3,
               data.frame(S1=sum(tempS$t020 > 0.05/nrow(tempS) & tempS$GWAS <= 0.05/nrow(tempS) & tempS$dist <= 5, na.rm=TRUE),
                          S2=sum(tempS$t020 > 0.05/nrow(tempS) & tempS$GWAS <= 0.05/nrow(tempS) & tempS$dist > 5, na.rm=TRUE),
                          S3=sum(tempS$t020 > 0.05/nrow(tempS) & tempS$GWAS <= 0.05/nrow(tempS), na.rm=TRUE),
                          U1=sum(tempU$t020 > 0.05/nrow(tempU) & tempU$GWAS <= 0.05/nrow(tempU) & tempU$dist <= 5, na.rm=TRUE),
                          U2=sum(tempU$t020 > 0.05/nrow(tempU) & tempU$GWAS <= 0.05/nrow(tempU) & tempU$dist > 5, na.rm=TRUE),
                          U3=sum(tempU$t020 > 0.05/nrow(tempU) & tempU$GWAS <= 0.05/nrow(tempU), na.rm=TRUE)))
  
}

### calculate the number of significant markers tagging QTLs for all simulations.
qc <- vector()
for(i in idx){
  load(paste("output/sim/sim", formatC(x=i, width=3, flag="0"), ".RData", sep=""))
  z <- rep(1, 19000)
  z[idxS] <- MS$t020
  z <- which(z <= 0.05/length(idxS))
  q <- rep(0, 32)
  for(j in 1:length(z)){
    temp <- abs(qtl - z[j])
    if(any(temp <= 50)) q[which.min(temp)] <- q[which.min(temp)] + 1
  }
  qc <- rbind(qc, q)
}

### Figure S3 - Count of significant markers within 5cM of QTLs.
#####
d1 <- data.frame(QTL=sort(rep(1:32, nrow(qc))), marker=c(qc))
d2 <- data.frame(QTL=1:32,
                 marker=colSums(qc)/nrow(qc))

ggplot() +
  annotate("rect", xmin=0.5, xmax=2.5, ymin=-Inf, ymax=Inf, fill="#2171B5", color=NA, alpha=0.75) +
  annotate("rect", xmin=2.5, xmax=6.5, ymin=-Inf, ymax=Inf, fill="#6BAED6", color=NA, alpha=0.75) +
  annotate("rect", xmin=6.5, xmax=14.5, ymin=-Inf, ymax=Inf, fill="#BDD7E7", color=NA, alpha=0.75) +
  annotate("rect", xmin=14.5, xmax=32.5, ymin=-Inf, ymax=Inf, fill="#EFF3FF", color=NA, alpha=0.75) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_jitter(data=d1, aes(x=QTL, y=marker), width=0.25, size=1) +
  geom_line(data=d2, aes(x=QTL, y=marker), color="#FF0000", size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  scale_y_continuous(breaks=c(0,10,20,30,40,50)) +
  scale_x_continuous(breaks=1:32, expand=c(0,0))

ggsave(filename="output/S3_count_sim_sig_markers.png",
       width=6.5,
       height=3,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/S3_count_sim_sig_markers.svg",
       width=6.5,
       height=3,
       units="in",
       scale=4/3)
#####


### Figure 2B - RALLY sig + GWAS sig
#####
df.plot <- rbind(data.frame(count=df1$S1, x=center.jitter(x=1, y=df1$S1, width=0.4)),
                 data.frame(count=df1$S2, x=center.jitter(x=2, y=df1$S2, width=0.4)),
                 data.frame(count=df1$S3, x=center.jitter(x=3, y=df1$S3, width=0.4)),
                 data.frame(count=df1$U1, x=center.jitter(x=4, y=df1$U1, width=0.4)),
                 data.frame(count=df1$U2, x=center.jitter(x=5, y=df1$U2, width=0.4)),
                 data.frame(count=df1$U3, x=center.jitter(x=6, y=df1$U3, width=0.4)))

temp <- data.frame(x=1:6,
                   count=sapply(1:6, FUN=function(x) median(df1[,x])))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=df.plot, aes(x=x, y=count), size=0.5, color="#999999") +
  geom_point(data=temp, aes(x=x, y=count), size=1, color="#FF0000") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=8)) +
  theme(axis.title=element_text(size=8), axis.text.y=element_text(size=8)) +
  theme(plot.subtitle=element_text(size=8)) +
  scale_x_continuous(breaks=1:6, labels=c("S, \u2264 5 cM", "S, > 5 cM", "S, total", "U, \u2264 5 cM", "U, > 5 cM", "U, total")) +
  xlab("population, marker") +
  ylab("count") +
  labs(subtitle="Sig RALLY + Sig GWAS")

ggsave(filename="output/2B_sigRALLY_sigGWAS.png",
       width=1.75,
       height=1.5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/2B_sigRALLY_sigGWAS.svg",
       width=1.75,
       height=1.5,
       units="in",
       scale=4/3)

#####

### Figure 2C - RALLY sig + GWAS NS
#####
df.plot <- rbind(data.frame(count=df2$S1, x=center.jitter(x=1, y=df2$S1, width=0.4)),
                 data.frame(count=df2$S2, x=center.jitter(x=2, y=df2$S2, width=0.4)),
                 data.frame(count=df2$S3, x=center.jitter(x=3, y=df2$S3, width=0.4)),
                 data.frame(count=df2$U1, x=center.jitter(x=4, y=df2$U1, width=0.4)),
                 data.frame(count=df2$U2, x=center.jitter(x=5, y=df2$U2, width=0.4)),
                 data.frame(count=df2$U3, x=center.jitter(x=6, y=df2$U3, width=0.4)))

temp <- data.frame(x=1:6,
                   count=sapply(1:6, FUN=function(x) median(df2[,x])))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=df.plot, aes(x=x, y=count), size=0.5, color="#999999") +
  geom_point(data=temp, aes(x=x, y=count), size=1, color="#FF0000") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=8)) +
  theme(axis.title=element_text(size=8), axis.text.y=element_text(size=8)) +
  theme(plot.subtitle=element_text(size=8)) +
  scale_x_continuous(breaks=1:6, labels=c("S, \u2264 5 cM", "S, > 5 cM", "S, total", "U, \u2264 5 cM", "U, > 5 cM", "U, total")) +
  xlab("population, marker") +
  ylab("count") +
  labs(subtitle="Sig RALLY + NS GWAS")

ggsave(filename="output/2C_sigRALLY_nsGWAS.png",
       width=1.75,
       height=1.5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/2C_sigRALLY_sigGWAS.svg",
       width=1.75,
       height=1.5,
       units="in",
       scale=4/3)

#####

### Figure 2D - RALLY NS + GWAS sig
#####
df.plot <- rbind(data.frame(count=df3$S1, x=center.jitter(x=1, y=df3$S1, width=0.4)),
                 data.frame(count=df3$S2, x=center.jitter(x=2, y=df3$S2, width=0.4)),
                 data.frame(count=df3$S3, x=center.jitter(x=3, y=df3$S3, width=0.4)),
                 data.frame(count=df3$U1, x=center.jitter(x=4, y=df3$U1, width=0.4)),
                 data.frame(count=df3$U2, x=center.jitter(x=5, y=df3$U2, width=0.4)),
                 data.frame(count=df3$U3, x=center.jitter(x=6, y=df3$U3, width=0.4)))

temp <- data.frame(x=1:6,
                   count=sapply(1:6, FUN=function(x) median(df3[,x])))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=df.plot, aes(x=x, y=count), size=0.5, color="#999999") +
  geom_point(data=temp, aes(x=x, y=count), size=1, color="#FF0000") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=8)) +
  theme(axis.title=element_text(size=8), axis.text.y=element_text(size=8)) +
  theme(plot.subtitle=element_text(size=8)) +
  scale_x_continuous(breaks=1:6, labels=c("S, \u2264 5 cM", "S, > 5 cM", "S, total", "U, \u2264 5 cM", "U, > 5 cM", "U, total")) +
  xlab("population, marker") +
  ylab("count") +
  labs(subtitle="NS RALLY + Sig GWAS")

ggsave(filename="output/2D_nsRALLY_sigGWAS.png",
       width=1.75,
       height=1.5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/2D_nsRALLY_sigGWAS.svg",
       width=1.75,
       height=1.5,
       units="in",
       scale=4/3)

#####

### Figure 2E,F,G,H - Manhattan plot for RALLY and GWAS in one simulation.
#####

load("output/sim/sim001.RData")

### adjust the positions for the manhattan plot.
### create the marker name.
glen <- seq(1.0, 2.8, 0.2)
marker <- data.frame(chr=unlist(lapply(1:10, FUN=function(x) rep(x, glen[x]*1000))),
                     pos=unlist(lapply(1:10, FUN=function(x) 0:(glen[x]*1000-1)))/10)
temp <- unname(c(tapply(marker$pos, marker$chr, max)))
temp <- cumsum(c(0, temp+20))
dfS <- data.frame(chr=MS$chr,
                  pos=unlist(lapply(1:10, FUN=function(x) MS$pos[MS$chr==x] + temp[x])),
                  GWAS=-log10(MS$GWAS),
                  RALLY=-log10(MS$t020))
dfS <- melt(dfS, id.vars=c("chr", "pos"))
dfS$variable <- factor(dfS$variable, levels=c("RALLY", "GWAS"))

anno <- data.frame(chr=marker$chr, 
                   pos=unlist(lapply(1:10, FUN=function(x) marker$pos[marker$chr==x] + temp[x])))
anno <- cbind(anno[qtl,], effect=meff[qtl])
anno$start <- anno$pos - 5
anno$end <- anno$pos + 5
anno$effect <- as.factor(anno$effect)
anno$effect <- factor(anno$effect, labels=c("a", "2a", "3a", "4a"))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_rect(data=anno, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill=effect), alpha=0.5) +
  geom_point(data=dfS, aes(x=pos, y=value), size=0.5, color="#555555", na.rm=TRUE) +
  facet_wrap(vars(variable), ncol=1) +
  annotate("segment", x=-Inf, xend=Inf, y=-log10(0.05/nrow(MS)), yend=-log10(0.05/nrow(MS)), color="#000000", linetype=2, size=0.2) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.title.x=element_blank(), axis.title.y=element_text(size=8)) +
  theme(axis.text=element_text(size=8)) +
  theme(strip.background=element_blank(), strip.text=element_text(hjust=0, size=8)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8), legend.text=element_text(size=8)) +
  scale_fill_manual(values=c("#BAE4B3", "#74C476", "#238B45", "#006D2C")) +
  scale_x_continuous(expand=c(0.02,0), breaks=(temp[-1] + temp[-11] - 50)/2, labels=1:10) +
  scale_y_continuous(expand=c(0.02,0)) +
  ylab(expression("-log"[10]*"P"))

ggsave(filename="output/2EF_S001.png",
       width=3.5,
       height=2.5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/2EF_S001.svg",
       width=3.5,
       height=2.5,
       units="in",
       scale=4/3)


dfU <- data.frame(chr=MU$chr,
                  pos=unlist(lapply(1:10, FUN=function(x) MU$pos[MU$chr==x] + temp[x])),
                  GWAS=-log10(MU$GWAS),
                  RALLY=-log10(MU$t020))
dfU <- melt(dfU, id.vars=c("chr", "pos"))
dfU$variable <- factor(dfU$variable, levels=c("RALLY", "GWAS"))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_rect(data=anno, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill=effect), alpha=0.5) +
  geom_point(data=dfU, aes(x=pos, y=value), size=0.5, color="#555555", na.rm=TRUE) +
  facet_wrap(vars(variable), ncol=1) +
  annotate("segment", x=-Inf, xend=Inf, y=-log10(0.05/nrow(MU)), yend=-log10(0.05/nrow(MU)), color="#000000", linetype=2, size=0.2) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.title.x=element_blank(), axis.title.y=element_text(size=8)) +
  theme(axis.text=element_text(size=8)) +
  theme(strip.background=element_blank(), strip.text=element_text(hjust=0, size=8)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_text(size=8), legend.text=element_text(size=8)) +
  scale_fill_manual(values=c("#BAE4B3", "#74C476", "#238B45", "#006D2C")) +
  scale_x_continuous(expand=c(0.02,0), breaks=(temp[-1] + temp[-11] - 50)/2, labels=1:10) +
  scale_y_continuous(expand=c(0.02,0), limits=c(0, max(dfS$value))) +
  ylab(expression("-log"[10]*"P"))

ggsave(filename="output/2GH_U001.png",
       width=3.5,
       height=2.5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/2GH_U001.svg",
       width=3.5,
       height=2.5,
       units="in",
       scale=4/3)

#####

### Figure I,J - histogram of p-values.
#####
ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_histogram(aes(x=MS$t020), breaks=seq(0,1,0.05), fill="#CCCCCC", color="#000000") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  scale_y_continuous(breaks=c(0, 1000, 2000)) +
  xlab("p-values")
ggsave(filename="output/2I_hist_S001.png",
       width=3.5,
       height=1,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/2I_hist_S001.svg",
       width=3.5,
       height=1,
       units="in",
       scale=4/3)

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_histogram(aes(x=MU$t020), breaks=seq(0,1,0.05), fill="#CCCCCC", color="#000000") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
  scale_y_continuous(breaks=c(0, 1000, 2000), limits=c(0, 2480)) +
  xlab("p-values")
ggsave(filename="output/2J_hist_U001.png",
       width=3.5,
       height=1,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/2J_hist_U001.svg",
       width=3.5,
       height=1,
       units="in",
       scale=4/3)
#####

# save the plot data in case we need it.
save(idx, df, df1, df2, df3, dfS, dfU, file="output/2AJ.RData")

### Table S1.
temp <- lapply(unique(df$t), FUN=function(x) df[df$t==x, ])
df0 <- vector()
for(i in 1:47){
  df0 <- rbind(df0, data.frame(t=temp[[i]]$t[1],
                               t(sapply(2:7, FUN=function(x) mean(temp[[i]][,x], na.rm=TRUE))),
                               t(sapply(2:7, FUN=function(x) min(temp[[i]][,x], na.rm=TRUE))),
                               t(sapply(2:7, FUN=function(x) max(temp[[i]][,x], na.rm=TRUE))),
                               nsim=nrow(temp[[i]])))
}
colnames(df0)[2:19] <- c(paste(colnames(df)[2:7], ".mean", sep=""),
                         paste(colnames(df)[2:7], ".min", sep=""),
                         paste(colnames(df)[2:7], ".max", sep=""))
write.csv(df0, "output/S1_sim_summary_null_threshold.csv", quote=FALSE, row.names=FALSE)
