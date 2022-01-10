library(AlphaSimR)
library(sommer)
library(rrBLUP)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(shadowtext)

# set working directory.
setwd()

# initial values.
n.founder <- 10
gmap <- lapply(seq(1.0, 2.8, 0.2), FUN=function(x) seq(0, x-0.01, 0.01)) # 10 chromosomes, markers at 1 cM spacing.
sel=c(0.10, 0.10, 0.10)

# function to create founders and simulate one generation of selection.
simsel <- function(n.founder, gmap, sel){
  
  # set up the founder haplotypes.
  founder <- runMacs(nInd=n.founder,
                     nChr=length(gmap),
                     segSites=sapply(1:length(gmap), FUN=function(x) length(gmap[[x]])),
                     inbred=TRUE,
                     species="GENERIC",
                     ploidy=2L,
                     manualGenLen=sapply(1:length(gmap), FUN=function(x) max(gmap[[x]])))
  SP <- SimParam$new(founder)
  founder <- newPop(founder, simParam=SP)
  
  # get the allele count/frequency in the founder.
  af <- colSums(pullSegSiteGeno(founder, simParam=SP)==2)
  
  # identify 15 QTLs for each trait (3 traits total).
  # 9 unique QTLs, 2+2 QTLs shared between any two traits, 2 QTLs shared among 3 traits.
  qtl1 <- qtl2 <- qtl3 <- vector()
  meff1 <- meff2 <- meff3 <- rep(0, length(af))
  
  # largest effect QTL (0,0,1).
  temp <- which(af==1)
  idx <- sample(temp, 1, replace=FALSE)
  qtl1 <- c(qtl1, idx[1])
  qtl2 <- c(qtl2, idx[1])
  qtl3 <- c(qtl3, idx[1])
  meff1[idx[1]] <- meff2[idx[1]] <- meff3[idx[1]] <- 4
  
  # second largest QTL (0,1+1,0).
  temp <- which(af==2)
  idx <- sample(temp, 3, replace=FALSE)
  qtl1 <- c(qtl1, idx[c(1,2)])
  qtl2 <- c(qtl2, idx[c(1,3)])
  qtl3 <- c(qtl3, idx[c(2,3)])
  meff1[idx[c(1,2)]] <- 3
  meff2[idx[c(1,3)]] <- 3
  meff3[idx[c(2,3)]] <- 3
  
  # third largest QTL (3,0,1).
  temp <- which(af==3)
  idx <- sample(temp, 10, replace=FALSE)
  qtl1 <- c(qtl1, idx[c(1,2,3,4)])
  qtl2 <- c(qtl2, idx[c(1,5,6,7)])
  qtl3 <- c(qtl3, idx[c(1,8,9,10)])
  meff1[idx[c(1,2,3,4)]] <- 2
  meff2[idx[c(1,5,6,7)]] <- 2
  meff3[idx[c(1,8,9,10)]] <- 2
  
  # smallest QTL (6, 1+1, 0).
  temp <- which(af==4 | af==5)
  idx <- sample(temp, 21, replace=FALSE)
  qtl1 <- c(qtl1, idx[c(1,2,4,5,6,7,8,9)])
  qtl2 <- c(qtl2, idx[c(1,3,10,11,12,13,14,15)])
  qtl3 <- c(qtl3, idx[c(2,3,16,17,18,19,20,21)])
  meff1[idx[c(1,2,4,5,6,7,8,9)]] <- 1
  meff2[idx[c(1,3,10,11,12,13,14,15)]] <- 1
  meff3[idx[c(2,3,16,17,18,19,20,21)]] <- 1
  
  # make random crosses among the founders and create 300 DH individuals.
  pop0 <- randCross(pop=founder,
                    nCrosses=300,
                    nProgeny=1,
                    simParam=SP)
  pop0 <- makeDH(pop=pop0,
                 nDH=1,
                 simParam=SP)
  
  # get pop0 marker data.
  geno0 <- pullSegSiteGeno(pop0, simParam=SP)
  
  # scale the qtl effects such that Vg=1 in the 300 DH.
  meff1 <- meff1/sd(c(geno0%*%meff1))
  meff2 <- meff2/sd(c(geno0%*%meff2))
  meff3 <- meff3/sd(c(geno0%*%meff3))
  
  # calculate the genetic value in pop0.
  gv0 <- cbind(trait1=c(geno0%*%meff1),
               trait2=c(geno0%*%meff2),
               trait3=c(geno0%*%meff3))
  
  # add residual effect such that Ve=1.
  pv0 <- gv0 + matrix(rnorm(nrow(gv0)*ncol(gv0), 0, 1), nrow=nrow(gv0), ncol=ncol(gv0))

  # select on trait 1.
  idx1 <- which(rank(-pv0[,1]) <= floor(300*sel[1]))
  s1 <- colSums(pv0[idx1,])/length(idx1) - colSums(pv0)/nrow(pv0)
  
  # select on trait 1,2.
  idx2 <- which(rank(-pv0[,1]) <= floor(300*sel[2]) | rank(-pv0[,2]) <= floor(300*sel[2]))
  s2 <- colSums(pv0[idx2,])/length(idx2) - colSums(pv0)/nrow(pv0)

  # select on trait 1,2,3.
  idx3 <- which(rank(-pv0[,1]) <= floor(300*sel[3]) | rank(-pv0[,2]) <= floor(300*sel[3]) | rank(-pv0[,3]) <= floor(300*sel[3]))
  s3 <- colSums(pv0[idx3,])/length(idx3) - colSums(pv0)/nrow(pv0)
  
  # extract the selected individuals, make crosses and DH.
  # select on trait 1.
  pop1 <- selectInd(pop=pop0,
                    nInd=length(idx1),
                    use="rand",
                    candidates=idx1,
                    simParam=SP)
  pop1 <- randCross(pop=pop1,
                    nCrosses=300,
                    nProgeny=1,
                    simParam=SP)
  pop1 <- makeDH(pop=pop1,
                 nDH=1,
                 simParam=SP)
  
  # select on trait 1,2.
  pop2 <- selectInd(pop=pop0,
                    nInd=length(idx2),
                    use="rand",
                    candidates=idx2,
                    simParam=SP)
  pop2 <- randCross(pop=pop2,
                    nCrosses=300,
                    nProgeny=1,
                    simParam=SP)
  pop2 <- makeDH(pop=pop2,
                 nDH=1,
                 simParam=SP)
  
  # select on trait 1,2,3.
  pop3 <- selectInd(pop=pop0,
                    nInd=length(idx3),
                    use="rand",
                    candidates=idx3,
                    simParam=SP)
  pop3 <- randCross(pop=pop3,
                    nCrosses=300,
                    nProgeny=1,
                    simParam=SP)
  pop3 <- makeDH(pop=pop3,
                 nDH=1,
                 simParam=SP)
  
  # get pop1, pop2, pop3 marker data.
  geno1 <- pullSegSiteGeno(pop1, simParam=SP)
  geno2 <- pullSegSiteGeno(pop2, simParam=SP)
  geno3 <- pullSegSiteGeno(pop3, simParam=SP)
  
  # calculate the genetic value in pop1.
  gv1 <- cbind(trait1=c(geno1%*%meff1),
               trait2=c(geno1%*%meff2),
               trait3=c(geno1%*%meff3))
  
  # add residual effect such that Ve=1.
  pv1 <- gv1 + matrix(rnorm(nrow(gv1)*ncol(gv1), 0, 1), nrow=nrow(gv1), ncol=ncol(gv1))
  
  # calculate the genetic value in pop2.
  gv2 <- cbind(trait1=c(geno2%*%meff1),
               trait2=c(geno2%*%meff2),
               trait3=c(geno2%*%meff3))
  
  # add residual effect such that Ve=1.
  pv2 <- gv2 + matrix(rnorm(nrow(gv2)*ncol(gv2), 0, 1), nrow=nrow(gv2), ncol=ncol(gv2))
  
  # calculate the genetic value in pop1.
  gv3 <- cbind(trait1=c(geno3%*%meff1),
               trait2=c(geno3%*%meff2),
               trait3=c(geno3%*%meff3))

  # add residual effect such that Ve=1.
  pv3 <- gv3 + matrix(rnorm(nrow(gv3)*ncol(gv3), 0, 1), nrow=nrow(gv3), ncol=ncol(gv3))

  # estimate G and P for pop0 using mixed model.
  df0 <- data.frame(ID=rownames(geno0), pv0)
  vcov0 <- calc.vcov.sim(geno=geno0, pheno=df0, trait.idx=2:4)

  # estimate G and P for pop0 from actual var/cov of the genetic and phenotypic values.
  vcov0.b <- list(G=var(gv0), E=var(pv0-gv0), P=var(pv0))
  
  # calculate Z1, Z2, Z3.
  Z1 <- colSums(pv1)/nrow(pv1) - colSums(pv0)/nrow(pv0)
  Z2 <- colSums(pv2)/nrow(pv2) - colSums(pv0)/nrow(pv0)
  Z3 <- colSums(pv3)/nrow(pv3) - colSums(pv0)/nrow(pv0)
  
  # estimate beta (B1, B2, B3).
  B1 <- solve(a=vcov0$G, b=Z1)
  B2 <- solve(a=vcov0$G, b=Z2)
  B3 <- solve(a=vcov0$G, b=Z3)
  
  B1.b <- solve(a=vcov0.b$G, b=Z1)
  B2.b <- solve(a=vcov0.b$G, b=Z2)
  B3.b <- solve(a=vcov0.b$G, b=Z3)
  
  # estimate selection differential (S1, S2, S3).
  S1 <- c(vcov0$P%*%B1)
  S2 <- c(vcov0$P%*%B2)
  S3 <- c(vcov0$P%*%B3)
  
  S1.b <- c(vcov0.b$P%*%B1.b)
  S2.b <- c(vcov0.b$P%*%B2.b)
  S3.b <- c(vcov0.b$P%*%B3.b)
  
  # estimate selection intensity (I1, I2, I3).
  I1 <- S1/sqrt(diag(vcov0$P))
  I2 <- S2/sqrt(diag(vcov0$P))
  I3 <- S3/sqrt(diag(vcov0$P))
  
  I1.b <- S1.b/sqrt(diag(vcov0.b$P))
  I2.b <- S2.b/sqrt(diag(vcov0.b$P))
  I3.b <- S3.b/sqrt(diag(vcov0.b$P))
  
  # partition Z into direct and indirect.
  Z1.part <- cbind(Z.direct=diag(vcov0$G)*B1, Z.indirect=Z1-diag(vcov0$G)*B1)
  Z2.part <- cbind(Z.direct=diag(vcov0$G)*B2, Z.indirect=Z2-diag(vcov0$G)*B2)
  Z3.part <- cbind(Z.direct=diag(vcov0$G)*B3, Z.indirect=Z3-diag(vcov0$G)*B3)
  
  Z1.part.b <- cbind(Z.direct=diag(vcov0.b$G)*B1.b, Z.indirect=Z1-diag(vcov0.b$G)*B1.b)
  Z2.part.b <- cbind(Z.direct=diag(vcov0.b$G)*B2.b, Z.indirect=Z2-diag(vcov0.b$G)*B2.b)
  Z3.part.b <- cbind(Z.direct=diag(vcov0.b$G)*B3.b, Z.indirect=Z3-diag(vcov0.b$G)*B3.b)
  
  # partition S into direct and indirect.
  S1.part <- cbind(S.direct=diag(vcov0$P)*B1, S.indirect=S1-diag(vcov0$P)*B1)
  S2.part <- cbind(S.direct=diag(vcov0$P)*B2, S.indirect=S2-diag(vcov0$P)*B2)
  S3.part <- cbind(S.direct=diag(vcov0$P)*B3, S.indirect=S3-diag(vcov0$P)*B3)
  
  S1.part.b <- cbind(S.direct=diag(vcov0.b$P)*B1.b, S.indirect=S1.b-diag(vcov0.b$P)*B1.b)
  S2.part.b <- cbind(S.direct=diag(vcov0.b$P)*B2.b, S.indirect=S2.b-diag(vcov0.b$P)*B2.b)
  S3.part.b <- cbind(S.direct=diag(vcov0.b$P)*B3.b, S.indirect=S3.b-diag(vcov0.b$P)*B3.b)
  
  # partition I into direct and indirect.
  I1.part <- S1.part/sqrt(diag(vcov0$P))
  I2.part <- S2.part/sqrt(diag(vcov0$P))
  I3.part <- S3.part/sqrt(diag(vcov0$P))
  colnames(I1.part) <- colnames(I2.part) <- colnames(I3.part) <- c("I.direct", "I.indirect")
  
  I1.part.b <- S1.part.b/sqrt(diag(vcov0.b$P))
  I2.part.b <- S2.part.b/sqrt(diag(vcov0.b$P))
  I3.part.b <- S3.part.b/sqrt(diag(vcov0.b$P))
  colnames(I1.part.b) <- colnames(I2.part.b) <- colnames(I3.part.b) <- c("I.direct", "I.indirect")
  
  # get the true selection intensity and beta.
  i1 <- s1/sqrt(diag(vcov0.b$P))
  i2 <- s2/sqrt(diag(vcov0.b$P))
  i3 <- s3/sqrt(diag(vcov0.b$P))
  
  b1 <- solve(a=vcov0.b$P, b=s1)
  b2 <- solve(a=vcov0.b$P, b=s2)
  b3 <- solve(a=vcov0.b$P, b=s3)
  
  z1 <- c(vcov0.b$G%*%b1)
  z2 <- c(vcov0.b$G%*%b2)
  z3 <- c(vcov0.b$G%*%b3)
  
  # get the true partitions.
  s1.part <- cbind(S.direct=diag(vcov0.b$P)*b1, S.indirect=s1-diag(vcov0.b$P)*b1)
  s2.part <- cbind(S.direct=diag(vcov0.b$P)*b2, S.indirect=s2-diag(vcov0.b$P)*b2)
  s3.part <- cbind(S.direct=diag(vcov0.b$P)*b3, S.indirect=s3-diag(vcov0.b$P)*b3)
  
  i1.part <- s1.part/sqrt(diag(vcov0.b$P))
  i2.part <- s2.part/sqrt(diag(vcov0.b$P))
  i3.part <- s3.part/sqrt(diag(vcov0.b$P))
  colnames(i1.part) <- colnames(i2.part) <- colnames(i3.part) <- c("I.direct", "I.indirect")
  
  z1.part <- cbind(Z.direct=diag(vcov0.b$G)*b1, Z.indirect=z1-diag(vcov0.b$G)*b1)
  z2.part <- cbind(Z.direct=diag(vcov0.b$G)*b2, Z.indirect=z2-diag(vcov0.b$G)*b2)
  z3.part <- cbind(Z.direct=diag(vcov0.b$G)*b3, Z.indirect=z3-diag(vcov0.b$G)*b3)
  
  # compile the outputs.
  out.mm1 <- data.frame(Z=Z1, Z1.part, B=B1, S=S1, S1.part, I=I1, I1.part)
  out.sim1 <- data.frame(Z=Z1, Z1.part.b, B=B1.b, S=S1.b, S1.part.b, I=I1.b, I1.part.b)  
  out.true1 <- data.frame(Z=z1, z1.part, B=b1, S=s1, s1.part, I=i1, i1.part)
  
  out.mm2 <- data.frame(Z=Z2, Z2.part, B=B2, S=S2, S2.part, I=I2, I2.part)
  out.sim2 <- data.frame(Z=Z2, Z2.part.b, B=B2.b, S=S2.b, S2.part.b, I=I2.b, I2.part.b)  
  out.true2 <- data.frame(Z=z2, z2.part, B=b2, S=s2, s2.part, I=i2, i2.part)

  out.mm3 <- data.frame(Z=Z3, Z3.part, B=B3, S=S3, S3.part, I=I3, I3.part)
  out.sim3 <- data.frame(Z=Z3, Z3.part.b, B=B3.b, S=S3.b, S3.part.b, I=I3.b, I3.part.b)  
  out.true3 <- data.frame(Z=z3, z3.part, B=b3, S=s3, s3.part, I=i3, i3.part)
  
  return(list(sel1=list(mm=out.mm1, sim=out.sim1, true=out.true1),
              sel2=list(mm=out.mm2, sim=out.sim2, true=out.true2),
              sel3=list(mm=out.mm3, sim=out.sim3, true=out.true3),
              vcov.mm=vcov0,
              vcov.sim=vcov0.b))
}

# for convenience, let's make a function that takes geno and pheno and calculate G and P.
# pheno first columns must be ID.
# geno is in 0/1/2 format and missing data will be imputed.
calc.vcov.sim <- function(geno, pheno, trait.idx){
  
  # get the number of traits.
  n.trait <- length(trait.idx)
  
  # get the trait names.
  trait <- colnames(pheno)[trait.idx]
  
  # calculate the kinship matrix.
  A <- rrBLUP::A.mat(geno - 1)
  
  # get the mean of A-diagonals (need this to correct for G).
  da <- mean(diag(A))
  
  # create a matrix of indices for all possible trait pairs.
  tp <- cbind(sort(rep(trait.idx, n.trait)),
              rep(trait.idx, n.trait))
  tp <- tp[tp[,2] >= tp[,1], ]
  
  # identifies the matrix position for var/cov.
  mat.idx <- matrix(NA, nrow=n.trait, ncol=n.trait)
  mat.idx[lower.tri(mat.idx, diag=TRUE)] <- 1:nrow(tp)
  mat.idx[upper.tri(mat.idx, diag=FALSE)] <- t(mat.idx)[upper.tri(t(mat.idx), diag=FALSE)]
  mat.idx <- c(mat.idx)
  
  # create empty vector to store var/cov for G and E.
  vcov.G <- vcov.E <- vector()
  
  # loop to fit the mixed model (univariate and bivariate).
  for(j in 1:nrow(tp)){
    
    # univariate
    if(tp[j,1]==tp[j,2]){
      temp <- pheno[, c(1, tp[j,1])]
      colnames(temp) <- c("ID", "x")
      mm <- tryCatch(summary(sommer::mmer(x~1,
                                          random=~sommer::vs(ID, Gu=A),
                                          rcov=~units,
                                          data=temp,
                                          date.warning=FALSE,
                                          verbose=FALSE))$varcomp, error=function(e){})
      if(!is.null(mm)){
        vcov.G <- c(vcov.G, mm[1,1])
        vcov.E <- c(vcov.E, mm[2,1])
      } else {
        vcov.G <- c(vcov.G, NA)
        vcov.E <- c(vcov.E, NA)
      }
      
    # bivariate.
    } else {
      temp <- pheno[, c(1, tp[j,1], tp[j,2])]
      colnames(temp) <- c("ID", "x", "y")
      mm <- tryCatch(summary(sommer::mmer(cbind(x,y)~1,
                                          random=~sommer::vs(ID, Gu=A),
                                          rcov=~units,
                                          data=temp,
                                          date.warning=FALSE,
                                          verbose=FALSE))$varcomp, error=function(e){})
      if(!is.null(mm)){
        vcov.G <- c(vcov.G, mm[2,1])
        vcov.E <- c(vcov.E, mm[5,1])
      } else {
        vcov.G <- c(vcov.G, NA)
        vcov.E <- c(vcov.E, NA)
      }
    }
    
    # update on the loop.
    message(j, "/", nrow(tp), " ... done")
    
  }
  
  # create G, E and P matrices.
  G <- matrix(vcov.G[mat.idx], nrow=n.trait, ncol=n.trait)
  rownames(G) <- colnames(G) <- trait
  G <- da*G
  
  E <- matrix(vcov.E[mat.idx], nrow=n.trait, ncol=n.trait)
  rownames(E) <- colnames(E) <- trait
  
  P <- G + E
  
  # return the results.
  return(list(G=G, E=E, P=P))
  
}


out <- list()
for(i in 1:100){
  out[[i]] <- simsel(n.founder=10,
                     gmap=lapply(seq(1.0, 2.8, 0.2), FUN=function(x) seq(0, x-0.01, 0.01)),
                     sel=c(0.10, 0.10, 0.10))
  if(i%%10==0) message(i)
}
saveRDS(out, "output/simsel_20211014.RDS")

### G, E and P.
#####
# plot G covariances.
Gc <- vector()
for(i in 1:100){
  Gc <- rbind(Gc, data.frame(mm=out[[i]]$vcov.mm$G[lower.tri(out[[i]]$vcov.mm$G, diag=FALSE)],
                             true=out[[i]]$vcov.sim$G[lower.tri(out[[i]]$vcov.sim$G, diag=FALSE)],
                             trait=c("1_2", "1_3", "2_3"),
                             sim=i))
}
Gc$trait <- as.factor(Gc$trait)
Gc$trait <- factor(Gc$trait, labels=c(paste("r(1-2) = ", formatC(x=cor(Gc$mm[Gc$trait=="1_2"], Gc$true[Gc$trait=="1_2"]), digits=3, flag="#"), sep=""),
                                      paste("r(1-3) = ", formatC(x=cor(Gc$mm[Gc$trait=="1_3"], Gc$true[Gc$trait=="1_3"]), digits=3, flag="#"), sep=""),
                                      paste("r(2-3) = ", formatC(x=cor(Gc$mm[Gc$trait=="2_3"], Gc$true[Gc$trait=="2_3"]), digits=3, flag="#"), sep="")))

temp <- ggplot() +
  annotate("rect", xmin=0.59, xmax=Inf, ymin=-Inf, ymax=0.17, fill="#F2F2F2", color=NA) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Gc, aes(x=true, y=mm, color=trait), size=0.5) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(legend.justification=c(1,0), legend.position=c(1,0), legend.title=element_blank()) +
  theme(legend.text=element_text(size=4), legend.key.size=unit(0.5, "lines")) +
  theme(legend.background=element_blank()) +
  theme(axis.text=element_text(size=7)) +
  scale_x_continuous(limits=range(c(Gc$mm, Gc$true))) +
  scale_y_continuous(limits=range(c(Gc$mm, Gc$true))) +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  ylab("estimated") +
  coord_fixed()
ggsave(plot=temp,
       filename="output/simsel/png/Gcov.png",
       width=1.7,
       height=1.7,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/Gcov.svg",
       width=1.7,
       height=1.7,
       units="in",
       scale=4/3)


# plot G variances.
Gv <- vector()
for(i in 1:100){
  Gv <- rbind(Gv, data.frame(mm=diag(out[[i]]$vcov.mm$G),
                             true=diag(out[[i]]$vcov.sim$G),
                             trait=1:3,
                             sim=i))
}
Gv$trait <- as.factor(Gv$trait)

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  annotate("segment", x=-Inf, xend=Inf, y=1, yend=1, color="#FF0000") +
  geom_boxplot(data=Gv, aes(x=trait, y=mm), outlier.shape=NA) +
  geom_jitter(data=Gv, aes(x=trait, y=mm, color=trait), size=0.5, height=0, width=0.25) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=7)) +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  ylab("estimated") +
  xlab("trait")
ggsave(plot=temp,
       filename="output/simsel/png/Gvar.png",
       width=1.7,
       height=1.7,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/Gvar.svg",
       width=1.7,
       height=1.7,
       units="in",
       scale=4/3)


# plot E (residual) covariances.
Ec <- vector()
for(i in 1:100){
  Ec <- rbind(Ec, data.frame(mm=out[[i]]$vcov.mm$E[lower.tri(out[[i]]$vcov.mm$E, diag=FALSE)],
                             true=out[[i]]$vcov.sim$E[lower.tri(out[[i]]$vcov.sim$E, diag=FALSE)],
                             trait=c("1_2", "1_3", "2_3"),
                             sim=i))
}

Ec$trait <- as.factor(Ec$trait)
Ec$trait <- factor(Ec$trait, labels=c(paste("r(1-2) = ", formatC(x=cor(Ec$mm[Ec$trait=="1_2"], Ec$true[Ec$trait=="1_2"]), digits=3, flag="#"), sep=""),
                                      paste("r(1-3) = ", formatC(x=cor(Ec$mm[Ec$trait=="1_3"], Ec$true[Ec$trait=="1_3"]), digits=3, flag="#"), sep=""),
                                      paste("r(2-3) = ", formatC(x=cor(Ec$mm[Ec$trait=="2_3"], Ec$true[Ec$trait=="2_3"]), digits=3, flag="#"), sep="")))

temp <- ggplot() +
  annotate("rect", xmin=0.07, xmax=Inf, ymin=-Inf, ymax=-0.13, fill="#F2F2F2", color=NA) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Ec, aes(x=true, y=mm, color=trait), size=0.5) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(legend.justification=c(1,0), legend.position=c(1,0), legend.title=element_blank()) +
  theme(legend.text=element_text(size=4), legend.key.size=unit(0.5, "lines")) +
  theme(legend.background=element_blank()) +
  theme(axis.text=element_text(size=7)) +
  scale_x_continuous(limits=range(c(Ec$mm, Ec$true))) +
  scale_y_continuous(limits=range(c(Ec$mm, Ec$true))) +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  ylab("estimated") +
  coord_fixed()
ggsave(plot=temp,
       filename="output/simsel/png/Ecov.png",
       width=1.7,
       height=1.7,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/Ecov.svg",
       width=1.7,
       height=1.7,
       units="in",
       scale=4/3)

# plot E (residual) variances.
Ev <- vector()
for(i in 1:100){
  Ev <- rbind(Ev, data.frame(mm=diag(out[[i]]$vcov.mm$E),
                             true=diag(out[[i]]$vcov.sim$E),
                             trait=1:3,
                             sim=i))
}
Ev$trait <- as.factor(Ev$trait)
Ev$trait <- factor(Ev$trait, labels=c(paste("r(1) = ", formatC(x=cor(Ev$mm[Ev$trait==1], Ev$true[Ev$trait==1]), digits=3, flag="#"), sep=""),
                                      paste("r(2) = ", formatC(x=cor(Ev$mm[Ev$trait==2], Ev$true[Ev$trait==2]), digits=3, flag="#"), sep=""),
                                      paste("r(3) = ", formatC(x=cor(Ev$mm[Ev$trait==3], Ev$true[Ev$trait==3]), digits=3, flag="#"), sep="")))

temp <- ggplot() +
  annotate("rect", xmin=1.15, xmax=Inf, ymin=-Inf, ymax=0.82, fill="#F2F2F2", color=NA) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Ev, aes(x=true, y=mm, color=trait), size=0.5) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(legend.justification=c(1,0), legend.position=c(1,0), legend.title=element_blank()) +
  theme(legend.text=element_text(size=4), legend.key.size=unit(0.5, "lines")) +
  theme(legend.background=element_blank()) +
  theme(axis.text=element_text(size=7)) +
  scale_x_continuous(limits=range(c(Ev$mm, Ev$true))) +
  scale_y_continuous(limits=range(c(Ev$mm, Ev$true))) +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  ylab("estimated") +
  coord_fixed()
ggsave(plot=temp,
       filename="output/simsel/png/Evar.png",
       width=1.7,
       height=1.7,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/Evar.svg",
       width=1.7,
       height=1.7,
       units="in",
       scale=4/3)

#####

### sel1
#####
# plot Z.
Z <- vector()
for(i in 1:100){
  Z <- rbind(Z, data.frame(estimated=c(out[[i]]$sel1$mm$Z, out[[i]]$sel1$mm$Z.direct, out[[i]]$sel1$mm$Z.indirect),
                           realized=c(out[[i]]$sel1$sim$Z, out[[i]]$sel1$sim$Z.direct, out[[i]]$sel1$sim$Z.indirect),
                           expected=c(out[[i]]$sel1$true$Z, out[[i]]$sel1$true$Z.direct, out[[i]]$sel1$true$Z.indirect),
                           type=c(rep("Z", 3), rep("direct", 3), rep("indirect", 3)),
                           trait=rep(1:3,3),
                           sim=i))
}
Z$trait <- as.factor(Z$trait)

# Z
ZZ <- Z[Z$type=="Z", ]

# estimated/realized Z.
round(sapply(1:3, FUN=function(x) cor(ZZ$realized[ZZ$trait==x], ZZ$expected[ZZ$trait==x])), 3)
# 0.144 0.857 0.721

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=ZZ, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_x_continuous(limits=range(c(ZZ$realized, ZZ$expected))) +
  scale_y_continuous(limits=range(c(ZZ$realized, ZZ$expected))) +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  ylab("realized (Z)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel1_ZZr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_ZZr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)


# Z - direct
Zd <- Z[Z$type=="direct", ]

# Z estimated.
round(sapply(1:3, FUN=function(x) cor(Zd$estimated[Zd$trait==x], Zd$expected[Zd$trait==x])), 3)
#0.230 0.346 0.314

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Zd, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("estimated (Z, direct)") +
  coord_fixed(xlim=c(-0.5, 1.5))

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Zde.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Zde.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# Z realized.
round(sapply(1:3, FUN=function(x) cor(Zd$realized[Zd$trait==x], Zd$expected[Zd$trait==x])), 3)
# 0.069 0.611 0.449

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Zd, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("realized (Z, direct)") +
  coord_fixed(xlim=c(-0.5, 1.5))

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Zdr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Zdr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# Z - indirect
Zi <- Z[Z$type=="indirect", ]

# Z estimated
round(sapply(1:3, FUN=function(x) cor(Zi$estimated[Zi$trait==x], Zi$expected[Zi$trait==x])), 3)
# 0.190 0.213 0.315

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Zi, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("estimated (Z, indirect)") +
  coord_fixed(xlim=c(-0.5, 1.5))

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Zie.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Zie.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# Z realized
round(sapply(1:3, FUN=function(x) cor(Zi$realized[Zi$trait==x], Zi$expected[Zi$trait==x])), 3)
# 0.550 0.756 0.843

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Zi, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("realized (Z, indirect)") +
  coord_fixed(xlim=c(-0.5, 1.5))

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Zir.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Zir.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# plot beta.
B <- vector()
for(i in 1:100){
  B <- rbind(B, data.frame(estimated=out[[i]]$sel1$mm$B,
                           realized=out[[i]]$sel1$sim$B,
                           expected=out[[i]]$sel1$true$B,
                           trait=1:3,
                           sim=i))
}
B$trait <- as.factor(B$trait)

# beta from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(B$estimated[B$trait==x], B$expected[B$trait==x])), 3)
# 0.402 0.360 0.331

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=B, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab(expression("estimated ("~beta~")")) +
  coord_fixed(xlim=c(-0.5, 1.5))

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Be.png",
       width=4.6,
       height=4.6,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Be.svg",
       width=4.6,
       height=4.6,
       units="in",
       scale=4/3)

# beta from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(B$realized[B$trait==x], B$expected[B$trait==x])), 3)
# 0.069 0.611 0.449

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=B, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(B$realized, B$expected)), breaks=c(0,1)) +
  scale_y_continuous(limits=range(c(B$realized, B$expected))) +
  ylab(expression("realized ("~beta~")")) +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Br.png",
       width=4.6,
       height=4.6,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Br.svg",
       width=4.6,
       height=4.6,
       units="in",
       scale=4/3)




# plot S.
S <- vector()
for(i in 1:100){
  S <- rbind(S, data.frame(estimated=c(out[[i]]$sel1$mm$S, out[[i]]$sel1$mm$S.direct, out[[i]]$sel1$mm$S.indirect),
                           realized=c(out[[i]]$sel1$sim$S, out[[i]]$sel1$sim$S.direct, out[[i]]$sel1$sim$S.indirect),
                           expected=c(out[[i]]$sel1$true$S, out[[i]]$sel1$true$S.direct, out[[i]]$sel1$true$S.indirect),
                           type=c(rep("S", 3), rep("direct", 3), rep("indirect", 3)),
                           trait=rep(1:3,3),
                           sim=i))
}
S$trait <- as.factor(S$trait)

# S.
SS <- S[S$type=="S", ]

# S from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(SS$estimated[SS$trait==x], SS$expected[SS$trait==x])), 3)
# 0.095 0.460 0.308

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=SS, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(SS$estimated, SS$expected))) +
  scale_y_continuous(limits=range(c(SS$estimated, SS$expected))) +
  ylab("estimated (S)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel1_SSe.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_SSe.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(SS$realized[SS$trait==x], SS$expected[SS$trait==x])), 3)
# 0.522 0.795 0.605

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=SS, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(SS$realized, SS$expected))) +
  scale_y_continuous(limits=range(c(SS$realized, SS$expected))) +
  ylab("realized (S)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel1_SSr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_SSr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S - direct
Sd <- S[S$type=="direct", ]

# S from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(Sd$estimated[Sd$trait==x], Sd$expected[Sd$trait==x])), 3)
# 0.083 0.344 0.317

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Sd, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,2)) +
  ylab("estimated (S, direct)") +
  coord_fixed(xlim=c(-1,3))

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Sde.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Sde.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(Sd$realized[Sd$trait==x], Sd$expected[Sd$trait==x])), 3)
# 0.441 0.618 0.451

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Sd, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,2)) +
  ylab("realized (S, direct)") +
  coord_fixed(xlim=c(-1,3))

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Sdr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Sdr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S - indirect
Si <- S[S$type=="indirect", ]

# S from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(Si$estimated[Si$trait==x], Si$expected[Si$trait==x])), 3)
# 0.294 0.367 0.385

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Si, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("estimated (S, indirect)") +
  coord_fixed(xlim=c(-0.5, 1.4))

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Sie.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Sie.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(Si$realized[Si$trait==x], Si$expected[Si$trait==x])), 3)
# 0.591 0.837 0.879

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Si, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("realized (S, indirect)") +
  coord_fixed(xlim=c(-0.5, 1.4))

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Sir.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Sir.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)


# plot I.
I <- vector()
for(i in 1:100){
  I <- rbind(I, data.frame(estimated=c(out[[i]]$sel1$mm$I, out[[i]]$sel1$mm$I.direct, out[[i]]$sel1$mm$I.indirect),
                           realized=c(out[[i]]$sel1$sim$I, out[[i]]$sel1$sim$I.direct, out[[i]]$sel1$sim$I.indirect),
                           expected=c(out[[i]]$sel1$true$I, out[[i]]$sel1$true$I.direct, out[[i]]$sel1$true$I.indirect),
                           type=c(rep("I", 3), rep("direct", 3), rep("indirect", 3)),
                           trait=rep(1:3,3),
                           sim=i))
}
I$trait <- as.factor(I$trait)

# I
II <- I[I$type=="I", ]

# I from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(II$estimated[II$trait==x], II$expected[II$trait==x])), 3)
# 0.289 0.476 0.325

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=II, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(II$estimated, II$expected))) +
  scale_y_continuous(limits=range(c(II$estimated, II$expected))) +
  ylab("estimated (i)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Ie.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Ie.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(II$realized[II$trait==x], II$expected[II$trait==x])), 3)
# 0.155 0.791 0.603

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=II, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(II$realized, II$expected))) +
  scale_y_continuous(limits=range(c(II$realized, II$expected))) +
  ylab("realized (i)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Ir.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Ir.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I - direct
Id <- I[I$type=="direct", ]

# I from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(Id$estimated[Id$trait==x], Id$expected[Id$trait==x])), 3)
# 0.270 0.352 0.323

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Id, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,2)) +
  ylab("estimated (i, direct)") +
  coord_fixed(xlim=c(-1,3))

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Ide.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Ide.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(Id$realized[Id$trait==x], Id$expected[Id$trait==x])), 3)
# 0.154 0.615 0.449

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Id, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,2)) +
  ylab("realized (i, direct)") +
  coord_fixed(xlim=c(-1,3))

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Idr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Idr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I - indirect
Ii <- I[I$type=="indirect", ]

# I from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(Ii$estimated[Ii$trait==x], Ii$expected[Ii$trait==x])), 3)
# 0.294 0.373 0.394

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Ii, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("estimated (i, indirect)") +
  coord_fixed(xlim=c(-0.5, 1.4))

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Iie.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Iie.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(Ii$realized[Ii$trait==x], Ii$expected[Ii$trait==x])), 3)
# 0.593 0.830 0.877

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Ii, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("realized (i, indirect)") +
  coord_fixed(xlim=c(-0.5, 1.4))

ggsave(plot=temp,
       filename="output/simsel/png/sel1_Iir.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel1_Iir.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

#####

### sel2
#####
# plot Z.
Z <- vector()
for(i in 1:100){
  Z <- rbind(Z, data.frame(estimated=c(out[[i]]$sel2$mm$Z, out[[i]]$sel2$mm$Z.direct, out[[i]]$sel2$mm$Z.indirect),
                           realized=c(out[[i]]$sel2$sim$Z, out[[i]]$sel2$sim$Z.direct, out[[i]]$sel2$sim$Z.indirect),
                           expected=c(out[[i]]$sel2$true$Z, out[[i]]$sel2$true$Z.direct, out[[i]]$sel2$true$Z.indirect),
                           type=c(rep("Z", 3), rep("direct", 3), rep("indirect", 3)),
                           trait=rep(1:3,3),
                           sim=i))
}
Z$trait <- as.factor(Z$trait)

# Z
ZZ <- Z[Z$type=="Z", ]

# Z estimated/realized
round(sapply(1:3, FUN=function(x) cor(ZZ$realized[ZZ$trait==x], ZZ$expected[ZZ$trait==x])), 3)
# 0.609 0.692 0.730

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=ZZ, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_x_continuous(limits=range(c(ZZ$realized, ZZ$expected))) +
  scale_y_continuous(limits=range(c(ZZ$realized, ZZ$expected))) +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  ylab("realized (Z)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel2_ZZr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_ZZr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)


# Z - direct
Zd <- Z[Z$type=="direct", ]

# Z estimated
round(sapply(1:3, FUN=function(x) cor(Zd$estimated[Zd$trait==x], Zd$expected[Zd$trait==x])), 3)
#0.268 0.247 0.209

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Zd, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("estimated (Z, direct)") +
  coord_fixed(xlim=c(-0.5,1))

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Zde.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Zde.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# Z realized
round(sapply(1:3, FUN=function(x) cor(Zd$realized[Zd$trait==x], Zd$expected[Zd$trait==x])), 3)
# 0.310 0.394 0.458

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Zd, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("realized (Z, direct)") +
  coord_fixed(xlim=c(-0.5,1))

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Zdr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Zdr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# Z - indirect
Zi <- Z[Z$type=="indirect", ]

# Z estimated
round(sapply(1:3, FUN=function(x) cor(Zi$estimated[Zi$trait==x], Zi$expected[Zi$trait==x])), 3)
# 0.174 0.263 0.335

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Zi, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("estimated (Z, indirect)") +
  coord_fixed(xlim=c(0,1))

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Zie.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Zie.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# Z realized
round(sapply(1:3, FUN=function(x) cor(Zi$realized[Zi$trait==x], Zi$expected[Zi$trait==x])), 3)
# 0.738 0.735 0.779

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Zi, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("realized (Z, indirect)") +
  coord_fixed(xlim=c(0,1))

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Zir.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Zir.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# plot beta.
B <- vector()
for(i in 1:100){
  B <- rbind(B, data.frame(estimated=out[[i]]$sel2$mm$B,
                           realized=out[[i]]$sel2$sim$B,
                           expected=out[[i]]$sel2$true$B,
                           trait=1:3,
                           sim=i))
}
B$trait <- as.factor(B$trait)

# beta from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(B$estimated[B$trait==x], B$expected[B$trait==x])), 3)
# 0.314 0.223 0.200

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=B, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab(expression("estimated ("~beta~")")) +
  coord_fixed(xlim=c(-0.5, 1.0))

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Be.png",
       width=4.6,
       height=4.6,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Be.svg",
       width=4.6,
       height=4.6,
       units="in",
       scale=4/3)

# beta from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(B$realized[B$trait==x], B$expected[B$trait==x])), 3)
# 0.310 0.394 0.458

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=B, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(B$realized, B$expected)), breaks=c(0,1)) +
  scale_y_continuous(limits=range(c(B$realized, B$expected))) +
  ylab(expression("realized ("~beta~")")) +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Br.png",
       width=4.6,
       height=4.6,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Br.svg",
       width=4.6,
       height=4.6,
       units="in",
       scale=4/3)

# plot S.
S <- vector()
for(i in 1:100){
  S <- rbind(S, data.frame(estimated=c(out[[i]]$sel2$mm$S, out[[i]]$sel2$mm$S.direct, out[[i]]$sel2$mm$S.indirect),
                           realized=c(out[[i]]$sel2$sim$S, out[[i]]$sel2$sim$S.direct, out[[i]]$sel2$sim$S.indirect),
                           expected=c(out[[i]]$sel2$true$S, out[[i]]$sel2$true$S.direct, out[[i]]$sel2$true$S.indirect),
                           type=c(rep("S", 3), rep("direct", 3), rep("indirect", 3)),
                           trait=rep(1:3,3),
                           sim=i))
}
S$trait <- as.factor(S$trait)

# S.
SS <- S[S$type=="S", ]

# S from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(SS$estimated[SS$trait==x], SS$expected[SS$trait==x])), 3)
# 0.304 0.182 0.292

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=SS, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(SS$estimated, SS$expected))) +
  scale_y_continuous(limits=range(c(SS$estimated, SS$expected))) +
  ylab("estimated (S)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel2_SSe.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_SSe.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(SS$realized[SS$trait==x], SS$expected[SS$trait==x])), 3)
# 0.637 0.614 0.627

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=SS, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(SS$realized, SS$expected))) +
  scale_y_continuous(limits=range(c(SS$realized, SS$expected))) +
  ylab("realized (S)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel2_SSr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_SSr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S - direct
Sd <- S[S$type=="direct", ]

# S from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(Sd$estimated[Sd$trait==x], Sd$expected[Sd$trait==x])), 3)
# 0.204 0.026 0.194

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Sd, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,2)) +
  ylab("estimated (S, direct)") +
  coord_fixed(xlim=c(-1,2))

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Sde.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Sde.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(Sd$realized[Sd$trait==x], Sd$expected[Sd$trait==x])), 3)
# 0.450 0.434 0.448

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Sd, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,2)) +
  ylab("realized (S, direct)") +
  coord_fixed(xlim=c(-1,2))

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Sdr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Sdr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S - indirect
Si <- S[S$type=="indirect", ]

# S from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(Si$estimated[Si$trait==x], Si$expected[Si$trait==x])), 3)
# 0.375 0.422 0.333

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Si, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("estimated (S, indirect)") +
  coord_fixed(xlim=c(0, 1))

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Sie.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Sie.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(Si$realized[Si$trait==x], Si$expected[Si$trait==x])), 3)
# 0.759 0.827 0.839

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Si, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("realized (S, indirect)") +
  coord_fixed(xlim=c(0,1))

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Sir.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Sir.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)


# plot I.
I <- vector()
for(i in 1:100){
  I <- rbind(I, data.frame(estimated=c(out[[i]]$sel2$mm$I, out[[i]]$sel2$mm$I.direct, out[[i]]$sel2$mm$I.indirect),
                           realized=c(out[[i]]$sel2$sim$I, out[[i]]$sel2$sim$I.direct, out[[i]]$sel2$sim$I.indirect),
                           expected=c(out[[i]]$sel2$true$I, out[[i]]$sel2$true$I.direct, out[[i]]$sel2$true$I.indirect),
                           type=c(rep("I", 3), rep("direct", 3), rep("indirect", 3)),
                           trait=rep(1:3,3),
                           sim=i))
}
I$trait <- as.factor(I$trait)

# I
II <- I[I$type=="I", ]

# I from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(II$estimated[II$trait==x], II$expected[II$trait==x])), 3)
# 0.358 0.268 0.321

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=II, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(II$estimated, II$expected))) +
  scale_y_continuous(limits=range(c(II$estimated, II$expected))) +
  ylab("estimated (i)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Ie.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Ie.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(II$realized[II$trait==x], II$expected[II$trait==x])), 3)
# 0.554 0.581 0.637

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=II, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(II$realized, II$expected))) +
  scale_y_continuous(limits=range(c(II$realized, II$expected))) +
  ylab("realized (i)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Ir.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Ir.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I - direct
Id <- I[I$type=="direct", ]

# I from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(Id$estimated[Id$trait==x], Id$expected[Id$trait==x])), 3)
# 0.257 0.112 0.196

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Id, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("estimated (i, direct)") +
  coord_fixed(xlim=c(-0.5,1.5))

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Ide.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Ide.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(Id$realized[Id$trait==x], Id$expected[Id$trait==x])), 3)
# 0.360 0.405 0.453

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Id, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("realized (i, direct)") +
  coord_fixed(xlim=c(-0.5,1.5))

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Idr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Idr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I - indirect
Ii <- I[I$type=="indirect", ]

# I from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(Ii$estimated[Ii$trait==x], Ii$expected[Ii$trait==x])), 3)
# 0.351 0.421 0.334

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Ii, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("estimated (i, indirect)") +
  coord_fixed(xlim=c(0,1))

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Iie.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Iie.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(Ii$realized[Ii$trait==x], Ii$expected[Ii$trait==x])), 3)
# 0.753 0.822 0.831

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Ii, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("realized (i, indirect)") +
  coord_fixed(xlim=c(0,1))

ggsave(plot=temp,
       filename="output/simsel/png/sel2_Iir.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel2_Iir.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

#####

### sel3
#####
# plot Z.
Z <- vector()
for(i in 1:100){
  Z <- rbind(Z, data.frame(estimated=c(out[[i]]$sel3$mm$Z, out[[i]]$sel3$mm$Z.direct, out[[i]]$sel3$mm$Z.indirect),
                           realized=c(out[[i]]$sel3$sim$Z, out[[i]]$sel3$sim$Z.direct, out[[i]]$sel3$sim$Z.indirect),
                           expected=c(out[[i]]$sel3$true$Z, out[[i]]$sel3$true$Z.direct, out[[i]]$sel3$true$Z.indirect),
                           type=c(rep("Z", 3), rep("direct", 3), rep("indirect", 3)),
                           trait=rep(1:3,3),
                           sim=i))
}
Z$trait <- as.factor(Z$trait)

# Z
ZZ <- Z[Z$type=="Z", ]

# Z estimated/realized
round(sapply(1:3, FUN=function(x) cor(ZZ$realized[ZZ$trait==x], ZZ$expected[ZZ$trait==x])), 3)
# 0.546 0.612 0.448

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=ZZ, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_x_continuous(limits=range(c(ZZ$realized, ZZ$expected))) +
  scale_y_continuous(limits=range(c(ZZ$realized, ZZ$expected))) +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  ylab("realized (Z)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel3_ZZr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_ZZr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)


# Z - direct
Zd <- Z[Z$type=="direct", ]

# Z estimated.
round(sapply(1:3, FUN=function(x) cor(Zd$estimated[Zd$trait==x], Zd$expected[Zd$trait==x])), 3)
#0.087 0.205 0.041

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Zd, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("estimated (Z, direct)") +
  coord_fixed(xlim=c(0,1))

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Zde.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Zde.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# Z realized
round(sapply(1:3, FUN=function(x) cor(Zd$realized[Zd$trait==x], Zd$expected[Zd$trait==x])), 3)
# 0.229 0.359 0.181

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Zd, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("realized (Z, direct)") +
  coord_fixed(xlim=c(0,1))

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Zdr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Zdr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# Z - indirect
Zi <- Z[Z$type=="indirect", ]

# Z estimated
round(sapply(1:3, FUN=function(x) cor(Zi$estimated[Zi$trait==x], Zi$expected[Zi$trait==x])), 3)
# 0.122 0.101 0.355

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Zi, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("estimated (Z, indirect)") +
  coord_fixed(xlim=c(0,1))

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Zie.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Zie.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# Z realized
round(sapply(1:3, FUN=function(x) cor(Zi$realized[Zi$trait==x], Zi$expected[Zi$trait==x])), 3)
# 0.731 0.543 0.677

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Zi, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("realized (Z, indirect)") +
  coord_fixed(xlim=c(0,1))

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Zir.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Zir.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# plot beta.
B <- vector()
for(i in 1:100){
  B <- rbind(B, data.frame(estimated=out[[i]]$sel3$mm$B,
                           realized=out[[i]]$sel3$sim$B,
                           expected=out[[i]]$sel3$true$B,
                           trait=1:3,
                           sim=i))
}
B$trait <- as.factor(B$trait)

# beta from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(B$estimated[B$trait==x], B$expected[B$trait==x])), 3)
# 0.121 0.155 0.084

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=B, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab(expression("estimated ("~beta~")")) +
  coord_fixed(xlim=c(0, 1.0))

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Be.png",
       width=4.6,
       height=4.6,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Be.svg",
       width=4.6,
       height=4.6,
       units="in",
       scale=4/3)

# beta from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(B$realized[B$trait==x], B$expected[B$trait==x])), 3)
# 0.229 0.359 0.181

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=B, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(B$realized, B$expected))) +
  scale_y_continuous(limits=range(c(B$realized, B$expected))) +
  ylab(expression("realized ("~beta~")")) +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Br.png",
       width=4.6,
       height=4.6,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Br.svg",
       width=4.6,
       height=4.6,
       units="in",
       scale=4/3)

# plot S.
S <- vector()
for(i in 1:100){
  S <- rbind(S, data.frame(estimated=c(out[[i]]$sel3$mm$S, out[[i]]$sel3$mm$S.direct, out[[i]]$sel3$mm$S.indirect),
                           realized=c(out[[i]]$sel3$sim$S, out[[i]]$sel3$sim$S.direct, out[[i]]$sel3$sim$S.indirect),
                           expected=c(out[[i]]$sel3$true$S, out[[i]]$sel3$true$S.direct, out[[i]]$sel3$true$S.indirect),
                           type=c(rep("S", 3), rep("direct", 3), rep("indirect", 3)),
                           trait=rep(1:3,3),
                           sim=i))
}
S$trait <- as.factor(S$trait)

# S.
SS <- S[S$type=="S", ]

# S from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(SS$estimated[SS$trait==x], SS$expected[SS$trait==x])), 3)
# 0.119 0.197 0.030

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=SS, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(SS$estimated, SS$expected))) +
  scale_y_continuous(limits=range(c(SS$estimated, SS$expected))) +
  ylab("estimated (S)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel3_SSe.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_SSe.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(SS$realized[SS$trait==x], SS$expected[SS$trait==x])), 3)
# 0.521 0.559 0.342

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=SS, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(SS$realized, SS$expected))) +
  scale_y_continuous(limits=range(c(SS$realized, SS$expected))) +
  ylab("realized (S)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel3_SSr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_SSr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S - direct
Sd <- S[S$type=="direct", ]

# S from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(Sd$estimated[Sd$trait==x], Sd$expected[Sd$trait==x])), 3)
# 0.052  0.073 -0.015

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Sd, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("estimated (S, direct)") +
  coord_fixed(xlim=c(0,1.5))

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Sde.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Sde.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(Sd$realized[Sd$trait==x], Sd$expected[Sd$trait==x])), 3)
# 0.313 0.409 0.209

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Sd, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("realized (S, direct)") +
  coord_fixed(xlim=c(0,1.5))

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Sdr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Sdr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S - indirect
Si <- S[S$type=="indirect", ]

# S from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(Si$estimated[Si$trait==x], Si$expected[Si$trait==x])), 3)
# 0.390 0.287 0.309

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Si, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("estimated (S, indirect)") +
  coord_fixed(xlim=c(0, 1))

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Sie.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Sie.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# S from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(Si$realized[Si$trait==x], Si$expected[Si$trait==x])), 3)
# 0.810 0.694 0.740

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Si, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("realized (S, indirect)") +
  coord_fixed(xlim=c(0,1))

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Sir.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Sir.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)


# plot I.
I <- vector()
for(i in 1:100){
  I <- rbind(I, data.frame(estimated=c(out[[i]]$sel3$mm$I, out[[i]]$sel3$mm$I.direct, out[[i]]$sel3$mm$I.indirect),
                           realized=c(out[[i]]$sel3$sim$I, out[[i]]$sel3$sim$I.direct, out[[i]]$sel3$sim$I.indirect),
                           expected=c(out[[i]]$sel3$true$I, out[[i]]$sel3$true$I.direct, out[[i]]$sel3$true$I.indirect),
                           type=c(rep("I", 3), rep("direct", 3), rep("indirect", 3)),
                           trait=rep(1:3,3),
                           sim=i))
}
I$trait <- as.factor(I$trait)

# I
II <- I[I$type=="I", ]

# I from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(II$estimated[II$trait==x], II$expected[II$trait==x])), 3)
# 0.130 0.253 0.086

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=II, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(II$estimated, II$expected))) +
  scale_y_continuous(limits=range(c(II$estimated, II$expected))) +
  ylab("estimated (i)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Ie.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Ie.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(II$realized[II$trait==x], II$expected[II$trait==x])), 3)
# 0.450 0.514 0.280

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=II, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(limits=range(c(II$realized, II$expected))) +
  scale_y_continuous(limits=range(c(II$realized, II$expected))) +
  ylab("realized (i)") +
  coord_fixed()

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Ir.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Ir.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I - direct
Id <- I[I$type=="direct", ]

# I from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(Id$estimated[Id$trait==x], Id$expected[Id$trait==x])), 3)
# 0.076 0.107 0.029

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Id, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("estimated (i, direct)") +
  coord_fixed(xlim=c(0,1))

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Ide.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Ide.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(Id$realized[Id$trait==x], Id$expected[Id$trait==x])), 3)
# 0.257 0.378 0.182

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Id, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,1)) +
  ylab("realized (i, direct)") +
  coord_fixed(xlim=c(0,1))

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Idr.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Idr.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I - indirect
Ii <- I[I$type=="indirect", ]

# I from estimated (mixed model).
round(sapply(1:3, FUN=function(x) cor(Ii$estimated[Ii$trait==x], Ii$expected[Ii$trait==x])), 3)
# 0.380 0.279 0.314

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Ii, aes(x=expected, y=estimated, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,0.5)) +
  ylab("estimated (i, indirect)") +
  coord_fixed(xlim=c(0,0.5))

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Iie.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Iie.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

# I from realized (calculated using true G).
round(sapply(1:3, FUN=function(x) cor(Ii$realized[Ii$trait==x], Ii$expected[Ii$trait==x])), 3)
# 0.798 0.682 0.735

temp <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=Ii, aes(x=expected, y=realized, color=trait), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=7), legend.position="none") +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  scale_x_continuous(breaks=c(0,0.5)) +
  ylab("realized (i, indirect)") +
  coord_fixed(xlim=c(0,0.5))

ggsave(plot=temp,
       filename="output/simsel/png/sel3_Iir.png",
       width=2.3,
       height=2.3,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=temp,
       filename="output/simsel/svg/sel3_Iir.svg",
       width=2.3,
       height=2.3,
       units="in",
       scale=4/3)

#####


### check Z from BLUE vs BLUP.
simsel2 <- function(n.founder, gmap, sel){
  
  # set up the founder haplotypes.
  founder <- runMacs(nInd=n.founder,
                     nChr=length(gmap),
                     segSites=sapply(1:length(gmap), FUN=function(x) length(gmap[[x]])),
                     inbred=TRUE,
                     species="GENERIC",
                     ploidy=2L,
                     manualGenLen=sapply(1:length(gmap), FUN=function(x) max(gmap[[x]])))
  SP <- SimParam$new(founder)
  founder <- newPop(founder, simParam=SP)
  
  # get the allele count/frequency in the founder.
  af <- colSums(pullSegSiteGeno(founder, simParam=SP)==2)
  
  # identify 15 QTLs for each trait (3 traits total).
  # 9 unique QTLs, 2+2 QTLs shared between any two traits, 2 QTLs shared among 3 traits.
  qtl1 <- qtl2 <- qtl3 <- vector()
  meff1 <- meff2 <- meff3 <- rep(0, length(af))
  
  # largest effect QTL (0,0,1).
  temp <- which(af==1)
  idx <- sample(temp, 1, replace=FALSE)
  qtl1 <- c(qtl1, idx[1])
  qtl2 <- c(qtl2, idx[1])
  qtl3 <- c(qtl3, idx[1])
  meff1[idx[1]] <- meff2[idx[1]] <- meff3[idx[1]] <- 4
  
  # second largest QTL (0,1+1,0).
  temp <- which(af==2)
  idx <- sample(temp, 3, replace=FALSE)
  qtl1 <- c(qtl1, idx[c(1,2)])
  qtl2 <- c(qtl2, idx[c(1,3)])
  qtl3 <- c(qtl3, idx[c(2,3)])
  meff1[idx[c(1,2)]] <- 3
  meff2[idx[c(1,3)]] <- 3
  meff3[idx[c(2,3)]] <- 3
  
  # third largest QTL (3,0,1).
  temp <- which(af==3)
  idx <- sample(temp, 10, replace=FALSE)
  qtl1 <- c(qtl1, idx[c(1,2,3,4)])
  qtl2 <- c(qtl2, idx[c(1,5,6,7)])
  qtl3 <- c(qtl3, idx[c(1,8,9,10)])
  meff1[idx[c(1,2,3,4)]] <- 2
  meff2[idx[c(1,5,6,7)]] <- 2
  meff3[idx[c(1,8,9,10)]] <- 2
  
  # smallest QTL (6, 1+1, 0).
  temp <- which(af==4 | af==5)
  idx <- sample(temp, 21, replace=FALSE)
  qtl1 <- c(qtl1, idx[c(1,2,4,5,6,7,8,9)])
  qtl2 <- c(qtl2, idx[c(1,3,10,11,12,13,14,15)])
  qtl3 <- c(qtl3, idx[c(2,3,16,17,18,19,20,21)])
  meff1[idx[c(1,2,4,5,6,7,8,9)]] <- 1
  meff2[idx[c(1,3,10,11,12,13,14,15)]] <- 1
  meff3[idx[c(2,3,16,17,18,19,20,21)]] <- 1
  
  # make random crosses among the founders and create 300 DH individuals.
  pop0 <- randCross(pop=founder,
                    nCrosses=300,
                    nProgeny=1,
                    simParam=SP)
  pop0 <- makeDH(pop=pop0,
                 nDH=1,
                 simParam=SP)
  
  # get pop0 marker data.
  geno0 <- pullSegSiteGeno(pop0, simParam=SP)
  
  # scale the qtl effects such that Vg=1 in the 300 DH.
  meff1 <- meff1/sd(c(geno0%*%meff1))
  meff2 <- meff2/sd(c(geno0%*%meff2))
  meff3 <- meff3/sd(c(geno0%*%meff3))
  
  # calculate the genetic value in pop0.
  gv0 <- cbind(trait1=c(geno0%*%meff1),
               trait2=c(geno0%*%meff2),
               trait3=c(geno0%*%meff3))
  
  # add residual effect such that Ve=1.
  pv0 <- gv0 + matrix(rnorm(nrow(gv0)*ncol(gv0), 0, 1), nrow=nrow(gv0), ncol=ncol(gv0))
  
  # select on trait 1.
  idx1 <- which(rank(-pv0[,1]) <= floor(300*sel[1]))
  s1 <- colSums(pv0[idx1,])/length(idx1) - colSums(pv0)/nrow(pv0)
  
  # extract the selected individuals, make crosses and DH.
  # select on trait 1.
  pop1 <- selectInd(pop=pop0,
                    nInd=length(idx1),
                    use="rand",
                    candidates=idx1,
                    simParam=SP)
  pop1 <- randCross(pop=pop1,
                    nCrosses=300,
                    nProgeny=1,
                    simParam=SP)
  pop1 <- makeDH(pop=pop1,
                 nDH=1,
                 simParam=SP)

  # get pop1, pop2, pop3 marker data.
  geno1 <- pullSegSiteGeno(pop1, simParam=SP)

  # calculate the genetic value in pop1.
  gv1 <- cbind(trait1=c(geno1%*%meff1),
               trait2=c(geno1%*%meff2),
               trait3=c(geno1%*%meff3))
  
  # add residual effect such that Ve=1.
  pv1 <- gv1 + matrix(rnorm(nrow(gv1)*ncol(gv1), 0, 1), nrow=nrow(gv1), ncol=ncol(gv1))

  # estimate G and P for pop0 using mixed model.
  df0 <- data.frame(ID=rownames(geno0), pv0)
  vcov0 <- calc.vcov.sim(geno=geno0, pheno=df0, trait.idx=2:4)
  
  # estimate Z from BLUP.
  df1 <- data.frame(ID=rownames(geno1), pv1)
  mm.a1 <- tryCatch(sommer::mmer(trait1~1,
                                 random=~sommer::vs(at(year), ID, Gu=A.mat(rbind(geno0,geno1)-1)),
                                 rcov=~units,
                                 data=rbind(cbind(df0, year="A"), cbind(df1, year="B")),
                                 date.warning=FALSE,
                                 verbose=FALSE), error=function(e){})
  mm.a2 <- tryCatch(sommer::mmer(trait2~1,
                                 random=~sommer::vs(at(year), ID, Gu=A.mat(rbind(geno0,geno1)-1)),
                                 rcov=~units,
                                 data=rbind(cbind(df0, year="A"), cbind(df1, year="B")),
                                 date.warning=FALSE,
                                 verbose=FALSE), error=function(e){})
  mm.a3 <- tryCatch(sommer::mmer(trait3~1,
                                 random=~sommer::vs(at(year), ID, Gu=A.mat(rbind(geno0,geno1)-1)),
                                 rcov=~units,
                                 data=rbind(cbind(df0, year="A"), cbind(df1, year="B")),
                                 date.warning=FALSE,
                                 verbose=FALSE), error=function(e){})
  
  mm.b1 <- tryCatch(sommer::mmer(trait1~year,
                                random=~sommer::vs(ID, Gu=A.mat(rbind(geno0,geno1)-1)),
                                rcov=~units,
                                data=rbind(cbind(df0, year="A"), cbind(df1, year="B")),
                                date.warning=FALSE,
                                verbose=FALSE), error=function(e){})
  mm.b2 <- tryCatch(sommer::mmer(trait2~year,
                                 random=~sommer::vs(ID, Gu=A.mat(rbind(geno0,geno1)-1)),
                                 rcov=~units,
                                 data=rbind(cbind(df0, year="A"), cbind(df1, year="B")),
                                 date.warning=FALSE,
                                 verbose=FALSE), error=function(e){})
  mm.b3 <- tryCatch(sommer::mmer(trait3~year,
                                 random=~sommer::vs(ID, Gu=A.mat(rbind(geno0,geno1)-1)),
                                 rcov=~units,
                                 data=rbind(cbind(df0, year="A"), cbind(df1, year="B")),
                                 date.warning=FALSE,
                                 verbose=FALSE), error=function(e){})
  
  mm.c1 <- tryCatch(sommer::mmer(trait1~1,
                                random=~sommer::vs(ID, Gu=A.mat(rbind(geno0,geno1)-1)),
                                rcov=~units,
                                data=rbind(cbind(df0, year="A"), cbind(df1, year="B")),
                                date.warning=FALSE,
                                verbose=FALSE), error=function(e){})
  mm.c2 <- tryCatch(sommer::mmer(trait2~1,
                                 random=~sommer::vs(ID, Gu=A.mat(rbind(geno0,geno1)-1)),
                                 rcov=~units,
                                 data=rbind(cbind(df0, year="A"), cbind(df1, year="B")),
                                 date.warning=FALSE,
                                 verbose=FALSE), error=function(e){})
  mm.c3 <- tryCatch(sommer::mmer(trait3~1,
                                 random=~sommer::vs(ID, Gu=A.mat(rbind(geno0,geno1)-1)),
                                 rcov=~units,
                                 data=rbind(cbind(df0, year="A"), cbind(df1, year="B")),
                                 date.warning=FALSE,
                                 verbose=FALSE), error=function(e){})
  
  blup.a <- data.frame(A1=mm.a1$U$`A:ID`$trait1[df0$ID],
                       B1=mm.a1$U$`B:ID`$trait1[df1$ID],
                       A2=mm.a2$U$`A:ID`$trait2[df0$ID],
                       B2=mm.a2$U$`B:ID`$trait2[df1$ID],
                       A3=mm.a3$U$`A:ID`$trait3[df0$ID],
                       B3=mm.a3$U$`B:ID`$trait3[df1$ID])
  blup.b <- data.frame(A1=mm.b1$U$`u:ID`$trait1[df0$ID],
                       B1=mm.b1$U$`u:ID`$trait1[df1$ID],
                       A2=mm.b2$U$`u:ID`$trait2[df0$ID],
                       B2=mm.b2$U$`u:ID`$trait2[df1$ID],
                       A3=mm.b3$U$`u:ID`$trait3[df0$ID],
                       B3=mm.b3$U$`u:ID`$trait3[df1$ID])
  blup.c <- data.frame(A1=mm.c1$U$`u:ID`$trait1[df0$ID],
                       B1=mm.c1$U$`u:ID`$trait1[df1$ID],
                       A2=mm.c2$U$`u:ID`$trait2[df0$ID],
                       B2=mm.c2$U$`u:ID`$trait2[df1$ID],
                       A3=mm.c3$U$`u:ID`$trait3[df0$ID],
                       B3=mm.c3$U$`u:ID`$trait3[df1$ID])
  Z1.a <- c(mean(blup.a$B1) - mean(blup.a$A1),
            mean(blup.a$B2) - mean(blup.a$A2),
            mean(blup.a$B3) - mean(blup.a$A3))
  Z1.b <- c(mean(blup.b$B1) - mean(blup.b$A1),
            mean(blup.b$B2) - mean(blup.b$A2),
            mean(blup.b$B3) - mean(blup.b$A3))
  Z1.c <- c(mean(blup.c$B1) - mean(blup.c$A1),
            mean(blup.c$B2) - mean(blup.c$A2),
            mean(blup.c$B3) - mean(blup.c$A3))
  
  # estimate G and P for pop0 from actual var/cov of the genetic and phenotypic values.
  vcov0.d <- list(G=var(gv0), P=var(pv0))
  
  # calculate Z1.
  Z1 <- colSums(pv1)/nrow(pv1) - colSums(pv0)/nrow(pv0)
  
  # estimate beta.
  B1 <- solve(a=vcov0$G, b=Z1)
  B1.a <- solve(a=vcov0$G, b=Z1.a)
  B1.b <- solve(a=vcov0$G, b=Z1.b)
  B1.c <- solve(a=vcov0$G, b=Z1.c)
  B1.d <- solve(a=vcov0.d$G, b=Z1)
  
  # get the true selection intensity and beta.
  b1 <- solve(a=vcov0.d$P, b=s1)
  z1 <- c(vcov0.d$G%*%b1)
  
  
  # compile the outputs.
  out <- data.frame(Z=Z1, Z.a=Z1.a, Z.b=Z1.b, Z.c=Z1.c, z=z1, B=B1, B.a=B1.a, B.b=B1.b, B.c=B1.c, B.d=B1.d, b=b1)

  return(out)
}


out2 <- list()
for(i in 1:100){
  out2[[i]] <- simsel2(n.founder=10,
                       gmap=lapply(seq(1.0, 2.8, 0.2), FUN=function(x) seq(0, x-0.01, 0.01)),
                       sel=0.10)
  if(i%%10==0) message(i)
}

dat.plot.1 <- dat.plot.2 <- dat.plot.3 <- vector()
for(i in 1:100){
  dat.plot.1 <- rbind(dat.plot.1, out2[[i]][1,])
  dat.plot.2 <- rbind(dat.plot.2, out2[[i]][2,])
  dat.plot.3 <- rbind(dat.plot.3, out2[[i]][3,])
}

# Z, trait 1.
dat.plot.Z1 <- dat.plot.1[,1:4] - dat.plot.1[,5]
colnames(dat.plot.Z1) <- c("original", "blup_heterogeneous", "blup_fixed", "blup_none")
dat.plot.Z1 <- melt(dat.plot.Z1)
dat.plot.Z1$variable <- as.factor(dat.plot.Z1$variable)
dat.plot.Z1$variable <- factor(dat.plot.Z1$variable, levels=c("original", "blup_heterogeneous", "blup_fixed", "blup_none", "expected"))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, color="#FF0000") +
  geom_boxplot(data=dat.plot.Z1, aes(x=variable, y=value), outlier.shape=NA) +
  geom_jitter(data=dat.plot.Z1, aes(x=variable, y=value), color="#999999", height=0, width=0.25, size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=6)) +
  xlab("method") +
  ylab("diff") +
  labs(title="Z: trait 1")

ggsave(filename="output/simsel/sel1_diff_Z_trait1.png",
       width=4,
       height=3,
       units="in",
       dpi=600)

# Z, trait 2.
dat.plot.Z2 <- dat.plot.2[,1:4] - dat.plot.2[,5]
colnames(dat.plot.Z2) <- c("original", "blup_heterogeneous", "blup_fixed", "blup_none")
dat.plot.Z2 <- melt(dat.plot.Z2)
dat.plot.Z2$variable <- as.factor(dat.plot.Z2$variable)
dat.plot.Z2$variable <- factor(dat.plot.Z2$variable, levels=c("original", "blup_heterogeneous", "blup_fixed", "blup_none", "expected"))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, color="#FF0000") +
  geom_boxplot(data=dat.plot.Z2, aes(x=variable, y=value), outlier.shape=NA) +
  geom_jitter(data=dat.plot.Z2, aes(x=variable, y=value), color="#999999", height=0, width=0.25, size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=6)) +
  xlab("method") +
  ylab("diff") +
  labs(title="Z: trait 2")

ggsave(filename="output/simsel/sel1_diff_Z_trait2.png",
       width=4,
       height=3,
       units="in",
       dpi=600)

# Z, trait 3.
dat.plot.Z3 <- dat.plot.3[,1:4] - dat.plot.3[,5]
colnames(dat.plot.Z3) <- c("original", "blup_heterogeneous", "blup_fixed", "blup_none")
dat.plot.Z3 <- melt(dat.plot.Z3)
dat.plot.Z3$variable <- as.factor(dat.plot.Z3$variable)
dat.plot.Z3$variable <- factor(dat.plot.Z3$variable, levels=c("original", "blup_heterogeneous", "blup_fixed", "blup_none", "expected"))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, color="#FF0000") +
  geom_boxplot(data=dat.plot.Z3, aes(x=variable, y=value), outlier.shape=NA) +
  geom_jitter(data=dat.plot.Z3, aes(x=variable, y=value), color="#999999", height=0, width=0.25, size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=6)) +
  xlab("method") +
  ylab("diff") +
  labs(title="Z: trait 3")

ggsave(filename="output/simsel/sel1_diff_Z_trait3.png",
       width=4,
       height=3,
       units="in",
       dpi=600)


# beta, trait 1.
dat.plot.B1 <- dat.plot.1[,6:9] - dat.plot.1[,11]
colnames(dat.plot.B1) <- c("original", "blup_heterogeneous", "blup_fixed", "blup_none")
dat.plot.B1 <- melt(dat.plot.B1)
dat.plot.B1$variable <- as.factor(dat.plot.B1$variable)
dat.plot.B1$variable <- factor(dat.plot.B1$variable, levels=c("original", "blup_heterogeneous", "blup_fixed", "blup_none"))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, color="#FF0000") +
  geom_boxplot(data=dat.plot.B1, aes(x=variable, y=value), outlier.shape=NA) +
  geom_jitter(data=dat.plot.B1, aes(x=variable, y=value), color="#999999", height=0, width=0.25, size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=6)) +
  xlab("method") +
  ylab("diff") +
  coord_cartesian(ylim=c(min(dat.plot.B1$value), 5)) +
  labs(title="beta: trait 1")

ggsave(filename="output/simsel/sel1_diff_beta_trait1.png",
       width=4,
       height=3,
       units="in",
       dpi=600)

# beta, trait 2.
dat.plot.B2 <- dat.plot.2[,6:9] - dat.plot.2[,11]
colnames(dat.plot.B2) <- c("original", "blup_heterogeneous", "blup_fixed", "blup_none")
dat.plot.B2 <- melt(dat.plot.B2)
dat.plot.B2$variable <- as.factor(dat.plot.B2$variable)
dat.plot.B2$variable <- factor(dat.plot.B2$variable, levels=c("original", "blup_heterogeneous", "blup_fixed", "blup_none"))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, color="#FF0000") +
  geom_boxplot(data=dat.plot.B2, aes(x=variable, y=value), outlier.shape=NA) +
  geom_jitter(data=dat.plot.B2, aes(x=variable, y=value), color="#999999", height=0, width=0.25, size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=6)) +
  xlab("method") +
  ylab("diff") +
  coord_cartesian(ylim=c(-3.2, max(dat.plot.B2$value))) +
  labs(title="beta: trait 2")

ggsave(filename="output/simsel/sel1_diff_beta_trait2.png",
       width=4,
       height=3,
       units="in",
       dpi=600)

# beta, trait 3.
dat.plot.B3 <- dat.plot.3[,6:9] - dat.plot.3[,11]
colnames(dat.plot.B3) <- c("original", "blup_heterogeneous", "blup_fixed", "blup_none")
dat.plot.B3 <- melt(dat.plot.B3)
dat.plot.B3$variable <- as.factor(dat.plot.B3$variable)
dat.plot.B3$variable <- factor(dat.plot.B3$variable, levels=c("original", "blup_heterogeneous", "blup_fixed", "blup_none"))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, color="#FF0000") +
  geom_boxplot(data=dat.plot.B3, aes(x=variable, y=value), outlier.shape=NA) +
  geom_jitter(data=dat.plot.B3, aes(x=variable, y=value), color="#999999", height=0, width=0.25, size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text=element_text(size=6)) +
  xlab("method") +
  ylab("diff") +
  coord_cartesian(ylim=c(min(dat.plot.B3$value), 1.9)) +
  labs(title="beta: trait 3")

ggsave(filename="output/simsel/sel1_diff_beta_trait3.png",
       width=4,
       height=3,
       units="in",
       dpi=600)
