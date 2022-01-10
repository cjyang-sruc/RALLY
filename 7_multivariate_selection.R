library(ggplot2)
library(sommer)
library(phytools)
library(ape)

# set the working directory.
setwd()

# load the required data.
load("analysis_TG2.RData")

# convert the geno from 0/1/2 to -1/0/1.
geno.dart2 <- geno.dart - 1
geno.gbs2 <- geno.gbs - 1

# create a matrix of trait pairs.
tpair <- t(combn(12,2))

# create a function to run mantel test and calculate the off-diagonal correlations.
mtest <- function(m1, m2, nperm){
  p <- ape::mantel.test(m1=m1, m2=m2, nperm=nperm)$p
  r <- cor(m1[lower.tri(m1, diag=FALSE)], m2[lower.tri(m1, diag=FALSE)], use="complete.obs")
  return(c(r=r, p=p))
}

### first, we check if Country makes any difference in the G-matrix (gbs data only).
# calculate the GRM for Country group.
A.DE <- rrBLUP::A.mat(X=geno.gbs2[pheno.gbs$Country=="DE", ])
A.FR <- rrBLUP::A.mat(X=geno.gbs2[pheno.gbs$Country=="FR", ])
A.UK <- rrBLUP::A.mat(X=geno.gbs2[pheno.gbs$Country=="UK", ])

# fit univariate mixed models for Country group.
uvar.DE <- list()
uvar.FR <- list()
uvar.UK <- list()
for(i in 1:12){
  temp <- pheno.gbs[pheno.gbs$Country=="DE", c(1:2, i+3)]
  colnames(temp)[3] <- "x"
  uvar.DE[[i]] <- tryCatch(summary(mmer(x~Year,
                                        random=~vs(Variety, Gu=A.DE),
                                        rcov=~units,
                                        data=temp,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
  temp <- pheno.gbs[pheno.gbs$Country=="FR", c(1:2, i+3)]
  colnames(temp)[3] <- "x"
  uvar.FR[[i]] <- tryCatch(summary(mmer(x~Year,
                                        random=~vs(Variety, Gu=A.FR),
                                        rcov=~units,
                                        data=temp,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
  temp <- pheno.gbs[pheno.gbs$Country=="UK", c(1:2, i+3)]
  colnames(temp)[3] <- "x"
  uvar.UK[[i]] <- tryCatch(summary(mmer(x~Year,
                                        random=~vs(Variety, Gu=A.UK),
                                        rcov=~units,
                                        data=temp,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
  message(i)
}

# check for convergence.
sum(sapply(1:12, FUN=function(x) is.null(uvar.DE[[x]]))) # 0
sum(sapply(1:12, FUN=function(x) is.null(uvar.FR[[x]]))) # 0
sum(sapply(1:12, FUN=function(x) is.null(uvar.UK[[x]]))) # 0

# fit bivariate mixed models for Country group.
mvar.DE <- list()
mvar.FR <- list()
mvar.UK <- list()
for(i in 1:66){
  temp <- pheno.gbs[pheno.gbs$Country=="DE", c(1:2, tpair[i,1]+3, tpair[i,2]+3)]
  colnames(temp)[3:4] <- c("x", "y")
  mvar.DE[[i]] <- tryCatch(summary(mmer(cbind(x,y)~Year,
                                        random=~vs(Variety, Gu=A.DE),
                                        rcov=~units,
                                        data=temp,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
  temp <- pheno.gbs[pheno.gbs$Country=="FR", c(1:2, tpair[i,1]+3, tpair[i,2]+3)]
  colnames(temp)[3:4] <- c("x", "y")
  mvar.FR[[i]] <- tryCatch(summary(mmer(cbind(x,y)~Year,
                                        random=~vs(Variety, Gu=A.FR),
                                        rcov=~units,
                                        data=temp,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
  temp <- pheno.gbs[pheno.gbs$Country=="UK", c(1:2, tpair[i,1]+3, tpair[i,2]+3)]
  colnames(temp)[3:4] <- c("x", "y")
  mvar.UK[[i]] <- tryCatch(summary(mmer(cbind(x,y)~Year,
                                        random=~vs(Variety, Gu=A.UK),
                                        rcov=~units,
                                        data=temp,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
  message(i)
}

# check for convergence.
sum(sapply(1:66, FUN=function(x) is.null(mvar.DE[[x]]))) # 10
sum(sapply(1:66, FUN=function(x) is.null(mvar.FR[[x]]))) # 0
sum(sapply(1:66, FUN=function(x) is.null(mvar.UK[[x]]))) # 12

# re-run the models that did not converge by using higher tolparinv.
redo <- which(sapply(1:66, FUN=function(x) is.null(mvar.DE[[x]])))
for(i in redo){
  temp <- pheno.gbs[pheno.gbs$Country=="DE", c(1:2, tpair[i,1]+3, tpair[i,2]+3)]
  colnames(temp)[3:4] <- c("x", "y")
  mvar.DE[[i]] <- tryCatch(summary(mmer(cbind(x,y)~Year,
                                        random=~vs(Variety, Gu=A.DE),
                                        rcov=~units,
                                        data=temp,
                                        tolparinv=1e-1,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
}

redo <- which(sapply(1:66, FUN=function(x) is.null(mvar.UK[[x]])))
for(i in redo){
  temp <- pheno.gbs[pheno.gbs$Country=="UK", c(1:2, tpair[i,1]+3, tpair[i,2]+3)]
  colnames(temp)[3:4] <- c("x", "y")
  mvar.UK[[i]] <- tryCatch(summary(mmer(cbind(x,y)~Year,
                                        random=~vs(Variety, Gu=A.UK),
                                        rcov=~units,
                                        data=temp,
                                        tolparinv=1e-1,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
}

# extract the Va and Ve for Country group.
Va.DE <- Ve.DE <- vector()
Va.FR <- Ve.FR <- vector()
Va.UK <- Ve.UK <- vector()
for(i in 1:12){
  Va.DE <- rbind(Va.DE, uvar.DE[[i]][1,])
  Ve.DE <- rbind(Ve.DE, uvar.DE[[i]][2,])
  Va.FR <- rbind(Va.FR, uvar.FR[[i]][1,])
  Ve.FR <- rbind(Ve.FR, uvar.FR[[i]][2,])
  Va.UK <- rbind(Va.UK, uvar.UK[[i]][1,])
  Ve.UK <- rbind(Ve.UK, uvar.UK[[i]][2,])
}

# extract the Ca and Ce for Country group.
Ca.DE <- Ce.DE <- vector()
Ca.FR <- Ce.FR <- vector()
Ca.UK <- Ce.UK <- vector()
for(i in 1:66){
  Ca.DE <- rbind(Ca.DE, mvar.DE[[i]][2,])
  Ce.DE <- rbind(Ce.DE, mvar.DE[[i]][5,])
  Ca.FR <- rbind(Ca.FR, mvar.FR[[i]][2,])
  Ce.FR <- rbind(Ce.FR, mvar.FR[[i]][5,])
  Ca.UK <- rbind(Ca.UK, mvar.UK[[i]][2,])
  Ce.UK <- rbind(Ce.UK, mvar.UK[[i]][5,])
}

# create the genetic (G) and phenotypic (P) variance-covariance matrices for Country group.
G.DE <- diag(Va.DE[,1])
G.FR <- diag(Va.FR[,1])
G.UK <- diag(Va.UK[,1])
P.DE <- diag(Va.DE[,1] + Ve.DE[,1])
P.FR <- diag(Va.FR[,1] + Ve.FR[,1])
P.UK <- diag(Va.UK[,1] + Ve.UK[,1])
rownames(G.DE) <- colnames(G.DE) <- rownames(P.DE) <- colnames(P.DE) <- colnames(pheno.gbs)[4:15]
rownames(G.FR) <- colnames(G.FR) <- rownames(P.FR) <- colnames(P.FR) <- colnames(pheno.gbs)[4:15]
rownames(G.UK) <- colnames(G.UK) <- rownames(P.UK) <- colnames(P.UK) <- colnames(pheno.gbs)[4:15]
for(i in 1:nrow(tpair)){
  G.DE[tpair[i,1], tpair[i,2]] <- G.DE[tpair[i,2], tpair[i,1]] <- Ca.DE[i,1]
  P.DE[tpair[i,1], tpair[i,2]] <- P.DE[tpair[i,2], tpair[i,1]] <- Ca.DE[i,1] + Ce.DE[i,1]
  G.FR[tpair[i,1], tpair[i,2]] <- G.FR[tpair[i,2], tpair[i,1]] <- Ca.FR[i,1]
  P.FR[tpair[i,1], tpair[i,2]] <- P.FR[tpair[i,2], tpair[i,1]] <- Ca.FR[i,1] + Ce.FR[i,1]
  G.UK[tpair[i,1], tpair[i,2]] <- G.UK[tpair[i,2], tpair[i,1]] <- Ca.UK[i,1]
  P.UK[tpair[i,1], tpair[i,2]] <- P.UK[tpair[i,2], tpair[i,1]] <- Ca.UK[i,1] + Ce.UK[i,1]
}

# exclude AWNS because DE has 0 variance.
G.DE <- G.DE[-7, -7]
G.FR <- G.FR[-7, -7]
G.UK <- G.UK[-7, -7]

g.DE <- cov2cor(G.DE)
g.FR <- cov2cor(G.FR)
g.UK <- cov2cor(G.UK)

# RS test.
skewers(X=G.DE, Y=G.FR, nsim=1000, method="unifcorrmat") # r=0.8874988, p=0.
skewers(X=G.DE, Y=G.UK, nsim=1000, method="unifcorrmat") # r=0.8910844, p=0.
skewers(X=G.FR, Y=G.UK, nsim=1000, method="unifcorrmat") # r=0.9390572, p=0.

skewers(X=g.DE, Y=g.FR, nsim=1000, method="unifcorrmat") # r=0.4561751, p=0.966.
skewers(X=g.DE, Y=g.UK, nsim=1000, method="unifcorrmat") # r=0.2707087, p=1.
skewers(X=g.FR, Y=g.UK, nsim=1000, method="unifcorrmat") # r=0.3676850, p=1.

# mantel test.
mtest(m1=G.DE, m2=G.FR, nperm=10000) #r=-0.48350355  p=0.00349965
mtest(m1=G.DE, m2=G.UK, nperm=10000) #r=-0.01255277  p=0.43245675
mtest(m1=G.FR, m2=G.UK, nperm=10000) #r=-0.20335360  p=0.11128890

mtest(m1=g.DE, m2=g.FR, nperm=10000) #r=0.22311130 p=0.10508950
mtest(m1=g.DE, m2=g.UK, nperm=10000) #r=0.01760614 p=0.91280872
mtest(m1=g.FR, m2=g.UK, nperm=10000) #r=0.10802260 p=0.42985700


### next, we check if Year makes any difference in the G-matrix.
# calculate the GRM for Year group.
y1 <- pheno.gbs$Year <= 1990
y2 <- pheno.gbs$Year > 1990 & pheno.gbs$Year < 2002
y3 <- pheno.gbs$Year >= 2002
A.y1 <- rrBLUP::A.mat(X=geno.gbs2[y1, ])
A.y2 <- rrBLUP::A.mat(X=geno.gbs2[y2, ])
A.y3 <- rrBLUP::A.mat(X=geno.gbs2[y3, ])

# fit univariate mixed models for Year group.
uvar.y1 <- list()
uvar.y2 <- list()
uvar.y3 <- list()
for(i in 1:12){
  temp <- pheno.gbs[y1, c(1:3, i+3)]
  colnames(temp)[4] <- "x"
  uvar.y1[[i]] <- tryCatch(summary(mmer(x~Year+Country,
                                        random=~vs(Variety, Gu=A.y1),
                                        rcov=~units,
                                        data=temp,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
  temp <- pheno.gbs[y2, c(1:3, i+3)]
  colnames(temp)[4] <- "x"
  uvar.y2[[i]] <- tryCatch(summary(mmer(x~Year+Country,
                                        random=~vs(Variety, Gu=A.y2),
                                        rcov=~units,
                                        data=temp,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
  temp <- pheno.gbs[y3, c(1:3, i+3)]
  colnames(temp)[4] <- "x"
  uvar.y3[[i]] <- tryCatch(summary(mmer(x~Year+Country,
                                        random=~vs(Variety, Gu=A.y3),
                                        rcov=~units,
                                        data=temp,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
  message(i)
}

# check for convergence.
sum(sapply(1:12, FUN=function(x) is.null(uvar.y1[[x]]))) # 0
sum(sapply(1:12, FUN=function(x) is.null(uvar.y2[[x]]))) # 0
sum(sapply(1:12, FUN=function(x) is.null(uvar.y3[[x]]))) # 0

# fit bivariate mixed models for Year group.
mvar.y1 <- list()
mvar.y2 <- list()
mvar.y3 <- list()
for(i in 1:66){
  temp <- pheno.gbs[y1, c(1:3, tpair[i,1]+3, tpair[i,2]+3)]
  colnames(temp)[4:5] <- c("x", "y")
  mvar.y1[[i]] <- tryCatch(summary(mmer(cbind(x,y)~Year+Country,
                                        random=~vs(Variety, Gu=A.y1),
                                        rcov=~units,
                                        data=temp,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
  temp <- pheno.gbs[y2, c(1:3, tpair[i,1]+3, tpair[i,2]+3)]
  colnames(temp)[4:5] <- c("x", "y")
  mvar.y2[[i]] <- tryCatch(summary(mmer(cbind(x,y)~Year+Country,
                                        random=~vs(Variety, Gu=A.y2),
                                        rcov=~units,
                                        data=temp,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
  temp <- pheno.gbs[y3, c(1:3, tpair[i,1]+3, tpair[i,2]+3)]
  colnames(temp)[4:5] <- c("x", "y")
  mvar.y3[[i]] <- tryCatch(summary(mmer(cbind(x,y)~Year+Country,
                                        random=~vs(Variety, Gu=A.y3),
                                        rcov=~units,
                                        data=temp,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
  message(i)
}

# check for convergence.
sum(sapply(1:66, FUN=function(x) is.null(mvar.y1[[x]]))) # 32
sum(sapply(1:65, FUN=function(x) is.null(mvar.y2[[x]]))) # 5
sum(sapply(1:66, FUN=function(x) is.null(mvar.y3[[x]]))) # 16

# re-run the models that did not converge by using higher tolparinv.
redo <- which(sapply(1:66, FUN=function(x) is.null(mvar.y1[[x]])))
for(i in redo){
  temp <- pheno.gbs[y1, c(1:3, tpair[i,1]+3, tpair[i,2]+3)]
  colnames(temp)[4:5] <- c("x", "y")
  mvar.y1[[i]] <- tryCatch(summary(mmer(cbind(x,y)~Year+Country,
                                        random=~vs(Variety, Gu=A.y1),
                                        rcov=~units,
                                        data=temp,
                                        tolparinv=1e-1,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
}

redo <- which(sapply(1:66, FUN=function(x) is.null(mvar.y2[[x]])))
for(i in redo){
  temp <- pheno.gbs[y2, c(1:3, tpair[i,1]+3, tpair[i,2]+3)]
  colnames(temp)[4:5] <- c("x", "y")
  mvar.y2[[i]] <- tryCatch(summary(mmer(cbind(x,y)~Year+Country,
                                        random=~vs(Variety, Gu=A.y2),
                                        rcov=~units,
                                        data=temp,
                                        tolparinv=1e-1,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
}

redo <- which(sapply(1:66, FUN=function(x) is.null(mvar.y3[[x]])))
for(i in redo){
  temp <- pheno.gbs[y3, c(1:3, tpair[i,1]+3, tpair[i,2]+3)]
  colnames(temp)[4:5] <- c("x", "y")
  mvar.y3[[i]] <- tryCatch(summary(mmer(cbind(x,y)~Year+Country,
                                        random=~vs(Variety, Gu=A.y3),
                                        rcov=~units,
                                        data=temp,
                                        tolparinv=1e-1,
                                        date.warning=FALSE,
                                        verbose=FALSE))$varcomp, error=function(e){})
}

# extract the Va and Ve for year group.
Va.y1 <- Ve.y1 <- vector()
Va.y2 <- Ve.y2 <- vector()
Va.y3 <- Ve.y3 <- vector()
for(i in 1:12){
  Va.y1 <- rbind(Va.y1, uvar.y1[[i]][1,])
  Ve.y1 <- rbind(Ve.y1, uvar.y1[[i]][2,])
  Va.y2 <- rbind(Va.y2, uvar.y2[[i]][1,])
  Ve.y2 <- rbind(Ve.y2, uvar.y2[[i]][2,])
  Va.y3 <- rbind(Va.y3, uvar.y3[[i]][1,])
  Ve.y3 <- rbind(Ve.y3, uvar.y3[[i]][2,])
}

# extract the Ca and Ce for year group.
Ca.y1 <- Ce.y1 <- vector()
Ca.y2 <- Ce.y2 <- vector()
Ca.y3 <- Ce.y3 <- vector()
for(i in 1:66){
  Ca.y1 <- rbind(Ca.y1, mvar.y1[[i]][2,])
  Ce.y1 <- rbind(Ce.y1, mvar.y1[[i]][5,])
  Ca.y2 <- rbind(Ca.y2, mvar.y2[[i]][2,])
  Ce.y2 <- rbind(Ce.y2, mvar.y2[[i]][5,])
  Ca.y3 <- rbind(Ca.y3, mvar.y3[[i]][2,])
  Ce.y3 <- rbind(Ce.y3, mvar.y3[[i]][5,])
}

# create the genetic (G) and phenotypic (P) variance-covariance matrices for year group.
G.y1 <- diag(Va.y1[,1])
G.y2 <- diag(Va.y2[,1])
G.y3 <- diag(Va.y3[,1])
P.y1 <- diag(Va.y1[,1] + Ve.y1[,1])
P.y2 <- diag(Va.y2[,1] + Ve.y2[,1])
P.y3 <- diag(Va.y3[,1] + Ve.y3[,1])
rownames(G.y1) <- colnames(G.y1) <- rownames(P.y1) <- colnames(P.y1) <- colnames(pheno.gbs)[4:15]
rownames(G.y2) <- colnames(G.y2) <- rownames(P.y2) <- colnames(P.y2) <- colnames(pheno.gbs)[4:15]
rownames(G.y3) <- colnames(G.y3) <- rownames(P.y3) <- colnames(P.y3) <- colnames(pheno.gbs)[4:15]
for(i in 1:nrow(tpair)){
  G.y1[tpair[i,1], tpair[i,2]] <- G.y1[tpair[i,2], tpair[i,1]] <- Ca.y1[i,1]
  P.y1[tpair[i,1], tpair[i,2]] <- P.y1[tpair[i,2], tpair[i,1]] <- Ca.y1[i,1] + Ce.y1[i,1]
  G.y2[tpair[i,1], tpair[i,2]] <- G.y2[tpair[i,2], tpair[i,1]] <- Ca.y2[i,1]
  P.y2[tpair[i,1], tpair[i,2]] <- P.y2[tpair[i,2], tpair[i,1]] <- Ca.y2[i,1] + Ce.y2[i,1]
  G.y3[tpair[i,1], tpair[i,2]] <- G.y3[tpair[i,2], tpair[i,1]] <- Ca.y3[i,1]
  P.y3[tpair[i,1], tpair[i,2]] <- P.y3[tpair[i,2], tpair[i,1]] <- Ca.y3[i,1] + Ce.y3[i,1]
}

# convert covariance to correlation.
g.y1 <- cov2cor(G.y1)
g.y2 <- cov2cor(G.y2)
g.y3 <- cov2cor(G.y3)

# RS test.
skewers(X=G.y1, Y=G.y2, nsim=1000, method="unifcorrmat") # r=0.8557984, p=0.
skewers(X=G.y1, Y=G.y3, nsim=1000, method="unifcorrmat") # r=0.9767594, p=0.
skewers(X=G.y2, Y=G.y3, nsim=1000, method="unifcorrmat") # r=0.8862480, p=0.

skewers(X=g.y1, Y=g.y2, nsim=1000, method="unifcorrmat") # r=0.4339233, p=0.982.
skewers(X=g.y1, Y=g.y3, nsim=1000, method="unifcorrmat") # r=0.4905931, p=0.912.
skewers(X=g.y2, Y=g.y3, nsim=1000, method="unifcorrmat") # r=0.4732028, p=0.950.

# mantel test.
mtest(m1=G.y1, m2=G.y2, nperm=10000) #r=-0.30336302  p=0.08069193
mtest(m1=G.y1, m2=G.y3, nperm=10000) #r= 0.43764525  p=0.01029897
mtest(m1=G.y2, m2=G.y3, nperm=10000) #r=-0.33100233  p=0.02929707

mtest(m1=g.y1, m2=g.y2, nperm=10000) #r=0.22321074 p=0.06809319 
mtest(m1=g.y1, m2=g.y3, nperm=10000) #r=0.43305613 p=0.00009999
mtest(m1=g.y2, m2=g.y3, nperm=10000) #r=0.43652084 p=0.00029997


### construct the G matrix.
# calculate the GRM.
Ad <- rrBLUP::A.mat(geno.dart2)
Ag <- rrBLUP::A.mat(geno.gbs2)

# fit univariate mixed models.
uvar.dart <- list()
uvar.gbs <- list()
for(i in 1:12){
  temp <- pheno.dart[, c(1:3, i+3)]
  colnames(temp)[4] <- "x"
  uvar.dart[[i]] <- tryCatch(summary(mmer(x~Year+Country,
                                         random=~vs(Variety, Gu=Ad),
                                         rcov=~units,
                                         data=temp,
                                         date.warning=FALSE,
                                         verbose=FALSE))$varcomp, error=function(e){})
  temp <- pheno.gbs[, c(1:3, i+3)]
  colnames(temp)[4] <- "x"
  uvar.gbs[[i]] <- tryCatch(summary(mmer(x~Year+Country,
                                         random=~vs(Variety, Gu=Ag),
                                         rcov=~units,
                                         data=temp,
                                         date.warning=FALSE,
                                         verbose=FALSE))$varcomp, error=function(e){})
  message(i)
}

# check for convergence.
sum(sapply(1:12, FUN=function(x) is.null(uvar.dart[[x]]))) # 0
sum(sapply(1:12, FUN=function(x) is.null(uvar.gbs[[x]]))) # 0

# fit bivariate mixed models for country group.
mvar.dart <- list()
mvar.gbs <- list()
for(i in 1:66){
  temp <- pheno.dart[, c(1:3, tpair[i,1]+3, tpair[i,2]+3)]
  colnames(temp)[4:5] <- c("x", "y")
  mvar.dart[[i]] <- tryCatch(summary(mmer(cbind(x,y)~Year+Country,
                                          random=~vs(Variety, Gu=Ad),
                                          rcov=~units,
                                          data=temp,
                                          date.warning=FALSE,
                                          verbose=FALSE))$varcomp, error=function(e){})
  temp <- pheno.gbs[, c(1:3, tpair[i,1]+3, tpair[i,2]+3)]
  colnames(temp)[4:5] <- c("x", "y")
  mvar.gbs[[i]] <- tryCatch(summary(mmer(cbind(x,y)~Year+Country,
                                         random=~vs(Variety, Gu=Ag),
                                         rcov=~units,
                                         data=temp,
                                         date.warning=FALSE,
                                         verbose=FALSE))$varcomp, error=function(e){})
  message(i)
}

# check for convergence.
sum(sapply(1:66, FUN=function(x) is.null(mvar.dart[[x]]))) # 0
sum(sapply(1:66, FUN=function(x) is.null(mvar.gbs[[x]]))) # 0

# extract the Va and Ve.
Va.dart <- Ve.dart <- vector()
Va.gbs <- Ve.gbs <- vector()
for(i in 1:12){
  Va.dart <- rbind(Va.dart, uvar.dart[[i]][1,])
  Ve.dart <- rbind(Ve.dart, uvar.dart[[i]][2,])
  Va.gbs <- rbind(Va.gbs, uvar.gbs[[i]][1,])
  Ve.gbs <- rbind(Ve.gbs, uvar.gbs[[i]][2,])
}

# extract the Ca and Ce.
Ca.dart <- Ce.dart <- vector()
Ca.gbs <- Ce.gbs <- vector()
for(i in 1:66){
  Ca.dart <- rbind(Ca.dart, mvar.dart[[i]][2,])
  Ce.dart <- rbind(Ce.dart, mvar.dart[[i]][5,])
  Ca.gbs <- rbind(Ca.gbs, mvar.gbs[[i]][2,])
  Ce.gbs <- rbind(Ce.gbs, mvar.gbs[[i]][5,])
}

# create the genetic (G) and phenotypic (P) variance-covariance matrices.
G.dart <- diag(Va.dart[,1])
P.dart <- diag(Va.dart[,1] + Ve.dart[,1])
rownames(G.dart) <- colnames(G.dart) <- rownames(P.dart) <- colnames(P.dart) <- colnames(pheno.dart)[4:15]
G.gbs <- diag(Va.gbs[,1])
P.gbs <- diag(Va.gbs[,1] + Ve.gbs[,1])
rownames(G.gbs) <- colnames(G.gbs) <- rownames(P.gbs) <- colnames(P.gbs) <- colnames(pheno.gbs)[4:15]
for(i in 1:nrow(tpair)){
  G.dart[tpair[i,1], tpair[i,2]] <- G.dart[tpair[i,2], tpair[i,1]] <- Ca.dart[i,1]
  P.dart[tpair[i,1], tpair[i,2]] <- P.dart[tpair[i,2], tpair[i,1]] <- Ca.dart[i,1] + Ce.dart[i,1]
  G.gbs[tpair[i,1], tpair[i,2]] <- G.gbs[tpair[i,2], tpair[i,1]] <- Ca.gbs[i,1]
  P.gbs[tpair[i,1], tpair[i,2]] <- P.gbs[tpair[i,2], tpair[i,1]] <- Ca.gbs[i,1] + Ce.gbs[i,1]
}

# create the genetic correlation matrix.
g.dart <- cov2cor(G.dart)
g.gbs <- cov2cor(G.gbs)

# compare G and g from DArT vs GBS.
skewers(X=G.dart, Y=G.gbs, nsim=1000, method="unifcorrmat") # r=0.9942361, p=0
skewers(X=g.dart, Y=g.gbs, nsim=1000, method="unifcorrmat") # r=0.9412567, p=0

mtest(m1=G.dart, m2=G.gbs, nperm=10000) #r=0.87017357 p=0.00009999
mtest(m1=g.dart, m2=g.gbs, nperm=10000) #r=0.89949989 p=0.00009999


# compare G and g to those from Country and Year groups.
# remove the traits that were previously excluded to compare.
# RS test.
skewers(X=G.DE, Y=G.gbs[-7,-7], nsim=1000, method="unifcorrmat") # r=0.9235453, p=0.
skewers(X=G.FR, Y=G.gbs[-7,-7], nsim=1000, method="unifcorrmat") # r=0.9793231, p=0.
skewers(X=G.UK, Y=G.gbs[-7,-7], nsim=1000, method="unifcorrmat") # r=0.9739909, p=0.

skewers(X=g.DE, Y=g.gbs[-7,-7], nsim=1000, method="unifcorrmat") # r=0.6249224, p=0.191.
skewers(X=g.FR, Y=g.gbs[-7,-7], nsim=1000, method="unifcorrmat") # r=0.8757477, p=0.
skewers(X=g.UK, Y=g.gbs[-7,-7], nsim=1000, method="unifcorrmat") # r=0.5096690, p=0.852.


skewers(X=G.y1, Y=G.gbs, nsim=1000, method="unifcorrmat") # r=0.9839012, p=0.
skewers(X=G.y2, Y=G.gbs, nsim=1000, method="unifcorrmat") # r=0.8751377, p=0.
skewers(X=G.y3, Y=G.gbs, nsim=1000, method="unifcorrmat") # r=0.9784891, p=0.

skewers(X=g.y1, Y=g.gbs, nsim=1000, method="unifcorrmat") # r=0.7613373, p=0.
skewers(X=g.y2, Y=g.gbs, nsim=1000, method="unifcorrmat") # r=0.7552862, p=0.
skewers(X=g.y3, Y=g.gbs, nsim=1000, method="unifcorrmat") # r=0.6598474, p=0.049.

# mantel test.
mtest(m1=G.DE, m2=G.gbs[-7,-7], nperm=10000) #r= 0.542879700 p=0.00289971
mtest(m1=G.FR, m2=G.gbs[-7,-7], nperm=10000) #r=-0.001044577 p=0.920407959
mtest(m1=G.UK, m2=G.gbs[-7,-7], nperm=10000) #r= 0.614620100 p=0.00189981

mtest(m1=g.DE, m2=g.gbs[-7,-7], nperm=10000) #r=0.44918196 p=0.00229977
mtest(m1=g.FR, m2=g.gbs[-7,-7], nperm=10000) #r=0.82488742 p=0.00009999
mtest(m1=g.UK, m2=g.gbs[-7,-7], nperm=10000) #r=0.28020946 p=0.03749625

mtest(m1=G.y1, m2=G.gbs, nperm=10000) #r=0.78819033 p=0.00009999
mtest(m1=G.y2, m2=G.gbs, nperm=10000) #r=0.19472470 p=0.1076892
mtest(m1=G.y3, m2=G.gbs, nperm=10000) #r=0.29957734 p=0.03809619

mtest(m1=g.y1, m2=g.gbs, nperm=10000) #r=0.70019733 p=0.00009999
mtest(m1=g.y2, m2=g.gbs, nperm=10000) #r=0.62352373 p=0.00009999
mtest(m1=g.y3, m2=g.gbs, nperm=10000) #r=0.69290139 p=0.00009999


### calculate selection gradient (Beta) from Z=GB.
# we will approximate Z from the slope of trait~year+country
Z.dart <- vector()
Z.gbs <- vector()
for(i in 1:12){
  temp <- data.frame(year=pheno.dart$Year, country=pheno.dart$Country, trait=pheno.dart[,i+3])
  Z.dart <- c(Z.dart, summary(lm(trait ~ year + country, data=temp))$coefficients[2,1])
  temp <- data.frame(year=pheno.gbs$Year, country=pheno.gbs$Country, trait=pheno.gbs[,i+3])
  Z.gbs <- c(Z.gbs, summary(lm(trait ~ year + country, data=temp))$coefficients[2,1])
}

B.dart <- solve(a=G.dart, b=Z.dart)
S.dart <- solve(a=G.dart%*%solve(P.dart), b=Z.dart)
I.dart <- S.dart/sqrt(diag(P.dart))

B.gbs <- solve(a=G.gbs, b=Z.gbs)
S.gbs <- solve(a=G.gbs%*%solve(P.gbs), b=Z.gbs)
I.gbs <- S.gbs/sqrt(diag(P.gbs))


round(cbind(B.dart,
            Z.dart, Z_dir=diag(G.dart)*B.dart, Z_ind=Z.dart-diag(G.dart)*B.dart,
            S.dart, S_dir=diag(P.dart)*B.dart, S_ind=S.dart-diag(P.dart)*B.dart,
            I.dart, I_dir=sqrt(diag(P.dart))*B.dart, I_ind=I.dart-sqrt(diag(P.dart))*B.dart), 3)
#     B.dart Z.dart  Z_dir  Z_ind S.dart  S_dir  S_ind I.dart  I_dir  I_ind
#FT   -0.578 -0.017 -1.383  1.366 -0.199 -2.824  2.626 -0.090 -1.277  1.188
#LODG  0.877 -0.033  0.085 -0.118  0.156  0.323 -0.167  0.257  0.532 -0.275
#YLD   0.233  0.365  1.221 -0.856  2.941  3.934 -0.993  0.716  0.957 -0.242
#HT   -0.024 -0.312 -0.249 -0.064 -0.071 -0.624  0.553 -0.014 -0.123  0.109
#PROT  1.176 -0.028  0.133 -0.161 -0.093  0.323 -0.416 -0.178  0.616 -0.794
#WK   -0.563 -0.004 -0.208  0.203 -0.410 -0.530  0.120 -0.422 -0.546  0.124
#AWNS  0.977  0.003  0.025 -0.022  0.001  0.040 -0.039  0.006  0.198 -0.192
#SPWT -0.210 -0.001 -0.100  0.099 -0.288 -0.622  0.335 -0.167 -0.362  0.194
#TGW   0.114 -0.017  0.077 -0.094  0.402  0.364  0.038  0.225  0.204  0.021
#EM2   0.023  0.101  4.211 -4.110 17.103 18.234 -1.131  0.604  0.644 -0.040
#TILL -2.433  0.001 -0.020  0.021 -0.118 -0.180  0.062 -0.435 -0.662  0.228
#MAT   2.029 -0.001  0.352 -0.353  0.340  0.874 -0.534  0.518  1.332 -0.814

round(cbind(B.gbs,
            Z.gbs, Z_dir=diag(G.gbs)*B.gbs, Z_ind=Z.gbs-diag(G.gbs)*B.gbs,
            S.gbs, S_dir=diag(P.gbs)*B.gbs, S_ind=S.gbs-diag(P.gbs)*B.gbs,
            I.gbs, I_dir=sqrt(diag(P.gbs))*B.gbs, I_ind=I.gbs-sqrt(diag(P.gbs))*B.gbs), 3)
#      B.gbs  Z.gbs  Z_dir  Z_ind  S.gbs  S_dir  S_ind  I.gbs  I_dir  I_ind
#FT   -0.155 -0.017 -0.368  0.351  0.721 -0.775  1.496  0.323 -0.347  0.670
#LODG -1.116 -0.033 -0.090  0.058 -0.335 -0.412  0.076 -0.552 -0.678  0.126
#YLD   0.071  0.352  0.341  0.011  2.761  1.092  1.669  0.702  0.278  0.424
#HT    0.071 -0.312  0.772 -1.084  0.509  1.753 -1.244  0.103  0.354 -0.251
#PROT -0.934 -0.028 -0.092  0.064 -0.252 -0.238 -0.014 -0.500 -0.472 -0.028
#WK   -0.317 -0.004 -0.094  0.089 -0.472 -0.290 -0.182 -0.493 -0.303 -0.190
#AWNS -0.591  0.003 -0.014  0.017 -0.031 -0.022 -0.010 -0.164 -0.113 -0.052
#SPWT  0.471 -0.001  0.266 -0.267  1.269  1.344 -0.075  0.752  0.796 -0.044
#TGW  -0.081 -0.017 -0.054  0.037 -0.061 -0.254  0.194 -0.034 -0.143  0.109
#EM2  -0.010  0.101 -2.529  2.630 -3.661 -7.747  4.086 -0.132 -0.279  0.147
#TILL  0.355  0.001  0.003 -0.002 -0.007  0.026 -0.033 -0.026  0.096 -0.122
#MAT   0.777 -0.001  0.125 -0.126  0.323  0.330 -0.006  0.496  0.506 -0.010


# display the covariance in lower half, and correlation in upper half.
temp <- G.dart
temp[upper.tri(temp, diag=FALSE)] <- g.dart[upper.tri(g.dart, diag=FALSE)]
round(temp, 3)
#         FT   LODG    YLD     HT   PROT     WK   AWNS   SPWT    TGW     EM2   TILL    MAT
#FT    2.394  0.226  0.074  0.409 -0.070 -0.510 -0.215 -0.254 -0.476  -0.039  0.059  0.867
#LODG  0.109  0.097 -0.198  0.596  0.128 -0.309  0.262  0.438 -0.085  -0.188  0.628  0.052
#YLD   0.262 -0.141  5.241 -0.314 -0.768  0.184  0.041 -0.420  0.054   0.033 -0.062 -0.009
#HT    2.034  0.597 -2.309 10.328  0.352 -0.491  0.006  0.455 -0.063  -0.112  0.086  0.198
#PROT -0.036  0.013 -0.591  0.381  0.113 -0.228  0.024  0.643  0.061  -0.181  0.151 -0.053
#WK   -0.479 -0.059  0.256 -0.959 -0.047  0.369  0.056 -0.134  0.121   0.275  0.224 -0.187
#AWNS -0.053  0.013  0.015  0.003  0.001  0.005  0.026  0.454  0.270   0.058  0.234 -0.403
#SPWT -0.271  0.094 -0.665  1.011  0.150 -0.056  0.051  0.477 -0.150   0.295 -0.135 -0.499
#TGW  -0.603 -0.022  0.100 -0.165  0.017  0.060  0.036 -0.085  0.669  -0.517  0.329 -0.433
#EM2  -0.826 -0.798  1.013 -4.916 -0.828  2.275  0.127  2.770 -5.755 185.220 -0.210 -0.130
#TILL  0.008  0.018 -0.013  0.025  0.005  0.012  0.003 -0.008  0.024  -0.257  0.008  0.176
#MAT   0.559  0.007 -0.009  0.266 -0.007 -0.047 -0.027 -0.144 -0.148  -0.738  0.007  0.174

temp <- G.gbs
temp[upper.tri(temp, diag=FALSE)] <- g.gbs[upper.tri(g.gbs, diag=FALSE)]
round(temp, 3)
#         FT   LODG    YLD     HT   PROT     WK   AWNS   SPWT    TGW     EM2   TILL    MAT
#FT    2.367  0.498 -0.077  0.388  0.053 -0.432 -0.284 -0.205 -0.492  -0.065 -0.009  0.808
#LODG  0.218  0.081 -0.205  0.723  0.048 -0.315  0.094  0.235 -0.522   0.072  0.415  0.102
#YLD  -0.260 -0.128  4.830 -0.461 -0.726  0.259  0.144 -0.214  0.190   0.138  0.081 -0.057
#HT    1.964  0.677 -3.333 10.826  0.395 -0.302 -0.116  0.305 -0.235  -0.042 -0.069 -0.033
#PROT  0.026  0.004 -0.500  0.408  0.098 -0.285 -0.079  0.498 -0.096  -0.202 -0.140 -0.047
#WK   -0.361 -0.049  0.309 -0.540 -0.049  0.295  0.020 -0.021  0.130   0.253 -0.111 -0.095
#AWNS -0.068  0.004  0.049 -0.059 -0.004  0.002  0.024  0.405  0.319  -0.093  0.380 -0.296
#SPWT -0.237  0.050 -0.353  0.755  0.117 -0.009  0.047  0.565 -0.293   0.310  0.023 -0.510
#TGW  -0.618 -0.122  0.342 -0.631 -0.025  0.058  0.041 -0.180  0.668  -0.529  0.301 -0.579
#EM2  -1.585  0.326  4.813 -2.199 -1.005  2.176 -0.228  3.686 -6.846 250.916 -0.034 -0.019
#TILL -0.001  0.011  0.017 -0.022 -0.004 -0.006  0.006  0.002  0.024  -0.051  0.009  0.302
#MAT   0.500  0.012 -0.050 -0.044 -0.006 -0.021 -0.019 -0.154 -0.190  -0.118  0.012  0.162


save.image("estimate_selection.RData")



