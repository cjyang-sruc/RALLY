library(ggplot2)
library(reshape2)
library(gridExtra)
library(sommer)

# set working directory.
setwd()


###################################
### PART 1 - read the GBS data. ###
###################################

# read the marker data in hapmap.
# sourced from Funmi.
# the figshare version seems bugged (https://doi.org/10.6084/m9.figshare.7350284.v1).
# because the part 1 and 2 are mixed together in figshare.
dat <- read.table("raw/TG_GBS_NoHetsAllimputed_344lines_42795SNPs.hmp.txt", as.is=TRUE, sep="\t", header=TRUE)

# create a data frame with marker information.
marker.gbs <- data.frame(chr=dat$chrom,
                         pos=dat$pos,
                         major=stringr::str_split_fixed(dat$alleles, "/", 2)[,1],
                         minor=stringr::str_split_fixed(dat$alleles, "/", 2)[,2])
keep <- which(!(marker.gbs$chr=="UN"))
marker.gbs <- marker.gbs[keep, ]

# we need to add this to the part2 positions.
# these positions are taken from doi.org/10.1186/s12864-018-5228-2
# http://crobiad.agwine.adelaide.edu.au/dawn/coord/
split.pos <- c(471304005, 438720154, 452179604,
               462376173, 453218924, 462216879,
               454103970, 448155269, 476235359,
               452555092, 451014251, 451004620,
               453230519, 451372872, 451901030,
               452440856, 452077197, 450509124,
               450046986, 453822637, 453812268)
temp <- unique(marker.gbs$chr)[seq(2,42,2)]
for(i in 1:length(temp)) marker.gbs$pos[marker.gbs$chr==temp[i]] <- marker.gbs$pos[marker.gbs$chr==temp[i]] + split.pos[i]
marker.gbs$chr <- gsub(pattern="_.*", replacement="", x=marker.gbs$chr)

# randomize marker coding (maf coding doesn't work with MLE delta control).
temp <- sample(c(0,1), nrow(marker.gbs), replace=TRUE)
marker.gbs$allele0 <- sapply(1:nrow(marker.gbs), FUN=function(x) marker.gbs[x, temp[x]+3])
marker.gbs$allele2 <- sapply(1:nrow(marker.gbs), FUN=function(x) marker.gbs[x, (1-temp[x])+3])

# create a matrix for genotype.
geno.gbs <- as.matrix(dat[keep, 12:355])
geno.gbs[geno.gbs==marker.gbs$allele0] <- 0
geno.gbs[geno.gbs==marker.gbs$allele2] <- 2
geno.gbs <- t(geno.gbs)
class(geno.gbs) <- "numeric"
rownames(geno.gbs) <- gsub("\\.", "_", rownames(geno.gbs))

# read the phenotypic data.
pheno.gbs <- read.table("raw/TG-344lines-phenotype_data.txt", as.is=TRUE, sep="\t", header=TRUE)
pheno.gbs$Taxa <- gsub("-", "_", pheno.gbs$Taxa)
rownames(pheno.gbs) <- pheno.gbs$Taxa
pheno.gbs <- pheno.gbs[rownames(geno.gbs), ]
rownames(pheno.gbs) <- NULL

# rename pheno.gbs column names and country names.
colnames(pheno.gbs)[c(1:2,9,11,13,14)] <- c("Variety", "Year", "WK", "SPWT", "EM2", "TILL")
pheno.gbs$Country[pheno.gbs$Country=="FRA"] <- "FR"
pheno.gbs$Country[pheno.gbs$Country=="DEU"] <- "DE"


####################################
### PART 2 - read the DArT data. ###
####################################

# read the marker data (sourced from Ian).
geno.dart <- read.csv("input/TRITICEAE_geno.csv", as.is=T, row.names=1)
geno.dart <- 2*as.matrix(geno.dart)

# read the marker information.
marker.dart <- read.csv("input/TRITICEAE_marker.csv", as.is=T)

# remove the marker names and replace them by 1,2,3,...,2012.
colnames(geno.dart) <- 1:ncol(geno.dart)
marker.dart <- marker.dart[,2:3]
colnames(marker.dart) <- c("chr", "pos")

# randomize marker coding (maf coding doesn't work with MLE delta control).
marker.dart$convert <- sample(c(0,1), nrow(marker.dart), replace=TRUE)
geno.dart[, marker.dart$convert==1] <- 2 - geno.dart[, marker.dart$convert==1]

# read the phenotype data.
pheno.dart <- read.csv("input/TRITICEAE_pheno.csv", as.is=T)
pheno.dart <- pheno.dart[, (1:15)]
colnames(pheno.dart)[c(4,8,12,14,15)] <- c("YLD", "EM2", "SPWT", "TILL", "WK")
pheno.dart <- pheno.dart[,colnames(pheno.gbs)]
pheno.dart$Country[pheno.dart$Country=="DEU"] <- "DE"
pheno.dart$Country[pheno.dart$Country=="GBR"] <- "UK"
pheno.dart$Country[pheno.dart$Country=="FRA"] <- "FR"

# replace the missing phenotype data in AWNS and WK with the values from pheno.gbs.
pheno.dart$Variety[is.na(pheno.dart$AWNS)] # "SPERBER"
pheno.dart$Variety[is.na(pheno.dart$WK)] # "MYTHOS"   "OPTIDOR"  "SOCRATES"
pheno.dart$AWNS[is.na(pheno.dart$AWNS)] <- 0.5
pheno.dart$WK[is.na(pheno.dart$WK)] <- 5.10

# read the major loci marker data.
mloci <- read.csv("input/TRITICEAE_cov.csv", as.is=T)


#########################################
### PART 3 - QC in GBS and DArT data. ###
#########################################

# identify varieties in common between GBS and DArT.
temp <- rownames(geno.dart)[rownames(geno.dart)%in%rownames(geno.gbs)]

# retain those 342 varieties in common.
geno.gbs <- geno.gbs[temp,]
geno.dart <- geno.dart[temp,]

rownames(pheno.gbs) <- pheno.gbs$Variety
pheno.gbs <- pheno.gbs[temp,]
rownames(pheno.gbs) <- NULL

rownames(pheno.dart) <- pheno.dart$Variety
pheno.dart <- pheno.dart[temp,]
rownames(pheno.dart) <- NULL

rownames(mloci) <- mloci$Variety
mloci <- mloci[temp,]
rownames(mloci) <- NULL

# identify varieties with contradictory year and country information.
temp <- which(!(pheno.dart$Year==pheno.gbs$Year) | !(pheno.dart$Country==pheno.gbs$Country))
cbind(pheno.dart[temp, 1:3], pheno.gbs[temp, 1:3])
#         Variety Year Country      Variety Year Country
#79     CAMP_REMY 1980      FR    CAMP_REMY 2006      FR
#80       CAMPARI 2003      DE      CAMPARI 1980      FR
#81       CAMPERO 2006      FR      CAMPERO 2003      DE
#116      DI_9714 2007      FR      DI_9714 2005      FR
#117      DINOSOR 2005      FR      DINOSOR 2007      FR
#213       MILVUS 2004      DE       MILVUS 2006      DE
#214 MIRONOVSKAJA 1963      DE MIRONOVSKAJA 2004      DE
#215      MITCHEL 2001      FR      MITCHEL 1963      DE
#216      MOISSON 1985      FR      MOISSON 2001      FR

# remove 9 varieties with contradictory information.
geno.gbs <- geno.gbs[-temp,]
pheno.gbs <- pheno.gbs[-temp,]
rownames(pheno.gbs) <- NULL

geno.dart <- geno.dart[-temp,]
pheno.dart <- pheno.dart[-temp,]
rownames(pheno.dart) <- NULL
mloci <- mloci[-temp,]
rownames(mloci) <- NULL

# calculate the number of markers that have r2 > 0.2 (linked) with each marker.
# we will discard markers with more linked markers in different chromosomes.
# gbs data.
qc.gbs <- vector()
for(i in 1:ncol(geno.gbs)){
  temp <- table(marker.gbs$chr[sapply(1:ncol(geno.gbs), FUN=function(x) cor(geno.gbs[,i], geno.gbs[,x], use="complete.obs")^2) > 0.2])
  qc.gbs <- rbind(qc.gbs, data.frame(chr=marker.gbs$chr[i],
                                     count=unname(temp[names(temp)==marker.gbs$chr[i]]),
                                     chr2=names(sort(temp, decreasing=TRUE))[1],
                                     count2=unname(sort(temp, decreasing=TRUE)[1])))
  if(i%%100 == 0) message(i)
}
saveRDS(qc.gbs, "output/qc_gbs.RDS")
table(qc.gbs$chr==qc.gbs$chr2)
#FALSE  TRUE 
# 3025 38836

# dart data.
qc.dart <- vector()
for(i in 1:ncol(geno.dart)){
  temp <- table(marker.dart$chr[sapply(1:ncol(geno.dart), FUN=function(x) cor(geno.dart[,i], geno.dart[,x], use="complete.obs")^2) > 0.2])
  qc.dart <- rbind(qc.dart, data.frame(chr=marker.dart$chr[i],
                                     count=unname(temp[names(temp)==marker.dart$chr[i]]),
                                     chr2=names(sort(temp, decreasing=TRUE))[1],
                                     count2=unname(sort(temp, decreasing=TRUE)[1])))
  if(i%%100 == 0) message(i)
}
saveRDS(qc.dart, "output/qc_dart.RDS")
table(qc.dart$chr==qc.dart$chr2)
#FALSE  TRUE 
#  102  1910

# remove 3,025 gbs markers that have higher LD with other markers in different chromosomes.
marker.gbs <- marker.gbs[qc.gbs$chr==qc.gbs$chr2, ]
geno.gbs <- geno.gbs[, qc.gbs$chr==qc.gbs$chr2]
rownames(marker.gbs) <- NULL
colnames(geno.gbs) <- 1:ncol(geno.gbs)

marker.dart <- marker.dart[qc.dart$chr==qc.dart$chr2, ]
geno.dart <- geno.dart[, qc.dart$chr==qc.dart$chr2]
rownames(marker.dart) <- NULL
colnames(geno.dart) <- 1:ncol(geno.dart)

# allele frequency in gbs.
af.gbs <- colSums(geno.gbs==2, na.rm=TRUE)/colSums(!is.na(geno.gbs))

# allele frequency in dart.
af.dart <- colSums(geno.dart==2, na.rm=TRUE)/colSums(!is.na(geno.dart))


##############################################
### PART 4 - Run RALLY and identify groups. ###
##############################################

# a function to find the groups within each chromosome
find.group <- function(rc, tsc, geno, min.r2=0.2, min.marker=10){
  
  # identify the significant markers.
  sig <- which(rc$score >= tsc)
  sig <- sig[order(-rc$score[sig])]
  
  # create an empty list to store the groups.
  g <- list()
  
  if(length(sig) > 0){
    
    # set the initial values for the while-loop.
    a <- 1:length(sig)
    b <- vector()
    
    # loop to find groups.
    while(length(a) > 0){
      
      # identify markers above LD threshold.
      temp <- which(sapply(1:ncol(geno), FUN=function(x) cor(geno[,sig[a[1]]], geno[,x], use="complete.obs")^2) >= min.r2)
      
      if(length(temp) >= min.marker){
        
        # set the most significant marker as the group name.
        temp <- list(temp)
        names(temp) <- temp[[1]][which.max(rc$score[temp[[1]]])]
        
        # remove significant markers that are assigned to group.
        a <- a[-which(sig[a]%in%temp[[1]])]
        if(length(b) > 0) b <- b[-which(sig[b]%in%temp[[1]])]
        
        # save the markers into group.
        g <- c(g, temp)

      } else {
        
        # move the marker that fails to form a group into b.
        a <- a[-1]
        b <- c(b, i)
        
      }
    }
    
  }
  
  return(g)
}

# a function to manually check the groups identified using find.group.
check.group <- function(group){
  df <- vector()
  for(i in 1:length(group)) df <- rbind(df, data.frame(idx=group[[i]], group=i, peak=c(group[[i]]==as.numeric(names(group)[i]))))
  plot1 <- ggplot() +
    geom_point(data=df, aes(x=idx, y=group), na.rm=TRUE) +
    geom_point(data=df[df$peak, ], aes(x=idx, y=group), color="#FF5555") +
    scale_y_continuous(breaks=1:10)
  print(plot1)
  return(c(unname(group), list(t(sapply(1:length(group), FUN=function(x) range(group[[x]]))))))
}

# a function to merge the groups after they have been manually checked.
merge.group <- function(gid, group, n){
  out <- rep(NA, n)
  for(i in 1:length(gid)) out[ unique(unlist(group[ gid[[i]] ])) ] <- i
  return(out)
}

# LL function for normal distribution pdf (to be used in GC + DC).
LL.norm <- function(p, x){
  a <- 1/(p[2]*sqrt(2*pi))
  L1 <- exp(-(x-p[1])^2/(2*p[2]^2))
  L2 <- exp(-(-x-p[1])^2/(2*p[2]^2))
  LL <- -( sum(log(L1/2 + L2/2)) + length(x)*log(a) )
  return(LL)
}

# RALLY in gbs.
rc.gbs <- vector()
afc.gbs <- vector()
for(i in 1:ncol(geno.gbs)){
  temp <- data.frame(allele=geno.gbs[,i]/2,
                     year=pheno.gbs$Year,
                     country=pheno.gbs$Country)
  temp.glm <- tryCatch(glm(allele ~ year + country, data=temp, family=binomial(link="logit")), error=function(e){})
  rc.gbs <- rbind(rc.gbs, summary(temp.glm)$coefficients[2,])
  
  t0 <- predict(object=temp.glm, newdata=data.frame(year=rep(1948,3), country=c("DE","FR","UK")), type="response")
  t1 <- predict(object=temp.glm, newdata=data.frame(year=rep(2007,3), country=c("DE","FR","UK")), type="response")
  afc.gbs <- c(afc.gbs, mean(abs(t1 - t0)))
  
  if(i%%100==0) message(i)
}

# arrange the results in a data frame.
rc.gbs <- data.frame(idx=1:nrow(marker.gbs), marker.gbs[,1:2], rc.gbs, afc=afc.gbs)
colnames(rc.gbs)[4:7] <- c("est", "se", "Z", "P")

# apply GC + DC correction with markers of afc < 0.20 as null
temp <- rc.gbs$Z[!is.na(rc.gbs$Z) & rc.gbs$afc < 0.20] # 59.4% of the markers.
temp <- nlm(f=LL.norm, p=c(mean(temp), sd(temp)), x=temp, print.level=0, gradtol=1e-36)
rc.gbs$Z.adj <- (rc.gbs$Z - temp$estimate[1])/temp$estimate[2]
rc.gbs$P.adj <- 2*pnorm(q=abs(rc.gbs$Z.adj), mean=0, sd=1, lower.tail=FALSE)

# convert the P-value into -log10p scores and get the Bonferroni threshold.
rc.gbs$score <- -log10(rc.gbs$P.adj)
tsc.gbs <- -log10(0.05/nrow(rc.gbs))

# create a vector to store the group information.
group.gbs <- vector()

# use find.group to identify the groups.
# use check.group to check the groups manually.
# use merge.group to merge the groups after they have been checked.
i <- rc.gbs$chr=="1A"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
check.group(group=temp)
temp <- merge.group(gid=list(c(1), c(2)), group=temp, n=sum(i))
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="1B"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
check.group(group=temp)
temp <- merge.group(gid=list(c(1,2,3)), group=temp, n=sum(i))
temp <- temp + max(group.gbs, na.rm=TRUE)
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="1D"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
temp <- rep(NA, sum(i))
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="2A"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
check.group(group=temp)
temp[[2]] <- temp[[2]][ temp[[2]] <= 496 ]
temp <- merge.group(gid=list(c(2), c(1)), group=temp, n=sum(i))
temp <- temp + max(group.gbs, na.rm=TRUE)
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="2B"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
check.group(group=temp)
temp <- merge.group(gid=list(c(1:8)), group=temp, n=sum(i))
temp <- temp + max(group.gbs, na.rm=TRUE)
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="2D"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
temp <- rep(NA, sum(i))
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="3A"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.gbs, na.rm=TRUE)
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="3B"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
check.group(group=temp)
temp[[1]] <- temp[[1]][ temp[[1]] >= 3281 ]
temp <- merge.group(gid=list(c(2), c(1)), group=temp, n=sum(i))
temp <- temp + max(group.gbs, na.rm=TRUE)
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="3D"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
temp <- rep(NA, sum(i))
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="4A"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.gbs, na.rm=TRUE)
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="4B"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.gbs, na.rm=TRUE)
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="4D"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
temp <- rep(NA, sum(i))
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="5A"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
check.group(group=temp)
temp <- merge.group(gid=list(c(1,2)), group=temp, n=sum(i))
temp <- temp + max(group.gbs, na.rm=TRUE)
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="5B"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.gbs, na.rm=TRUE)
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="5D"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.gbs, na.rm=TRUE)
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="6A"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
check.group(group=temp)
temp <- merge.group(gid=list(c(1), c(2), c(3)), group=temp, n=sum(i))
temp <- temp + max(group.gbs, na.rm=TRUE)
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="6B"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
temp <- rep(NA, sum(i))
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="6D"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
temp <- rep(NA, sum(i))
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="7A"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
check.group(group=temp)
temp <- merge.group(gid=list(c(2), c(1)), group=temp, n=sum(i))
temp <- temp + max(group.gbs, na.rm=TRUE)
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="7B"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
check.group(group=temp)
temp <- merge.group(gid=list(c(3), c(2), c(1)), group=temp, n=sum(i))
temp <- temp + max(group.gbs, na.rm=TRUE)
group.gbs <- c(group.gbs, temp)

i <- rc.gbs$chr=="7D"
temp <- find.group(rc=rc.gbs[i, ], tsc=tsc.gbs, geno=geno.gbs[, i], min.r2=0.2, min.marker=10)
temp <- rep(NA, sum(i))
group.gbs <- c(group.gbs, temp)

group.gbs[is.na(group.gbs)] <- 0
rc.gbs$group <- group.gbs


# RALLY in dart.
rc.dart <- vector()
afc.dart <- vector()
for(i in 1:ncol(geno.dart)){
  temp <- data.frame(allele=geno.dart[,i]/2,
                     year=pheno.dart$Year,
                     country=pheno.dart$Country)
  temp.glm <- tryCatch(glm(allele ~ year + country, data=temp, family=binomial(link="logit")), error=function(e){})
  rc.dart <- rbind(rc.dart, summary(temp.glm)$coefficients[2,])
  
  t0 <- predict(object=temp.glm, newdata=data.frame(year=rep(1948,3), country=c("DE","FR","UK")), type="response")
  t1 <- predict(object=temp.glm, newdata=data.frame(year=rep(2007,3), country=c("DE","FR","UK")), type="response")
  afc.dart <- c(afc.dart, mean(abs(t1 - t0)))
  
  if(i%%100==0) message(i)
}

# arrange the results in a data frame.
rc.dart <- data.frame(idx=1:nrow(marker.dart), marker.dart[,1:2], rc.dart, afc=afc.dart)
colnames(rc.dart)[4:7] <- c("est", "se", "Z", "P")

# apply GC + DC correction with markers of afc < 0.20 as null
temp <- rc.dart$Z[!is.na(rc.dart$Z) & rc.dart$afc < 0.20] # 65.8% of the markers.
temp <- nlm(f=LL.norm, p=c(0, 1), x=temp, print.level=0, gradtol=1e-36) # using the observed mean/sd leads to wrong convergence
rc.dart$Z.adj <- (rc.dart$Z - temp$estimate[1])/temp$estimate[2]
rc.dart$P.adj <- 2*pnorm(q=abs(rc.dart$Z.adj), mean=0, sd=1, lower.tail=FALSE)

# convert the P-value into -log10p scores and get the Bonferroni threshold.
rc.dart$score <- -log10(rc.dart$P)
tsc.dart <- -log10(0.05/nrow(rc.dart))

# create a vector to store the group information.
group.dart <- vector()

# use find.group to identify the groups.
# use check.group to check the groups manually.
# use merge.group to merge the groups after they have been checked.
i <- rc.dart$chr=="1A"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(2,3), c(1)), group=temp, n=sum(i))
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="1B"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="1D"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.dart, na.rm=TRUE)
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="2A"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(1,3), c(2)), group=temp, n=sum(i))
temp <- temp + max(group.dart, na.rm=TRUE)
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="2B"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(2), c(3), c(1), c(4)), group=temp, n=sum(i))
temp <- temp + max(group.dart, na.rm=TRUE)
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="2D"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="3A"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="3B"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.dart, na.rm=TRUE)
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="3D"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="4A"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="4B"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="4D"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="5A"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="5B"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="5D"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="6A"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(2), c(3), c(1)), group=temp, n=sum(i))
temp <- temp + max(group.dart, na.rm=TRUE)
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="6B"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.dart, na.rm=TRUE)
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="6D"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="7A"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.dart, na.rm=TRUE)
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="7B"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.dart <- c(group.dart, temp)

i <- rc.dart$chr=="7D"
temp <- find.group(rc=rc.dart[i, ], tsc=tsc.dart, geno=geno.dart[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.dart <- c(group.dart, temp)

group.dart[is.na(group.dart)] <- 0
rc.dart$group <- group.dart

# RALLY in mloci.
rc.mloci <- vector()
afc.mloci <- vector()
for(i in 2:ncol(mloci)){
  temp <- data.frame(allele=mloci[,i],
                     year=pheno.gbs$Year,
                     country=pheno.gbs$Country)
  temp.glm <- tryCatch(glm(allele ~ year + country, data=temp, family=binomial(link="logit")), error=function(e){})
  rc.mloci <- rbind(rc.mloci, summary(temp.glm)$coefficients[2,])
  
  t0 <- predict(object=temp.glm, newdata=data.frame(year=rep(1948,3), country=c("DE","FR","UK")), type="response")
  t1 <- predict(object=temp.glm, newdata=data.frame(year=rep(2007,3), country=c("DE","FR","UK")), type="response")
  afc.mloci <- c(afc.mloci, mean(abs(t1 - t0)))
}

# arrange the results in a data frame.
rc.mloci <- data.frame(locus=colnames(mloci)[2:14], rc.mloci, afc=afc.mloci)
colnames(rc.mloci)[2:5] <- c("est", "se", "Z", "P")
rc.mloci$score <- -log10(rc.mloci$P)


#################################
### PART 5 - Manhattan plots. ###
#################################

# create the manhattan plot for rc.gbs and highlight the groups.
plot.gbs <- list()
for(i in 1:21){
  chr <- unique(rc.gbs$chr)[i]
  n.group <- length(unique(rc.gbs$group[rc.gbs$chr==chr]))
  color.opt <- hcl(h=seq(15, 375, length=n.group), l=65, c=100)
  color.opt <- color.opt[-n.group]
  
  plot.gbs[[i]] <- ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
    annotate("segment", x=-Inf, xend=Inf, y=tsc.gbs, yend=tsc.gbs, color="#AAAAAA", linetype=2) +
    geom_point(data=rc.gbs[rc.gbs$chr==chr, ], aes(x=pos/1e6, y=score, color=as.factor(group)), size=1, alpha=0.8, stroke=0) +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    theme(axis.title=element_text(size=8), axis.text=element_text(size=8)) +
    theme(legend.title=element_text(size=8), legend.text=element_text(size=8)) +
    theme(legend.key.size=unit(0.75, "lines"), plot.subtitle=element_text(size=8)) +
    theme(legend.justification="top") +
    guides(color=guide_legend(override.aes=list(size=2))) +
    scale_color_manual(values=c("#555555", color.opt), name="group") +
    scale_y_continuous(limits=c(0, max(rc.gbs$score))) +
    scale_x_continuous(limits=c(0, max(rc.gbs$pos)/1e6)) +
    xlab("Pos (Mb)") +
    ylab(expression("-log"[10]*"P")) +
    labs(subtitle=chr)
}

ggsave(plot=arrangeGrob(grobs=plot.gbs, layout_matrix=matrix(1:21, nrow=7, ncol=3, byrow=TRUE)),
       filename="output/S5_rc_TG_gbs.png",
       width=7,
       height=9,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=arrangeGrob(grobs=plot.gbs, layout_matrix=matrix(1:21, nrow=7, ncol=3, byrow=TRUE)),
       filename="output/S5_rc_TG_gbs.svg",
       width=7,
       height=9,
       units="in",
       scale=4/3)

# create the manhattan plot for rc.dart and highlight the groups.
plot.dart <- list()
for(i in 1:21){
  chr <- unique(rc.dart$chr)[i]
  n.group <- length(unique(rc.dart$group[rc.dart$chr==chr]))
  color.opt <- hcl(h=seq(15, 375, length=n.group), l=65, c=100)
  color.opt <- color.opt[-n.group]
  
  plot.dart[[i]] <- ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
    annotate("segment", x=-Inf, xend=Inf, y=tsc.dart, yend=tsc.dart, color="#AAAAAA", linetype=2) +
    geom_point(data=rc.dart[rc.dart$chr==chr, ], aes(x=pos, y=score, color=as.factor(group)), size=1, alpha=0.8, stroke=0) +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    theme(axis.title=element_text(size=8), axis.text=element_text(size=8)) +
    theme(legend.title=element_text(size=8), legend.text=element_text(size=8)) +
    theme(legend.key.size=unit(0.75, "lines"), plot.subtitle=element_text(size=8)) +
    theme(legend.justification="top") +
    guides(color=guide_legend(override.aes=list(size=2))) +
    scale_color_manual(values=c("#555555", color.opt), name="group") +
    scale_y_continuous(limits=c(0, max(rc.dart$score))) +
    scale_x_continuous(limits=c(0, max(rc.dart$pos))) +
    xlab("Pos (cM)") +
    ylab(expression("-log"[10]*"P")) +
    labs(subtitle=chr)
}

ggsave(plot=arrangeGrob(grobs=plot.dart, layout_matrix=matrix(1:21, nrow=7, ncol=3, byrow=TRUE)),
       filename="output/S0_rc_TG_dart.png",
       width=7,
       height=9,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=arrangeGrob(grobs=plot.dart, layout_matrix=matrix(1:21, nrow=7, ncol=3, byrow=TRUE)),
       filename="output/S0_rc_TG_dart.svg",
       width=7,
       height=9,
       units="in",
       scale=4/3)


##################################################
### PART 6 - Peak, range and h2 of each group. ###
##################################################

# identify the chr, pos range and peak pos for each group.
group.gbs <- vector()
for(i in 1:22){
  temp <- which(rc.gbs$group==i)
  group.gbs <- rbind(group.gbs, 
                     c(i,
                       length(temp),
                       rc.gbs$chr[temp[1]],
                       rc.gbs$pos[temp][which.max(rc.gbs$score[temp])],
                       range(rc.gbs$pos[temp])))
}
group.gbs <- data.frame(group.gbs)
colnames(group.gbs) <- c("group", "count", "chr", "peak", "start", "end")

# add a row for background h2 in group.gbs.
group.gbs <- rbind(group.gbs, data.frame(group=0, count=NA, chr=NA, peak=NA, start=NA, end=NA))

# add columns for the h2.
temp <- matrix(NA, nrow=23, ncol=12)
colnames(temp) <- colnames(pheno.gbs)[4:15]
group.gbs <- cbind(group.gbs, temp)

# identify groups with h2 > 0.001 for each traits.
h2init.gbs <- replicate(12, vector())
for(i in 1:22){
  
  # calculate the kinship matrix for group i and the rest.
  A1 <- rrBLUP::A.mat(geno.gbs[, rc.gbs$group==i, drop=FALSE] - 1)
  A2 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group==i), drop=FALSE] - 1)
  
  # get h2 for group i in each trait.
  for(j in 1:12){
    dat <- data.frame(ID1=pheno.gbs$Variety,
                      ID2=pheno.gbs$Variety,
                      year=pheno.gbs$Year,
                      country=pheno.gbs$Country,
                      trait=pheno.gbs[, j+3])
    mm <- mmer(trait ~ year + country,
               random=~vs(ID1, Gu=A1) + vs(ID2, Gu=A2),
               rcov=~units,
               data=dat,
               date.warning=FALSE,
               verbose=FALSE)
    h2init.gbs[[j]] <- rbind(h2init.gbs[[j]], 
                             c(unlist(vpredict(mm, h1 ~ V1/(V1+V2+V3))),
                               unlist(vpredict(mm, h2 ~ V2/(V1+V2+V3)))))
  }
  message(i)
}


# FT
g <- which(h2init.gbs[[1]][,1] > 0.001)
length(g) # 8
A1 <- lapply(g, FUN=function(x) rrBLUP::A.mat(geno.gbs[, rc.gbs$group==x, drop=FALSE] - 1))
A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group%in%g), drop=FALSE] - 1)
temp <- data.frame(replicate(length(g) + 1, pheno.gbs$Variety),
                   year=pheno.gbs$Year,
                   country=pheno.gbs$Country,
                   trait=pheno.gbs$FT)
mm <- mmer(trait ~ year + country,
           random=~
             vs(X1, Gu=A1[[1]]) + vs(X2, Gu=A1[[2]]) + vs(X3, Gu=A1[[3]]) + vs(X4, Gu=A1[[4]]) +
             vs(X5, Gu=A1[[5]]) + vs(X6, Gu=A1[[6]]) + vs(X7, Gu=A1[[7]]) + vs(X8, Gu=A1[[8]]) +
             vs(X9, Gu=A0),
           rcov=~units,
           data=temp,
           date.warning=FALSE,
           verbose=FALSE)
temp <- summary(mm)$varcomp
group.gbs$FT[c(g,23)] <- temp[-nrow(temp),1]/sum(temp[,1])

# LODG
g <- which(h2init.gbs[[2]][,1] > 0.001)
length(g) # 8
A1 <- lapply(g, FUN=function(x) rrBLUP::A.mat(geno.gbs[, rc.gbs$group==x, drop=FALSE] - 1))
A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group%in%g), drop=FALSE] - 1)
temp <- data.frame(replicate(length(g) + 1, pheno.gbs$Variety),
                   year=pheno.gbs$Year,
                   country=pheno.gbs$Country,
                   trait=pheno.gbs$LODG)
mm <- mmer(trait ~ year + country,
           random=~
             vs(X1, Gu=A1[[1]]) + vs(X2, Gu=A1[[2]]) + vs(X3, Gu=A1[[3]]) + vs(X4, Gu=A1[[4]]) +
             vs(X5, Gu=A1[[5]]) + vs(X6, Gu=A1[[6]]) + vs(X7, Gu=A1[[7]]) + vs(X8, Gu=A1[[8]]) +
             vs(X9, Gu=A0),
           rcov=~units,
           data=temp,
           date.warning=FALSE,
           verbose=FALSE)
temp <- summary(mm)$varcomp
group.gbs$LODG[c(g,23)] <- temp[-nrow(temp),1]/sum(temp[,1])

# YLD
g <- which(h2init.gbs[[3]][,1] > 0.001)
length(g) # 15
A1 <- lapply(g, FUN=function(x) rrBLUP::A.mat(geno.gbs[, rc.gbs$group==x, drop=FALSE] - 1))
A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group%in%g), drop=FALSE] - 1)
temp <- data.frame(replicate(length(g) + 1, pheno.gbs$Variety),
                   year=pheno.gbs$Year,
                   country=pheno.gbs$Country,
                   trait=pheno.gbs$YLD)
mm <- mmer(trait ~ year + country,
           random=~
             vs(X1, Gu=A1[[1]]) + vs(X2, Gu=A1[[2]]) + vs(X3, Gu=A1[[3]]) + vs(X4, Gu=A1[[4]]) +
             vs(X5, Gu=A1[[5]]) + vs(X6, Gu=A1[[6]]) + vs(X7, Gu=A1[[7]]) + vs(X8, Gu=A1[[8]]) +
             vs(X9, Gu=A1[[9]]) + vs(X10, Gu=A1[[10]]) + vs(X11, Gu=A1[[11]]) + vs(X12, Gu=A1[[12]]) +
             vs(X13, Gu=A1[[13]]) + vs(X14, Gu=A1[[14]]) + vs(X15, Gu=A1[[15]]) +
             vs(X16, Gu=A0),
           rcov=~units,
           data=temp,
           date.warning=FALSE,
           verbose=FALSE)
temp <- summary(mm)$varcomp
group.gbs$YLD[c(g,23)] <- temp[-nrow(temp),1]/sum(temp[,1])

# HT
g <- which(h2init.gbs[[4]][,1] > 0.001)
length(g) # 9
A1 <- lapply(g, FUN=function(x) rrBLUP::A.mat(geno.gbs[, rc.gbs$group==x, drop=FALSE] - 1))
A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group%in%g), drop=FALSE] - 1)
temp <- data.frame(replicate(length(g) + 1, pheno.gbs$Variety),
                   year=pheno.gbs$Year,
                   country=pheno.gbs$Country,
                   trait=pheno.gbs$HT)
mm <- mmer(trait ~ year + country,
           random=~
             vs(X1, Gu=A1[[1]]) + vs(X2, Gu=A1[[2]]) + vs(X3, Gu=A1[[3]]) + vs(X4, Gu=A1[[4]]) +
             vs(X5, Gu=A1[[5]]) + vs(X6, Gu=A1[[6]]) + vs(X7, Gu=A1[[7]]) + vs(X8, Gu=A1[[8]]) +
             vs(X9, Gu=A1[[9]]) +
             vs(X10, Gu=A0),
           rcov=~units,
           data=temp,
           date.warning=FALSE,
           verbose=FALSE)
temp <- summary(mm)$varcomp
group.gbs$HT[c(g,23)] <- temp[-nrow(temp),1]/sum(temp[,1])

# PROT
g <- which(h2init.gbs[[5]][,1] > 0.001)
length(g) # 11
A1 <- lapply(g, FUN=function(x) rrBLUP::A.mat(geno.gbs[, rc.gbs$group==x, drop=FALSE] - 1))
A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group%in%g), drop=FALSE] - 1)
temp <- data.frame(replicate(length(g) + 1, pheno.gbs$Variety),
                   year=pheno.gbs$Year,
                   country=pheno.gbs$Country,
                   trait=pheno.gbs$PROT)
mm <- mmer(trait ~ year + country,
           random=~
             vs(X1, Gu=A1[[1]]) + vs(X2, Gu=A1[[2]]) + vs(X3, Gu=A1[[3]]) + vs(X4, Gu=A1[[4]]) +
             vs(X5, Gu=A1[[5]]) + vs(X6, Gu=A1[[6]]) + vs(X7, Gu=A1[[7]]) + vs(X8, Gu=A1[[8]]) +
             vs(X9, Gu=A1[[9]]) + vs(X10, Gu=A1[[10]]) + vs(X11, Gu=A1[[11]]) + 
             vs(X12, Gu=A0),
           rcov=~units,
           data=temp,
           date.warning=FALSE,
           verbose=FALSE)
temp <- summary(mm)$varcomp
group.gbs$PROT[c(g,23)] <- temp[-nrow(temp),1]/sum(temp[,1])

# WK
g <- which(h2init.gbs[[6]][,1] > 0.001)
length(g) # 10
A1 <- lapply(g, FUN=function(x) rrBLUP::A.mat(geno.gbs[, rc.gbs$group==x, drop=FALSE] - 1))
A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group%in%g), drop=FALSE] - 1)
temp <- data.frame(replicate(length(g) + 1, pheno.gbs$Variety),
                   year=pheno.gbs$Year,
                   country=pheno.gbs$Country,
                   trait=pheno.gbs$WK)
mm <- mmer(trait ~ year + country,
           random=~
             vs(X1, Gu=A1[[1]]) + vs(X2, Gu=A1[[2]]) + vs(X3, Gu=A1[[3]]) + vs(X4, Gu=A1[[4]]) +
             vs(X5, Gu=A1[[5]]) + vs(X6, Gu=A1[[6]]) + vs(X7, Gu=A1[[7]]) + vs(X8, Gu=A1[[8]]) +
             vs(X9, Gu=A1[[9]]) + vs(X10, Gu=A1[[10]]) +
             vs(X11, Gu=A0),
           rcov=~units,
           data=temp,
           date.warning=FALSE,
           verbose=FALSE)
temp <- summary(mm)$varcomp
group.gbs$WK[c(g,23)] <- temp[-nrow(temp),1]/sum(temp[,1])

# AWNS
g <- which(h2init.gbs[[7]][,1] > 0.001)
length(g) # 9
A1 <- lapply(g, FUN=function(x) rrBLUP::A.mat(geno.gbs[, rc.gbs$group==x, drop=FALSE] - 1))
A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group%in%g), drop=FALSE] - 1)
temp <- data.frame(replicate(length(g) + 1, pheno.gbs$Variety),
                   year=pheno.gbs$Year,
                   country=pheno.gbs$Country,
                   trait=pheno.gbs$AWNS)
mm <- mmer(trait ~ year + country,
           random=~
             vs(X1, Gu=A1[[1]]) + vs(X2, Gu=A1[[2]]) + vs(X3, Gu=A1[[3]]) + vs(X4, Gu=A1[[4]]) +
             vs(X5, Gu=A1[[5]]) + vs(X6, Gu=A1[[6]]) + vs(X7, Gu=A1[[7]]) + vs(X8, Gu=A1[[8]]) +
             vs(X9, Gu=A1[[9]]) +
             vs(X10, Gu=A0),
           rcov=~units,
           data=temp,
           date.warning=FALSE,
           verbose=FALSE)
temp <- summary(mm)$varcomp
group.gbs$AWNS[c(g,23)] <- temp[-nrow(temp),1]/sum(temp[,1])

# SPWT
g <- which(h2init.gbs[[8]][,1] > 0.001)
length(g) # 7
A1 <- lapply(g, FUN=function(x) rrBLUP::A.mat(geno.gbs[, rc.gbs$group==x, drop=FALSE] - 1))
A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group%in%g), drop=FALSE] - 1)
temp <- data.frame(replicate(length(g) + 1, pheno.gbs$Variety),
                   year=pheno.gbs$Year,
                   country=pheno.gbs$Country,
                   trait=pheno.gbs$SPWT)
mm <- mmer(trait ~ year + country,
           random=~
             vs(X1, Gu=A1[[1]]) + vs(X2, Gu=A1[[2]]) + vs(X3, Gu=A1[[3]]) + vs(X4, Gu=A1[[4]]) +
             vs(X5, Gu=A1[[5]]) + vs(X6, Gu=A1[[6]]) + vs(X7, Gu=A1[[7]]) +
             vs(X8, Gu=A0),
           rcov=~units,
           data=temp,
           date.warning=FALSE,
           verbose=FALSE)
temp <- summary(mm)$varcomp
group.gbs$SPWT[c(g,23)] <- temp[-nrow(temp),1]/sum(temp[,1])

# TGW
g <- which(h2init.gbs[[9]][,1] > 0.001)
length(g) # 13
A1 <- lapply(g, FUN=function(x) rrBLUP::A.mat(geno.gbs[, rc.gbs$group==x, drop=FALSE] - 1))
A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group%in%g), drop=FALSE] - 1)
temp <- data.frame(replicate(length(g) + 1, pheno.gbs$Variety),
                   year=pheno.gbs$Year,
                   country=pheno.gbs$Country,
                   trait=pheno.gbs$TGW)
mm <- mmer(trait ~ year + country,
           random=~
             vs(X1, Gu=A1[[1]]) + vs(X2, Gu=A1[[2]]) + vs(X3, Gu=A1[[3]]) + vs(X4, Gu=A1[[4]]) +
             vs(X5, Gu=A1[[5]]) + vs(X6, Gu=A1[[6]]) + vs(X7, Gu=A1[[7]]) + vs(X8, Gu=A1[[8]]) +
             vs(X9, Gu=A1[[9]]) + vs(X10, Gu=A1[[10]]) + vs(X11, Gu=A1[[11]]) + vs(X12, Gu=A1[[12]]) +
             vs(X13, Gu=A1[[13]]) +
             vs(X14, Gu=A0),
           rcov=~units,
           data=temp,
           date.warning=FALSE,
           verbose=FALSE)
temp <- summary(mm)$varcomp
group.gbs$TGW[c(g,23)] <- temp[-nrow(temp),1]/sum(temp[,1])

# EM2
g <- which(h2init.gbs[[10]][,1] > 0.001)
length(g) # 9
A1 <- lapply(g, FUN=function(x) rrBLUP::A.mat(geno.gbs[, rc.gbs$group==x, drop=FALSE] - 1))
A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group%in%g), drop=FALSE] - 1)
temp <- data.frame(replicate(length(g) + 1, pheno.gbs$Variety),
                   year=pheno.gbs$Year,
                   country=pheno.gbs$Country,
                   trait=pheno.gbs$EM2)
mm <- mmer(trait ~ year + country,
           random=~
             vs(X1, Gu=A1[[1]]) + vs(X2, Gu=A1[[2]]) + vs(X3, Gu=A1[[3]]) + vs(X4, Gu=A1[[4]]) +
             vs(X5, Gu=A1[[5]]) + vs(X6, Gu=A1[[6]]) + vs(X7, Gu=A1[[7]]) + vs(X8, Gu=A1[[8]]) +
             vs(X9, Gu=A1[[9]]) +
             vs(X10, Gu=A0),
           rcov=~units,
           data=temp,
           date.warning=FALSE,
           verbose=FALSE)
temp <- summary(mm)$varcomp
group.gbs$EM2[c(g,23)] <- temp[-nrow(temp),1]/sum(temp[,1])

# TILL
g <- which(h2init.gbs[[11]][,1] > 0.001)
length(g) # 8
A1 <- lapply(g, FUN=function(x) rrBLUP::A.mat(geno.gbs[, rc.gbs$group==x, drop=FALSE] - 1))
A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group%in%g), drop=FALSE] - 1)
temp <- data.frame(replicate(length(g) + 1, pheno.gbs$Variety),
                   year=pheno.gbs$Year,
                   country=pheno.gbs$Country,
                   trait=pheno.gbs$TILL)
mm <- mmer(trait ~ year + country,
           random=~
             vs(X1, Gu=A1[[1]]) + vs(X2, Gu=A1[[2]]) + vs(X3, Gu=A1[[3]]) + vs(X4, Gu=A1[[4]]) +
             vs(X5, Gu=A1[[5]]) + vs(X6, Gu=A1[[6]]) + vs(X7, Gu=A1[[7]]) + vs(X8, Gu=A1[[8]]) +
             vs(X9, Gu=A0),
           rcov=~units,
           data=temp,
           date.warning=FALSE,
           verbose=FALSE)
temp <- summary(mm)$varcomp
group.gbs$TILL[c(g,23)] <- temp[-nrow(temp),1]/sum(temp[,1])

# MAT
g <- which(h2init.gbs[[12]][,1] > 0.001)
length(g) # 6
A1 <- lapply(g, FUN=function(x) rrBLUP::A.mat(geno.gbs[, rc.gbs$group==x, drop=FALSE] - 1))
A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group%in%g), drop=FALSE] - 1)
temp <- data.frame(replicate(length(g) + 1, pheno.gbs$Variety),
                   year=pheno.gbs$Year,
                   country=pheno.gbs$Country,
                   trait=pheno.gbs$MAT)
mm <- mmer(trait ~ year + country,
           random=~
             vs(X1, Gu=A1[[1]]) + vs(X2, Gu=A1[[2]]) + vs(X3, Gu=A1[[3]]) + vs(X4, Gu=A1[[4]]) +
             vs(X5, Gu=A1[[5]]) + vs(X6, Gu=A1[[6]]) +
             vs(X7, Gu=A0),
           rcov=~units,
           data=temp,
           date.warning=FALSE,
           verbose=FALSE)
temp <- summary(mm)$varcomp
group.gbs$MAT[c(g,23)] <- temp[-nrow(temp),1]/sum(temp[,1])


# export the group.gbs h2 results.
write.csv(rc.gbs, "output/rc_TG_gbs.csv", quote=FALSE, row.names=FALSE)
write.csv(group.gbs, "output/rc_TG_gbs_h2.csv", quote=FALSE, row.names=FALSE)

# check the total h2 from the 22 groups.
h2total.gbs <- vector()
for(i in 1:12){
  temp <- which(group.gbs[1:22, 6+i] > 0.001)
  A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group%in%temp), drop=FALSE] - 1)
  A1 <- rrBLUP::A.mat(geno.gbs[, rc.gbs$group%in%temp, drop=FALSE] - 1)
  
  temp <- data.frame(ID0=pheno.gbs$Variety,
                     ID1=pheno.gbs$Variety,
                     year=pheno.gbs$Year,
                     country=pheno.gbs$Country,
                     trait=pheno.gbs[, i+3])
  mm <- mmer(trait ~ year + country,
             random=~vs(ID0, Gu=A0) + vs(ID1, Gu=A1),
             rcov=~units,
             data=temp,
             date.warning=FALSE,
             verbose=FALSE)
  h2total.gbs <- rbind(h2total.gbs, c(vpredict(mm, h0 ~ V1/(V1+V2+V3))$Estimate,
                                      vpredict(mm, h1 ~ V2/(V1+V2+V3))$Estimate))
}
h2total.gbs <- data.frame(trait=colnames(pheno.gbs)[4:15], h2total.gbs, colSums(group.gbs[23,7:18], na.rm=TRUE), colSums(group.gbs[-23,7:18], na.rm=TRUE))
colnames(h2total.gbs)[2:5] <- c("bg_total", "qtl_total", "bg", "qtl")

h2total.gbs$total <- NA
h2total.gbs <- h2total.gbs[,c(1,6,2:5)]
for(i in 1:12){
  A <- rrBLUP::A.mat(geno.gbs - 1)

  temp <- data.frame(ID=pheno.gbs$Variety,
                     year=pheno.gbs$Year,
                     country=pheno.gbs$Country,
                     trait=pheno.gbs[, i+3])
  mm <- mmer(trait ~ year + country,
             random=~vs(ID, Gu=A),
             rcov=~units,
             data=temp,
             date.warning=FALSE,
             verbose=FALSE)
  h2total.gbs$total[i] <- unname(c(vpredict(mm, h ~ V1/(V1+V2))$Estimate))
}

h2total.gbs
#     trait     total   bg_total  qtl_total         bg        qtl
#FT      FT 0.4748070 0.30378017 0.17762299 0.28359105 0.26038031
#LODG  LODG 0.2197346 0.11399181 0.21162018 0.12329775 0.19369669
#YLD    YLD 0.3125302 0.17394044 0.12226243 0.15967530 0.15255068
#HT      HT 0.4402769 0.27667331 0.31499249 0.22254937 0.40217027
#PROT  PROT 0.3855805 0.26330616 0.12212762 0.20475306 0.27484583
#WK      WK 0.3225396 0.27563015 0.06077987 0.21957710 0.13254774
#AWNS  AWNS 0.6631252 0.60436316 0.05574612 0.58344423 0.07417284
#SPWT  SPWT 0.1981242 0.16503363 0.05316908 0.15583548 0.12943033
#TGW    TGW 0.2117209 0.10618366 0.14710490 0.05016476 0.25402972
#EM2    EM2 0.3263935 0.29816782 0.02406646 0.21618382 0.15355681
#TILL  TILL 0.1268249 0.05706365 0.07692826 0.05137391 0.06844361
#MAT    MAT 0.3805576 0.24037949 0.13754340 0.17920244 0.19737515


save("geno.dart", "geno.gbs",
     "marker.dart", "marker.gbs",
     "pheno.dart", "pheno.gbs",
     "rc.dart", "rc.gbs",
     "tsc.dart", "tsc.gbs",
     "group.gbs", "h2total.gbs",
     "mloci", "rc.mloci",
     file="analysis_TG2.RData")
save("h2init.gbs", file="h2init.RData")


# plot the h2 for each group.
df <- data.frame(chr=group.gbs$chr,
                 marker=paste(group.gbs$group,
                              " (",
                              group.gbs$chr,
                              ", ",
                              round(as.numeric(group.gbs$peak)/1e6),
                              " Mb",
                              ")",
                              sep=""),
                 group.gbs[,7:18])
df$chr[23] <- "8"
df$marker[23] <- "others"
temp <- df$marker

df <- melt(df, id.vars=c("chr", "marker"))
df$marker <- as.factor(df$marker)
df$marker <- factor(df$marker, levels=temp)
df <- df[!is.na(df$value), ]



for(i in c(1,2,4,5,6)) group.gbs[,i] <- as.numeric(group.gbs[,i])

df <- melt(group.gbs[1:22, c(1,7:18)], id.vars="group")
df <- df[!is.na(df$value), ]
colnames(df)[2] <- "trait"

anno <- data.frame(marker=paste(group.gbs$chr[1:22], ", ", round(group.gbs$peak[1:22]/1e6), " Mb", sep=""),
                   group=1:22,
                   value=0.01 + rowSums(group.gbs[1:22, 7:18], na.rm=TRUE))

anno2 <- data.frame(a=0.5 + c(0, which(!(group.gbs$chr[-c(22,23)]==group.gbs$chr[-c(1,23)]))),
                    b=0.5 + c(which(!(group.gbs$chr[-c(22,23)]==group.gbs$chr[-c(1,23)])), 22),
                    d=rep(c("#FFFFFF", "#EEEEEE"), 7))

ggplot() +
  annotate("rect", xmin=anno2$a, xmax=anno2$b, ymin=-Inf, ymax=Inf, fill=anno2$d, color=NA) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_col(data=df, aes(x=group, y=value, fill=trait), width=0.8) +
  geom_text(data=anno, aes(x=group, y=value, label=marker), angle=90, hjust=0, size=2, color="#555555") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(legend.position="bottom", axis.text=element_text(size=8)) +
  theme(legend.title=element_text(size=8), legend.text=element_text(size=8)) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0), breaks=1:37) +
  scale_fill_brewer(palette="Paired") +
  ylab(expression("h"^"2")) +
  coord_cartesian(ylim=c(0, 0.5))

ggsave(filename="output/5A_groups_h2.png",
       height=3,
       width=1.75*3,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/5A_groups_h2.svg",
       height=3,
       width=1.75*3,
       units="in",
       scale=4/3)


df <- data.frame(trait=colnames(group.gbs)[7:18],
                 groups=colSums(group.gbs[1:22, 7:18], na.rm=TRUE),
                 others=colSums(group.gbs[23, 7:18], na.rm=TRUE))
df <- melt(df, id.vars="trait")
df$trait <- as.factor(df$trait)
df$trait <- factor(df$trait, levels=colnames(group.gbs)[7:18])
colnames(df)[2] <- "estimate"

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_col(data=df, aes(x=trait, y=value, fill=estimate), width=0.8, color="#555555") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(legend.position="bottom") +
  theme(axis.text.y=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8)) +
  theme(legend.title=element_text(size=8), legend.text=element_text(size=8)) +
  scale_fill_manual(values=c("#AAAAAA", "#EEEEEE")) +
  scale_y_continuous(expand=c(0,0)) +
  ylab(expression("h"^"2")) +
  coord_cartesian(ylim=c(0, 0.8))

ggsave(filename="output/5B_total_h2.png",
       height=3,
       width=1.75,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/5B_total_h2.svg",
       height=3,
       width=1.75,
       units="in",
       scale=4/3)


### LRT to test if the QSLs h2 are significantly different from 0.
load("analysis_TG2.RData")
load("h2init.RData")

### test individual QSL.
LRT.ind <- matrix(NA, nrow=12, ncol=22)
for(j in 1:22){
  
  # calculate the GRM for QSLs and non-QSL.
  A1 <- rrBLUP::A.mat(geno.gbs[, (rc.gbs$group==j), drop=FALSE] - 1)
  A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group==j), drop=FALSE] - 1)
  
  for(i in 1:12){
    
    # create a data.frame for use in mmer.
    temp <- data.frame(X1=pheno.gbs$Variety,
                       X0=pheno.gbs$Variety,
                       year=pheno.gbs$Year,
                       country=pheno.gbs$Country,
                       trait=pheno.gbs[,i+3])
    
    # mixed model fit with mmer (with QSL).
    mm1 <- mmer(trait ~ year + country,
                random=~vs(X1, Gu=A1) + vs(X0, Gu=A0),
                rcov=~units,
                data=temp,
                date.warning=FALSE,
                verbose=FALSE)
    
    # mixed model fit with mmer (without QSL).
    mm0 <- mmer(trait ~ year + country,
                random=~vs(X0, Gu=A0),
                rcov=~units,
                data=temp,
                date.warning=FALSE,
                verbose=FALSE)
    
    # get the likelihood and run LRT.
    LL1 <- summary(mm1)$logo$logLik
    LL0 <- summary(mm0)$logo$logLik
    LRT.ind[i,j] <- 2*(LL1 - LL0)
  }
}

### insufficient power to test individual QSL, so we will combine them all.
LRT <- vector()
for(i in 1:12){
  
  # identify the groups with h2 > 0.001.
  g <- which(h2init.gbs[[i]][,1] > 0.001)
  
  # calculate the GRM for QSLs and non-QSL.
  A1 <- rrBLUP::A.mat(geno.gbs[, (rc.gbs$group%in%g), drop=FALSE] - 1)
  A0 <- rrBLUP::A.mat(geno.gbs[, !(rc.gbs$group%in%g), drop=FALSE] - 1)
  
  # create a data.frame for use in mmer.
  temp <- data.frame(X1=pheno.gbs$Variety,
                     X0=pheno.gbs$Variety,
                     year=pheno.gbs$Year,
                     country=pheno.gbs$Country,
                     trait=pheno.gbs[,i+3])
  
  # mixed model fit with mmer (with QSL).
  mm1 <- mmer(trait ~ year + country,
              random=~vs(X1, Gu=A1) + vs(X0, Gu=A0),
              rcov=~units,
              data=temp,
              date.warning=FALSE,
              verbose=FALSE)
  
  # mixed model fit with mmer (without QSL).
  mm0 <- mmer(trait ~ year + country,
              random=~vs(X0, Gu=A0),
              rcov=~units,
              data=temp,
              date.warning=FALSE,
              verbose=FALSE)
  
  # get the likelihood and run LRT.
  LL1 <- summary(mm1)$logo$logLik
  LL0 <- summary(mm0)$logo$logLik
  LRT <- c(LRT, 2*(LL1 - LL0))
}

data.frame(trait=colnames(pheno.gbs)[4:15],
           chisq=LRT,
           p=pchisq(q=LRT, df=1, lower.tail=FALSE))
#   trait     chisq            p
#1     FT  5.438278 1.970009e-02
#2   LODG  9.668368 1.874678e-03
#3    YLD 15.506166 8.223656e-05
#4     HT 27.408367 1.647170e-07
#5   PROT 18.527613 1.674605e-05
#6     WK  4.192930 4.059288e-02
#7   AWNS  2.843832 9.172521e-02
#8   SPWT  4.181742 4.086169e-02
#9    TGW  8.067585 4.506413e-03
#10   EM2  2.422253 1.196228e-01
#11  TILL 10.498969 1.194412e-03
#12   MAT  7.358731 6.673789e-03

df <- data.frame(QSL=c(1:22,0),
                 pchisq(q=t(cbind(LRT.ind, LRT)), df=1, lower.tail=FALSE))
colnames(df)[2:13] <- colnames(pheno.gbs)[4:15]
rownames(df) <- NULL
write.csv(df, "output/S5_LRT_local_h2.csv", row.names=FALSE)
