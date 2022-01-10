library(ggplot2)
library(reshape2)
library(gridExtra)

# set working directory.
setwd()

# read in the wheat refseq v1.0 physical position.
gff <- read.csv("input/wheat_refseq_v1_pmap.csv", as.is=T)

# read in the WAGTAIL marker data.
geno <- read.csv("input/WAGTAIL_geno.csv", as.is=TRUE, row.names=1)
geno <- as.matrix(geno)
marker <- read.csv("input/WAGTAIL_marker.csv", as.is=TRUE)

# get the physical positions for WAGTAIL markers.
marker$chr_p <- NA
marker$pos_p <- NA
for(i in 1:nrow(marker)){
  temp <- which(gff$marker==marker$Name[i])
  if(length(temp)==1){
    marker$chr_p[i] <- gff$chr[temp]
    marker$pos_p[i] <- gff$pos[temp]
  }
}

# keep only markers with matching chromosomes (17,841 to 5,592 markers).
temp <- which(marker$Chr==marker$chr_p)
marker <- marker[temp, ]
geno <- geno[, temp]

# remove the marker genetic positions and reorder.
marker <- marker[, c(1,4,5)]
colnames(marker) <- c("name", "chr", "pos")
temp <- order(marker$chr, marker$pos)
marker <- marker[temp, ]
geno <- geno[, temp]

# check for fixed markers (none).
af <- colSums(geno==2, na.rm=TRUE)/colSums(!is.na(geno))
temp <- which(af==0 | af==1)

# read the phenotypic data.
pheno <- read.csv("input/WAGTAIL_pheno.csv", as.is=TRUE)

# calculate the number of markers that have r2 > 0.2 (linked) with each marker.
# we will discard markers with more linked markers in different chromosomes.
qc <- vector()
for(i in 1:ncol(geno)){
  temp <- table(marker$chr[sapply(1:ncol(geno), FUN=function(x) cor(geno[,i], geno[,x], use="complete.obs")^2) > 0.2])
  qc <- rbind(qc, data.frame(chr=marker$chr[i],
                             count=unname(temp[names(temp)==marker$chr[i]]),
                             chr2=names(sort(temp, decreasing=TRUE))[1],
                             count2=unname(sort(temp, decreasing=TRUE)[1])))
  if(i%%100 == 0) message(i)
}
saveRDS(qc, "output/WT_qc.RDS")
table(qc$chr==qc$chr2)

# remove 319 markers that have higher LD with other markers in different chromosomes.
marker <- marker[qc$chr==qc$chr2, ]
geno <- geno[, qc$chr==qc$chr2]
rownames(marker) <- NULL

# LL function for normal distribution pdf (to be used in GC + DC).
LL.norm <- function(p, x){
  a <- 1/(p[2]*sqrt(2*pi))
  L1 <- exp(-(x-p[1])^2/(2*p[2]^2))
  L2 <- exp(-(-x-p[1])^2/(2*p[2]^2))
  LL <- -( sum(log(L1/2 + L2/2)) + length(x)*log(a) )
  return(LL)
}

# calculate the allele regression coefficients.
rc <- vector()
afc <- vector()
for(i in 1:ncol(geno)){
  temp <- data.frame(allele=geno[,i]/2,
                     year=pheno$Year,
                     country=pheno$Country)
  temp.glm <- tryCatch(glm(allele ~ year + country, data=temp, family=binomial(link="logit")), error=function(e){})
  rc <- rbind(rc, summary(temp.glm)$coefficients[2,])
  
  # skip the other countries.
  t0 <- predict(object=temp.glm, newdata=data.frame(year=rep(1916,3), country=c("DE","FR","UK")), type="response")
  t1 <- predict(object=temp.glm, newdata=data.frame(year=rep(2010,3), country=c("DE","FR","UK")), type="response")
  afc <- c(afc, mean(abs(t1 - t0)))
  
  if(i%%100==0) message(i)
}

# arrange the results in a data frame.
rc <- data.frame(idx=1:nrow(marker), marker[,2:3], rc, afc=afc)
colnames(rc)[4:7] <- c("est", "se", "Z", "P")

# apply GC + DC correction with markers of afc < 0.20 as null
temp <- rc$Z[!is.na(rc$Z) & rc$afc < 0.20] # 46.1% of the markers.
temp <- nlm(f=LL.norm, p=c(0, 1), x=temp, print.level=0, gradtol=1e-36)
rc$Z.adj <- (rc$Z - temp$estimate[1])/temp$estimate[2]
rc$P.adj <- 2*pnorm(q=abs(rc$Z.adj), mean=0, sd=1, lower.tail=FALSE)

# convert the P-value into -log10p scores and get the Bonferroni threshold.
rc$score <- -log10(rc$P.adj)
tsc <- -log10(0.05/nrow(rc))

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

# create a vector to store the group information.
group.wt <- vector()

# use find.group to identify the groups.
# use check.group to check the groups manually.
# use merge.group to merge the groups after they have been checked.
i <- rc$chr=="1A"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp[[1]] <- temp[[1]][ temp[[1]] >= 58 ]
temp <- merge.group(gid=list(c(2), c(1)), group=temp, n=sum(i))
group.wt <- c(group.wt, temp)

i <- rc$chr=="1B"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.wt, na.rm=TRUE)
group.wt <- c(group.wt, temp)

i <- rc$chr=="1D"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.wt, na.rm=TRUE)
group.wt <- c(group.wt, temp)

i <- rc$chr=="2A"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp[[1]] <- temp[[1]][ temp[[1]] <= 162 ]
temp <- merge.group(gid=list(c(2), c(4), c(1), c(3)), group=temp, n=sum(i))
temp <- temp + max(group.wt, na.rm=TRUE)
group.wt <- c(group.wt, temp)

i <- rc$chr=="2B"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(1:5), c(6)), group=temp, n=sum(i))
temp <- temp + max(group.wt, na.rm=TRUE)
group.wt <- c(group.wt, temp)

i <- rc$chr=="2D"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.wt, na.rm=TRUE)
group.wt <- c(group.wt, temp)

i <- rc$chr=="3A"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(2), c(1,3)), group=temp, n=sum(i))
temp <- temp + max(group.wt, na.rm=TRUE)
group.wt <- c(group.wt, temp)

i <- rc$chr=="3B"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.wt <- c(group.wt, temp)

i <- rc$chr=="3D"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.wt <- c(group.wt, temp)

i <- rc$chr=="4A"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp[[1]] <- temp[[1]][ temp[[1]] >= 42 ]
temp <- merge.group(gid=list(c(1,2)), group=temp, n=sum(i))
temp <- temp + max(group.wt, na.rm=TRUE)
group.wt <- c(group.wt, temp)

i <- rc$chr=="4B"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.wt, na.rm=TRUE)
group.wt <- c(group.wt, temp)

i <- rc$chr=="4D"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.wt <- c(group.wt, temp)

i <- rc$chr=="5A"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.wt, na.rm=TRUE)
group.wt <- c(group.wt, temp)

i <- rc$chr=="5B"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.wt, na.rm=TRUE)
group.wt <- c(group.wt, temp)

i <- rc$chr=="5D"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.wt <- c(group.wt, temp)

i <- rc$chr=="6A"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(1,2)), group=temp, n=sum(i))
temp <- temp + max(group.wt, na.rm=TRUE)
group.wt <- c(group.wt, temp)

i <- rc$chr=="6B"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.wt <- c(group.wt, temp)

i <- rc$chr=="6D"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.wt <- c(group.wt, temp)

i <- rc$chr=="7A"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
check.group(group=temp)
temp <- merge.group(gid=list(c(1)), group=temp, n=sum(i))
temp <- temp + max(group.wt, na.rm=TRUE)
group.wt <- c(group.wt, temp)

i <- rc$chr=="7B"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.wt <- c(group.wt, temp)

i <- rc$chr=="7D"
temp <- find.group(rc=rc[i, ], tsc=tsc, geno=geno[, i], min.r2=0.2, min.marker=5)
temp <- rep(NA, sum(i))
group.wt <- c(group.wt, temp)

group.wt[is.na(group.wt)] <- 0
rc$group <- group.wt


# create the manhattan plot for rc and highlight the groups.
plot.wt <- list()
for(i in 1:21){
  chr <- unique(rc$chr)[i]
  n.group <- length(unique(rc$group[rc$chr==chr]))
  color.opt <- hcl(h=seq(15, 375, length=n.group), l=65, c=100)
  color.opt <- color.opt[-n.group]
  
  plot.wt[[i]] <- ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
    annotate("segment", x=-Inf, xend=Inf, y=tsc, yend=tsc, color="#AAAAAA", linetype=2) +
    geom_point(data=rc[rc$chr==chr, ], aes(x=pos/1e6, y=score, color=as.factor(group)), size=1, alpha=0.8, stroke=0) +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    theme(axis.title=element_text(size=8), axis.text=element_text(size=8)) +
    theme(legend.title=element_text(size=8), legend.text=element_text(size=8)) +
    theme(legend.key.size=unit(0.75, "lines"), plot.subtitle=element_text(size=8)) +
    theme(legend.justification="top") +
    guides(color=guide_legend(override.aes=list(size=2))) +
    scale_color_manual(values=c("#555555", color.opt), name="group") +
    scale_y_continuous(limits=c(0, max(rc$score))) +
    scale_x_continuous(limits=c(0, max(rc$pos)/1e6)) +
    xlab("Pos (Mb)") +
    ylab(expression("-log"[10]*"P")) +
    labs(subtitle=chr)
}

ggsave(plot=arrangeGrob(grobs=plot.wt, layout_matrix=matrix(1:21, nrow=7, ncol=3, byrow=TRUE)),
       filename="output/S6_rc_WT.png",
       width=7,
       height=9,
       units="in",
       dpi=600,
       scale=4/3)
ggsave(plot=arrangeGrob(grobs=plot.wt, layout_matrix=matrix(1:21, nrow=7, ncol=3, byrow=TRUE)),
       filename="output/S6_rc_WT.svg",
       width=7,
       height=9,
       units="in",
       scale=4/3)

# identify the chr, pos range and peak pos for each group.
group.wt <- vector()
for(i in 1:19){
  temp <- which(rc$group==i)
  group.wt <- rbind(group.wt, 
                    c(i,
                      length(temp),
                      rc$chr[temp[1]],
                      rc$pos[temp][which.max(rc$score[temp])],
                      range(rc$pos[temp])))
}
group.wt <- data.frame(group.wt)
colnames(group.wt) <- c("group", "count", "chr", "peak", "start", "end")

# export the rc results.
write.csv(rc, "output/rc_WT.csv", quote=FALSE, row.names=FALSE)
write.csv(group.wt, "output/rc_WT_group.csv", quote=FALSE, row.names=FALSE)


save(geno, marker, pheno, rc, group.wt, tsc, file="analysis_WT.RData")
