library(ggplot2)
library(ggforce)
library(reshape2)
library(gridExtra)
library(rrBLUP)
library(glmnet)

# set working directory.
setwd()

# load the RALLY results for TG.
load("analysis_TG2.RData")

### PLOT 1: plot the variety count per year, colored by country.
ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_bar(data=pheno.gbs, aes(x=Year, fill=Country)) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  scale_x_continuous(expand=c(0,0), breaks=seq(1950, 2000, 10)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values=c("#882255", "#44AA99", "#DDCC77"))

ggsave(filename="output/S2A_variety_count_per_year.png",
       height=2.5,
       width=5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/S2A_variety_count_per_year.svg",
       height=2.5,
       width=5,
       units="in",
       scale=4/3)

### PLOT 2: plot the PCA, colored by country.
PC <- eigen(A.mat(geno.gbs-1))
head(PC$values/sum(PC$values)*100) # PC1 = 7.42%, PC2 = 5.11%
PC <- data.frame(Variety=pheno.gbs$Variety, Year=pheno.gbs$Year, Country=pheno.gbs$Country, PC1=PC$vectors[,1], PC2=PC$vectors[,2])

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=PC, aes(x=PC1, y=PC2, color=Country), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77")) +
  xlab("PC1 (7.42 %)") +
  ylab("PC2 (5.11 %)") +
  coord_fixed()

ggsave(filename="output/S2B_PCA.png",
       height=2.5,
       width=2.5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/S2B_PCA.svg",
       height=2.5,
       width=2.5,
       units="in",
       scale=4/3)


# load the RALLY results for WAGTAIL.
load("analysis_WT.RData")

### PLOT 1: plot the variety count per year, colored by country.
pheno$Country <- as.factor(pheno$Country)
pheno$Country <- factor(pheno$Country, levels=c("DE", "FR", "UK", "AU", "BE", "CA", "CH", "DK", "NL", "SE", "US"))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_bar(data=pheno, aes(x=Year, fill=Country)) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  scale_x_continuous(expand=c(0,0), breaks=seq(1920, 2000, 20)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values=c("#882255", "#44AA99", "#DDCC77", "#1F78B4", "#E31A1C", "#FF7F00", "#CAB2D6", "#B2DF8A", "#A6CEE3", "#FB9A99", "#FDBF6F"))

ggsave(filename="output/S2C_variety_count_per_year_WT.png",
       height=2.5,
       width=5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/S2C_variety_count_per_year_WT.svg",
       height=2.5,
       width=5,
       units="in",
       scale=4/3)

### PLOT 2: plot the PCA, colored by country.
PC <- eigen(A.mat(geno-1))
head(PC$values/sum(PC$values)*100) # PC1 = 7.95%, PC2 = 4.93%
PC <- data.frame(Variety=pheno$Name, Year=pheno$Year, Country=pheno$Country, PC1=PC$vectors[,1], PC2=PC$vectors[,2])

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=PC, aes(x=PC1, y=PC2, color=Country), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  scale_color_manual(values=c("#882255", "#44AA99", "#DDCC77", "#1F78B4", "#E31A1C", "#FF7F00", "#CAB2D6", "#B2DF8A", "#A6CEE3", "#FB9A99", "#FDBF6F")) +
  xlab("PC1 (7.95 %)") +
  ylab("PC2 (4.93 %)") +
  coord_fixed()

ggsave(filename="output/S2D_PCA_WT.png",
       height=2.5,
       width=2.5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/S2D_PCA_WT.svg",
       height=2.5,
       width=2.5,
       units="in",
       scale=4/3)






# convert geno.gbs from 0/1/2 to -1/0/1 format.
geno.gbs2 <- geno.gbs - 1

# predict using rrBLUP.
# create 5-fold cv.
set.seed(54715)
cv <- c(sample(1:nrow(pheno.gbs), nrow(pheno.gbs), replace=FALSE), 0, 0)
cv <- lapply(1:5, FUN=function(x) sort(cv[(67*(x-1)+1):(67*x)]))
cv[[5]] <- cv[[5]][-c(1,2)]

# rrBLUP.
rr <- list(mean=vector(),
           meff=vector(),
           pred=matrix(NA, nrow=nrow(pheno.gbs), ncol=12))

for(i in 1:12){
  
  mean.0 <- vector()
  meff.0 <- vector()
  
  for(j in 1:5){
    temp <- mixed.solve(y=pheno.gbs[-cv[[j]], i+3], Z=geno.gbs2[-cv[[j]], ])
    mean.0 <- c(mean.0, c(temp$beta))
    meff.0 <- cbind(meff.0, c(temp$u))
    rr$pred[cv[[j]], i] <- c(geno.gbs2[cv[[j]], ]%*%meff.0[,j] + mean.0[j])
  }
  
  rr$mean <- c(rr$mean, mean(mean.0))
  rr$meff <- cbind(rr$meff, rowSums(meff.0)/rowSums(!(meff.0==0)))
  message(i)
  
}

# glmnet LASSO.
# glmnet (alpha = 1).
lasso <- list(mean=vector(),
              meff=vector(),
              pred=matrix(NA, nrow=nrow(pheno.gbs), ncol=12))


mean.0 <- vector()
meff.0 <- replicate(5, list(matrix(0, nrow=ncol(geno.gbs2), ncol=12)))

for(j in 1:5){
  temp <- cv.glmnet(y=as.matrix(pheno.gbs[-cv[[j]], 4:15]), x=geno.gbs2[-cv[[j]], ], alpha=1, family="mgaussian")
  temp <- coef.glmnet(temp, s="lambda.min")
  mean.0 <- cbind(mean.0, sapply(1:12, FUN=function(y) temp[[y]]@x[1]))
  meff.0[[j]][temp[[1]]@i[-1], ] <- sapply(1:12, FUN=function(y) temp[[y]]@x[-1])
  lasso$pred[cv[[j]], ] <- geno.gbs2[cv[[j]], ]%*%meff.0[[j]] + matrix(rep(mean.0[,j], length(cv[[j]])), ncol=12, byrow=TRUE)
  
  message(j)
}

lasso$mean <- rowSums(mean.0)/5
temp <- !(meff.0[[1]]==0)
for(i in 2:5) temp <- temp + !(meff.0[[i]]==0)
temp2 <- meff.0[[1]]
for(i in 2:5) temp2 <- temp2 + meff.0[[i]]
lasso$meff <- temp2/temp


### Figure S6: plot the observed/predicted trait over year (rrBLUP).
dat.plot <- data.frame(year=rep(pheno.gbs$Year, 2),
                       type=c(rep("observed", nrow(pheno.gbs)), rep("predicted", nrow(pheno.gbs))),
                       rbind(as.matrix(pheno.gbs[,4:15]), rr$pred))
dat.plot <- melt(dat.plot, id.vars=c("year", "type"))
colnames(dat.plot)[3] <- "trait"
dat.plot$trait <- factor(dat.plot$trait, labels=paste(levels(dat.plot$trait),
                                                      " (r = ",
                                                      formatC(x=sapply(1:12, FUN=function(x) cor(rr$pred[,x], pheno.gbs[,x+3])), digits=3, format="f"),
                                                      ")",
                                                      sep=""))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_jitter(data=dat.plot, aes(x=year, y=value, color=type), stroke=0, alpha=0.7, width=0.2) +
  geom_smooth(data=dat.plot[dat.plot$type=="observed", ], aes(x=year, y=value), method="lm", se=FALSE, color="#999999") +
  facet_wrap(vars(trait), ncol=3, scales="free_y") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(strip.background=element_blank(), strip.text=element_text(hjust=0)) +
  guides(color=guide_legend(override.aes=list(size=3))) +
  scale_x_continuous(breaks=seq(1950, 2000, 10)) +
  scale_color_manual(values=c("#44AA99", "#882255")) +
  ylab("trait value")

ggsave(filename="output/S7_obs_pred_year_rr.png",
       height=6,
       width=7,
       units="in",
       dpi=600,
       scale=4/3)



### Figure S7: plot the observed/predicted trait over year (LASSO).
dat.plot <- data.frame(year=rep(pheno.gbs$Year, 2),
                       type=c(rep("observed", nrow(pheno.gbs)), rep("predicted", nrow(pheno.gbs))),
                       rbind(as.matrix(pheno.gbs[,4:15]), lasso$pred))
dat.plot <- melt(dat.plot, id.vars=c("year", "type"))
colnames(dat.plot)[3] <- "trait"
dat.plot$trait <- factor(dat.plot$trait, labels=paste(levels(dat.plot$trait),
                                                      " (r = ",
                                                      formatC(x=sapply(1:12, FUN=function(x) cor(lasso$pred[,x], pheno.gbs[,x+3])), digits=3, format="f"),
                                                      ")",
                                                      sep=""))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_jitter(data=dat.plot, aes(x=year, y=value, color=type), stroke=0, alpha=0.7, width=0.2) +
  geom_smooth(data=dat.plot[dat.plot$type=="observed", ], aes(x=year, y=value), method="lm", se=FALSE, color="#999999") +
  facet_wrap(vars(trait), ncol=3, scales="free_y") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(strip.background=element_blank(), strip.text=element_text(hjust=0)) +
  guides(color=guide_legend(override.aes=list(size=3))) +
  scale_x_continuous(breaks=seq(1950, 2000, 10)) +
  scale_color_manual(values=c("#44AA99", "#882255")) +
  ylab("trait value")

ggsave(filename="output/S8_obs_pred_year_lasso.png",
       height=6,
       width=7,
       units="in",
       dpi=600,
       scale=4/3)


# identify the directions of allele effects and regression over year coefficients.
rr$dir <- sapply(1:12, FUN=function(x) sign(rc.gbs$Z.adj)*sign(rr$meff[,x]))
lasso$dir <- sapply(1:12, FUN=function(x) sign(rc.gbs$Z.adj)*sign(lasso$meff[,x]))

### PLOT 6: Plot the increasing/decreasing effects w.r.t. allele reg over year (rrBLUP).
dat.plot <- data.frame(trait=rep(colnames(pheno.gbs)[4:15], 6),
                       sig=c(rep("A", 12*2), rep("B", 12*2), rep("C", 12*2)),
                       direction=rep(c(rep("positive", 12), rep("negative", 12)), 3),
                       count=c(colSums(rr$dir[rc.gbs$P.adj <= 0.05/nrow(rc.gbs), ] > 0),
                               colSums(rr$dir[rc.gbs$P.adj <= 0.05/nrow(rc.gbs), ] < 0),
                               colSums(rr$dir[rc.gbs$P.adj > 0.05/nrow(rc.gbs) & rc.gbs$P.adj <= 0.05, ] > 0),
                               colSums(rr$dir[rc.gbs$P.adj > 0.05/nrow(rc.gbs) & rc.gbs$P.adj <= 0.05, ] < 0),
                               colSums(rr$dir[rc.gbs$P.adj > 0.05, ] > 0),
                               colSums(rr$dir[rc.gbs$P.adj > 0.05, ] < 0)))
dat.plot$trait <- as.factor(dat.plot$trait)
dat.plot$trait <- factor(dat.plot$trait, levels=colnames(pheno.gbs)[4:15])
dat.plot$trait <- factor(dat.plot$trait, labels=c("FT", "LODG", "YLD", "HT", "PROT", "WK", "AWNS", "SPWT", "TGW", "EM2", "TILL", "MAT"))
dat.plot$direction <- as.factor(dat.plot$direction)
dat.plot$direction <- factor(dat.plot$direction, levels=c("positive", "negative"))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_col(data=dat.plot, aes(x=trait, y=count, fill=direction)) +
  facet_wrap(vars(sig), nrow=3, scales="free_y") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(strip.background=element_blank(), strip.text=element_text(hjust=0)) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=8)) +
  theme(axis.text.y=element_text(size=8)) +
  scale_y_continuous(n.breaks=3) +
  scale_fill_manual(values=c("#44AA99", "#882255"))

ggsave(filename="output/6_allele_effect_direction_rr.png",
       height=3.5,
       width=3.5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/6_allele_effect_direction_rr.svg",
       height=3.5,
       width=3.5,
       units="in",
       scale=4/3)


### PLOT S8: Plot the increasing/decreasing effects w.r.t. allele reg over year (LASSO).
dat.plot <- data.frame(trait=rep(colnames(pheno.gbs)[4:15], 6),
                       sig=c(rep("A", 12*2), rep("B", 12*2), rep("C", 12*2)),
                       direction=rep(c(rep("positive", 12), rep("negative", 12)), 3),
                       count=c(colSums(lasso$dir[rc.gbs$P.adj <= 0.05/nrow(rc.gbs), ] > 0, na.rm=TRUE),
                               colSums(lasso$dir[rc.gbs$P.adj <= 0.05/nrow(rc.gbs), ] < 0, na.rm=TRUE),
                               colSums(lasso$dir[rc.gbs$P.adj > 0.05/nrow(rc.gbs) & rc.gbs$P.adj <= 0.05, ] > 0, na.rm=TRUE),
                               colSums(lasso$dir[rc.gbs$P.adj > 0.05/nrow(rc.gbs) & rc.gbs$P.adj <= 0.05, ] < 0, na.rm=TRUE),
                               colSums(lasso$dir[rc.gbs$P.adj > 0.05, ] > 0, na.rm=TRUE),
                               colSums(lasso$dir[rc.gbs$P.adj > 0.05, ] < 0, na.rm=TRUE)))
dat.plot$trait <- as.factor(dat.plot$trait)
dat.plot$trait <- factor(dat.plot$trait, levels=colnames(pheno.gbs)[4:15])
dat.plot$trait <- factor(dat.plot$trait, labels=c("FT", "LODG", "YLD", "HT", "PROT", "WK", "AWNS", "SPWT", "TGW", "EM2", "TILL", "MAT"))
dat.plot$direction <- as.factor(dat.plot$direction)
dat.plot$direction <- factor(dat.plot$direction, levels=c("positive", "negative"))

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_col(data=dat.plot, aes(x=trait, y=count, fill=direction)) +
  facet_wrap(vars(sig), nrow=3, scales="free_y") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(strip.background=element_blank(), strip.text=element_text(hjust=0)) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=8)) +
  theme(axis.text.y=element_text(size=8)) +
  scale_y_continuous(n.breaks=3) +
  scale_fill_manual(values=c("#44AA99", "#882255"))

ggsave(filename="output/S9_allele_effect_direction_lasso.png",
       height=3.5,
       width=3.5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/S9_allele_effect_direction_lasso.svg",
       height=3.5,
       width=3.5,
       units="in",
       scale=4/3)


### Plot allele effects for pairs of traits.
dat.plot <- vector()
dat.sig <- vector()
for(i in 1:11){
  for(j in (i+1):12){
    temp <- data.frame(trait1=sign(rc.gbs$Z.adj)*sign(rr$meff[,i]),
                       trait2=sign(rc.gbs$Z.adj)*sign(rr$meff[,j]))
    colnames(temp) <- colnames(pheno.gbs)[c(3+i, 3+j)]
    dat.plot <- rbind(dat.plot,
                      data.frame(trait1=colnames(pheno.gbs)[3+i],
                                 trait2=colnames(pheno.gbs)[3+j],
                                 x=c(4*i-1, 4*i-3, 4*i-1, 4*i-3),
                                 y=c(4*j-1, 4*j-1, 4*j-3, 4*j-3),
                                 count=c(table(temp)),
                                 count2=sqrt(c(table(temp))/nrow(temp)),
                                 type=c("same", "opposite", "opposite", "same"),
                                 model="RR"))
    dat.sig <- rbind(dat.sig,
                     data.frame(trait1=colnames(pheno.gbs)[3+i],
                                trait2=colnames(pheno.gbs)[3+j],
                                x1=4*i-4,
                                x2=4*i,
                                y1=4*j-4,
                                y2=4*j,
                                model="RR",
                                P=chisq.test(as.matrix(table(temp)))$p.value))
    
    temp <- data.frame(trait1=sign(rc.gbs$Z.adj)*sign(lasso$meff[,i]),
                       trait2=sign(rc.gbs$Z.adj)*sign(lasso$meff[,j]))
    colnames(temp) <- colnames(pheno.gbs)[c(3+i, 3+j)]
    dat.plot <- rbind(dat.plot,
                      data.frame(trait1=colnames(pheno.gbs)[3+i],
                                 trait2=colnames(pheno.gbs)[3+j],
                                 x=c(4*j-1, 4*j-1, 4*j-3, 4*j-3),
                                 y=c(4*i-1, 4*i-3, 4*i-1, 4*i-3),
                                 count=c(table(temp)),
                                 count2=sqrt(c(table(temp))/sum(rowSums(!is.na(temp))==2)),
                                 type=c("same", "opposite", "opposite", "same"),
                                 model="LASSO"))
    dat.sig <- rbind(dat.sig,
                     data.frame(trait1=colnames(pheno.gbs)[3+i],
                                trait2=colnames(pheno.gbs)[3+j],
                                x1=4*j-4,
                                x2=4*j,
                                y1=4*i-4,
                                y2=4*i,
                                model="LASSO",
                                P=chisq.test(as.matrix(table(temp)))$p.value))
  }
}

dat.sig2 <- dat.sig[dat.sig$model=="LASSO" & dat.sig$P <= 0.05/66, ]

ggplot() +
  geom_rect(data=dat.sig2, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="#FAFAD2", color=NA) +
  annotate("segment", x=-Inf, xend=Inf, y=seq(0, 48, 4), yend=seq(0, 48, 4), color="#DDDDDD") +
  annotate("segment", y=-Inf, yend=Inf, x=seq(0, 48, 4), xend=seq(0, 48, 4), color="#DDDDDD") +
  geom_circle(data=dat.plot, aes(x0=x, y0=y, r=count2, fill=type), color=NA) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank()) +
  scale_fill_manual(values=c("#44AA99", "#BBBBBB")) +
  scale_y_reverse(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  coord_fixed()

ggsave(filename="output/7_allele_effect_direction_trait_pairs.png",
       height=6,
       width=6.5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/7_allele_effect_direction_trait_pairs.svg",
       height=6,
       width=6.5,
       units="in",
       scale=4/3)

save(rr, lasso, file="allele_effect.RData")
write.csv(dat.plot, "output/source_data1_figure7.csv", quote=FALSE, row.names=FALSE)
write.csv(cbind(dat.sig, matrix(dat.plot$count, ncol=4, byrow=TRUE)), "output/source_data2_figure7.csv", quote=FALSE, row.names=FALSE)
write.csv(rc.mloci, "output/rc_mloci.csv", quote=FALSE, row.names=FALSE)

# plot the allele regressions at major loci.
plot.mloci <- list()
mloci.names <- c("vrn_A1", "vrn_B1", "vrn_D1", "vrn_B3", "1BL_1RS", "Ppd_B1a", "Ppd_A1a", "Ppd_D1a", "Rht_B1b", "Rht_D1b", "Rht8_165", "Rht8_174", "Rht8_192")

for(i in c(1:5,7,6,8:13)){
  temp <- data.frame(allele=mloci[,i+1],
                     year=pheno.gbs$Year,
                     country=pheno.gbs$Country)
  temp.glm <- tryCatch(glm(allele ~ year + country, data=temp, family=binomial(link="logit")),
                       error=function(e){})
  
  temp2 <- rbind(data.frame(fit=predict(temp.glm, newdata=data.frame(year=1948:2007, country="DE"), type="response"),
                            year=c(1948:2007),
                            country="DE"),
                 data.frame(fit=predict(temp.glm, newdata=data.frame(year=1948:2007, country="FR"), type="response"),
                            year=c(1948:2007),
                            country="FR"),
                 data.frame(fit=predict(temp.glm, newdata=data.frame(year=1948:2007, country="UK"), type="response"),
                            year=c(1948:2007),
                            country="UK"))
  
  temp3 <- data.frame(freq=tapply(temp$allele, temp$year, mean, na.rm=TRUE),
                      year=as.numeric(names(tapply(temp$allele, temp$year, mean, na.rm=TRUE))))
  
  plot.mloci[[i]] <- ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
    geom_col(data=temp3, aes(x=year, y=freq), color="#CCCCCC", fill="#EEEEEE", na.rm=TRUE) +
    geom_jitter(data=temp, aes(x=year, y=allele), width=0.4, height=0.02, size=0.5, na.rm=TRUE, color="#555555") +
    geom_line(data=temp2, aes(x=year, y=fit, color=country)) +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) +
    theme(legend.position="none") +
    scale_x_continuous(breaks=seq(1950,2000,10)) +
    scale_color_manual(values=c("#FF5555", "#55FF55", "#5555FF")) +
    ylab("allele") +
    labs(subtitle=mloci.names[i])
  
}

plot.mloci[[14]] <- ggplot() + geom_blank() + theme_void()
plot.mloci[[15]] <- ggplot() + geom_blank() + theme_void()

ggsave(plot=arrangeGrob(grobs=plot.mloci, layout_matrix=matrix(1:15, ncol=3, byrow=TRUE)),
       filename="output/mloci.png",
       height=7,
       width=7,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(plot=arrangeGrob(grobs=plot.mloci, layout_matrix=matrix(1:15, ncol=3, byrow=TRUE)),
       filename="output/mloci.svg",
       height=7,
       width=7,
       units="in",
       scale=4/3)






