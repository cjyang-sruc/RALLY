library(ggplot2)

# set working directory.
setwd()

# load the RALLY results for TG.
load("analysis_TG2.RData")

# load the GWAS peaks from Ladejobi et al (2019).
peak <- read.csv("input/GWAS_peaks.csv", as.is=TRUE)

# change group.gbs columns from characters to numerics.
for(i in c(1,2,4,5,6)) group.gbs[,i] <- as.numeric(group.gbs[,i])
group.gbs$score <- NA
for(i in 1:22) group.gbs$score[i] <- rc.gbs$score[rc.gbs$chr==group.gbs$chr[i] & rc.gbs$pos==group.gbs$peak[i]]

# remove the last row and h2 columns in group.gbs
group.gbs <- group.gbs[1:22, c(1:6,19)]

# change chr from characters to factors.
rc.gbs$chr <- as.factor(rc.gbs$chr)
group.gbs$chr <- as.factor(group.gbs$chr)
group.gbs$chr <- factor(group.gbs$chr, levels=levels(rc.gbs$chr))
peak$chr <- as.factor(peak$chr)
peak$chr <- factor(peak$chr, levels=levels(rc.gbs$chr))

### manhattan plot.

plot1 <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  annotate("segment", x=-Inf, xend=Inf, y=tsc.gbs, yend=tsc.gbs, color="#AAAAAA", linetype=2) +
  geom_point(data=rc.gbs, aes(x=pos/1e6, y=score, color=chr), size=0.5, show.legend=FALSE) +
  geom_point(data=group.gbs, aes(x=peak/1e6, y=score), color="#FF7777", size=1.2) +
  geom_segment(data=group.gbs, aes(x=start/1e6, xend=end/1e6, y=score, yend=score), color="#FF7777", size=0.8) +
  geom_point(data=peak, aes(x=pos/1e6, y=score), color="#0077FF", size=1.2) +
  facet_grid(cols=vars(chr), scales="free_x", space="free_x") +
  theme(panel.spacing.x=unit(0.2, "lines")) +
  theme(strip.background=element_blank()) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4, size=6)) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.ticks=element_line(color="#DDDDDD")) +
  scale_y_continuous(expand=c(0.01, 0), breaks=seq(0, 12, 2)) +
  scale_x_continuous(breaks=seq(0, 800, 200), expand=c(0.1, 0)) +
  scale_color_manual(values=c("#CCCCCC", "#999999")[c(1:21)%%2 + 1]) +
  ylab(expression("-log"[10]*"p")) +
  xlab("position (Mb)")

ggsave(plot=plot1,
       filename="output/3_manhattan_TG_with_GWAS_anno.png",
       height=2.75,
       width=7,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(plot=plot1,
       filename="output/3_manhattan_TG_with_GWAS_anno.svg",
       height=2.75,
       width=7,
       units="in",
       scale=4/3)


### bar plot of the unique/overlapping hits from RALLY and GWAS.
dat.plot <- data.frame(method=c(rep("AT", 4), rep("AW", 4), rep("GT", 4)),
                       overlap=c("unique", "AT-AW", "AT-GT", "AT-AW-GT",
                                 "unique", "AT-AW", "AW-GT", "AT-AW-GT",
                                 "unique", "AT-GT", "AW-GT", "AT-AW-GT"),
                       count=c(7,3,5,7, 8,4,1,6, 31,8,2,9))
dat.plot$overlap <- as.factor(dat.plot$overlap)
dat.plot$overlap <- factor(dat.plot$overlap, levels=c("unique", "AT-AW", "AT-GT", "AW-GT", "AT-AW-GT"))

plot1 <- ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_col(data=dat.plot, aes(x=method, y=count, fill=overlap)) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  scale_fill_manual(values=c("#CCCCCC", "#66CCCC", "#6666CC", "#CCCC66", "#CC66CC"))

ggsave(plot=plot1,
       filename="output/4_barplot_RALLY_GWAS_overlap.png",
       height=3,
       width=3.5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(plot=plot1,
       filename="output/4_barplot_RALLY_GWAS_overlap.svg",
       height=3,
       width=3.5,
       units="in",
       scale=4/3)

