library(ggplot2)

### set working directory.
setwd("")

### calculate the probability of haplotype combinations
### given one generation of selfing of each combination.
### x0 is the index of haplotype combination
### 1 = QM:QM, 2 = QM:Qm, 3 = QM:qM, 4 = QM:qm, 5 = Qm:Qm
### 6 = Qm:qM, 7 = Qm:qm, 8 = qM:qM, 9 = qM:qm, 10 = qm:qm
### p0 is the probability.
### r is the recombination fraction between QTL (Q) and marker (M).
hcalc <- function(x0, p0, r){
  if(x0==1 | x0==5 | x0==8 | x0==10){
    x <- x0
    p <- p0*1
  } else if(x0==2){
    x <- c(1,5,2)
    p <- p0*c(1/4, 1/4, 1/2)
  } else if(x0==3){
    x <- c(1,8,3)
    p <- p0*c(1/4, 1/4, 1/2)
  } else if(x0==7){
    x <- c(5,10,7)
    p <- p0*c(1/4, 1/4, 1/2)
  } else if(x0==9){
    x <- c(8,10,9)
    p <- p0*c(1/4, 1/4, 1/2)
  } else if(x0==4){
    x <- c(1,10,4,2,3,7,9,5,8,6)
    p <- p0*c((1-r)^2/4,(1-r)^2/4,(1-r)^2/2,r*(1-r)/2,r*(1-r)/2,r*(1-r)/2,r*(1-r)/2,r^2/4,r^2/4,r^2/2)
  } else if(x0==6){
    x <- c(5,8,6,2,3,7,9,1,10,4)
    p <- p0*c((1-r)^2/4,(1-r)^2/4,(1-r)^2/2,r*(1-r)/2,r*(1-r)/2,r*(1-r)/2,r*(1-r)/2,r^2/4,r^2/4,r^2/2)
  }
  return(list(x=x, p=p))
}

### calculate the haplotype frequency.
### a = c(QM, Qm, qM, qm)
hprob <- function(x, p){
  a <- rep(0, 4)
  a[1] <- sum(p[x==1]) + sum(p[x==2]/2) + sum(p[x==3]/2) + sum(p[x==4]/2)
  a[2] <- sum(p[x==5]) + sum(p[x==2]/2) + sum(p[x==6]/2) + sum(p[x==7]/2)
  a[3] <- sum(p[x==8]) + sum(p[x==3]/2) + sum(p[x==6]/2) + sum(p[x==9]/2)
  a[4] <- sum(p[x==10]) + sum(p[x==4]/2) + sum(p[x==7]/2) + sum(p[x==9]/2)
  return(a/sum(a)) #dividing by the sum to fix the rounding issue so a always sum to 1.
}


### calculate the marker frequency given initial haplotype frequency (a0),
### QTL frequency across all generations (q) and recombination fraction (r).
hsel <- function(a0, q, r){
  
  n.gen <- length(q) - 1
  m <- a0[1] + a0[3]
  a <- a0
  
  for(h in 1:n.gen){
    
    if(h > 1) a0 <- c(a.star, b.star, c.star, d.star)
    
    x0 <- 1:10
    p0 <- c(a0[1]^2, 2*a0[1]*a0[2], 2*a0[1]*a0[3], 2*a0[1]*a0[4], a0[2]^2,
            2*a0[2]*a0[3], 2*a0[2]*a0[4], a0[3]^2, 2*a0[3]*a0[4], a0[4]^2)
    
    x <- replicate(6, vector())
    p <- replicate(6, vector())
    x[[1]] <- x0
    p[[1]] <- p0
    for(i in 2:6){
      for(j in 1:length(x[[i-1]])){
        temp <- hcalc(x0=x[[i-1]][j], p0=p[[i-1]][j], r=r)
        x[[i]] <- c(x[[i]], temp$x)
        p[[i]] <- c(p[[i]], temp$p)
      }
    }
    
    # m* = a' x q* / q + c' x t* / t
    a.prime <- hprob(x=x[[6]], p=p[[6]])[1]
    b.prime <- hprob(x=x[[6]], p=p[[6]])[2]
    c.prime <- hprob(x=x[[6]], p=p[[6]])[3]
    d.prime <- hprob(x=x[[6]], p=p[[6]])[4]
    
    q.ratio <- q[h+1]/q[h]
    t.ratio <- (1-q[h+1])/(1-q[h])
    
    a.star <- a.prime*q.ratio
    b.star <- b.prime*q.ratio
    c.star <- c.prime*t.ratio
    d.star <- d.prime*t.ratio
    
    m <- c(m, a.star + c.star)
    a <- rbind(a, c(a.star, b.star, c.star, d.star))
  }
  
  return(list(m=m, a=a))
}

### derive the QTL frequency under a logistic model.
qlogit <- function(q0, qF, start, end){
  b0 <- log(q0/(1-q0))
  b1 <- 1/(end-start)*log(qF/(1-qF)*(1-q0)/q0)
  return(1/(1 + exp(-b0-b1*c(start:end))))
}

### derive the QTL frequency under a linear model.
qlinear <- function(q0, qF, start, end){
  b0 <- q0
  b1 <- (qF-q0)/(end-start)
  return(b0 + b1*c(start:end))
}

### create a large matrix for all possible initial QTL-marker haplotype frequency.
init <- vector()
for(q0 in seq(1/32, 16/32, 1/32)){
  for(aa in seq(0, q0, 1/32)){
    for(bb in seq(0, q0, 1/32)){
      for(cc in seq(0, 1-q0, 1/32)){
        for(dd in seq(0, 1-q0, 1/32)){
          temp <- c(aa, bb, cc, dd)
          if(sum(temp[1:2])==q0 & sum(temp)==1) init <- rbind(init, c(q0, temp))
        }
      }
    }
  }
}

### QTL increases using logistic model.
df <- vector()
Z <- list()
for(r in seq(0.01, 0.10, 0.01)){
  z <- matrix(NA, nrow=100, ncol=nrow(init))
  for(val in 1:nrow(init)){
    q <- qlogit(q0=init[val, 1], qF=31/32, start=0, end=50)
    m <- hsel(a0=init[val, 2:5], q=q, r=r)$m[-1]
    temp <- vector()
    for(i in 1:100){
      allele <- c(sapply(1:length(m), FUN=function(j) rbinom(n=8, size=1, prob=m[j])))
      year <- sort(rep(1:length(m), 8))
      lm1 <- glm(allele ~ year, family=binomial)
      temp <- c(temp, summary(lm1)$coefficients[2,4] < 0.05/19000)
      z[i, val] <- summary(lm1)$coefficients[2,3]
    }
    df <- rbind(df, c(r, init[val, ], count=sum(temp, na.rm=TRUE)))
  }
  Z <- c(Z, list(z))
  message(r)
}
colnames(df) <- c("r", "q0", "aa", "bb", "cc", "dd", "count")
df <- data.frame(df)

save(df, Z, file="output/RALLY_precision_logistic.RData")

### plot the proportion of marker-qtl haplotypes that are significant in X out of 100 simulations.
### logistic QTL
#####
dat.plot <- vector()
for(i in unique(df$r)){
  for(j in unique(df$q0)){
    dat.plot <- rbind(dat.plot,
                      data.frame(r=i, q0=j, sig=sum(df$count[df$r==i & df$q0==j] > 0)/sum(df$r==i & df$q0==j), type="a"),
                      data.frame(r=i, q0=j, sig=sum(df$count[df$r==i & df$q0==j] > 9)/sum(df$r==i & df$q0==j), type="b"),
                      data.frame(r=i, q0=j, sig=sum(df$count[df$r==i & df$q0==j] > 49)/sum(df$r==i & df$q0==j), type="c"),
                      data.frame(r=i, q0=j, sig=sum(df$count[df$r==i & df$q0==j] == 100)/sum(df$r==i & df$q0==j), type="d"))
  }
}
dat.plot$q0 <- as.factor(dat.plot$q0)
dat.plot$q0 <- factor(dat.plot$q0, labels=paste(1:16, 32, sep="/"))
dat.plot$r <- 100*dat.plot$r

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=dat.plot, aes(x=r, y=sig, color=q0), size=1) +
  geom_line(data=dat.plot, aes(x=r, y=sig, color=q0)) +
  facet_wrap(vars(type), nrow=2, scales="free") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8), legend.key.size=unit(1, "lines")) +
  theme(strip.text=element_blank(), strip.background=element_blank()) +
  #guides(color=guide_legend(override.aes=list(size = 0.5), ncol=1)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) +
  scale_x_continuous(breaks=1:10) +
  ylab("proportion") +
  xlab("r (cM)")

ggsave(filename="output/S4_sig_q0_r_logistic.png",
       width=6.5,
       height=3,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/S4_sig_q0_r_logistic.svg",
       width=6.5,
       height=3,
       units="in",
       scale=4/3)

#####


### QTL increases using linear model.
df2 <- vector()
Z2 <- list()
for(r in seq(0.01, 0.10, 0.01)){
  z <- matrix(NA, nrow=100, ncol=nrow(init))
  for(val in 1:nrow(init)){
    q <- qlinear(q0=init[val, 1], qF=31/32, start=0, end=50)
    m <- hsel(a0=init[val, 2:5], q=q, r=r)$m[-1]
    temp <- vector()
    for(i in 1:100){
      allele <- c(sapply(1:length(m), FUN=function(j) rbinom(n=8, size=1, prob=m[j])))
      year <- sort(rep(1:length(m), 8))
      lm1 <- glm(allele ~ year, family=binomial)
      temp <- c(temp, summary(lm1)$coefficients[2,4] < 0.05/19000)
      z[i, val] <- summary(lm1)$coefficients[2,3]
    }
    df2 <- rbind(df2, c(r, init[val, ], count=sum(temp, na.rm=TRUE)))
  }
  Z2 <- c(Z2, list(z))
  message(r)
}
colnames(df2) <- c("r", "q0", "aa", "bb", "cc", "dd", "count")
df2 <- data.frame(df2)

save(df2, Z2, file="output/RALLY_precision_linear.RData")


### plot the proportion of marker-qtl haplotypes that are significant in X out of 100 simulations.
### linear QTL
#####
dat.plot <- vector()
for(i in unique(df2$r)){
  for(j in unique(df2$q0)){
    dat.plot <- rbind(dat.plot,
                      data.frame(r=i, q0=j, sig=sum(df2$count[df2$r==i & df2$q0==j] > 0)/sum(df2$r==i & df2$q0==j), type="a"),
                      data.frame(r=i, q0=j, sig=sum(df2$count[df2$r==i & df2$q0==j] > 9)/sum(df2$r==i & df2$q0==j), type="b"),
                      data.frame(r=i, q0=j, sig=sum(df2$count[df2$r==i & df2$q0==j] > 49)/sum(df2$r==i & df2$q0==j), type="c"),
                      data.frame(r=i, q0=j, sig=sum(df2$count[df2$r==i & df2$q0==j] == 100)/sum(df2$r==i & df2$q0==j), type="d"))
  }
}
dat.plot$q0 <- as.factor(dat.plot$q0)
dat.plot$q0 <- factor(dat.plot$q0, labels=paste(1:16, 32, sep="/"))
dat.plot$r <- 100*dat.plot$r

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_point(data=dat.plot, aes(x=r, y=sig, color=q0), size=1) +
  geom_line(data=dat.plot, aes(x=r, y=sig, color=q0)) +
  facet_wrap(vars(type), nrow=2, scales="free") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8), legend.key.size=unit(1, "lines")) +
  theme(strip.text=element_blank(), strip.background=element_blank()) +
  #guides(color=guide_legend(override.aes=list(size = 0.5), ncol=1)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) +
  scale_x_continuous(breaks=1:10) +
  ylab("proportion") +
  xlab("r (cM)")

ggsave(filename="output/S4_sig_q0_r_linear.png",
       width=6.5,
       height=3,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/S4_sig_q0_r_linear.svg",
       width=6.5,
       height=3,
       units="in",
       scale=4/3)

#####



### plot the change in marker frequency using two different initial
### QTL-marker haplotype frequencies as examples.
q <- qlogit(q0=1/32, qF=31/32, start=0, end=50)
m33 <- hsel(a0=init[33, 2:5], q=q, r=0.01)$m
m54 <- hsel(a0=init[54, 2:5], q=q, r=0.01)$m

set.seed(12304)
z33 <- colSums(sapply(1:length(m33), FUN=function(j) rbinom(n=8, size=1, prob=m33[j]))==1)/8
z54 <- colSums(sapply(1:length(m54), FUN=function(j) rbinom(n=8, size=1, prob=m54[j]))==1)/8

ggplot() +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  geom_line(aes(x=0:50, y=z33), color="#BB55BB", alpha=0.5) +
  geom_line(aes(x=0:50, y=z54), color="#55BBBB", alpha=0.5) +
  geom_point(aes(x=0:50, y=q), shape=15) +
  geom_point(aes(x=0:50, y=m33), color="#BB55BB", shape=16) +
  geom_point(aes(x=0:50, y=m54), color="#55BBBB", shape=17) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  xlab("generation") +
  ylab("allele frequency")

ggsave(filename="output/S1C_sim_qtl_marker_example.png",
       height=2.5,
       width=5,
       units="in",
       dpi=600,
       scale=4/3)

ggsave(filename="output/S1C_sim_qtl_marker_example.svg",
       height=2.5,
       width=5,
       units="in",
       scale=4/3)