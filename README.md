# `bmcp`
R package for the Normal data application of the Bayesian multipartition change point model introduced in [Pedroso et al (2021)](https://arxiv.org/abs/2107.11456).


# Install
To install `bmcp`, consider the code
```R
install.packages("devtools")
devtools::install_github("rcpedroso/bmcp")
```

If installation fails, please report the problem to ricardocunhap@gmail.com. Any other problems or comments can also be reported!


# Available functions in the `bmcp` package

### `bmcp`
Executes the partially colapsed Gibbs sampler scheme given in Algorithm 1 to sample from the posterior distribution described in section 3.2.

**Input:**

- burn: number of initial samples to be discarded as a burn-in period.
- ns: number of final samples to be generated.
- X: vector of observations.
- alpha1,beta1: parameter values for the Beta prior of the cohesion parameter p1.
- alpha2,beta2: parameter values for the Beta prior of the cohesion parameter p2.
- a,d: parameter values for the Inverse-Gamma prior of the cluster variance.
- mu0,s02: parameter values for the Normal prior of the cluster mean.

**Output:** an object of class `list` having the following objects:

- p1: vector of samples of the parameter p1.
- p2: vector of samples of the parameter p2.
- b1: vector of samples of the number of blocks in rho1.
- b2: vector of samples of the number of blocks in rho2.
- mu: matrix of samples of the mean vector.
- s2: matrix of samples of the variance vector.
- u: matrix of samples of the partition rho1.
- v: matrix of samples of the partition rho2.
- alpha1,beta1,alpha2,beta2,mu0,s02,a,d: parameter values.




### `u_to_index`
**Input:** any matrix M of 0/1 elements.

**Output:** if M has m rows and (n-1) columns, the function returns a matrix with m rows
and n columns where the elements indicate the positions of the 0's in each row
and column n indicates the number of 0's in each respective row. If M is the matrix where each line is a sample of some specific partition, the function returns the end points
and the number of end points of each sample.


# Code example
### Case study 1 of [Pedroso et al (2021)](https://arxiv.org/abs/2107.11456)

```R
### loading required packages
require(bmcp)
require(bcp)
require(MCMCvis)
require(ggplot2)
require(calibrate)
require(dplyr)


### data set
data(RealInt) # available in bcp package
Y <- as.vector(RealInt)
Yt <- as.vector(time(RealInt))
Ytq <- paste0(Yt,"/",((Yt-floor(Yt))+.25)/.25)
n = length(Y)


### Figure 10
plot(y=Y, x=Yt, type="l", xaxt="n", cex.lab=1.2,
     xlab="Quarters", ylab="US ex-post interest rate")
axis(1, at=Yt[seq(1,n,12)], labels=Ytq[seq(1,n,12)], cex.axis=1)


### MCMC (run time: approx. 30s using a computer with an Intel® Core i7-7500U/2.9GHz with 16Gb of RAM)
# settings
burn = 3e4
ns = 2e4
# run
set.seed(1000)
bmcp.est <- bmcp(burn=burn, ns=ns, X=Y,
                 alpha1=1, beta1=1, alpha2=1, beta2=1, a=0.1, d=2.1, mu0=0, s02=100)


### posterior distributions of partitions rho1 and rho2
# posterior samples of the end points
eps.U = u_to_index(bmcp.est$u[,-n])
eps.V = u_to_index(bmcp.est$v[,-n])
colnames(eps.U) <- c(paste0('t',1:(n-1)),'N1')
colnames(eps.V) <- c(paste0('t',1:(n-1)),'N2')
max1 = max(eps.U[,'N1'])
max2 = max(eps.V[,'N2'])
top = 5

# Table 4, rho1 posterior distribution (top most probable partitions for the mean)
top.U = as_tibble(eps.U) %>% group_by(across(starts_with("t"))) %>%
  summarize(p=n()/ns) %>% arrange(desc(p)) %>% head(top)
top.U[1:top,c(1:max1,n)]

# Table 4, rho2 posterior distribution (top most probable partitions for the variance)
top.V = as_tibble(eps.V) %>% group_by(across(starts_with("t"))) %>%
  summarize(p=n()/ns) %>% arrange(desc(p)) %>% head(top)
top.V[1:top,c(1:max2,n)]

# end points of the posterior modes of rho1 and rho2 distributions
eps_mu = c(47,79)
eps_s2 = c(51)


### mean and variance posterior estimates
# HPD interval
pHPD = 0.9
bmcp_mu.hpd <- MCMCsummary(bmcp.est$mu, HPD=T, hpd_prob=pHPD, Rhat=F, n.eff=F)
bmcp_s2.hpd <- MCMCsummary(bmcp.est$s2, HPD=T, hpd_prob=pHPD, Rhat=F, n.eff=F)
# moving variance
w = 2
mvar <- array(NA_real_,n)
for (i in (1+w):(n-w)) mvar[i] <- var(Y[(i-w):(i+w)]) 
# plot settings
xby=20;
ylim_mu = c(-6.5,12)
ylim_s2 = c(0,22)

# Figure 11(a)
bmcp_PE_mu <- ggplot() + ylab(expression(mu)) + xlab("") +
  theme(panel.background = element_blank()) +
  theme(axis.text = element_text(size=12, colour="black")) +
  theme(axis.line = element_line(color="black", size=.2)) +
  theme(axis.title.y = element_text(angle=90, vjust=0.5, size=18, face="bold")) +
  theme(axis.text.y = element_text(angle=90, hjust=0.5)) +
  theme(plot.margin=unit(c(.1,.2,0,0),"cm")) +
  scale_x_continuous(limits=c(1,n) , expand=c(0.01,1), breaks=seq(0,n,by=xby)) +
  scale_y_continuous(limits=ylim_mu, expand=c(0.01,0)) +
  # data
  geom_line(aes(x=i, y=obs), data=data.frame(i=1:n, obs=Y), col='black', size=.2) +
  # posterior estimates (mu)
  geom_point(aes(x=i, y=mu), data=data.frame(i=1:n, mu=bmcp_mu.hpd[[1]]), size=1, pch=20) +
  geom_line(aes(x=i, y=q05), data=data.frame(i=1:n, q05=bmcp_mu.hpd[[3]]), size=.5, lty='longdash') +
  geom_line(aes(x=i, y=q95), data=data.frame(i=1:n, q95=bmcp_mu.hpd[[4]]), size=.5, lty='longdash') +
  # mode(rho1|X)
  geom_vline(xintercept=eps_mu, col='gray10', size=.4, lty='dotted') +
  geom_vline(xintercept=eps_s2, col='gray10', size=.4, lty='dashed')

bmcp_PE_mu


# Figure 11(e)
bmcp_PE_s2 <- ggplot() + ylab(expression(sigma^2)) + xlab("") +
  theme(panel.background = element_blank()) +
  theme(axis.text = element_text(size=12, colour="black")) +
  theme(axis.line = element_line(color="black", size=.2)) +
  theme(axis.title.y = element_text(angle=90, vjust=0.5, size=18, face="bold")) +
  theme(axis.text.y = element_text(angle=90, hjust=0.5)) +
  theme(plot.margin=unit(c(.1,.2,0,0),"cm")) +
  scale_x_continuous(limits=c(1,n) , expand=c(0.01,1), breaks=seq(0,n,by=xby)) +
  scale_y_continuous(limits=ylim_s2, expand=c(0.01,0)) +
  # data
  geom_line(aes(x=i, y=obs), data=data.frame(i=1:n, obs=mvar), col='black', size=.2) +
  # posterior estimates (s2)
  geom_point(aes(x=i, y=s2), data=data.frame(i=1:n, s2=bmcp_s2.hpd[[1]]), size=1, pch=20) +
  geom_line(aes(x=i, y=q05), data=data.frame(i=1:n, q05=bmcp_s2.hpd[[3]]), size=.5, lty='longdash') +
  geom_line(aes(x=i, y=q95), data=data.frame(i=1:n, q95=bmcp_s2.hpd[[4]]), size=.5, lty='longdash') +
  # mode(rho2|X)
  geom_vline(xintercept=eps_mu, col='gray10', size=.4, lty='dotted') +
  geom_vline(xintercept=eps_s2, col='gray10', size=.4, lty='dashed')

bmcp_PE_s2


### posterior samples of the number of changes
N1 = rowSums(1-bmcp.est$u[,-n])
N2 = rowSums(1-bmcp.est$v[,-n])
tab.U = table(N1)/ns
tab.V = table(N2)/ns

# Figure 12(a)
barplot(tab.U, ylim=c(0,1))
title(xlab=expression(N[1]*" | "*X), family = "serif", line=2.5, cex.lab=1.5)

# Figure 12(b)
barplot(tab.V, ylim=c(0,1))
title(xlab=expression(N[2]*" | "*X), family = "serif", line=2.5, cex.lab=1.5)


### posterior probabilities of a change
prob.U = ifelse(colMeans(1-bmcp.est$u[,-n])==0,NA,colMeans(1-bmcp.est$u[,-n]))
prob.V = ifelse(colMeans(1-bmcp.est$v[,-n])==0,NA,colMeans(1-bmcp.est$v[,-n]))
th.mu = 0.2
th.s2 = 0.2

# Figure 13(a)
plot(prob.U, type='p', lwd=1.5, pch=20, col='black', cex=1.5,
     bty="n", yaxt="n", xlab='', ylab='probability of a change', ylim=c(0,1))
axis(1, lwd=2, lwd.ticks=2, las=1)
axis(2, lwd=2, lwd.ticks=2, las=1)
abline(h=seq(.2,.8,.2),col='gray70',lwd=1.5, lty='dashed')
abline(v=eps_mu,col='black',lwd=2.5, lty='dotted')
abline(v=eps_s2,col='black',lwd=1.5, lty='longdash')
leg <- rbind(which(prob.U>th.mu),prob.U[which(prob.U>th.mu)])
textxy(X=leg[1,1]-7, Y=leg[2,1], labs=round(leg[1,1],2), cex=.95, offset=.75, lwd=5)
textxy(X=leg[1,2]-7, Y=leg[2,2], labs=round(leg[1,2],2), cex=.95, offset=.75, lwd=5)
textxy(X=leg[1,3]+1, Y=leg[2,3], labs=round(leg[1,3],2), cex=.95, offset=.75, lwd=5)
textxy(X=leg[1,4]+1, Y=leg[2,4], labs=round(leg[1,4],2), cex=.95, offset=.75, lwd=5)

# Figure 13(b)
plot(prob.V, type= "p", lwd=1.5, pch=20, col='black', cex=1.5,
     bty="n", yaxt="n", xlab='', ylab='', ylim=c(0,1))
axis(1, lwd=2, lwd.ticks=2, las=1)
axis(2, lwd=2, lwd.ticks=2, las=1)
abline(h=seq(.2,.8,.2),col='gray70',lwd=1.5, lty='dashed')
abline(v=eps_mu,col='black',lwd=2.5, lty='dotted')
abline(v=eps_s2,col='black',lwd=1.5, lty='longdash')
leg <- rbind(which(prob.V>th.s2),prob.V[which(prob.V>th.s2)])
textxy(X=leg[1,]+2, Y=leg[2,], labs=round(leg[1,],2), cex=.95, offset=.75, lwd=5)


```



