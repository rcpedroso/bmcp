# `bmcp`
R package for the Normal data application of the Bayesian multipartition change point model introduced in [Pedroso et al (2021)](https://arxiv.org/abs/2107.11456).

# Install
To install `bmcp`, consider the code

```R
install.packages("devtools")
devtools::install_github("rcpedroso/bmcp")
```
If installation fails, please report the problem to ricardocunhap@gmail.com. Any other problems or comments can also be reported!

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


### MCMC run (run time: ~30 seconds using a computer with an Intel® Core i7-7500U/2.9GHz with 16Gb of RAM)
burn = 3e4
ns = 2e4
niter = burn + ns

set.seed(1000)
bmcp.est <- bmcp(burn=burn, N=niter, X=Y,
                 alpha1=1, beta1=1, alpha2=1, beta2=1,
                 a=0.1, d=2.1,
                 mu0=0, s02=100)


### posterior distributions of partitions rho1 and rho2

# function u_to_index
# input: any matrix M of 0/1 elements
# output: if M has m rows and (n-1) columns, the function returns a matrix with m rows
# and n columns where the elements indicate the positions of the 0's in each row
# and column n indicates the number of 0's in each respective row.
# If M is the matrix where each line is a sample of some specific partition, the function
# returns the end points and the number of end points of each sample.

# posterior samples of the end points
eps.U = u_to_index(bmcp.est$u[,-n])
eps.V = u_to_index(bmcp.est$v[,-n])

colnames(eps.U) <- c(paste0('t',1:(n-1)),'N1')
colnames(eps.V) <- c(paste0('t',1:(n-1)),'N2')
max1 = max(eps.U[,'N1'])
max2 = max(eps.V[,'N2'])

top = 5

# rho1 posterior distribution
top.U = as_tibble(eps.U) %>% group_by(across(starts_with("t"))) %>%
  summarize(p=n()/ns) %>% arrange(desc(p)) %>% head(top)
top.U[1:top,c(1:max1,n)]

# rho2 posterior distribution
top.V = as_tibble(eps.V) %>% group_by(across(starts_with("t"))) %>%
  summarize(p=n()/ns) %>% arrange(desc(p)) %>% head(top)
top.V[1:top,c(1:max2,n)]

# end points of the posterior modes of rho1 and rho2
eps_mu = c(47,79)
eps_s2 = c(51)

```



