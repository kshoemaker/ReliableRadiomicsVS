### plots for paper

library(R.matlab)
library(ggplot2)
setwd("~/Documents/GitHub/StatinMedicineProject/MatLabCode")

## plotting sensitivity plots
{
### plotting sens plots ####
alpha_0 <- c( -2.7500,   -2.5000,   -2.2500,   -2.0000,   -1.7500,   -1.5000,
              -1.2500,   -1.0000,   -0.7500,   -0.5000,   -0.2500,         0)

avg_ppi_real <- c(     0.7157,
                       0.7229,
                       0.9997,
                       0.9498,
                       0.9819,
                       1.0000,
                       0.9086,
                       1.0000,
                       0.9922,
                       1.0000,
                       1.0000,
                       1.0000)
avg_ppi_false <- c(    0.0006,
                       0.0014,
                       0.0020,
                       0.0038,
                       0.0050,
                       0.0077,
                       0.0112,
                       0.0164,
                       0.0224,
                       0.0298,
                       0.0411,
                       0.0662)

alpha0df <- data.frame(alpha_0, avg_ppi_false, avg_ppi_real)

ggplot(alpha0df, aes(x = alpha_0, y = avg_ppi_real)) + geom_line(colour = "black") +
  geom_line(aes(x = alpha_0, y = avg_ppi_false), colour = "red") +
  labs( title = "Sensitivity of PPI", subtitle = "Alpha 0") + xlab("alpha_0") + ylab("Average PPI") 
  

alpha_1 <- seq(from = 5, to = 7, by = 0.1)

avg_ppi_real <- c(    0.8717,    0.9716  ,  0.9486  ,  0.9167 ,   0.8810 ,
                      0.9428 ,   0.9904,    0.9828,    0.9646, 0.9870,
                      0.9819,    0.9854  ,  0.9994  ,  0.9933 ,   1.0000  ,
                      0.9959 ,   1.0000 ,   0.9792,
                      0.9167,    0.9838,    0.9501 )


avg_ppi_false <- c(    0.0042,
                       0.0045,
                       0.0046,
                       0.0040,
                       0.0039,
                       0.0042,
                       0.0052,
                       0.0042,
                       0.0051,
                       0.0048,
                       0.0050,
                       0.0051,
                       0.0050,
                       0.0059,
                       0.0060,
                       0.0055,
                       0.0053,
                       0.0050,
                       0.0056,
                       0.0054,
                       0.0060)


alpha1df <- data.frame(alpha_1, avg_ppi_false, avg_ppi_real)

ggplot(alpha1df, aes(x = alpha_1, y = avg_ppi_real)) + geom_line(colour = "black") +
  geom_line(aes(x = alpha_1, y = avg_ppi_false), colour = "red") +
  labs( title = "Sensitivity of PPI", subtitle = "Alpha 1" )  + xlab("alpha_1") + ylab("Average PPI") 


c <- seq(0.3,0.7, by = 0.05)

ppi_real <- c(    1.0000,
                  1.0000,
                  0.9965,
                  1.0000,
                  0.9819,
                  0.9657,
                  0.9199,
                  0.7268,
                  0.7362)

ppi_noise <- c( 0.0691,
                0.0335,
                0.0155,
                0.0086,
                0.0050,
                0.0037,
                0.0020,
                0.0015,
                0.0010)


c_df <- data.frame(c, ppi_real, ppi_noise)

ggplot(c_df, aes(x = c, y = ppi_real)) + geom_line(colour = "black") +
  geom_line(aes(x = c, y = ppi_noise), colour = "red") +
  labs( title = "Sensitivity of PPI", subtitle = "c") + xlab("c") + ylab("Average PPI") 
}


#### Chains #####
simdata <-readMat("simulation_two.mat")
simdata <-readMat("gammabi.mat")
gamma_chosen <- flatten(simdata$GammaBI)
step1 <- lapply(gamma_chosen, str_split, pattern = "[\\s]")
step2 <- lapply(step1, unlist)
step3 <- lapply(step2, as.numeric)
step4 <- lapply(step3, na.omit)
step5 <- lapply(step4, length)
chain1 <- unlist(step5)

gamma_chosen <- flatten(simdata$GammaBI2)
step1 <- lapply(gamma_chosen, str_split, pattern = "[\\s]")
step2 <- lapply(step1, unlist)
step3 <- lapply(step2, as.numeric)
step4 <- lapply(step3, na.omit)
step5 <- lapply(step4, length)
chain2 <- unlist(step5)

gamma_chosen <- flatten(simdata$GammaBI3)
step1 <- lapply(gamma_chosen, str_split, pattern = "[\\s]")
step2 <- lapply(step1, unlist)
step3 <- lapply(step2, as.numeric)
step4 <- lapply(step3, na.omit)
step5 <- lapply(step4, length)
chain3 <- unlist(step5)

plot(chain1, type = "l", xlab = "Iteration", ylab = "Variables Selected", main = "")
points(chain1, type = "l", col = 3)
points(chain3, type = "l", col = 4)

df <- data.frame(chain1, chain2, chain3, iteration = 1:80001)
df %>% ggplot(aes(x = iteration, y = chain1)) + 
  geom_line() + 
  geom_line(aes(x = iteration, y = chain2), color = "lightblue") + 
  geom_line(aes(x = iteration, y = chain3), color = "coral") + 
  labs(title  = "Number of Variables Selected", subtitle = "Three Simulation Chains")+ 
  xlab("Iterations") +
  ylab("Num. Variables Selected") 
ggsave("numvarsim.png")




informdata <- readMat('predictingHPVwprior.mat')
bi = informdata$bi

gamma_chosen <- flatten(informdata$GammaBI1)
step1 <- lapply(gamma_chosen, str_split, pattern = "[\\s]")
step2 <- lapply(step1, unlist)
step3 <- lapply(step2, as.numeric)
step4 <- lapply(step3, na.omit)
step5 <- lapply(step4, length)
chain1 <- unlist(step5)

gamma_chosen <- flatten(informdata$GammaBI2)
step1 <- lapply(gamma_chosen, str_split, pattern = "[\\s]")
step2 <- lapply(step1, unlist)
step3 <- lapply(step2, as.numeric)
step4 <- lapply(step3, na.omit)
step5 <- lapply(step4, length)
chain2 <- unlist(step5)

gamma_chosen <- flatten(informdata$GammaBI3)
step1 <- lapply(gamma_chosen, str_split, pattern = "[\\s]")
step2 <- lapply(step1, unlist)
step3 <- lapply(step2, as.numeric)
step4 <- lapply(step3, na.omit)
step5 <- lapply(step4, length)
chain3 <- unlist(step5)

plot(chain3, type = "l", xlab = "Iteration", ylab = "Variables Selected", main = "")
points(chain1, type = "l", col = 3)
points(chain2, type = "l", col = 4)

df <- data.frame(chain1, chain2, chain3, iteration = 1:100001)
df %>% ggplot(aes(x = iteration, y = chain3)) + 
  geom_line(color = "darkorchid") + 
  geom_line(aes(x = iteration, y = chain2), color = "turquoise") + 
  geom_line(aes(x = iteration, y = chain1), color = "coral") + 
  labs(title  = "Number of Variables Selected", subtitle = "Three Simulation Chains, Informative Prior")+ 
  xlab("Iterations") +
  ylab("Num. Variables Selected") 
ggsave("numvarwith.png")

df2 <- data.frame(MargGam = informdata$MargGam, Variable = 1:162)
df2 %>% ggplot(aes(x = Variable, y = MargGam)) + 
  geom_linerange(aes(x = Variable, ymin = 0, ymax = MargGam)) +
  labs(title = "Posterior Probability of Inclusion", subtitle = "Informative Prior") +
  ylab("PPI")
ggsave("ppiwith.png")



noninformdata <- readMat('predictingHPVwoprior.mat')
bi = noninformdata$bi

gamma_chosen <- flatten(noninformdata$GammaBI1)
step1 <- lapply(gamma_chosen, str_split, pattern = "[\\s]")
step2 <- lapply(step1, unlist)
step3 <- lapply(step2, as.numeric)
step4 <- lapply(step3, na.omit)
step5 <- lapply(step4, length)
chain1 <- unlist(step5)

gamma_chosen <- flatten(noninformdata$GammaBI2)
step1 <- lapply(gamma_chosen, str_split, pattern = "[\\s]")
step2 <- lapply(step1, unlist)
step3 <- lapply(step2, as.numeric)
step4 <- lapply(step3, na.omit)
step5 <- lapply(step4, length)
chain2 <- unlist(step5)

gamma_chosen <- flatten(noninformdata$GammaBI3)
step1 <- lapply(gamma_chosen, str_split, pattern = "[\\s]")
step2 <- lapply(step1, unlist)
step3 <- lapply(step2, as.numeric)
step4 <- lapply(step3, na.omit)
step5 <- lapply(step4, length)
chain3 <- unlist(step5)

plot(chain2, type = "l", xlab = "Iteration", ylab = "Variables Selected", main = "")
points(chain1, type = "l", col = 3)
points(chain3, type = "l", col = 4)

df <- data.frame(chain1, chain2, chain3, iteration = 1:100001)
df %>% ggplot(aes(x = iteration, y = chain1)) + 
  geom_line(color = "coral") + 
  geom_line(aes(x = iteration, y = chain2), color = "turquoise") + 
  geom_line(aes(x = iteration, y = chain3), color = "darkorchid") + 
  labs(title  = "Number of Variables Selected", subtitle = "Three Simulation Chains, Noninformative Prior")+ 
  xlab("Iterations") +
  ylab("Num. Variables Selected") 
ggsave("numvarwithout.png")


df2 <- data.frame(MargGam = noninformdata$MargGam, Variable = 1:162)
df2 %>% ggplot(aes(x = Variable, y = MargGam)) + 
  geom_linerange(aes(x = Variable, ymin = 0, ymax = MargGam)) +
  labs(title = "Posterior Probability of Inclusion", subtitle = "Non-Informative Prior") +
  ylab("PPI")
ggsave("ppiwithout.png")

chains = data.frame(betachain1 = simdata$Beta01.gam[,1],
           betachain2 = simdata$Beta01.gam[,2],
           betachain3 = simdata$Beta01.gam[,3], 
           Iteration = 1:51068)
ggplot(chains,aes(x = Iteration,  y = betachain1)) +  
  geom_line(colour = "slateblue", alpha = 0.5 ) +
    geom_line(aes(x = Iteration,  y = betachain3), colour = "turquoise", alpha = 0.5) + 
  geom_line(aes(x = Iteration,  y = betachain2), colour = "salmon", alpha = 0.5) + 
  geom_hline(aes(yintercept = -1), colour = "gray75", lty = "dashed") + 
  labs(title = "Class Specific Mean Mixing", subtitle = "Class 1, true mean = -1") +
  ylab("Parameter Value")  + ylim(-2,0)

chains = data.frame(betachain1 = simdata$Beta02.gam[,1],
                    betachain2 = simdata$Beta02.gam[,2],
                    betachain3 = simdata$Beta02.gam[,3], 
                    Iteration = 1:51068)
ggplot(chains,aes(x = Iteration,  y = betachain1)) + geom_line(colour = "slateblue", alpha = 0.5 ) +
  geom_line(aes(x = Iteration,  y = betachain3), colour = "turquoise", alpha = 0.5) + 
  geom_line(aes(x = Iteration,  y = betachain2), colour = "salmon", alpha = 0.5) + 
  geom_hline(aes(yintercept = 1), colour = "gray75", lty = "dashed") +
  labs(title = "Class Specific Mean Mixing", subtitle = "Class 2, true mean = 1") +
  ylab("Parameter Value") + ylim(0,2)

selected_chains = data.frame(beta1 = simdata$Beta01.gam, beta2 = simdata$Beta02.gam, Iteration = 1:41250)

ggplot(selected_chains, aes(x = Iteration, y = beta1.1)) + geom_line(colour = "slateblue") +
  geom_line(aes(x = Iteration, y = beta1.2), colour = "salmon") + 
  geom_line(aes(x = Iteration, y = beta1.3), colour = "gray25") 


## trying a gelman?

gelmans1 <- vector(length = length(informdata$ROI.sel))
for (k in 1:length(informdata$ROI.sel)){
i <- informdata$ROI.sel[k]
mu1 <- informdata$mu.1.mat[i,100001:250000]
mu2 <- informdata$mu.1.mat[i,350001:500000]
mu3 <- informdata$mu.1.mat[i,600001:750000]

gd <- gelman.diag(mcmc.list(list(mcmc(data = mu1, start = 1), mcmc(data = mu2, start = 1), mcmc(data = mu3, start = 1))), autoburnin = F)
gelmans1[k] <- gd$psrf[1]

}


gelmans2 <- vector(length = length(informdata$ROI.sel))
for (k in 1:length(informdata$ROI.sel)){
  i <- informdata$ROI.sel[k]
  mu1 <- informdata$mu.2.mat[i,100001:250000]
  mu2 <- informdata$mu.2.mat[i,350001:500000]
  mu3 <- informdata$mu.2.mat[i,600001:750000]
  
  gd <- gelman.diag(mcmc.list(list(mcmc(data = mu1, start = 1), mcmc(data = mu2, start = 1), mcmc(data = mu3, start = 1))), autoburnin = F)
  gelmans2[k] <- gd$psrf[1]
  
}


gd_all <- c(gelmans1, gelmans2)
summary(gd_all)


#######
# redo with SIM ppi 
df2 <- data.frame(MargGam = informdata$MargGam, Variable = 1:104)
df2 %>% ggplot(aes(x = Variable, y = MargGam)) + 
  geom_linerange(aes(x = Variable, ymin = 0, ymax = MargGam)) +
  labs(title = "Case Study: Posterior Probability of Inclusion", subtitle = "Proposed Model - Reliability Prior ") +
  ylab("PPI")
ggsave("ppiwith.png")

informMarg <- readMat("InformativeMargGam.mat")
df2 <- data.frame(MargGam = informMarg$MargGam, Variable = 1:104)
df2 %>% ggplot(aes(x = Variable, y = MargGam)) + 
  geom_linerange(aes(x = Variable, ymin = 0, ymax = MargGam)) +
  labs(title = "Simulation: Posterior Probability of Inclusion", subtitle = "Proposed Model - Reliability Prior ") +
  ylab("PPI")
ggsave("PPISimwith.png")

