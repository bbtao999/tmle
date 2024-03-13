#------------------------------------------------------------------------------#
#  this is to create a simulation exercise to test the tmle                    #
# replicate https://migariane.github.io/TMLE.nb.html                           #
#------------------------------------------------------------------------------#

#why tmle?
#Evidence shows that TMLE provides the less unbiased ATE estimate compared with 
#other double-robust estimators (Neugebauer and Laan, 2005), (Laan and Rose, 2011) 
#such as the combination of regression adjustment with inverse probability of treatment weighting (IPTW-RA) 
#and the augmented inverse probability of treatment weighting (AIPTW).

#the example is just for demonstration purpose because technically
#we should fit the PPS and response model in a better way for the theoretical 
#properties of the estimators to hold
library(tidyverse)

#----------------------------------------------------------------#
#0. create simu data #
#----------------------------------------------------------------#

simu_data<-function(rseed,n=1000){
  set.seed(rseed)
  #define covariates w1 - w4
  # w1 <- rbinom(n, size=1, prob=0.5) #gender
  # w2 <- rbinom(n, size=1, prob=0.65) #age
  # w3 <- round(runif(n, min=0, max=4), digits=3) #disease condition scale 1-4
  # w4 <-rnorm(n, mean=0,sd=1) #standardized income
  
  w1 <- rbinom(n, size=1, prob=0.5)
  w2 <- rbinom(n, size=1, prob=0.65)
  w3 <- round(runif(n, min=0, max=4), digits=3)
  w4 <- round(runif(n, min=0, max=5), digits=3)
  
  #Treatment
  # A  <- rbinom(n, size=1, prob= plogis(-0.4 + 0.2*w2 + 0.15*w3 + 0.2*w4 + 0.15*w2*w4)) #using logistic regression to make sure the prob is within 0 and 1
  A  <- rbinom(n, size=1, prob= plogis(-0.4 + 0.2*w2 + 0.15*w3 + 0.2*w4 + 0.15*w2*w4))
  
  # counterfactual
  # Y.1 <- rbinom(n, size=1, prob= plogis(-1 + 1 -0.1*w1 + 0.3*w2 + 0.25*w3 - 0.2*exp(w4) + 0.15*w2*w4))
  # Y.0 <- rbinom(n, size=1, prob= plogis(-1 + 0 -0.1*w1 + 0.3*w2 + 0.25*w3  - 0.2*exp(w4) + 0.15*w2*w4))
  
  Y.1 <- rbinom(n, size=1, prob= plogis(-1 + 1 -0.1*w1 + 0.3*w2 + 0.25*w3 + 0.2*w4 + 0.15*w2*w4))
  Y.0 <- rbinom(n, size=1, prob= plogis(-1 + 0 -0.1*w1 + 0.3*w2 + 0.25*w3 + 0.2*w4 + 0.15*w2*w4))
  
  
  # Observed outcome
  Y <- Y.1*A + Y.0*(1 - A)
  
  # return data.frame
  data.frame(w1, w2, w3, w4, A, Y, Y.1, Y.0)
}


ObsData <- simu_data(7777,n=10000)
True_Psi <- mean(ObsData$Y.1-ObsData$Y.0);
cat(" True_Psi:", True_Psi)


#----------------------------------------------------------------#
#1. Stat Tests                                                   #
#----------------------------------------------------------------#

#a. ols
Bias_Psi_ols <- lm(data=ObsData, Y~ A + w1 + w2 + w3 + w4)
cat("\n")
cat("\n Bias_Psi_ols:",summary(Bias_Psi_ols)$coef[2, 1])

#----------------------------------------------------------------#
#2. TMLE Tests  (manual steps)                                   #
#----------------------------------------------------------------#

#------------------------------------------------#
#step 1:
#use super learner library to estimate Q_n^0(A,W)
#------------------------------------------------#

#for simplicity, we use the logistic (sigmoid) to estimate Q_n^0
#use all data (without train/test split)
m  <- glm(Y ~ A + w1 + w2 + w3 + w4, family=binomial, data=ObsData)

Y  <- ObsData$Y
A  <- ObsData$A
w1 <- ObsData$w1
w2 <- ObsData$w2
w3 <- ObsData$w3
w4 <- ObsData$w4

Q  <- as.data.frame(cbind(Q_AW = predict(m),
            Q_1W = predict(m, newdata=data.frame(A = 1, w1, w2, w3, w4)),
            Q_0W = predict(m, newdata=data.frame(A = 0, w1, w2, w3, w4))))

#now we have the log odds as the dependent variable, then convert back to P
Q_A1 <- exp(Q$Q_1W)/(1+exp(Q$Q_1W))
Q_A0 <- exp(Q$Q_0W)/(1+exp(Q$Q_0W))

psi_hat<-mean(Q_A1 - Q_A0)
cat("\n Q0:", psi_hat)

#------------------------------------------------#
#step 2: estimate the PPS g^0(A,W)
#use logistic (sigmoid)
#------------------------------------------------#

g <- glm(A ~ w2 + w3 + w4, data=ObsData, family = binomial)
g1W = predict(g, type ="response")
cat("\n Propensity score = g1W","\n")
summary(g1W)

#------------------------------------------------#
#step 3: HAW and epsilon (initial guess)                     
#------------------------------------------------#
ObsData_reg<-cbind(ObsData,g1W)
ObsData_reg<-ObsData_reg %>% 
  mutate(h=A/g1W - (1 - A) / (1 - g1W)) #or the H(A,W) function)

#offset function ensures the coef of this variable in the offset function to be fixed at 1
epsilon <- coef(glm(Y ~ -1 + h + offset(Q[,"Q_AW"]), data=ObsData_reg, family = binomial))
cat("\n Epsilon:",epsilon)

#------------------------------------------------#
#step 4:update the Q*_n bar                   
#------------------------------------------------#

#Afterwards, the estimated probability of the potential outcomes is updated by 
#the substitution parameters (ϵ0^,ϵ1^)
#The substitution update is performed by setting A = 0 and A = 1 for each subject 
#in the initial estimate probability of the potential outcomes bar(Q0(0,W)),bar(Q0(1,W))
#, as well as in the clever covariates H0(0,W) and H1(1,W)

Qstar<-Q %>% mutate(Q1Wstar=plogis(Q_1W+epsilon/g1W),
                    Q0Wstar=plogis(Q_0W-1*epsilon/(1-g1W)))%>%
             select(Q1Wstar,Q0Wstar)
  
Psi_new <- mean(Qstar[,"Q1Wstar"] - Qstar[,"Q0Wstar"]);
cat("TMLE_Psi:", Psi_new)

#------------------------------------------------#
#step 5: inferences              
#------------------------------------------------#

#the idea is to compute the IC (influence curve)
#b/c the variance of the psi is the variance of IC /n

Q  <- as.data.frame(Q)
Qstar <- as.data.frame(Qstar)

#compute IC
IC <- ObsData_reg$h*(ObsData_reg$Y-plogis(Q$Q_AW)) 
      + plogis(Qstar$Q_1W - Qstar$Q_0W) - Psi_new

summary(IC)

n <- nrow(ObsData_reg)
varHat.IC <- var(IC)/n
cat("TMLE_Psi_Variance:", varHat.IC)
cat("\n point estimate of ATE tmle:", Psi_new,
    "95%CI:", c(Psi_new-1.96*sqrt(varHat.IC),  Psi_new+1.96*sqrt(varHat.IC)))


#----------------------------------------------------------------#
#3. AIPW (also doubly robust)                                    #
#----------------------------------------------------------------#
#instead of having a Q*, it just uses the Q
#the second part is the augmented part

AIPW<-mean(plogis(Q$Q_1W)-plogis(Q$Q_0W)+
  ObsData_reg$A/g1W*(ObsData_reg$Y-plogis(Q$Q_1W))-
  (1-ObsData_reg$A)/(1-g1W)*(ObsData_reg$Y-plogis(Q$Q_0W)))
AIPW

#----------------------------------------------------------------#
#4. TMLE Tests  R-TMLE                                           #
#----------------------------------------------------------------#
library(tmle)
w <- subset(ObsData_reg, select=c(w1,w2,w3,w4))
Y<-ObsData_reg$Y
A<-ObsData_reg$A

tmle <- tmle(Y, A, W=w)
cat("TMLER_Psi:", tmle$estimates[[2]][[1]],
    ";","95%CI(", tmle$estimates[[2]][[3]],")")

#----------------------------------------------------------------#
#5. TMLE Tests  R-TMLE                                           #
#----------------------------------------------------------------#

#In addition to the default algorithms implemented in the R-tmle package, 
#we can even decrease more the bias of our estimation by calling more efficient 
# machine learing algorithms, such as generalized additive models, Random Forest,
# Recursive Partitioning and Regression Trees (it is highly recomended to include the the highly adaptive Lasso [HAL9001] in your SL library but for computing efficiency we did not include it here)
# library(hal9001)

SL.TMLER.Psi <- tmle(Y=Y, A=A, W=w, family="binomial", 
                     Q.SL.library = c("SL.glm", "SL.step", "SL.glm.interaction",
                                      "SL.gam", "SL.ranger"),
                     g.SL.library = c("SL.glm", "SL.step", "SL.glm.interaction", "SL.gam", "SL.ranger"))

cat("SL.TMLER.Psi:", SL.TMLER.Psi$estimates[[2]][[1]],";","95%CI(", SL.TMLER.Psi$estimates[[2]][[3]],")")

#----------------------------------------------------------------#
#6. Result Comparison                                            #
#----------------------------------------------------------------#

cat(" True_Psi:", True_Psi)

cat("\n Bias_Psi_ols:",summary(Bias_Psi_ols)$coef[2, 1])

cat("\n point estimate of ATE tmle:", Psi_new,
    "95%CI:", c(Psi_new-1.96*sqrt(varHat.IC),  Psi_new+1.96*sqrt(varHat.IC)))

cat("\n TMLER_Psi:", tmle$estimates[[2]][[1]],
    ";","95%CI(", tmle$estimates[[2]][[3]],")")

cat("\n SL.TMLER.Psi:", SL.TMLER.Psi$estimates[[2]][[1]],";",
    "95%CI(", SL.TMLER.Psi$estimates[[2]][[3]],")")

