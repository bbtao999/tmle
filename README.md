This is to replicate an example from the epidemiology field in which the author uses a simulated data set that compares the ols to the tmle methods to estimate the ATE, 
without knowing the true population parameter.

The author's example walk throughs the under-the-hood steps of the tmle manually, and show that it corrects the bias (at least partially) the estimation bias of OLS. 
However, if you manipulate the simulation example a bit, one may find that the tmle is not always working. A good bias correction requires specifying the propensity score matching model or the 
reponse function good enough (ideally using ML).
