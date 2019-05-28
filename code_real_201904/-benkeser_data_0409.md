

for propensity score -> use lasso regression

for outcome regression -> use lasso regression



all_data$Y == ic50.censored

W = {

country of virus

hxb (amino acid) = (which position, which amino acid) binarized multinomial



}



find the most common variant of amino-acid -> code as 1

other -> code as 0

=======================

## fit full blown methods

* ~~lasso Q+G~~ TMLE
* ~~lasso Q+G~~ one-step
* ~~lasso Q + HAL reduced G~~ TMLE
  * ~~take Q_hat and G_hat; just do the targeting~~
* ~~lasso Q + HAL reduced G one-step~~
  * ~~evaluate the EIF; add to the psi_n~~

~~not truncate the proensity scores~~

