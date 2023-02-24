# Spatial illness-death model with multivariate Leroux random effects via INLA (`R` code)

This is the `R` code associated to the paper:

[Llopis-Cardona, F., Armero, C., & Sanfélix-Gimeno, G. (2022). A Bayesian multivariate spatial approach for illness-death survival models. arXiv preprint arXiv:2210.07101.](https://arxiv.org/abs/2210.07101)

* `data_preparation.R`: function to prepare data according to INLA input formatting.
* `multiLeroux3d function.R`: INLA `rgeneric` function defining Leroux-multivariate lattent effects.
* `multiLeroux PREV2FO.R`: main `R` script containing INLA calls.
* `Transition probabilities INLA.R`: `R` functions for the calculation of posterior outcomes (transition probabilities and cumulative incidences).

