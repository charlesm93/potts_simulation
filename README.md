# potts_simulation
Code for the paper "Simulating Ising and Potts models at critical and cold temperatures using auxiliary Gaussian variables" by Charles Margossian, Chenyang Zhong and Sumit Mukherjee. The paper is still in preparation but a (now outdated) preprint can be found at: https://arxiv.org/abs/2110.10801.

 The R code defining the algorithms and used to run the numerical experiments is located under the `scripts` directory and is organised as follows.

### Algorithms

* `algorithms.r`: specifies (specialized) sampling algorithms for the Ising model (q = 2).
* `algorithms_potts.r`: specifies sampling algorithms for the Potts model.
* `algorithm_tempering.r`: specifies the tempering algorithm, which acts as a wrapper around a sampling algorithm.
* `tools.r`: defines a series of handy functions to support the above files.
* `potts_unit_test.r`: runs small Ising and Potts models, and checks that the algorithms produce results in agreement with one another.

### Numerical experiments

All experiments can be run using the functions in `run_simulations.r`. The results are saved in the `deliv_v3` directory and can then be analyzed using `evaluate_simulations.r` to generate figures.
