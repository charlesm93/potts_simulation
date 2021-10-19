# potts_simulation
Code for the paper "Simulating Ising and Potts models at critical and cold temperatures using auxiliary Gaussian variables." The R code defining the algorithms and used to run the numerical experiments is located under the `scripts` directory and is organised as follows.

### Algorithms

* `algorithms.r`: specifies (specialized) sampling algorithms for the Ising model (q = 2).
* `algorithms_potts.r`: specifies sampling algorithms for the Potts model.
* `algorithm_tempering.r`: specifies the tempering algorithm, which acts as a wrapper around a sampling algorithm.
* `tools.r`: defines a series of handy functions to support the above files.
* `potts_unit_test.r`: runs small Ising and Potts models, and checks that the algorithms produce results in agreement with one another.

### Numerical experiments

The `performance` files can be used to reproduce our numerical experiments, and run further experiments.
At the top of the file adjust your the directories to your setup.
Specify which algorithms would you like to run, the type of graphs, and the temperatures.
The model then runs all sampling algorithms on the specified Potts models and prints out MCMC summaries in an output file.
These summaries include estimated mean, Monte Carlo error, estimated variance, Rhat, effective sample size and effective sample size per second for various quantities of interest.
Be sure to specify the location of the output files.

* `performance_test_v2.r`: run numerical experiments for the Ising model.
* `performance_potts_test.r`: run numerical experiments for the Potts model. Also covers the Ising model case.
* `performance_test_temp.r`: run numerical experiments with the tempering algorithm.
* `make_plots_potts.r`: the file used to create the plots that appear in the paper. To use this file, outputs from the numerical experiments need to be created.

