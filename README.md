GA_ANN
======

R scripts used to implement the GA_ANN input variable selection (IVS) algorithm described in Humphrey et. al. This wrapper IVS algorithm is a combination of a genetic algorithm (GA) search procedure with an artificial neural network (ANN) model. In this implementation, a simple 1-hidden node multilayer perceptron was utilised. The model training process is performed by means of a simulated annealing algorithm, which is used each time a new combination of inputs is evaluated. The GA adopted is a relatively simple variant outlined in Goldberg (1989). The user is referred to the original publication for further details regarding the GA_ANN algorithm.

Contents:
* `GA_Search.R`: implements a Genetic Algorithm (GA) search to the maximum of a user-provided objective (fitness) function.
* `GA_ANN_run.R`: run a k-fold cross-validation for an ensemble of Extra-Trees.


Humphrey, G.B., S. Galelli, H.R. Maier, A. Castelletti, G.C. Dandy and M.S. Gibbs (*submitted*), An evaluation framework for input variable selection algorithms for environmental data-driven models, Environmental Modelling & Software.
Goldberg, D.E., 1989. Genetic Algorithms in Search, Optimization and Machine Learning. Addison-Wesley Pub. Co., Reading, MA.

Copyright 2014 Greer Humphrey.
