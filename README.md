GA_ANN
======

R scripts used to implement the GA_ANN input variable selection (IVS) algorithm described in Galelli et. al. (2014). This wrapper IVS algorithm is a combination of a genetic algorithm (GA) search procedure with an artificial neural network (ANN) model. In this implementation, a simple 1-hidden node multilayer perceptron was utilised. The model training process is performed by means of a simulated annealing algorithm, which is used each time a new combination of inputs is evaluated. The GA adopted is a relatively simple variant outlined in Goldberg (1989). Further details of this implementation of the GA_ANN algorithm can be found in:

Galelli S., Humphrey G.B., Maier H.R., Castelletti A., Dandy G.C. and Gibbs M.S. (2014)  An evaluation framework for input variable selection algorithms for environmental data-driven models, *Environmental Modelling and Software*, 62, 33-51, DOI: 10.1016/j.envsoft.2014.08.015.

Contents:
* `GA_Search.R`: implements a Genetic Algorithm (GA) search to the maximum of a user-provided objective (fitness) function.
* `GA_ANN_run.R`: run the GA_ANN IVS algorithm to select the optimal inputs for a given set of input data.
* `inp_dat.csv`: an example input data file. Column 1 contains an array of data labels or IDs (e.g. dates on which data were recorded); columns 2 to P+1 contain the P candidate input variables; and column P+2 contains the response variable, while the rows are data points. The first row contains the variable names.

To run the GA_ANN algorithm, the following command should be used:

`R --args [`*filename*`] [`*out_dir*`] < GA_ANN_run.R`

where *filename* is the name of the name of the input data file (including path) and *out_dir* is the name of the output directory (i.e. the directory to which results will be written. This name should *NOT* include the whole path).  


Goldberg, D.E., 1989. Genetic Algorithms in Search, Optimization and Machine Learning. Addison-Wesley Pub. Co., Reading, MA.

Copyright 2014 Greer Humphrey.
