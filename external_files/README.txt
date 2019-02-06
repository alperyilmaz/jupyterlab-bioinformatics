This archive contains code used in the paper:

A Nonparametric Significance Test for Sampled Networks

The code is in python and relies on the Networkx library. Networkx can be
installed via pip, or from the libraries website (https://networkx.github.io/).

The code is in Python 2.7, but the use of dictionary comprehensions may prevent
use in earlier versions of python. 

Please note that the multicore version of the code assumes that it can run on
20 virtual cores, if this is not the case, either do not run it or change the
constants in the code.

The files in this folder are as follows:

sampling.py
    Contains code that performs graph sampling.

redundantSimple.py
    Contains code to compute the smallest seed list that generates the same network.

networkmetricsets.py
    Contains wrappers to compute the statistics used in the paper.

random_bio_seeds.py
    Contains code to perform the Monte Carlo sampling of seeds over degrees,
    also contains wrappers around the sampling techniques code to facilitate
    its use.

If you have nay questions please email, ande.elliott@gmail.com
