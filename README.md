
GWMCMC
======

GWMCMC is an implementation of the Goodman and Weare 2010 Affine
invariant ensemble Markov Chain Monte Carlo (MCMC) sampler. MCMC sampling
enables bayesian inference. The problem with many traditional MCMC samplers
is that they can have slow convergence for badly scaled problems, and that
it is difficult to optimize the random walk for high-dimensional problems.
This is where the GW-algorithm really excels as it is affine invariant. It
can achieve much better convergence on badly scaled problems. It is much
simpler to get to work straight out of the box, and for that reason it
truly deserves to be called the MCMC hammer. 


Authors: [Aslak Grinsted](http://www.glaciology.net) 


	
## Licensing

The majority of the code is licensed under a very permissive MIT license, but some routines and example data are licensed under other terms. See licensing details in LICENSE.txt and individual files. 


## Acknowledgements

This software has been developed at [Centre for Ice and Climate](http://www.iceandclimate.nbi.ku.dk), Niels Bohr Institute, University of Copenhagen. It is partly inspired by emcee for python, but not modelled after it. 

