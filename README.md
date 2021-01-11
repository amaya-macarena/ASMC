# Adaptive Sequential Monte Carlo (ASMC) python codes for posterior inference and evidence computation 

Python 3.7 implementation of the Adaptive Sequential Monte Carlo method presented in Amaya et al. (under review), a particle approach to infer the posterior probability density function and compute de evidence (marginal likelihood) introduced by in Zhou et al (2016) [algorithm 4]. 

This implementation (referred to as ASMC-DREAM) uses the code presented by Laloy et al. (2018) for GAN-based probabilistic inversion using DREAMzs MCMC sampler (ter Braak and Vrugt, 2008; Vrugt, 2009; Laloy and Vrugt, 2012). The associated synthetic cross-hole ground penetrating radar (GPR) tomography data first-arrival times are calculated using the time2d algorithm by Podvin & Lecomte (1991).

## Test cases
CM1: binary channelized training image (CM1) (Zahner et al., 2016) 
CM2: tri-categorical training image representing braided-river aquifer deposits (Pirot et al., 2015).

## Codes and files

run_asmc.py : control the user-defined parameters and run the inversion. 

asmc.py : main asmc code.

asmc_func.py : contain the auxiliar functions called by asmc.py.

Eikonal_solver.py : forward solver (function for times2d)

.pht : SGAN

gen_from_z.py, generator.py, generate.py and generators.py: contain the SGAN generators

Z_trumodel_vector.npy: Latent space parameter values for the reference model. 

datatruemodel_sigma1.npy : first arrival times obtained with times2d for the reference model plus sigma=1ns gaussian random noise. 

noise_vector_sigma1.npy : the noise that was added to the first arrival times. 

forward_setup_0 folder: to use for forward solver parallel computation, contains the times2d codes.


## GAN-based representation
Pytorch and SGAN specifications


## Performing the inversion

Modify the user-defined parameters and run run_asmc.py. 


## References:

Amaya, M., Linde, N., Laloy, E. under review. Adaptive sequential Monte Carlo for posterior inference and model selection among complex geological priors 
encoded with deep generative neural networks [submitted to Geophysical Journal International on December 2020].

Laloy, E., Hérault, R., Jacques, D., and Linde, N. 2018. Training-image based geostatistical inversion using
a spatial generative adversarial neural network. Water Resources Research, 54, 381–406. https://doi.org/10.1002/2017WR022148.

Podvin, P. & Lecomte, I., 1991.  Finite difference computation of traveltimes in very contrasted velocity models: 
a massively parallel approach and its associated tools,Geophysical Journal International,105(1), 271–284

Ter Braak, C. J., & Vrugt, J. A. 2008. Differential evolution Markov chain with snooker updater and fewer chains. 
Statistics and Computing, 18(4), 435-446.

Vrugt, J. A., ter Braak, C., Diks, C., Robinson, B. A., Hyman, J. M., & Higdon, D. 2009. Accelerating Markov chain Monte Carlo simulation by
differential evolution with self-adaptive randomized subspace sampling. International Journal of Nonlin ear Sciences and Numerical Simu-
lation, 10(3), 273–290.
          
Zhou,  Y.,  Johansen,  A.  M.,  &  Aston,  J.  A.,  2016.   Toward  automatic  model  comparison:  an adaptive sequential 
Monte Carlo approach, Journal of Computational and Graphical Statistics,69925(3), 701–726.      


## Citation :

Amaya, M., Linde, N., Laloy, E. under review. Adaptive sequential Monte Carlo for posterior inference and model selection among complex geological priors 
encoded with deep generative neural networks [submitted to Geophysical Journal International on December 2020].

## License


## Contact
Macarena Amaya (macarena.amaya@unil.ch)
