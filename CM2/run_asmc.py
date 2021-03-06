# -*- coding: utf-8 -*-
"""
Python 3.7 implementation of the Adaptive Sequential Monte Carlo method presented in Amaya et al. (under review), a particle approach 
to infer the posterior probability density function and compute de evidence (marginal likelihood) introduced by in Zhou et al (2016) [algorithm 4]. 

This implementation (reffered to as ASMC-DREAM) uses the code presented by Laloy et al. (2018) for GAN-based probabilistic inversion using 
DREAMzs MCMC sampler (ter Braak and Vrugt, 2008; Vrugt, 2009; Laloy and Vrugt, 2012). The associated synthetic cross-hole ground 
penetrating radar (GPR) tomography data first-arrival times are calculated using the time2d algorithm by Podvin & Lecomte (1991).

In case you have a question or if you find a bug, please write me an email to macarena.amaya@unil.ch. 

===                               

References:

Amaya, M., Linde, N., Laloy, E. (under review). Adaptive sequential Monte Carlo for posterior inference and model selection among complex geological priors 
encoded with deep generative neural networks [submitted to Geophysical Journal International on December 2020].

Laloy, E., & Vrugt, J. A. (2012). High‐dimensional posterior exploration of hydrologic models using multiple‐try DREAM (ZS) and high‐performance computing. 
Water Resources Research, 48(1).

Laloy, E., Hérault, R., Jacques, D., and Linde, N. (2018). Training-image based geostatistical inversion using
a spatial generative adversarial neural network. Water Resources Research, 54, 381–406. https://doi.org/10.1002/2017WR022148.

Podvin, P. & Lecomte, I., (1991).  Finite difference computation of traveltimes in very contrasted velocity models: 
a massively parallel approach and its associated tools,Geophysical Journal International,105(1), 271–284

Ter Braak, C. J., & Vrugt, J. A. (2008). Differential evolution Markov chain with snooker updater and fewer chains. 
Statistics and Computing, 18(4), 435-446.

Vrugt, J. A., ter Braak, C., Diks, C., Robinson, B. A., Hyman, J. M., & Higdon, D. (2009). Accelerating Markov chain Monte Carlo simulation by
differential evolution with self-adaptive randomized subspace sampling. International Journal of Nonlin ear Sciences and Numerical Simu-
lation, 10(3), 273–290.
          
Zhou,  Y.,  Johansen,  A.  M.,  &  Aston,  J.  A.,  (2016).   Toward  automatic  model  comparison:  an adaptive sequential 
Monte Carlo approach, Journal of Computational and Graphical Statistics,69925(3), 701–726.                                                                                                                                                                                                                             
"""


import os
import time
import numpy as np
import shutil

work_dir=os.getcwd()

import asmc

#% Set random seed and case study

CaseStudy=2   # Tricategorical braided aquifer (CM2)
Restart=False
 
    
if  CaseStudy==2: 

    # User defined parameters:
    
    seq=40 # Number of particles (N)
    thin=60  # Thinning parameter, rate for saving MCMC steps  
    steps=60 # Iterations per intermediate distribution (K)
    jr_scale=10 # Starting value for de proposal scale
    CESSf_div=0.999996 # targeted CESS (CESS_op) 
    ESSf_div=0.5 # ESS treshold (ESS*) 
    AR_min=25.0 # Min acceptance rate before decreasing jr_scale
    jr_factor=0.2 # Fraction of the jr_scale that is decrease when the accepatnce rate gets lower than AR_min [0,1]
    
    

    ndraw=seq*1000000 # Set a high number of iterations to stop in case the beta sequence becomes too long due to a bad choice of CESS  
    it_b=steps 
    
    #Decide if to run forward solver in parallel
    
    
    DoParallel=True
    parallel_jobs=seq
    MakeNewDir=True
    
    if MakeNewDir==True:
        src_dir=work_dir+'/forward_setup_0'
        for i in range(1,parallel_jobs+1):
            dst_dir=work_dir+'/forward_setup_'+str(i)
            if os.path.exists(dst_dir):
                shutil.rmtree(dst_dir)
            shutil.copytree(src_dir,dst_dir)
 
#% Run the DREAMzs algorithm
if __name__ == '__main__':
    
    #start_time = time.time()

    q=asmc.Sampler(main_dir=work_dir,CaseStudy=CaseStudy,seq=seq,ndraw=ndraw,parallel_jobs=seq,steps=steps,
                   parallelUpdate = 1,pCR=False,thin=thin,nCR=3,DEpairs=1,pJumpRate_one=0.2,BoundHandling='Fold',
                   lik_sigma_est=False,DoParallel=DoParallel,jr_scale=jr_scale,it_b=it_b,jr_factor=jr_factor,CESSf_div=CESSf_div,ESSf_div=ESSf_div,AR_min=AR_min)
    
    print("Iterating")
    
    if Restart:
        tmpFilePath=work_dir+'/out_tmp.pkl'
    else:
        tmpFilePath=None 
    
    Sequences, Z, OutDiag, fx, MCMCPar, MCMCVar, beta_run, jr_seq, weig_seq, CESS_ev, increment, ESS_ev, evid_cont, evid_ev, weights_unnorm, new_weight_ev, weig_cont, eve_seq = q.sample(RestartFilePath=tmpFilePath)
    
    #end_time = time.time()
    
    #print("This sampling run took %5.4f seconds." % (end_time - start_time))
    
    np.save('Sequences_states',Sequences) # Evolution of the states for every particle (latent parameters) and its likelihood
    np.save('AR',OutDiag.AR) # Acceptance Rate
    np.save('beta_seq',beta_run) # Sequence that defines the intermediate distributions (resulting from the adaptive procedure) 
    np.save('jr_ev',jr_seq) # Proposal scale evolution
    np.save('weig_ev',weig_seq) # Weights evolution	 	
    np.save('evid_ev',evid_ev) # Evidence evolution
    np.save('eve_ev',eve_seq) # Eve index evolution
