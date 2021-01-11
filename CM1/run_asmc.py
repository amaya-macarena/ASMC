# -*- coding: utf-8 -*-
"""
Python 3.7 implementation of the Adaptive Sequential Monte Carlo method presented in Amaya et al. (2021) [under review].
This is a particle approach to infer the posterior probability density function and compute de evidence (marginal likelihood) 
introduced by in Zhou et al (2016) [algorithm 4]. 

This implementation (referred to as ASMC-DREAM) uses the code presented by Laloy et al. (2018a) for GAN-based probabilistic inversion using 
DREAMzs MCMC sampler (ter Braak and Vrugt, 2008; Vrugt, 2009; Laloy and Vrugt, 2012). The associated synthetic cross-hole ground 
penetrating radar (GPR) tomography data first-arrival times are calculated using the time2d algorithm by Podvin & Lecomte (1991).


Please write an email if you have any question and/or if you find a problem: macarena.amaya@unil.ch. 


===

Copyright (C) 2018  Eric Laloy

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

===                               

References:

Amaya, M., Linde, N., Laloy, E. Adaptive sequential Monte Carlo for posterior inference 
and model selection among complex geological priors encoded with deep generative neural networks 
[Submitted to Geophysical Journal International on December 4th, 2020].

Laloy, E., Hérault, R., Jacques, D., and Linde, N. 2018a. Training-image based geostatistical inversion using
a spatial generative adversarial neural network. Water Resources Research, 54, 381–406. https://doi.org/10.1002/2017WR022148.

Podvin, P. & Lecomte, I., 1991.  Finite difference computation of traveltimes in very contrasted velocity models: 
a massively parallel approach and its associated tools,Geophysical Journal International,105(1), 271–284.
          
Zhou,  Y.,  Johansen,  A.  M.,  &  Aston,  J.  A.,  2016.   Toward  automatic  model  comparison:  an adaptive sequential 
Monte Carlo approach, Journal of Computational and Graphical Statistics,69925(3), 701–726.                       


                                                                                                                                                                                                       
"""

import os
import time
import numpy as np
import shutil

work_dir=os.getcwd()

import asmc

#% Set random seed and case study

#rng_seed=12345

CaseStudy=1
Restart=False
    

if  CaseStudy==1: # Channelized conceptual model (CM1)
    
    # User defined parameters:
    
    seq=1 # Number of particles (N)
    thin=5  # Thinning parameter, rate for saving MCMC steps  
    steps=5 # Iterations per intermediate distribution (K)
    jr_scale=10 # Starting value for de proposal scale
    CESSf_div=0.99 # targeted CESS (CESS_op) 
    ESSf_div=0.5 # ESS treshold (ESS*) 
    AR_min=25.0 # Min acceptance rate before decreasing jr_scale
    jr_factor=0.2 # Fraction of the jr_scale that is decrease when the accepatnce rate gets lower than AR_min [0,1]
    
    
    
    ndraw=seq*3000000 # Set a high number of iterations to stop in case the beta sequence becomes too long due to a bad choice of CESS  
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
            
     
#% Run the ASMC-DREAM algorithm
            
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
    
    Sequences, Z, OutDiag, fx, MCMCPar, MCMCVar, beta_run, jr_seq, weig_seq, CESS_ev, increment, ESS_ev, evid_cont, evid_ev, weights_unnorm, new_weight_ev, weig_cont = q.sample(RestartFilePath=tmpFilePath)
    
    
    np.save('Sequences_states',Sequences) # Evolution of the states for every particle (latent parameters) and its likelihood
    np.save('AR',OutDiag.AR) # Acceptance Rate
    np.save('beta_seq',beta_run) # Sequence that defines the intermediate distributions (resulting from the adaptive procedure) 
    np.save('jr_ev',jr_seq) # Proposal scale evolution
    np.save('weig_ev',weig_seq) # Weights evolution	 	
    np.save('evid_ev',evid_ev) # Evidence evolution
