# -*- coding: utf-8 -*-
"""
Main asmc code. For control of the user-defined parameters and for running the inversion go to run_asmc.py 
"""
from __future__ import print_function

import numpy as np
import numpy.matlib as matlib
try:
    import cPickle as pickle
except:
    import pickle

import time

from scipy.stats import triang

from asmc_func import* # This imports both all ASMC, MCMC and inverse problem-related functions

import sys

from attrdict import AttrDict

MCMCPar=AttrDict()

MCMCVar=AttrDict()

Measurement=AttrDict()

OutDiag=AttrDict()

Extra=AttrDict()

class Sampler:

    

    
    def __init__(self, main_dir=None,CaseStudy=0,seq = 3,ndraw=10000,thin = 1,  nCR = 3, 
                 DEpairs = 3, parallelUpdate = 1.0, pCR=True,k=10,pJumpRate_one=0.2,
                 steps=100,savemodout=False, saveout=True,save_tmp_out=True,Prior='LHS',
                 DoParallel=True,eps=5e-2,BoundHandling='Reflect',
                 lik_sigma_est=False,parallel_jobs=4,jr_scale=1.0,rng_seed=123,it_b=10,jr_factor=0.2,CESSf_div=0.999993,ESSf_div=0.5,AR_min=25.0):
        
        
        
        
        self.CaseStudy=CaseStudy
        MCMCPar.seq = seq
        MCMCPar.ndraw=ndraw
        MCMCPar.thin=thin
        MCMCPar.nCR=nCR
        MCMCPar.DEpairs=DEpairs
        MCMCPar.parallelUpdate=parallelUpdate
        MCMCPar.Do_pCR=pCR
        MCMCPar.k=k
        MCMCPar.pJumpRate_one=pJumpRate_one
        MCMCPar.steps=steps
        MCMCPar.savemodout=savemodout
        MCMCPar.saveout=saveout  
        MCMCPar.save_tmp_out=save_tmp_out  
        MCMCPar.Prior=Prior
        MCMCPar.DoParallel=DoParallel
        MCMCPar.eps = eps
        MCMCPar.BoundHandling = BoundHandling
        MCMCPar.jr_scale=jr_scale
        MCMCPar.lik_sigma_est=lik_sigma_est
        Extra.n_jobs=parallel_jobs
        Extra.main_dir=main_dir
        np.random.seed(seed=None) 
        MCMCPar.rng_seed=rng_seed
        MCMCPar.it_b=it_b
        MCMCPar.jr_factor=jr_factor
        MCMCPar.AR_min=AR_min
        MCMCPar.ESSf_div=ESSf_div
        MCMCPar.CESSf_div=CESSf_div



        if self.CaseStudy==1:   
            
            ModelName='nonlinear_gpr_tomo'
            MCMCPar.lik=2
           
            MCMCPar.savemodout=True
            MCMCPar.lik_sigma_est=False
            
            MCMCPar.lb=np.ones((1,15))*-1
            MCMCPar.ub=np.ones((1,15))
            MCMCPar.n=MCMCPar.ub.shape[1]

            # Forward model:
            nx= 60     # Here x is the horizontal axis (number of columns) and not the number of rows
            ny = 125    # Here y is the vertical axis (number of rows) and not the number of columns
            finefac = 1     # not working for channelized models
            spacing = 0.1/finefac
            nx = np.int(60*finefac) 
            ny = np.int(125*finefac)
    	    # The x-axis is varying the fastest 
            sourcex = 0.1
            sourcez = np.arange(0.5,12+0.5,0.5)      #sources positions in meters
            receiverx = 5.9      
            receiverz = np.arange(0.5,12+0.5,0.5)    #receivers positions in meters
            xs = np.float32(sourcex/spacing)        #sources positions in model domain coordinates
            ys = sourcez/spacing       # divided by the spacing to get the domain coordinate                   
            rx = receiverx/spacing      #receivers positions in model domain coordinates
            rz = receiverz/spacing     # divided by the spacing to get receiverzthe domain coordinate  
            nsource = len(sourcez); nreceiver = len(receiverz)
            ndata=nsource*nreceiver
            data=np.zeros((ndata,4))
            x = np.arange(0,(nx/10)+0.1,0.1)                      
            y = np.arange(0,(ny/10)+0.1,0.1) 

            
            # Neural network (SGAN):
            DNN=AttrDict()
            DNN.nx=nx
            DNN.ny=ny
            DNN.zx=5
            DNN.zy=3
            DNN.nc=1
            DNN.nz=1
            DNN.depth=5
            DNN.threshold=True
            DNN.filtering=False
            DNN.cuda=False
            
            DNN.gpath=Extra.main_dir+'/netG.pth'

            from generator import Generator as Generator
            
            DNN.npx=(DNN.zx-1)*2**DNN.depth + 1
            DNN.npy=(DNN.zy-1)*2**DNN.depth + 1
            DNN.netG = Generator(cuda=DNN.cuda, gpath=DNN.gpath)
            for param in DNN.netG.parameters():
                param.requires_grad = False
            DNN.netG.eval()
            if DNN.cuda:
                DNN.netG.cuda()
            self.DNN=DNN
            
            # Load measurements

            Measurement.Sigma=1 # Measurement error is 1 ns
                   
            Measurement.MeasData=np.load('datatruemodel_sigma1.npy')        

            
            #Filter travel times with source-receiver angle > 45 degrees
            source=[]
            for i in range(1,25) :   
                s=np.ones(24)*(0.5*i)
                source=np.concatenate((source,s))
    
            rec=np.linspace(0.5, 12, 24)
            receiver=[]
    
            for ii in range(1,25) :
                receiver=np.concatenate((receiver,rec))
    
            dist=np.abs(source-receiver)
            ind=np.where(dist > 6)[0]
            
            Measurement.MeasData=np.delete(Measurement.MeasData,ind)
            Measurement.N=len(Measurement.MeasData) 
            
            
#           Define the beta list to choose from in the binary search, to crate an adaptative beta sequence.
            
            betanum=500
#
            beta_list=np.zeros(betanum-1)

            for i in range(1,betanum):
    
                beta_list[i-1]=i*2*10**-5
            
            MCMCPar.betainc_seq=beta_list
            

        MCMCPar.m0=10*MCMCPar.n
        
        self.MCMCPar=MCMCPar
        self.Measurement=Measurement
        self.Extra=Extra
        self.ModelName=ModelName
       
        
    def _init_sampling(self):
        
        start_time = time.time()
        
        Iter=self.MCMCPar.seq
        iteration=2
        iloc=0
        T=0
        

        Zinit=np.random.uniform(-1.,1.,(self.MCMCPar.m0+self.MCMCPar.seq,self.MCMCPar.lb.shape[1]))
            
            
        self.MCMCPar.CR=np.cumsum((1.0/self.MCMCPar.nCR)*np.ones((1,self.MCMCPar.nCR)))
        Nelem=np.floor(self.MCMCPar.ndraw/self.MCMCPar.seq)++self.MCMCPar.seq*2
        OutDiag.CR=np.zeros((np.int(np.floor(Nelem/self.MCMCPar.steps))+2,self.MCMCPar.nCR+1))
        OutDiag.AR=np.zeros((np.int(np.floor(Nelem/self.MCMCPar.steps))+2,2))
        OutDiag.AR[0,:] = np.array([self.MCMCPar.seq,-1])
        OutDiag.R_stat = np.zeros((np.int(np.floor(Nelem/self.MCMCPar.steps))+2,self.MCMCPar.n+1))
        pCR = (1.0/self.MCMCPar.nCR) * np.ones((1,self.MCMCPar.nCR))
        
        # Calculate the actual CR values based on pCR
        CR,lCR = GenCR(self.MCMCPar,pCR)  
        
        if self.MCMCPar.savemodout:
            self.fx = np.zeros((self.Measurement.N,np.int(np.floor(self.MCMCPar.ndraw/self.MCMCPar.thin))))
            MCMCVar.m_func = self.MCMCPar.seq     
        
        self.Sequences = np.zeros((np.int(np.floor(Nelem/self.MCMCPar.thin)),self.MCMCPar.n+2,self.MCMCPar.seq))
           
        self.MCMCPar.Table_JumpRate=np.zeros((self.MCMCPar.n,self.MCMCPar.DEpairs))
        for zz in range(0,self.MCMCPar.DEpairs):
            self.MCMCPar.Table_JumpRate[:,zz] = 2.38/np.sqrt(2 * (zz+1) * np.linspace(1,self.MCMCPar.n,self.MCMCPar.n).T)
        
        
        self.Z = np.zeros((np.floor(self.MCMCPar.m0 + self.MCMCPar.seq * (self.MCMCPar.ndraw - self.MCMCPar.m0) / (self.MCMCPar.seq * self.MCMCPar.k)).astype('int64')+self.MCMCPar.seq*100,self.MCMCPar.n+2))
        self.Z[:self.MCMCPar.m0,:self.MCMCPar.n] = Zinit[:self.MCMCPar.m0,:self.MCMCPar.n]

        X = Zinit[self.MCMCPar.m0:(self.MCMCPar.m0+self.MCMCPar.seq),:self.MCMCPar.n]
        
        ###X[0,:]=Extra.z_true.cpu().numpy().flatten()
        
        del Zinit
        
        # Run forward model, if any this is done in parallel
        if  self.CaseStudy > 0:
            if self.MCMCPar.lik_sigma_est==True: # The inferred sigma must always occupy the last position in the parameter vector
                fx0 = RunFoward(X[:,:-1],self.MCMCPar,self.Measurement,self.ModelName,self.Extra,DNN=self.DNN)
            else:
                fx0 = RunFoward(X,self.MCMCPar,self.Measurement,self.ModelName,self.Extra,DNN=self.DNN) #the one we use   
        
        # Compute likelihood from simulated data    
        of,log_p = CompLikelihood(X,fx0,self.MCMCPar,self.Measurement,self.Extra)

        X = np.concatenate((X,of,log_p),axis=1)
        Xfx = fx0
        
        if self.MCMCPar.savemodout==True:
            self.fx=fx0
        else:
            self.fx=None

        self.Sequences[0,:self.MCMCPar.n+2,:self.MCMCPar.seq] = np.reshape(X.T,(1,self.MCMCPar.n+2,self.MCMCPar.seq))

        # Store N_CR
        OutDiag.CR[0,:MCMCPar.nCR+1] = np.concatenate((np.array([Iter]).reshape((1,1)),pCR),axis=1)
        delta_tot = np.zeros((1,self.MCMCPar.nCR))

        # Compute the R-statistic of Gelman and Rubin  ##(This is only meaninful for mcmc without asmc)
        OutDiag.R_stat[0,:self.MCMCPar.n+1] = np.concatenate((np.array([Iter]).reshape((1,1)),GelmanRubin(self.Sequences[:1,:self.MCMCPar.n,:self.MCMCPar.seq],self.MCMCPar)),axis=1)
      
        self.OutDiag=OutDiag
        
        # Also return the necessary variable parameters
        MCMCVar.m=self.MCMCPar.m0
        MCMCVar.Iter=Iter
        MCMCVar.iteration=iteration
        MCMCVar.iloc=iloc; MCMCVar.T=T; MCMCVar.X=X
        MCMCVar.Xfx=Xfx; MCMCVar.CR=CR; MCMCVar.pCR=pCR
        MCMCVar.lCR=lCR; MCMCVar.delta_tot=delta_tot
        self.MCMCVar=MCMCVar
        
        if self.MCMCPar.save_tmp_out==True:
            with open('out_tmp'+'.pkl','wb') as f:
                 pickle.dump({'Sequences':self.Sequences,'Z':self.Z,
                 'OutDiag':self.OutDiag,'fx':self.fx,'MCMCPar':self.MCMCPar,
                 'MCMCVar':self.MCMCVar,'Measurement':self.Measurement,
                 'ModelName':self.ModelName,'Extra':self.Extra},f, protocol=pickle.HIGHEST_PROTOCOL)
    
        end_time = time.time()
    
        print("init_sampling took %5.4f seconds." % (end_time - start_time))
      
    def sample(self,RestartFilePath=None):
        start_time1b = time.time()
        
        if not(RestartFilePath is None):
            print('This is a restart')
            with open(RestartFilePath, 'rb') as fin:
                tmp_obj = pickle.load(fin)
            self.Sequences=tmp_obj['Sequences']
            self.Z=tmp_obj['Z']
            self.OutDiag=tmp_obj['OutDiag']
            self.fx=tmp_obj['fx']
            self.MCMCPar=tmp_obj['MCMCPar']
            self.MCMCVar=tmp_obj['MCMCVar']
            self.Measurement=tmp_obj['Measurement']
            self.ModelName=tmp_obj['ModelName']
            self.Extra=tmp_obj['Extra']
            del tmp_obj
            
            self.ndim=self.MCMCPar.n
#                
            self.MCMCPar.ndraw = 2 * self.MCMCPar.ndraw
            
            # Reset rng
            np.random.seed(np.floor(time.time()).astype('int'))
            
            # Extend Sequences, Z, OutDiag.AR,OutDiag.Rstat and OutDiag.CR
            self.Sequences=np.concatenate((self.Sequences,np.zeros((self.Sequences.shape))),axis=0)
            self.Z=np.concatenate((self.Z,np.zeros((self.Z.shape))),axis=0)
            self.OutDiag.AR=np.concatenate((self.OutDiag.AR,np.zeros((self.OutDiag.AR.shape))),axis=0)
            self.OutDiag.R_stat=np.concatenate((self.OutDiag.R_stat,np.zeros((self.OutDiag.R_stat.shape))),axis=0)
            self.OutDiag.CR=np.concatenate((self.OutDiag.CR,np.zeros((self.OutDiag.CR.shape))),axis=0)
      
            
        else:
            self._init_sampling()
        

        #Initialize some variables and arrays to store the results:

        prov_AR=30     # (no real meaning, just to not change the jr_scale on the first loop),   
        beta_run=[]
        increment=[]
        jr_seq=MCMCPar.jr_scale
        
        weig_seq=[]
        weig_cont=[]
        weights_unnorm=[]
        weig_unn=np.ones(self.MCMCPar.seq)
        norm_weight=np.ones(self.MCMCPar.seq)/self.MCMCPar.seq
        new_weight_ev=[]
        
        CESSf=self.MCMCPar.CESSf_div*self.MCMCPar.seq
        CESS_ev=[]
        ESS_ev=[]
        beta=0.
        beta_run=np.append(beta_run,beta)
        ind=0
        
        evid_cont=[]
        evid_ev=[]
        evid_evolution=1.
        
        eve_prev=np.arange(0,self.MCMCPar.seq,dtype=int)
        anc_prev=np.arange(0,self.MCMCPar.seq,dtype=int)
        
        eve_seq=[]
        eve_seq=np.append(eve_seq,eve_prev)
             
        end_time1b = time.time()
        print("The initialization of SAMPLE (before the main loop) took %5.4f seconds." % (end_time1b - start_time1b))



            
        # Main sampling loop  
       
        while self.MCMCVar.Iter < self.MCMCPar.ndraw:
            start_time1c = time.time()
            print('Iter =',self.MCMCVar.Iter)
                    
            #Calculate if it is necessary to modify jr_scale
            
            #Increase jr_scale if acceptance rate is too low
            if (prov_AR < self.MCMCPar.AR_min):
                self.MCMCPar.jr_scale = self.MCMCPar.jr_scale*(1-self.MCMCPar.jr_factor)
                
            
            # Initialize totaccept
            totaccept = 0

            jr_seq=np.append(jr_seq,self.MCMCPar.jr_scale)
                       
         
            # Loop a number of K mcmc steps for each intermediate distribution

            for gen_number in range(0,self.MCMCPar.steps):
                
                # Update T
                self.MCMCVar.T = self.MCMCVar.T + 1
                
                # Define the current locations and associated log-densities
                xold = np.array(self.MCMCVar.X[:self.MCMCPar.seq,:self.MCMCPar.n])
                log_p_xold = np.array(self.MCMCVar.X[:self.MCMCPar.seq,self.MCMCPar.n + 2-1])
                

                # Without replacement draw rows from Z for proposal creation
                R=np.random.permutation(self.MCMCVar.m)
                R=R[0:2 * self.MCMCPar.DEpairs * self.MCMCPar.seq]
                Zoff = np.array(self.Z[R,:self.MCMCPar.n])
             
                
                # Determine to do parallel direction or snooker update
                if (np.random.rand(1) <= self.MCMCPar.parallelUpdate):
                    Update = 'Parallel_Direction_Update'
                else:
                    Update = 'Snooker_Update'

                # Generate candidate points (proposal) in each chain using either snooker or parallel direction update
                xnew,self.MCMCVar.CR[:,gen_number] ,alfa_s = DreamzsProp(xold,Zoff,self.MCMCVar.CR[:,gen_number],self.MCMCPar,Update)
                              
                    
                # Get simulated data (done in parallel if DoParalel='True')
                if  self.CaseStudy > 0:
                    if self.MCMCPar.lik_sigma_est==True: # The inferred sigma must always occupy the last position in the parameter vector
                        fx_new = RunFoward(xnew[:,:-1],self.MCMCPar,self.Measurement,self.ModelName,self.Extra,DNN=self.DNN)
                    else:
                        fx_new = RunFoward(xnew,self.MCMCPar,self.Measurement,self.ModelName,self.Extra,DNN=self.DNN) #the one we are using
                else:
                    fx_new = RunFoward(xnew,self.MCMCPar,self.Measurement,self.ModelName,self.Extra)
                
              
                # Compute the likelihood of each proposal in each chain
                of_xnew,log_p_xnew = CompLikelihood(xnew,fx_new,self.MCMCPar,self.Measurement,self.Extra)
     
                # Calculate the Metropolis ratio
                accept = Metrop(self.MCMCPar,xnew,log_p_xnew,xold,log_p_xold,alfa_s,Extra,beta)

                # And update X and the model simulation
                idx_X= np.argwhere(accept==1);idx_X=idx_X[:,0]
                         
                if not(idx_X.size==0):
                     
                    self.MCMCVar.X[idx_X,:] = np.concatenate((xnew[idx_X,:],of_xnew[idx_X,:],log_p_xnew[idx_X,:]),axis=1)
                    self.MCMCVar.Xfx[idx_X,:] = fx_new[idx_X,:]
                
                                  
                # Check whether to add the current points to the chains or not?
                if self.MCMCVar.T == self.MCMCPar.thin:
                    # Store the current sample in Sequences
                    self.MCMCVar.iloc = self.MCMCVar.iloc + 1
                    self.Sequences[self.MCMCVar.iloc,:self.MCMCPar.n+2,:self.MCMCPar.seq] = np.reshape(self.MCMCVar.X.T,(1,self.MCMCPar.n+2,self.MCMCPar.seq))
                   
                   # Check whether to store the simulation results of the function evaluations
                    if self.MCMCPar.savemodout==True:
                        self.fx=np.append(self.fx,self.MCMCVar.Xfx,axis=0)
                        # Update m_func
                        self.MCMCVar.m_func = self.MCMCVar.m_func + self.MCMCPar.seq
                    else:
                        self.MCMCVar.m_func=None
                    # And set the T to 0
                    self.MCMCVar.T = 0

                # Compute squared jumping distance for each CR value
                if (self.MCMCPar.Do_pCR==True and self.MCMCVar.Iter < 0.1 * self.MCMCPar.ndraw):
                   
                    # Calculate the standard deviation of each dimension of X
                    r = matlib.repmat(np.std(self.MCMCVar.X[:,:self.MCMCPar.n],axis=0),self.MCMCPar.seq,1)
                    # Compute the Euclidean distance between new X and old X
                    delta_normX = np.sum(np.power((xold[:,:self.MCMCPar.n] - self.MCMCVar.X[:,:self.MCMCPar.n])/r,2),axis=1)
                                        
                    # Use this information to update delta_tot which will be used to update the pCR values
                    self.MCMCVar.delta_tot = CalcDelta(self.MCMCPar.nCR,self.MCMCVar.delta_tot,delta_normX,self.MCMCVar.CR[:,gen_number])

                # Check whether to append X to Z
                if np.mod((gen_number+1),self.MCMCPar.k) == 0:
                   
                    ## Append X to Z
                    self.Z[self.MCMCVar.m + 0 : self.MCMCVar.m + self.MCMCPar.seq,:self.MCMCPar.n+2] = np.array(self.MCMCVar.X[:,:self.MCMCPar.n+2])
                    # Update MCMCPar.m
                    self.MCMCVar.m = self.MCMCVar.m + self.MCMCPar.seq

                # Compute number of accepted moves
                totaccept = totaccept + np.sum(accept)

                # Update total number of MCMC iterations
                self.MCMCVar.Iter = self.MCMCVar.Iter + self.MCMCPar.seq
                
                
             
            curr_log_lik=np.array(self.MCMCVar.X[:self.MCMCPar.seq,self.MCMCPar.n + 2-1])    

            # Store acceptance rate
            self.OutDiag.AR[self.MCMCVar.iteration-1,:] = np.concatenate((np.array([self.MCMCVar.Iter]).reshape((1,1)), np.array([100 * totaccept/(self.MCMCPar.steps * self.MCMCPar.seq)]).reshape((1,1))),axis=1)
            
            prov_AR=100 * totaccept/(self.MCMCPar.steps * self.MCMCPar.seq)
            
            # Store probability of individual crossover values
            self.OutDiag.CR[self.MCMCVar.iteration-1,:self.MCMCPar.nCR+1] = np.concatenate((np.array([self.MCMCVar.Iter]).reshape((1,1)), self.MCMCVar.pCR),axis=1)
            
            # Is pCR updating required?
            if (self.MCMCPar.Do_pCR==True): #and self.MCMCVar.Iter < 0.1 * self.MCMCPar.ndraw):

                # Update pCR values
                self.MCMCVar.pCR = AdaptpCR(self.MCMCPar.seq,self.MCMCVar.delta_tot,self.MCMCVar.lCR,self.MCMCVar.pCR)

            # Generate CR values from current pCR values
            self.MCMCVar.CR,lCRnew = GenCR(MCMCPar,self.MCMCVar.pCR); self.MCMCVar.lCR = self.MCMCVar.lCR + lCRnew

            # Calculate Gelman and Rubin Convergence Diagnostic
            start_idx = np.maximum(1,np.floor(0.5*self.MCMCVar.iloc)).astype('int64')-1; end_idx = self.MCMCVar.iloc
            
            current_R_stat = GelmanRubin(self.Sequences[start_idx:end_idx,:self.MCMCPar.n,:self.MCMCPar.seq],self.MCMCPar)
            
            self.OutDiag.R_stat[self.MCMCVar.iteration-1,:self.MCMCPar.n+1] = np.concatenate((np.array([self.MCMCVar.Iter]).reshape((1,1)),np.array([current_R_stat]).reshape((1,self.MCMCPar.n))),axis=1)

            # Update number of complete generation loops
            self.MCMCVar.iteration = self.MCMCVar.iteration + 1

            if self.MCMCPar.save_tmp_out==True:
                with open('out_tmp'+'.pkl','wb') as f:
                    pickle.dump({'Sequences':self.Sequences,'Z':self.Z,
                    'OutDiag':self.OutDiag,'fx':self.fx,'MCMCPar':self.MCMCPar,
                    'MCMCVar':self.MCMCVar,'Measurement':self.Measurement,
                    'ModelName':self.ModelName,'Extra':self.Extra},f, protocol=pickle.HIGHEST_PROTOCOL)
            
            if beta>=1.:
                break
           
            #Binary search test for next beta
            
            next_beta,incr,CESS_found  = binary_search(curr_log_lik,CESSf,beta,norm_weight,self.MCMCPar.betainc_seq,self.MCMCPar.seq)
            
            CESS_ev=np.append(CESS_ev,CESS_found)
                 
            #Calculate importance weights for current beta 
    
            contribution = np.exp((next_beta - beta) * curr_log_lik)
            
            weig_cont=np.append(weig_cont,contribution)
            
            new_weight = np.multiply(norm_weight,contribution)
            
            new_weight_ev=np.append(new_weight_ev,new_weight)
            
            norm_weight = new_weight / np.sum(new_weight)
            
            weig_seq=np.append(weig_seq,norm_weight)
            
            weig_unn=np.multiply(weig_unn,contribution)
            
            weights_unnorm=np.append(weights_unnorm,weig_unn)
            
        
            evid=np.sum(new_weight)
            
            evid_cont=np.append(evid_cont, evid)
            
            evid_evolution=evid_evolution*evid
            
            evid_ev=np.append(evid_ev, evid_evolution)
            
                           
            ESS=(np.sum(norm_weight*contribution))**2 / np.sum(norm_weight**2*contribution**2)
            
            ESS_ev=np.append(ESS_ev,ESS)
            
            print('ESS=',ESS)
            
            
            if (ESS/self.MCMCPar.seq < MCMCPar.ESSf_div):
                
                print('Resample')
                
                Xres, ind, eve = resampling(norm_weight,self.MCMCPar.seq,self.MCMCVar.X,self.MCMCPar.n, anc_prev,eve_prev)
                
                eve_seq=np.append(eve_seq,eve)
                
                anc_prev=ind
                
                eve_prev=eve
                
                self.MCMCVar.X = Xres
                
                norm_weight= (1 / self.MCMCPar.seq) * np.ones(self.MCMCPar.seq)
              
            else:
                
                eve_seq=np.append(eve_seq,eve_prev)
                

            
            beta=next_beta
            
            increment=np.append(increment,incr)
            beta_run=np.append(beta_run,beta)   
                       

        
        start_time22 = time.time()
        # Remove zeros from pre-allocated variavbles if needed
        self.Sequences,self.Z,self.OutDiag,self.fx = Dreamzs_finalize(self.MCMCPar,self.Sequences,self.Z,self.OutDiag,self.fx,self.MCMCVar.iteration,self.MCMCVar.iloc,self.MCMCVar.pCR,self.MCMCVar.m,self.MCMCVar.m_func)
        
       
        if self.MCMCPar.saveout==True:
            with open('dreamzs_out'+'.pkl','wb') as f:
                pickle.dump({'Sequences':self.Sequences,'Z':self.Z,'OutDiag':self.OutDiag,'fx':self.fx,'MCMCPar':self.MCMCPar,'Measurement':self.Measurement,'Extra':self.Extra},f
                , protocol=pickle.HIGHEST_PROTOCOL)
        
        
        end_time22 = time.time()
        print("This saving took %5.4f seconds." % (end_time22 - start_time22))


        self.Sequences=self.Sequences[1:,:,:]        
        
        eve_seq=eve_seq.reshape(int(eve_seq.shape[0]/self.MCMCPar.seq),self.MCMCPar.seq)
        return self.Sequences, self.Z, self.OutDiag,  self.fx, self.MCMCPar, self.MCMCVar, beta_run, jr_seq, weig_seq, CESS_ev, increment, ESS_ev, evid_cont, evid_ev, weights_unnorm, new_weight_ev, weig_cont, eve_seq




        
