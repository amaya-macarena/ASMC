#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(Shiran Levy)s
"""
from ctypes import *
import numpy as np
# import time

def time2d_py(nx, ny, sourcez, nsource, nreceiver, xs, ys, rz, rx, s_model, spacing):
    
    s_model = s_model*spacing
    s_model=np.append(s_model,s_model[-1,:].reshape([1,nx]), axis = 0) # adds a ghost node to the bottom of the model (to avoid the boundry)
    ny = len(s_model[:,0]) 
    s_model = np.append(s_model, s_model[:,-1].reshape((ny,1)), axis = 1)
    nx = len(s_model[0,:]) 
    s_model = s_model.flatten('F')
    
    so_file = './time_2d.so'
    
    dll = CDLL(so_file)
    time_2d = dll.time_2d 
    # time_2d.restype = c_double
    time_2d.argtypes = [POINTER(nx*ny*c_double), POINTER(c_double), c_int, c_int, c_double,c_double, c_double, c_int] # specify c types argument for the function
    LP_c_double = POINTER(c_double)     # define a double array for c
    myVecType = nx*ny*c_double
    t = myVecType()
    hs = myVecType()
    
    for i in range(0,nx*ny):
        hs[i]=s_model[i]
    
    travel_time = np.zeros((len(sourcez),ny,nx))
    
    # tic = time.process_time()
    for jj in range(0,nsource):
        tt = time_2d(hs,t,nx,ny,xs,ys[jj],0.001,0)
        travel_time[jj,:,:]=np.reshape(t,(nx,ny)).T
    # toc = time.process_time()
    #calaculate first arrival time to the receivers:
    first_arrivals = np.ones((nsource,nreceiver))*-200
    
    for ii in range(0,nsource):
        first_arrivals[ii,:] = travel_time[ii,rz.astype(int),np.int(rx)]
    
    first_arrivals = first_arrivals.flatten()
    first_arrivals = first_arrivals.T
    
    return first_arrivals
