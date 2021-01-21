# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 11:54:17 2018

@author: elaloy
"""

import torch
import torch.nn as nn
import torch.nn.parallel


class G(nn.Module):
    def __init__(self, gpath, nc = 1, nz = 1, ngf = 64, gfs = 5, ngpu = 1,cuda=True):
        super(G, self).__init__()
        self.ngpu = ngpu
   
        if gfs < 5:
            out_pad=1
        else:
            out_pad=0

        self.main = nn.Sequential(
             
                nn.ConvTranspose2d(     nz, ngf * 8, gfs, 2, gfs//2, bias=False), 
                nn.ReLU(True),
                nn.BatchNorm2d(ngf * 8),
               
                nn.ConvTranspose2d(ngf * 8, ngf * 4, gfs, 2, gfs//2, bias=False),
                nn.ReLU(True),
                nn.BatchNorm2d(ngf * 4),
                
                nn.ConvTranspose2d(ngf * 4, ngf * 2, gfs, 2, gfs//2, bias=False),
                nn.ReLU(True),
                nn.BatchNorm2d(ngf * 2),
                
                nn.ConvTranspose2d(ngf * 2,     ngf, gfs, 2, gfs//2, out_pad,bias=False),
                nn.ReLU(True),
                nn.BatchNorm2d(ngf),
                
                nn.ConvTranspose2d(    ngf,      nc, gfs, 2, 2, bias=False),
                nn.Tanh()
            )
        
        if cuda:
            self.load_state_dict(torch.load(gpath))
        else:
            self.load_state_dict(torch.load(gpath,
                                            map_location=lambda storage,
                                            loc: storage))

    def forward(self, input):
        
        output = self.main(input)
        return output
    

class G_dilat(nn.Module):
    def __init__(self, gpath,nc = 1, nz = 1, ngf = 64, gfs = 5, ngpu = 1,cuda=True):
        super(G_dilat, self).__init__()
        self.ngpu = ngpu

        ngf2=ngf
        self.main = nn.Sequential(

                nn.ConvTranspose2d(     nz, ngf * 8, gfs, 2, gfs//2, bias=False), 
                nn.ReLU(True),
                nn.BatchNorm2d(ngf * 8),
               
                nn.ConvTranspose2d(ngf * 8, ngf * 4, gfs, 2, gfs//2, bias=False),
                nn.ReLU(True),
                nn.BatchNorm2d(ngf * 4),
                
                nn.ConvTranspose2d(ngf * 4, ngf * 2, gfs, 2, gfs//2, bias=False),
                nn.ReLU(True),
                nn.BatchNorm2d(ngf * 2),
                
                nn.ConvTranspose2d(ngf * 2,     ngf, gfs, 2, gfs//2, bias=False),
                nn.ReLU(True),
                nn.BatchNorm2d(ngf),
                
                nn.ConvTranspose2d(    ngf,      nc, gfs, 2, 2, bias=False),

                nn.ReLU(True),
#                
#                ### Start dilations ###
        
                nn.ConvTranspose2d(     nc, ngf2, gfs, 1, 6, output_padding=0,bias=False,dilation=3), #4, dilation 2 / 6, dilation 3
                nn.ReLU(True),
                nn.BatchNorm2d(ngf2),
               
                nn.ConvTranspose2d(    ngf2,  nc, gfs, 1, 10, output_padding=0, bias=False,dilation=5), #6, dilation 3 / 10, dilation 5
                
                nn.Tanh()
                
            )
        if cuda:
            self.load_state_dict(torch.load(gpath))
        else:
            self.load_state_dict(torch.load(gpath,
                                            map_location=lambda storage,
                                            loc: storage))


    def forward(self, input):
        
        output = self.main(input)
        return output
    




