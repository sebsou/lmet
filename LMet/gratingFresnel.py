#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2019, LTM/CNRS UGA (Laboratoire des Technologies de la 
# microélectronique, Centre National de la Recherche Scientifique, 
# Université Grenoble Alpes) 

# sebastien.soulan_at_univ-grenoble-alpes.fr

# This software is a computer program whose purpose is to develop 
# advanced metrology techniques.

# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 

# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 

# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 

# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

# Credit to : Mohammed EL KODADI

import numpy as np
import pandas as pd
from math import sin, cos, sqrt, exp, pi

from LMet.grating import *
from LMet.SVGScene import *

class gratingFresnel(grating) :
    
    
    def build(self, params = None):
        
        if params == None :
            params = self.parameters
        
        if hasattr(self, 'instructions') :
            
            # Deploie les paramètres           
            for p in params :
                #print("%s = params['%s']" % (p,p))
                exec("global %s" % p)
                exec("%s = params['%s']" % (p,p))
                
            # Déploie les fonctions de constructions            
            exec(self.buildinstructions)
             
            # Construit.            
            exec(self.instructions)   

        # Finalise          
        self.ec = self.parameters['ec']
        self.mat = self.parameters['mat']
     
     
   
     
    def rprs(self, e, mat, context) :
        
        def singlerprs(e, mat, PHI0, LAMBDA) :
                    
            N0 = mat[0](LAMBDA) # air
            N1 = mat[1](LAMBDA) # res
            N2 = mat[2](LAMBDA) # sub
            
            # limite milieu 0-1 air reSINe
            
            SIN_PHI1 = (N0*sin(PHI0)) / N1
            COS_PHI1 = sqrt(1-(SIN_PHI1)**2)
  
            R01P=( N1 * cos( PHI0 ) - N0* COS_PHI1 ) / ( N1 * cos( PHI0 ) + N0* COS_PHI1 ) 
            R01S=( N0 * cos( PHI0 ) - N1* COS_PHI1 ) / ( N0 * cos( PHI0 ) + N1* COS_PHI1 )
            
            # limite milieu 1-2 resine substrat#
            
            SIN_PHI2 = ( N1 * SIN_PHI1) / N2
            COS_PHI2 = sqrt( 1 - (SIN_PHI2)**2 )
  
            R12P=( N2* COS_PHI1 - N1* COS_PHI2 ) / ( N2* COS_PHI1 + N1* COS_PHI2 ) 
            R12S=( N1* COS_PHI1 - N2* COS_PHI2 ) / ( N1* COS_PHI1 + N2* COS_PHI2 ) 
            
            #
            
            BETA = (( 4 * pi) / LAMBDA)* e * N1* COS_PHI1
            
            RP = (    R01P + R12P * exp( -1.0j * BETA) ) /  (1 + R01P * R12P * exp( -1.0j * BETA) )
            RS = (    R01S + R12S * exp( -1.0j * BETA) ) /  (1 + R01S * R12S * exp( -1.0j * BETA) ) 

            return  RP / RS
       
        
        
        d = []
        
        for idx, row in context.iterrows():
            LAMBDA = row['Hv']
            PHI0 = row['angle']
            rprs = singlerprs(e, mat, PHI0, LAMBDA) 
            d.append ([PHI0, LAMBDA, np.real(rprs), np.imag(rprs)] )
         
                  
        d = np.array(d)
        return pd.DataFrame(d, columns = ["angle","Hv", "rp", "rs"])
            
    def compute(self, context) :
        
        ctx = context
        
        if isinstance(context, dict) : # traditionnel
            
            import itertools
            ctx = pd.DataFrame(list(itertools.product(context['angle'], context['Hv'])))
            ctx.columns = ['angle', 'Hv']

        rp_rs = self.rprs(self.ec, self.mat, ctx)
                     
        return rp_rs
        
    def display(self):
 
        print ("ec : ", self.ec)
      
    def showSVG(self) :
        # Vu sur http://nbviewer.ipython.org/gist/rpmuller/5666810       
   
        sw = 400
        y0 = 50
        e0 = 100 # epaisseur ciel et substrat
        
        #http://www.december.com/html/spec/colorsvg.html
        colors = ( "aliceblue", "lightgreen", "lightsteelblue", "lightsalmon",
                  "lightgrey", "limegreen", "plum", "peachpuff")
        
        scene = SVGScene()
            
        newnl = np.hstack(self.mat)
        matdict = {}
        i=0
        for m in newnl :
            if m  not in matdict :
                matdict[m]= colors[i]
                i +=1

        y=y0
        for m in matdict :
            scene.rectangle((sw+25,y), 50,  50, matdict[m]) #(0,255,255)) 
            scene.text((sw+100,y+30), m.name, size=16) 
            y += 75
        
   
        y=y0
        
        ec2 = [0, self.ec, 0]
        
        for i, e in enumerate(reversed(ec2)) :
         
         #print "couche ", i, e
    
            if e == 0 : e=e0
        
            f = self.mat[i]
            scene.rectangle((0,y), e,  sw, matdict[f]) #(0,255,255)) 
            
            y += e
        
        return scene
    
    

        

        