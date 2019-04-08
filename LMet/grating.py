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


import numpy as np
from LMet.material import *
from LMet.signature import *
from LMet.mmfe import *

import copy
from scipy.optimize import minimize


class grating(object) :

    buildinstructions =   """ 
    """

    def ip(self, sig, ptooptimize, *args, **kwargs):

        g = copy.deepcopy(self)
       
        g.build() # je suppose que instr est déjà entré
        
        def f(pvect) :
            for i, p in enumerate(ptooptimize) :
                g.parameters[p] = pvect[i] 
               
            g.build()

            tempsig = g.signature(sig.context)
            
            return sig.distance(tempsig, ["Is", "Ic"]) 
            
            # Remarque : le choix de Is Ic est forcé ici
            # car il est plus long et plus compliqué de faire
            # des distance de fonction modulo 2pi.
        
        res = minimize(f, *args, **kwargs)
     
        return res
    
                
    def signature(self, context = None, ellipso = "JobinYvon") :

        if context == None : 
            context = self.context
        
            
        rp_rs = self.compute(context)
                     
        sig = signature(rp_rs, ellipso)
        sig.context = context
        sig.name = "Computed from "+ self.__class__.__name__ + "."
        sig.isic
    
        return sig
        
    
    def library(self, context = None) :
        
        if context == None : 
            context = self.context
        
        l = []
        
        for p in self.libparameters :
            self.build(p)
            s = np.hstack(self.compute(context)[['rp', 'rs']].values)
            l.append(s)
        
        return np.array(l)
    
    # uses dataframe
    # trouver un moyen de melanger les deux
      
    
    def librarydf(self, df, context = None, progress_bar = lambda x : None) :
        
        if context == None : 
            context = self.context
        
        l = []
            
        p = copy.deepcopy(self.parameters)
         
        for i in range(len(df)) :
            for k in df.ix[i].keys():
                p[k] = df.ix[i][k]
            self.build(p)
            s = np.hstack(self.compute(context)[['rp', 'rs']].values)
            l.append(s)
            progress_bar(i+1)  
        
        return np.array(l)

        
