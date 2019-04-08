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



"""
N.B :
- J'ai choisi de créer les fonctions interp1d dans le constructeur.
Ca devrait être plus rapide (quoique) mais ça empeche de passe d'autres méthode d'interp (cubic, etc)
"""


from scipy.interpolate import interp1d
from scipy.constants import c, h, physical_constants
import numpy as np
from pylab import *


class material(object) :


    name = ""

    def __init__(self, model, p, name = "NoName") :

        self.name = name
        self.model = model

        self.wlmin = 400 # by default, for plotting
        self.wlmax = 800

        getattr(self, self.model+'_init')(p)


    def __call__(self, wl) :
        return getattr(self, self.model+'_call')(wl)


    def plot(self, others = []) :

        fig, ax = plt.subplots()
        wl = linspace( self.wlmin, self.wlmax, 100)

        ax.plot(wl, real(self.__call__(wl)), label=self.name+'-n')
        ax.plot(wl, -imag(self.__call__(wl)), label=self.name+'-k')

        for other in others :
              ax.plot(wl, real(other(wl)), label=other.name+'-n')
              ax.plot(wl, -imag(other(wl)), label=other.name+'-k')

        ax.legend(loc=2); # upper left corner
        ax.set_xlabel('$\lambda$')
        ax.set_ylabel('n, k')
        ax.set_title(self.name)


    # Modèle constant

    def constant_init(self, p):
        self.nk = p

    def constant_call(self, wl) :
        return np.ones_like(wl)*self.nk

    # Modèle file

    def file_init(self, p):

        data  = genfromtxt(p, dtype=float) #, delimiter='\t')

        wl = data[:, 0]

        e = physical_constants['elementary charge'][0]

        if max(wl) <= 10 :
            wl = h*c  / e / np.array(wl) # en m
            wl = wl * 10**9 # en nm

        data[:, 0] = wl

        points = sorted(data, key=lambda point: point[0])

        wl, n, k = zip(*points)

        self.wlmin = min(wl)
        self.wlmax = max(wl)
        self.n = interp1d(wl, n)
        self.k = interp1d(wl, k)

    def file_call(self, wl):
        return self.n(wl)- 1.0j*self.k(wl)


    # Modèle Cauchy

    def Cauchy_init(self, p):
        self.ni = p['ni']
        self.ki = p['ki']

    def Cauchy_call(self, wl):
        def Cauchy(wl) :
             l2i = np.power(wl, np.dot(range(len(self.ni)),2))
             n = sum(self.ni/l2i, axis=0)
             k = sum(self.ki/l2i, axis=0)
             return n - 1.0j*k

        return np.vectorize(Cauchy)(wl)
