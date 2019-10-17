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
from LMet.grating import *

import copy
from scipy.optimize import minimize
from LMet.SVGScene import *

class gratingMMFE1D(grating) :

    #with open ("buildinstructions.py", "r") as f:
    #    buildinstructions1D = f.read()

    buildinstructions =   """

def initsquare(dd, cd, h, sub, res, air, nlayer=1) :
    global d, ec, eta, nulam

    n = nlayer
    d = dd

    ec = np.zeros(n+2)
    ec[1:n+1] = h/n

    eta = np.zeros( (2, n+2))
    eta[0,1:n+1] = 0.5-cd/d/2
    eta[1,1:n+1] = 0.5+cd/d/2

    nulam = np.zeros((3, n+2))
    nulam = np.array([[air]*(n+2)]*3)
    nulam[1, 1:n+1] = res
    nulam[:, n+1] = sub





def addlayer(i, e, cd, air, mat) :
    global ec, eta, nulam

    ec = np.insert(ec, i, e)
    eta = np.insert(eta, i, np.array([0.5-cd/d/2, 0.5+cd/d/2]), 1) # d'ou la necessité de bien faire les eta, nulam...
    nulam = np.insert(nulam, i, np.array([air, mat, air]),1)

    """

    def build(self, params = None):

        if params == None :
            params = self.parameters

        #instr = self.instructions

        global d, ec, eta, nulam
        d = None
        ec = None
        eta = None
        nulam = None

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

        self.ec = ec
        self.eta = eta
        self.nulam = nulam
        self.d = d

    def display(self):

        print ("ec : ", self.ec)
        print ("eta : ", self.eta)

        def showname(mat):
            return mat.name
        shownamevec = np.vectorize(showname)
        print ("nulam : ", shownamevec(self.nulam))

    def compute(self, context):

        ctx = context

        # if isinstance(context, dict): # traditionnel

        #     import itertools
        #     ctx = pd.DataFrame(list(itertools.product(context['angle'], context['𝜆'], [context['M']])))
        #     ctx.columns = ['angle', '𝜆', 'M']


        rp_rs = rprs(self.ec,
                     self.eta,
                     self.nulam,
                     self.d,
                     ctx)

        return rp_rs

    # TODO :
    """
    def _repr_svg_(self) :
        self.showSVG()
    """

    def showSVG(self):
        # Vu sur http://nbviewer.ipython.org/gist/rpmuller/5666810

        sw = 400
        y0 = 50
        e0 = 100  # epaisseur ciel et substrat

        #http://www.december.com/html/spec/colorsvg.html
        colors = ("aliceblue", "lightgreen", "lightsteelblue", "lightsalmon",
                  "lightgrey", "limegreen", "plum", "peachpuff")

        scene = SVGScene()

        newnl = np.hstack(self.nulam)
        matdict = {}
        i = 0
        for m in newnl:
            if m not in matdict:
                matdict[m] = colors[i]
                i += 1

        y = y0
        for m in matdict:
            scene.rectangle((sw+25, y), 50,  50, matdict[m])
            scene.text((sw+100, y+30), m.name, size=16)
            y += 75

        y = y0
        for i in range(len(self.ec)):

            #for i, e in enumerate(reversed(self.ec)):
            #print "couche ", i, e
            e = self.ec[i]

            if e == 0:
                e = e0

            eta = self.eta[:, i]
            #print "eta = ", eta

            if not any(eta):
                f = self.nulam[:, i][0]
                scene.rectangle((0, y), e,  sw, matdict[f])
            else:
                eta2 = eta[:]
                eta2 = np.insert(eta2, 0, 0)
                eta2 = np.append(eta2, 1)

                x = 0
                for j, f in enumerate(self.nulam[:, i]):
                    w = sw*(eta2[j+1]-eta2[j])
                    scene.rectangle((x, y), e, w, matdict[f])
                    #print "rectangle :", (x,y), w, e
                    x += w

            y += e

        scene.bbox()

        return scene
