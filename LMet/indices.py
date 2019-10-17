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


import itertools
import numpy as np
import pandas as pd

from sklearn.neighbors import NearestNeighbors
from multiprocessing import Manager, Process

"""
def nklib(g, wl, mat, ns, ks) :

    lib = []
    for n, k in itertools.product(ns, ks):

        mat.nk = n + 1.j*k
        g.context["𝜆"] = np.array([wl])
        lib.append(g.signature().data[['Is', 'Ic']].values[0])

    return np.array(lib)
"""
def nklib_keys(ns, ks) :
    return list(itertools.product(ns, ks))


"""

def nklibset(g, wls, mat, ns, ks, progress_bar = lambda x : None) :

    def worker(procnum, wl, return_dict):
        '''worker function'''
        #print str(procnum) + ' represent!'
        libset = nklib(g, wl, mat, ns, ks)
        return_dict[procnum] = libset



    manager = Manager()
    return_dict = manager.dict()

    jobs = []
    for i, wl in enumerate(wls) :
        p = Process(target=worker, args=(i, wl, return_dict))
        jobs.append(p)
        p.start()


    for i, proc in enumerate(jobs):
        proc.join()
        progress_bar(i+1)

    return return_dict.values()
"""

# TODO : a fusionner avec celui précedent.
# je pense qu'elle est générique presque.
def nklib_newctx(g, wl, angles, mat, ns, ks) :

    lib = []

    for n, k in itertools.product(ns, ks):

        ct = {"𝜆" : [wl],
              'angle' : angles}

        mat.__init__("constant",  n + 1.j*k)

        lib.append(g.signature(ct).data[['Is', 'Ic']].values)

    n_nk = len(ns)*len(ks)
    sigsize = 2*len(angles)

    return np.array(lib).reshape(n_nk, sigsize)


"""
def nklibset(g, wls, mat, ns, ks, progress_bar = lambda x : None) :

    libset = [] # ensemble de nklib indéxé par wl

    for i, wl in enumerate(wls) :
        libset.append(nklib(g, wl, mat, ns, ks))
        progress_bar(i+1)

    return libset
"""
# ici, on considère que le context ne fait pas partie de g.
# on parallelise wl par wl; donc on envoie une serie d'angle.

def nklibset_newctx(g, ctx, mat, ns, ks, progress_bar = lambda x : None) :

    groups = ctx.groupby(['𝜆']).groups
    wls = np.sort(groups.keys())

    libset = [] # ensemble de nklib indéxé par wl

    for i, wl in enumerate(wls):
        libset.append(nklib_newctx(g, wl, groups[wl], mat, ns, ks))
        progress_bar(i+1)

    return libset


    """
    def worker(procnum, row, return_dict):
        '''worker function'''
        #print str(procnum) + ' represent!'
        libset = nklib_newctx(g, row, mat, ns, ks)
        return_dict[procnum] = libset


    manager = Manager()
    return_dict = manager.dict()

    jobs = []
    for i, row in ctx.iterrows() :
        p = Process(target=worker, args=(i, row, return_dict))
        jobs.append(p)
        p.start()


    for i, proc in enumerate(jobs):
        proc.join()
        progress_bar(i+1)

    return return_dict.values()
    """

# exemple : cut (libset, 'n', [1.5, 1.6], [500, 600])
def cut(wls, keys, libset, n_or_k, nkrange, wlrange) :


    if n_or_k == 'n' : inork = 0
    if n_or_k == 'k' : inork = 1

    BIG = -10

    for iwl, wl in enumerate(wls) :
        if wl >= wlrange[0] and wl <= wlrange[1] :

            for i, key in enumerate(keys) :
                if key[inork] >= nkrange[0] and key[inork] <= nkrange[1]:
                    libset[iwl][i] = [BIG, BIG]

def determine_library(sig, wls, keys, libset, progress_bar = lambda x : None, nn = 2) :

    pts = []

    for i, wl in enumerate(np.sort(wls)) :


        nklib = libset[i] #(g, wl, mat, ns, ks)
        nbrs = NearestNeighbors(n_neighbors=nn, algorithm='ball_tree').fit(nklib)

        refsig = np.hstack(sig.data.loc[sig.data['𝜆'] == wl][['Is', 'Ic']].values)

        distances, indices = nbrs.kneighbors(refsig) #sig.data[['Is', 'Ic']].values[i])

        for inn in indices[0] :
            pts.append( [wl, keys[inn][0], keys[inn][1]] )

        progress_bar(i+1)

    return  pd.DataFrame(pts, columns = ['𝜆', 'n', 'k'])