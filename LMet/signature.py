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



import os
import numpy as np
from numpy import cos, arccos,  sin, arcsin, arctan, mod, pi, angle
from scipy.constants import c, h, physical_constants
import math
import pandas as pd
import linecache


class signature():

    ellipso = 'JobinYvon'  # default
    modulator = math.radians(0)
    analyser = math.radians(45)

    def __init__(self, df=None, ellipso="JobinYvon"):

        if type(df) is str:
            self.openfile(df)
        else:
            self.data = df
            self.ellipso = ellipso

        """
    def __call__(self, 𝜆) :
        if  not hasattr(self, "f1") : self.f1 = interp1d(self.data[:, 0], self.data[:, 1])
        if  not hasattr(self, "f2") : self.f2 = interp1d(self.data[:, 0], self.data[:, 2])
        return  self.f1(𝜆), self.f2(𝜆)
"""
    def distance(self, sig, columns, *args, **kwargs):

        m = self.data[columns].values - sig.data[columns].values
        return np.linalg.norm(m)  # , args, kwargs)

    def openfile(self, fname):

        if os.path.splitext(fname)[1] == '.spe':
            self.openJobinYvon(fname)
        elif os.path.splitext(fname)[1] == '.dat':
            self.openWoolam(fname)

    def openJobinYvon(self, f):

        n = 0  
        datablock = 0
      

        with open(f,  encoding = "ISO-8859-1") as search:
      
            for line in search:
                line = line.rstrip()  # remove '\n' at end of line
                if line == "## ELLIPSOMETRIC CONFIGURATION:":
                    self.modulator = float(next(search).split(' ')[0]) ;  n += 1
                    self.analyser = float(next(search).split(' ')[0]) ;  n += 1
                if line == "# INCIDENCE ANGLE:" :
                    self.incidence = float(next(search).split(' ')[0]) ;  n += 1
                if line == '# DATA:': 
                    datablock = n + 1
                    break
                if line == '# NUMBER OF POINTS:' : 
                    number_of_points = int(next(search).split(' ')[0]) ;  n += 1
                
                n += 1
        
        
        self.data = pd.read_csv(f, header=datablock,
                                delim_whitespace=True,
                                skip_blank_lines=False,
                                encoding = "ISO-8859-1",
                                nrows = number_of_points)
        
        self.ellipso = 'JobinYvon'


        """
        with open(fname, encoding = "ISO-8859-1") as search:


            for i, line in enumerate(search):
                line = line.decode(errors='ignore').rstrip()  # remove '\n' at end of line
                if line == "## ELLIPSOMETRIC CONFIGURATION:":
                    elconf = i
                #if line == "# INCIDENCE ANGLE:" :
                #    incang = i
                if line == "# NUMBER OF POINTS:" :
                    inpoints = i
                if line == "# DATA:":
                    toskip = i
                    break


        #print (linecache.getline(fname, elconf+2).decode(errors='ignore').split(' '))
        #print (linecache.getline(fname, elconf+3)) #.split(' '))
        self.modulator = float(linecache.getline(fname, elconf+2).split(' ')[0])
        self.analyser = float(linecache.getline(fname, elconf+3).split(' ')[0])
        npoints = int(linecache.getline(fname, inpoints+2).split(' ')[0])

        self.data = pd.read_csv(fname, header=toskip+1,
                                delim_whitespace=True,
                                skip_blank_lines=False,
                                encoding = "ISO-8859-1",
                                nrows = npoints)
        self.ellipso = 'JobinYvon'
        """

    def openNanometrix(self, fname): # Nanometrix Atlas XP+
                # TODO : analyser

                self.data = pd.read_csv(fname,header=None, skiprows=2,
                                        delim_whitespace=True,
                                        skip_blank_lines=False)

                self.data.columns = ['𝜆', 'Is', 'Ic']

                self.ellipso = 'Nivea'


    def openWoolam(self, fname):

        # TODO : analyser


        with open(fname) as search:
            for i, line in enumerate(search):
                line = line.rstrip()  # remove '\n' at end of line
                if line == "nm":
                    toskip = i

        self.data = pd.read_csv(fname, header=toskip+1,
                                delim_whitespace=True,
                                skip_blank_lines=False)

        self.data.columns = ['𝜆', 'angle', 'Psi',
                             'Delta', 'ErrPsi', 'ErrDelta']

        self.data[['angle', 'Psi', 'Delta', 'ErrPsi', 'ErrDelta']] = np.radians(
            self.data[['angle', 'Psi', 'Delta', 'ErrPsi', 'ErrDelta']])

        #self.data['angle'] = np.radians(self.data['angle'])
        #self.data['angle'] = np.radians(self.data['angle'])

        self.ellipso = 'Woolam'

    @property
    def isic(self):
        if ('Is' not in self.data.columns) or ('Ic' not in self.data.columns):
            self.computesignals()

        return self.data[['𝜆', 'Is', 'Ic']]

    def eV(self):
        e = physical_constants['elementary charge'][0]
        self.data['eV'] = h * c / e / self.data['𝜆'] / 10**(-9)

    def 𝜆(self):
        e = physical_constants['elementary charge'][0]
        self.data['𝜆'] = h * c / e / self.data['eV'] / 10**(-9)


    def cmmo(self):
        self.data['cmmo'] = 1 / (self.data['𝜆'] * 10**(-7))

    def compute_isic(self):

        # TODO : verifier que polarize et analiser est pareil
        P_M = self.analyser - self.modulator  # polariseur-modulateur

        psi = self.data['Psi']
        delta = self.data['Delta']

        Is = sin(2*P_M)*sin(2*self.analyser)*sin(2*psi)*sin(delta)
        # normalisé, (cf. Techniques de l'ingénieur, § R6 490 2.3 Ellipsométrie à modulation de phase)

        Ic = sin(2*P_M)*( (cos(2*psi)-cos(2*self.analyser))*sin(2*self.modulator) +
                         sin(2*self.analyser)*cos(2*self.modulator)*sin(2*psi)*cos(delta) )

        self.data['Is'] = Is
        self.data['Ic'] = Ic

    def computesignals(self):

        rprs = self.data['rp']+1.0j*self.data['rs']
        tan_psi, cos_del = abs(rprs), cos(angle(rprs))
        psi, delta = arctan(tan_psi), angle(rprs)  # en RADIAN !

        #  ICI:
        #          0 < tan_psi < inf            0 < psi < pi/2
        #         -1 < cos_del < 1            -pi < delta < pi
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #
        #     NOTE IMPORTANTE :
        #
        #  les donnees ellipso sont fournies de maniere variable :
        #
        # - MMFE :       0 < psi < pi/2  et    -pi < delta < pi
        # - KLA 1280SE : 0 < psi < pi/2  et      0 < delta < pi   ,
        #                car psi obtenu à partir de tan(psi) et
        #                delta obtenu à partir de cos(delta)
        #                [typiquement une determination de psi et delta
        #                par une approche mesure tan(psi)-cos(delta)]
        # - Jobin-Yvon : 0 < psi < pi/4  et      0 < delta < 2*pi ,
        #                car psi obtenu à partir de sin(2*psi)  et delta
        #                obtenu à partir d'une combinaison
        #                cos(delta) & sin(delta)
        #                [typiquement une determination de psi et delta par
        #                une approche mesure Is-Ic]
        # - Woolam :     0 < psi < pi/2  et   0 < delta < 2*pi
        #
        # Il en resulte le switch ci-dessous:

        if self.ellipso == 'JobinYvon':
            psi = arcsin(sin(2*psi))/2.
            # car mesures experimentales psi determinees à partir de sin(2*psi)
            # avec Jobin-Yvon, avec 0 < psi_simul < pi/4

            delta = mod(delta, 2*pi)
            # en radian, avec transfo en :  0 < delta_simul < 2*pi

        elif self.ellipso == 'Woolam':
            # psi, comme standard
            delta = mod(delta, 2*pi)

        else:  # ellipso=='KLA 1280SE' ou autre
            delta = arccos(cos_del)
            # transfo en radian pour Matlab, avec 0 < delta_simul < pi

        if self.modulator == 0:
            S1, S2 = cos(2*psi), sin(2*psi)*cos(delta)
        #else:
            #S1, S2 = NaN*rprs, NaN*rprs  # [NaN*rprs]*2 est un alias de liste dangereux !
            #alpha, beta = (tan_psi**2-1)/(tan_psi**2+1), (2*tan_psi*cos_del)/(tan_psi**2+1)

        self.data['Psi'] = psi
        self.data['Delta'] = delta

        self.compute_isic()  # nécessite psi-delta

        self.data['S1'] = S1
        self.data['S2'] = S2

        #self.data['alpha'] = alpha
        #self.data['beta'] = beta





"""
    @Isic.setter
    def isic(self, value) : self._isic = value

    @isic.deleter
    def isic(self): del self._isic

    def tanpsicosdelta(self, 𝜆) :

        l=self._rprs[0]
        tanpsi = abs(self._rprs[1])
        cosdel = cos(angle(self._rprs[1]))

        f1 = interp1d(l, _isic[:, 1])
        f2 = interp1d(l, _isic[:, 2])

        return ( f1(𝜆), f2(𝜆) )

    @property
    def tanpsicosdelta(self):
        tan_psi = abs(self.rprs[1])
        cos_del = cos(angle(self.rprs[1]))
        return [ self._rprs[0], tanpsi, cosdelta ]

    @isic.setter
    def psidelta(self, value) : self._isic = value

    @isic.deleter
    def psidelta(self): del self._isic
     """
