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


# Credits to :
# David FUARD 
# Jean-Hervé TORTAI

"""mmfe.py: 1D Modal method by Fourier expansion source code."""

import numpy as np
import pandas as pd


# reorganiser ce bordel ...
#from math import *
from operator import itemgetter
from random import random, randrange

from numpy import pi, NaN, Inf,imag, real, conj, angle, floor, ceil, \
     sin,cos,tan,arcsin,arccos,arctan, deg2rad,rad2deg, sqrt,exp,log,log10,mod,\
     array,ndarray,matrix, arange,linspace,hstack,vstack, zeros,ones,diag,diagflat,eye,meshgrid, diff,argsort,prod,dot, nansum,nan_to_num,isnan,\
     transpose,tril,triu, polyfit,\
     set_printoptions

from numpy import max as npmax
from numpy import min as npmin
from numpy import mean as npmean
from numpy import sum as npsum
from numpy.linalg import eig, inv,lstsq
from scipy.linalg import toeplitz

# tri 1D
def tri1D(cos_carre_teta):
    """

    """

    epsi = 1e-06
    rp = sqrt(cos_carre_teta)
    # index =( rp.imag>epsi or (rp.real<0 and abs(rp.imag)<epsi) ) ;   rp(index) =-rp(index)
    
    A = rp.imag > epsi
    B = (rp.real < 0) & (abs(rp.imag) < epsi)


    #rp[A+B-A*B] *=-1 # rp*(1-2*(A+B-A*B)) # mais c'est XOR en fait!!!
    rp[A^B] *= -1 
    
    return rp
    # print 'Code testÈ... rÈsultat ok'

    

# cren 1D
def cren1D(etaC,nulamC,M):
    " fonction crÈneau "
    eta1, eta2 =vstack([0,etaC.reshape(-1,1)]), vstack([etaC.reshape(-1,1),1]) # remise sous forme colonne : .reshape(-1,1) ; remise sous forme ligne : .reshape(1,-1)
    epr,  m_   =nulamC.reshape(-1,1)**2,  -2*pi*(arange(2*M)+1)

    #print (eta1, eta2, epr)
    epr0 = npsum( (eta2-eta1) * epr )
    eprm = dot(epr,ones((1,2*M)))
    m =  m_.reshape(1,-1)

    e1m,  e2m  =1j*dot(eta1,m),  1j*dot(eta2,m)
    epr_m =hstack([0., -1j* npsum( (exp(e2m) -exp(e1m)) *eprm ,0) /m_ ])
    epr_p =hstack([0.,  1j* npsum( (exp(-e2m)-exp(-e1m))*eprm ,0) /m_ ])
    return diagflat(epr0*ones(2*M+1)) + triu( toeplitz(epr_m.real)+1j*toeplitz(epr_m.imag) ) + tril( toeplitz(epr_p.real)+1j*toeplitz(epr_p.imag) )
    # print 'Code testÈ... rÈsultat ok'


# <codecell>

# casc 1D
def casc1D(S1,S2,nmax):
    " cascadage des matrices S entre les couches "
    #
    GA11 =S1[:nmax,:nmax]
    GA12 =S1[:nmax,nmax:]
    GA21 =S1[nmax:,:nmax]
    GA22 =S1[nmax:,nmax:]
    #
    GB11 =S2[:nmax,:nmax]
    GB12 =S2[:nmax,nmax:]
    GB21 =S2[nmax:,:nmax]
    GB22 =S2[nmax:,nmax:]
    #
    ID =eye(nmax)
    U1 =ID-dot(GA22,GB11)
    U2 =ID-dot(GB11,GA22)
    #
    S11 =GA11+ dot(GA12, dot( dot(inv(U2),GB11) ,GA21) )
    S12 =dot(GA12, dot(inv(U2),GB12) )
    S21 =dot(GB21, dot(inv(U1),GA21) )
    S22 =GB22+ dot(GB21, dot( dot(inv(U1),GA22) ,GB12) )
    #
    return  vstack( (hstack((S11,S12)),hstack((S21,S22))) )
    # print 'Code testÈ... rÈsultat ok'

# <codecell>

# rp & rs 1D
def S_rp_1D(M,nmax,Ncouches, ec,eta,nulam, nu1,alpham,k):
    " mat S (te & tm modes) computation "
    epsi =1e-06
    # -
    S =zeros((2*nmax,2*nmax),dtype='complex')
    S_te =array(S)
    S_te[:nmax,nmax:] =eye(nmax,dtype='complex')
    S_te[nmax:,:nmax] =eye(nmax,dtype='complex')
    S_tm =array(S_te)
    # -
    Fp_te_homogene =eye(nmax,dtype='complex')
    Fp_tm_homogene =eye(nmax,dtype='complex')
    #
    rp_te_air =tri1D(nu1**2-alpham**2)+0j # superstrat forcement homogËne,  cos_carre_teta =nu1**2-alpham*alpham
    rp_te =zeros((rp_te_air.size,Ncouches),dtype='complex')
    rp_te[:,0:1] =rp_te_air
    rp_tm_air =array(rp_te_air,dtype='complex')
    rp_tm =array(rp_te)
    #
    Gp_te_air =diagflat(rp_te_air)
    Gp_tm_air =array(Gp_te_air)
    rp_te1 =array(rp_te_air)
    Fp_te1 =array(Fp_te_homogene)
    Gp_te1 =diagflat(rp_te_air)
    rp_tm1 =array(rp_te1)
    Fp_tm1 =array(Fp_te1)
    Gp_tm1 =diagflat(rp_te_air/(nu1**2))
    #
    for j in range(ec.size-1):
        phip_ligne_te =dot( ones((rp_te1.size,1)) , exp(-1j*k*ec[j]*rp_te1.reshape(1,-1)) )
        phip_col_te =dot( exp(-1j*k*ec[j]*rp_te1.reshape(-1,1)) , ones((1,rp_te1.size)) )
        S_te[:nmax,nmax:] =S_te[:nmax,nmax:]*phip_ligne_te
        S_te[nmax:,:nmax] =S_te[nmax:,:nmax]*phip_col_te
        S_te[nmax:,nmax:] =S_te[nmax:,nmax:]*phip_ligne_te*phip_col_te
        #
        phip_ligne_tm =dot( ones((rp_tm1.size,1)) , exp(-1j*k*ec[j]*rp_tm1.reshape(1,-1)) )
        phip_col_tm =dot( exp(-1j*k*ec[j]*rp_tm1.reshape(-1,1)) , ones((1,rp_tm1.size)) )
        S_tm[:nmax,nmax:] =S_tm[:nmax,nmax:]*phip_ligne_tm
        S_tm[nmax:,:nmax] =S_tm[nmax:,:nmax]*phip_col_tm
        S_tm[nmax:,nmax:] =S_tm[nmax:,nmax:]*phip_ligne_tm*phip_col_tm
        #
        if npmax(nulam[:,j+1])-npmin(nulam[:,j+1])>1e-8: # si couche hÈtÈrogËne
            # ---
            epr =cren1D(eta[:,j+1],nulam[:,j+1],M)

            # "awad_te", calcul de valeurs propres, mÈthode Moharam et Gaylord, formulation de D.M. Pai et K.A. Awada [Vol.8, No. 5/May 1991/J. Opt. Soc.Am. A pp755 761.]
            rp_te2,Fp_te2 =eig(epr-diagflat(alpham**2))
            rp_te2 =sqrt(rp_te2)
            #A,B =rp_te2.imag>epsi, (rp_te2.real<0)*(abs(rp_te2.imag)<epsi)
            #rp_te2[A+B-A*B] *=-1
            
            A = rp_te2.imag>epsi
            B = (rp_te2.real<0) & (abs(rp_te2.imag)<epsi)
            rp_te2[A^B] *= -1
            


            rp_te[:,j+1:j+2] =rp_te2.reshape(-1,1)
            Gp_te2 =Fp_te2* dot( ones((rp_te2.size,1)) , rp_te2.reshape(1,-1) )
            # ---
            eprm_inv =cren1D( eta[:,j+1], 1./nulam[:,j+1], M)



            # "awad_tm"
            C =dot( inv(eprm_inv) , eye(nmax)-dot(diagflat(alpham),dot(inv(epr),diagflat(alpham))) )
            rp_tm2,Fp_tm2 =eig(C)
            rp_tm2 =sqrt(rp_tm2)
            A =abs(rp_tm2.imag)<epsi
            rp_tm2[A] =rp_tm2.real[A]
            #A,B =rp_tm2.imag>epsi, (rp_tm2.real<0)*(abs(rp_tm2.imag)<epsi)
            #rp_tm2[A+B-A*B] *=-1
            
            A = rp_tm2.imag>epsi
            B = (rp_tm2.real<0) & (abs(rp_tm2.imag)<epsi)
            rp_tm2[A^B] *= -1
            



            rp_tm[:,j+1:j+2] =rp_tm2.reshape(-1,1)
            Gp_tm2 = dot(eprm_inv, Fp_tm2)* dot(ones((rp_tm2.size, 1)), rp_tm2.reshape(1, -1))
            # ---
        else:
            rp_te2 = tri1D(nulam[0, j+1]**2-alpham**2)
            Fp_te2 = array(Fp_te_homogene)
            rp_te[:,j+1:j+2] =rp_te2.reshape(-1,1)
            Gp_te2 = diagflat(rp_te2[:nmax])
            # ---
            rp_tm2 = array(rp_te2)
            Fp_tm2 = array(Fp_tm_homogene)
            rp_tm[:,j+1:j+2] =rp_tm2.reshape(-1,1)
            Gp_tm2 = diagflat(rp_tm2[:nmax]/(nulam[0,j+1]**2))
            # ---
        # ---
        SA_te =vstack( (hstack((Fp_te1,-Fp_te2)) , hstack((Gp_te1,Gp_te2))) )
        SB_te =vstack( (hstack((-Fp_te1,Fp_te2)) , hstack((Gp_te1,Gp_te2))) )
        S_te2 =dot(inv(SA_te),SB_te) # interface (j/j+1)     # S_te2=SA_te\SB_te, o˘ S_te2 est solution de dot(SA_te,S_te2)=SB_te, soit S_te2 =dot( inv(SA_te) , SB_te )
        S_te =casc1D(S_te,S_te2,nmax)
        #
        rp_te1 =array(rp_te2)
        Fp_te1 =array(Fp_te2)
        Gp_te1 =array(Gp_te2)
        # ---
        SA_tm =vstack( (hstack((Fp_tm1,-Fp_tm2)) , hstack((Gp_tm1,Gp_tm2))) )
        SB_tm =vstack( (hstack((-Fp_tm1,Fp_tm2)) , hstack((Gp_tm1,Gp_tm2))) )
        S_tm2 =dot(inv(SA_tm),SB_tm) # interface (j/j+1)     # X=A\B, o˘ X est solution de dot(A,X)=B, soit X =dot(inv(A),B)
        S_tm =casc1D(S_tm,S_tm2,nmax)
        #
        rp_tm1 =array(rp_tm2)
        Fp_tm1 =array(Fp_tm2)
        Gp_tm1 =array(Gp_tm2)
    return ( S_te, rp_te, S_tm, rp_tm )
    # print 'Code testÈ... rÈsultat ok'

# <codecell>

# efficacitÈ 1D
def efficacite_ref_1D(rp, S, M, nmax) :
    " calcul des efficacitÈs de diffraction "
    #
    # ****************  efficacites de diffraction  ************************
    #  etx1 : coeff de reflexion    sur les differents ordres non evanescents
    #  etxS : coeff de transmission sur les differents ordres non evanescents
    #
    m =arange(-M,M+1)
    #
    # directions de diffraction (constantes de propagation rÈelles)
    or1 =abs(rp[:,0].imag)<1e-25
    orS =abs(rp[:,-1].imag)<1e-25
    I1 =m*or1
##    IS =m*orS
##    ordr1 =arange(npmin(I1),npmax(I1)+1)           PARTIE DE CODE ERRONEE POUR LE CALCUL DE L'EFFICACITE (cf. pb avec prÈsence de modes Èvanescents)
##    ordrS =arange(npmin(IS),npmax(IS)+1)           PARTIE DE CODE ERRONEE POUR LE CALCUL DE L'EFFICACITE (cf. pb avec prÈsence de modes Èvanescents)
##    nm1 =M+npmin(I1)                               PARTIE DE CODE ERRONEE POUR LE CALCUL DE L'EFFICACITE
##    nmS =M+npmin(IS)                               PARTIE DE CODE ERRONEE POUR LE CALCUL DE L'EFFICACITE
    ordr1new,ordrSnew =m[or1],m[orS]
    if ordr1new.size==0:  ordr1new =arange(0,1)
    if ordrSnew.size==0:  ordrSnew =arange(0,1)
    #
    costeta1 =rp[or1,0]
    costeta0 =costeta1[abs(npmin(I1))]
    costetaS =rp[orS,-1]
    #
    # Efficacites
##    S11 =S[nm1:nm1+ordr1.size,M]                   PARTIE DE CODE ERRONEE POUR LE CALCUL DE L'EFFICACITE
##    S21 =S[nmax+nmS:nmax+nmS+ordrS.size,M]         PARTIE DE CODE ERRONEE POUR LE CALCUL DE L'EFFICACITE
    S11 =S[M+ordr1new,M]
    S21 =S[nmax+M+ordrSnew,M]
    #
    e1 =(costeta1 * S11 * conj(S11) /costeta0 ).real
    eS =(costetaS * S21 * conj(S21) /costeta0 ).real
    #
##    return ( e1 , eS , ordr1 , ordrS , costeta1 , costetaS , S11 )
    return ( e1 , eS , ordr1new , ordrSnew , costeta1 , costetaS , S11 )
    # print 'Code testÈ... rÈsultat ok'

# <codecell>



def rprs(ec, eta, nl, pitch, context) : #angle, wl, M) :
    """
    angle : de taille quelconque
    wl : idem
    nl : tableau (iwl) de tableau de complexes (3d)
         ou tableau (2D) de materials
    """


    def mat2cx(nl, wl):
        def ope (f) : return f(wl)
        ope = np.vectorize(ope)
        return ope(nl)


    """
    if ndim == 2 : # utilisation de nulam_mat
        nulam = []
        for iL,lambda0 in enumerate(wl):
            nulam.append( mat2cx(nl, lambda0))
    else :
        nulam = nl
    """


    #pitch =float(G.lay[0].param['pitch'].currentVal)

    #

    #d = pd.DataFrame()

    def singlerprs(ang, wl, M, ec, eta, nulam, pitch) :

        theta = ang
        lambda0 = wl

        sinteta = sin(theta)

        amp_te =0j*wl
        amp_tm =array(amp_te)

        nmax = 2*M+1 # a verifier (seb : OK)
        Ncouches = len(ec) # est-on vraiment obligé d'avoir cette variable?


        k = 2*pi/lambda0

        nu1 = nulam[0,0] # nu1 =G.nulam[0,0]

        alpham =nu1*sinteta + arange(-M,M+1).reshape(-1,1)*lambda0/pitch  # for lambda0 in G.spectre:

        # matrice de passage S entre les couches
        S_te_lite, rp_te, S_tm_lite, rp_tm = \
        S_rp_1D(M, nmax, Ncouches, ec, eta, nulam, nu1, alpham,k)

        #if '(EXTENDED)' not in G.optim_mthd:  G.infolog +=' . S_rp_1D :  %.1f ms\n' % (1000*(toc-tic2))

        # efficacitÈs de diffraction
        ete1,eteS,ordr1_te,ordrS_te,costeta1_te,costetaS_te,S11_te = \
        efficacite_ref_1D(rp_te, S_te_lite, M, nmax) # mode TE

        etm1,etmS,ordr1_tm,ordrS_tm,costeta1_tm,costetaS_tm,S11_tm = \
        efficacite_ref_1D(rp_tm,S_tm_lite, M, nmax) # mode TM

        # amplitude des modes diffractÈs
        amp_te =S11_te[ (ordr1_te==0) ]  # remplissage de la matrice d'amplitude en TE
        amp_tm =S11_tm[ (ordr1_tm==0) ]  # remplissage de la matrice d'amplitude en TM


        # comp MMFE-S4, seb
        #print ("(ordr1_te, tm==0)", ordr1_te, ordr1_tm)
        #print ("amp_te, tm", amp_te, amp_tm)
        #print ("S11", S11_te, S11_tm)


        #if '(EXTENDED)' not in G.optim_mthd:  G.infolog +=' . rprs %s (1pt) :  %.1f ms\n' % (G.method,1000*(toc-tic1))


        #rp_rs +=[amp_tm/amp_te] #????? a verifier, marche vraiment si plusieurs angles?

        rprs     = amp_tm/amp_te # Ici rprps équivalent à rp/rs, non? rpors=abs(amp_tm)/abs(amp_te) ??
        rp_polar = abs(amp_tm) # !!!!modif JH
        rs_polar = abs(amp_te) # !!!!modif JH
        rpors    = abs(amp_tm)/abs(amp_te) # !!!!modif JH

        
        # @TODO : enlever ces [0], d'ou viennent-ils maintenant????

        return (rprs[0], rp_polar[0], rs_polar[0], rpors[0]) # !!!!modif JH

        """
        d = d.append(pd.DataFrame({'angle' : np.repeat(theta,len(wl)),
                                   'Hv' : wl,
                                  'rp' : np.real(rprs),
                                  'rs' : np.imag(rprs)}))
        """



    d = []

    # Si context est un dictionnaire, on le transforme en DataFrame

    ctx = context
    if isinstance(context, dict): # traditionnel

        import itertools
        ctx = pd.DataFrame(list(itertools.product(context['angle'], context['Hv'], [context['M']])))
        ctx.columns = ['angle', 'Hv', 'M']


    for idx, row in ctx.iterrows():
        wl = row['Hv']
        ang = row['angle']
        M = int(row['M'])

        nulam = mat2cx(nl, wl)
        (rprs,rp_polar,rs_polar,rpors)  = singlerprs(ang, wl, M, ec, eta, nulam, pitch) # !!!!modif JH

       # d.append ([ang, wl, M, np.real(rprs), np.imag(rprs), np.real(rprs)/np.imag(rprs) ] )
        d.append ([ang, wl, M, np.real(rprs), np.imag(rprs),rp_polar,rs_polar, rpors ] ) # !!!!modif JH

    d = np.array(d)
    return pd.DataFrame(d, columns = ["angle","Hv", "M", "rp", "rs", "rp_polar", "rs_polar", "rpors"]) # !!!!modif JH


import sys

def main(args):
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
   

