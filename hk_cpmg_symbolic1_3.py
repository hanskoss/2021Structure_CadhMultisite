# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 13:57:09 2018
Python 2.7; Sympy 1.1.1; Numpy 1.14.3; Scipy 1.1.1
This was tested in Spyder. UMP activated (except for hk_cpmg_symbolic1_2).
Performs symbolic pre-calculations
@author: hanskoss
"""
from __future__ import division
import sympy as smp
from hkimports2 import pickle

pickle.settings['recurse'] = True 

def save_obj(obj, name ):
    f = open(name + '.dill', 'wb+')
    pickle.dump(obj, f)#, pickle.HIGHEST_PROTOCOL)
    f.close()
    
# works 
#fileW = open('file_where_I_dump.dill', 'wb')
#pickle.dump([f_lbda, f_lbda], fileW)
#fileR = open('file_where_I_dump.dill', 'rb')
#f_lbda_loaded = pickle.load(fileR)
    #with open(name + '.dat', 'wb+') as f:
       

def load_obj(name ):
    f = open(name + '.dill', 'rb')
    loaded_obj = pickle.load(f)
    f.close()
    return loaded_obj

def getnhmat():
    return nhmatnewsc

print 'symbolic calculations ...'

dwb = smp.Symbol('dwb');dwc = smp.Symbol('dwc');dwd = smp.Symbol('dwd');pb = smp.Symbol('pb');
pc = smp.Symbol('pc');pd = smp.Symbol('pd');k12 = smp.Symbol('k12');k13 = smp.Symbol('k13');
k23 = smp.Symbol('k23');k14 = smp.Symbol('k14');k24 = smp.Symbol('k24');k34 = smp.Symbol('k34');
t = smp.Symbol('t');mode = smp.Symbol('mode');tt = smp.Symbol('tt');ttt = smp.Symbol('ttt');
typ = smp.Symbol('typ');pa=smp.Symbol('pa');dwt=smp.Symbol('dwt');dwtsq=smp.Symbol('dwtsq')

popmod=4 #can be set to 3 to limit extent of smp.Symbolic calculations.
if popmod ==4:
    scn=8
elif popmod ==3:
    scn=2
else:
    scn=1
#nhmatnewsc=load_obj('/home/hanskoss/scripts/website/CPMG/CPMGapprox/nhmatnewsc')
pypath='C:\\Users\\Hans\\Documents\\GitHub\\RDfitting\\'
nhmatnewsc=load_obj(pypath+'nhmatnewsc')

print 'finished symbolic calculations'      
