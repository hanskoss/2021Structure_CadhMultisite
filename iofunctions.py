#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 12:50:34 2019
@author: hanskoss
Some simple functions for data input and output.
"""
from __future__ import division
from hkimports2 import os
from hkimports2 import np
from hkimports2 import datetime
from hkimports2 import pickle
from hkimports2 import flatten
#savstatdir=windowsss
def savstatus2(savstatdir,name,comment1,comment2,datver,scriptver,precalc,resnam,poscoll,resultcoll,setparameters,allsets,dataset):
    status={}
    cwdx=os.getcwd()
    os.chdir(savstatdir)
    status['savdat1']=datetime.now()
    status['savdat2']=str(status['savdat1'])
    status['name']=name
    status['precalc']=precalc
    status['comment1']=comment1
    status['comment2']=comment2
    status['poscoll']=poscoll
    status['resnam']=resnam
    status['resultcoll']=resultcoll
    status['datver']=datver
    status['scriptver']=scriptver
    status['setparameters']=setparameters
    status['allsets']=allsets
    status['dataset']=dataset
    save_obj(status, name)
    os.chdir(cwdx)
    return status

def savstatus2b(savstatdir,name,resnam,poscoll,allsets,allconditions):
    status={}
    cwdx=os.getcwd()
    os.chdir(savstatdir)
    status['savdat1']=datetime.now()
    status['savdat2']=str(status['savdat1'])
    status['name']=name
    status['poscoll']=poscoll
    status['resnam']=resnam
    status['allsets']=allsets
    for i in allsets:
        print i.cost
    print 'printall', len(allsets), len(allsets[0])
    status['allconditions']=allconditions
    save_obj(status, name+'_part')
    os.chdir(cwdx)
    return status

def savss(savstatdir,spinsystems,name):
    cwdx=os.getcwd()
    os.chdir(savstatdir)
    save_obj(spinsystems, name+'_ss.dat')
    os.chdir(cwdx)

def loadss(savstatdir,name):
    cwdx=os.getcwd()
    os.chdir(savstatdir)
    ss=load_obj(name+'_ss.dat')
    os.chdir(cwdx)
    return ss
    
def save_obj(obj, name ):
    f = open(name + '.dill', 'wb+')
    pickle.dump(obj, f)#, pickle.HIGHEST_PROTOCOL)
    f.close()

def load_obj(name):
    with open(name + '.dill', 'rb') as f:
        loaded_obj = pickle.load(f)
    return loaded_obj

def loadstatus2(savstatdir,name):
    status={}
    cwdx=os.getcwd()
    os.chdir(savstatdir)#'/home/hanskoss/scripts/relaxproc/savstat')
    status=load_obj(name)
    os.chdir(cwdx)
    try:
        cond=status['allconditions']
    except:
        cond=0
    alls=status['allsets']
    return status['name'],status['resnam'],status['poscoll'],alls,cond

def loadeverything(savstatdir,filenames,fieldnumber,**kwargs):
    keylist=[key for key, value in kwargs.items()]
    valuelist=[value for key, value in kwargs.items()]
    if "progcycle" in keylist:
        pp=valuelist[keylist.index("progcycle")]
    else:
        pp=-1
    if "selproc" in keylist:
        sellist=valuelist[keylist.index("selproc")]
    else:
        sellist=np.arange(120)
    resultcoll=[];resnam=[];poscoll=[]
    allsetcoll=[]
    explen=0
    for filename in filenames:
        for selproc in sellist:
            try:
                name,resnamx,poscollx,allsets,cond=loadstatus2(savstatdir,filename+str(selproc)+'_part')
                resultcoll.append(allsets[pp])
                allsetcoll.append(allsets)
                if len(allsets) > explen:
                    explen=len(allsets)
                resnam.append(resnamx)
                poscoll.append(poscollx)
            except:
                _=0
        os.chdir(savstatdir)

    resultcoll=flatten(resultcoll,levels=1)
    resnam=flatten(resnam,levels=1)
    poscoll=flatten(poscoll,levels=1)
    return poscoll,resnam,allsetcoll,resultcoll,0,0,0,0,0,0,0,0,0,0,0,0,0,0,cond


