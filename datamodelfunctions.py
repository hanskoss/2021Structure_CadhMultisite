# -*- coding: utf-8 -*-
"""
This collection establishes and deals with parameter axes objects, property/parameter
objects, and contains a function for running parallelized fitting procedures
(might be moved to a different file in future).

List of classes and functions
Mapping (class)
propertyaxis (class)
GeneratePropertyAxesCollection (function)
params (class)
sitemaker, listpermut, coupllists (small helper functions)
parammake (function)
parallelmultifit4
resc2param

"""

from __future__ import division
import globalfitfunctions as hkfit2
from hkimports2 import flatten
from hkimports2 import np
import itertools
import multiprocessing

setprotoncpmg=0

class Mapping(dict):
    """this mapping enables us to define dictionary-like objects"""
    def __setitem__(self, key, item):
        self.__dict__[key] = item
    def __getitem__(self, key):
        return self.__dict__[key]
    def __repr__(self):
        return repr(self.__dict__)
    def __len__(self):
        return len(self.__dict__)
    def __delitem__(self, key):
        del self.__dict__[key]
    def clear(self):
        return self.__dict__.clear()
    def copy(self):
        return self.__dict__.copy()
    def has_key(self, k):
        return k in self.__dict__
    def update(self, *args, **kwargs):
        return self.__dict__.update(*args, **kwargs)
    def keys(self):
        return self.__dict__.keys()
    def values(self):
        return self.__dict__.values()
    def items(self):
        return self.__dict__.items()
    def pop(self, *args):
        return self.__dict__.pop(*args)
    def __cmp__(self, dict_):
        return self.__cmp__(self.__dict__, dict_)
    def __contains__(self, item):
        return item in self.__dict__
    def __iter__(self):
        return iter(self.__dict__)
    def __unicode__(self):
        return unicode(repr(self.__dict__))

class propertyaxis(Mapping):
    """this class defines property dimensions as objects, see paper"""
    def __init__(self):#,praxes):
        pass
    def addaxis(self,praxes,*args):
        if praxes[0] not in self.keys():
            self[praxes[0]]={}
        for i in praxes[1]+['num']:
            if i not in self[praxes[0]].keys():
                self[praxes[0]][i]=[]
    def addprop(self,pr1,proplists):
        for i,j in proplists:
            #self[pr1]['num']=range(0,i+1) if len(self[pr1]['num'])<i+1 else self[pr1]['num']
            if i not in self[pr1]['num']:
                self[pr1]['num'].append(i)
            for k,l in j:
                self.addaxis([pr1,[k]])
                self[pr1][k].append(l)
    def addproplist(self,pr1,propsets):
        for i,j in propsets:
            if len(np.shape(j)) == 1:
                m=range(len(j))
            for k,l in enumerate(j):
                self.addprop(pr1,[[m[k],[[i,l]]]])
    def fastaddlist(self,pr1,pr2,pr3):
        self.addaxis([pr1,[pr2]])
        self.addproplist(pr1,[[pr2,pr3]])
    def fastaddlists(self,pr1,*args):
        for i,j in [[i,i+1] for i in np.arange(0,len(args),2)]:
            self.fastaddlist(pr1,args[i],args[j])
    def display(self,*args):
        dispthose=self.keys() if len(args) == 0 else args
        for i in dispthose:
            for j in self[i].keys():
                pass
                #print 'property key', j, ' ', self[i][j]
    def clear(self,*args):
        delthose=self.keys() if len(args) == 0 else args
        for i in delthose:
            del self[i]
    def getactforfit(self,pr1):
        actset=[]
        for i in np.arange(len(self[pr1]['num'])):
            if 'fitact' in self[pr1]:
                actset.append(self[pr1][i].fitact[i])
            else:
                actset.append(1)
        return actset

    def getnforfit(self,pr1):
        return np.arange(len(self[pr1]['num']))

def GeneratePropertyAxesCollection(reslist,resnolist):
    """This generates a specific set of property axes as required for our
    specifi work"""
    PropAxesColl=propertyaxis()
    PropAxesColl.fastaddlists('sites','name',['A','B','C'])
    PropAxesColl.fastaddlists('sites2','name',['A','B','C'])
    PropAxesColl.fastaddlists('residues','name',reslist,'seqno',resnolist)
    if setprotoncpmg == 0:
        PropAxesColl.fastaddlists('B1field','rounded',[50,70,80,90],'field',[50.12,70.434,80.351,90.2341])
        PropAxesColl.fastaddlists('B0field','rounded',[500,700,800,900],'field',[501.2,704.34,803.51,902.341])
    else:
        PropAxesColl.fastaddlists('B1field','rounded',[500,600,800,900],'field',[501.2,602.34,803.51,902.341])
        PropAxesColl.fastaddlists('B0field','rounded',[500,600,800,900],'field',[501.2,602.34,803.51,902.341])
    PropAxesColl.fastaddlists('temp','value',[285,298])
    PropAxesColl.fastaddlists('conc','value',['2.475','9.9'])
    PropAxesColl.fastaddlists('w1field','rounded',[25,50,2500],'field',[25,50,2500],'type',['CEST','CEST','HH'],'posno',[0,1,0])
    PropAxesColl.fastaddlists('TR','name',['T','X']) #'S'
    PropAxesColl.fastaddlists('type','name',['cpmg','Rex','cest'])
    PropAxesColl.display()
    return PropAxesColl


class params(Mapping):
    """
    parameter objects. also see paper.
    #The system of par, parbds and bounds works the following way.
    #single value for par provided:
    #When a single par value is provided, then this is the starting value and no
    #random generation of a starting parameter is supposed to happen.
    #In that case, parbds have to be set tightly.
    #seperate from that the question whether this value is supposed to be
    #restrained during the fitting process. This is achieved by setting bounds
    #to [].
    #For bouth parbds and bounds, one-sided bound settings are possible. The
    #determination whether the provided bound is lower upper determines the
    #relationship to par, which is a required setting.
    #no value for par provided:
    #with parbds, and bounds provided (must be two): random generation in the
    given range. bonds has separate range
    """
    
    def __init__(self):#,praxes):
        pass
    def checkbounds(self,pname):
            if len(self[pname][-1]['parbds']) == 0 and len(self[pname][-1]['bounds']) == 2:
                self[pname][-1]['parbds']=self[pname][-1]['bounds']
            if len(self[pname][-1]['bounds']) > 0 and len(self[pname][-1]['parbds']) > 0:
                if 'bounds' in self[pname][-1] and 'parbds' in self[pname][-1]:
                    if self[pname][-1]['bounds'][1] < self[pname][-1]['parbds'][1] or self[pname][-1]['bounds'][0] > self[pname][-1]['parbds'][0]:
                        print 'check boundaries for bounds and parbds' , pname, l, self[pname][-1]['bounds'], self[pname][-1]['parbds']
            if len(self[pname][-1]['par']) > 0 and len(self[pname][-1]['bounds']) > 0:
                if 'bounds' in self[pname][-1] and 'par' in self[pname][-1]:
                    if self[pname][-1]['bounds'][1] < self[pname][-1]['par'][0] or self[pname][-1]['bounds'][0] > self[pname][-1]['par'][0]:
                        print 'check boundaries for bounds and par' , pname, l, self[pname][-1]['bounds'], self[pname][-1]['par']
            if len(self[pname][-1]['parbds']) > 0 and len(self[pname][-1]['par']) > 0:
                if 'parbds' in self[pname][-1] and 'par' in self[pname][-1]:
                    if self[pname][-1]['parbds'][1] < self[pname][-1]['par'][0] or self[pname][-1]['parbds'][0] > self[pname][-1]['par'][0]:
                        print 'check boundaries for parbds and par' , pname, l, self[pname][-1]['parbds'], self[pname][-1]['par']

    def addpar(self,pname,axnames,axislist,**kwargs):
        self[pname]=[] if pname not in self.keys() else self[pname]
        for l,i in enumerate(axislist):
            self[pname].append({})
            self[pname][-1]['axnames']=axnames
            for k,j in enumerate(axnames):
                self[pname][-1][j]=i[k]
            #print kwargs
            if 'parbdslist' in kwargs.keys():
                self[pname][-1]['parbds']=kwargs['parbdslist'][l]
            else:
                self[pname][-1]['parbds']=kwargs['parbds'] if 'parbds' in kwargs else []
            if 'parlist' in kwargs.keys():
                self[pname][-1]['par']=kwargs['parlist'][l]
            else:
                self[pname][-1]['par']=kwargs['par'] if 'par' in kwargs else []
            if 'boundlist' in kwargs.keys():
                self[pname][-1]['bounds']=kwargs['boundlist'][l]
            else:
                self[pname][-1]['bounds']=kwargs['bounds'] if 'bounds' in kwargs.keys() else []
            if 'parbndtol' in kwargs.keys():
                self[pname][-1]['parbndtol']=kwargs['parbndtol']
            else:
                if pname == 'p':
                    self[pname][-1]['parbndtol']=0.001
                else:
                    self[pname][-1]['parbndtol']=1
            if 'bndtol' in kwargs.keys():
                self[pname][-1]['bndtol']=kwargs['bndtol']
            else:
                self[pname][-1]['bndtol']=self[pname][-1]['parbndtol']
            if 'fitactive' in kwargs.keys():
                self[pname][-1]['fitactive']=kwargs['fitactive'][l]
            else:
                self[pname][-1]['fitactive']=1

            self.checkbounds(pname)
    def generatemissing(self,pname):
        for j,i in enumerate(self[pname]):
            if 'par' in i and len(i['par']) > 0:
                #if 'parbds' not in self[pname][j] or len(self[pname][j]['parbds']) == 0:
                if 'parbds' not in self[pname][j]:
                    self[pname][j]['parbds']=[self[pname][j]['par'][0]-self[pname][j]['parbndtol'],self[pname][j]['par'][0]+self[pname][j]['parbndtol']]
                if 'bounds' not in self[pname][j] or len(self[pname][j]['bounds']) == 0:
                 #   print 'gothere'
                    self[pname][j]['bounds']=[]
                    for parx in self[pname][j]['par']:
                     #   print self[pname][j]['bndtol'], '!eee!'
                        self[pname][j]['bounds'].append([parx-self[pname][j]['bndtol'],parx+self[pname][j]['bndtol']])
                else:
                    newbnd=[]
                    for numb,bndx in enumerate(self[pname][j]['bounds']):
                        if len(bndx) == 0:
                            newbnd.append([self[pname][j]['par'][numb]-self[pname][j]['bndtol'],self[pname][j]['par'][numb]+self[pname][j]['bndtol']])
                        elif len(bndx) == 1: #might not work with dwb
                            if bndx[0] < self[pname][j]['par'][numb]:
                                newbnd.append([bndx[0],parx+self[pname][j]['bndtol']])
                            else:
                                newbnd.append([parx-self[pname][j]['bndtol'],bndx[0]])
                        else:
                            if pname == 'dw' and (not self[pname][j]['par'][numb] > bndx[0] or not self[pname][j]['par'][numb] < bndx[1]):
                               # print [-bndx[1],-bndx[0]]
                                newbnd.append([-bndx[1],-bndx[0]])
                            else:
                                newbnd.append(bndx)
                    self[pname][j]['bounds']=newbnd
#                    if len(self[pname][-1]['parbds']) == 0 and len(self[pname][-1]['bounds']) == 2:
#                        self[pname][-1]['parbds']=self[pname][-1]['bounds']
            else:
                if 'parbds' not in self[pname][j] or len(self[pname][j]['parbds']) == 0 :
                    self[pname][j]['parbds']=self[pname][j]['bounds']
                if 'bounds' not in self[pname][j] or len(self[pname][j]['bounds']) == 0: #at this stage not equipped to fix randomized parameter in place
                    self[pname][j]['bounds']=self[pname][j]['parbds']
#    else:
#        if parbds exist but not bnds, bnds=parbds:
 #       if bnds exist but not parbds, bnds=parbds

    def randomizepar(self,pname,numb):
        for j,i in enumerate(self[pname]):
            bdscoll=i['parbds']
            #print bdscoll, j, pname
            if pname == 'dwsign':
                self[pname][j]['par']=np.random.choice(bdscoll,numb)#np.array([np.random.choice(bdscoll) for i in np.arange(numb)])
            elif pname == 'dw':
                #print np.random.random(numb)*(bdscoll[1]-bdscoll[0])+bdscoll[0], '!!!!', bdscoll[1]-bdscoll[0], bdscoll[0]
#                print
                #rnddat=np.random.random(numb)*
                #print rnddat, self['dwsign'][j]['par'],bdscoll[0]
                self[pname][j]['par']=(self['dwsign'][j]['par']*np.random.random(numb))*(bdscoll[1]-bdscoll[0])+bdscoll[0]*self['dwsign'][j]['par']

            else:
                self[pname][j]['par']=(np.random.random(numb)*(bdscoll[1]-bdscoll[0])+bdscoll[0])
            if len(np.shape(self[pname][j]['bounds'])) == 1 or len(self[pname][j]['bounds']) != len(self[pname][j]['par']):
                bndsx=self[pname][j]['bounds'] if len(np.shape(self[pname][j]['bounds'])) == 1 else self[pname][j]['bounds'][0]
                self[pname][j]['bounds']=[bndsx for i in np.arange(len(self[pname][j]['par']))]
            try:
                self.generatemissing(pname)
            except:
                pass

    def setonepar(self,pname,datlist,**kwargs):
        #print self[pname]
        inactivecount=0
        for j,i in enumerate(self[pname]):
            if self[pname][j]['fitactive'] == 1:
                if pname == 'dw':
                    self[pname][j]['par']=[datlist[j-inactivecount]]#(self['dwsign'][j]['par']*np.random.random(numb))*(bdscoll[1]-bdscoll[0])+bdscoll[0]*self['dwsign'][j]['par']
                    self['dwsign'][j]['par']=0
                else:
                    self[pname][j]['par']=[datlist[j-inactivecount]]
                if len(np.shape(self[pname][j]['bounds'])) == 1 or len(self[pname][j]['bounds']) != len(self[pname][j]['par']):
                    bndsx=self[pname][j]['bounds'] if len(np.shape(self[pname][j]['bounds'])) == 1 else self[pname][j]['bounds'][0]
                    self[pname][j]['bounds']=[bndsx for i in np.arange(len(self[pname][j]['par']))]
                if 'bounds' in kwargs.keys():
                    self[pname][j]['bounds']=kwargs['bounds']
                if 'bndtol' in kwargs.keys():
                    self[pname][j]['bndtol']=kwargs['bndtol']
                if 'resetbounds' in kwargs.keys():
                    self[pname][j]['bounds']=[]
        #        print self[pname][j]['bounds'], 'ye now1',self[pname][j]['parbds']
                self.generatemissing(pname)
       #         print self[pname][j]['bounds'], 'ye now2',self[pname][j]['parbds']
                if 'resetparbounds' in kwargs.keys(): #resetparbounds
                    self[pname][j]['parbds']=self[pname][j]['bounds']
            else:
                inactivecount+=1
            

    def getparandbnds(self,PropAxesColl,pname,**kwargs):
  #      print kwargs['inclfilt']
        plist=[]
        blist=[]
        nlist=[]
        mapppar2=[]
        countformap=0
        mapppar=[]
        for j,i in enumerate(self[pname]):
            plistx=[]
            blistx=[]
            if 'par' in i and 'bounds' in i:
                if len(i['par']) > 0 and len(i['bounds'][0])== 2 and i['fitactive'] == 1:
  #                  print i['fitactive'], 'ii', i
                    nottake=0
                    if 'inclfilt' in kwargs:
                        for k in kwargs['inclfilt']:

                            if k[0] in self[pname][j] and PropAxesColl[k[0]][k[1]][self[pname][j][k[0]]] not in k[2]:
                            #    print pname, k[2]
                            #    if pname == 'R2mult' and k[2] == [50]:
                            #        nottake=2
                            #    if nottake != 2:
                                nottake=1
                 #   print nottake, 'nottake'
                    if nottake == 0:
                        mapppar.append(countformap)
                        for numb in range(len(i['par'])):
                            plistx.append(i['par'][numb])
                            blistx.append(i['bounds'][numb])
                    elif nottake == 2:
                        mapppar.append('X')
                    else:
                        nlist.append(j)
                    countformap+=1


                else:
                    nlist.append(j)
            else:
                nlist.append(j)
            if plistx != []:
                plist.append(plistx)
                blist.append(blistx)
        mapppar.append(countformap)
        plist=np.transpose(plist)
        if pname == 'R2mult' or pname == 'R20500':
            try:
                return [plist,np.transpose(blist)[0],np.transpose(blist)[1],flatten(nlist),mapppar]
            except:
                return [[[1]],[[0.999]],[[1.001]],flatten(nlist),mapppar]
        else:
#            print 'a', plist, blist, nlist, mapppar
            return [plist,np.transpose(blist)[0],np.transpose(blist)[1],flatten(nlist),mapppar]

    def setallparandbnds(self,paralist,pnameorder):
        pcollect=[];b1collect=[];b2collect=[];ncollect=[];mapcollect=[]
        for pname in pnameorder:
            if 'inclfilt' in kwargs:
                a,b1,b2,c,d=self.getparandbnds(pname,inclfilt=kwargs['inclfilt'])
            else:
                a,b1,b2,c,d=self.getparandbnds(pname,kwargs)
            pcollect.append(a)
            b1collect.append(b1)
            b2collect.append(b2)
            ncollect.append(c)
            mapcollect.append(d)
        pcollect=[flatten([[k for k in l[j]] for l in pcollect]) for j,i in enumerate(pcollect[0])]
        b1collect=[flatten([[k for k in l[j]] for l in b1collect]) for j,i in enumerate(b1collect[0])]
        b2collect=[flatten([[k for k in l[j]] for l in b2collect]) for j,i in enumerate(b2collect[0])]
        return pcollect,b1collect,b2collect,ncollect,mapcollect

    def getallparandbnds(self,PropAxesColl,pnameorder,**kwargs):
        #global pcollect
        pcollect=[];b1collect=[];b2collect=[];ncollect=[];mapcollect=[]
        maxlen=0
        for pname in pnameorder:
            if 'inclfilt' in kwargs:
                a,b1,b2,c,d=self.getparandbnds(PropAxesColl,pname,inclfilt=kwargs['inclfilt'])
            else:
                a,b1,b2,c,d=self.getparandbnds(PropAxesColl,pname,kwargs)
            if len(a) > maxlen:
                maxlen=len(a)
            if len(a) < maxlen:
                pcollect.append(np.array(a*maxlen))
                b1collect.append(np.array(b1*maxlen))
                b2collect.append(np.array(b2*maxlen))
                ncollect.append(c)
                mapcollect.append(d)
                
            else:
                pcollect.append(a)
                b1collect.append(b1)
                b2collect.append(b2)
                ncollect.append(c)
                mapcollect.append(d)
     #   return pcollect
     #   print pcollect, 'pcoll'
#        print pcollect[0]
 #       print [[l for l in pcollect] for j,i in enumerate(pcollect[0])]
        pcollect=[flatten([[k for k in l[j]] for l in pcollect]) for j,i in enumerate(pcollect[0])]
        b1collect=[flatten([[k for k in l[j]] for l in b1collect]) for j,i in enumerate(b1collect[0])]
        b2collect=[flatten([[k for k in l[j]] for l in b2collect]) for j,i in enumerate(b2collect[0])]
        return pcollect,b1collect,b2collect,ncollect,mapcollect

def sitemaker(n,named):
    """some funny helper function"""
    if n == 3:
        if named == 'triang':
            return [[0,1],[0,2],[1,2]]
        
def listpermut(*lists):
    """helper function"""
    step1=list(itertools.product(*lists))
    return [flatten(i) for i in step1]

def coupllists(a,b):
    """helper function"""
    nl=[]
    for j,i in enumerate(a):
        nl.append(i+b[j])
    return nl

def parammake(PropAxesColl,parbdslistb,parbdslista,k12min,k12max,k13min,k13max,k23min,k23max,pmin,pmax):
    """create initial parameters and boundaries as needed (requires editing)
    (requires optimization and simplification in an update"""
    R20500TRmin=2
    R20500TRmax=11
    R20500nTRmin=2
    R20500nTRmax=12
    R2multmin=1
    R2multmax=3
    
    paramsxx=[]
    for u in np.arange(100):
        paramsx=params() #100,400 k23
    #    paramsx.addpar('k',['sites','sites2','conc'],listpermut(sitemaker(3,'triang'),PropAxesColl.getnforfit('conc')),parbdslist=[[2000,12000],[2000,12000],[1,3500],[1,3500],[1,3500],[1,3500]])#,bounds=[20,400])
        paramsx.addpar('k',['sites','sites2','conc'],listpermut(sitemaker(3,'triang'),PropAxesColl.getnforfit('conc')),parbdslist=[[k12min,k12max],[k12min,k12max],[k13min,k13max],[k13min,k13max],[k23min,k23max],[k23min,k23max]])#,bounds=[20,400])
        paramsx.generatemissing('k')
        
        
        paramsx.addpar('p',['sites','conc'],listpermut(PropAxesColl.getnforfit('sites'),PropAxesColl.getnforfit('conc')),parbds=[pmin,pmax],fitactive=flatten([np.product(i) for i in listpermut([0,1,1],PropAxesColl.getactforfit('conc'))]))#,bounds=[20,400])
        if parbdslistb != 0:
            paramsx.addpar('dw',['residues','sites'],listpermut(PropAxesColl.getnforfit('residues'),PropAxesColl.getnforfit('sites')),boundlist=flatten(parbdslista,levels=1),parbdslist=flatten(parbdslistb,levels=1),fitactive=flatten([np.product(i) for i in listpermut(PropAxesColl.getactforfit('residues'),PropAxesColl.getactforfit('conc'),[0,1,1])]))#,bounds=[20,400])
        else:
            paramsx.addpar('dw',['residues','sites'],listpermut(PropAxesColl.getnforfit('residues'),PropAxesColl.getnforfit('sites')),parbdslist=flatten(parbdslista,levels=1),fitactive=flatten([np.product(i) for i in listpermut(PropAxesColl.getactforfit('residues'),PropAxesColl.getactforfit('conc'),[0,1,1])]))#,bounds=[20,400])
        paramsx.addpar('dwsign',['residues','sites'],listpermut(PropAxesColl.getnforfit('residues'),PropAxesColl.getnforfit('sites')),parbds=[1],fitactive=flatten([np.product(i) for i in listpermut(PropAxesColl.getactforfit('residues'),[0,1,1])]))#,bounds=[20,400])
        paramsx.addpar('R20500',['residues','TR'],listpermut(PropAxesColl.getnforfit('residues'),PropAxesColl.getnforfit('TR')),parbdslist=[[R20500TRmin,R20500TRmax] if i[1] == 0 else [R20500nTRmin,R20500nTRmax] for i in listpermut(PropAxesColl.getnforfit('residues'),PropAxesColl.getnforfit('TR'))])
        paramsx.addpar('R2mult',['residues','B1field'],listpermut(PropAxesColl.getnforfit('residues'),PropAxesColl.getnforfit('B1field')),parbds=[R2multmin,R2multmax],fitactive=[1 if i[1] > 0 else 0 for i in listpermut(PropAxesColl.getnforfit('residues'),PropAxesColl.getnforfit('B1field'))])
        for p in ['p','k','dw','dwsign','R20500','R2mult']:
            paramsx.generatemissing(p)
        for p in ['p','k','dwsign','dw','R20500','R2mult']:
            paramsx.randomizepar(p,20)
        for p in ['p','k','dw','R20500','R2mult']:
            paramsx.generatemissing(p)
        paramsxx.append(paramsx)
    
    #'/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/'
    #spinsystems,setlabels=prepro.launch('C:\\Users\\Hans\\Desktop\\TRANSFER\\2020Feb\\','test2.dat')
    #spinsystems,setlabels=prepro.launch('/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/','combo10.dat')
    #poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr, intmin,ac,oc,rateconstpre,cond=hkio.loadeverything(['testxx'],0,decoupl=0)
    
    prxx=[]
    for paramsx in paramsxx:
        prx=paramsx
        for i in ['p','k']:
            for j in np.arange(len(paramsx[i])/2):
                for which in ['par','bounds','parbds']:
#                    print i,j*2,which,1+j*2
                    prx[i][int(1+j*2)][which]=paramsx[i][int(j*2)][which]
        prxx.append(prx)
    paramsxx=prxx
    return paramsxx



def parallelmultifit4(path2020,savstatdir,setparameters3,outernum,numrep,paramsxx,PropAxesColl):
    """parallel fitting engine"""
    processes=[]
    pn=0
    for outnum in np.arange(outernum):
        pn0=list([pn])[0]
        for runnum in np.arange(numrep):
            print 'run ', str(outnum*numrep+runnum)
            #p=parallelfit(setparameters,runnum)
            p=multiprocessing.Process(target=hkfit2.parallelfit3, args=(path2020,savstatdir,setparameters3,outnum*numrep+runnum,paramsxx[outnum*numrep+runnum],PropAxesColl))
            processes.append(p)
            processes[pn].start()
            pn+=1
        pn=list([pn0])[0]
        for runnum in np.arange(numrep):
            processes[pn].join()
            pn+=1
        print 'really done with ', str(outnum*numrep+numrep), 'runs'

def resc2param(PropAxesColl,paramsxx,resultcoll,resc2paramtype):
    """transfers a result to a property/parameter collection using a defined
    property axis system. simplifies analysis a lot, e.g. easy modification
    of bounds etc
    This function can be operated in different modes.
    Mode 2 is the most frequently used - The parameters but no other conditions
    of a given parameter set are modified on an input (typically a result
    from an earlier stage).
    Mode 1 and 3 are variations of this, however, the boundaries are then
    locked to a very small range to prevent further calculations.
    Mode 4 and 5 can be used to transfer residue-independent residue from one
    calculation with residue set A to another calculation with residue set B.
    """
    whichlist=['p','k','dw','R20500','R2mult']
    lengthofoutput = len(resultcoll)
    
    if resc2paramtype == 4 or resc2paramtype == 5:
        whichlist=['p','k']#,'dw']
        lengthofoutput = len(paramsxx)
    
    for pn in np.arange(lengthofoutput):
        posit=0
        if resc2paramtype == 4  or resc2paramtype == 5:
            paralist=resultcoll[0].x
        else:
            paralist=resultcoll[pn].x
        for i in whichlist:
            a,b1,b2,c,e=paramsxx[pn].getallparandbnds(PropAxesColl,[i],inclfilt=[])
            posit0=posit
            posit+=len(a[0])
            prlst=paralist[posit0:posit]
            prlst=[x for x in prlst]
            bnds=[]
            if resc2paramtype == 1:
                bndtl=0.1
                if i == 'k':
                    print 'k',prlst
                if i == 'p':
                    bndtl=0.0001
                if i == 'dwf':
                    bndtl=1 #0000
                    paramsxx[pn].setonepar(i,prlst,bounds=bnds,bndtol=bndtl,resetbounds=1,resetparbounds=1)
            elif resc2paramtype == 2: #lets parameter bounds untouched, and sets all parameters
                paramsxx[pn].setonepar(i,prlst,resetparbounds=0)
            elif resc2paramtype == 3:
                bndtl=0.1
                if i == 'k':
                    bndtl=1
                if i == 'p':
                    bndtl=0.0001
                if i == 'dw':
                    bndtl=1 #0000
                paramsxx[pn].setonepar(i,prlst,bounds=bnds,bndtol=bndtl,resetbounds=1,resetparbounds=1)
            elif resc2paramtype == 4:
                #paramsxx can have originated from a different set of residues than resultcoll
                if i == 'k':
                    bndtl=1
                    paramsxx[pn].setonepar(i,prlst,bounds=bnds,bndtol=bndtl,resetbounds=1,resetparbounds=1)
                if i == 'p':
                    bndtl=0.0001
                    paramsxx[pn].setonepar(i,prlst,bounds=bnds,bndtol=bndtl,resetbounds=1,resetparbounds=1)
            elif resc2paramtype == 5:
                #paramsxx can have originated from a different set of residues than resultcoll
                if i == 'k':
                    bndtl=1
                    paramsxx[pn].setonepar(i,prlst,resetparbounds=0)
                if i == 'p':
                    bndtl=0.0001
                    paramsxx[pn].setonepar(i,prlst,resetparbounds=0)
    return paramsxx