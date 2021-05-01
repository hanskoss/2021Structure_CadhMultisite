# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 00:57:19 2021

@author: Hans
"""

from __future__ import division
import globalfitfunctions as hkfit2
import iofunctions as hkio2
import datamodelfunctions as mainfuncts
from hkimports2 import flatten
from hkimports2 import np
import copy
from hkimports2 import multiprocessing
import multiprocessing

linux='/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/'
linux='/home/hanskoss/scripts/Github/2021Structure_CadhMultisite/exptl_data/'
windows='C:\\Users\\Hans\\Desktop\\TRANSFER\\2020Feb\\'
linuxss='/home/hanskoss/scripts/relaxproc/savstat'
#linuxss='/home/hanskoss/scripts/Cadh11_multis/savstat'
windowsss='C:\\Users\\Hans\\Desktop\\TRANSFER\\relaxproc\\savstat'
savstatdir=linuxss
path2020=linux#linux#windows

def readoutresults(reslall,resnall,pickthese,DeltaOmegaParametersBoundaries,savstatdir,namresults,conditions):
    """
    reads global fitting results from file.
    reslall: list of all residues possibly used for the results output, with a letter
        in front indicating whether this is referring to the folded (A), random
        coil (B) or dimer(C) peak.
    resnall: same list as reslall, without the letter
    pickthese: list pointing to the residues in reslall, resnall actually used
        for the output
    DeltaOmegaParametersBoundaries: boundaries for delta omegas, this is a dummy
        variable in most case. format: [residue1, residue2, ...]; residue1=
        [site A, site B, site C]; site A: [min, max boundary]
    namresults: name of results file. typically has an err_ suffix when referring
        to the error calculations based on a corresponding file without the suffix
        results files have the .dill file extension.
    DEL? conditions: list of fitting conditions, not important (dummy)
        format will me explained later. subject to removal 
    
    First, a collection of Property Axes is triggered. Containing informations
        about possible magnetic fields, temperatures, sites, concentrations,
        B1 fields, experiment types, TROSY type, residues
    Then a ParameterCollection object based on the params class, using the
        PropertyAxesColl object, is generated.
    All relaxation dispersion fitting data are loaded into a variety
        of variables.
    The resultcoll variable contains the actual results. costs are extracted.
    The resulting parameters are transferred to the ParameterCollection object.
    The ParameterCollection object is re-sorted in order of increasing fitting cost.
    Important variables including a list of fitting conditions saved in the
        results file are returned. Also resultcoll, cost list, the resorted
        Parameter Collection and Property Axes Collection.
    """
    
    PropAxesColl=mainfuncts.GeneratePropertyAxesCollection([x for y,x in \
        enumerate(reslall) if y in pickthese],[x for y,x in enumerate(resnall)\
                 if y in pickthese])
    ParameterColl=mainfuncts.parammake(PropAxesColl,0,[x for y,x in enumerate\
        (DeltaOmegaParametersBoundaries) if y in pickthese],1,10000,100,900,1,900,0.005,0.5) #2 or 900
    
    #DEL?#setparameters2=['combo10l.dat','C:\\Users\\Hans\\Desktop\\TRANSFER\\2020Feb\\',[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults]
    #DEL?#namresults=namresults+str(resnall[0])+'_'+str(nn)+'_'
    
    #A variety of pre-calculated data is loaded with this command. Not all of
    #them are in use, so no need to understand every sirngle variable here.
    poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,\
        results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr,\
        intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namresults],0,decoupl=0)
    costcoll=[i.cost for i in resultcoll]
    #resultcoll=[resultcoll[np.argsort(costcoll)[0]]]
    ParameterColl2=mainfuncts.resc2param(PropAxesColl,ParameterColl,resultcoll,2)
    ParameterColl2=[ParameterColl2[i] for i in np.argsort(costcoll)]
    costcoll=[costcoll[i] for i in np.argsort(costcoll)]
    return PropAxesColl,ParameterColl2, costcoll, resultcoll,cond
# attempt to set up a completely new run from scratch.

def produceindivplots(path2020,savstatdir,datafile,namresults,reslallx,resnallx,parbdslistaxallfx,namout,mode,k121,k122,k131,k132,k231,k232,errruns):
    #print 'a'
    for pickthis in np.arange(len(reslallx)):
     #   print 'b'
        reslall=[reslallx[pickthis]]
        resnall=[resnallx[pickthis]]
        pickthese=[0]
      #  print 'c'
        PropAxesColl=mainfuncts.GeneratePropertyAxesCollection([x for y,x in enumerate(reslall) if y in pickthese],[x for y,x in enumerate(resnall) if y in pickthese])
        parbdslistaxallf=[parbdslistaxallfx[pickthis]]
        #paramsxx=mainfuncts.parammake(PropAxesColl,0,[x for y,x in enumerate(parbdslistaxallf) if y in pickthese],1,8000,100,900,k231,k232,0.005,0.1)
        #paramsxx=mainfuncts.parammake(PropAxesColl,0,[x for y,x in enumerate(parbdslistaxallf) if y in pickthese],3000,4000,600,700,k231,k232,0.01,0.02)
        paramsxx=mainfuncts.parammake(PropAxesColl,0,[x for y,x in enumerate(parbdslistaxallf) if y in pickthese],k121,k122,k131,k132,k231,k232,0.01,0.02)
        ##conditions=[0,[1,0,0,5,5]] #individual nitrogen fits
        conditions=[0,[1,0,0,10,5]] #individualproton fits
        #setparameters2=['combo10l.dat','C:\\Users\\Hans\\Desktop\\TRANSFER\\2020Feb\\',[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults]
        if namresults == []:
            ParameterColl6=paramsxx
        else:
            poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr, intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namresults],0,decoupl=0)
            costcoll=[i.cost for i in resultcoll]
            resultcoll=[resultcoll[np.argsort(costcoll)[0]]]
        
            ParameterColl6=mainfuncts.resc2param(PropAxesColl,paramsxx,resultcoll,mode)
        for nn in [0]:#np.arange(3):
            namresultsOUT=namout+str(nn)+'_'
            setparameters3=[datafile,path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresultsOUT,0]
            hkfit2.parallelfit3(path2020,savstatdir,setparameters3,0,ParameterColl6[nn],PropAxesColl)
            dataname=datafile
            
            setparameters2=[dataname,path2020,[x for y,x in enumerate(reslall) if y in pickthese],conditions,namresultsOUT] #,'A34' #'/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/
            selset=0
            ss=hkfit2.evaluaterdfit(path2020,savstatdir,PropAxesColl,setparameters2,1,ParameterColl6[selset],0)
            hkio2.savss(savstatdir,ss,namresultsOUT)

            for errc in np.arange(errruns):
                setparameters3=[datafile,path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresultsOUT+'err_'+str(errc)+'_',namresultsOUT]
                hkfit2.parallelfit3(path2020,savstatdir,setparameters3,0,ParameterColl6[nn],PropAxesColl)
        setparameters3=['combo10l.dat',path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,'FINAL_stage4_23_','FINAL_stage3_23x_'+str(runnum)]
        p=multiprocessing.Process(target=hkfit2.parallelfit3, args=(path2020,savstatdir,setparameters3,outnum*numrep+runnum,ParameterColl6[runnum],PropAxesColl))            

def parallelindiv(path2020,savstatdir,datafile,outnum,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,k121,k122,k131,k132,k23min,k23max,errruns):
    processes=[]
    pn=0
    outnum=1
    numrep=len(reslall)
    for outnum in np.arange(1):
        pn0=list([pn])[0]
        for runnum in np.arange(numrep):
            print 'run ', str(outnum*numrep+runnum)
            #p=parallelfit(setparameters,runnum)
            namout=namoutfrag+str(runnum)+'_'
            #produceindivplots(path2020,savstatdir,'FINAL_stage4_',[reslall[runnum]],[resnall[runnum]],[DeltaOmegaParametersBoundaries[runnum]],namout,5,1,2)
            p=multiprocessing.Process(target=produceindivplots, args=(path2020,savstatdir,datafile,namin,[reslall[runnum]],[resnall[runnum]],[DeltaOmegaParametersBoundaries[runnum]],namout,fittype,k121,k122,k131,k132,k23min,k23max,errruns))
            processes.append(p)
            processes[pn].start()
            pn+=1
        pn=list([pn0])[0]
        for runnum in np.arange(numrep):
            processes[pn].join()
            pn+=1
        print 'really done with ', str(outnum*numrep+numrep), 'runs'


#%%
conditions=[0,[1,0,0,1,2]]
processes=[]
""" The parallelization for this process is explicitly coded because each of 
the 16 results from stage 5 uses a slightly different calculated result for 
resampling."""
pn=0
outnum=0
pn0=list([pn])[0]

reslallx=['A31','A32','A33','A34','A44','A45','A52','A59','A72','A81','A82','A83','A84']
resnallx=[31,32,33,34,44,45,52,59,72,81,82,83,84]
DeltaOmegaParametersBoundariesx=[[[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]]]
DeltaOmegaParametersBoundaries=[[[-1, 1], [-2645, -680.9999336824217], [-771.6725093375005, -571.6725093375005]],\
                  [[-1, 1], [309.97784129346957, 1291], [100.1194520336694, 300.1194520336694]],\
                  [[-1, 1], [-1891, -884.260978369682], [-627.1592472723664, -427.15924727236643]],\
                  [[-1, 1], [-352, 269.9601625592079], [59.27874753699962, 259.2787475369996]],\
                  [[-1, 1], [-319.28316722056434, 1571], [-268.07520696705285, -68.07520696705285]],\
                  [[-1, 1], [-3914, 158.74778262212615], [-220.95131716320975, -20.95131716320975]],\
                  [[-1, 1], [-3022, -199.07962062174886], [-192.67698328090054, 7.323016719099456]],\
                  [[-1, 1], [-30.884961621023805, 3311], [-140.52654523130875, 59.47345476869125]],\
                  [[-1, 1], [445.69464392854763, 2099], [-161.88937527571977, 38.110624724280214]],\
                  [[-1, 1], [-393.11059457992536, -193.11059457992536], [-326.822989589179, -126.822989589179]],\
                  [[-1, 1], [-96.54424808105091, 103.45575191894909], [-193.93362034233354, 6.066379657666474]],\
                  [[-1, 1], [-491.4424446372856, -291.4424446372856], [-559.9291644855452, -359.9291644855453]],\
                  [[-1, 1], [-48.16372121576645, 151.83627878423354], [-50.04867680792152, 149.9513231920785]]]



dataname='combo_hcpmg285.dat'
for basnam in ['FINAL_indiv_proton_omrelax_']:#'FINAL_indiv_proton_lock2_']:#,'FINAL_indiv_proton_omrelax_']:#
    for pickthis in np.arange(13):
        print pickthis, 'step pic'
        reslall=[reslallx[pickthis]]
        resnall=[resnallx[pickthis]]
        DeltaOmegaParametersBoundaries=[DeltaOmegaParametersBoundariesx[pickthis]]
        namresults=basnam+str(pickthis)+'_2_'
        pickthese=[0]
        PropAxesColl=mainfuncts.GeneratePropertyAxesCollection([x for y,x in enumerate(reslall) if y in pickthese],[x for y,x in enumerate(resnall) if y in pickthese])
        paramsxx=mainfuncts.parammake(PropAxesColl,0,[x for y,x in enumerate(DeltaOmegaParametersBoundaries) if y in pickthese],2930,3110,664,674,1,2,0.0154,0.0163)
        conditions=[0,[1,0,0,5,10]]
        setparameters3=[dataname,path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults,0]
        hkfit2.parallelfit3(path2020,savstatdir,setparameters3,0,paramsxx[0],PropAxesColl)
        
        poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr, intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namresults],0,decoupl=0)
        costcoll=[i.cost for i in resultcoll]
        resultcoll=[resultcoll[np.argsort(costcoll)[0]]]
        setparameters2=[dataname,path2020,[x for y,x in enumerate(reslall) if y in pickthese],conditions,namresults] #,'A34' #'/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/
        selset=0
        
        ParameterColl6=mainfuncts.resc2param(PropAxesColl,paramsxx,resultcoll,2)
        ss=hkfit2.evaluaterdfit(path2020,savstatdir,PropAxesColl,setparameters2,1,ParameterColl6[selset],0)
        hkio2.savss(savstatdir,ss,namresults)
        processes=[]
        pn=0
        pn0=list([pn])[0]
        
        numrep=10
        outnum=0
        print 'goparallel'
        for runnum in np.arange(numrep):
            print 'run ', str(outnum*numrep+runnum)
            setparameters3=[dataname,path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults+'err',namresults]
            p=multiprocessing.Process(target=hkfit2.parallelfit3, args=(path2020,savstatdir,setparameters3,outnum*numrep+runnum,ParameterColl6[runnum],PropAxesColl))
            processes.append(p)
            processes[pn].start()
            pn+=1
        pn=list([pn0])[0]
        for runnum in np.arange(numrep):
            processes[pn].join()
            pn+=1
        print 'really done with ', str(outnum*numrep+numrep), 'runs'
        
#%%global protonfit
DeltaOmegaParametersBoundariesx=[[[-1, 1], [-2645, -680.9999336824217], [-771.6725093375005, -571.6725093375005]],\
          [[-1, 1], [309.97784129346957, 1291], [100.1194520336694, 300.1194520336694]],\
          [[-1, 1], [-1891, -884.260978369682], [-627.1592472723664, -427.15924727236643]],\
          [[-1, 1], [-352, 269.9601625592079], [59.27874753699962, 259.2787475369996]],\
          [[-1, 1], [-319.28316722056434, 1571], [-268.07520696705285, -68.07520696705285]],\
          [[-1, 1], [-3914, 158.74778262212615], [-220.95131716320975, -20.95131716320975]],\
          [[-1, 1], [-3022, -199.07962062174886], [-192.67698328090054, 7.323016719099456]],\
          [[-1, 1], [-30.884961621023805, 3311], [-140.52654523130875, 59.47345476869125]],\
          [[-1, 1], [445.69464392854763, 2099], [-161.88937527571977, 38.110624724280214]],\
          [[-1, 1], [-393.11059457992536, -193.11059457992536], [-326.822989589179, -126.822989589179]],\
          [[-1, 1], [-96.54424808105091, 103.45575191894909], [-193.93362034233354, 6.066379657666474]],\
          [[-1, 1], [-491.4424446372856, -291.4424446372856], [-559.9291644855452, -359.9291644855453]],\
          [[-1, 1], [-48.16372121576645, 151.83627878423354], [-50.04867680792152, 149.9513231920785]]]  
reslallx=['A31','A32','A33','A34','A44','A45','A52','A59','A72','A81','A82','A83','A84']
resnallx=[31,32,33,34,44,45,52,59,72,81,82,83,84]

namin='FINAL_stage3_'
datafile='combo_hcpmg285.dat'

  #  print 'c'
PropAxesColl=mainfuncts.GeneratePropertyAxesCollection(reslallx,resnallx)
paramsxx=mainfuncts.parammake(PropAxesColl,0,DeltaOmegaParametersBoundariesx,1000,8000,100,900,1,2,0.005,0.1)
conditions=[0,[1,0,0,10,5]] #individualproton fits
poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr, intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namin],0,decoupl=0)
costcoll=[i.cost for i in resultcoll]
resultcoll=[resultcoll[np.argsort(costcoll)[0]] for i in np.arange(100)]
ParameterColl6=mainfuncts.resc2param(PropAxesColl,paramsxx,resultcoll,5)

namresultsOUT='protonglobal1b'
setparameters3=[datafile,path2020,[reslallx],conditions,namresultsOUT,0]
print 'a'
#%%
mainfuncts.parallelmultifit4(path2020,savstatdir,setparameters3,1,10,ParameterColl6,PropAxesColl)
#%%
import globalfitfunctions as hkfit2
DeltaOmegaParametersBoundariesx=[[[-1, 1], [-2645, -680.9999336824217], [-771.6725093375005, -571.6725093375005]],\
          [[-1, 1], [309.97784129346957, 1291], [100.1194520336694, 300.1194520336694]],\
          [[-1, 1], [-1891, -884.260978369682], [-627.1592472723664, -427.15924727236643]],\
          [[-1, 1], [-352, 269.9601625592079], [59.27874753699962, 259.2787475369996]],\
          [[-1, 1], [-319.28316722056434, 1571], [-268.07520696705285, -68.07520696705285]],\
          [[-1, 1], [-3914, 158.74778262212615], [-220.95131716320975, -20.95131716320975]],\
          [[-1, 1], [-3022, -199.07962062174886], [-192.67698328090054, 7.323016719099456]],\
          [[-1, 1], [-30.884961621023805, 3311], [-140.52654523130875, 59.47345476869125]],\
          [[-1, 1], [445.69464392854763, 2099], [-161.88937527571977, 38.110624724280214]],\
          [[-1, 1], [-393.11059457992536, -193.11059457992536], [-326.822989589179, -126.822989589179]],\
          [[-1, 1], [-96.54424808105091, 103.45575191894909], [-193.93362034233354, 6.066379657666474]],\
          [[-1, 1], [-491.4424446372856, -291.4424446372856], [-559.9291644855452, -359.9291644855453]],\
          [[-1, 1], [-48.16372121576645, 151.83627878423354], [-50.04867680792152, 149.9513231920785]]]  
reslallx=['A31','A32','A33','A34','A44','A45','A52','A59','A72','A81','A82','A83','A84']
resnallx=[31,32,33,34,44,45,52,59,72,81,82,83,84]
namin='protonglobal1'
datafile='combo_hcpmg285.dat'
PropAxesColl=mainfuncts.GeneratePropertyAxesCollection(reslallx,resnallx)
paramsxx=mainfuncts.parammake(PropAxesColl,0,DeltaOmegaParametersBoundariesx,1000,8000,100,900,1,2,0.005,0.1)
conditions=[0,[1,0,0,5,5]] #individualproton fits
poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr, intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namin],0,decoupl=0)
costcoll=[i.cost for i in resultcoll]
#resultcoll=[resultcoll[np.argsort(costcoll)[0]]]
resultcoll=[resultcoll[i] for i in np.argsort(costcoll)]
ParameterColl6=mainfuncts.resc2param(PropAxesColl,paramsxx,resultcoll,2)
namresultsOUT='protonglobal2_'
setparameters3=[datafile,path2020,[reslallx],conditions,namresultsOUT,0]
mainfuncts.parallelmultifit4(path2020,savstatdir,setparameters3,1,5,ParameterColl6,PropAxesColl)

#%%
import globalfitfunctions as hkfit2
DeltaOmegaParametersBoundariesx=[[[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]]]
reslallx=['A31','A32','A33','A34','A44','A45','A52','A59','A72','A81','A82','A83','A84']
resnallx=[31,32,33,34,44,45,52,59,72,81,82,83,84]
namin='protonglobal1'
datafile='combo_hcpmg285.dat'
PropAxesColl=mainfuncts.GeneratePropertyAxesCollection(reslallx,resnallx)
paramsxx=mainfuncts.parammake(PropAxesColl,0,DeltaOmegaParametersBoundariesx,1000,8000,100,900,1,2,0.005,0.1)
conditions=[0,[1,0,0,5,5]] #individualproton fits
poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr, intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namin],0,decoupl=0)
costcoll=[i.cost for i in resultcoll]
#resultcoll=[resultcoll[np.argsort(costcoll)[0]]]
resultcoll=[resultcoll[i] for i in np.argsort(costcoll)]
ParameterColl6=mainfuncts.resc2param(PropAxesColl,paramsxx,resultcoll,2)

namresultsOUT='protonglobal3_'
setparameters3=[datafile,path2020,[reslallx],conditions,namresultsOUT,0]
mainfuncts.parallelmultifit4(path2020,savstatdir,setparameters3,1,5,ParameterColl6,PropAxesColl)

#%%
import globalfitfunctions as hkfit2
DeltaOmegaParametersBoundariesx=[[[-1, 1], [-2645, -680.9999336824217], [-771.6725093375005, -571.6725093375005]],\
          [[-1, 1], [309.97784129346957, 1291], [100.1194520336694, 300.1194520336694]],\
          [[-1, 1], [-1891, -884.260978369682], [-627.1592472723664, -427.15924727236643]],\
          [[-1, 1], [-352, 269.9601625592079], [59.27874753699962, 259.2787475369996]],\
          [[-1, 1], [-319.28316722056434, 1571], [-268.07520696705285, -68.07520696705285]],\
          [[-1, 1], [-3914, 158.74778262212615], [-220.95131716320975, -20.95131716320975]],\
          [[-1, 1], [-3022, -199.07962062174886], [-192.67698328090054, 7.323016719099456]],\
          [[-1, 1], [-30.884961621023805, 3311], [-140.52654523130875, 59.47345476869125]],\
          [[-1, 1], [445.69464392854763, 2099], [-161.88937527571977, 38.110624724280214]],\
          [[-1, 1], [-393.11059457992536, -193.11059457992536], [-326.822989589179, -126.822989589179]],\
          [[-1, 1], [-96.54424808105091, 103.45575191894909], [-193.93362034233354, 6.066379657666474]],\
          [[-1, 1], [-491.4424446372856, -291.4424446372856], [-559.9291644855452, -359.9291644855453]],\
          [[-1, 1], [-48.16372121576645, 151.83627878423354], [-50.04867680792152, 149.9513231920785]]]  
reslallx=['A31','A32','A33','A34','A44','A45','A52','A59','A72','A81','A82','A83','A84']
resnallx=[31,32,33,34,44,45,52,59,72,81,82,83,84]
namin='protonglobal1'
datafile='combo_hcpmg285.dat'
PropAxesColl=mainfuncts.GeneratePropertyAxesCollection(reslallx,resnallx)
paramsxx=mainfuncts.parammake(PropAxesColl,0,DeltaOmegaParametersBoundariesx,1000,8000,100,900,1,2,0.005,0.1)
conditions=[0,[1,0,0,5,5]] #individualproton fits
poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr, intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namin],0,decoupl=0)
costcoll=[i.cost for i in resultcoll]
#resultcoll=[resultcoll[np.argsort(costcoll)[0]]]
resultcoll=[resultcoll[i] for i in np.argsort(costcoll)]
ParameterColl6=mainfuncts.resc2param(PropAxesColl,paramsxx,resultcoll,2)


namin='FINAL_stage3_'
poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr, intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namin],0,decoupl=0)
costcoll=[i.cost for i in resultcoll]
resultcoll=[resultcoll[np.argsort(costcoll)[0]] for i in np.arange(100)]
ParameterColl7=mainfuncts.resc2param(PropAxesColl,ParameterColl6,resultcoll,5) #4 for protonglobal4, 5 for protonglobal6


namresultsOUT='protonglobal4_'
setparameters3=[datafile,path2020,[reslallx],conditions,namresultsOUT,0]
#%%
mainfuncts.parallelmultifit4(path2020,savstatdir,setparameters3,1,3,ParameterColl6,PropAxesColl)

#%%
import globalfitfunctions as hkfit2
DeltaOmegaParametersBoundariesx=[[[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]]]
reslallx=['A31','A32','A33','A34','A44','A45','A52','A59','A72','A81','A82','A83','A84']
resnallx=[31,32,33,34,44,45,52,59,72,81,82,83,84]
namin='protonglobal1'
datafile='combo_hcpmg285.dat'
PropAxesColl=mainfuncts.GeneratePropertyAxesCollection(reslallx,resnallx)
paramsxx=mainfuncts.parammake(PropAxesColl,0,DeltaOmegaParametersBoundariesx,1000,8000,100,900,1,2,0.005,0.1)
conditions=[0,[1,0,0,5,5]] #individualproton fits
poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr, intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namin],0,decoupl=0)
costcoll=[i.cost for i in resultcoll]
#resultcoll=[resultcoll[np.argsort(costcoll)[0]]]
resultcoll=[resultcoll[i] for i in np.argsort(costcoll)]
ParameterColl6=mainfuncts.resc2param(PropAxesColl,paramsxx,resultcoll,2)

namin='FINAL_stage3_'
poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr, intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namin],0,decoupl=0)
costcoll=[i.cost for i in resultcoll]
resultcoll=[resultcoll[np.argsort(costcoll)[0]] for i in np.arange(100)]
ParameterColl6=mainfuncts.resc2param(PropAxesColl,ParameterColl6,resultcoll,5)


namresultsOUT='protonglobal5_'
setparameters3=[datafile,path2020,[reslallx],conditions,namresultsOUT,0]
mainfuncts.parallelmultifit4(path2020,savstatdir,setparameters3,1,5,ParameterColl6,PropAxesColl)


#%%




reslall=['A31','A32','A33','A34','A44','A45','A52','A59','A72','A81','A82','A83','A84']
resnall=[31,32,33,34,44,45,52,59,72,81,82,83,84]
DeltaOmegaParametersBoundaries=[[[-1, 1], [-2645, -680.9999336824217], [-771.6725093375005, -571.6725093375005]],\
                  [[-1, 1], [309.97784129346957, 1291], [100.1194520336694, 300.1194520336694]],\
                  [[-1, 1], [-1891, -884.260978369682], [-627.1592472723664, -427.15924727236643]],\
                  [[-1, 1], [-352, 269.9601625592079], [59.27874753699962, 259.2787475369996]],\
                  [[-1, 1], [-319.28316722056434, 1571], [-268.07520696705285, -68.07520696705285]],\
                  [[-1, 1], [-3914, 158.74778262212615], [-220.95131716320975, -20.95131716320975]],\
                  [[-1, 1], [-3022, -199.07962062174886], [-192.67698328090054, 7.323016719099456]],\
                  [[-1, 1], [-30.884961621023805, 3311], [-140.52654523130875, 59.47345476869125]],\
                  [[-1, 1], [445.69464392854763, 2099], [-161.88937527571977, 38.110624724280214]],\
                  [[-1, 1], [-393.11059457992536, -193.11059457992536], [-326.822989589179, -126.822989589179]],\
                  [[-1, 1], [-96.54424808105091, 103.45575191894909], [-193.93362034233354, 6.066379657666474]],\
                  [[-1, 1], [-491.4424446372856, -291.4424446372856], [-559.9291644855452, -359.9291644855453]],\
                  [[-1, 1], [-48.16372121576645, 151.83627878423354], [-50.04867680792152, 149.9513231920785]]]

k23min=1
k23max=2
namin='FINAL_stage3_'
fittype=4
namoutfrag='FINAL_indiv_proton_lock2_'
parallelindiv(path2020,savstatdir,'combo_hcpmg285.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,3000,4000,600,700,k23min,k23max,20)

DeltaOmegaParametersBoundaries=[[[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]]]
k23min=1
k23max=2
namin='FINAL_stage3_'
fittype=4
namoutfrag='FINAL_indiv_proton_omrelax_'
parallelindiv(path2020,savstatdir,'combo_hcpmg285.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,3000,4000,600,700,k23min,k23max,20)


#%%

reslall=['A31','A32','A33','A34','A44','A45','A52','A59','A72','A81','A82','A83','A84']
resnall=[31,32,33,34,44,45,52,59,72,81,82,83,84]
DeltaOmegaParametersBoundaries=[[[-1, 1], [-2645, -680.9999336824217], [-771.6725093375005, -571.6725093375005]],\
                  [[-1, 1], [309.97784129346957, 1291], [100.1194520336694, 300.1194520336694]],\
                  [[-1, 1], [-1891, -884.260978369682], [-627.1592472723664, -427.15924727236643]],\
                  [[-1, 1], [-352, 269.9601625592079], [59.27874753699962, 259.2787475369996]],\
                  [[-1, 1], [-319.28316722056434, 1571], [-268.07520696705285, -68.07520696705285]],\
                  [[-1, 1], [-3914, 158.74778262212615], [-220.95131716320975, -20.95131716320975]],\
                  [[-1, 1], [-3022, -199.07962062174886], [-192.67698328090054, 7.323016719099456]],\
                  [[-1, 1], [-30.884961621023805, 3311], [-140.52654523130875, 59.47345476869125]],\
                  [[-1, 1], [445.69464392854763, 2099], [-161.88937527571977, 38.110624724280214]],\
                  [[-1, 1], [-393.11059457992536, -193.11059457992536], [-326.822989589179, -126.822989589179]],\
                  [[-1, 1], [-96.54424808105091, 103.45575191894909], [-193.93362034233354, 6.066379657666474]],\
                  [[-1, 1], [-491.4424446372856, -291.4424446372856], [-559.9291644855452, -359.9291644855453]],\
                  [[-1, 1], [-48.16372121576645, 151.83627878423354], [-50.04867680792152, 149.9513231920785]]]

#reslall=[reslall[0]]
#resnall=[resnall[0]]
#DeltaOmegaParametersBoundaries=[DeltaOmegaParametersBoundaries[0]]
k23min=1
k23max=2
namin='FINAL_stage3_'
fittype=4
namoutfrag='FINAL_indiv_proton_lock_'
parallelindiv(path2020,savstatdir,'combo_hcpmg285.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,3000,4000,600,700,k23min,k23max)
fittype=5
namoutfrag='FINAL_indiv_proton_relax_'
parallelindiv(path2020,savstatdir,'combo_hcpmg285.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,3000,4000,600,700,k23min,k23max)
k23min=800
k23max=1100
namin='FINAL_stage3_23_'
fittype=4
namoutfrag='FINAL_indiv_proton_23lock_'
parallelindiv(path2020,savstatdir,'combo_hcpmg285.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,2500,4500,150,500,k23min,k23max)
fittype=5
namoutfrag='FINAL_indiv_proton_23relax_'
parallelindiv(path2020,savstatdir,'combo_hcpmg285.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,2500,4500,150,500,k23min,k23max)

reslall=['A31','A32','A33','A34','A44','A45','A52','A59','A72','A81','A82','A83','A84']
resnall=[31,32,33,34,44,45,52,59,72,81,82,83,84]
DeltaOmegaParametersBoundariesx=[[[-1, 1], [-2645, -680.9999336824217], [-771.6725093375005, -571.6725093375005]],\
                  [[-1, 1], [309.97784129346957, 1291], [100.1194520336694, 300.1194520336694]],\
                  [[-1, 1], [-1891, -884.260978369682], [-627.1592472723664, -427.15924727236643]],\
                  [[-1, 1], [-352, 269.9601625592079], [59.27874753699962, 259.2787475369996]],\
                  [[-1, 1], [-319.28316722056434, 1571], [-268.07520696705285, -68.07520696705285]],\
                  [[-1, 1], [-3914, 158.74778262212615], [-220.95131716320975, -20.95131716320975]],\
                  [[-1, 1], [-3022, -199.07962062174886], [-192.67698328090054, 7.323016719099456]],\
                  [[-1, 1], [-30.884961621023805, 3311], [-140.52654523130875, 59.47345476869125]],\
                  [[-1, 1], [445.69464392854763, 2099], [-161.88937527571977, 38.110624724280214]],\
                  [[-1, 1], [-393.11059457992536, -193.11059457992536], [-326.822989589179, -126.822989589179]],\
                  [[-1, 1], [-96.54424808105091, 103.45575191894909], [-193.93362034233354, 6.066379657666474]],\
                  [[-1, 1], [-491.4424446372856, -291.4424446372856], [-559.9291644855452, -359.9291644855453]],\
                  [[-1, 1], [-48.16372121576645, 151.83627878423354], [-50.04867680792152, 149.9513231920785]]]
for rs in np.arange(10):
    DO2=[]    
    for lx in reslall:
        DO2.append([[-1,1]])
        j=np.abs(np.random.normal(0,np.average([i[1][1]-i[1][0] for i in DeltaOmegaParametersBoundariesx])))
        k=np.random.normal(0,np.std(flatten([i[1] for i in DeltaOmegaParametersBoundariesx])));
        DO2[-1].append([k-j,k+j])
        #print [i for i in DeltaOmegaParametersBoundaries]
        #print [i[2][1]-i[2][0] for i in DeltaOmegaParametersBoundaries]
        j=np.abs(np.random.normal(0,np.average([i[2][1]-i[2][0] for i in DeltaOmegaParametersBoundariesx])))
        k=np.random.normal(0,np.std(flatten([i[2] for i in DeltaOmegaParametersBoundariesx])));
        #print j
        DO2[-1].append([k-j,k+j])
    DeltaOmegaParametersBoundaries=DO2

    k23min=1
    k23max=2
    namin='FINAL_stage3_'
    fittype=5
    namoutfrag='FINAL_indiv_proton_relax_s_'+str(rs)+'_'
    parallelindiv(path2020,savstatdir,'combo_hcpmg285.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,3000,4000,600,700,k23min,k23max)
    fittype=4
    namoutfrag='FINAL_indiv_proton_lock_s_'+str(rs)+'_'
    parallelindiv(path2020,savstatdir,'combo_hcpmg285.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,3000,4000,600,700,k23min,k23max)
    k23min=800
    k23max=1100
    namin='FINAL_stage3_23_'
    fittype=4
    namoutfrag='FINAL_indiv_proton_23lock_s_'+str(rs)+'_'
    parallelindiv(path2020,savstatdir,'combo_hcpmg285.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,2500,4500,150,500,k23min,k23max)
    fittype=5
    namoutfrag='FINAL_indiv_proton_23relax_s_'+str(rs)+'_'
    parallelindiv(path2020,savstatdir,'combo_hcpmg285.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,2500,4500,150,500,k23min,k23max)
#%%
DeltaOmegaParametersBoundariesx=[[[-1, 1], [-2645, -680.9999336824217], [-771.6725093375005, -571.6725093375005]],\
                  [[-1, 1], [309.97784129346957, 1291], [100.1194520336694, 300.1194520336694]],\
                  [[-1, 1], [-1891, -884.260978369682], [-627.1592472723664, -427.15924727236643]],\
                  [[-1, 1], [-352, 269.9601625592079], [59.27874753699962, 259.2787475369996]],\
                  [[-1, 1], [-319.28316722056434, 1571], [-268.07520696705285, -68.07520696705285]],\
                  [[-1, 1], [-3914, 158.74778262212615], [-220.95131716320975, -20.95131716320975]],\
                  [[-1, 1], [-3022, -199.07962062174886], [-192.67698328090054, 7.323016719099456]],\
                  [[-1, 1], [-30.884961621023805, 3311], [-140.52654523130875, 59.47345476869125]],\
                  [[-1, 1], [445.69464392854763, 2099], [-161.88937527571977, 38.110624724280214]],\
                  [[-1, 1], [-393.11059457992536, -193.11059457992536], [-326.822989589179, -126.822989589179]],\
                  [[-1, 1], [-96.54424808105091, 103.45575191894909], [-193.93362034233354, 6.066379657666474]],\
                  [[-1, 1], [-491.4424446372856, -291.4424446372856], [-559.9291644855452, -359.9291644855453]],\
                  [[-1, 1], [-48.16372121576645, 151.83627878423354], [-50.04867680792152, 149.9513231920785]]]  
    
#%%
namin='FINAL_stage3_'
namoutfrag='FINAL_stage3_indiv_restrGmore_'
fittype=4
k23min=1
k23max=2
reslall=['A81','A82','A83','A84']
resnall=[81,82,83,84]

DeltaOmegaParametersBoundaries=[[[-1, 1], [-1347, -1345], [-1116.0, -916.0]],\
           [[-1, 1], [-773.0, -771], [-1122.0, -922.0]],\
           [[-1, 1], [-1372.0, -1370], [-1287.0, -1087.0]],\
           [[-1, 1], [654.0, 656], [405.0, 605.0]]]
parallelindiv(path2020,savstatdir,'combo10l.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,k23min,k23max)

reslall=['A29','A30','A31','A32','A33','A34','A52','A53','A54','A78','A81','A82','A83','A84','A86']
resnall=[29,30,31,32,33,34,52,53,54,78,81,82,83,84,86]
DeltaOmegaParametersBoundaries=[[[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]]]

reslall=['A29','A30','A31','A32','A33','A34','A52','A53','A54','A78','A81','A82','A83','A84','A86']
resnall=[29,30,31,32,33,34,52,53,54,78,81,82,83,84,86]

namin='FINAL_stage3_23_'
namoutfrag='FINAL_stage3_23_indiv_'
fittype=5
k23min=1
k23max=900
parallelindiv(path2020,savstatdir,'combo10l.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,k23min,k23max)

namin='FINAL_stage3_'
namoutfrag='FINAL_stage3_indiv_'
fittype=5
k23min=1
k23max=2
parallelindiv(path2020,savstatdir,'combo10l.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,k23min,k23max)

namin='FINAL_stage3_'
namoutfrag='FINAL_stage3_indiv_restrkp_'
fittype=4
k23min=1
k23max=2
parallelindiv(path2020,savstatdir,'combo10l.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,k23min,k23max)


DeltaOmegaParametersBoundaries=[[[-1, 1], [3199.0, 13980], [1440.0, 1640.0]],\
           [[-1, 1], [-22500.0, -6856.0], [-6307.0, -6107.0]],\
           [[-1, 1], [-1896.0, -136.0], [-736.0, -536.0]],\
           [[-1, 1], [2641.0, 7887.0], [1569.0, 1769.0]],\
           [[-1, 1], [-389.0, -123.0], [-38.0, 162.0]],\
           [[-1,-1,],[-269,868],[-386,-186]],\
       #out    [[-1, 1], [-389.0, -5341], [-13000, 13000]],\
           [[-1, 1], [-5341, -588.0], [-518.0, -318.0]],\
           [[-1, 1], [-26000,26000], [-2151,-1951]],\
           [[-1, 1], [-2465.0, 5540.0], [-3519.0, -3319.0]],\
           [[-1, 1], [7279,12491], [7151,7351]],\
           [[-1, 1], [-1346.0, 15376.0], [-1116.0, -916.0]],\
           [[-1, 1], [-772.0, 22607.0], [-1122.0, -922.0]],\
           [[-1, 1], [-1371.0, 8777.0], [-1287.0, -1087.0]],\
           [[-1, 1], [655.0, 5131.0], [405.0, 605.0]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]]]


namin='FINAL_stage3_'
namoutfrag='FINAL_stage3_indiv_restrall_'
fittype=4
k23min=1
k23max=2

parallelindiv(path2020,savstatdir,'combo10l.dat',1,reslall,resnall,DeltaOmegaParametersBoundaries,namin,namoutfrag,fittype,k23min,k23max)



#%%
         
#%%
"""
These calculations do not restrict delta omega based on field dependent shifts (control)
"""
"""
Stage 1:"Fitting of data using very wide ranges of parameter boundaries.
Established initial parameters used in Stage 2."
Stage 2: "More restricted boundaries, very similar to stage 3."
The following section shows how stage 2 is executed. The only difference to
stage 1 is pointed out.
"""

""" select residues to be fitted; define delta O boundaries, for each residue
one line, in the order deltaO_A, deltaO_B, deltaO_C 
The boundaries (lower upper) are printed in units of s-1 at 600 MHz. These
are nitrogen experiments so actual deltaOs will be much smaller.
"""
reslall=['A14','A30','A32','A37','A38','A43','A45','A50','A53','A54','A55','A73','A77','A78','A86']
resnall=[14,30,32,37,38,43,45,50,53,54,55,73,77,78,86]
pickthese=[1,2,3,4,6,7,8,9,11,12,13,14]

DeltaOmegaParametersBoundaries=[[[-1, 1], [-26000, 26000], [-13000, 13000]],\
                  [[-1, 1], [-22500,-6855], [-6827,-6206]],\
                  [[-1, 1], [2640, 8675], [1669, 1836]],\
                  [[-1, 1], [-26000, 26000], [-13000, 13000]],\
                  [[-1, 1], [-26000, 26000], [-24000, -100]],\
                  [[-1, 1], [-790,-718], [773, 850]],\
                  [[-1, 1], [-26000, 26000], [-471,-428]],\
                  [[-1, 1], [-26000, 26000], [-13000, 13000]],\
                  [[-1, 1], [-26000, 26000], [-2255,-2050]],\
                  [[-1, 1], [-2711,-2465], [-3760,-3418]],\
                  [[-1, 1], [-133,-121], [2037, 2241]],\
                  [[-1, 1], [-26000, 26000], [1000, 13000]],\
                  [[-1, 1], [-26000, 26000], [1000, 13000]],\
                  [[-1, 1], [7279, 13742], [7251, 7976]],\
                  [[-1, 1], [-26000, 26000],[1000, 13000]]]

namresults='FINALstage2a'
conditions=[0,[3,3,3,1,1]]
PropAxesColl=mainfuncts.GeneratePropertyAxesCollection([x for y,x in \
    enumerate(reslall) if y in pickthese],[x for y,x in enumerate(resnall)\
             if y in pickthese])

"""Based on the defined property axes, parameter sets with certain boundaries
    are created. The first (commented out) set belongs to stage 1, the other
    to stage 2"""
#ParameterColl=mainfuncts.parammake(PropAxesColl,0,[x for y,x in enumerate\     run N
#    (DeltaOmegaParametersBoundaries) if y in pickthese],1,10000,100,900,1,2,0.005,0.5)
ParameterColl=mainfuncts.parammake(PropAxesColl,0,[x for y,x in enumerate\
    (DeltaOmegaParametersBoundaries) if y in pickthese],1000,8000,100,900,1,2,0.005,0.1)

"""compile setting and start global fit (stage 3a): 5 x 20 parallel calculations."""
setparameters3=['combo10l.dat',path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults,0]

mainfuncts.parallelmultifit4(path2020,savstatdir,setparameters3,3,30,ParameterColl,PropAxesColl)

"""
combine results form the stage 2a calculations, gives 100 runs, with 3 attempts
each"   
"""
reslalmall=[reslall[i] for i in pickthese]
namresultsprev='FINALstage2a'
poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,\
    results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr,\
    intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namresultsprev],0,decoupl=0)


"""
from the 90 runs, the results of the last step from each of the three 
precalculation attempts is collected. The best 90 (lowest cost) of the 270
collected results are then selected to continue the fit (stage 3b), with 10
documented steps of 5 undocumented fitting substeps.
"""

resultcoll=[]
allcostcoll=[]
for i in allsetcoll:
    allcostcoll.append([i[j].cost for j in [2,5,8]])
    
allcostcoll=flatten(allcostcoll)
sortedcosts=np.argsort(allcostcoll)

for sel in np.arange(90):
    pos2=sortedcosts[sel]%3*3+2
    pos1=int(np.floor(sortedcosts[sel]/3))
    resultcoll.append(allsetcoll[pos1][pos2])
namresults='FINAL_stage2b'
    
ParameterColl2=mainfuncts.resc2param(PropAxesColl,ParameterColl,resultcoll,2)
conditions=[0,[1,0,0,10,5]]
setparameters3=['combo10l.dat',path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults,0]

mainfuncts.parallelmultifit4(path2020,savstatdir,setparameters3,3,30,ParameterColl2,PropAxesColl)


#%%
""" with restrictions of delta omega based on field-dependent chemical shifts"""
reslall=['A14','A30','A32','A37','A38','A43','A45','A50','A53','A54','A55','A73','A77','A78','A86']
resnall=[14,30,32,37,38,43,45,50,53,54,55,73,77,78,86]
pickthese=[1,2,3,4,6,7,8,9,11,12,13,14]

DeltaOmegaParametersBoundaries=[[[-1, 1], [-26000, 26000], [-13000, 13000]],\
                  [[-1, 1], [-22500,-6855], [-6827,-6206]],\
                  [[-1, 1], [2640, 8675], [1669, 1836]],\
                  [[-1, 1], [-26000, 26000], [-13000, 13000]],\
                  [[-1, 1], [-26000, 0], [-24000, -100]],\
                  [[-1, 1], [-790,-718], [773, 850]],\
                  [[-1, 1], [-26000, 0], [-471,-428]],\
                  [[-1, 1], [-0, 26000], [-13000, 13000]],\
                  [[-1, 1], [-26000, 0], [-2255,-2050]],\
                  [[-1, 1], [-2711,-2465], [-3760,-3418]],\
                  [[-1, 1], [-133,-121], [2037, 2241]],\
                  [[-1, 1], [-26000, 0], [1000, 13000]],\
                  [[-1, 1], [0, 26000], [1000, 13000]],\
                  [[-1, 1], [7279, 13742], [7251, 7976]],\
                  [[-1, 1], [-26000, 26000],[1000, 13000]]]

namresults='FINALstage2c'
conditions=[0,[3,3,3,1,1]]
PropAxesColl=mainfuncts.GeneratePropertyAxesCollection([x for y,x in \
    enumerate(reslall) if y in pickthese],[x for y,x in enumerate(resnall)\
             if y in pickthese])

"""Based on the defined property axes, parameter sets with certain boundaries
    are created. The first (commented out) set belongs to stage 1, the other
    to stage 2"""
#ParameterColl=mainfuncts.parammake(PropAxesColl,0,[x for y,x in enumerate\     run N
#    (DeltaOmegaParametersBoundaries) if y in pickthese],1,10000,100,900,1,2,0.005,0.5)
ParameterColl=mainfuncts.parammake(PropAxesColl,0,[x for y,x in enumerate\
    (DeltaOmegaParametersBoundaries) if y in pickthese],1000,8000,100,900,1,2,0.005,0.1)

"""compile setting and start global fit (stage 3a): 5 x 20 parallel calculations."""
setparameters3=['combo10l.dat',path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults,0]

mainfuncts.parallelmultifit4(path2020,savstatdir,setparameters3,1,20,ParameterColl,PropAxesColl)

"""
combine results form the stage 2c calculations, gives 100 runs, with 3 attempts
each"   
"""
reslalmall=[reslall[i] for i in pickthese]
namresultsprev='FINALstage2c'
poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,\
    results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr,\
    intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namresultsprev],0,decoupl=0)


"""
from the 20 runs, the results of the last step from each of the three 
precalculation attempts is collected. The best 20 (lowest cost) of the 60
collected results are then selected to continue the fit (stage 3b), with 10
documented steps of 5 undocumented fitting substeps.
"""

resultcoll=[]
allcostcoll=[]
for i in allsetcoll:
    allcostcoll.append([i[j].cost for j in [2,5,8]])
    
allcostcoll=flatten(allcostcoll)
sortedcosts=np.argsort(allcostcoll)

for sel in np.arange(20):
    pos2=sortedcosts[sel]%3*3+2
    pos1=int(np.floor(sortedcosts[sel]/3))
    resultcoll.append(allsetcoll[pos1][pos2])
namresults='FINAL_stage2d'
    
ParameterColl2=mainfuncts.resc2param(PropAxesColl,ParameterColl,resultcoll,2)
conditions=[0,[1,0,0,10,5]]
setparameters3=['combo10l.dat',path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults,0]

mainfuncts.parallelmultifit4(path2020,savstatdir,setparameters3,1,20,ParameterColl2,PropAxesColl)

#%%

"""Stage 3: extended fitting procedure / minimization to account for all 16 
different possible combinations of ambiguous dwb and dwc"""

""
namresultsprev='FINAL_stage2d'

""" read in results from stage 2d"""
conditions=[0,[1,0,0,10,5]]
PropAxesColl,ParameterColl2, costcoll, resultcoll,cond=readoutresults(reslall,resnall,pickthese,DeltaOmegaParametersBoundaries,savstatdir,namresultsprev,conditions)
reslalmall=[reslall[i] for i in pickthese]

""" Filter out all results which match the final conditions for dwb and dwc
signs. Additional limitations for the sign emerge from field-dependent chemical
shift, and for residue 45 from the fact that a positive sign would give an
unusually large Gly N shift.
For each residue, the sign conditions are defined for dwb and dwc. -1 means
negative; 1 means positive; 0 means ambiguous. The following loop than checks
for each result whether all residues match the condition. If not, a large
number is added to the cost.
Then, the result with the lowest cost is selected and 16 copies of this results
are created; a corresponding parameter set is created.
"""

resbcond=[['A30',-1],['A32',1],['A37',0],['A38',-1],['A45',-1],['A50',1],['A53',-1],['A54',-1],['A73',-1],['A77',1],['A78',1],['A86',0]]
resbcondc=[['A30',-1],['A32',1],['A37',0],['A38',-1],['A45',-1],['A50',0],['A53',-1],['A54',-1],['A73',1],['A77',1],['A78',1],['A86',1]]
for dnwn,pc2 in enumerate(ParameterColl2):
    dnw=0
    for rbc in resbcond:
        dwb=flatten(pc2.getparandbnds(PropAxesColl,'dw',inclfilt=[['residues','name',rbc[0]],['sites','name','B']])[0])[0]
        if (rbc[1] < 0 and dwb > 0) or (rbc[1] > 0 and dwb < 0):
            dnw+=1
    if dnw == 0:
        for rbcn,rbc in enumerate(resbcond):
            if (rbc[1] == 0):
                dwb=flatten(pc2.getparandbnds(PropAxesColl,'dw',inclfilt=[['residues','name',rbc[0]],['sites','name','B']])[0])[0]

    else:
        costcoll[dnwn]=costcoll[dnwn]+dnw
        
costcoll2=[costcoll[np.argsort(costcoll)[0]] for i in np.arange(16)]
resultcoll2=[resultcoll[np.argsort(costcoll)[0]] for i in np.arange(16)]
ParameterColl3=[ParameterColl2[np.argsort(costcoll)[0]] for i in np.arange(16)]
ParameterColl4=mainfuncts.resc2param(PropAxesColl,ParameterColl3,resultcoll2,2)


"""The parameter set with the 16 copies of the most suited result from stage 2d
is now modified (into a copy of the set) by testing all 16 combinations of
ambiguous dwb and dwc signs."""

cntgo=0
pcn=0
ParameterColl5=copy.deepcopy(ParameterColl4)
ParameterColl6=[]
for s1var in [1,-1]:
    for s2var in [1,-1]:
        for s3var in [1,-1]:
            for s4var in [1,-1]:
                posb=0;posc=0
                for rbp,rb in enumerate(resbcond):
                    if rb[1] == 0:
                        svar = s1var if posb == 0 else s2var
                        posb+=1
                        ParameterColl5[pcn]['dw'][rbp*3+1]['par']=[ParameterColl4[pcn]['dw'][rbp*3+1]['par'][0]*svar]
                for rcp,rc in enumerate(resbcondc):
                    if rc[1] == 0:
                        svar = s3var if posc == 0 else s4var
                        posc+=1
                        ParameterColl5[pcn]['dw'][rcp*3+2]['par']=[ParameterColl4[pcn]['dw'][rcp*3+2]['par'][0]*svar]
                fragm=copy.deepcopy(ParameterColl5[pcn])
                ParameterColl6.append(fragm)
                pcn+=1

""" With the 16 sets of parameters, a 5-step, 5-substep minimization is
performed to see how the altering the sign of the 16 ambiguous dwb and dwc
signs affects the results. From previous test we know that the results will be
similar, so it is enough to perform a simple minimization.
The averages of the resulting parameter correspond to the final results.
For ambiguous dwb, dwc, there are two sets of results. Other parameters
are similar enough to be merged because the difference between them is much
smaller than the arising from the subsequent error calculation"""
    
conditions=[0,[1,0,0,5,5]]
setparameters3=['combo10l.dat',path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,'FINAL_stage3_23_',0]
#%%
mainfuncts.parallelmultifit4(path2020,savstatdir,setparameters3,1,16,ParameterColl6,PropAxesColl)

#%%
import globalfitfunctions as hkfit2
""" Stage 4 - error calculation. In this final stage, errors are obtained by
resampling. For each of the 16 outputs from the previous calculation, 5
parameter sets from randomly resampled data are calculated. This gives a
total of 80 paarametersets to then report the error by calculating the standard
deviation for each parameter """
dataname=cond[-1][-1]
setparameters2=[dataname,path2020,[x for y,x in enumerate(reslall) if y in pickthese],conditions,'FINAL_stage3_23_'] #,'A34' #'/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/
for selset in np.arange(16):
    ss=hkfit2.evaluaterdfit(path2020,savstatdir,PropAxesColl,setparameters2,1,ParameterColl6[selset],0)
    hkio2.savss(savstatdir,ss,'FINAL_stage3_23x_'+str(selset))

conditions=[0,[1,0,0,5,5]]
processes=[]
""" The parallelization for this process is explicitly coded because each of 
the 16 results from stage 5 uses a slightly different calculated result for 
resampling."""
pn=0
outnum=5
numrep=16
for outnum in np.arange(outnum):
    pn0=list([pn])[0]
    for runnum in np.arange(numrep):
        print 'run ', str(outnum*numrep+runnum)
        setparameters3=['combo10l.dat',path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,'FINAL_stage4_23_','FINAL_stage3_23x_'+str(runnum)]
        p=multiprocessing.Process(target=hkfit2.parallelfit3, args=(path2020,savstatdir,setparameters3,outnum*numrep+runnum,ParameterColl6[runnum],PropAxesColl))
        processes.append(p)
        processes[pn].start()
        pn+=1
    pn=list([pn0])[0]
    for runnum in np.arange(numrep):
        processes[pn].join()
        pn+=1
    print 'really done with ', str(outnum*numrep+numrep), 'runs'


