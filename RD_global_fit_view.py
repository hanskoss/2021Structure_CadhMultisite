# -*- coding: utf-8 -*-
"""
@author: Hans Koss
Python 2.7 based
This script contains cells to be evaluated in Spyder. Future adjustment to run
it directly as a python script are desired.

Viewing global fitting results
"""

from __future__ import division
import plotRD as procall2
import globalfitfunctions as hkfit2
from hkimports2 import flatten
from hkimports2 import np
import iofunctions as hkio2
import datamodelfunctions as mainfuncts
#import hkwatch


linux='/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/'
windows='C:\\Users\\Hans\\Desktop\\TRANSFER\\2020Feb\\'
linuxss='/home/hanskoss/scripts/relaxproc/savstat'
windowsss='C:\\Users\\Hans\\Desktop\\TRANSFER\\relaxproc\\savstat'
savstatdir=windowsss
path2020=windows#linux#windows

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
        (DeltaOmegaParametersBoundaries) if y in pickthese],1,10000,100,900,1,2,0.014,0.0165)
    #DEL?#setparameters2=['combo10l.dat','C:\\Users\\Hans\\Desktop\\TRANSFER\\2020Feb\\',[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults]
    #DEL?#namresults=namresults+str(resnall[0])+'_'+str(nn)+'_'
    
    #A variety of pre-calculated data is loaded with this command. Not all of
    #them are in use, so no need to understand every single variable here.
    poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,\
        results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr,\
        intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namresults],0,decoupl=0)
    costcoll=[i.cost for i in resultcoll]
    #resultcoll=[resultcoll[np.argsort(costcoll)[0]]]
    ParameterColl2=mainfuncts.resc2param(PropAxesColl,ParameterColl,resultcoll,2)
    ParameterColl2=[ParameterColl2[i] for i in np.argsort(costcoll)]
    costcoll=[costcoll[i] for i in np.argsort(costcoll)]
    return PropAxesColl,ParameterColl2, costcoll, resultcoll,cond
#%%

#probably should rewrite this as a function.
reslall=['A14','A30','A32','A37','A38','A43','A45','A50','A53','A54','A55','A73','A77','A78','A86']
resnall=[14,30,32,37,38,43,45,50,53,54,55,73,77,78,86]
pickthese=[1,2,3,4,6,7,8,9,11,12,13,14]
#pickthese=[1]

DeltaOmegaParametersBoundaries=[[[-1, 1], [-26000, 26000], [-13000, 13000]],\
                  [[-1, 1], [-22500,-6855], [-6827,-6206]],\
                  [[-1, 1], [2640, 8675], [1669, 1836]],\
                  [[-1, 1], [-26000, 26000], [-13000, 13000]],\
                  [[-1, 1], [-26000, 26000], [-24000, -100]],\
                  [[-1, 1], [-26000, 26000], [773, 850]],\
                  [[-1, 1], [-26000, -100], [-471,-428]],\
                  [[-1, 1], [-26000, 26000], [-13000, 13000]],\
                  [[-1, 1], [-26000, 26000], [-2255,-2050]],\
                  [[-1, 1], [-26000, 26000], [-3760,-3418]],\
                  [[-1, 1], [-26000, 26000], [2037, 2241]],\
                  [[-1, 1], [-26000, 26000], [1000, 13000]],\
                  [[-1, 1], [-26000, 26000], [1000, 13000]],\
                  [[-1, 1], [7279, 13742], [7251, 7976]],\
                  [[-1, 1], [-26000, 26000],[1000, 13000]]]  

namresults='multix13_k13_err_'
namresults='multix13_k13_ext_0'
namresults='multix13_NEW_M2_'
namresults='multix13_NEW_pfr_N2_'
namresults='multix13_NEW_pfr_ST2_test_'
namresults='multix13_NEW_pfr_STU_'
namresults='multix13_NEW_pfr_STU_errb_'
namresults='multix13_NEW_pfr_STU_err3f_'
namresults='multix13_NEW_pfr_STU_errb3_'
namresults='multix13_NEW_pfr_STU_err3g_'
#namresults='multix13_NEW_pfr_S2_'
#namresults='multix13_NEW_pfr_Sdeb2_'

#namresults='multix13_NEW_pfr_M2_err_'
conditions=[0,[1,3,3,3,3]]

PropAxesColl,ParameterColl2, costcoll, resultcoll,cond=readoutresults(reslall,resnall,pickthese,DeltaOmegaParametersBoundaries,savstatdir,namresults,conditions)
reslalmall=[reslall[i] for i in pickthese]

#%%
alldat=[]
for pn in np.arange(94):
    pnl=[]
    for re in resultcoll:
        pnl.append(re.x[pn])
    alldat.append(pnl)
    print np.average(pnl), np.std(pnl)
import csv
#with open('/home/hanskoss/scripts/relaxproc/savstat/fitresults_final16.csv','w') as f:
with open('/home/hanskoss/scripts/relaxproc/savstat/fitresults_errg.csv','w') as f:
    write = csv.writer(f)
    write.writerows(alldat)

with open('/home/hanskoss/scripts/relaxproc/savstat/fitresults_errgcost.csv','w') as f:
    write = csv.writer(f)
    write.writerows([costcoll])

#%%
#with open('/home/hanskoss/scripts/relaxproc/savstat/fitresults_final16cost.csv','w') as f:
with open('/home/hanskoss/scripts/relaxproc/savstat/fitresults_all100cost.csv','w') as f:
    write = csv.writer(f)
    write.writerows([costcoll])


#%%
newpc=[]
#resbcond=[['A30',-1],['A32',1],['A37',0],['A38',-1],['A45',0],['A50',1],['A53',-1],['A54',-1],['A73',-1],['A77',1],['A78',1],['A86',0]]
resbcond=[['A30',-1],['A32',1],['A37',0],['A38',-1],['A45',-1],['A50',1],['A53',-1],['A54',-1],['A73',-1],['A77',1],['A78',1],['A86',0]]
for dnwn,pc2 in enumerate(ParameterColl2):
    dnw=0
    for rbc in resbcond:
        dwb=flatten(pc2.getparandbnds(PropAxesColl,'dw',inclfilt=[['residues','name',rbc[0]],['sites','name','B']])[0])[0]
        if dnwn == 4 or dnwn == 8 or dnwn == 95:
            print dwb
        if (rbc[1] < 0 and dwb > 0) or (rbc[1] > 0 and dwb < 0):
            #print 'fail', dwb
            dnw+=1
        #else:
        #    print 'pass', dwb
    if dnw == 0:
        print 'all pass', dnwn
    else:
        print 'fail ', dnw, dnwn
        costcoll[dnwn]=costcoll[dnwn]+dnw

ParameterColl2=[ParameterColl2[i] for i in np.argsort(costcoll)]
costcoll2=[costcoll[i] for i in np.argsort(costcoll)]

#costcoll=[costcoll[i] for i in np.argsort(costcoll)]
        
        #ParameterColl2[2]['dw'][19]['par']=[-10258.957]
#%%
import plotRD as procall2
import globalfitfunctions as hkfit2
#selnam=[reslalmall[0],reslalmall[1]] #7 symbolic
selnam=['A30','A32','A37']#reslalmall[0],reslalmall[1]] #7 symbolic
selnam=['A38','A45','A38']
selnam=['A32','A38','A77','A32','A38','A77']
#selnam=['A73','A50','A86']
#selnam=['A53','A50','A73']
#selnam=['A30']
#selnam=['A53','A54','A73']
#selnam=['A77','A78','A86']
#ParameterColl2[1]['k'][0]['par']=[10000]#ParameterColl2[0]['k'][0]['par']
#ParameterColl2[1]['k'][1]['par']=[10000]#ParameterColl2[0]['k'][1]['par']
#ParameterColl2[1]['k'][2]['par']=[1]#ParameterColl2[0]['k'][2]['par']
#ParameterColl2[1]['k'][3]['par']=[1]#ParameterColl2[0]['k'][3]['par']
#ParameterColl2[1]['k'][4]['par']=[1]#ParameterColl2[0]['k'][4]['par']
#ParameterColl2[1]['k'][5]['par']=[1]#ParameterColl2[0]['k'][5]['par']

#selnam=['A53','A54','A73']
#selnam=['A77','A78','A86']

#selnam

dataname=cond[-1][-1]
dataname='combo10l.dat'
selset=0#2
setparameters2=[dataname,path2020,selnam,conditions,namresults] #,'A34' #'/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/
setparameters2=[dataname,path2020,selnam,conditions,namresults] #,'A34' #'/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/
mustevaluate = 1
if mustevaluate == 1:
    ss=hkfit2.evaluaterdfit(path2020,savstatdir,PropAxesColl,setparameters2,1,ParameterColl2[selset],0)
    hkio2.savss(savstatdir,ss,namresults)
else:
    ss=hkfit2.evaluaterdfit(path2020,savstatdir,PropAxesColl,setparameters2,1,ParameterColl2[selset],namresults) #k12_

#%%

spinsystems=ss #['x','labxx']
selnam=['A32','A38','A77']
labelstring1=[[1,'160 '+r'$\mu$'+'M\n500 MHz'],[4,'640 '+r'$\mu$'+'M\n500 MHz'],[10,'160 '+r'$\mu$'+'M\n700 MHz'],[6,'160 '+r'$\mu$'+'M\n900 MHz'],[8,'640 '+r'$\mu$'+'M\n900 MHz'],[9,'640 '+r'$\mu$'+'M\n900 MHz']]
labelstring2=[[0,'160 '+r'$\mu$'+'M, 500 MHz'],[2,'640 '+r'$\mu$'+'M, 500 MHz'],[3,'640 '+r'$\mu$'+'M, 500 MHz'],[5,'160 '+r'$\mu$'+'M, 900 MHz'],[7,'640 '+r'$\mu$'+'M, 900 MHz']]
#labelstring3=[[11,'160 '+r'$\mu$'+'M, 700 MHz, B'+r'$_1$'+' = 15 Hz'],[12,'160 '+r'$\mu$'+'M, 700 MHz, B'+r'$_1$'+' = 25 Hz'],[13,'160 '+r'$\mu$'+'M, 700 MHz, B'+r'$_1$'+' = 50 Hz'],[14,'640 '+r'$\mu$'+'M, 800 MHz, B'+r'$_1$'+' = 35 Hz']]
labelstring3=[[11,'160 '+r'$\mu$'+'M, 15 Hz'],[12,'160 '+r'$\mu$'+'M, 25 Hz'],[13,'160 '+r'$\mu$'+'M, 50 Hz'],[14,'640 '+r'$\mu$'+'M, 35 Hz']]

labelstring4=[[1,'160 '+r'$\mu$'+'M\n500 MHz'],[4,'640 '+r'$\mu$'+'M\n500 MHz'],[9,'160 '+r'$\mu$'+'M\n700 MHz'],[6,'160 '+r'$\mu$'+'M\n900 MHz'],['x','640 '+r'$\mu$'+'M\n900 MHz'],[8,'640 '+r'$\mu$'+'M\n900 MHz']]
labelstring5=[[10,'160 '+r'$\mu$'+'M, 15 Hz'],[11,'160 '+r'$\mu$'+'M, 25 Hz'],[12,'160 '+r'$\mu$'+'M, 50 Hz'],[13,'640 '+r'$\mu$'+'M, 35 Hz']]
labelstring6=[[0,'160 '+r'$\mu$'+'M, 500 MHz']]

#special plot when the last line is residue 73

procall2.plotelements(spinsystems,selnam,0,[[1,2,5],[1,2,5],[1,2,5]],[[labelstring3,labelstring2,labelstring1],[labelstring3,labelstring2,labelstring1],[labelstring3,labelstring2,labelstring1]],[0,0,0,0,0,0,0,0,0],1,104)
#procall2.plotelements(spinsystems,selnam,0,[[1,2,5,1,2,5],[1,2,5,1,2,5],[1,2,5,1,2,5],[1,2,5,1,2,5],[1,2,5,1,2,5],[1,2,5,1,2,5]],[[labelstring3,labelstring2,labelstring1],[labelstring3,labelstring2,labelstring1],[labelstring3,labelstring2,labelstring1,labelstring3,labelstring2,labelstring1]],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],104)

#procall2.plotelements(spinsystems,selnam,0,[[1,2,5],[1,2,5],[1,2,5]],[[labelstring3,labelstring2,labelstring1],[labelstring3,labelstring2,labelstring1],[labelstring5,labelstring2,labelstring4]],[0,0,0,0,0,0,0,0,0])
#procall2.plotelements(spinsystems,selnam,0,[[2,2,2],[2,2,2],[2,2,2]],[[labelstring6,labelstring6,labelstring6],[labelstring6,labelstring6,labelstring6],[labelstring6,labelstring6,labelstring6]],[0,0,0,0,0,0,0,0,0])

#procall2.plotelements(spinsystems,selnam,v,[1,2,5],[labelstring3,labelstring2,labelstring1],0)
#procall2.plotelements(spinsystems,selnam,v,2,labelstring2,1)
#procall2.plotelements(spinsystems,selnam,v,1,labelstring3,1)
#procall2.plotelements(spinsystems,selnam,v,2,[[0,'160 '+r'$\mu$'+'M\n500 MHz'],[2,'640'],[3,'640'],[5,'160'],[7,'640']],1)



#%%
import hkwatch
hkwatch.parameterplotter2(PropAxesColl,ParameterColl2,[2],cost=costcoll2,lenresults=3,resnam=reslalmall,conc=['2.475','9.9'],figno=106) #32
#%%
pbrat_collect=[]
for mm in np.arange(80):
    pb_lo=[];pb_hi=[];pc_lo=[];pc_hi=[];k12_lo=[];k12_hi=[];k13_lo=[];k13_hi=[];k23_lo=[];k23_hi=[];pb_rat=[];pc_rat=[];k12r=[];k13r=[];
    for n in np.arange(mm):
        pb1=flatten(ParameterColl2[n].getparandbnds(PropAxesColl,'p',inclfilt=[['sites','name','B'],['conc','value','2.475']]))[0]
        pb2=flatten(ParameterColl2[n].getparandbnds(PropAxesColl,'p',inclfilt=[['sites','name','B'],['conc','value','9.9']]))[0]
        pb_lo.append(pb1)
        pb_hi.append(pb2)
        pb_rat.append(pb1/pb2)
        pc1=flatten(ParameterColl2[n].getparandbnds(PropAxesColl,'p',inclfilt=[['sites','name','C'],['conc','value','2.475']]))[0]
        pc2=flatten(ParameterColl2[n].getparandbnds(PropAxesColl,'p',inclfilt=[['sites','name','C'],['conc','value','9.9']]))[0]
        pc_lo.append(pc1)
        pc_hi.append(pc2)
        pc_rat.append(pc1/pc2)
      #  print pc1, pc2, 'test'
        k121=flatten(ParameterColl2[n].getparandbnds(PropAxesColl,'k',inclfilt=[['sites','name','A'],['sites2','name','B'],['conc','value','2.475']]))[0]
        k122=flatten(ParameterColl2[n].getparandbnds(PropAxesColl,'k',inclfilt=[['sites','name','A'],['sites2','name','B'],['conc','value','9.9']]))[0]
        k131=flatten(ParameterColl2[n].getparandbnds(PropAxesColl,'k',inclfilt=[['sites','name','A'],['sites2','name','C'],['conc','value','2.475']]))[0]
        k132=flatten(ParameterColl2[n].getparandbnds(PropAxesColl,'k',inclfilt=[['sites','name','A'],['sites2','name','C'],['conc','value','9.9']]))[0]
        k12_lo.append(k121)
        k12_hi.append(k122)
        k13_lo.append(k131)
        k13_hi.append(k132)
        k12r.append(k121/k122)
        k13r.append(k131/k132)
        k23_lo.append(flatten(ParameterColl2[n].getparandbnds(PropAxesColl,'k',inclfilt=[['sites','name','B'],['sites2','name','C'],['conc','value','2.475']]))[0])
        k23_hi.append(flatten(ParameterColl2[n].getparandbnds(PropAxesColl,'k',inclfilt=[['sites','name','B'],['sites2','name','C'],['conc','value','9.9']]))[0])
    

    if mm > 1:
        #print np.percentile(np.sort(pb_lo),100-68-(100-68)/2), np.percentile(np.sort(pb_lo),100-(100-68)/2), '!'
        #print np.percentile(np.sort(pb_rat),100-68-(100-68)/2)#-np.percentile(np.sort(pb_rat),100-(100-68)/2)
      #  print np.percentile(np.sort(k12_lo),100-68-(100-68)/2), np.percentile(np.sort(k12_lo),100-(100-68)/2)
        print np.percentile(np.sort(k12_hi),100-68-(100-68)/2), np.percentile(np.sort(k12_hi),100-(100-68)/2)
        #pbrat_collect.append(np.percentile(np.sort(pb_rat),100-68-(100-68)/2)-np.percentile(np.sort(pb_rat),100-(100-68)/2))
        

print np.percentile(np.sort(pb_lo),100-68-(100-68)/2), np.percentile(np.sort(pb_lo),100-(100-68)/2)
print np.percentile(np.sort(pb_hi),100-68-(100-68)/2), np.percentile(np.sort(pb_hi),100-(100-68)/2)
print np.percentile(np.sort(pc_lo),100-68-(100-68)/2), np.percentile(np.sort(pc_lo),100-(100-68)/2)
print np.percentile(np.sort(pc_hi),100-68-(100-68)/2), np.percentile(np.sort(pc_hi),100-(100-68)/2)
print np.percentile(np.sort(pb_rat),100-68-(100-68)/2), np.percentile(np.sort(pb_rat),100-(100-68)/2)
print np.percentile(np.sort(pc_rat),100-68-(100-68)/2), np.percentile(np.sort(pc_rat),100-(100-68)/2)
print np.percentile(np.sort(k12_lo),100-68-(100-68)/2), np.percentile(np.sort(k12_lo),100-(100-68)/2)
print np.percentile(np.sort(k12_hi),100-68-(100-68)/2), np.percentile(np.sort(k12_hi),100-(100-68)/2)
print np.percentile(np.sort(k13_lo),100-68-(100-68)/2), np.percentile(np.sort(k13_lo),100-(100-68)/2)
print np.percentile(np.sort(k13_hi),100-68-(100-68)/2), np.percentile(np.sort(k13_hi),100-(100-68)/2)
print np.percentile(np.sort(k23_lo),100-68-(100-68)/2), np.percentile(np.sort(k23_lo),100-(100-68)/2)
print np.percentile(np.sort(k23_hi),100-68-(100-68)/2), np.percentile(np.sort(k23_hi),100-(100-68)/2)
print np.percentile(np.sort(k12r),100-68-(100-68)/2), np.percentile(np.sort(k12r),100-(100-68)/2)
print np.percentile(np.sort(k13r),100-68-(100-68)/2), np.percentile(np.sort(k13r),100-(100-68)/2)

#print np.percentile(np.sort(pc_lo),100-68-(100-68)/2), np.percentile(np.sort(pb_lo),100-(100-68)/2)
#%% attempt to set up a completely new run from scratch.
reslall=['A14','A30','A32','A37','A38','A43','A45','A50','A53','A54','A55','A73','A77','A78','A86']
resnall=[14,30,32,37,38,43,45,50,53,54,55,73,77,78,86]
pickthese=[1,2,3,4,6,7,8,9,11,12,13,14]
DeltaOmegaParametersBoundaries=[[[-1, 1], [-26000, 26000], [-13000, 13000]],\
                  [[-1, 1], [-22500,-6855], [-6827,-6206]],\
                  [[-1, 1], [2640, 8675], [1669, 1836]],\
                  [[-1, 1], [-26000, 26000], [-13000, 13000]],\
                  [[-1, 1], [-26000, 26000], [-24000, -100]],\
                  [[-1, 1], [-26000, 26000], [773, 850]],\
                  [[-1, 1], [-26000, -100], [-471,-428]],\
                  [[-1, 1], [-26000, 26000], [-13000, 13000]],\
                  [[-1, 1], [-26000, 26000], [-2255,-2050]],\ 
                  [[-1, 1], [-26000, 26000], [-3760,-3418]],\
                  [[-1, 1], [-26000, 26000], [2037, 2241]],\
                  [[-1, 1], [-26000, 26000], [1000, 13000]],\
                  [[-1, 1], [-26000, 26000], [1000, 13000]],\
                  [[-1, 1], [7279, 13742], [7251, 7976]],\
                  [[-1, 1], [-26000, 26000],[1000, 13000]]]  

namresults='multix13_NEW'
conditions=[0,[1,3,3,3,3]]
PropAxesColl=mainfuncts.GeneratePropertyAxesCollection([x for y,x in \
    enumerate(reslall) if y in pickthese],[x for y,x in enumerate(resnall)\
             if y in pickthese])
ParameterColl=mainfuncts.parammake(PropAxesColl,0,[x for y,x in enumerate\
    (DeltaOmegaParametersBoundaries) if y in pickthese],1,10000,100,900,1,0.005,0.5)
#%%
setparameters3=['combo10l.dat',path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults,0]
hkfit2.parallelfit3(setparameters3,0,ParameterColl[0],PropAxesColl)
#%%
import hkio2
import hkwatch
costdev=[]
pbx=[]

pcx=[]
dwbx=[]
dwcx=[]
#for nrp,namresults in enumerate(['multix13_NEW_pfr_S2_ext_']):#,'multix13_NEW_pfr_M2_ext_']):
for nrp,namresults in enumerate(['multix13_NEW_pfr_STU_err_']):#,'multix13_NEW_pfr_M2_ext_']):
#for nrp,namresults in enumerate(['multix13_NEW_pfr_M2_','multix13_NEW_pfr_M2_ext_']):
    resnam=['A30']
    selproc=0
    inc=0
    decoupl=0
    magfnum=0
    poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,\
        results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr,\
        intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namresults],0,decoupl=0)
    allsetcoll=[allsetcoll[cs] for cs in np.argsort([rc.cost for rc in resultcoll])]
    for jx,j in enumerate(np.arange(len(allsetcoll))):
        if nrp == 0:
            costdev.append([])
            pbx.append([])
            pcx.append([])
            dwbx.append([[],[],[],[],[],[],[],[],[],[],[],[]])
            dwcx.append([[],[],[],[],[],[],[],[],[],[],[],[]])

        for i in np.arange(len(allsetcoll[j])):
            costdev[jx].append(allsetcoll[j][i].cost)
            #print here
            pbx[jx].append(allsetcoll[j][i].x[0])
            
            pcx[jx].append(allsetcoll[j][i].x[2])
            dwblist=allsetcoll[j][i].x[10:10+2*11+1:2]
            dwclist=allsetcoll[j][i].x[11:11+2*11+1:2]
            for dwpos,dw in enumerate(dwblist):
                dwbx[jx][dwpos].append(dwblist[dwpos])
                dwcx[jx][dwpos].append(dwclist[dwpos])
        #    pcx[jx].append(allsetcoll[j][i].x[10])
        
#%%
from matplotlib import pyplot as plt
plt.figure(3)
plt.clf()
for cdp in np.argsort([costdev[i][-1] for i in np.arange(len(costdev))]):
    cdl=costdev[cdp]
    plt.plot(np.arange(len(cdl)),cdl,linestyle='solid')
plt.figure(5)
plt.clf()

#for cdp in np.argsort([pcx[i][-1] for i in np.arange(len(pcx))]):
for cdp in np.argsort([costdev[i][-1] for i in np.arange(len(costdev))]):
    cdl=pcx[cdp]
    plt.plot(np.arange(len(cdl)),cdl,linestyle='solid')
plt.figure(4)
plt.clf()
#for cdp in np.argsort([pbx[i][-1] for i in np.arange(len(pbx))]):
for cdp in np.argsort([costdev[i][-1] for i in np.arange(len(costdev))]):
    cdl=pbx[cdp]
    plt.plot(np.arange(len(cdl)),cdl,linestyle='solid')

for figno in np.arange(12):
    plt.figure(10+figno)
    plt.clf()
    for cdp in np.argsort([costdev[i][-1] for i in np.arange(len(costdev))]):
        cdl=dwbx[cdp][figno]
        plt.plot(np.arange(len(cdl)),cdl,linestyle='solid')

for figno in np.arange(12):
    plt.figure(30+figno)
    plt.clf()
    for cdp in np.argsort([costdev[i][-1] for i in np.arange(len(costdev))]):
        cdl=dwcx[cdp][figno]
        plt.plot(np.arange(len(cdl)),cdl,linestyle='solid')

#%% Figure out why extension of run goes up in chi square in the beginning.
namresults='multix13_NEW_pfr_R2_'
poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,\
    results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr,\
    intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namresults],0,decoupl=0)
c1=allsetcoll[0][-1].cost
x1=allsetcoll[0][-1].x
dataname=cond[-1][-1]
selnam=['A30']#,'A32','A37','A38','A45','A50','A53','A54','A73','A77','A78','A86']
selset=0
setparameters2=[dataname,path2020,selnam,conditions,namresults] #,'A34' #'/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/
mustevaluate = 1
PropAxesColl2,ParameterColl3, costcoll, resultcoll,cond=readoutresults(reslall,resnall,pickthese,DeltaOmegaParametersBoundaries,savstatdir,namresults,conditions)
#%%
import hk_fitengines1_1b as hkfit2 #4.4669783020706575, 3.8352619428032675
namresults='multix13_NEW_pfr_M2_'
setparameters2=[dataname,path2020,selnam,conditions,namresults] #,'A34' #'/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/
#5.25043918658713 5.3670662370553135
ss=hkfit2.evaluaterdfit(path2020,savstatdir,PropAxesColl2,setparameters2,1,ParameterColl3[selset],0)
ss1a=hkfit2.prepro.spinsystems1
ssa=hkfit2.prepro.spinsystems
#5.169550500484723
#5.25043918658713 5.310617455770056 4.912346663061861
#5.25043918658713 5.25043918658713 5.316890190913362 4.83710583621936
#%%
namresults='multix13_NEW_pfr_R2_'
poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,\
    results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr,\
    intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything(savstatdir,[namresults],0,decoupl=0)
c2=allsetcoll[0][0].cost
x2=allsetcoll[0][0].x  
    
dataname=cond[-1][-1]
selnam=['A30']#,'A32','A37','A38','A45','A50','A53','A54','A73','A77','A78','A86']
selset=0
setparameters2=[dataname,path2020,selnam,conditions,namresults] #,'A34' #'/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/
mustevaluate = 1
PropAxesColl,ParameterColl2, costcoll, resultcoll,cond=readoutresults(reslall,resnall,pickthese,DeltaOmegaParametersBoundaries,savstatdir,namresults,conditions)
#%%
import hk_fitengines1_1b as hkfit2
namresults='multix13_NEW_pfr_M2_ext_'
setparameters2=[dataname,path2020,selnam,conditions,namresults]
ss=hkfit2.evaluaterdfit(path2020,savstatdir,PropAxesColl,setparameters2,1,ParameterColl2[selset],0)
ss1b=hkfit2.prepro.spinsystems1
ssb=hkfit2.prepro.spinsystems
