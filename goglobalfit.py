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
from hkimports2 import copy
from hkimports2 import multiprocessing


linux='/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/'
windows='C:\\Users\\Hans\\Desktop\\TRANSFER\\2020Feb\\'
path2020=linux#windows

def readoutresults(reslall,resnall,pickthese,DeltaOmegaParametersBoundaries,namresults,conditions):
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
        (DeltaOmegaParametersBoundaries) if y in pickthese],1,10000,100,900,1,2,0.005,0.5)
    #DEL?#setparameters2=['combo10l.dat','C:\\Users\\Hans\\Desktop\\TRANSFER\\2020Feb\\',[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults]
    #DEL?#namresults=namresults+str(resnall[0])+'_'+str(nn)+'_'
    
    #A variety of pre-calculated data is loaded with this command. Not all of
    #them are in use, so no need to understand every single variable here.
    poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,\
        results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr,\
        intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything([namresults],0,decoupl=0)
    costcoll=[i.cost for i in resultcoll]
    #resultcoll=[resultcoll[np.argsort(costcoll)[0]]]
    ParameterColl2=mainfuncts.resc2param(PropAxesColl,ParameterColl,resultcoll,2)
    ParameterColl2=[ParameterColl2[i] for i in np.argsort(costcoll)]
    costcoll=[costcoll[i] for i in np.argsort(costcoll)]
    return PropAxesColl,ParameterColl2, costcoll, resultcoll,cond
# attempt to set up a completely new run from scratch.
#%%

"""
Step 1:"Fitting of data using very wide ranges of parameter boundaries.
Established initial parameters used in Step 2."
Step 2: "More restricted boundaries (including deltaO), very similar to step 3."
Step 3: "Based on the results from Step 2, boundaries were restricted a little
more"
The following section shows how step 3 is executed. The only difference to
step 2 is pointed out. Step 1 was similar, but hardly any boundaries,
also excluding deltaO boundaries.
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


"""Generate output file name; then define conditions of the run, in this case:
    3 attempts of 3 documented pre-run steps with three undocumented fitting
    steps each. The property axes for all residues are created. The property
    axes are important to define and describe parameters, as outlined elsewhere
    and in the paper.
    Note: For step 3 in the paper, this multi-run was split into three (S, T, U)"""
namresults='multix13_NEW_pfr_U2_'
conditions=[0,[3,3,3,1,1]]
PropAxesColl=mainfuncts.GeneratePropertyAxesCollection([x for y,x in \
    enumerate(reslall) if y in pickthese],[x for y,x in enumerate(resnall)\
             if y in pickthese])

"""Based on the defined property axes, parameter sets with certain boundaries
    are created. The first (commented out) set belongs to step 2, the other
    to step 3"""
#ParameterColl=mainfuncts.parammake(PropAxesColl,0,[x for y,x in enumerate\     run N
#    (DeltaOmegaParametersBoundaries) if y in pickthese],1,10000,100,900,1,2,0.005,0.5)
ParameterColl=mainfuncts.parammake(PropAxesColl,0,[x for y,x in enumerate\
    (DeltaOmegaParametersBoundaries) if y in pickthese],1000,8000,100,900,1,2,0.005,0.1)

"""compile setting and start global fit (step 3): 5 x 20 parallel calculations."""
setparameters3=['combo10l.dat',path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults,0]
mainfuncts.parallelmultifit4(setparameters3,5,20,ParameterColl,PropAxesColl)


#%%
"""
combine results form the step 3 calculations, gives 100 runs, with 3 attempts
each"   
"""
reslalmall=[reslall[i] for i in pickthese]
namresultsprev='multix13_NEW_pfr_S2_'
poscoll,resnam,allsetcoll1,resultcoll,relaxrat0,relaxrat,lookatratio,\
    results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr,\
    intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything([namresultsprev],0,decoupl=0)
namresultsprev='multix13_NEW_pfr_T2_'
poscoll,resnam,allsetcoll2,resultcoll,relaxrat0,relaxrat,lookatratio,\
    results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr,\
    intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything([namresultsprev],0,decoupl=0)
namresultsprev='multix13_NEW_pfr_U2_'
poscoll,resnam,allsetcoll3,resultcoll,relaxrat0,relaxrat,lookatratio,\
    results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr,\
    intmin,ac,oc,rateconstpre,cond=hkio2.loadeverything([namresultsprev],0,decoupl=0)

allsetcoll=allsetcoll1+allsetcoll2+allsetcoll3

"""
from the 100 runs, the results of the last step from each of the three 
precalculation attempts is collected. The best 100 (lowest cost) of the 300
collected results are then selected to continue the fit (step 4), with 10
documented steps of 5 undocumented fitting substeps.
"""

resultcoll=[]
allcostcoll=[]
for i in allsetcoll:
    allcostcoll.append([i[j].cost for j in [2,5,8]])
    
allcostcoll=flatten(allcostcoll)
sortedcosts=np.argsort(allcostcoll)

for sel in np.arange(100):
    pos2=sortedcosts[sel]%3*3+2
    pos1=int(np.floor(sortedcosts[sel]/3))
    resultcoll.append(allsetcoll[pos1][pos2])
namresults='multix13_NEW_pfr_STU_'
    
ParameterColl2=mainfuncts.resc2param(PropAxesColl,ParameterColl,resultcoll,2)
conditions=[0,[1,1,1,10,5]]
setparameters3=['combo10l.dat',path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults,0]

mainfuncts.parallelmultifit4(setparameters3,5,20,ParameterColl2,PropAxesColl)
#%%
"""Stage 5: extended fitting procedure / minimization to account for all 16 
different possible combinations of ambiguous dwb and dwc"""

""
namresultsprev='multix13_NEW_pfr_STU_'

""" read in results from stage 4"""
conditions=[0,[1,1,1,10,5]]
PropAxesColl,ParameterColl2, costcoll, resultcoll,cond=readoutresults(reslall,resnall,pickthese,DeltaOmegaParametersBoundaries,namresultsprev,conditions)
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


"""The parameter set with the 16 copies of the most suited result from stage 4
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
    
conditions=[0,[1,1,1,5,5]]
setparameters3=['combo10l.dat',path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,'multix13_NEW_pfr_STU_errb_',0]
mainfuncts.parallelmultifit4(setparameters3,1,16,ParameterColl6,PropAxesColl)


#%%
""" Stage 6 - error calculation. In this final stage, errors are obtained by
resampling. For each of the 16 outputs from the previous calculation, 5
parameter sets from randomly resampled data are calculated. This gives a
total of 80 paarametersets to then report the error by calculating the standard
deviation for each parameter """
dataname=cond[-1][-1]
setparameters2=[dataname,path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,'multix13_NEW_pfr_STU_errb_'] #,'A34' #'/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/
for selset in np.arange(16):
    ss=hkfit2.evaluaterdfit(PropAxesColl,setparameters2,1,ParameterColl6[selset],0)
    hkio2.savss(ss,'multix13_NEW_pfr_STU_err2b_'+str(selset))

conditions=[0,[1,1,1,5,5]]
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
        setparameters3=['combo10l.dat',path2020,[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,'multix13_NEW_pfr_STU_err3b_','multix13_NEW_pfr_STU_err2b_'+str(runnum)]
        p=multiprocessing.Process(target=hkfit2.parallelfit3, args=(setparameters3,outnum*numrep+runnum,ParameterColl6[runnum],PropAxesColl))
        processes.append(p)
        processes[pn].start()
        pn+=1
    pn=list([pn0])[0]
    for runnum in np.arange(numrep):
        processes[pn].join()
        pn+=1
    print 'really done with ', str(outnum*numrep+numrep), 'runs'