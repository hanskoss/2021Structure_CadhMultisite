# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 16:41:51 2020

@author: Hans
"""
import hk_fitengines1_1 as hkfit2
import preprocess as prepro
import procall2
from __future__ import division
import hk_fitengines1_1 as hkfit2
import multiprocessing
from hkimports2 import flatten
from hkimports2 import np
from hkimports2 import optimize
from hkimports2 import time
from hkimports2 import os
import preprocess as prepro
import hkio
import hk_mainfuncts as mainfuncts

def produceindivplots(reslallx,resnallx,parbdslistaxallfx,namout,mode,k231,k232):
    for pickthis in np.arange(len(reslallx)):

        reslall=[reslallx[pickthis]]
        resnall=[resnallx[pickthis]]
        pickthese=[0]
        #global praxs1
        praxs1=mainfuncts.praxismake([x for y,x in enumerate(reslall) if y in pickthese],[x for y,x in enumerate(resnall) if y in pickthese])

        parbdslistaxallf=[parbdslistaxallfx[pickthis]]
    
    
        namresults='multix13_k13_ext_0'#######'multindiv_'+str(resnall[0])
        paramsxx=mainfuncts.parammake(praxs1,0,[x for y,x in enumerate(parbdslistaxallf) if y in pickthese],1,10000,100,900,k231,k232)
        conditions=[0,[1,3,3,3,3]]
        #setparameters2=['combo10l.dat','C:\\Users\\Hans\\Desktop\\TRANSFER\\2020Feb\\',[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults]
        poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr, intmin,ac,oc,rateconstpre,cond=hkio.loadeverything([namresults],0,decoupl=0)
        costcoll=[i.cost for i in resultcoll]
        resultcoll=[resultcoll[np.argsort(costcoll)[0]]]
        paramsxy=mainfuncts.resc2param(praxs1,paramsxx,resultcoll,mode)
        for nn in [0]:#np.arange(3):
            namresultsOUT=namout+str(resnall[0])+'_'+str(nn)+'_'
            setparameters3=['combo10l.dat','C:\\Users\\Hans\\Desktop\\TRANSFER\\2020Feb\\',[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresultsOUT,0]
            hkfit2.parallelfit3(praxs1,setparameters3,0,paramsxy[nn])
            


#%%
reslall=['A29','A30','A31','A32','A33','A34','A52','A53','A54','A81','A82','A83','A84']
resnall=[29,30,31,32,33,34,52,53,54,81,82,83,84]
reslall=['A78','A86']
resnall=[78,86]
reslall=['A86']
resnall=[86]
parbdslistaxallf=[[[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1,-1,],[-26000,26000],[-13000, 13000]],\
       #out    [[-1, 1], [-389.0, -5341], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]]]
parbdslistaxallf=[[[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]]]
parbdslistaxallf=[[[-1, 1], [-26000,26000], [-13000, 13000]]]
namout='indiv_free2_'
produceindivplots(reslall,resnall,parbdslistaxallf,namout,5,1,5)
#%%
reslall=['A29','A30','A31','A32','A33','A34','A52','A53','A54','A81','A82','A83','A84']
resnall=[29,30,31,32,33,34,52,53,54,81,82,83,84]
reslall=['A78','A86']
resnall=[78,86]
reslall=['A86']
resnall=[86]

parbdslistaxallf=[[[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1,-1,],[-26000,26000],[-13000, 13000]],\
       #out    [[-1, 1], [-389.0, -5341], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]]]
parbdslistaxallf=[[[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]]]
namout='indiv_free_'
produceindivplots(reslall,resnall,parbdslistaxallf,namout,5,1,10000)

#%%
reslall=['A29','A30','A31','A32','A33','A34','A52','A53','A54','A81','A82','A83','A84']
resnall=[29,30,31,32,33,34,52,53,54,81,82,83,84]
reslall=['A78','A86']
resnall=[78,86]


parbdslistaxallf=[[[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1,-1,],[-26000,26000],[-13000, 13000]],\
       #out    [[-1, 1], [-389.0, -5341], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]]]
parbdslistaxallf=[[[-1, 1], [-26000,26000], [-13000, 13000]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]]]
namout='indiv_fix_1_'
produceindivplots(reslall,resnall,parbdslistaxallf,namout,4,1,5)
#%%
reslall=['A29','A30','A31','A32','A33','A34','A52','A53','A54','A81','A82','A83','A84']
resnall=[29,30,31,32,33,34,52,53,54,81,82,83,84]
reslall=['A78','A86']
resnall=[78,86]
parbdslistaxallf=[[[-1, 1], [3199.0, 13980], [1440.0, 1640.0]],\
           [[-1, 1], [-22500.0, -6856.0], [-6307.0, -6107.0]],\
           [[-1, 1], [-1896.0, -136.0], [-736.0, -536.0]],\
           [[-1, 1], [2641.0, 7887.0], [1569.0, 1769.0]],\
           [[-1, 1], [-389.0, -123.0], [-38.0, 162.0]],\
           [[-1,-1,],[-269,868],[-386,-186]],\
       #out    [[-1, 1], [-389.0, -5341], [-13000, 13000]],\
           [[-1, 1], [-5341, -588.0], [-518.0, -318.0]],\
           [[-1, 1], [-26000,26000], [-2151,-1951]],\
           [[-1, 1], [-2465.0, 5540.0], [-3519.0, -3319.0]],\
           [[-1, 1], [-1346.0, 15376.0], [-1116.0, -916.0]],\
           [[-1, 1], [-772.0, 22607.0], [-1122.0, -922.0]],\
           [[-1, 1], [-1371.0, 8777.0], [-1287.0, -1087.0]],\
           [[-1, 1], [655.0, 5131.0], [405.0, 605.0]]]
parbdslistaxallf=[[[-1, 1], [7279,12491], [7151,7351]],\
           [[-1, 1], [-26000,26000], [-13000, 13000]]]
namout='indiv_fix_0_'
produceindivplots(reslall,resnall,parbdslistaxallf,namout,4,1,5)

#%%
reslall=['A81','A82','A83','A84'] #A52 is done
resnall=[81,82,83,84]

#parbdslistaxallf=[[[-1, 1], [-589, -587.0], [-518.0, -318.0]],\
parbdslistaxallf=[[[-1, 1], [-1347, -1345], [-1116.0, -916.0]],\
           [[-1, 1], [-773.0, -771], [-1122.0, -922.0]],\
           [[-1, 1], [-1372.0, -1370], [-1287.0, -1087.0]],\
           [[-1, 1], [654.0, 656], [405.0, 605.0]]]
namout='indiv_fix_x_'
produceindivplots(reslall,resnall,parbdslistaxallf,namout,4,1,5)


#%%

def readoutresults(namresultsx):
    reslall=['A29','A30','A31','A32','A33','A34','A52','A53','A54','A78','A81','A82','A83','A84','A86']
    resnall=[29,30,31,32,33,34,52,53,54,78,81,82,83,84,86]
    freevsfix='free' #'fix_0
    pickthis=1
    list1=[]
    for pickthis in np.arange(15):
        list2=[]
        for nn in np.arange(3):
            reslall=['A29','A30','A31','A32','A33','A34','A52','A53','A54','A78','A81','A82','A83','A84','A86']
            resnall=[29,30,31,32,33,34,52,53,54,78,81,82,83,84,86]
            reslall=[reslall[pickthis]]
            resnall=[resnall[pickthis]]
            pickthese=[0]
            
            praxs1=mainfuncts.praxismake([x for y,x in enumerate(reslall) if y in pickthese],[x for y,x in enumerate(resnall) if y in pickthese])
            parbdslistaxallf=[[[-1, 1], [-26000,26000], [-13000, 13000]],\
                       [[-1, 1], [-26000,26000], [-13000, 13000]],\
                       [[-1, 1], [-26000,26000], [-13000, 13000]],\
                       [[-1, 1], [-26000,26000], [-13000, 13000]],\
                       [[-1, 1], [-26000,26000], [-13000, 13000]],\
                       [[-1,-1,],[-26000,26000],[-13000, 13000]],\
                   #out    [[-1, 1], [-389.0, -5341], [-13000, 13000]],\
                       [[-1, 1], [-26000,26000], [-13000, 13000]],\
                       [[-1, 1], [-26000,26000], [-13000, 13000]],\
                       [[-1, 1], [-26000,26000], [-13000, 13000]],\
                       [[-1, 1], [-26000,26000], [-13000, 13000]],\
                       [[-1, 1], [-26000,26000], [-13000, 13000]],\
                       [[-1, 1], [-26000,26000], [-13000, 13000]],\
                       [[-1, 1], [-26000,26000], [-13000, 13000]],\
                       [[-1, 1], [-26000,26000], [-13000, 13000]],\
                       [[-1, 1], [-26000,26000], [-13000, 13000]]]
            parbdslistaxallf=[parbdslistaxallf[pickthis]]
            
            
            paramsxx=mainfuncts.parammake(praxs1,0,[x for y,x in enumerate(parbdslistaxallf) if y in pickthese],1,10000,100,900,1,10000)
            conditions=[0,[1,3,3,3,3]]
            #setparameters2=['combo10l.dat','C:\\Users\\Hans\\Desktop\\TRANSFER\\2020Feb\\',[[x for y,x in enumerate(reslall) if y in pickthese]],conditions,namresults]
            namresults=namresultsx+str(resnall[0])+'_'+str(nn)+'_'
            
            poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr, intmin,ac,oc,rateconstpre,cond=hkio.loadeverything([namresults],0,decoupl=0)
            costcoll=[i.cost for i in resultcoll]
            resultcoll=[resultcoll[np.argsort(costcoll)[0]]]
            paramsxy=mainfuncts.resc2param(praxs1,paramsxx,resultcoll,5)
            list2.append(costcoll[0])
        list1.append(list2)
    return list1

listx=[]
listx.append(readoutresults('indiv_free_'))
listx.append(readoutresults('indiv_free2_'))
listx.append(readoutresults('indiv_fix_1_'))
listx.append(readoutresults('indiv_fix_0_'))

#%%
dwblist=[]
dwclist=[]
for pickthis in np.arange(15):
    reslall=['A29','A30','A31','A32','A33','A34','A52','A53','A54','A78','A81','A82','A83','A84','A86']
    resnall=[29,30,31,32,33,34,52,53,54,78,81,82,83,84,86]
    reslall=[reslall[pickthis]]
    resnall=[resnall[pickthis]]
    pickthese=[0]
    
    praxs1=mainfuncts.praxismake([x for y,x in enumerate(reslall) if y in pickthese],[x for y,x in enumerate(resnall) if y in pickthese])
    parbdslistaxallf=[[[-1, 1], [-26000,26000], [-13000, 13000]],\
               [[-1, 1], [-26000,26000], [-13000, 13000]],\
               [[-1, 1], [-26000,26000], [-13000, 13000]],\
               [[-1, 1], [-26000,26000], [-13000, 13000]],\
               [[-1, 1], [-26000,26000], [-13000, 13000]],\
               [[-1,-1,],[-26000,26000],[-13000, 13000]],\
           #out    [[-1, 1], [-389.0, -5341], [-13000, 13000]],\
               [[-1, 1], [-26000,26000], [-13000, 13000]],\
               [[-1, 1], [-26000,26000], [-13000, 13000]],\
               [[-1, 1], [-26000,26000], [-13000, 13000]],\
               [[-1, 1], [-26000,26000], [-13000, 13000]],\
               [[-1, 1], [-26000,26000], [-13000, 13000]],\
               [[-1, 1], [-26000,26000], [-13000, 13000]],\
               [[-1, 1], [-26000,26000], [-13000, 13000]],\
               [[-1, 1], [-26000,26000], [-13000, 13000]],\
               [[-1, 1], [-26000,26000], [-13000, 13000]]]
    parbdslistaxallf=[parbdslistaxallf[pickthis]]
    paramsxx=mainfuncts.parammake(praxs1,0,[x for y,x in enumerate(parbdslistaxallf) if y in pickthese],1,10000,100,900,1,10000)
    resultcollx=[]
    costcollx=[]
    for nn in np.arange(3):
        namresults='indiv_fix_0_'+str(resnall[0])+'_'+str(nn)+'_'
        
        poscoll,resnam,allsetcoll,resultcoll,relaxrat0,relaxrat,lookatratio,results,relaxrat1,relaxrat2,relaxrat_1,relaxrat_2,intdiffcorr, intcorr, intmin,ac,oc,rateconstpre,cond=hkio.loadeverything([namresults],0,decoupl=0)
        costcoll=[i.cost for i in resultcoll]
        resultcollx.append(resultcoll[np.argsort(costcoll)[0]])
        costcollx.append(resultcoll[np.argsort(costcoll)[0]].cost)
    resultcoll=[resultcollx[np.argsort(costcoll)[0]]]
    paramsxy=mainfuncts.resc2param(praxs1,paramsxx,resultcoll,2)
    dwblist.append(paramsxy[0]['dw'][1]['par'][0])
    dwclist.append(paramsxy[0]['dw'][2]['par'][0])
reslall=['A29','A30','A31','A32','A33','A34','A52','A53','A54','A78','A81','A82','A83','A84','A86']
for j,i in enumerate(reslall):
    print i, int(np.round(dwblist[j],0)), int(np.round(dwclist[j],0))

#%%
from matplotlib import pyplot as plt
reslall=['A29','A30','A31','A32','A33','A34','A52','A53','A54','A78','A81','A82','A83','A84','A86']
fullnamlist=['Ser','Gly', 'Trp', 'Val', 'Trp', 'Asn', 'Gln', 'Phe', 'Phe', 'Val', 'Ile', 'Glu', 'Glu', 'Tyr', 'Thr', 'Gly', 'Pro', 'Asp', 'Pro', 'Val', 'Leu', 'Val', 'Gly', 'Arg', 'Leu', 'His', 'Ser', 'Asp', 'Ile', 'Asp', 'Ser', 'Gly', 'Asp', 'Gly', 'Asn', 'Ile', 'Lys', 'Tyr', 'Ile', 'Leu', 'Ser', 'Gly', 'Glu', 'Gly', 'Ala', 'Gly', 'Thr', 'Ile', 'Phe', 'Val', 'Ile', 'Asp', 'Asp', 'Lys', 'Ser', 'Gly', 'Asn', 'Ile', 'His', 'Ala', 'Thr', 'Lys', 'Thr', 'Leu', 'Asp', 'Arg', 'Glu', 'Glu', 'Arg', 'Ala', 'Gln', 'Tyr', 'Thr', 'Leu', 'Met', 'Ala', 'Gln', 'Ala', 'Val', 'Asp', 'Arg', 'Asp', 'Thr', 'Asn', 'Arg', 'Pro', 'Leu', 'Glu', 'Pro', 'Pro', 'Ser', 'Glu', 'Phe', 'Ile', 'Val', 'Lys', 'Val', 'Gln', 'Asp']
plt.rcParams.update({'font.size': 18})
lab=[];y1=[];y2=[];y3=[];y4=[]
for j,i in enumerate(listx[0]):
    lab.append(reslall[j])
    y1.append(np.min(listx[0][j]))
    y2.append(np.min(listx[1][j]))
    y3.append(np.min(listx[2][j]))
    y4.append(np.min(listx[3][j]))
plt.figure(700)
plt.clf()
plt.bar(np.arange(len(lab))-0.2,y1,width=0.15,label='free fit; '+r'$\overline{\chi^2}$'+' = '+str(np.round(np.average(y1),2)),color='blue')
plt.bar(np.arange(len(lab))-0.2+0.35/2,y2,width=0.15,label='free fit, but k'+r'$_2$$_3$'+'=0; '+r'$\overline{\chi^2}$'+' = '+str(np.round(np.average(y2),2)),color='orange')
plt.bar(np.arange(len(lab))+0.35-0.2,y3,width=0.15,label='restrained k, p, but '+r'$\Delta\omega$'+' free; '+r'$\overline{\chi^2}$'+' = '+str(np.round(np.average(y3),2)),color='cyan')
plt.bar(np.arange(len(lab))+0.35-0.2+0.35/2,y4,width=0.15,label='restrained k, p, '+r'$\Delta\omega$; '+r'$\overline{\chi^2}$'+' = '+str(np.round(np.average(y4),2)),color='red')
ind=np.arange(len(lab))
plt.xticks(ind,tuple(['','','','']))
plt.xticks(ind,tuple([fullnamlist[int(xx[1]+xx[2])]+xx[1]+xx[2] for xx in reslall]))
plt.setp(plt.gca().get_xticklabels(), rotation=70)
plt.legend()
plt.ylabel(r'$\chi^2$')
plt.tight_layout()
