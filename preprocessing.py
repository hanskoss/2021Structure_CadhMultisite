#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 01:23:31 2020

@author: hanskoss
"""
from __future__ import division

import numpy as np

from sympy import flatten
import sympy as smp

import re
import csv
import os
import os.path

class spinsystem:
    def __init__(self):
        self.name = []
        self.ssnumlist = []
        self.orig = []
        self.orifiletype = []
        self.duplstat=0
        self.datasets= []
        self.nucleus = []
        self.seqno = 0

def addsystem(spinsystems,spinlistpath,filetype):
    #whichres=[]
    os.chdir(spinlistpath)
    relpath='./'
    if filetype ==4: #peaklist format
        with open(relpath+'peaks.reslist') as csvfile:
            filetext = csv.reader(csvfile,delimiter=' ')
            for row in filetext:
                name=row[0]
                if len([j for j,i in enumerate(spinsystems) if name in i.name]) == 0:
                    spinsystems.append(spinsystem())
                    spinsystems[-1].name.append(name)
                    namtoint=name
                    #print namtoint, spinlistpath
                    letnumtoint=0
                    for letnum, letter in enumerate(['A','B','C','D','E','F']):
                        if letter in namtoint:
                            namtoint=re.sub(letter,'',namtoint)
                            letnumtoint+=1000*(letnum+1)
                    letnumtoint+=int(namtoint)
                    spinsystems[-1].ssnumlist.append(letnumtoint)
                    spinsystems[-1].orig=spinlistpath
                    spinsystems[-1].orifiletype=filetype
                    spinsystems[-1].seqno=int(namtoint)
    else: #reslist.csv, simple residue list as created by CEST data processing
        with open(relpath+'reslist.csv') as csvfile:
            filetext = csv.reader(csvfile,delimiter=',')
            for row in filetext:
                name=row[0]
                try:
                    offset=row[1]
                except:
                    offset=0
                if len([j for j,i in enumerate(spinsystems) if name in i.name]) == 0:
                    spinsystems.append(spinsystem())
                    spinsystems[-1].name.append(name)
                    spinsystems[-1].offset=offset
                    namtoint=name
                    letnumtoint=0
                    for letnum, letter in enumerate(['A','B','C','D','E','F']):
                        if letter in namtoint:
                            namtoint=re.sub(letter,'',namtoint)
                            letnumtoint+=1000*(letnum+1)
                            spinsystems[-1].letter=letter
                    letnumtoint+=int(namtoint)
                    spinsystems[-1].ssnumlist.append(letnumtoint)
                    spinsystems[-1].orig=spinlistpath
                    spinsystems[-1].orifiletype=filetype
                    spinsystems[-1].seqno=int(namtoint)
                    #pass
                    #print 'new ss', name
    return spinsystems

def precalcerr(pathlist,errpos,actualpos,merge,slt):
    intdat=[];hshift=[];nshift=[];sloffr=[];subexptype=[];slhard=[];offsSL=[]
    slhard=[];sllen=[];nam=[];group=[];offs=[];firstrr=1;resl=[]
    offse=[];intdat=[];nshiftx=[];b0field=[]
    
    with open('resall_1r.dat','rb') as csvfile:
        filetext = csv.reader(csvfile, delimiter=' ')
        for row in filetext:
            resl.append(row[0])
            intdat.append(row[1::]) #-1 instead of second :
    ctr=0
    hshiftl=[]
    nshiftl=[]
    with open('res_1r.dat','rb') as csvfile:
        filetext = csv.reader(csvfile, delimiter=' ')
        for row in filetext:
            if ctr > 0:
                hshiftl.append(row[1])
                nshiftl.append(row[2])
                if firstrr == 1:
                    hshift.append(row[1])
                    nshift.append(row[2])
                    
                nshiftx.append(row[2]) 
            ctr+=1
    ctr=0
    with open('res_1.dat.headerdat2','rb') as csvfile:
        filetext = csv.reader(csvfile, delimiter=' ')
        for row in filetext:
            if ctr > 0:
                offse.append(row[2])
                nam.append(row[0])
                group.append(row[0].split('@'))
                offsSL.append(row[1])
                offs.append(row[2])
                slhard.append(row[3])
                sllen.append(row[4])
                subexptype.append(row[5])
                b0field.append(row[6])
                sloffr.append(row[7])
            ctr+=1
    nshift=np.array(nshift).astype('float')
    offsSL=np.array(offsSL).astype('float')
    sloffr=np.array(sloffr).astype('float')
    b0field=np.array(b0field).astype('float')

    i=0
    resoecoll=[];resangcoll=[];reserrfac=[];resomegacoll=[]
    intwitherrcoll=[];intonlycoll=[]
    t1=[];t2=[];tr=[];hsl=[];t1err=[];t2err=[];t2errcollx=[];trerrcollx=[]
    trerr=[];hslerr=[];t2errn=[];hsllen=[]
    for res in resl:
        j=0
        resoe=[];resang=[];resomega=[];intwitherr=[];intonly=[]
        errfac=[];errfac2=[];t1errcoll=[];t2errcoll=[];trerrcoll=[];hslerrcoll=[]
        
        k=0
        for sepos in subexptype:
            if sepos == 'HH_offr_or':
                intwitherr.append(intdat[i][j])
                k+=1
            j+=1
        intwitherrcoll.append(intwitherr)
        j=0
        k=0
        cnthard=0
        for sepos in subexptype:
            if sepos == 'HH_offr_or':
                bigomega=np.round((offsSL[j]-nshift[i])*b0field[j])
                omegae=np.sqrt(bigomega**2+(1000*sloffr[j])**2)  ### effective field
                flipang=np.arctan((1000*sloffr[j])/bigomega) ###57.2958*
                
                if k in actualpos:
                    resoe.append(omegae)
                    resang.append(flipang)
                    resomega.append(bigomega)
                    f=np.array(smp.flatten([[intwitherr[c] for d,c in enumerate(e)] for e in merge if k in e])).astype(int)
                    if len(f) == 0:
                        f=np.array([intwitherr[k]]).astype(int)
                    errfac.append(1/np.sqrt(len(f)))
                    intonly.append(np.average(f))
                k+=1
                errfac2.append([1])
            if sepos == 'HH_offr_t1r':
                t1errcoll.append(intdat[i][j])
            if sepos == 'HH_offr_t2r':
                t2errcoll.append(intdat[i][j])
            if sepos == 'HH_offr_hslr':
                if np.abs((offsSL[j]-nshift[i])) < 4: #### OFFSET cutoff for HARD spinlock
                    hslerrcoll.append(intdat[i][j])
                    cnthard+=1
            if sepos == 'HH_offr_trr':
                trerrcoll.append(intdat[i][j])
            if sepos == 'HH_offr_t1':
                t1x=float(intdat[i][j])
            if sepos == 'HH_offr_t2':
                t2x=float(intdat[i][j])
            if sepos == 'HH_offr_hsl':
                hslx=float(intdat[i][j])
            if sepos == 'HH_offr_tr':
                trx=float(intdat[i][j])
            j+=1
        i+=1
        #print hslerrcoll, 'hslerrcoll'
        #print t2x, hslx
        hslx=np.average(np.array(hslerrcoll).astype('int')) ###if average values are obtained from individual
   #     print hslx, 'hslx', res
        t2x=np.average(np.array(t2errcoll).astype('int'))
        t1x=np.average(np.array(t1errcoll).astype('int'))
        trx=np.average(np.array(trerrcoll).astype('int'))
        intonlycoll.append(intonly)
        t1err.append(np.average((np.array(t1errcoll).astype(float)-t1x)**2)/(len(t1errcoll)-1))
        t2err.append(np.average((np.array(t2errcoll).astype(float)-t2x)**2)/(len(t2errcoll)-1))
        t2errn.append(len(t2errcoll))
        t2errcollx.append(t2errcoll)
        trerrcollx.append(trerrcoll)
        trerr.append(np.average((np.array(trerrcoll).astype(float)-trx)**2)/(len(trerrcoll)-1))
       # print len(hslerrcoll)
        if len(hslerrcoll) > 1:
            hslerr.append(np.sum((np.array(hslerrcoll).astype(float)-hslx)**2)/(len(hslerrcoll)-1))
        hsllen.append(len(hslerrcoll))       
        reserrfac.append(errfac);resoecoll.append(resoe)
        resangcoll.append(resang);resomegacoll.append(resomega)
        t1.append(t1x);t2.append(t2x);tr.append(trx);hsl.append(hslx)
    if os.path.exists("errset") == 1:
        errext=[]
        with open('errset','rb') as csvfile:
            filetext = csv.reader(csvfile, delimiter=' ')
            for row in filetext:
                errext.append(int(row[0]))
        trerr=np.sqrt(errext[0])
        hslerr=np.sqrt(errext[1])
        t1err=np.sqrt(errext[2])
        t2err=np.sqrt(errext[3])
    else:
        t1err=np.sqrt(np.average(t1err))
        t2err=np.sqrt(np.average(t2err))
        trerr=np.sqrt(np.average(trerr))
        hslerr=np.sqrt(np.average(hslerr))

    t2errall=t2err
    intwitherrcoll=np.array(intwitherrcoll).astype('int')
    intonlycoll=np.array(intonlycoll).astype('int')
    #DO NOT DELETE. preferred error estimation########errestimate=np.sqrt(np.average(np.array(smp.flatten([[np.sqrt(np.sum((m-np.average(m))**2)/(len(m)-1)) for m in n] for n in [[k for b,k in enumerate(intwitherrcoll[:,l])] for l in [smp.flatten([x for x in merge if j in x]) for j in errpos]]]))**2))
    ##reserrfac=1
    errestimate=hslerr#t2err #np.max([t2err,hslerr,trerr,hslerr])
    reserrfac=errestimate*np.array(reserrfac)
    resangcoll=np.array(resangcoll).astype('float')
    resoecoll=np.array(resoecoll).astype('float')
    resomegacoll=np.array(resomegacoll).astype('float')
    return [nshift,hsllen,resl,hsl,t2errn,t2errall,t2err,trerr,t1err,reserrfac,resangcoll,hslerr,t2,t1,intonlycoll,tr,resoecoll, b0field,resomegacoll]
    
def processrelaxdat(pathlist,errpos,actualpos,merge,posy,slt,precalcstuff):
    nshift,hsllen,resl,hsl,t2errn,t2errall,t2err,trerr,t1err,reserrfac,resangcoll,hslerr,t2,t1,intonlycoll,tr,resoecoll, b0field,resomegacoll=precalcstuff
 #   print hsl, 'hsl' 
    #print hsllen[12]
    
    numbtest=10000
    #slt=70000
    if posy in resl:
        posx=resl.index(posy)
        
        r1rhardmat=np.array(np.random.normal(hsl[posx],hslerr,[numbtest]))
        
       # if t2errn[posx] > 1000: #np.std(np.array(t2errcollx[posx]).astype('int'))
       #     r2mat=np.array(np.random.normal(t2[posx],np.sqrt(t2errall[posx]),[numbtest]))
       # else:
        r2mat=np.array(np.random.normal(t2[posx],t2err,[numbtest]))
        trmat=np.array(np.random.normal(tr[posx],trerr,[numbtest]))
     #   print r2mat, r1rhardmat, 'r2', 'rhard'
        try:
            r1mat=np.array(np.random.normal(t1[posx],t1err,[numbtest]))
            r1roffmat=np.array(np.random.normal(intonlycoll[posx],reserrfac[posx],[numbtest,len(intonlycoll[posx])]))
        
            sinfcx=np.sin(resangcoll[posx])**2
            cosfcx=np.cos(resangcoll[posx])**2
            
            datcomb=[((r1roffmat[z])/((r1rhardmat[z]**sinfcx)*(r1mat[z]**cosfcx))) for z in np.arange(numbtest)]
            ########datcomb=[(r1roffmat[z]/trmat[z]) for z in np.arange(numbtest)]
            #datcomb=[(r1roffmat[z]) for z in np.arange(numbtest)]
            r1rmatresf1=np.array([[-np.log(y)*(1000000/slt) if y > 0 else 0 for y in x] for x in datcomb])
            #r1rmatresf1=np.array([[y if y > 0 else 0 for y in x] for x in datcomb])
            r1rmatresf=r1rmatresf1
            #r1rmatresf=r1rmatresf1/sinfcx
            #####datcomb=[(r1roffmat[z])/(trmat[z]) for z in np.arange(numbtest)]
            r1rav=np.average(r1rmatresf,axis=0) #
            confid=655
            selarr=np.array([list(z)[int((1000-confid)*numbtest/2000):int(numbtest-(1000-confid)*numbtest/2000)] for z in np.array([np.sort(x) for x in np.transpose(r1rmatresf)])])
            err1=r1rav-np.min(selarr,axis=1)
            err2=np.max(selarr,axis=1)-r1rav
        except:
            r1rav=0
            err1=0
            err2=0
            
        datcomb=[((r2mat[z]/r1rhardmat[z])) for z in np.arange(numbtest)]
        
        r1rmatresf=np.array([-np.log(y)*(1000000/slt) if y > 0 else 0 for y in datcomb])
        #r1rav2=np.average(r1rmatresf)
        confid=655
        selarr=np.array([z for z in np.array(np.sort(r1rmatresf))][int((1000-confid)*numbtest/2000):int(numbtest-(1000-confid)*numbtest/2000)])
   #     print selarr
        #try:
        
        r1rav2=-np.log(t2[posx]/hsl[posx])*(1000000/slt)
   #     print t2[posx],hsl[posx],t2[posx]/hsl[posx], r1rav2, 'rexjustcalct', slt
        #    r1rav2=np.average(selarr)#average(selarr)
        
        errcomb=-(t2[posx]/hsl[posx])*np.sqrt(np.average((t2err/t2[posx])**2+(hslerr/hsl[posx])**2))
        errmin=np.min(selarr)
        errmax=np.max(selarr)
        err12=r1rav2-np.min(selarr)
        err22=np.max(selarr)-r1rav2
    #    errmin=-np.log(t2[posx]/hsl[posx]-errcomb)*(1000000/slt)
   #     errmax=-np.log(t2[posx]/hsl[posx]+errcomb)*(1000000/slt)
    #    erralt=[errmax-r1rav2,r1rav2-errmin]
        #np.average([-np.log(t2[posx]/hsl[posx])*(1000000/slt)+np.log(t2[posx]/hsl[posx]+errcomb)*(1000000/slt),-np.log(t2[posx]/hsl[posx])*(1000000/slt)-np.log(t2[posx]/hsl[posx]-errcomb)*(1000000/slt)])
      #  print -np.log(t2[posx]/hsl[posx])*(1000000/slt), r1rav2, err12, err22, erralt, errmin, errmax, errcomb
     #   print err12,err22,r1rav2
        datcomb=[((r1rhardmat[z]/trmat[z])) for z in np.arange(numbtest)]
        r1rmatresf=np.array([-np.log(y)*(1000000/slt) if y > 0 else 0 for y in datcomb])
        r1rav4=np.average(r1rmatresf)
        confid=655
        selarr=np.array([z for z in np.array(np.sort(r1rmatresf))][int((1000-confid)*numbtest/2000):int(numbtest-(1000-confid)*numbtest/2000)])
        err14=r1rav4-np.min(selarr)
        err24=np.max(selarr)-r1rav4
    
    
        datcomb=[((r2mat[z])/trmat[z]) for z in np.arange(numbtest)]
        r1rmatresf=np.array([-np.log(y)*(1000000/slt) if y > 0 else 0 for y in datcomb])
        r1rav3=np.average(r1rmatresf)
        confid=655
        selarr=np.array([z for z in np.array(np.sort(r1rmatresf))][int((1000-confid)*numbtest/2000):int(numbtest-(1000-confid)*numbtest/2000)])
        err13=r1rav3-np.min(selarr)
        err23=np.max(selarr)-r1rav3
    else:
        err13=0;err23=0;err24=0;err14=0;err12=0;err22=0;err1=0;err2=0;r1rav=0;r1rav2=0;r1rav4=0;r1rav3=0
    return resl, resomegacoll, resangcoll, resoecoll, b0field, nshift, r1rav, err1, err2, r1rav2,err12,err22,r1rav4,err14,err24,r1rav3,err13,err23

                        
def adddata(spinsystems,spinlistpath,filetype,dataposition,setselect,setlabel,pathlist):
    """
    This function reads in data and connects them to the individual spin systems
    input variables:
    spinsystems: spinsystems (still without data sets)
    spinlistpath: the datasets are physically stored in numbered folders. This
    is the folder of the dataset.
    filetype: number definining the type of data set which are being read in.
        See below.
    dataposition:
        
    setselect: The list of data sets being read in was supplied elsewhere,this
    stores the position of the dataset in that list
    setlabel: Each dataset has a unique label
    pathlist: useless variable.
        
    The settings.dat files is read in and all lines are saved in the expcond
    variables, which is later added to any of the data sets.
    The settings.dat file contains the following information.
    
    The next step depends on the file type which is being read in (variable fileypte)
    Some details about file types here.
    5: peak lists. mostly uses the addshifts function for dataset dict object
        see details for shiftlisb.csv file there
    4: multi- titration data (with pseudodimension).
        files which are being read in:
        peaks.properties
            first row, column2: type of pseudodimension, for example pressure
            second row, column2: temperature
        peaks.pseudodim
            first column: numbers counting up, referring to pseudodimension
            second column: pseudodimension value, z.B. pressure in psi
        peaks.peaklist: columns:
            1: spin system name
            2: number of pseudodimension position
            3: proton shift
            4: nitrogen shift
            5: heights
            6: volumes
            7: peak width 1H dimension
            8: peak width 15N dimension
    3: R1rho-experiment, including Rex experiment
            r1rsetting.dat. at this moment it contains the spin lock duration
                in line 2, in us.

    res_1.dat.headerdat2 contains line-by-line information about the different
    sub-experiments in the pseudo-3D-r1rho cube.
    Each line contains the following information (column by column):
    name@pseudodimposition@1
    spin lock position (15N, in ppm)
    15N offset (center of spectrum)
    spinlock power (high power experiment) in kHz
    spinlock duration in us
    type of subexperiment:
        HH_offr_or off-resonance (low power) experiment
        HH_offr_t1r T1 experiment
        HH_offr_t2r T2 experiment
        HH_offr_trr T=0 (regular correlation experiment)
        HH_offr_hslr on-resonance (high power) experiment

    4: CPMG experiments

    cpmgsetting.dat contains the following information (line by line):
    row1: space-separated list of positions of TROSY experiments (T=0) in pseudo-
    3D list
    row2: pi pulse length (in relaxation dispersion dimension)
    row3: field (B0) in MHz
    row4: nitrotag, 2 for nitrogen (probably not needed anymore)
    row5: number of pi pulses per CPMG loops (2 loops with 2 pulses means 4)
    
    vclist: number of CPMG loop repetitions for each subexperiment
    vdlist: list of tau_cps - half of the period between the pi pulses
    
    CPMGpeaks.dat: starting from line two, eah column contains the following information:
    column1: protein type (A: folded; B: unfolded; C: dimer; D: other)
    column2: sequence number of residue
    column3: pseudodimension position
    column4: peak height
        
        
        
        
    """
    print spinlistpath,filetype
    os.chdir(spinlistpath)
    relpath='./'
    expcond={}
    try:
        with open('settings.dat','rb') as csvfile:
            filetext = csv.reader(csvfile, delimiter=' ')
            for row in filetext:
                expcond[row[0]] = row[1]
    except:
        print 'no experimental conditions defined'
    if filetype == 5:
        for nf,l in enumerate([k.name[0] for k in spinsystems]):
            spinsystems[nf].datasets.append(RDset(spinlistpath,filetype,dataposition,l,'peaklist',expcond['nucleus'],setselect))
            spinsystems[nf].datasets[-1].addexpcond(expcond)
            field=float(expcond['B0'])/(500*9.868831346500839)
            spinsystems[nf].datasets[-1].addshifts(setlabel,spinsystems,expcond['TR'],field,nf,len(spinsystems[nf].datasets),filetype)

    elif filetype ==4:
        with open(relpath+'peaks.properties') as csvfile:
            filetext = csv.reader(csvfile,delimiter=' ')
            for rnum, row in enumerate(filetext):
                if rnum == 0:
                    psdimprop=row[1]
                elif rnum == 1:
                    temp=row[1]
        with open(relpath+'peaks.pseudodim') as csvfile:
            filetext = csv.reader(csvfile,delimiter=' ')
            psdim=[];psdimx=[]
            for rnum, row in enumerate(filetext):
                psdimx.append(row[0])
                psdim.append(row[1])
            
        for nf,l in enumerate([k.name[0] for k in spinsystems]):
            spinsystems[nf].datasets.append(titpeaklist(spinlistpath,filetype,dataposition,l,'titpeaklist',expcond['nucleus'],psdimprop,temp,psdimx,psdim,setselect))
            spinsystems[nf].datasets[-1].addexpcond(expcond)
            xdatx=[];xdat=[]; hshifts=[]; nshifts=[]; heights=[]; volumes=[]; pw1=[]; pw2=[]
            datth=0
            for xx, x in enumerate(psdimx):
                with open(relpath+'peaks.peaklist') as csvfile:
                    filetext = csv.reader(csvfile,delimiter=' ')
                    for row in filetext:
                        if row[0] == l and int(row[1]) == int(x):
                            datth=1
                            xdatx.append(x)
                            xdat.append(psdim[xx])
                            hshifts.append(row[2])
                            nshifts.append(row[3])
                            heights.append(row[4])
                            volumes.append(row[5])
                            pw1.append(row[6])
                            pw2.append(row[7])
            hshifts=np.array(hshifts).astype(float)
            nshifts=np.array(nshifts).astype(float)
            heights=np.array(heights).astype(int)
            volumes=np.array(volumes).astype(int)
            pw1=np.array(pw1).astype(float)
            pw2=np.array(pw2).astype(float)
            xdat=np.array(xdat).astype(int)
            xdatx=np.array(xdatx).astype(int)
            if datth ==1:
                if spinsystems[nf].datasets[-1].psdimprop == 'pressure':
                    diffpos=[]
                    shdiff=[]
                    shdirect=[]
                    shdirect2=[]
                    xdatargs=np.argsort(xdat)
                    xdatx=xdat[xdatargs]
                    field=float(expcond['B0'])/(500*9.868831346500839)
                    if spinsystems[nf].datasets[-1].expcond['TR'][0] == 'T' or spinsystems[nf].datasets[-1].expcond['TR'][0] == 'S':
                        trshiftn=-92/(2*field*500)
                        trshifth=92/(2*field*500*9.868831346500839)
                    else:
                        trshifth=0
                        trshiftn=0
                    hshiftsx=hshifts[xdatargs]
                    nshiftsx=nshifts[xdatargs]
                    hshifts=hshifts+trshifth-(float(spinsystems[nf].datasets[-1].expcond['temp'])-298)*0.0046
                    nshifts=nshifts+trshiftn
                    
                    for xx,x in enumerate(xdat):
                        if xx > 0:
                            diffpos.append((prev-x)/2+prev)
                            shdiff.append((np.sqrt(((hshiftsx[xx-1]-hshiftsx[xx])**2+(0.14*(nshiftsx[xx-1]-nshiftsx[xx]))**2))/2)/(x-prev))
                            ang1=(hshiftsx[xx-1]-hshiftsx[xx])
                            ang2=(0.14*(nshiftsx[xx-1]-nshiftsx[xx]))
                            shdirect.append(np.arctan2(ang1,ang2)*57.2958)
                        prev=x
                shdirect2=[i+360 if i<0 else i for i in shdirect]
                spinsystems[nf].datasets[-1].addtitdat(setlabel,xdat,xdatx,hshifts,nshifts,heights,volumes,pw1,pw2,diffpos,shdiff,shdirect,shdirect2)
                
                            
                        
        
    if filetype ==3:
        with open(relpath+'r1rsetting.dat') as csvfile:
            filetext = csv.reader(csvfile,delimiter=',')
            for rnum, row in enumerate(filetext):
                if rnum == 1:
                    slt=np.int(row[0]) #np.array(row).astype('int')
                elif rnum == 2:
                    errpos=list(np.array(row).astype('int'))
                elif rnum == 3:
                    actualpos=list(np.array(row).astype('int'))
                elif rnum == 4:
                    merge=[list(np.array(row).astype('int'))]
                elif rnum > 4:
                    merge.append(list(np.array(row).astype('int')))
        precalcerrstuff=precalcerr(pathlist,errpos,actualpos,merge,slt)
        for nf,l in enumerate([k.name[0] for k in spinsystems]):
#            try:
            
            resl, resomegacoll, resangcoll, resoecoll, b0field, nshift, r1rav, err1, err2,rex, errex1, errex2,r20,err201,err202,r2f,err211,err212=processrelaxdat(pathlist,errpos,actualpos,merge,l,slt,precalcerrstuff)
  #          print rex, 'rexrexrex' 
            if errex1 != 0 and errex2 != 0 and rex != 0:
            #if len(resangcoll) > 0:
               # print rex, 'rexhere', pathlist,errpos,actualpos,merge,l,slt,precalcerrstuff
                spinsystems[nf].datasets.append(RDset(spinlistpath,filetype,dataposition,l,'Rex',expcond['nucleus'],setselect))
                spinsystems[nf].datasets.append(RDset(spinlistpath,filetype,dataposition,l,'R2',expcond['nucleus'],setselect))
                spinsystems[nf].datasets.append(RDset(spinlistpath,filetype,dataposition,l,'R20',expcond['nucleus'],setselect))
                spinsystems[nf].datasets.append(RDset(spinlistpath,filetype,dataposition,l,'r1rho',expcond['nucleus'],setselect))
                for i in [-4,-3,-2,-1]:
                    spinsystems[nf].datasets[i].addexpcond(expcond)
                for i in [-4,-3,-2,-1]:
                    spinsystems[nf].datasets[i].addshifts(setlabel,spinsystems,expcond['TR'],float(expcond['B0'])/(500*9.868831346500839),nf,len(spinsystems[nf].datasets),filetype)
                    
                spinsystems[nf].datasets[-1].addrelax(setlabel,rex, errex1, errex2,r20,err201,err202,r2f,err211,err212)
                
                try:
                    if len(r1rav)>0:
                        ng=resl.index(spinsystems[nf].name[0])
                        spinsystems[nf].datasets[-1].addr1rho(setlabel,resomegacoll[ng], resangcoll[ng], resoecoll[ng], b0field[ng], nshift[ng], r1rav, err1, err2)
                except:
                    trash=0
                try:
                    if r1rav != 0 and err1 != 0 and err2 != 0:
                        ng=resl.index(spinsystems[nf].name[0])
                        spinsystems[nf].datasets[-1].addr1rho(setlabel,resomegacoll[ng], resangcoll[ng], resoecoll[ng], b0field[ng], nshift[ng], r1rav, err1, err2)
                except:
                    trash=0
                try:
                    titlfrag=expcond['B0']
                except:
                    titlfrag=''
                #print rex, errex1,errex2, 'add',spinsystems[nf].name
                spinsystems[nf].datasets[-4].adddatapoint('Rex','Rex ' + setlabel,'s-1',rex,errex1,errex2)
                spinsystems[nf].datasets[-3].adddatapoint('R2','R2 ' + setlabel,'s-1',r2f,err211,err212) #titlfrag
                spinsystems[nf].datasets[-2].adddatapoint('R20','R20 ' + setlabel,'s-1',r20,err201,err202) #titlfrag
            else:
                pass
                #print spinsystems[nf].name, 'not there'
            
            
#except:
#                    'no r1r there'    
#                spinsystems[nf].datasets[-1].addr1rho(resomegacoll[ng], resangcoll[ng], resoecoll[ng], b0field[ng], nshift[ng], r1rav, err1, err2,rex, errex1, errex2,r20,err201,err202,r2f,err211,err212)
 #           except:
            
        
    if filetype ==2: #standard nitrogen-CPMG data
        with open(relpath+'cpmgsetting.dat') as csvfile:
            filetext = csv.reader(csvfile,delimiter=' ')
            for rnum, row in enumerate(filetext):
                if rnum == 1:
                    #trpos=np.array(row).astype('int')
                    trpos=np.array(row[:][0:-1]).astype('int')
                elif rnum == 2:
                    p2=np.float(row[0])
                elif rnum == 3:
                    field=np.float(row[0])/500
                elif rnum == 4:
                    nitrotag=np.int(row[0])
                elif rnum == 5:
                    p4pulsenumber=np.int(row[0])
        vc=[]
        vd=[]
        with open('vclist','rb') as csvfile:
            filetext = csv.reader(csvfile, delimiter=' ')
            for row in filetext:
                vc.append(row[0])
        with open('vdlist','rb') as csvfile:
            filetext = csv.reader(csvfile, delimiter=' ')
            for row in filetext:
                vd.append(row[0])
        vc=np.array(vc).astype('int')
        vd=np.array(vd).astype('float')
        vcf=np.array([i[0] if j+1 not in trpos else 0 for j,i in enumerate(zip(vc,vd))])
        vdf=np.array([i[1] if j+1 not in trpos else 0 for j,i in enumerate(zip(vc,vd))])
        vcfs=np.sort(np.array(list(set(vcf))))
        vdfs=np.array([vdf[j] for j in [list(vcf).index(i) for i in vcfs]])
        fdata=1/(4*vdfs)
        decaytimes=vcfs*(p4pulsenumber*p2/1000000+p4pulsenumber*2*vdfs)
        label=[];intens=[];peakno=[]
        
        with open('CPMGpeaks.dat','rb') as csvfile:
            filetext = csv.reader(csvfile, delimiter=' ')
            for i,row in enumerate(filetext):
                if i>0:
                    try:
                        trash=int(row[0])
                        digis=1
                    except:
                        digis=0
                    if digis == 1:
                        label.append(row[0])
                    else:
                        label.append(row[0]+row[1])
                    intens.append(row[3])
                    peakno.append(row[2])
        intens0=np.array(intens).astype('int')
        peakno=np.array(peakno).astype('int')

        for nf,l in enumerate([k.name[0] for k in spinsystems]):
            intens=[j for i,j in enumerate(intens0) if label[i]==l]
            spinsystems[nf].datasets.append(RDset(spinlistpath,filetype,dataposition,l,'cpmg',expcond['nucleus'],setselect))
            spinsystems[nf].datasets[-1].addexpcond(expcond)
            spinsystems[nf].datasets[-1].addshifts(setlabel,spinsystems,expcond['TR'],field,nf,len(spinsystems[nf].datasets),filetype)
            if len(intens) > 0:
                errcnt=[]
                errget=[]
                dataint=[]
                trint=[]
                pntcnt=[]
                #now collect CPMG data, going through each future data point of the CPMG curve.
                for i in vcfs:  #leaves option open to create a supra-dataset vcfs from vclist and vdlist and then match whichever is in the actual data
                        #but careful. right now,  it assumes that intensity positions map to vcf positions!
                    if np.sum(i == vcf) > 1: #if more than one point
                        a=[intens[k] for k,j in enumerate(vcf) if j == i]
                        errget.append((a-np.average(a))**2)
                        errcnt.append(len(a))
                    a=[intens[k] for k,j in enumerate(vcf) if (j == i)] # and k not in trpos
                    if i == 0:
                        trint.append(np.average(a))
                        pntcntr=np.sum([1 for k,j in enumerate(vcf) if (j == i)])
                        
                    else:
                        dataint.append(np.average(a))
                        pntcnt.append(np.sum([1 for k,j in enumerate(vcf) if (j == i)]))
                dataint=flatten(dataint)
                trint=flatten(trint)[0]
                err=np.sqrt(np.average(flatten(errget)))
                #err=np.sqrt(np.average(flatten(errget)/np.sqrt(np.average(errcnt))))
                j=[dataint/trint]
                #print (err/(np.array(flatten(np.sqrt(np.array(pntcnt))))*np.array(dataint)))**2
                #print np.array((err/np.array(dataint))**2)
                #print np.array((err/np.array(dataint))**2)/np.array(flatten(np.sqrt(np.array(pntcnt))))
                #print '!',l, pntcnt
                j.append((dataint/trint)*np.sqrt((np.array(err/np.array(dataint))/np.array(flatten(np.sqrt(np.array(pntcnt)))))**2+(np.array((err/np.array(trint)))/np.sqrt(pntcntr))**2))
                #print (dataint/trint)*np.sqrt((np.array(err/np.array(dataint))/np.array(flatten(np.sqrt(np.array(pntcnt)))))**2+(np.array((err/np.array(trint)))/np.sqrt(pntcntr))**2)
                decaytimex=decaytimes[1:]
                ydata=([-np.log(j[0][i])*(1/decaytimex[i]) for i in np.arange(len(j[0]))])
                if np.min([(j[0][i]-j[1][i]) for i in np.arange(len(j[0]))]) < 0:
                    ydataerr=zip([-np.log(j[0][i])*(1/decaytimex[i])+np.log(j[0][i]+j[1][i])*(1/decaytimex[i]) for i in np.arange(len(j[0]))],[-np.log(j[0][i])*(1/decaytimex[i])+np.log(j[0][i]+j[1][i])*(1/decaytimex[i]) for i in np.arange(len(j[0]))])
                else:
                    ydataerr=zip([-np.log(j[0][i]-j[1][i])*(1/decaytimex[i])+np.log(j[0][i])*(1/decaytimex[i]) for i in np.arange(len(j[0]))],[-np.log(j[0][i])*(1/decaytimex[i])+np.log(j[0][i]+j[1][i])*(1/decaytimex[i]) for i in np.arange(len(j[0]))])
                if fdata[-1] > 3500:
                    spinsystems[nf].datasets[-1].addcpmgdata(setlabel,field,nitrotag,p2,fdata[1:-1],ydata[:-1],ydataerr[:-1],vcfs[1:-1],vdfs[1:-1],decaytimex) #[1,;]
                else:
                    spinsystems[nf].datasets[-1].addcpmgdata(setlabel,field,nitrotag,p2,fdata[1:],ydata,ydataerr,vcfs[1:],vdfs[1:],decaytimex) #[1,;]
             #   datasets[-1].addss(nf,l)
                #print np.sqrt((err/dataint)**2+(err/trint)**2)
                #print (dataint/trint)*np.sqrt((err/dataint)**2+(err/trint)**2)
    elif filetype ==1:
        with open(relpath+'expinfo','rb') as csvfile:
            filetext = csv.reader(csvfile, delimiter=' ')
            for i,row in enumerate(filetext):
                if i == 0:
                    field1=np.float(row[1])
                if i == 1:
                    field2=np.float(row[1])
        field=float(expcond['B0'])/(500*9.868831346500839)
        
        cutoffx=90
    #    datasets.append(dataset(dataposition,spinlistpath,filetype,'cest',expcond['nucleus']))
      #  datasets[-1].addcestset()
        ng=-1
        for nf,name in enumerate([w.name[0] for w in spinsystems]):
            spinsystems[nf].datasets.append(RDset(spinlistpath,filetype,dataposition,name,'cest',expcond['nucleus'],setselect))
            spinsystems[nf].datasets[-1].addexpcond(expcond)
            spinsystems[nf].datasets[-1].addshifts(setlabel,spinsystems,expcond['TR'],field,nf,len(spinsystems[nf].datasets),filetype)
            openable=0
            x=[];y=[]
            try:
                
                with open(relpath+name+'.csv','rb') as csvfile:
                    filetext = csv.reader(csvfile, delimiter=',')
                    for row in filetext:
                        x.append(float(row[0]))
                        y.append(int(row[1]))
                        
                x1xa=(np.array(x))[np.array([np.array(x),np.array(y)])[0,:].argsort()]
                y1bxa=(np.array(y))[np.array([np.array(x),np.array(y)])[0,:].argsort()]
                if cutoffx > 0:
                    x1x=[x1xa[k] for k in [i for i,j in enumerate(x1xa) if j>=cutoffx]]
                    y1bx=[y1bxa[k] for k in [i for i,j in enumerate(x1xa) if j>=cutoffx]]
                else:
                    x1x=[x1xa[k] for k in [i for i,j in enumerate(x1xa) if j<-cutoffx]]
                    y1bx=[y1bxa[k] for k in [i for i,j in enumerate(x1xa) if j<-cutoffx]]
                ax=[x1x];bx=[y1bx]
                openable=1
                ng+=1
            except:
                pass#rint 'cant open'

            if openable == 1:
                v3=np.array(smp.flatten(bx))
                v3max=v3
                v3min=v3
                #print v3
                spinsystems[nf].datasets[-1].addcestdata(setlabel,x1x,v3,v3min,v3max,field1,field2)
                #print v3
            #    datasets[-1].addss(nf,name)
                
#        for nf,name in enumerate([w.name[0] for w in spinsystems]):
            #this makes use of a very particular format of expinfo. rows 3++
            #contain: spinsystemname signalfrom signalto zeroregion1from
            # zeroregion1to zeroregion2from zeroregion2to etc
                signalrange=[]
                nullranges=[]
                zerolist=[]
                signlist=[]
                yzeros=[]
                ysign=[]
                xsign=[]
                #print relpath
                #print os.pwd()
                with open(relpath+'expinfo','rb') as csvfile:
                    filetext = csv.reader(csvfile, delimiter=' ')
                    #print name, row[0]
                    for i,row in enumerate(filetext):
                        if i > 2 and row[0] == name:
                            signalrange=[row[1],row[2]]
                            nullranges=[list(m) for m in zip(row[3::2],row[4::2])]
                    #print signalrange
                    if len(signalrange) > 0:
                        #print nf
                   #     print signalrange[0], signalrange[1]
                        for i,j in enumerate(spinsystems[nf].datasets[-1].fdata):
                            for k in nullranges:
                                if j >= np.float(k[0]) and j <= np.float(k[1]):
                                    zerolist.append(i)
                                    yzeros.append(spinsystems[nf].datasets[-1].y[i])
                                    
                            if j >=np.float(signalrange[0]) and j <= np.float(signalrange[1]):
                                signlist.append(i)
                                xsign.append(j)
                                ysign.append(spinsystems[nf].datasets[-1].y[i])
                        #print yzeros, np.average(yzeros), np.std(yzeros)
                        #print len(spinsystems[nf].datasets), nf
                        if expcond['TR'] == 'T':
                            #print expcond['TR'], 'move', spinsystems[nf].datasets[-1].nshiftcorr
                            spinsystems[nf].datasets[-1].xorig=np.array(spinsystems[nf].datasets[-1].fdata)+spinsystems[nf].datasets[-1].nshiftcorr
                            spinsystems[nf].datasets[-1].fdata=np.array(xsign)+spinsystems[nf].datasets[-1].nshiftcorr
                        else:
                            #print expcond['TR'], 'not move'#, spinsystems[nf].datasets[-1].nshiftcorr
                            spinsystems[nf].datasets[-1].xorig=np.array(spinsystems[nf].datasets[-1].fdata)
                            spinsystems[nf].datasets[-1].fdata=np.array(xsign)
                        spinsystems[nf].datasets[-1].yorig=spinsystems[nf].datasets[-1].y/np.average(yzeros)
                        spinsystems[nf].datasets[-1].y=np.array(ysign)/np.average(yzeros)
                     #   print spinsystems[nf].datasets[-1].y, spinsystems[nf].name
                        spinsystems[nf].datasets[-1].ymax=(np.array(ysign)+np.std(yzeros))/np.average(yzeros)
                        spinsystems[nf].datasets[-1].ymin=(np.array(ysign)-np.std(yzeros))/np.average(yzeros)
                        spinsystems[nf].datasets[-1].nullranges=nullranges
                        spinsystems[nf].datasets[-1].signalrange=signalrange
                        spinsystems[nf].datasets[-1].intensav=np.average(yzeros)
                        spinsystems[nf].datasets[-1].intensstd=np.std(yzeros)
                        deletethese=[j for j,i in enumerate(abs(spinsystems[nf].datasets[-1].xorig-spinsystems[nf].datasets[-1].nshiftav) < 0.5) if i == 1]
                        spinsystems[nf].datasets[-1].xorig=np.array([b for c,b in enumerate(spinsystems[nf].datasets[-1].xorig) if c not in deletethese])
                        spinsystems[nf].datasets[-1].yorig=np.array([b for c,b in enumerate(spinsystems[nf].datasets[-1].yorig) if c not in deletethese])
                        deletethese=[j for j,i in enumerate(abs(spinsystems[nf].datasets[-1].fdata-spinsystems[nf].datasets[-1].nshiftav) < 0.5) if i == 1]
                        spinsystems[nf].datasets[-1].fdata=np.array([b for c,b in enumerate(spinsystems[nf].datasets[-1].fdata) if c not in deletethese])
                        spinsystems[nf].datasets[-1].y=np.array([b for c,b in enumerate(spinsystems[nf].datasets[-1].y) if c not in deletethese])
                        spinsystems[nf].datasets[-1].ymax=np.array([b for c,b in enumerate(spinsystems[nf].datasets[-1].ymax) if c not in deletethese])
                        spinsystems[nf].datasets[-1].ymin=np.array([b for c,b in enumerate(spinsystems[nf].datasets[-1].ymin) if c not in deletethese])
    return spinsystems
                        #print spinsystems[42].datasets[-1].fdata, nf

class RDsetx:
    def __init__(self,spinlistpath,filetype,dataposition,l,datatype,nucleus,setselect):
        self.name=l
        self.setselect=setselect
        self.dataposition=dataposition
        self.spinlistpath=spinlistpath
        self.filetype=filetype
        self.datatype=datatype
        self.datathere=0
        self.nucleus=nucleus
    def addcestdata(self,fdata,y,ymin,ymax,field,field2):
        self.fdata=fdata
        self.y=y
        self.ymin=ymin
        self.ymax=ymax
        self.datathere=1
        self.field=field
        self.field2=field2
    def addexpcond(self,expcond):
        self.expcond=expcond

        
class RDset:
    def __init__(self,spinlistpath,filetype,dataposition,l,datatype,nucleus,setselect):
        self.name=l
        self.setselect=setselect
        self.dataposition=dataposition
        self.spinlistpath=spinlistpath
        self.filetype=filetype
        self.datatype=datatype
        self.datathere=0
        self.nucleus=nucleus
    def adddatapoint(self,gendatatype,xlabel,ylabel,yval,yerr1,yerr2):
        self.gendatatype=gendatatype
        self.xlabel=xlabel
        self.ylabel=ylabel
        self.yval=yval
        self.yerr1=yerr1
        self.yerr2=yerr2
        self.datathere=1
    def addcpmgdata(self,xlabel,field,nitrotag,p2,fdata,rcpmg,rcpmgerr,vcfs,vdfs,decaytimes):
        self.xlabel=xlabel
        self.field=field
        self.nitrotag=nitrotag
        self.p2=p2
        self.decaytimes=decaytimes
        self.fdata=fdata
        self.rcpmg=rcpmg
        self.rcpmgerr=rcpmgerr
        self.vcfs=vcfs
        self.vdfs=vdfs
        self.datathere=1
    def addcestdata(self,xlabel,fdata,y,ymin,ymax,field,field2):
        self.xlabel=xlabel
        self.fdata=fdata
        self.y=y
        self.ymin=ymin
        self.ymax=ymax
        self.datathere=1
        self.field=field
        self.field2=field2
    def addrelax(self,xlabel,rex, errex1, errex2,r20,err201,err202,r2f,err211,err212):
        self.xlabel=xlabel
        self.relaxthere=1
        self.datathere=1
        self.r20=r20
        self.r20err1=err201
        self.r20err2=err202
        self.rex=rex
        self.rexerr1=errex1
        self.rexerr2=errex2
        self.r2=r2f
        self.r2err1=err211
        self.r2err2=err212
    def addr1rho(self,xlabel,resomegacoll, resangcoll, resoecoll, b0field, nshift, r1rav, err1, err2):
        self.xlabel=xlabel
        self.resomegacoll=resomegacoll
        self.resangcoll=resangcoll
        self.resoecoll=resoecoll
        self.b0field=b0field
        self.nshift=nshift
        self.r1rav=r1rav
        self.err1=err1
        self.err2=err2
        self.datathere=1
        self.r1rthere=1
    def addtitdat(self,xdat,xdatx,hshifts,nshifts,heights,volumes,pw1,pw2,diffpos,shdiff,shdirect,shdirect2):
        self.xdat=xdat
        self.xdatx=xdatx
        self.hshifts=hshifts
        self.nshifts=nshifts
        self.heights=heights
        self.volumes=volumes
        self.pw1=pw1
        self.pw2=pw2
        self.diffpos=diffpos
        self.shdiff=shdiff
        self.shdirect=shdirect
        self.shdirect2=shdirect2
        self.datathere=1
    def addexpcond(self,expcond):
        self.expcond=expcond
    def addshifts(self,xlabel,spinsystems,trosy,field,ss,dsn,filetype):
        """
        adds shifts to a specific dataset of a specific spinsystem
        input variables:
        xlabel: "name of shift data set"
        spinsystems: spinsystems object which we are working on
        trosy: TROSY status (T - TROSY, X - not TROSY, S - treated as TROSY)
        field: B0 field relative to 500 MHz
        ss: number of corrent spin system of interest
        dsn: total number of datasets (giving the current, last dataset by
        subtracting one)
        filetype: see above. relevant here when it is 5 (peak lists)
        
        The following file is opened:
        shiftlistb.csv format: one line one peak. Columns:
        column1: residue name format YX.
            Y refers to the state: A - folded; B - unfolded; C - dimer
                D- other/unkown
            X number refering to sequence position of residue
        column2: either has some peak type such as 'TR' in there, or 'ref'.
        'ref' means reference.
        column3: 1H shift
        column4: 15N shift
        column5: peak intensity from CCPN Analysis fit (peak height)
        column6: Y from column1
        
        The file is read in. reference data are read in separately. There can
        be multiple lines supplying shifts or peak intensities per residue.
        Average intensities (minus reference intensities) and shifts are determined.
        Then a shift correction is calculated: one for
        TROSY file, and another one shifts the spectrum to match the 298K
        reference peak (It was determined that the proton shift is 0.0046 ppm/K
        off.). The shift correction is applied. shift corrections are stored
        in a variable for the record.
        For error determination (standard deviation), reference intensities/
        shifts are subtracted first
        
        """
        self.xlabel=xlabel
        hdat=[];ndat=[];intens=[]
       # print 'a'
        with open('shiftlistb.csv','rb') as csvfile:
            filetext = csv.reader(csvfile, delimiter=' ')
            found =0
            for row in filetext:
        #        print self.name, row[0]
                if self.name == row[0]:
                    if row[1] == 'ref':
                        found=1
                        refh=float(row[2])
                        refn=float(row[3])
                        refi=int(row[4])
                    else:
                        if filetype == 5:
                            found=1
                            refi=0;refn=0;refh=0
                        hdat.append(row[2])
                        ndat.append(row[3])
                        intens.append(row[4])
        hdat=np.array(hdat).astype('float')
        ndat=np.array(ndat).astype('float')
        intens=np.array(intens).astype('int')
        if trosy[0] == 'T' or trosy[0] == 'S':
            trshiftn=-92/(2*field*500)
            trshifth=92/(2*field*500*9.868831346500839)
        else:
            trshifth=0
            trshiftn=0
        if found == 1:
            self.hshiftav=np.average(hdat)#+trshifth-(float(spinsystems[ss].datasets[dsn-1].expcond['temp'])-298)*0.0046
            self.nshiftav=np.average(ndat)#+trshiftn
            
            self.hshiftav=np.average(hdat)+trshifth-(float(spinsystems[ss].datasets[dsn-1].expcond['temp'])-298)*0.0046
            self.nshiftav=np.average(ndat)+trshiftn

            self.nshiftcorr=trshiftn
            self.hshiftcorr=trshifth-(float(spinsystems[ss].datasets[dsn-1].expcond['temp'])-298)*0.0046
            #print field,trosy[0], np.average(hdat), np.average(hdat)+trshifth,np.average(ndat),np.average(ndat)+trshiftn
            self.intensav=np.average(np.array(intens)-refi)
            self.hshiftstd=np.std(np.array(hdat)-np.array(refh))
            self.nshiftstd=np.std(np.array(ndat)-np.array(refn))
            self.intensstd=np.std(np.array(intens)-refi)
            if filetype == 5:
                self.datathere=1
        


class shiftlist:
    def __init__(self,name):
        self.name=name
    def addshift(self,protein,state,hshift,nshift):
        self.protein=protein
        self.state=state
        self.hshift=hshift
        self.nshift=nshift
class peakwitherr:
    def __init__(self,name):
        self.name=name
    def addpeak(self,state,hshift,nshift,intens,hshifterr,nshifterr,intenserr):
        self.state=state
        self.hshift=hshift
        self.nshift=nshift
        self.intens=intens
        self.hshifterr=hshifterr
        self.nshifterr=nshifterr
        self.intenserr=intenserr
    def addexpcond(self,expcond):
        self.expcond=expcond
        
class titpeaklist:
    def __init__(self,spinlistpath,filetype,dataposition,name,datatype,nucleus,psdimprop,temp,psdimx,psdim,setselect):
        self.name=name
        self.setselect=setselect
        self.dataposition=dataposition
        self.spinlistpath=spinlistpath
        self.filetype=filetype
        self.datatype=datatype
        self.datathere=0
        self.nucleus=nucleus
        self.psdimprop=psdimprop
        self.temp=temp
        self.psdimx=psdimx
        self.psdim=psdim
    def addtitdat(self,xlabel,xdat,xdatx,hshifts,nshifts,heights,volumes,pw1,pw2,diffpos,shdiff,shdirect,shdirect2):
        self.xlabel=xlabel
        self.xdat=xdat
        self.xdatx=xdatx
        self.hshifts=hshifts
        self.nshifts=nshifts
        self.heights=heights
        self.volumes=volumes
        self.pw1=pw1
        self.pw2=pw2
        self.diffpos=diffpos
        self.shdiff=shdiff
        self.shdirect=shdirect
        self.shdirect2=shdirect2
        self.datathere=1
    def addexpcond(self,expcond):
        self.expcond=expcond

class dataset:
    def __init__(self,dataposition,spinlistpath,filetype,datatype,nucleus):
        self.dataposition=dataposition
        self.spinlistpath=spinlistpath
        self.filetype=filetype
        self.datatype=datatype
        self.datathere=0
        self.nucleus=nucleus
    def addcpmgset(self,field,nitrotag,p2,fdata,vcfs,vdfs,decaytimes):
        self.field=field
        self.nitrotag=nitrotag
        self.p2=p2
        self.decaytimes=decaytimes
        self.fdata=fdata
        self.vcfs=vcfs
        self.vdfs=vdfs
        self.datathere=1
        self.spinnum=[]
        self.spinnam=[]
    def addcestset(self):
        self.spinnum=[]
        self.spinnam=[]
        self.datathere=1
    def addss(self,nf,l):
        self.spinnum.append(nf)
        self.spinnam.append(l)
        
        # This was all done to determine decaytimes, fdata, vcfs, vdfs
        # Intensity data are now being loaded.
        
        #decaytime=vcf*(16*p2/1000000+vdf*32)
        #tdata=decaytime/(vc*32)

#        vcfs=np.sort(np.array(list(set(vcf))))
        #print vcf
        
#class CPMG_set:\
#    def __init__(self)
        



def launch(pathprefix,collection):
    os.chdir(pathprefix)
    pathlist=[]
    filetypelist=[]
    setlabels=[]
    global spinsystems1, spinsystems
    spinsystems=[]
    j=[]
    with open(collection,'rb') as csvfile:
        filetext = csv.reader(csvfile, delimiter=' ')
        for i,row in enumerate(filetext):
            pathlist.append(row[0])
            filetypelist.append(int(row[1]))
            setlabels.append(row[2])
            j.append([i])
    setselectlist=[j]
    superselector=0
    
    cnt=0
 #   print 'should start subst'
    for sspos1, setselecta in enumerate(setselectlist[superselector]):
        for sspos2, setselect in enumerate(setselecta):
            #print filetypelist[setselect],'filetypelist'
   #         print 'moreparameters',pathprefix+pathlist[setselect],filetypelist[setselect]
            spinsystems1=addsystem(spinsystems,pathprefix+pathlist[setselect],filetypelist[setselect])
            cnt+=1
    cnt=0
    for sspos1, setselecta in enumerate(setselectlist[superselector]):
        for sspos2, setselect in enumerate(setselecta):
     #       print pathprefix+pathlist[setselect],filetypelist[setselect],cnt,setselect,setlabels[setselect],pathlist, 'someparameters'
            spinsystems=adddata(spinsystems1,pathprefix+pathlist[setselect],filetypelist[setselect],cnt,setselect,setlabels[setselect],pathlist)
            cnt+=1
    spss=[spinsystems[j] for j in list(np.argsort([i.seqno for i in spinsystems]).astype(int))]
 #   try:
 #       print 'try'
  #      print spinsystems1[5].datasets[1]
 #   except:
 #       print 'x'

    for u in np.arange(len(spss)):
        try:
            xfb=spss[u].letter
        except:
            spss[u].letter='X'
 #   print spss[0].datasets[7].yval
    return [spss,setlabels]

def launchthis(pathprefix,collection,kickoutlist):
    spinsystems=[]
    datasets=[]
    #collection='concsel.dat'
    spinsystems,setlabels=launch(pathprefix,collection)
    kickoutlist=['C77','C38','A42','D38']
    
    ls=len(spinsystems)
    for i in np.arange(ls):
      #  print spinsystems[ls-i-1].name[0]
        if str(spinsystems[ls-i-1].name[0]) in kickoutlist:
            del spinsystems[ls-i-1]
    return spinsystems


def passdatatofitn(spinsystems,selres,prepro):
    """The experimental data and properties and properties, which are stored in
    the spinsystems object, are converted here to lists which are used by the
    fitting routine. The input variables are:
    spinsystems: spinsystem object
    selres: list of residues/states which are selected for further processing
    prepro: switch to select whether re-sampled data (=1, for error calculation)
        or original data (=0) are used.
    The output variables are:
    resnaml: list of residues
    timedat: lists of x-axis data lists for dispersion experiments:
        CPMG: tau_cp [insert precise definition of tau_cp here]
        Rex: two power levels converted to tau_cp; in this version of the tool,
            Rex experiments are treated like CPMG experiments which is true
            for high power levels. Convert Spinlock data to pseudo-CPMG-data:
            no spin lock: calculate pseudo CPMG frequency directly from duration\
                of tSL period (during which no spin lock is applied)
            spinlock in kHz *2*pi*1000: s-1. then divide by sqrt(3) \
                for comparability with CPMG-tau (1/2 period) according 
                to Torchia, Ishima, J Biomol NMR, 1999.
        CEST: CEST offset from main peak in ppm (0 for 0 offset)
    y: lists of y-axis data lists for relaxation dispersion 
        CPMG/CEST: either takes the original experimental data, or the resampled
            data, depending on the prepro setting
        Rex: takes the Rex result twice (is used twice because experiment is
            interpreted as pseudo-CPMG with small and large tau_cp)
            alternative, resampled data are used
    field: B0 field, relative to 500 MHz (possibly incorrect for Rex experiment
        but not used therein), for protons relative to 50 MHz
    field2: B1 field (CEST experiments), in Hz
    field3: B0 field in MHz, rounded to 10^2 (e.g. 500 rather than 501)
    tr: TROSY-type of experiments. Each dataset has a expcond dictrionary with
        a 'TR' entry: 'T' (full TROSY), 'S' (was supposed to mean semi-TROSY),
        and 'X' (not TROSY). These entries are modified here: 'S' is changed to
        'T'. The use of 'S' needs to be checked in later versions of this tool.
        also, expcond['tr'] is saved in the tr variable.
    equationtype:
        3 - CPMG experiment
        6 - Rex experiment
        > 10 - CEST experiments, set number to B1 field
    poscoll: numbers counting the spin systems, to record the order of how the 
    lists were put together in this script
    expcnd: collection of experimental conditions. From the datasets, and,
    in addition the following dictionary entries: B1field (probably a typo,
    supposed to mean B0field), type (data type as define by the entry datatype
    in the datasets[i] dictionary object), w1type (same as equationtype)
    
    In the first line of the script, there is the parameter "setprotoncpmg",
    which has to be set to 1 for proton experiments.
    future improvements: 
    - Include setprotoncpmg as argument of the function or
    eliminate it.
    - check whether field3 is calculated calculated correctly for HCPMG expts
    - eliminate duplicate output variables, why the same output variable in
    expcond and as a list
    - possibly eliminate need for this function alltogether, but check for 
    computational disadvantages
    - consider treating Rex as a on-Resonance R1rho + Hahn-echo experiment
    rather than pseudo-CPMG
    - check the meaningfulness of allowing semi-TROSY experiments as experiment
    type
    - implement handing over B1 field via appropriate variable, not equationtype
    poscoll
    - eliminate poscoll
    - check naming of 'B1field' dict entry in expcnd
    """
    setprotoncpmg=0
    timedat=[]; y=[]; err=[]; datasetno=[]; resnaml=[]; equationtype=[]
    field=[]; field2=[]; field3=[]; tr=[]; expcnd=[]; poscoll=[]
    for rr,i in enumerate(spinsystems):
        #print i
        if selres==[] or i.name[0] in selres:
            poscoll.append(rr)
            resnaml.append(i.name[0])
     #       print i.name[0], 'name'
            y.append([]);err.append([]);timedat.append([]);datasetno.append([])
            equationtype.append([]);field.append([]);field2.append([]);tr.append([])
            field3.append([])
            expcnd.append([])
            #print rr 
            removelist=[]
            for l,j in enumerate(spinsystems[rr].datasets):
   #             print j.datatype, j.datathere, 'whatever'
                if j.datatype in ['cpmg'] and j.datathere == 1:
                    timedat[-1].append(j.fdata)
                    if prepro != 0:
                        y[-1].append(j.reshufy)
                    else:
                        y[-1].append(j.rcpmg)
                    err[-1].append([np.sqrt(np.average(np.array(xx)**2)) for xx in j.rcpmgerr])
                    field[-1].append([j.field for i in np.arange(len(j.fdata))])
                    field2[-1].append([0 for i in np.arange(len(j.fdata))])
                    if setprotoncpmg == 0:
                        field3[-1].append([min([50,70,80,90], key=lambda x: \
                        abs(x-500*j.field)) for i in np.arange(len(j.fdata))])
                    else:
                        field3[-1].append([min([500,600,800,900], key=lambda x:\
                        abs(x-500*j.field)) for i in np.arange(len(j.fdata))])
                    
                    
                    equationtype[-1].append([3 for i in np.arange(len(j.fdata))])
                    tr[-1].append([j.expcond['TR'] for i in np.arange(len(j.fdata))])
                    expcnd[-1].append(j.expcond)
                    if setprotoncpmg == 0:
                        expcnd[-1][-1]['B1field']=min([50,70,80,90], key=lambda x: abs(x-500*j.field))
                    else:
                        expcnd[-1][-1]['B1field']=min([500,600,800,900], key=lambda x: abs(x-500*j.field))
                    expcnd[-1][-1]['type']=j.datatype
                    expcnd[-1][-1]['w1type']=3
                elif j.datatype in ['Rex'] and j.datathere == 1:
                    tSL=j.expcond['tSL']
                    SLHz=j.expcond['SLHz']
                    timedat[-1].append([int(round(1/(int(tSL)/1000000),0)), int(round(float(SLHz)*3628,0))])
                    if prepro != 0 and prepro != 1:
                        y[-1].append(j.reshufy)
                    else:
                        y[-1].append([j.yval,j.yval])
                    err[-1].append([np.average([j.yerr1,j.yerr2]),np.average([j.yerr1,j.yerr2])])
                    field[-1].append([float(j.expcond['B1'])/500,float(j.expcond['B1'])/500])
                    field2[-1].append([0,0])
                    field3[-1].append([min([50,70,80,90], key=lambda x: abs(x-float(j.expcond['B1']))) for i in np.arange(2)])
                    equationtype[-1].append([6 for i in np.arange(2)])
                    tr[-1].append([j.expcond['TR'],j.expcond['TR']])
                    expcnd[-1].append(j.expcond)
                    expcnd[-1][-1]['B1field']=min([50,70,80,90], key=lambda x: abs(x-float(j.expcond['B1'])))
                    expcnd[-1][-1]['type']=j.datatype
                    expcnd[-1][-1]['w1type']=6
                elif j.datatype in ['cest'] and j.datathere == 1:
                    timedat[-1].append(j.fdata-float(j.nshiftav))
                    if prepro != 0:
                        y[-1].append(j.reshufy)
                    else:
                        y[-1].append(j.y)
                    err[-1].append((j.ymax-j.ymin)/2)
                    field[-1].append([j.field/500 for i in np.arange(len(j.y))])
                    field2[-1].append([j.field2 for i in np.arange(len(j.y))])#not used for now. shall substitute equationtype
                    field3[-1].append([min([50,70,80,90], key=lambda x: abs(x-j.field)) for i in np.arange(len(j.y))])
                    tr[-1].append(['T' if x == 'S' else x for x in [j.expcond['TR'] for i in np.arange(len(j.y))]])
                    equationtype[-1].append([j.field2 for i in np.arange(len(j.y))])
                    expcnd[-1].append(j.expcond)
                    expcnd[-1][-1]['B1field']=min([50,70,80,90], key=lambda x: abs(x-j.field))
                    expcnd[-1][-1]['type']=j.datatype
                    expcnd[-1][-1]['w1type']=j.field2
                    if expcnd[-1][-1]['TR'] == 'S':
                        expcnd[-1][-1]['TR'] = 'T'
                else:
                    removelist.append(l)
            for j in removelist[::-1]:
                del spinsystems[rr].datasets[j]
    return resnaml,timedat,y,err,field,field2,field3,tr,equationtype,poscoll,expcnd