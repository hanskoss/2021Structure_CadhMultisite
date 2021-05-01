#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This scripts contains functions related to fitting relaxation dispersion data

"""
from __future__ import division
from hkimports2 import flatten
from hkimports2 import np
from hkimports2 import optimize
from hkimports2 import time
from hkimports2 import os
#from hkimports2 import scipyexpm
import preprocessing as prepro    #used preprocess - passdatatofit previously
import iofunctions as hkio
import RDmath as hkRDmath
#print 1
fineres=1
#[  2.76800000e-03  -1.47694090e+03   2.07130011e-01   8.57390488e+00]
#linux='/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/'
#windows='C:\\Users\\Hans\\Desktop\\TRANSFER\\2020Feb\\'
#path2020=windows#linux#windows

def multifunctg2(params,x,m,dwbset,g):
    """performs fitting for a variety of relaxation dispersion
    experiments. Beause r20 of a higher field >= r20 of a lower field, r20
    is constructed by combining 1-2 parameters: R20 @ 500, and a scaling factor
    >= 1 for R20 @ B0 > 500.
    input variables: params (parameters); x: x values (1/(4*tau_cp) or offset);
    m=magnetic field B0 relative to 500 MHz, to be multiplied with delta_omegas (dwx, dcx)
    delta_omega refers to the chemical shift of the minor state from the major state
    in s-1 at 500 MHz. dwbset: not used anymore.
    g: type of experiment. 3: CPMG; 6: Rex. Can chose whether to calculate via pseudo-
    CPMG data or r1rho/free precession; 5: R1rho (not used); > 10: CEST (g is than
    equal to the B1 field in Hz). Hard-coded: R1=1.8 (arbitrary, virtually doesn't really change
    CEST curves); duration of CEST experiment: 400 ms.
    Possible Improvements: might reconsider hard-coded R1 at revision stage.
    might test more extensively whether it better to use a psuedo-CPMG model or a r1rho/free
    model to calculate Rex.
    """
    pbx,pcx,kexx,k13x,k23x=params[0:5]
    try:
        r20=params[7]*params[8]
    except:
        r20=params[7]

    result=[]
    g=g[0]
    if g > 1000:
        calcchoice = 0
        g=g-1000
    else:
        calcchoice = 1
    if g == 3:
        tx=1/(4*np.array(x)); #tau_cp
        for b,a in enumerate(tx):
            dwx=m[b]*params[5]
            dcx=m[b]*params[6]
            """ for step 5 and 6: exact calculation """
            result.append(hkRDmath.nEVapprox(a,dwx,pbx,kexx,2,0,0,0,dcx,pcx,k13x,k23x,0,0,0,0,0)+r20) 
            """ for steps 1-4: Exp0Log2Lamda2 approximation from
            Koss, H., Rance, M., & Palmer, A. G., 3rd. (2018). General \
            Expressions for Carr-Purcell-Meiboom-Gill Relaxation Dispersion \
            for N-Site Chemical Exchange. Biochemistry, 57(31), 4753-4763."""
            #result.append(hkRDmath.nEVapprox(a,dwx,pbx,kexx,2,0,2,2,dcx,pcx,k13x,k23x,0,0,0,0,0)+r20) 
    elif g == 6:
        #
        tx=1/(np.array(x))   #tx is exact duration of tau in tau-pi-tau-tau-pi-tau
        #print x,tx, 'tx2' 
        wx=np.array(x)*(np.sqrt(3))/(2*np.pi)
        for b,a in enumerate(tx):
            dwx=m[b]*params[5]
            dcx=m[b]*params[6]
            """Two distinct ways to calculate 15N-Rex:
            The low-power data point is always calculated as CPMG-like.
            The high-power point is calculated either as originating from a 
            CPMG-like experiment, or from an R1rho-like experiment.
            The R1rho-like calculation is a little bit more accurate but
            less stable (eigenvalue calculation can give extreme results).
            For Step 1-3, the CPMG-like calculation is preferred, then R1rho-
            like
            """
     #       rex1old=hkRDmath.r1req(0,dwx,pbx,kexx,2,4*wx[1]*np.pi*2,0,0,dcx,pcx,k13x,k23x,0,0,0,0,0)
            rex1=hkRDmath.r1req(0,dwx,pbx,kexx,2,wx[1]*np.pi*2,0,0,dcx,pcx,k13x,k23x,0,0,0,0,0)
            cpmg1=hkRDmath.nEVapprox(tx[1],dwx,pbx,kexx,2,0,0,0,dcx,pcx,k13x,k23x,0,0,0,0,0)
            cpmg2=hkRDmath.nEVapprox(tx[0],dwx,pbx,kexx,2,0,0,0,dcx,pcx,k13x,k23x,0,0,0,0,0)
     #       viarexold=cpmg2-rex1old
            viarex=cpmg2-rex1
            viacpmg=cpmg2-cpmg1#hkRDmath.nEVapprox(tx[0],dwx,pbx,kexx,2,0,0,0,dcx,pcx,k13x,k23x,0,0,0,0,0)-hkRDmath.nEVapprox(tx[1],dwx,pbx,kexx,2,0,0,0,dcx,pcx,k13x,k23x,0,0,0,0,0)
            if calcchoice == 1:
                result.append(viarex)
            else:
                result.append(viacpmg)
    elif g > 10:
        tx=np.array(x)
        dwx=params[5]
        dcx=params[6]
        for b,a in enumerate(tx):
            """R1 has a very small influence on the CEST calculation. We have
            simulated an average R1 for the C11 monomer at 700 and 800 MHz
            which is selected here.
            """
            if m[b] < 0.15:
                r1set=0.966
            else:
                r1set=1.144 #was 1.8 in test phase
            datapointx=hkRDmath.cestfunction([a],dwx,dcx,kexx,k13x,k23x,pbx,pcx,g,r1set,r20,m[b],1)
            result.append(datapointx)
    return result


def errfunctg3(parax,x,exp_data,m,dwbset,err,g,fl):
    """The main error function for the global relaxation dispersion data 
    fitting process"""
    value=[]
    parnocoll=[]
    cestchi=[];rdchi=[]
    for l,i in enumerate(fl):
        a=[parax[j] for j in i]
        parnocoll.append([j for j in i])
        val=np.array((flatten([multifunctg2(a,x[l],m[l],dwbset[l],g[l])])))
        value.append(val)
        if g[l][0] > 10:
            cestchi.append((np.array(exp_data[l])-np.array(val))/(np.array(err[l])))
        else:
            rdchi.append((np.array(exp_data[l])-np.array(val))/(np.array(err[l])))#        print np.average(((np.array(exp_data[l])-np.array(val))/(np.array(err[l])))**2), g[l][0]
    paramno=len(set(flatten(parnocoll)))
    x1=((len(flatten(cestchi)))/(len(flatten(cestchi))+len(flatten(rdchi))))
    x2=(np.sqrt(2)*np.array(flatten(cestchi))*(1/np.sqrt(len(flatten(exp_data))-paramno)))
    x3=((len(flatten(rdchi)))/(len(flatten(cestchi))+len(flatten(rdchi))))
    x4=(np.sqrt(2)*np.array(flatten(rdchi))*(1/np.sqrt(len(flatten(exp_data))-paramno)))
    return np.array(flatten([x4,x2]))

def printrd(praxs1,x,exp_data,m,err,g,dwbset,par,fl):
    """This function recalculates the theoretical value and chi square from
    parameters and experimental data."""
    value=[];valuealt=[]
    chicoll=[]
    chicoll2=[]
    a,b1,b2,c,e=par.getallparandbnds(praxs1,['p','k','dw','R20500','R2mult'],inclfilt=[])
    par1=a[0]
    parnocoll=[]
    for l,i in enumerate(fl):
        thisset=[j for j in i]
        parnocoll.append(thisset)
    paramno=len(set(flatten(parnocoll)))
    cestchi=[];rdchi=[]
    xalt=[list(np.arange(np.min(xi),np.max(xi),20)) for xi in x]
    for l,i in enumerate(fl):
        a=[par1[j] for j in i]
        thisset=[j for j in i]
        val=np.array((flatten([multifunctg2(a,x[l],m[l],dwbset[l],g[l])])))
        """ comment in the following line if a higher resolution CPMG curve is required"""
        if fineres == 1:
            #print len([m[l][0] for xa in xalt]), len([g[l][0] for xa in xalt]), len(xalt[l])
            valalt=np.array((flatten([multifunctg2(a,xalt[l],[m[l][0] for xa in xalt[l]],dwbset[l],[g[l][0] for xa in xalt[l]])])))
        valuealt.append(valalt)
        value.append(val)
        if g[l][0] > 10:
            cestchi.append((np.array(exp_data[l])-np.array(val))/(np.array(err[l])))
        else:
            rdchi.append((np.array(exp_data[l])-np.array(val))/(np.array(err[l])))
        chicoll2.append((np.array(exp_data[l])-np.array(val))/(np.array(err[l])))
        chicoll.append(np.sum(((np.array(exp_data[l])-np.array(val))/\
            (np.array(err[l])*np.sqrt(len(exp_data[l])-paramno*\
            (len(exp_data[l])/len(flatten(exp_data))))))**2))
    paramno=len(set(flatten(parnocoll)))
    err0=np.array(flatten(exp_data))-np.array(flatten(value))
    x1=((len(flatten(cestchi)))/(len(flatten(cestchi))+len(flatten(rdchi))))
    x2=np.sum((np.sqrt(2)*np.array(flatten(cestchi))*(1/np.sqrt(len(flatten(exp_data))-paramno)))**2)/2
    x3=((len(flatten(rdchi)))/(len(flatten(cestchi))+len(flatten(rdchi))))
    x4=np.sum((np.sqrt(2)*np.array(flatten(rdchi))*(1/np.sqrt(len(flatten(exp_data))-paramno)))**2)/2
    if fineres == 1:
        value=valuealt
        #print value, 'value'
    return value, chicoll,chicoll2, (x4*x1+x3*x2)*2,(np.array(flatten(err0))/\
            (np.array(flatten(err))*np.sqrt(len(flatten(exp_data))-paramno)))


def fitcpmg4(praxs1,timedat,rawdata,field,err,mode,equationtype,conditions,savstatdir,filenamsav,resnam,poscoll,moreconditions,paramsx,fl):
    """ Fitting multiple types of relaxation dispersion data (function title somewhat misleading).
    input arguments: praxs1: property axis collection; timedat: x axis data; rawdata: \
    y axis data; field: B0 field relative to 500 MHz, err: y error, mode: not used here, \
    was used in the parent script as precalc; equationtype: type of relaxation dispersion, \
    see multifunctg2 function for details;
    conditions: list of conditions: pos 0 - reshuffle, \
    this is not used anymore, but not entirely deleted yet (check). pos 1 - details \
    about many steps and attempts should be made for fitting. pos 1/0: number of initial \
    attempts of the "precalculation" round to get to a small chi square at the beginning \
    of the fitting procedure. High number can be useful when starting in a Monte Carlo \
    manner at the beginning of the project (large boundaries, little idea about the system);
    pos 1/1 precalcatt: Number of documented major calculation steps during the precalculation; \
    pos 1/2 precalclen: Number of undocumented steps within each major precalculation step; \
    pos 1/3 maincalcatt: Number of documented major calculation steps during the main calculation;\
    pos 1/4 maincalclen: Number of undocumented steps within each major main calculation step.
    filenamsav: Name under which progess of fit is saved to disk
    resnam: list of residue names included in fit (not used here)
    poscoll: (not used anymore, has to do with being able to undo some sorting action)
    moreconditions: package of conditions used by parent script, this is supposed to be saved with\
    the status and therefore an input argument.
    paramsx: collection of parameter objects, can have different bounds etc
    fl: list of pointer lists. Each list member point to certain parameter list positions, \
    there is one pointer list for each data point.
    """
    allconditions=[timedat,rawdata,field,err,mode,equationtype,conditions,\
                   filenamsav,resnam,poscoll,moreconditions]
    reshuffle=conditions[0] #0 or 1
    numbattempts,precalcatt,precalclen,maincalcatt,maincalclen=conditions[1][0:5]#1, 5, 5, 5, 10

    numdat=np.shape(rawdata)[0]
    dwbsetp=[]
    setdwb=0
    for i in np.arange(len(timedat)):
        for j in np.arange(len(flatten(timedat[i]))):
            dwbsetp.append(setdwb)
        setdwb+=1
    par6=np.array(field)
    gpar=np.array(equationtype)
    par7=np.array(np.array(dwbsetp).astype('int'))
    par2=np.array(timedat)
    par3=np.array(rawdata)
    errvalpar=np.array(err)
 #   print len(boundsl)
    if reshuffle == 1:
        par7x=par7;par6x=par6;par2x=par2;par3x=par3;errvalparx=errvalpar;gparx=gpar
        par7=[];par6=[];par2=[];par3=[];errvalpar=[];gpar=[]
        for i,j in enumerate(par2x):
            if gparx[i] == 6:
                par2.append(par2x[i])
                par3.append(par3x[i])
                par6.append(par6x[i])
                par7.append(par7x[i])
                errvalpar.append(errvalparx[i])
                gpar.append(gparx[i])
            else:
                foundrnd=0
                np.random.seed()
                while foundrnd == 0:
                    np.random.seed()
                    n=np.random.randint(len(par2x))
                    if gparx[n] == 3:
                        par2.append(par2x[n])
                        par3.append(par3x[n])
                        par6.append(par6x[n])
                        par7.append(par7x[n])
                        errvalpar.append(errvalparx[n])
                        gpar.append(gparx[n])
                        foundrnd=1
    allrescoll=[]
    par1coll=[]
    costcoll=[]
    for u in np.arange(numbattempts):
        np.random.seed()
        evalmode=1
        a,b1,b2,c,e=paramsx.getallparandbnds(praxs1,['p','k','dw','R20500','R2mult'],inclfilt=[])
        par1=a[u];boundsl=b1[u];boundsh=b2[u]
        if evalmode == 1:
            for i in zip(par1,boundsl,boundsh):
               if i[0] <= i[1] or i[0] >= i[2]:
                   print 'problem!', i[0], i[1], i[2]
        if precalcatt > 0:
            try:
                for k in np.arange(precalcatt):
                    if evalmode ==1:
                        for i in zip(par1,boundsl,boundsh):
                            print i[0],i[1],i[2]
                 #   print par1,par2, par3, par6, par7, errvalpar, gpar, fl
               #     print 'here'
                #    print gpar, 'gparold'
                    gparn=[]
                    for ggg in gpar:
                        gparn.append(list([gggg+1000 for gggg in ggg]))
                    #gpar=np.array([ggg+1000 for ggg in gpar])
                 #   print gparn, 'gparnew'
                    res=optimize.least_squares(errfunctg3,par1,max_nfev=precalclen,\
                        bounds=(boundsl,boundsh),args=(par2,par3,par6,par7,\
                        errvalpar,gparn,fl),method='trf',jac='3-point',x_scale='jac') #,
                    par1=res.x
                    allrescoll.append(res)
                    print 'attempt ', u, ' precalculation step ', u, k, par1,res.cost,filenamsav, allrescoll[-1].cost
                    hkio.savstatus2b(savstatdir,filenamsav,resnam,poscoll,allrescoll,allconditions)
            except:
                print 'well this one didnt work'
            par1coll.append(res.x)
            costcoll.append(res.cost)
    if precalcatt > 0:
        try:
            par1=par1coll[np.argmin(costcoll)]
        except:
            print "unfeasable result"
#    try:
    for k in np.arange(maincalcatt):
        print len(par1), len(par2), len(flatten(par2)), 'well'
        res=optimize.least_squares(errfunctg3,par1,max_nfev=maincalclen,\
            bounds=(boundsl,boundsh),args=(par2,par3,par6,par7,errvalpar,gpar,\
            fl),method='trf',jac='3-point',x_scale='jac')
        allrescoll.append(res)
        hkio.savstatus2b(savstatdir,filenamsav,resnam,poscoll,allrescoll,allconditions)
        par1=res.x
            
        print 'mainalculation step ', u, k, par1, res.cost,filenamsav, allrescoll[-1].cost
    for i in res.x:
        print i
    print 'final cost', res.cost
#    except:
#        res=[0]
#        print 'late fitting error'
    return res, allrescoll #xtol=1e-9

def reshuffle(ss,reslalmall,shuffletype):
    """ Resampling for error calculation by adding or subtracting residuals
    from calculated experimental data.
    """
    signrnd=1
    heteroskedacity=1
    withoutreplacement=0
    allresid={}
    alldsref={}
    dplcoll=[]
    datatypes=[i[0] for i in shuffletype]
    for dt in datatypes:
        allresid[dt]=[]
        alldsref[dt]=[]
    dsnum=0
    """calculating residuals"""
    for spinsyst in ss:
        selnam = spinsyst.name[0]
        if selnam in reslalmall:
            dpl=[]
            resid={}
            dsref={}
            if spinsyst.datasets[-1].setselect > dsnum:
                dsnum=spinsyst.datasets[-1].setselect
            for dt in datatypes:
                resid[dt]=[]
                dsref[dt]=[]
         #   setparameters2=[dataname,'/home/hanskoss/data/Cadherin/nmrCad/procandcoll/TSnewsort/2020Feb/',[selnam],conditions,namresults]
            for ds in spinsyst.datasets:
                
                if signrnd != 1:
                    signrnd=np.random.choice([-1,1],size=len(ds.fit))
                if ds.datatype == 'cpmg':
                    if heteroskedacity == 1:
                        errlist=np.array([np.sqrt(np.average(np.array(xx)**2)) for xx in ds.rcpmgerr])
                    else:
                        errlist=1
                    ds.resid=signrnd*(ds.rcpmg-ds.fit)/errlist
                elif ds.datatype == 'cest':
                    if heteroskedacity == 1:
                        errlist=(ds.ymax-ds.ymin)/2
                    else:
                        errlist=1
                    ds.resid=signrnd*(ds.y-ds.fit)/errlist
                elif ds.datatype == 'Rex':
                    if heteroskedacity == 1:
                        errlist=[np.average([ds.yerr1,ds.yerr2]),np.average([ds.yerr1,ds.yerr2])]
                    else:
                        errlist=1
                    ds.resid=signrnd*(ds.yval-ds.fit)/np.array(errlist)
             #   print ds.resid
             #   print np.average(ds.resid), np.std(ds.resid)
                dsref[ds.datatype].append(ds.setselect)
                resid[ds.datatype].append(ds.resid)
            for dt in datatypes:
                allresid[dt].append(resid[dt])
                alldsref[dt].append(dsref[dt])
            dplcoll.append(dpl)
    
    dspool={}
    sspool={}
    """for each data point, a random residual from a set of eligible residuals
    is selected. This ways, the residuals for a certain group of data points
    are mixed (with replacement).
    method "dataset": include all points belonging to a certain datatype and
    dataset (across different residues)
    method "spinsyst": include all points belonging to a certain data type
    and residue (across different datasets)
    method "any": include all points belonging to a certain data type (any 
    dataset and any spin system)
    method "each": include all points belonging to a certain dataset, data type
    and residue
    """
    for dt,method in shuffletype:
        dspool[dt]=[[] for i in np.arange(dsnum+1)]
        sspool[dt]=[]
        for xx,x in enumerate(allresid[dt]):
            
            for z,y in enumerate(x):
                if 'y' != []:
                    dspool[dt][alldsref[dt][xx][z]].append(y)
            sspool[dt].append(flatten(x))
        dspool[dt]=[flatten(i) for i in dspool[dt]]
        for q,x in enumerate(allresid[dt]):
            for z,y in enumerate(x):
                if method == 'dataset':
                    pool=dspool[dt][alldsref[dt][q][z]]
                elif method == 'spinsyst':
                    pool=sspool[dt][q]
                elif method == 'any':
                    pool=flatten(allresid[dt])
                elif method == 'each':
                    pool=allresid[dt][q][z]
                np.random.seed()
                if pool != []:
                    
                    if withoutreplacement == 0:
                        """with replacement works always"""
                        chosenwere=np.random.randint(len(pool),size=len(allresid[dt][q][z]))
                        """without replacement is not implemented for the 'any' method at this time."""
                    else:
                        chosenwere=np.random.choice(np.arange(len(pool)),size=len(allresid[dt][q][z]),replace=False)
                        allresid[dt][q][z]=np.array([pool[i] for i in chosenwere])
                        for delthis in np.sort(chosenwere)[::-1]:
                            try:
                                del(pool[delthis])
                            except:
                                np.delete(pool,delthis)
    
    #z=0
    """calculate resampled data points, used as experimental data for error
    estimation"""
    for dt in datatypes:
        q=0
        for ssn,spinsyst in enumerate(ss):
            selnam = spinsyst.name[0]
            if selnam in reslalmall:
                z=0
                for dsn,ds in enumerate(spinsyst.datasets):
                    if ds.datatype == dt:
                        if heteroskedacity == 1:
                            
                            if ds.datatype == 'cpmg':
                                errlist=np.array([np.sqrt(np.average(np.array(xx)**2)) for xx in ds.rcpmgerr])
                            elif ds.datatype == 'cest':
                                errlist=(ds.ymax-ds.ymin)/2
                            elif ds.datatype == 'Rex':
                                errlist=np.array([np.average([ds.yerr1,ds.yerr2]),np.average([ds.yerr1,ds.yerr2])])
                        else:
                            errlist = 1
                        print ds.fit, allresid[dt][q][z], errlist, ds.fit+allresid[dt][q][z]*errlist, 'reshuffled'
                        ss[ssn].datasets[dsn].reshufy=ds.fit+allresid[dt][q][z]*errlist
                        if ds.datatype == 'Rex':
                            ss[ssn].datasets[dsn].reshufy=[ss[ssn].datasets[dsn].reshufy[0] for i in np.arange(len(ss[ssn].datasets[dsn].reshufy))]
                        z+=1
                q+=1
    return ss

def runfit4b(praxs1,ctd,selresidues,precalc,resnam,conditions,path2020,savstatdir,files,filenamsav,paramsx,drawonly,cond):
    """ This function prepares the global fit
    input arguments
    praxs1: property axis collections
    ctd
    precalc: triggers different preparatory routines depending on this switch
        0: regular fit, no pre-calculated data
        1: uses re-sampled data supplied by cond[1].
        2: loads a previously performed calculation and calls a dedicated \
            re-sampling function.
    
    
    """
    global spinsystems

    moreconditions=[ctd,selresidues,precalc,resnam,files]
    print 'gethere'
    if precalc != 0 and precalc != 1:
        spinsystems=hkio.loadss(savstatdir,precalc)
        reslalmall=selresidues[0]#[reslall[i] for i in pickthese]
        shuffletype=[['cpmg','dataset'],['Rex','dataset'],['cest','each']]
        spinsystems=reshuffle(spinsystems,reslalmall,shuffletype)
    elif precalc == 1:
        print 'and'
        reslalmall=selresidues
   #     print reslalmall, 'reslalmall'
        spinsystems,setlabels=prepro.launch(path2020,files)
        nr=0
        for ssn,spins in enumerate(spinsystems):
            if spins.name[0] in reslalmall:
       #         print ''
        #        print spins.name[0]
                #dsn=0
                for dsn,ds in enumerate(spins.datasets):
                    if ds.datatype in ['Rex', 'cpmg', 'cest']:
                        if ds.datatype == 'Rex':
                            spinsystems[ssn].datasets[dsn].reshufy=cond[1][nr][0]
                        else:
                            spinsystems[ssn].datasets[dsn].reshufy=cond[1][nr]
          #              print 'reshuf'
          #              print spinsystems[ssn].datasets[dsn].reshufy

                        nr+=1
                    #    dsn+=1
        
    elif precalc == 0:
 #       print 'path'
  #      print path2020,files,'!!!!!!'
        spinsystems,setlabels=prepro.launch(path2020,files)
        print spinsystems[35].name, spinsystems[35].datasets[0].rcpmg
        
    #    print vars(spinsystems[55].datasets[5])
            #print i.xlabel
  #  print 'p2'
    if len(ctd) == 0:
        resultcolll=[]
        poscolll=[]
        resnaml=[]
        allresultcoll=[]
    else:
        resultcolll=ctd[0]
        poscolll=ctd[1]
        resnaml=resnam
        allresultcoll=[]
    for selectdatn in selresidues:
      #  global expcnd
        resnaml,timedat,rawdata,errd,field,field2,field3,tr,equationtype,poscoll,expcnd=prepro.passdatatofitn(spinsystems,selectdatn,precalc)
   #     print rawdata, 'rawdata'
   #ss[     print selectdatn, precalc, rawdata, 'step'
        filters=[[['residues','name'],[[rn] for rn in resnaml]],[['conc','value'],[['2.475'],['9.9']]],[['TR','name'],[['T'],['X']]],[['B1field','rounded'],[[50],[70],[80],[90]]],[['type','name'],[['cpmg'],['Rex'],['cest']]]]#[['w1type','numb'],[[3],[6],[15],[25],[50]]]]
        filters2=[[['conc','value'],[['2.475'],['9.9']]],[['TR','name'],[['T'],['X']]],[['B1field','rounded'],[[50],[70],[80],[90]]],[['type','name'],[['cpmg'],['Rex'],['cest']]]]#[['w1type','numb'],[[3],[6],[15],[25],[50]]]]
        selset=[]
        q=0
        explist=[]
        for j,i in enumerate(expcnd):
            for l,k in enumerate(i):
                selsetx=[j]
              #  print j,i,l,k
                for n,m in enumerate(filters2):
                    selsetx.append([p for p,o in enumerate(m[1]) if k[m[0][0]] in o][0])
                selset.append(selsetx)
                explist.append(q)
                q+=1
        filt=[]
        aa=selset
       # print aa
       # print len(aa), 'aa'
        bb=set(tuple(ix) for ix in selset)
        bb=[list(b) for b in aa]
        seldatasets=list(np.arange(len(aa)))
       # print seldatasets, 'seldatasets', aa
        inclfx=[]
        for l,i in enumerate(bb):
            inclf=[]
            for k,j in enumerate(i):
                inclf.append([filters[k][0][0],filters[k][0][1],filters[k][1][j]])
                    #inclf.append([filters[k-1][0][0],filters[k-1][0][1],filters[k-1][1][j]])
                inclfx.append(inclf)
            a,b1,b2,c,e=paramsx.getallparandbnds(praxs1,['p','k','dw','R20500','R2mult'],inclfilt=inclf)
            f=flatten([np.array(e[m][:-1])+np.sum([k for k in [0]+[e[j][-1] for j,i in enumerate(e) if j < len(e)-1]][0:(l+1)]) for m,l in enumerate(np.arange(len(e)))])
            filt.append(f)
        timedat=[flatten(timedat,levels=1)[i] for i in seldatasets]
        rawdata=[flatten(rawdata,levels=1)[i] for i in seldatasets]
        errd=[flatten(errd,levels=1)[i] for i in seldatasets]
        field=[flatten(field,levels=1)[i] for i in seldatasets]
        field2=[flatten(field2,levels=1)[i] for i in seldatasets]
        field3=[flatten(field3,levels=1)[i] for i in seldatasets]
        tr=[flatten(tr,levels=1)[i] for i in seldatasets]
        
        equationtype=[flatten(equationtype,levels=1)[i] for i in seldatasets]
        #paramsx.getallparandbnds()
        a,b1,b2,c,e=paramsx.getallparandbnds(praxs1,['p','k','dw','R20500','R2mult'],inclfilt=[])
        f=flatten([np.array(e[m][:-1])+np.sum([k for k in [0]+[e[j][-1] for j,i in enumerate(e) if j < len(e)-1]][0:(l+1)]) for m,l in enumerate(np.arange(len(e)))])
        poscolll.append(poscoll)
        if drawonly == 0:
            sojetzt,allrescoll=fitcpmg4(path2020,savstatdir,praxs1,timedat,rawdata,field,errd,precalc,equationtype,conditions,filenamsav,resnaml, poscolll,moreconditions,paramsx,filt)
        else:
            dwbsetp=[]
            setdwb=0
            for i in np.arange(len(timedat)):
                for j in np.arange(len(flatten(timedat[i]))):
                    dwbsetp.append(setdwb)
                setdwb+=1
            par7=np.array(np.array(dwbsetp).astype('int'))
       #     print timedat, 'timedat'
        #    print rawdata, 'rawdata'
        #    print field, 'field'
         #   print errd, 'errd'
            fittedcurve,chsq0,chsq1,chsq2,chsq3=printrd(praxs1,timedat,rawdata,field,errd,equationtype,par7,paramsx,filt)
      #      print fittedcurve, '!!!!!'
       #     print "overall chi quare", chsq2
            for j,i in enumerate(seldatasets):
                print i, 'datasetno', chsq0[j], equationtype[j][0]
       # print chsq0, chsq1, chsq2, chsq3
        if drawonly == 0:
            try:
                resultcolll.append(sojetzt)
                allresultcoll.append(allrescoll)
            except:
                print 'ugh1'
        else:
            for j in poscoll:
                for k,i in enumerate(seldatasets): #fittedcurve:
                    spinsystems[j].datasets[i].fit=fittedcurve[k]
    if drawonly == 0:
        return resultcolll, resnaml, poscolll, spinsystems, allresultcoll
    else:
        return spinsystems#else:

def runfit4(praxs1,ctd,selresidues,precalc,resnam,conditions,path2020,savstatdir,files,filenamsav,paramsx,drawonly,cond):
    """ This function prepares the global fit and converts data structures where
    appropriate.
    input arguments:
    praxs1: property axis collections
    ctd not used
    precalc: triggers different preparatory routines depending on this switch
        0: regular fit, no pre-calculated data
        (spinsystems object) containing pre-calcualted theoretical data for
        resampling
    
    """

    moreconditions=[ctd,selresidues,precalc,resnam,files]
    if precalc != 0 and precalc != 1:
        spinsystems=hkio.loadss(savstatdir,precalc)
        reslalmall=selresidues[0]#[reslall[i] for i in pickthese]
        shuffletype=[['cpmg','dataset'],['Rex','dataset'],['cest','each']]
        spinsystems=reshuffle(spinsystems,reslalmall,shuffletype)
    elif precalc == 0:
        spinsystems,setlabels=prepro.launch(path2020,files)
        print spinsystems[35].name, spinsystems[35].datasets[0].rcpmg
        
   # print spinsystems[57].datasets[13].xlabel
    resultcolll=[];
    poscolll=[]
    resnaml=[]
    allresultcoll=[]
    """This is a loop allowing to test various residue combination sets."""
    for selectdatn in selresidues:
        """The experimental data are stored in the spinsystems object; 
        the parameters are stored in the parameters object, using certain
        property axes. For the fitting engine, all data are flattened form.
        For each data point, there is a corresponding pointer set in a list of
        equal lengths which selects the appropriate parameters.
        All flattened lists, including pointer lists, are generated here.
        """
        resnaml,timedat,rawdata,errd,field,field2,field3,tr,equationtype,poscoll,expcnd=prepro.passdatatofitn(spinsystems,selectdatn,precalc)
        """hard-coded filters, has to be modified for other experimental
        combinations"""
        setprotoncpmg=1
        if setprotoncpmg == 1:
            
            filters=[[['residues','name'],[[rn] for rn in resnaml]],[['conc','value'],\
                    [['2.475'],['9.9']]],[['TR','name'],[['T'],['X']]],\
                    [['B1field','rounded'],[[500],[600],[800],[900]]],\
                    [['type','name'],[['cpmg'],['Rex'],['cest']]]]
            filters2=[[['conc','value'],[['2.475'],['9.9']]],\
                      [['TR','name'],[['T'],['X']]],\
                      [['B1field','rounded'],[[500],[600],[800],[900]]],\
                      [['type','name'],[['cpmg'],['Rex'],['cest']]]]
        else:
            filters=[[['residues','name'],[[rn] for rn in resnaml]],[['conc','value'],\
                    [['2.475'],['9.9']]],[['TR','name'],[['T'],['X']]],\
                    [['B1field','rounded'],[[50],[70],[80],[90]]],\
                    [['type','name'],[['cpmg'],['Rex'],['cest']]]]
            filters2=[[['conc','value'],[['2.475'],['9.9']]],\
                      [['TR','name'],[['T'],['X']]],\
                      [['B1field','rounded'],[[50],[70],[80],[90]]],\
                      [['type','name'],[['cpmg'],['Rex'],['cest']]]]
            
            
        selset=[]
        q=0
        explist=[]
      #  print len(expcnd)
        for j,i in enumerate(expcnd):
           # print i, len(i)
            for l,k in enumerate(i):
                selsetx=[j]
                for n,m in enumerate(filters2):
                    selsetx.append([p for p,o in enumerate(m[1]) if k[m[0][0]] in o][0])
                selset.append(selsetx)
                explist.append(q)
                q+=1
        filt=[]
        
        aa=selset
       # print len(aa), 'aa'
        bb=set(tuple(ix) for ix in selset)
        bb=[list(b) for b in aa]
        seldatasets=list(np.arange(len(aa)))
        inclfx=[]
        for l,i in enumerate(bb):
            inclf=[]
            for k,j in enumerate(i):
                inclf.append([filters[k][0][0],filters[k][0][1],filters[k][1][j]])
                inclfx.append(inclf)
            a,b1,b2,c,e=paramsx.getallparandbnds(praxs1,['p','k','dw','R20500','R2mult'],inclfilt=inclf)
            f=flatten([np.array(e[m][:-1])+np.sum([k for k in [0]+[e[j][-1] for \
                j,i in enumerate(e) if j < len(e)-1]][0:(l+1)]) for m,l in \
                enumerate(np.arange(len(e)))])
            filt.append(f)
        timedat=[flatten(timedat,levels=1)[i] for i in seldatasets]
        rawdata=[flatten(rawdata,levels=1)[i] for i in seldatasets]
        errd=[flatten(errd,levels=1)[i] for i in seldatasets]
        field=[flatten(field,levels=1)[i] for i in seldatasets]
        field2=[flatten(field2,levels=1)[i] for i in seldatasets]
        field3=[flatten(field3,levels=1)[i] for i in seldatasets]
        tr=[flatten(tr,levels=1)[i] for i in seldatasets]
        
        equationtype=[flatten(equationtype,levels=1)[i] for i in seldatasets]
        #paramsx.getallparandbnds()
        a,b1,b2,c,e=paramsx.getallparandbnds(praxs1,['p','k','dw','R20500','R2mult'],inclfilt=[])
        f=flatten([np.array(e[m][:-1])+np.sum([k for k in [0]+[e[j][-1] for \
            j,i in enumerate(e) if j < len(e)-1]][0:(l+1)]) for m,l \
            in enumerate(np.arange(len(e)))])
        poscolll.append(poscoll)
        if drawonly == 0:
            sojetzt,allrescoll=fitcpmg4(praxs1,timedat,rawdata,field,errd,\
                precalc,equationtype,conditions,savstatdir,filenamsav,resnaml, poscolll,\
                moreconditions,paramsx,filt)
        else:
            dwbsetp=[]
            setdwb=0
            for i in np.arange(len(timedat)):
                for j in np.arange(len(flatten(timedat[i]))):
                    dwbsetp.append(setdwb)
                setdwb+=1
            par7=np.array(np.array(dwbsetp).astype('int'))
            fittedcurve,chsq0,chsq1,chsq2,chsq3=printrd(praxs1,timedat,\
                            rawdata,field,errd,equationtype,par7,paramsx,filt)
       #     print fittedcurve, 'fittedcurve'
            for j,i in enumerate(seldatasets):
                print i, 'datasetno', chsq0[j], equationtype[j][0]
        if drawonly == 0:
            try:
                resultcolll.append(sojetzt)
                allresultcoll.append(allrescoll)
            except:
                print 'ugh1'
        else:
            for j in poscoll:
                for k,i in enumerate(seldatasets): #fittedcurve:
                    spinsystems[j].datasets[i].fit=fittedcurve[k]
    if drawonly == 0:
        return resultcolll, resnaml, poscolll, spinsystems, allresultcoll
    else:
        return spinsystems#else:


def parallelfit3(path2020,savstatdir,setparameters3,runnum,paramsx,praxis):
    time.sleep(runnum/10)
    databasis,datapath,selresidues,conditions,filenamsav,precalc2=setparameters3
    os.chdir(datapath)
    resultcoll=[];poscoll=[];resnam=[]#;precalc2=[]
    resultcoll, resnam, poscoll, dataset,allsets = runfit4(praxis,[],\
    selresidues,precalc2,resnam,conditions,path2020,savstatdir,databasis,filenamsav+str(runnum),paramsx,0,0)
    comment1='fitengines1_1';comment2='x';datver='x';scriptver='x'
    hkio.savstatus2(savstatdir,filenamsav+str(runnum),comment1,comment2,datver,scriptver,precalc2,resnam,poscoll,resultcoll,setparameters3,allsets,dataset)
    allsav=[filenamsav,comment1,comment2,datver,scriptver,precalc2,resnam,poscoll,resultcoll,setparameters3,allsets,dataset]
    return allsav


def evaluaterdfit(path2020,savstatdir,praxs1,setparameters2,runnum,paramsx,precalc2):
    """Reads in an experimental data set and recalculates the theoretical
    result from provided parameters. This theoretical result is obtained alongside
    chi square (which is printed during the process for inspection).
    The reuturned spin system contains dataset with the object 'fit' which are
    the calculated theoretical results, for example used to calculate residuals.
    """
    time.sleep(runnum/10)
    databasis,datapath,selresidues,conditions,filenamsav=setparameters2
    os.chdir(datapath)
    resnam=[]
    temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11\
    ,temp12,temp13,temp14,temp15,temp16,temp17,temp18,cond=hkio.\
    loadeverything(savstatdir,[filenamsav],0,decoupl=0)
    print 'pointa'
    ss=runfit4(praxs1,[],selresidues,precalc2,resnam,conditions,path2020,savstatdir,\
    databasis,filenamsav+str(runnum),paramsx,1,cond)
    return ss
    #return spinsystems,fittedcurve

def evaluaterdfitold(path2020,savstatdir,praxs1,setparameters2,runnum,paramsx,precalc2):
    time.sleep(runnum/10)
    databasis,datapath,selresidues,conditions,filenamsav=setparameters2
    #selresidues=[]
    os.chdir(datapath)
    resultcoll=[];poscoll=[];resnam=[]
    temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11,temp12,temp13,temp14,temp15,temp16,temp17,temp18,cond=hkio.loadeverything(savstatdir,[filenamsav],0,decoupl=0)
   # print cond
    ss=runfit4(praxs1,[],selresidues,precalc2,resnam,conditions,path2020,savstatdir,databasis,filenamsav+str(runnum),paramsx,1,cond)
    
    return ss   
   