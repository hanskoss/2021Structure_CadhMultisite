#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 17:49:04 2020

@author: hanskoss
"""
from __future__ import division
import numpy as np
from sympy import flatten
import matplotlib.pyplot as plt
#from random import *
import scipy as scp


def twolines(params,x):
    return np.array([params[0]*(i)+params[1] if i <= params[4] else params[2]*(i)+params[3] for i in x])

def twolineserr(params,x,exp_data):
    value=twolines(params,x)
    return value-exp_data


def oneline(params,x):
    return np.array([params[0]*(i)+params[1] for i in x])

def onelineerr(params,x,exp_data):
    value=oneline(params,x)
 #   print np.abs(value-exp_data)
    return np.abs(value-exp_data)

 #   print ulist, v
def plotall(spinsystems,ulist,v,titl,showtype,labelstring,legendon,ax,showxaxis):
    ploton = 1
    showtypes2=['cest','cpmg','Rex','titpeaklist']
    colorlist=['black', 'blue','cyan','orange','red','magenta','black','cyan','grey','blue','red','orange','green','magenta','black','cyan','grey','yellow']
    linestylelist=['solid','solid','solid','solid','solid','solid','solid','solid','solid','solid','solid','solid','solid']
    legendlist=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    datasetdeactive=[[],[],[],[],[],[],[],[],[],[],\
                     [],[],[],[],[],[],[],[],[],[],\
                     [],[],[],[],[],[],[],[],[],[],\
                     [],[],[],[],[],[],[],[],[],[],\
                     [],[],[],[],[],[],[],[],[],[],\
                     [],[],[],[],[],[],[],[],[],[],\
                     [],[],[],[],[],[],[],[],[],[],[],[]]
    u=ulist
    if showtype == 5:
        legendpos=0
        x=[]
        y=[]
        yerr1=[]
        yerr2=[]
        fitin=0;yfit=[]
  #      setplist=[]
        #for setpos in np.arange(len(spinsystems[u].datasets)):
        for poslabel in labelstring:
            setpos=poslabel[0];label=poslabel[1]
            
            if setpos != 'x' and spinsystems[u].datasets[setpos].datatype in ['R2','Rex','R20'] and spinsystems[u].datasets[setpos].datatype in showtypes2:
                legendentry=spinsystems[u].datasets[setpos].xlabel
                colorentry='black'#colorlist[legendpos]
                legendpos+=1
                #print spinsystems[u].name, datasetdeactive[setpos]
                if spinsystems[u].datasets[setpos].datathere == 1 and spinsystems[u].name[0] not in datasetdeactive[setpos]:
          #          setplist.append(setpos)
    #                x.append(spinsystems[u].datasets[setpos].xlabel+str(setpos))
                    x.append(label)
                    
                 #   print spinsystems[u].datasets[setpos].yval
                 #   print spinsystems[u].datasets[setpos].reshufy
                    if 'reshufy' in vars(spinsystems[u].datasets[setpos]):
                        y.append(spinsystems[u].datasets[setpos].reshufy)
                    else:
                        y.append(spinsystems[u].datasets[setpos].yval)
                    yerr1.append(spinsystems[u].datasets[setpos].yerr1)
                    yerr2.append(spinsystems[u].datasets[setpos].yerr2)
                    #print 'HHerr', np.average([spinsystems[u].datasets[setpos].yerr1,spinsystems[u].datasets[setpos].yerr2])
                    if 'fit' in vars(spinsystems[u].datasets[setpos]):
                        yfit.append(spinsystems[u].datasets[setpos].fit)
                        fitin=1
                    else:
                        yfit.append([])
                    nshift=round(spinsystems[u].datasets[setpos].nshiftav,1)
            else:
                x.append(label)
                y.append('null')
    #    print y, x
        xpos = [i for i, _ in enumerate(x) if y[i] != 'null']
        y = [_ for i, _ in enumerate(y) if y[i] != 'null']
        #print np.array(y), np.array(yerr1), np.array(yerr2), 'aaa'
        #print 'general SN', np.average([np.array(yerr1),np.array(yerr2)],axis=0)
        if len(yerr1) > 0:
            spinsystems[u].HHrange=np.average([np.array(yerr1),np.array(yerr2)],axis=0)
            
        legendentry=spinsystems[u].datasets[setpos].xlabel
        if fitin == 1:
        #spinsystems[4].datasets[0].fit != []:o
         #   print xpos, y, x
            _=ax.errorbar(xpos,y,yerr=np.array([yerr1,yerr2]),marker='o',color=colorentry,label='calculated',linestyle='None',mfc='none',ms=12,mew=2)   
            _=ax.plot(xpos,yfit,marker='o',color=colorentry,linestyle='None',ms=12,label='experimental')#########,linestyle=linestyleentry)
        else:
            _=ax.bar(xpos,y,yerr=[yerr1,yerr2],capsize=2,linewidth=1)
        ind=np.arange(len(x))
        if showxaxis == 1:
            _=ax.set_xticks(ind)#,tuple(x))
            _=ax.set_xticklabels(x, rotation=90)
        else:
            _=ax.set_xticks([])
        try:
            _=ax.set_ylim([0,np.max(y)*1.1])
        except:
            pass
        _=ax.set_ylabel(r'$R_{ex}$'+' / '+r'$s^{-1}$')

#            _=ax.set_ylabel(r'$R_{ex}$'+' / '+r'$s^{-1}$')
#            _=ax.set_ylim([0,np.max(y)*1.1])
        if legendon == 1:
            plt.legend()
    elif showtype == 1:
        colorlist=['grey', 'cyan','blue','red','red','magenta','black','cyan','grey','blue','red','orange','green','magenta','black','cyan','grey','yellow']
        linestylelist=['dashed','dashed','dashed','solid','solid','solid','solid','solid','solid','solid','solid','solid','solid']
        legendpos=0
        cestrange=[]
        for poslabel in labelstring:
       # for setpos in np.arange(len(spinsystems[u].datasets)):
            setpos=poslabel[0];label=poslabel[1]
            if spinsystems[u].datasets[setpos].datatype in ['cest']:
                
                colorentry=colorlist[legendpos]
                linestyleentry=linestylelist[legendpos]
                legendpos+=1
                if spinsystems[u].datasets[setpos].datathere == 1:
          ####          spinsystems[u].datasets[setpos].nshiftav=spinsystems[u].datasets[setpos].nshiftav+spinsystems[u].datasets[setpos].nshiftcorr
              #      print setpos, spinsystems[u].datasets[setpos].xlabel 
                    legendentry=label#spinsystems[u].datasets[setpos].xlabel
                    x=spinsystems[u].datasets[setpos].fdata
                    if 'reshufy' in vars(spinsystems[u].datasets[setpos]):
                        y=spinsystems[u].datasets[setpos].reshufy
                    else:
                        y=spinsystems[u].datasets[setpos].y
                    ymin=spinsystems[u].datasets[setpos].ymin
                    ymax=spinsystems[u].datasets[setpos].ymax
                    cestrange.append(1/np.average((ymax-ymin)/2))
                    #print legendentry, 'SN CEST', np.average((ymax-ymin)/2)
                    if ploton == 1:
                        if 'fit' in vars(spinsystems[u].datasets[setpos]):
                        #spinsystems[4].datasets[0].fit != []:
                            #plt.errorbar(x,y,yerr=np.transpose(np.array(yerr)),marker='x',color=colorentry,label=legendentry,linestyle='None') 
                            _=ax.errorbar(x,y,yerr=[y-ymin,ymax-y],marker='o',label=None,color=colorentry,linestyle='None')
                            ax.plot(x,spinsystems[u].datasets[setpos].fit,color=colorentry,label=legendentry,linestyle=linestyleentry,marker='None')
                        else:
                        #    ax.errorbar(x,y,yerr=np.transpose(np.array(yerr)),marker='x',color=colorentry,label=legendentry,linestyle=linestyleentry)
                            _=ax.errorbar(x,y,yerr=[y-ymin,ymax-y],marker='o',label=None,color=colorentry)
                            x2=[spinsystems[u].datasets[setpos].nshiftav-(i-spinsystems[u].datasets[setpos].nshiftav) for i in x if i > spinsystems[u].datasets[setpos].nshiftav]
                            y2=[y[j] for j,i in enumerate(x) if i > spinsystems[u].datasets[setpos].nshiftav]
                            _=ax.errorbar(x2,y2,marker='o',label=None,linestyle='dotted',color=colorentry)
                        try:
                            _=ax.errorbar(spinsystems[u].datasets[setpos].xorig,spinsystems[u].datasets[setpos].yorig,marker='None',color=colorentry,linestyle='None')
                            xnull=[];ynull=[]
                            for k in spinsystems[u].datasets[setpos].nullranges:
                                xnull.append([l for mm,l in enumerate(spinsystems[u].datasets[setpos].xorig) if l <= np.float(k[1]) and l >= np.float(k[0])])
                                ynull.append([spinsystems[u].datasets[setpos].yorig[mm] for mm,l in enumerate(spinsystems[u].datasets[setpos].xorig) if l <= np.float(k[1]) and l >= np.float(k[0])])
                            xnull=flatten(xnull,levels=1)
                            ynull=flatten(ynull,levels=1)
                            _=ax.errorbar(xnull,ynull,marker='o',label=None,color=colorentry,linestyle='None')#,linestyle='None')
                            xoth=[l for mm,l in enumerate(spinsystems[u].datasets[setpos].xorig) if l not in flatten(xnull,levels=1) and l not in x]
                            yoth=[spinsystems[u].datasets[setpos].yorig[mm] for mm,l in enumerate(spinsystems[u].datasets[setpos].xorig) if l not in flatten(xnull,levels=1) and l not in x]
                            _=ax.errorbar(xoth,yoth,marker='o',label=None,color=colorentry,linestyle='None')
                            _=ax.plot([spinsystems[u].datasets[setpos].nshiftav,spinsystems[u].datasets[setpos].nshiftav],[0,1],color=colorentry,linestyle='dotted')
                        except:
                            yorig=0
                    
        titl=spinsystems[u].name[0]
        spinsystems[u].cestrange=cestrange
        if showxaxis == 1:
            _=ax.set_xlabel(r'$\Omega_N$'' / ppm')
        _=ax.set_xlim([np.min(x),np.max(x)]) #105,132
        _=ax.set_ylabel('I (CEST) / I (B'+r'$_1$'+'=0 Hz)')
        _=ax.set_ylim([0,1.1])
        try:
            ax.plot([np.min(x),np.max(x)],[1,1],color='red')
        except:
            print 'noplot'
        if legendon == 1:
            _=ax.legend()
    elif showtype == 2:
        colorlist=['grey', 'blue','cyan','orange','red','magenta','black','cyan','grey','blue','red','orange','green','magenta','black','cyan','grey','yellow']
        linestylelist=['dashed','dashed','dashed','solid','solid','solid','solid','solid','solid','solid','solid','solid','solid']
        legendpos=0
      #  colorlist=['blue','red','orange','green','magenta']
        cpmgrange=[]
        for poslabel in labelstring:
            setpos=poslabel[0];label=poslabel[1]
            if spinsystems[u].datasets[setpos].datatype in ['cpmg']:
                
                colorentry=colorlist[legendpos]
                linestyleentry=linestylelist[legendpos]
                legendpos+=1
                if spinsystems[u].datasets[setpos].datathere == 1:
                  #  print setpos, 'setpos', spinsystems[u].datasets[setpos].xlabel
                    legendentry=label#spinsystems[u].datasets[setpos].xlabel
                    x=spinsystems[u].datasets[setpos].fdata
                    if 'reshufy' in vars(spinsystems[u].datasets[setpos]):
                        y=spinsystems[u].datasets[setpos].reshufy
                    else:
                        y=spinsystems[u].datasets[setpos].rcpmg
                    yerr=spinsystems[u].datasets[setpos].rcpmgerr
                    if int(float(spinsystems[u].datasets[setpos].expcond['B0'])) > 890:
                        cpmgrange.append(np.average(yerr)/(np.abs(np.average(y[0:2])-(np.average(y[-3:])))))
                    #ymax=spinsystems[u].datasets[setpos].ymax
                    if 'fit' in vars(spinsystems[u].datasets[setpos]):
                    #spinsystems[4].datasets[0].fit != []:
                        ax.errorbar(x,y,yerr=np.transpose(np.array(yerr)),marker='o',color=colorentry,label=legendentry,linestyle='None')   
                        ax.plot(x,spinsystems[u].datasets[setpos].fit,color=colorentry,linestyle=linestyleentry)
                    else:
                        ax.errorbar(x,y,yerr=np.transpose(np.array(yerr)),marker='o',color=colorentry,label=legendentry,linestyle=linestyleentry)
        #print cpmgrange, 'cpmgrange'
        spinsystems[u].cpmgrange=cpmgrange
        if showxaxis == 1:
            ax.set_xlabel(r'$1/\tau_{cp}$'+' / '+r'$s^{-1}$')
        ax.set_ylabel(r'$R_{cpmg}$'+' / '+r'$s^{-1}$')
        if legendon == 1:
            ax.legend()
    elif showtype == 3:
        cutoff=100
        legendpos=0
     #   colorlist=['blue','red','orange','green','magenta']
     #   print 'a'
        #legendlist=['500 MHz 283K','900 MHz 283K','900 MHz 298K']
        for setpos in np.arange(len(spinsystems[u].datasets)):
            if spinsystems[u].datasets[setpos].datatype in ['r1rho']:
                #legendentry=legendlist[legendpos]
                legendentry=spinsystems[u].datasets[setpos].xlabel
                colorentry=colorlist[legendpos]
                legendpos+=1
                if spinsystems[u].datasets[setpos].datathere == 1:
                    if spinsystems[u].datasets[setpos].r1rthere == 1:
                        fig2=plt.figure(12)
              #          plt.title('R1r off resonance', titl)
                        resangcollx=spinsystems[u].datasets[setpos].resangcoll
                        xminmax=[[i for i,j in enumerate(np.abs(1/(np.tan(resangcollx))) < cutoff) if j>0][0],[i for i,j in enumerate(np.abs(1/(np.tan(resangcollx))) < cutoff) if j>0][-1]]
                        nshift=spinsystems[u].datasets[setpos].nshift
                        fullxrange=(spinsystems[u].datasets[setpos].resomegacoll/spinsystems[u].datasets[setpos].b0field+nshift)
                        cutxrange=fullxrange[xminmax[0]:xminmax[1]]
                        r1rav=spinsystems[u].datasets[setpos].r1rav
                        err1=spinsystems[u].datasets[setpos].err1
                        err2=spinsystems[u].datasets[setpos].err2
                        plt.errorbar(cutxrange,r1rav[xminmax[0]:xminmax[1]],yerr=[err1[xminmax[0]:xminmax[1]],err2[xminmax[0]:xminmax[1]]],marker='o',color=colorentry,label=legendentry)
                        plt.legend() #,label=legendentry errorbar
                        trash=plt.errorbar([nshift,nshift],[0,1.1*np.max(r1rav)],color='orange')
    elif showtype == 7:
        ii=0
        legendpos=0
     #   colorlist=['blue','rot','orange','green']
        legendlist=[['1H lw 283K','15N lw 283K'],['1H lw 298K','15N lw 298K'],['1H lw HC','15N lw HC']]
      #  csps=[]
        for setpos in np.arange(len(spinsystems[u].datasets)):
            
            if spinsystems[u].datasets[setpos].datatype in ['titpeaklist']:
                typ=3
                legendentry=legendlist[legendpos]
                colorentry=colorlist[legendpos]
                if spinsystems[u].datasets[setpos].datathere == 1:
                    av=[[],[],[],[],[],[]]
                    xdatarg=np.argsort(spinsystems[u].datasets[setpos].xdat)
                    xa=spinsystems[u].datasets[setpos].xdat[xdatarg]
                    ya=[spinsystems[u].datasets[setpos].hshifts[xdatarg],spinsystems[u].datasets[setpos].nshifts[xdatarg],spinsystems[u].datasets[setpos].heights[xdatarg],\
                    spinsystems[u].datasets[setpos].volumes[xdatarg],spinsystems[u].datasets[setpos].pw1[xdatarg],spinsystems[u].datasets[setpos].pw2[xdatarg]]
                    xmatchmatrix=[xa==i for i in np.sort(list(set(xa)))]
                    xold=np.sort(list(set(xa)))
                    for h,g in enumerate(ya):
                        for b,a in enumerate(xmatchmatrix):
                            
                            av[h].append(np.average([ya[h][j] for j,i in enumerate(a) if a[j] == 1]))
                    xnew=np.arange(np.min(xold),np.max(xold),500)
                    ynew=[scp.interpolate.griddata(xold,np.array(y),xnew,method='linear') for y in av]
                    fig3=plt.figure(13)
                    plt.title(titl)
                    x=xold##old######xnew
                    y=av[typ] #####ynew[typ]
                    y0=av[0]
                    y1=av[1]
                    if legendpos == 0:
                        plt.errorbar(x,np.array(y)/(10**8)-shift1[zz])
                    else:
                        plt.errorbar(x,np.array(y)/(10**8)-shift2[zz])#-y[0])
                        plt.ylim([0,np.max(y)/(10**8)])
                    plt.xlim([0,30000])
                    
                legendpos+=1
        zz+=1
                        
    elif showtype == 4:
    #    print 'step1'
        plt.figure(166)
        plt.clf()
        
        truncplotmod = 4
        ii=0
        legendpos=0
     #   colorlist=['blue','rot','orange','green']
        legendlist=[['1H lw 283K','15N lw 283K'],['1H lw 298K','15N lw 298K'],['1H lw HC','15N lw HC']]
        legendlist=['283K, low conc', '298K, low conc','283K, high conc']
      #  csps=[]
        for setpos in np.arange(len(spinsystems[u].datasets)):
            if spinsystems[u].datasets[setpos].datatype in ['titpeaklist']:
                typ=3
             #   legendentry=spinsystems[u].datasets[setpos].xlabel
                legendentry=legendlist[legendpos]
                colorentry=colorlist[legendpos]
                legendpos+=1
                if spinsystems[u].datasets[setpos].datathere == 1:
                    av=[[],[],[],[],[],[]]
                    xdatarg=np.argsort(spinsystems[u].datasets[setpos].xdat)
                    xa=spinsystems[u].datasets[setpos].xdat[xdatarg]
                    ya=[spinsystems[u].datasets[setpos].hshifts[xdatarg],spinsystems[u].datasets[setpos].nshifts[xdatarg],spinsystems[u].datasets[setpos].heights[xdatarg],\
                    spinsystems[u].datasets[setpos].volumes[xdatarg],spinsystems[u].datasets[setpos].pw1[xdatarg],spinsystems[u].datasets[setpos].pw2[xdatarg]]
                    xmatchmatrix=[xa==i for i in np.sort(list(set(xa)))]
                    xold=np.sort(list(set(xa)))
                   #
                    for h,g in enumerate(ya):
                        for b,a in enumerate(xmatchmatrix):
                            
                            av[h].append(np.average([ya[h][j] for j,i in enumerate(a) if a[j] == 1]))
                    xnew=np.arange(np.min(xold),np.max(xold),500)
                    ynew=[scp.interpolate.griddata(xold,np.array(y),xnew,method='linear') for y in av]
                    x=xold##old######xnew
                    if truncplotmod < 2:
                       rpres=14200
                    elif truncplotmod == 4:
                        rpres=30000
                    else:
                       rpres=30000
                    refy=spinsystems[u].datasets[setpos].volumes[np.argmin(np.abs(spinsystems[u].datasets[setpos].xdat-rpres))] if spinsystems[u].letter == 'B' else spinsystems[u].datasets[setpos].volumes[np.argmin(np.abs(spinsystems[u].datasets[setpos].xdat-0))]
                    if truncplotmod != 1:# and truncplotmod != 4:
                        y=av[typ]/refy #####ynew[typ]
                    else:
                        y=av[typ]
                    y0=av[0]
                    y1=av[1]
                   #
                        #print newrp
#                        print 'gottilhere'
               #     y=yprev
           #         print 'hello'
                    if (truncplotmod > 0 and legendpos < 3) or truncplotmod < 2:
           #             print 'hullu'
                        plt.figure(13)
                        plt.title(titl)
                        if truncplotmod == 2:
                            legendentry=legendentry+'; I(14.2 kpsi)/I(30 kpsi): '+str(np.round(y[np.argmin(np.abs(x-14200))],3))
                        if truncplotmod == 1:
                            legendentry=legendentry+'; I(14.2 kpsi): '+str(np.round(y[np.argmin(np.abs(x-14200))]/10000000,1))+' 1e7'
                        plt.errorbar(x,y,label=legendentry,color=colorentry)
                        if len(spinsystems[u].datasets[setpos].slopes[0]) > 3:
                            plt.plot([spinsystems[u].datasets[setpos].slopes[0][4],spinsystems[u].datasets[setpos].slopes[0][4]],[0,1],linestyle='dotted',color=colorentry)
                            plt.plot(x,oneline([spinsystems[u].datasets[setpos].slopes[0][2],spinsystems[u].datasets[setpos].slopes[0][3]],x),linestyle='dashed',color=colorentry)
                        plt.plot(x,oneline([spinsystems[u].datasets[setpos].slopes[0][0],spinsystems[u].datasets[setpos].slopes[0][1]],x),linestyle='dashed',color=colorentry)
                        plt.plot([14200,14200],[0,1],linestyle='dotted',color='grey')
                        if truncplotmod == 2:
                            plt.ylabel('I(p) / I(30 kpsi)')
                        elif truncplotmod == 1:
                            plt.ylabel('I(p) / a.u.')
                        elif truncplotmod == 0:
                            plt.ylabel('I(p) / I(0 kpsi)')
                        plt.xlabel('pressure / psi')
                        plt.legend()
                        plt.draw()
                        plt.figure(166)
                        plt.plot(x,np.sqrt((np.array(y0)-y0[0])**2+(0.16*np.array(y1)-0.16*y1[0])**2))
                        plt.draw()
             #           print 'gothere'
                        plt.figure(140+ii)
                        plt.title(legendentry)
                        plt.scatter(y0,y1,c=y)#,yerr=np.transpose(np.array(yerr)),marker='x')
                        
                        
                        xnew2=list(x[0::5])
                        if xnew2[-1] != x[-1]:
                            xnew2.append(x[-1])
                        poslab=[h for h,i in enumerate(x) if i in xnew2]
                        for p,q in enumerate(xnew2):
                            plt.annotate(q, (y0[poslab[p]],y1[poslab[p]]))
                        ii+=1
                        facwide=4
                        xspac=0.7/10
                        yspac=4/10#xspac/0.14
                        xlims=y0; ylims=y1
                        xspac=0.7*facwide
                        yspac=4*facwide
                        xlimpress=[np.average(xlims)+xspac,np.average(xlims)-xspac]
                        ylimpress=[np.average(ylims)+yspac,np.average(ylims)-yspac]
                        try:
                            plt.xlim(xlimpress)
                            plt.ylim(ylimpress)
                        except:
                            pass
                        plt.xlabel('proton shift / ppm')
                        plt.ylabel('nitrogen shift / ppm')
                        plt.legend()
    elif showtype == 99:
        ii=0
        legendpos=0
        av=[]
        std=[]
        label=[]
        nslist=[]
        rglist=[]
        legendentrylist=[]
        conclist=[]
      #  print 'OK5'
        for setpos in np.arange(len(spinsystems[u].datasets)):
            if spinsystems[u].datasets[setpos].datatype in ['cpmg','Rex','cest','peaklist']:
    #                    legendentry=legendlist[setpos]
    #                    legendpos+=1
                #try:
                if spinsystems[u].datasets[setpos].datathere == 1:
                    av.append([spinsystems[u].datasets[setpos].hshiftav,spinsystems[u].datasets[setpos].nshiftav,spinsystems[u].datasets[setpos].intensav])
                    std.append([spinsystems[u].datasets[setpos].hshiftstd,spinsystems[u].datasets[setpos].nshiftstd,spinsystems[u].datasets[setpos].intensstd])
                    nslist.append(spinsystems[u].datasets[setpos].expcond['NS'])
                    rglist.append(spinsystems[u].datasets[setpos].expcond['RG'])
                    conclist.append(spinsystems[u].datasets[setpos].expcond['conc'])
                    legendentrylist.append(spinsystems[u].datasets[setpos].xlabel)
                    plt.figure(160)
                    plt.clf()
                    plt.scatter([i[0] for i in av],[i[1] for i in av],c=[i[2] for i in av])
    #                            plt.xlim([np.average([i[0] for i in av])+0.1,np.average([i[0] for i in av])-0.1])
    #                            plt.ylim([np.average([i[1] for i in av])+0.3,np.average([i[1] for i in av])-0.3])
                    label=legendentrylist
                    for p,q in enumerate(label):
                        plt.annotate(q, (av[p][0],av[p][1]))
                    #xlims.append(np.average(y))
                    plt.xlabel('proton shift /ppm')
                    plt.ylabel('nitrogen shift /ppm')
                    try:
                        plt.xlim(xlimpress)
                        plt.ylim(ylimpress)
                    except:
                        pass
                    plt.title(titl)
                    plt.legend()
                else:
                    pass
    _=plt.tight_layout()
    return ax

def plotelements(spinsystems,selnam,v,showtypes,labelstrings,legendon,segmentswitch,figno):
    ulist=[]
    #selnam=['A29']
    for sn in selnam:
        for u in np.arange(len(spinsystems)):
            if spinsystems[u].name[0] in [sn]:
                ulist.append(u)
    fig=plt.figure(figno)
    _=plt.clf()
    plt.rcParams.update({'font.size': 14})
    axs=[]
   # print ulist
    for lsnum, labelstring in enumerate(labelstrings):
        showtype=showtypes[lsnum]
        axs.append(fig.add_subplot(len(selnam),3,1+3*lsnum))
        if lsnum == len(labelstrings)-1 and segmentswitch == 0:
            _=plotall(spinsystems,ulist[lsnum],v,spinsystems[ulist[lsnum]].name,showtype[0],labelstring[0],legendon[0+3*lsnum],axs[0+3*lsnum],1)
        else:
            _=plotall(spinsystems,ulist[lsnum],v,spinsystems[ulist[lsnum]].name,showtype[0],labelstring[0],legendon[0+3*lsnum],axs[0+3*lsnum],0)
   #     _=plt.tight_layout()
        axs.append(fig.add_subplot(len(selnam),3,2+3*lsnum))
        if lsnum == len(labelstrings)-1 and segmentswitch == 0:
            _=plotall(spinsystems,ulist[lsnum],v,spinsystems[ulist[lsnum]].name,showtype[1],labelstring[1],legendon[1+3*lsnum],axs[1+3*lsnum],1)
        else:
            _=plotall(spinsystems,ulist[lsnum],v,spinsystems[ulist[lsnum]].name,showtype[1],labelstring[1],legendon[1+3*lsnum],axs[1+3*lsnum],0)

   #     _=plt.tight_layout()
        axs.append(fig.add_subplot(len(selnam),3,3+3*lsnum))
        
        if lsnum == len(labelstrings)-1 and segmentswitch == 0:
            _=plotall(spinsystems,ulist[lsnum],v,spinsystems[ulist[lsnum]].name,showtype[2],labelstring[2],legendon[2+3*lsnum],axs[2+3*lsnum],1)
        else:
            _=plotall(spinsystems,ulist[lsnum],v,spinsystems[ulist[lsnum]].name,showtype[2],labelstring[2],legendon[2+3*lsnum],axs[2+3*lsnum],0)
        #  plt.show()
        _=plt.tight_layout()

    


