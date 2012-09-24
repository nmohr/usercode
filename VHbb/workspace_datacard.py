#!/usr/bin/env python
import os
import sys
import ROOT 
from ROOT import TFile
from array import array
from math import sqrt
from copy import copy
#suppres the EvalInstace conversion warning bug
import warnings
warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )
from optparse import OptionParser
from BetterConfigParser import BetterConfigParser
from samplesclass import sample
from mvainfos import mvainfo
import pickle
from progbar import progbar
from printcolor import printc
from gethistofromtree import getHistoFromTree, orderandadd

#CONFIGURE
argv = sys.argv
parser = OptionParser()
parser.add_option("-P", "--path", dest="path", default="",
                      help="path to samples")
parser.add_option("-V", "--var", dest="variable", default="",
                      help="variable for shape analysis")
parser.add_option("-C", "--config", dest="config", default=[], action="append",
                      help="configuration file")
(opts, args) = parser.parse_args(argv)
if opts.config =="":
        opts.config = "config"
print opts.config
config = BetterConfigParser()
config.read(opts.config)
anaTag = config.get("Analysis","tag")



# -------------------- parsing configuration and options: (an ugly spaghetti code section) ----------------------------------------------------------------------
#get locations:
Wdir=config.get('Directories','Wdir')
#systematics
systematics=config.get('systematics','systematics')
systematics=systematics.split(' ')
weightF=config.get('Weights','weightF')
path=opts.path
var=opts.variable
plot=config.get('Limit',var)
infofile = open(path+'/samples.info','r')
info = pickle.load(infofile)
infofile.close()
options = plot.split(',')
if len(options) < 12:
    print "You have to choose option[11]: either Mjj or BDT"
    sys.exit("You have to choose option[11]: either Mjj or BDT")
name=options[1]
title = options[2]
nBins=int(options[3])
xMin=float(options[4])
xMax=float(options[5])
SIG=options[9]
data=options[10]
anType=options[11]
RCut=options[7]
setup=eval(config.get('LimitGeneral','setup'))
ROOToutname = options[6]
outpath=config.get('Directories','limits')
systematicsnaming=eval(config.get('LimitGeneral','systematicsnaming7TeV'))
if anaTag =='8TeV':
    systematicsnaming=eval(config.get('LimitGeneral','systematicsnaming8TeV'))
    MC_rescale_factor=1.0
elif anaTag =='7TeV':
    MC_rescale_factor=2.0
    printc('red','', 'I  RESCALE by 2.0! (from training)')  
else:
    print "What is your Analysis Tag in config? (anaTag)"
    sys.exit("What is your Analysis Tag in config? (anaTag)")
scaling=eval(config.get('LimitGeneral','scaling'))
rescaleSqrtN=eval(config.get('LimitGeneral','rescaleSqrtN'))
if 'RTight' in RCut:
    Datacradbin=options[10]
elif 'RMed' in RCut:
    Datacradbin=options[10]
else:
    Datacradbin=options[10]
blind=eval(config.get('LimitGeneral','blind'))
BKGlist = eval(config.get('LimitGeneral','BKG'))
#Groups for adding samples together
Group = eval(config.get('LimitGeneral','Group'))
#naming for DC
Dict= eval(config.get('LimitGeneral','Dict'))
weightF_sys = eval(config.get('LimitGeneral','weightF_sys'))
binstat = eval(config.get('LimitGeneral','binstat'))
bdt = False
mjj = False
#print str(anType)
#print len(options)
if str(anType) == 'BDT':
    bdt = True
    systematics = eval(config.get('LimitGeneral','sys_BDT'))
elif str(anType) == 'Mjj':
    mjj = True
    systematics = eval(config.get('LimitGeneral','sys_Mjj'))
sys_cut_suffix=eval(config.get('LimitGeneral','sys_cut_suffix'))




# -------------------- generate the Workspace with all Histograms: ----------------------------------------------------------------------

ROOT.gROOT.SetStyle("Plain")
datas = []
datatyps =[]
histos = []
typs = []
if blind: 
    printc('red','', 'I AM BLINDED!')  
counter=0

if weightF_sys:
    #weightF_sys_function=config.get('Weights','weightF_sys')
    weightF_sys_histos = []
    weightF_sys_Ups = []
    weightF_sys_Downs = []

for job in info:
    if eval(job.active):
        if job.subsamples:
            for subsample in range(0,len(job.subnames)):
                if job.subnames[subsample] in BKGlist:
                    hTemp, typ = getHistoFromTree(job,path,config,options,MC_rescale_factor,subsample)
                    histos.append(hTemp)
                    typs.append(Group[job.subnames[subsample]])
                    if weightF_sys:
                        hTempW, _ = getHistoFromTree(job,path,config,options,MC_rescale_factor,subsample,'weightF_sys')
                        weightF_sys_histos.append(hTempW)
                    if counter == 0:
                        hDummy = copy(hTemp)
                    else:
                        hDummy.Add(hTemp)
                    counter += 1
                    
                elif job.subnames[subsample] == SIG:
                    hTemp, typ = getHistoFromTree(job,path,config,options,MC_rescale_factor,subsample)
                    histos.append(hTemp)
                    typs.append(Group[job.subnames[subsample]])
                    if weightF_sys:
                        hTempW, _ = getHistoFromTree(job,path,config,options,MC_rescale_factor,subsample,'weightF_sys')
                        weightF_sys_histos.append(hTempW)
    
        else:
            if job.name in BKGlist:
                #print job.getpath()
                hTemp, typ = getHistoFromTree(job,path,config,options,MC_rescale_factor)
                histos.append(hTemp)
                typs.append(Group[job.name])                        
                if weightF_sys:
                    hTempW, _ = getHistoFromTree(job,path,config,options,MC_rescale_factor,-1,'weightF_sys')
                    weightF_sys_histos.append(hTempW)

                if counter == 0:
                    hDummy = copy(hTemp)
                else:
                    hDummy.Add(hTemp)
                counter += 1
                
            elif job.name == SIG:
                hTemp, typ = getHistoFromTree(job,path,config,options,MC_rescale_factor)
                histos.append(hTemp)
                typs.append(Group[job.name])                                        
                if weightF_sys:
                    hTempW, _ = getHistoFromTree(job,path,config,options,MC_rescale_factor,-1,'weightF_sys')
                    weightF_sys_histos.append(hTempW)

            elif job.name in data:
                #print 'DATA'
                hTemp, typ = getHistoFromTree(job,path,config,options)
                datas.append(hTemp)
                datatyps.append(typ)

MC_integral=0
MC_entries=0
for histo in histos:
    MC_integral+=histo.Integral()
printc('green','', 'MC integral = %s'%MC_integral)

#order and add together
typs2=copy(typs)
histos, typs = orderandadd(histos,typs,setup)
if weightF_sys:
	weightF_sys_histos,_=orderandadd(weightF_sys_histos,typs2,setup)

#SYSTEMATICS:
UD = ['Up','Down']
systhistosarray=[]
Coco=0 #iterates over (all systematics) * (up,down)
nominalShape = options[0]

outfile = ROOT.TFile(outpath+'vhbb_TH_'+ROOToutname+'.root', 'RECREATE')

for sys in systematics:
    for Q in UD: # Q = 'Up' and 'Down'
        #options[7] ist der CutString name
        new_cut=sys_cut_suffix[sys]
        new_options = copy(options)
        if not new_cut == 'nominal':
            old_str,new_str=new_cut.split('>')
            new_options[7]=[options[7],old_str,new_str.replace('?',Q)]
        ff=options[0].split('.')
        if bdt == True:
            ff[1]='%s_%s'%(sys,Q.lower())
            new_options[0]=nominalShape.replace('.nominal','.%s_%s'%(sys,Q.lower()))
        elif mjj == True:
            ff[0]='H_%s'%(sys)
            ff[1]='mass_%s'%(Q.lower())
            new_options[0]='.'.join(ff)

        print '\n'
        printc('blue','','\t--> doing systematic %s %s'%(sys,Q.lower())) 

        systhistosarray.append([])
        typsX = []

        for job in info:
            if eval(job.active):
                if job.subsamples:
                    for subsample in range(0,len(job.subnames)):
                        if job.subnames[subsample] in BKGlist:
                            hTemp, typ = getHistoFromTree(job,path,config,new_options,MC_rescale_factor,subsample)
                            systhistosarray[Coco].append(hTemp)
                            typsX.append(Group[job.subnames[subsample]])
                        elif job.subnames[subsample] == SIG:
                            hTemp, typ = getHistoFromTree(job,path,config,new_options,MC_rescale_factor,subsample)
                            systhistosarray[Coco].append(hTemp)
                            typsX.append(Group[job.subnames[subsample]])
                            
                else:
                    if job.name in BKGlist:
                        hTemp, typ = getHistoFromTree(job,path,config,new_options,MC_rescale_factor)
                        systhistosarray[Coco].append(hTemp)
                        typsX.append(Group[job.name])
                    elif job.name == SIG:
                        hTemp, typ = getHistoFromTree(job,path,config,new_options,MC_rescale_factor)
                        systhistosarray[Coco].append(hTemp)
                        typsX.append(Group[job.name])

        MC_integral=0
        for histoX in systhistosarray[Coco]:
            MC_integral+=histoX.Integral()
        printc('green','', 'MC integral = %s'%MC_integral)  
        systhistosarray[Coco], typsX = orderandadd(systhistosarray[Coco],typsX,setup)

        if scaling: #rescaling after the sys has been propagated through the BDT with a scaling
            for i in range(0,len(systhistosarray[Coco])):
                for bin in range(0,histos[i].GetSize()):
                    A=systhistosarray[Coco][i].GetBinContent(bin)
                    B=histos[i].GetBinContent(bin)
                    systhistosarray[Coco][i].SetBinContent(bin,B+((A-B)/4.))
        Coco+=1

nToys = 500
for nToy in range(nToys):
	WS = ROOT.RooWorkspace('%s_Toy%s'%(Datacradbin,nToy),'%s_Toy%s'%(Datacradbin,nToy)) #Zee
	print 'WS initialized'
	disc= ROOT.RooRealVar(name,name,xMin,xMax)
	obs = ROOT.RooArgList(disc)
	statUps=[]
	statDowns=[]
	for i in range(0,len(histos)):
	    newname=Dict[typs[i]]
	    histos[i].SetName(newname)
	    #histos[i].SetDirectory(outfile)
	    outfile.cd()
	    histos[i].Write()
	    errorsum=0
	    total=0
	    for j in range(histos[i].GetNbinsX()+1):
	        errorsum=errorsum+(histos[i].GetBinError(j))**2
	    errorsum=sqrt(errorsum)
	    total=histos[i].Integral()
	
	    if binstat: #treating statistics in single bins
	        for bin in range(0,nBins):
	            statUps.append(histos[i].Clone())
	            statDowns.append(histos[i].Clone())
	            statUps[i*nBins+bin].SetName('%sCMS_vhbb_stats_%s_%s_%sUp'%(newname,newname,bin,options[10]))
	            statDowns[i*nBins+bin].SetName('%sCMS_vhbb_stats_%s_%s_%sDown'%(newname,newname,bin,options[10]))
	            #shift up and down with statistical error
	            if rescaleSqrtN:
	                statUps[i*nBins+bin].SetBinContent(bin,statUps[i*nBins+bin].GetBinContent(bin)+statUps[i*nBins+bin].GetBinError(bin)/total*errorsum)
	                statDowns[i*nBins+bin].SetBinContent(bin,statDowns[i*nBins+bin].GetBinContent(bin)-statDowns[i*nBins+bin].GetBinError(bin)/total*errorsum)
	            else:
	                statUps[i*nBins+bin].SetBinContent(bin,statUps[i*nBins+bin].GetBinContent(bin)+statUps[i*nBins+bin].GetBinError(bin))
	                statDowns[i*nBins+bin].SetBinContent(bin,statDowns[i*nBins+bin].GetBinContent(bin)-statDowns[i*nBins+bin].GetBinError(bin))
	            statUps[i*nBins+bin].Write()
	            statDowns[i*nBins+bin].Write()
	            histPdf = ROOT.RooDataHist(newname,newname,obs,histos[i])
	            #UP stats of MCs
	            RooStatsUp = ROOT.RooDataHist('%sCMS_vhbb_stats_%s_%s_%sUp'%(newname,newname,bin,options[10]),'%sCMS_vhbb_stats_%s_%s_%sUp'%(newname,newname,bin,options[10]),obs, statUps[i*nBins+bin])
	            #DOWN stats of MCs
	            RooStatsDown = ROOT.RooDataHist('%sCMS_vhbb_stats_%s_%s_%sDown'%(newname,newname,bin,options[10]),'%sCMS_vhbb_stats_%s_%s_%sDown'%(newname,newname,bin,options[10]),obs, statDowns[i*nBins+bin])
	            getattr(WS,'import')(histPdf)
	            getattr(WS,'import')(RooStatsUp)
	            getattr(WS,'import')(RooStatsDown)
	
	    else:
	        statUps.append(histos[i].Clone())
	        statDowns.append(histos[i].Clone())
	        statUps[i].SetName('%sCMS_vhbb_stats_%s_%sUp'%(newname,newname,options[10]))
	        statDowns[i].SetName('%sCMS_vhbb_stats_%s_%sDown'%(newname,newname,options[10]))
	        #shift up and down with statistical error
	        for j in range(histos[i].GetNbinsX()+1):
	            if rescaleSqrtN:
	                statUps[i].SetBinContent(j,statUps[i].GetBinContent(j)+statUps[i].GetBinError(j)/total*errorsum)
	                statDowns[i].SetBinContent(j,statDowns[i].GetBinContent(j)-statDowns[i].GetBinError(j)/total*errorsum)
	            else:
	                statUps[i].SetBinContent(j,statUps[i].GetBinContent(j)+statUps[i].GetBinError(j))
	                statDowns[i].SetBinContent(j,statDowns[i].GetBinContent(j)-statDowns[i].GetBinError(j))
	        statUps[i].Write()
	        statDowns[i].Write()
	        histPdf = ROOT.RooDataHist(newname,newname,obs,histos[i])
	        #UP stats of MCs
	        RooStatsUp = ROOT.RooDataHist('%sCMS_vhbb_stats_%s_%sUp'%(newname,newname,options[10]),'%sCMS_vhbb_stats_%s_%sUp'%(newname,newname,options[10]),obs, statUps[i])
	        #DOWN stats of MCs
	        RooStatsDown = ROOT.RooDataHist('%sCMS_vhbb_stats_%s_%sDown'%(newname,newname,options[10]),'%sCMS_vhbb_stats_%s_%sDown'%(newname,newname,options[10]),obs, statDowns[i])
	        getattr(WS,'import')(histPdf)
	        getattr(WS,'import')(RooStatsUp)
	        getattr(WS,'import')(RooStatsDown)
	        
	    #And now WeightF sys
	    if weightF_sys:
	        weightF_sys_Downs.append(weightF_sys_histos[i].Clone())
	        weightF_sys_Ups.append(weightF_sys_histos[i].Clone())
	        weightF_sys_Downs[i].SetName('%sCMS_vhbb_weightF_%sDown'%(newname,options[10]))
	        weightF_sys_Ups[i].SetName('%sCMS_vhbb_weightF_%sUp'%(newname,options[10]))
	        for j in range(histos[i].GetNbinsX()+1):
	            weightF_sys_Ups[i].SetBinContent(j,2*histos[i].GetBinContent(j)-weightF_sys_Downs[i].GetBinContent(j))
	        weightF_sys_Ups[i].Write()
	        weightF_sys_Downs[i].Write()    
	        RooWeightFUp = ROOT.RooDataHist('%sCMS_vhbb_weightF_%sUp'%(newname,options[10]),'%sCMS_vhbb_weightF_%s_%sUp'%(newname,newname,options[10]),obs, weightF_sys_Ups[i])
	        RooWeightFDown = ROOT.RooDataHist('%sCMS_vhbb_weightF_%sDown'%(newname,options[10]),'%sCMS_vhbb_weightF_%s_%sDown'%(newname,newname,options[10]),obs, weightF_sys_Downs[i])
	        getattr(WS,'import')(RooWeightFUp)
	        getattr(WS,'import')(RooWeightFDown)
	
	
	#HISTOGRAMM of DATA    
	d1 = ROOT.TH1F('d1','d1',nBins,xMin,xMax)
	for i in range(0,len(datas)):
	    d1.Add(datas[i],1)
	printc('green','','\nDATA integral = %s\n'%d1.Integral())
	flow = d1.GetEntries()-d1.Integral()
	if flow > 0:
	    printc('red','','U/O flow: %s'%flow)
	d1.SetName(Dict['Data'])
	outfile.cd()
	d1.Write()
	
	if blind:
	    hDummy.SetName(Dict['Data'])
	    histPdf = ROOT.RooDataHist(Dict['Data'],Dict['Data'],obs,hDummy)
	    rooDummy = ROOT.RooDataHist('data_obs','data_obs',obs,hDummy)
	    toy = ROOT.RooHistPdf('data_obs','data_obs',ROOT.RooArgSet(obs),rooDummy)
	    rooDataSet = toy.generate(ROOT.RooArgSet(obs),int(hDummy.Integral()))
	    histPdf = ROOT.RooDataHist('data_obs','data_obs',ROOT.RooArgSet(obs),rooDataSet.reduce(ROOT.RooArgSet(obs)))
	else:
	    histPdf = ROOT.RooDataHist(Dict['Data'],Dict['Data'],obs,d1)
	#ROOT.RooAbsData.plotOn(histPdf,frame)
	getattr(WS,'import')(histPdf)
        Coco = 0
        for sys in systematics:
	        for Q in UD: # Q = 'Up' and 'Down'
	            # finaly lpop over histos
	            for i in range(0,len(systhistosarray[Coco])):
	                systhistosarray[Coco][i].SetName('%s%s%s'%(Dict[typs[i]],systematicsnaming[sys],Q))
	                outfile.cd()
	                systhistosarray[Coco][i].Write()            
	                histPdf = ROOT.RooDataHist('%s%s%s'%(Dict[typs[i]],systematicsnaming[sys],Q),'%s%s%s'%(Dict[typs[i]],systematicsnaming[sys],Q),obs,systhistosarray[Coco][i])
	                getattr(WS,'import')(histPdf)
                Coco+=1
	
	    
	WS.writeToFile(outpath+'vhbb_WS_'+ROOToutname+'_Toy%s.root'%nToy)
	
	
	
	# -------------------- write DATAcard: ----------------------------------------------------------------------
	columns=len(setup)
	
	if '8TeV' in options[10]:
	    pier = open(Wdir+'/pier8TeV.txt','r')
	else:
	    pier = open(Wdir+'/pier.txt','r')
	scalefactors=pier.readlines()
	pier.close()
	f = open(outpath+'vhbb_DC_'+ROOToutname+'_Toy'+str(nToy)+'.txt','w')
	f.write('imax\t1\tnumber of channels\n')
	f.write('jmax\t%s\tnumber of backgrounds (\'*\' = automatic)\n'%(columns-1))
	f.write('kmax\t*\tnumber of nuisance parameters (sources of systematical uncertainties)\n\n')
	if bdt==True:
	    f.write('shapes * * vhbb_WS_%s_Toy%s.root $CHANNEL:$PROCESS $CHANNEL:$PROCESS$SYSTEMATIC\n\n'%(ROOToutname,nToy))
	else:
	    f.write('shapes * * vhbb_TH_%s.root $PROCESS $PROCESS$SYSTEMATIC\n\n'%ROOToutname)
	f.write('bin\t%s_Toy%s\n\n'%(Datacradbin,nToy))
	if blind:
	    f.write('observation\t%s\n\n'%(int(hDummy.Integral())))
	else:
	    f.write('observation\t%s\n\n'%(int(d1.Integral())))
	
	f.write('bin')
	for c in range(0,columns): f.write('\t%s_Toy%s'%(Datacradbin,nToy))
	f.write('\n')
	
	f.write('process')
	for c in setup: f.write('\t%s'%Dict[c])
	f.write('\n')
	
	f.write('process')
	for c in range(0,columns): f.write('\t%s'%c)
	f.write('\n')

	f.write('rate')
	for c in range(0,columns): f.write('\t%s'%histos[c].Integral())
	f.write('\n')

	InUse=eval(config.get('Datacard','InUse'))
	#Parse from config
	for item in InUse:
	    f.write(item)
	    what=eval(config.get('Datacard',item))
	    f.write('\t%s'%what['type'])
	    for c in setup:
	        if c in what:
	            if item == 'CMS_eff_e' and 'Zmm' in options[10]: f.write('\t-')
	            elif item == 'CMS_eff_m' and 'Zee' in options[10]: f.write('\t-')
	            elif item == 'CMS_trigger_e' and 'Zmm' in options[10]: f.write('\t-')
	            elif item == 'CMS_trigger_m' and 'Zee' in options[10]: f.write('\t-')
	            else:
	                f.write('\t%s'%what[c])
	        else:
	            f.write('\t-')
	    f.write('\n')
	
	#Write shape stats and sys
	if binstat:
	    for c in setup:
	        for bin in range(0,nBins):
	            f.write('CMS_vhbb_stats_%s_%s_%s\tshape'%(Dict[c], bin, options[10]))
	            for it in range(0,columns):
	                if it == setup.index(c):
	                    f.write('\t1.0')
	                else:
	                    f.write('\t-')
	            f.write('\n')
	
	else:
	    for c in setup:
	        f.write('CMS_vhbb_stats_%s_%s\tshape'%(Dict[c], options[10]))
	        for it in range(0,columns):
	            if it == setup.index(c):
	                f.write('\t1.0')
	            else:
	                f.write('\t-')
	        f.write('\n')
	    
	if weightF_sys:
	    f.write('CMS_vhbb_weightF_%s\tshape'%(options[10]))
	    for it in range(0,columns): f.write('\t1.0')
	    f.write('\n')
	
	    
	    
	if scaling: sys_factor=0.25
	else: sys_factor=1.0
	for sys in systematics:
	    f.write('%s\tshape'%systematicsnaming[sys])
	    for c in range(0,columns): f.write('\t%s'%sys_factor)
	    f.write('\n')
	f.close()
	outfile.Close()
