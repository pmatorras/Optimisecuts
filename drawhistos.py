from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TLegend
from ROOT import gROOT,gDirectory, TTree
from collections import OrderedDict
import os, optparse

optim="/home/pablinux/cernbox/www/susy/optimisation"
hfolder='../Histograms/'
hfilenm="../Output/allpostcuts.root"
hfile  = TFile(hfilenm,"READ","Example");
os.system("mkdir -p "+hfolder)

#set x and y ranges
def setRanges(hT2tt,httbar):
    xMin=-999
    xMax= 999
    yMax= 999
    yMin= 0
    #    histo=
    yMax=1.1*max(hT2tt.GetMaximum(),httbar.GetMaximum())
    xMin=min(hT2tt.GetXaxis().GetXmin(),httbar.GetXaxis().GetXmin())
    xMax= max(hT2tt.GetXaxis().GetXmax(),httbar.GetXaxis().GetXmax())
    yMin=0.9*min(hT2tt.GetMinimum(),httbar.GetMinimum())
    hT2tt.GetYaxis().SetRangeUser(yMin,0.4)#yMax)
    #hT2tt.GetXaxis().SetRangeUser(xMin,xMax)
    #print xMax, xMin, "-----------------------------"
#Scale histograms
def ScaleToMax(histos,varControl):
    for histo in histos:
        if varControl is True: continue
        norm=1
        maxBin=histo.GetMaximumBin()
        norm=histo.GetBinContent(maxBin)
        histo.Scale(1/norm)

#Scale histograms
def OLDScaleToInt(histos,varControl):
    for histo in histos:
        norm=histo.GetSumOfWeights()
        if(norm<0.0001): norm=1.0
        if varControl is False: histo.Scale(1/norm)
        if norm<10: histo.Scale(0.1)
        print "SCALING:\n",histo.GetName(), norm, histo.GetEntries()
#Scale histograms
def ScaleToInt(histos,varControl):
    for sample in histos:
        histo=gDirectory.Get(histos[sample])
        #histo=histos[sample]
        print sample
        norm=histo.GetSumOfWeights()
        if(norm<0.0001): norm=1.0
        if varControl is False: histo.Scale(1/norm)
        if norm<10: histo.Scale(0.1)
        print "NORM", norm
#check if histogram exists
def OLDhistoexists(histo, hnm):
    cont=False
    histo_type=str(type(histo)) #'ROOT.TH1F'>"
    if (bool('ROOT.TH1F' not in histo_type)):
        print "-----------------------------------"
        print "|  "+hnm+" not a histogram\t  |"
        print "| perhaps variable doesn't exist? |"
        print "-----------------------------------"
        cont=True
    return cont
#check if histogram exists
def histoexists(histo, hnm):
    cont=False
    histo_type=str(type(histo)) #'ROOT.TH1F'>"
    if (bool('ROOT.TH1F' not in histo_type)):
        print "-----------------------------------"
        print "|  "+hnm+" not a histogram\t  |"
        print "| perhaps variable doesn't exist? |"
        print "-----------------------------------"
        cont=True
    return cont

#Scale x axis
def findXrange(hdict, binWOri,nBinFin,xMinOri):
    xMinBin=999999
    xMaxBin=0
    for samp in hdict:
        histo=gDirectory.Get(hdict[samp])
        xMinBin=min(xMinBin,histo.FindFirstBinAbove(0.1))
        xMaxBin=max(xMaxBin,histo.FindLastBinAbove(0.1))
    if(xMinBin>xMaxBin):
        xMinBin=0
        xMaxBin=nBinFin
    nBin=xMaxBin-xMinBin
    xLen=nBin*binWOri
    widthNew=xLen/nBinFin
    scaleFactor=int(widthNew/binWOri)
    while(scaleFactor>1 and nBinOri%scaleFactor>0): scaleFactor-=1
    if(scaleFactor<1): scaleFactor=1
    xMax=xMinOri+xMaxBin*binWOri
    xMin=xMinOri+xMinBin*binWOri   
    if(xMax>200): xMin=0
    return xMin, xMax, scaleFactor
#Scale y axis
def findYrange(hdict, scaleFactor):
    yMax=0
    for samp in hdict:
        histo=gDirectory.Get(hdict[samp])                
        histo.Rebin(scaleFactor)
        ymaxBin=histo.GetMaximumBin()
        yMax=max(yMax,histo.GetBinContent(ymaxBin))
    return yMax
            

#Optional parameters
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('--test'   , dest='test'   , help='run only few plots', default=False, action='store_true')
parser.add_option('--rmold'  , dest='rmold'  , help='remove old images' , default=False, action='store_true')
parser.add_option('--nobatch', dest='nobatch', help='remove old images' , default=False, action='store_true')
parser.add_option('--ISR'    , dest='ISR'   , help='do only ISR histograms'  , default=False, action='store_true')
parser.add_option('--btag'   , dest='btag'  , help='do only btag histograms' , default=False, action='store_true')
parser.add_option('--bveto'  , dest='bveto' , help='do only bveto histograms', default=False, action='store_true')
parser.add_option('--lowMT2' , dest='lowMT2', help='only draw histos for low MT2', default=False, action='store_true')
parser.add_option('--var'    , dest='var', help='only draw histos for sel var', default=False)

parser.add_option('--all', dest='allSamples', help='do all histograms', default=False, action='store_true')
(opt, args) = parser.parse_args()
if opt.var is not False:
    print "opt", opt.var
#Apply options
allsamples=opt.allSamples

if opt.rmold is True:
    print "removing histograms in", hfolder
    os.system('rm -r '+hfolder)
gROOT.SetBatch(not(opt.nobatch))

#Define variables to over which make histograms
allvars=["evid_", "isSF_", "btagW_", "bvetoW_", "nSV_", "nLepton_", "ISRcut_","nbCleanjets_","nbjets_" ,\
         "mll_", "mt2ll_", "detall_","dRll_", "MET_sumEt_","MET_pt_", "ptmiss_","lep1_pt_","lep2_pt_" ,\
         "jet1_pt_","jet2_pt_","dphill_", "dphijj_","detajj_","dRjj_","Dnjetstot_","Dnbjetstot_",\
         "SV_x_", "SV_y_", "SV_z_", "SV_chi2_","SV_eta_", "SV_phi_", "SV_pt_", "SV_mass_",\
         "PV_x_", "PV_y_", "PV_z_", "PV_chi2_","PV_npvs_","susyMstop_","susyMLSP_"]
           
controlvars = ["evid_","btagW_","bvetoW_", "susyMstop_","susyMLSP_","isSF_", "nLepton_"]
flavours= {'df':'-1','sf':'1'}
bjets   = {'btag': 'btagW_','bveto':'bvetoW_'}
samples = {'T2tt':'hT2tt', 'ttbar':'httbar', 'WW': 'hWW'}
massrange='&& susyMstop>=400 && susyMstop<700 '
deltamass=' && susyMstop-susyMLSP'
mStops = {'mS-400to700_dm-1to125'   :[massrange + deltamass +'<=125', 6],\
          'mS-400to700_dm-125to200' :[massrange + deltamass +'>125'+ deltamass +'<=200', 30],\
          'mS-400to700_dm-200to700' :[massrange + deltamass +'>200', 46],\
        }
testvars=["njets_", "nbjets_", "nCleanjets_", "nbCleanjets_",\
          "dphil1jmin_","dphil1jmax_","dphil1bmin_","dphil1bmax_",\
          "dphil2jmin_","dphil2jmax_","dphil2bmin_","dphil2bmax_",\
          "bjet1_pt_","bjet1_eta_", "bjet1_phi_", "bjet2_pt_", "bjet2_eta_", "bjet2_phi_"]


bkgs={'ttbar': 4,'WW': 3}
ISRcut={"ISR": "&& ISRcut>0", "Inclusive":"", "lowMT2": " && mt2ll<20", "lowMT2_ISR": " && mt2ll<20 && ISRcut>0"}
print "Creating histograms:"
c1 = TCanvas( 'c1', 'Dynamic Filling Example', 200,10, 1600, 900 )
if  opt.test is True:
    variables=testvars
else:
    variables=allvars
for var in variables:
    #Save in different folder if control variable
    cfolder=''
    isControl=bool(var in controlvars)        
    if(isControl is True): cfolder='Control/'

    if opt.var not in var: continue
    print "var", var
    btagVeto=False
    if(var in bjets.values()): btagVeto=True 
    for bjet in bjets:
        if allsamples is False and opt.btag  is True and "btag"  not in bjet: continue  
        if allsamples is False and opt.bveto is True and "bveto" not in bjet: continue
        if btagVeto is True and bjet not in var:
            print "Omitting histogram since var="+var+" and bjet="+bjet
            continue
        for isr in ISRcut: 
            if allsamples is False and opt.ISR is False and "ISR"     in isr: continue
            if allsamples is False and opt.ISR is True  and "ISR" not in isr: continue
            if allsamples is False and opt.lowMT2 is False and "lowMT"     in isr: continue
            if allsamples is False and opt.lowMT2 is True  and "lowMT" not in isr: continue

            samplefolder=hfolder+bjet+'/'+isr
            isrnm=isr
            if isrnm in "Inclusive": isrnm="All"
            os.system('mkdir -p '+samplefolder)
            for fl in flavours:
                histnm = var+bjet+'_'+fl+'_'+isrnm
                nBinOri= 1000000
                nBinFin= 10#40
                xMinOri=-5000.
                xMaxOri= 5000.
                lenOri= (xMaxOri-xMinOri)
                binWOri=lenOri/nBinOri
                xRanOri= "("+str(nBinOri)+","+str(xMinOri)+","+str(xMaxOri)+")"
                dtest=OrderedDict({})
                for sample in samples:
                    #make histo for each variable and condition
                    fillvar   = '>-999'
                    if var in ['jet1_pt', 'jet2_pt', 'mt2ll']: fillvar='>0'
                    tree      = hfile.Get(sample)
                    histo     = samples[sample]+histnm
                    treevar   = var+sample+">>"+histo+xRanOri
                    flcond    = "isSF_"+sample+"=="+flavours[fl]
                    isrcond   = ISRcut[isr]
                    number ='-999'
                    if "dphi" in var: number='-1'
                    onlyfill  = var+sample+">"+number
                    bjetcond  = bjets[bjet]+sample
                    condition = bjetcond+'*('+flcond+'&&'+onlyfill+isrcond+')'
                    tree.Draw(treevar, condition)
                    dtest[sample]=histo
                    #divide in different mass ranges
                    for mstop in mStops:
                        if sample not in 'T2tt': continue
                        splitcond    = condition.split(')')
                        mStopcond    = splitcond[0]+mStops[mstop][0]+')'
                        treemStop = var+sample+'>>'+histo+'_'+mstop+xRanOri
                        tree.Draw(treemStop,mStopcond)
                        dtest[sample+'_'+mstop]=histo+'_'+mstop
            
                #Set proper ranges and rebin histogram
                xMin, xMax, scaleFactor= findXrange(dtest, binWOri, nBinFin, xMinOri)
                ScaleToInt(dtest, isControl)
                print "scaleFactor", scaleFactor
                yMax=findYrange(dtest, scaleFactor)
                legend = TLegend(0.7,0.7,0.9,0.9);

            
                for idx,sam in enumerate(dtest):
                    isSame='same'
                    hist=gDirectory.Get(dtest[sam])
                    if idx is 0:
                        hist.GetXaxis().SetRangeUser(xMin,xMax)
                        hist.GetYaxis().SetRangeUser(0,1.1*yMax)
                        hist.SetTitle(histnm)
                        hist.SetStats(False)

                        isSame=''
                        hist.Draw('hist')
                    if sam in "T2tt":
                        hist.Draw('hist same')
                        hist.SetLineColor(2)
                        #print "histo", hist, xMin, xMax
                    for ms in mStops:
                        if ms in sam :
                            #print "colors", type( hist), hist.GetName()
                            hist.SetLineColor(mStops[ms][1])
                            hist.Draw('same')
                    isBkg=''
                    if sam in bkgs:
                        isBkg='hist '
                        hist.SetLineColor(bkgs[sam])
                    hist.Draw(isBkg+'same')
                    legend.AddEntry(hist,sam,"f");
                legend.Draw()
                #save to image
                outputfolder=samplefolder+'/'+cfolder
                os.system('mkdir -p '+outputfolder)
                os.system("cp "+optim+'/index.php '+outputfolder)
                c1.SaveAs(outputfolder+histnm+'.png')

print "\nFinished drawing histograms.\nCopying to "+optim
os.system("mkdir -p "+optim )
#os.system('cp -r '+hfolder+' '+optim)
