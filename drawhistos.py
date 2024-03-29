from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TLegend
from ROOT import gROOT,gDirectory, TTree
from collections import OrderedDict
import os, optparse
wloc    = os.environ['WWW']
optim   = wloc+"/susy/optimisation"
hfolder = 'Histograms/'
hfilenm = "Output/allcuts.root"
hfile   = TFile(hfilenm,"READ","Example");
os.system("mkdir -p "+hfolder)

#set x and y ranges
def setRanges(hT2tt,httbar):
    xMin = -999
    xMax =  999
    yMax =  999
    yMin =  0
    #    histo=
    yMax = 1.1*max(hT2tt.GetMaximum(),httbar.GetMaximum())
    xMin = min(hT2tt.GetXaxis().GetXmin(),httbar.GetXaxis().GetXmin())
    xMax = max(hT2tt.GetXaxis().GetXmax(),httbar.GetXaxis().GetXmax())
    yMin = 0.9*min(hT2tt.GetMinimum(),httbar.GetMinimum())
    hT2tt.GetYaxis().SetRangeUser(yMin,0.4)#yMax)
    #hT2tt.GetXaxis().SetRangeUser(xMin,xMax)
    #print xMax, xMin, "-----------------------------"
#Scale histograms
def ScaleToMax(histos,varControl):
    for histo in histos:
        if varControl is True: continue
        norm   = 1
        maxBin = histo.GetMaximumBin()
        norm   = histo.GetBinContent(maxBin)
        histo.Scale(1/norm)

#Scale histograms
def ScaleToInt(histos,varControl):
    for sample in histos:
        histo = gDirectory.Get(histos[sample])
        norm  = histo.GetSumOfWeights()
        if(norm<0.0001): norm = 1.0
        if varControl is False: histo.Scale(1/norm)
        if norm<10: histo.Scale(0.1)
        #print "NORM", histo.GetName(), histo.GetEntries()
#check if histogram exists
def OLDhistoexists(histo, hnm):
    cont       = False
    histo_type = str(type(histo)) #'ROOT.TH1F'>"
    if (bool('ROOT.TH1F' not in histo_type)):
        print "-----------------------------------"
        print "|  "+hnm+" not a histogram\t  |"
        print "| perhaps variable doesn't exist? |"
        print "-----------------------------------"
        cont = True
    return cont
#check if histogram exists
def histoexists(histo, hnm):
    cont       = False
    histo_type = str(type(histo)) #'ROOT.TH1F'>"
    if (bool('ROOT.TH1F' not in histo_type)):
        print "-----------------------------------"
        print "|  "+hnm+" not a histogram\t  |"
        print "| perhaps variable doesn't exist? |"
        print "-----------------------------------"
        cont = True
    return cont

#Scale x axis
def findXrange(hdict, binWOri,nBinFin,xMinOri):
    xMinBin = 999999
    xMaxBin = 0
    for samp in hdict:
        histo   = gDirectory.Get(hdict[samp])
        xMinBin = min(xMinBin,histo.FindFirstBinAbove(0.1))
        xMaxBin = max(xMaxBin,histo.FindLastBinAbove(0.1))
    if(xMinBin>xMaxBin):
        xMinBin = 0
        xMaxBin = nBinFin
    nBin        = xMaxBin-xMinBin
    xLen        = nBin*binWOri
    widthNew    = xLen/nBinFin
    scaleFactor = int(widthNew/binWOri)
    while(scaleFactor>1 and nBinOri%scaleFactor>0): scaleFactor-=1
    if(scaleFactor<1): scaleFactor=1
    xMax=xMinOri+xMaxBin*binWOri
    xMin=xMinOri+xMinBin*binWOri   
    if(xMax>200 and xMin>0): xMin=0
    if(xMin==xMax):
        xMin=-0
        xMax=2*xMax
    return xMin, xMax, scaleFactor
#Scale y axis
def findYrange(hdict, scaleFactor):
    yMax=0
    for samp in hdict:
        histo   = gDirectory.Get(hdict[samp])                
        histo.Rebin(scaleFactor)
        ymaxBin = histo.GetMaximumBin()
        yMax    = max(yMax,histo.GetBinContent(ymaxBin))
    return yMax
            

#Optional parameters
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('--cp'     , dest='copy'   , help='copy to webpage'   , default=False, action='store_true')
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
allvars=["evid", "isSF", "btagW", "bvetoW", "nSV", "nLepton", "ISRcut","nbCleanjets","nbjets" ,\
         "mll", "mt2ll", "detall","dRll", "MET_sumEt","MET_pt", "ptmiss","lep1_pt","lep2_pt" ,\
         "jet1_pt","jet2_pt","dphill", "dphijj","detajj","dRjj","Dnjetstot","Dnbjetstot",\
         "SV_x", "SV_y", "SV_z", "SV_chi2","SV_eta", "SV_phi", "SV_pt", "SV_mass",\
         "PV_x", "PV_y", "PV_z", "PV_chi2","PV_npvs","susyMstop","susyMLSP","ptlepsum",\
         "njets", "nbjets", "nCleanjets", "nbCleanjets","dphill:ptlepsum",\
         "dphil1jmin","dphil1jmax","dphil1bmin","dphil1bmax","dphil1MET","dphilepsumMET",\
         "dphil2jmin","dphil2jmax","dphil2bmin","dphil2bmax","dphil2MET","MET_phi",\
         "bjet1_pt","bjet1_eta", "bjet1_phi", "bjet2_pt", "bjet2_eta", "bjet2_phi"]


controlvars = ["evid","btagW","bvetoW", "susyMstop","susyMLSP","isSF","ISRcut", "nLepton"]
flavours= {'df':'-1','sf':'1'}
bjets   = {'btag': 'btagW','bveto':'bvetoW'}
samples = {'T2tt':'hT2tt', 'ttbar':'httbar', 'WW': 'hWW', 'ZZ':'hZZ','ttZ':'httZ'}
massrange='&& susyMstop>=400 && susyMstop<700 '
deltamass=' && susyMstop-susyMLSP'
mStops = {'mS-400to700_dm-1to125'   :[massrange + deltamass +'<=125', 6],\
          'mS-400to700_dm-125to200' :[massrange + deltamass +'>125'+ deltamass +'<=200', 30],\
          'mS-400to700_dm-200to700' :[massrange + deltamass +'>200', 46],\
        }
testvars=["ptlepsum"]
ptsumcond={"smalldphi":"&& dphill<1.2","largedphi": "&& dphill>1.2", "":""}

bkgs={'ttbar': 4,'WW': 3, 'ZZ' : 8 , 'ttZ': 7  }
ISRcut={"ISR": "ISRcut>0", "Inclusive":"1==1", "lowMT2": "mt2ll<20", "lowMT2_ISR": " mt2ll<20 && ISRcut>0",\
        "highMT2": "mt2ll>100", "highMT2_ISR": " mt2ll>100 && ISRcut>0"}
print "Creating histograms:"
c1 = TCanvas( 'c1', 'Dynamic Filling Example', 200,10, 1600, 900 )
if  opt.test is True:
    variables = testvars
else:
    variables = allvars
for var in variables:
    #Save in different folder if control variable
    cfolder   = ''
    isControl = bool(var in controlvars)        
    if(isControl is True): cfolder = 'Control/'

    if opt.var is not False and opt.var not in var: continue
    btagVeto = False
    if(var in bjets.values()): btagVeto = True 
    print "VARIABLE", var
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

            samplefolder = hfolder+bjet+'/'+isr
            isrnm        = isr
            if isrnm in "Inclusive": isrnm = "All"
            os.system('mkdir -p '+samplefolder)
            for fl in flavours:
                varsplit = var.split(':')
                varnm    = var
                ptcond1  = ""#smalldphi"#"smalldphi"
                histnm   = varnm+'_'+bjet+'_'+fl+'_'+isrnm
                nBinOri  = 500000
                nBinFin  =     12
                xMinOri  =      0.
                xMaxOri  =   5000.
                if "evid" in var:
                    xminOri = 0
                    xMaxOri = 1000000
                    nBinOri = 1000000
                if "SV" in var:
                    xMinOri =     -5
                    xMaxOri =   1000.
                    nBinOri = 100000
                if "deta" in var or "dphi" in var or "dR" in var:
                    xMinOri =    0.
                    xMaxOri =    5.
                    nBinOri = 5000
                if "SV_chi2" in var or "SV_x" in var or "SV_y" in var or "SV_z" in var:
                    xMinOri =  -20.
                    xMaxOri =   20.
                    nBinOri = 10000

                if var in "PV_y":
                    xMinOri =  0.15
                    xMaxOri =  0.20
                    nBinOri = 10000

                if var in "PV_x":
                    xMinOri =  0.09
                    xMaxOri =  0.12
                    nBinOri = 10000
                if(xMaxOri-xMinOri<0.000001): xMaxOri=xMinOri+1
                lenOri  = (xMaxOri-xMinOri)
                binWOri = lenOri/nBinOri
                xRanOri = ""
                is2D    = bool(len(var.split(':'))>1)
                if(is2D is False):
                    xRanOri = "("+str(nBinOri)+","+str(xMinOri)+","+str(xMaxOri)+")"
                elif is2D is True and "ptlepsum"  in var :
                    xRanOri = "(30,40,300,30,0,3.15)"
                dtest=OrderedDict({})
                for sample in samples:
                    #make histo for each variable and condition
                    fillvar   = '>-999'
                    if var in ['jet1_pt', 'jet2_pt', 'mt2ll']: fillvar='>0'
                    tree      = hfile.Get(sample)
                    histo     = samples[sample]+histnm
                    
                    onlyfill = ' && '
                    number   = '-999'
                    if is2D is False:
                        if "dphi" in var: number='-1'
                        onlyfill += var+">"+number+ " && "
                    elif is2D is True:
                        for v in var.split(':'):
                            n='-999'
                            if "dphi" in v: n='-1'
                            onlyfill+=v+">"+n+" && "        

                    treevar   = var+">>"+histo+xRanOri
                    flcond    = "isSF=="+flavours[fl]
                    isrcond   = ISRcut[isr]
                    
                    bjetcond  = bjets[bjet]
                    condition = bjetcond+'*('+flcond+onlyfill+isrcond+ptsumcond[ptcond1]+')'
                    tree.Draw(treevar, condition)
                    dtest[sample]=histo
                    #divide in different mass ranges
                    for mstop in mStops:
                        if sample not in 'T2tt': continue
                        splitcond    = condition[:-1]
                        mStopcond    = splitcond+mStops[mstop][0]+')'
                        treemStop = var+'>>'+histo+'_'+mstop+xRanOri
                        tree.Draw(treemStop,mStopcond)
                        dtest[sample+'_'+mstop]=histo+'_'+mstop
            
                #Set proper ranges and rebin histogram
                if(is2D is False):
                    xMin, xMax, scaleFactor = findXrange(dtest, binWOri, nBinFin, xMinOri)
                    ScaleToInt(dtest, isControl)
                    yMax   = findYrange(dtest, scaleFactor)
                    legend = TLegend(0.7,0.7,0.9,0.9);
                    for idx,sam in enumerate(dtest):
                        isSame = 'same'
                        hist   = gDirectory.Get(dtest[sam])
                        if idx is 0:
                            hist.GetXaxis().SetRangeUser(xMin,xMax)
                            hist.GetYaxis().SetRangeUser(0,1.1*yMax)
                            hist.SetTitle(histnm)
                            hist.SetStats(False)

                            isSame = ''
                            hist.Draw('hist')
                        if sam in "T2tt":
                            hist.Draw('hist same')
                            hist.SetLineColor(2)
                        for ms in mStops:
                            if ms in sam :
                                hist.SetLineColor(mStops[ms][1])
                                hist.Draw('same')
                        isBkg = ''
                        lab   = "le"
                        if sam in bkgs:
                            lab   = "f"
                            isBkg = 'hist '
                            hist.SetLineColor(bkgs[sam])
                        hist.Draw(isBkg+'same')
                        legend.AddEntry(hist,sam,lab);
                    legend.Draw()
                    #save to image
                    outputfolder = samplefolder+'/'+cfolder
                    os.system('mkdir -p '+outputfolder)
                    os.system("cp "+optim+'/index.php '+outputfolder)
                    c1.SaveAs(outputfolder+histnm+'.png')
                elif is2D is True:
            
                    outputfolder = samplefolder+'/'+varnm+'/'+cfolder
                    os.system('mkdir -p '+outputfolder)
                    os.system("cp "+optim+'/index.php '+outputfolder)
                    
                    for sam in dtest:
                        hist = gDirectory.Get(dtest[sam])
                        hist.SetStats(False)
                        hist.SetTitle(histnm+'_'+sam)
                        #xMax=hist.FindLastBinAbove(1)
                        #hist.GetXaxis().SetRangeUser(0,xMax)
                        hist.Draw('colz')
                        c1.SaveAs(outputfolder+histnm+sam+'.png')
               
                
print "\nFinished drawing histograms.\nCopying to "+optim
os.system("mkdir -p "+optim )
if opt.copy is True: os.system('cp -r '+hfolder+' '+optim)
