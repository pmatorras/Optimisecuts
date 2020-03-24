from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TLegend
from ROOT import gROOT,gDirectory, TTree
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
    yMax=1.1*max(hT2tt.GetMaximum(),httbar.GetMaximum())
    xMin=min(hT2tt.GetXaxis().GetXmin(),httbar.GetXaxis().GetXmin())
    xMax= max(hT2tt.GetXaxis().GetXmax(),httbar.GetXaxis().GetXmax())
    yMin=0.9*min(hT2tt.GetMinimum(),httbar.GetMinimum())
    hT2tt.GetYaxis().SetRangeUser(yMin,yMax)
    hT2tt.GetXaxis().SetRangeUser(xMin,xMax)
    print xMax, xMin, "-----------------------------"
#Scale histograms
def ScaleToMax(histos,varControl):
    for histo in histos:
        if varControl is True: continue
        norm=1
        maxBin=histo.GetMaximumBin()
        norm=histo.GetBinContent(maxBin)
        #if maxE>norm: norm=maxE
        #print "SCALING:\n",histo.GetName(), norm
        histo.Scale(1/norm)

#Scale histograms
def Scalehisto(histo,varControl):
    norm=histo.GetSumOfWeights()
    if(norm<0.0001): norm=1.0
    if varControl is False: histo.Scale(1/norm)

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
            

#Optional parameters
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('--test' , dest='test' , help='run only few plits', default=False, action='store_true')
parser.add_option('--rmold' , dest='rmold' , help='remove old images', default=False, action='store_true')
parser.add_option('--nobatch' , dest='nobatch' , help='remove old images', default=False, action='store_true')
(opt, args) = parser.parse_args()

#Apply options
if opt.rmold is True:
    print "removing histograms in", hfolder
    os.system('rm -r '+hfolder)
gROOT.SetBatch(not(opt.nobatch))
c1 = TCanvas( 'c1', 'Dynamic Filling Example', 200,\
 10, 1600, 900 )

#Define variables to over which make histograms
variables=["evid_",  "nSV_","Dnjetstot_", "isSF_", \
           "btagW_", "bvetoW_","MET_sumEt_","MET_pt_", \
           "PV_x_", "PV_y_", "PV_z_", "PV_npvs_", "PV_chi2_",\
           "SV_eta_", "SV_phi_", "SV_pt_", "SV_mass_", \
           "SV_x_", "SV_y_", "SV_z_", "SV_chi2_","mll_", \
           "mt2ll_", "nLepton_", "ptmiss_", "susyMstop_","susyMLSP__",\
           "lep1_pt_","lep2_pt_","jet1_pt_","jet2_pt_","dphill_",\
           "detall_","dRll_","dphijj_","detajj_","dRjj_"]
controlvars=["evid_","btagW_","bvetoW_", "susyMstop_","susyMLSP__","isSF_", "nLepton_"]
flavours={'df':'-1','sf':'1'}
bjets={'btag': 'btagW_','bveto':'bvetoW_'}
samples={'T2tt':'hT2tt', 'ttbar':'httbar', 'WW': 'hWW'}
mStops={'mS-400to700_dm-1to200' :'&& susyMstop>=400 && susyMstop<600 && susyMstop-susyMLSP<=200',\
        'mS-400to700_dm-200to700' :'&& susyMstop>=400 && susyMstop<600 && susyMstop-susyMLSP>200',\
        }
        #'mS-700to1200':'&& susyMstop>=700 && susyMstop<=1200'       }

print "Creating histograms:"
for var in variables:
    #Save in different folder if control variable
    cfolder=''
    isControl=bool(var in controlvars)        
    if(isControl is True): cfolder='Control/'

    if opt.test==True and "ptmiss" not in var: continue
    btagVeto=False
    if(var in bjets.values()): btagVeto=True 
    for bjet in bjets:
        
        if btagVeto is True and bjet not in var:
            print "Omitting histogram since var="+var+" and bjet="+bjet
            continue
        os.system('mkdir -p '+hfolder+bjet)
        for fl in flavours:
            histnm = var+bjet+'_'+fl
            nBinOri= 1000000
            nBinFin= 50
            xMinOri=-5000.
            xMaxOri= 5000.
            lenOri= (xMaxOri-xMinOri)
            binWOri=lenOri/nBinOri
            xRanOri= "("+str(nBinOri)+","+str(xMinOri)+","+str(xMaxOri)+")"
            for sample in samples:
                #make histo for each variable and condition
                tree      = hfile.Get(sample)
                histo     = samples[sample]+histnm
                treevar   = var+sample+">>"+histo+xRanOri
                flcond    = "isSF_"+sample+"=="+flavours[fl]
                onlyfill  = var+sample+'>-999'
                bjetcond  = bjets[bjet]+sample
                condition = bjetcond+'*('+flcond+'&&'+onlyfill+')'
                tree.Draw(treevar, condition)

                #divide in different mass ranges
                for mstop in mStops:
                    if sample not in 'T2tt': continue
                    splitcond    = condition.split(')')
                    mStopcond    = splitcond[0]+mStops[mstop]+')'
                    treemStop = var+sample+'>>'+histo+'_'+mstop+xRanOri
                    tree.Draw(treemStop,mStopcond)                
            histosig= samples['T2tt']+histnm
            print "HISTOGRAM:\t", histnm

            #Get  histograms and check its existance
            hT2tt      = gDirectory.Get(samples["T2tt"]+histnm)
            httbar     = gDirectory.Get(samples["ttbar"]+histnm)
            hWW        = gDirectory.Get(samples["WW"]+histnm)
            cont_T2tt  = histoexists(hT2tt,"hT2tt")
            cont_ttbar = histoexists(httbar,"httbar")
            cont_WW    = histoexists(hWW,"WW")

            if(True in [cont_T2tt,cont_ttbar,cont_WW]): continue
            
            #Set proper ranges
            
            #Draw histogram and legend
            hT2tt.SetTitle(histnm)
            hT2tt.SetStats(False)
            hT2tt.SetLineColor(2)
            hWW.SetLineColor(3)
            xBinMin=hT2tt.FindFirstBinAbove(0.0001)
            xBinMax=hT2tt.FindLastBinAbove(0.0001)
            nBin=xBinMax-xBinMin
            xLen=nBin*binWOri
            widthNew=xLen/nBinFin
            scaleFactor=int(widthNew/binWOri)
            while(nBinOri%scaleFactor>0): scaleFactor-=1
            hT2tt.Rebin(scaleFactor)
            httbar.Rebin(scaleFactor)
            hWW.Rebin(scaleFactor)
            xMax=xMinOri+xBinMax*binWOri
            xMin=xMinOri+xBinMin*binWOri
            if(xMin>0): xMin=0
            hT2tt.GetXaxis().SetRangeUser(xMin,xMax)
            Scalehisto(httbar,isControl)
            Scalehisto(hT2tt,isControl)
            Scalehisto(hWW,isControl)
            ScaleToMax([hT2tt,httbar,hWW], isControl)
            #setRanges(hT2tt,httbar)
            
            hT2tt.Draw('hist')
            httbar.Draw("hist same")
            hWW.Draw("hist same")
            legend = TLegend(0.7,0.7,0.9,0.9);
            legend.AddEntry(hT2tt,"T2tt","f");
            legend.AddEntry(httbar,"ttbar","f");
            legend.AddEntry(hWW,'WW','f')
            #draw masspoint ranges
            for idx,mstop in enumerate(mStops):
                hmStop = gDirectory.Get(histosig+'_'+mstop)
                normmStop  = hmStop.GetSumOfWeights();
                hmStop.SetLineColor(idx+6)
                hmStop.Rebin(scaleFactor)
                ScaleToMax([hmStop],isControl)
                #hmStop.Draw('same')
                legend.AddEntry(hmStop,"T2tt-"+mstop,"f");
                
            legend.Draw()

            #save to image
            outputfolder=hfolder+bjet+'/'+cfolder
            os.system('mkdir -p '+outputfolder)
            os.system("cp "+optim+'/index.php '+outputfolder)
            c1.SaveAs(outputfolder+histnm+'.png')

print "\nFinished drawing histograms.\nCopying to "+optim
os.system("mkdir -p "+optim )
os.system('cp -r '+hfolder+' '+optim)
