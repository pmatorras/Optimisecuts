from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TLegend
from ROOT import gROOT,gDirectory, TTree
from definebranches import hfilenm, defVectors
import os

hfolder='../Histograms/'
hfilenm="../Output/allpostcuts.root"
hfile  = TFile(hfilenm,"READ","Example");
os.system("mkdir -p "+hfolder)


gROOT.SetBatch(1)
c1 = TCanvas( 'c1', 'Dynamic Filling Example', 200,\
 10, 1600, 900 )
doLogy=False

variables=["evid_",  "nSV_","Dnjetstot_", "isSF_", \
           "btagW_", " bvetoW_","MET_sumEt_","MET_pt_", \
           "PV_x_", "PV_y_", "PV_z_", "PV_npvs_", "PV_chi2_",\
           "SV_eta_", "SV_phi_", "SV_pt_", "SV_mass_", \
           "SV_x_", "SV_y_", "SV_z_", "SV_chi2_"," mll_", \
           "mt2ll_", "nLepton_", "ptmiss_", "susyMstop_", "susyMLSP__"]
  
flavours={'df':'-1','sf':'1'}
bjets={'btag': 'btagW_','bveto':'bvetoW_'}
samples={'T2tt':'hT2tt', 'ttbar':'httbar'}

for var in variables:
    for bjet in bjets:
        os.system('mkdir -p '+hfolder+bjet)
        for fl in flavours:
            normT2tt=1
            normttbar=1
            nttbar=1
            for sample in samples:
                #print "sample",sample
                tree=hfile.Get(sample)
                hterm='_'+var+'_'+bjet+'_'+fl
                histo=var+sample+">>"+samples[sample]+hterm
                flcond="isSF_"+sample+"=="+flavours[fl]
                bjetcond=bjets[bjet]+sample
                condition=bjetcond+'*('+flcond+')'
                #print "condition", condition, "histo", histo
                tree.Draw(histo, condition)

            histnm = var+bjet+'_'+fl
            print "HISTOGRAM:\t", histnm

            #Get and normalise histograms
            hT2tt = gDirectory.Get(samples["T2tt"]+hterm);
            hT2tt.SetTitle(histnm)
            hT2tt.SetStats(False)
            httbar = gDirectory.Get(samples["ttbar"]+hterm)
            normT2tt  = hT2tt.GetSumOfWeights();
            normttbar = httbar.GetSumOfWeights();
            nttbar=httbar.GetEntries()
            if normttbar<0.01: normttbar=1.0
            hT2tt.Scale(1/normT2tt);
            httbar.Scale(1/normttbar);
            
            #Make legend
            legend = TLegend(0.7,0.7,0.9,0.9);
            legend.AddEntry(hT2tt,"T2tt","f");
            legend.AddEntry(httbar,"ttbar","f");

            #set x and y ranges
            xMin=-999
            xMax=-999
            yMax= 999
            yMin= 0

            yMax=max(hT2tt.GetMaximum(),httbar.GetMaximum())
            xMin=min(hT2tt.GetXaxis().GetXmin(),httbar.GetXaxis().GetXmin())
            xMax= max(hT2tt.GetXaxis().GetXmax(),httbar.GetXaxis().GetXmax())
            yMin=min(hT2tt.GetMinimum(),httbar.GetMinimum())
            hT2tt.GetYaxis().SetRangeUser(0.9*yMin,1.1*yMax)
            hT2tt.GetXaxis().SetRangeUser(xMin,xMax)
            hT2tt.SetLineColor(2)

            
            hT2tt.Draw('hist')
            legend.Draw()
            httbar.Draw("hist same")
            c1.SaveAs(hfolder+bjet+'/'+histnm+'.png')


optim="/home/pablinux/cernbox/www/susy/optimisation"
os.system("mkdir -p "+optim )
os.system('cp -r '+hfolder+' '+optim)
