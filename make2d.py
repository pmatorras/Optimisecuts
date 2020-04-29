import os, sys
from ROOT import TH1D, TH2D, TFile, TTree, TCanvas, gROOT, gStyle
import numpy as np


def getall(d, basepath="/"):
    "Generator function to recurse into a ROOT file/dir and yield (path, obj) pairs"
    for key in d.GetListOfKeys():
        kname = key.GetName()
        if key.IsFolder():
            # TODO: -> "yield from" in Py3
            for i in getall(d.Get(kname), basepath+kname+"/"):
                yield i
        else:
            yield basepath+kname, d.Get(kname)



inpfile=TFile ('../lxplus/plots_HighPtMissOptimisationRegion_SM-T2tt_mS-400to700.root', "READ")
wloc=os.environ['WWW']
optim=wloc+"/susy/optimisation/"
regions=["VR1_Tag_sf", "VR1_Tag_em", "VR1_Veto_em", "VR1_Veto_sf"]
dmass={"dm_1to125": [1,125], "dm_125to200" : [125,200] , "dm_200to700": [200,700]}
gROOT.SetBatch(True)
folder="../Histograms/significance/"

for dm in dmass:
    foldm=folder+dm
    os.system("mkdir -p "+foldm)
    for reg in regions:
        dmmin= dmass[dm][0]
        dmmax= dmass[dm][1]
        nMT2=100
        nptmiss=100
        signal1D= TH1D("signal1D"+reg+dm, reg, 10000, 0, 10000)
        bkg1D   = TH1D("bkgs1D"+reg+dm  , reg, 10000, 0, 10000)
        signal2D= TH2D("signal2D"+reg+dm, reg, nMT2, 0, 1000, nptmiss, 0, 2000)
        bkg2D   = TH2D("bkgs2D"+reg+dm  , reg, nMT2, 0, 1000, nptmiss, 0, 2000)
        signif2D= TH2D("signif2D"+reg+dm, reg, nMT2, 0, 1000, nptmiss, 0, 2000)
        mSmin=400
        mSmax=700
        signif2D.SetTitle("m_{T2ll} vs p_{T}^{miss} ("+reg+'-'+dm+")")
        signif2D.SetYTitle("m_{T2ll} [GeV]")
        signif2D.SetXTitle("p_{T}^{miss} [GeV]")
        allhistos=getall(inpfile)# folder)
        nsig=1
        for ihis, histloc in enumerate(allhistos):
            histnm=histloc[0]
            hsplit = histnm.split('mS')
            isSig =False
            if(reg not in histnm): continue
            #if ihis>3: break # continue

            if len(hsplit) >1:
                isSig=True
                mass=hsplit[1].replace('-','_').split('_')
                mS=int(mass[1])
                mX=int(mass[3])
                dm_i=mS-mX
                if(mS<mSmin or mS>mSmax):
                    #print "mS not in ", mSmin, mSmax
                    continue
                
                if(dm_i<=dmmin or dm_i>dmmax):
                    #print "histo not in mass range"
                    continue
    
            histo=inpfile.Get(histnm)
            if(isSig is True):
                nsig+=1
                signal1D.Add(histo)
            else: bkg1D.Add(histo)
    

        print nsig
        nbin=signal1D.GetNbinsX()
        y=0
        x=0
        n=0
        for i in range(1,nbin):
            if (i % nMT2 is 1):
                x+=1
                y=0
            y+=1
            isig=signal1D.GetBinContent(i)/nsig
            ibkg=bkg1D.GetBinContent(i)
            #if ientry<0.01: continue
            if(isig>0):signal2D.SetBinContent(x,y,isig)
            if(ibkg>0): bkg2D.SetBinContent(x,y, ibkg)
            if(ibkg+isig>0 and isig>0): signif2D.SetBinContent(x,y,isig/np.sqrt(ibkg+isig))
            #print nbin, signal1D.GetBinContent(nbin)
            #break
        #print signal1D.GetEntries()

        #exit()
        c1 = TCanvas( 'c1', 'Dynamic Filling Example', 200,10, 1600, 900 )
        #signal2D.Draw('colz')
        #bkg2D.Draw('colz')
        gStyle.SetPaintTextFormat("4.1f");
        signif2D.Draw('colz')
        signif2D.GetYaxis().SetRangeUser(0,400);
        signif2D.GetXaxis().SetRangeUser(0,800);
        c1.SaveAs(foldm+"/mt2vsptmiss"+reg+'-'+dm+'.png')
        #os.system("display "+regions[0]+".png")
        print signal1D.GetEntries()

        #exit()
        os.system("cp "+optim+'/index.php '+foldm)
cpweb= 'cp -r '+folder+" "+ optim
os.system(cpweb)
print cpweb    
