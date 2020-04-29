import os, sys
from ROOT import TH1D, TH2D, TFile, TTree, TCanvas, gROOT, gStyle
from array import array

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

def addhisto(sig, bkg, isSig, nsig, count=False):
    histo=inpfile.Get(histnm)
    if(isSig is True):
        if count is True: nsig+=1
        sig.Add(histo)
    else: bkg.Add(histo)
    return nsig
            

inpfile=TFile ('../rootfiles/plots_HighPtMissOptimisationRegion_SM-T2tt_mS-400to700.root', "READ")
wloc=os.environ['WWW']
optim=wloc+"/susy/optimisation/"
regions=["VR1_Tag_sf", "VR1_Tag_em", "VR1_Veto_em", "VR1_Veto_sf"]
dmass={"dm_1to125": [1,125], "dm_125to200" : [125,200] , "dm_200to700": [200,700]}
gROOT.SetBatch(True)#False)#True)
folder="../Histograms/significance/"
hminMT2   =    0
hmaxMT2   = 1000
hminPtm   =    0
hmaxPtm   = 2000
nMT2 =  100
nPtm =  100
        
for dm in dmass:
    print "MASS:\t", dm
    foldm=folder+dm
    os.system("mkdir -p "+foldm)
    for reg in regions:
        print "REGION:", reg
        dmmin= dmass[dm][0]
        dmmax= dmass[dm][1]
        sigunrol  = TH1D("sigunrol" +reg+dm, reg, 10000  , 0, 10000)
        bkgunrol  = TH1D("bkgunrol" +reg+dm, reg, 10000  , 0, 10000)
        sigPtm = TH1D("sigPtmiss"+reg+dm, reg, nPtm, 0, hmaxPtm)
        bkgPtm = TH1D("bkgPtmiss"+reg+dm, reg, nPtm, 0, hmaxPtm)
        sigMT2 = TH1D("sigMT2"   +reg+dm, reg, nMT2   , 0, hmaxMT2)
        bkgMT2 = TH1D("bkgMT2"   +reg+dm, reg, nMT2   , 0, hmaxMT2)

        signal2D = TH2D("signal2D"+reg+dm, reg, nPtm, hminPtm, hmaxPtm,nMT2, hminMT2,hmaxMT2)
        bkg2D    = TH2D("bkgs2D"+reg+dm  , reg, nPtm, hminPtm, hmaxPtm,nMT2, hminMT2,hmaxMT2)
        signif2D = TH2D("signif2D"+reg+dm, reg, nPtm, hminPtm, hmaxPtm,nMT2, hminMT2,hmaxMT2)
        mSmin=400
        mSmax=700
        signif2D.SetTitle("m_{T2ll} vs p_{T}^{miss} ("+reg+'-'+dm+")")
        signif2D.SetYTitle("m_{T2ll} [GeV]")
        signif2D.SetXTitle("p_{T}^{miss} [GeV]")
        allhistos = getall(inpfile)# folder)
        nsig      = 1

        #For a given region and dm, add all histograms
        for ihis, histloc in enumerate(allhistos):
            histnm = histloc[0]
            hsplit = histnm.split('mS')
            isSig  = False
            if(reg not in histnm): continue
            #if ihis>3: break # continue

            if len(hsplit) >1:
                isSig = True
                mass  = hsplit[1].replace('-','_').split('_')
                mS    = int(mass[1])
                mX    = int(mass[3])
                dm_i  = mS-mX
                if(mS<mSmin    or mS>mSmax  ): continue
                if(dm_i<=dmmin or dm_i>dmmax): continue
    
            if 'ptmissmt2' in histnm: nsig=addhisto(sigunrol,bkgunrol,isSig, nsig,True)
            elif 'mt2'     in histnm: nsig=addhisto(sigMT2,bkgMT2,isSig, nsig)
            elif 'ptmiss'  in histnm: nsig=addhisto(sigPtm,bkgPtm,isSig, nsig)
        print "nsig", nsig
        #if(ibkg>0): print 10*x,20*y, ibkg, isig


        
        #recover the 2D histogram from the unrolled 
        nbin = sigunrol.GetNbinsX()
        x = 0
        y = 0
        n = 0
        
        for i in range(1,nbin):
            if (i % nPtm is 1):
                y+=1
                x=0
            x+=1
            isig=sigunrol.GetBinContent(i)/nsig
            ibkg=bkgunrol.GetBinContent(i)
            #if ientry<0.01: continue
            if(isig>0):signal2D.SetBinContent(x,y,isig)
            if(ibkg>=0):
                bkg2D.SetBinContent(x,y, ibkg)
                #if(x<10 and y<5):print x, y, ibkg
            if(ibkg+isig>0 and isig>=0): signif2D.SetBinContent(x,y,isig/np.sqrt(ibkg+isig))



        
        #Make profiles (currently unused)
        bkgprofX=bkg2D.ProfileX()
        bkgprofY=bkg2D.ProfileY()
        sigprofY=signal2D.ProfileY()

        #Get desired variable binning
        wsumMT2=0
        wYmax=500
        wsumPtm=0
        binsMT2=[0]
        rangMT2=[0]
        binsPtm=[5]
        rangPtm=[100]
        totWPtm=bkgPtm.GetSumOfWeights()
        totWMT2=bkgMT2.GetSumOfWeights()
        wmaxPtm=0.05*totWPtm
        wmaxMT2=0.005*totWMT2#0.01
        
        #print wYmax/totWY
        #exit()
        print "tot2\t", int(totWPtm), '\t',int(totWMT2)
        for i in range(1,99):
            wbkgMT2=bkgMT2.GetBinContent(i)
            wbkgPtm=bkgPtm.GetBinContent(i)
            #print i,"\tX:", int(wbkgMT2), "\tY:", int(wbkgPtm)
            if(wbkgPtm>0): wsumPtm+=wbkgPtm
            if(wbkgMT2>0): wsumMT2+=wbkgMT2
            if wsumPtm>wmaxPtm:
                wsumPtm=0
                locPtm=i*(hmaxPtm-hminPtm)/nPtm
                binsPtm.append(i)
                rangPtm.append(locPtm)
                print "Y:",locPtm, "\tnew bin", i
            if wsumMT2>wmaxMT2:
                wsumMT2=0
                locMT2=i*(hmaxMT2-hminMT2)/nMT2
                binsMT2.append(i)
                rangMT2.append(locMT2)
                print "X:",locMT2, "\tnew bin"
        #print  binsPtm, binsMT2
        print "Ptm:\t",rangPtm, "\nMT2:\t",rangMT2

        #Fill signal and bakground for variable binwidth

        nvarMT2=len(binsMT2)
        nvarPtm=len(binsPtm)
        vartitle="s/#sqrt{s+b} ("+reg+"-"+dm+"); p_{T}^{miss}; m_{T2}^{ll}"
        bkgvar = TH2D("bkgvar"+reg+dm,"",len(rangPtm)-1,array('d',rangPtm),len(rangMT2)-1,array('d', rangMT2))
        sigvar = TH2D("sigvar"+reg+dm,"",len(rangPtm)-1,array('d',rangPtm),len(rangMT2)-1,array('d', rangMT2))
        signifvar = TH2D("signifvar"+reg+dm,vartitle,len(rangPtm)-1,array('d',rangPtm),len(rangMT2)-1,array('d', rangMT2))
        for ibin in range(1, nPtm):
            for jbin in range(1, nMT2):
                bkg=bkg2D.GetBinContent(ibin, jbin)
                sig=signal2D.GetBinContent(ibin, jbin)
                locMT2=(jbin-0.5)*(hmaxMT2-hminMT2)/nMT2
                locPtm=(ibin-0.5)*(hmaxPtm-hminPtm)/nPtm
                bkgvar.Fill(locPtm,locMT2,bkg)
                sigvar.Fill(locPtm,locMT2,sig)
                if locPtm<100 or locPtm>200 or locMT2>100: continue

                
        #Fill significance plot
        for xbin in range (1,nvarPtm+1): #One is due to empty underflow bin
            for ybin in range(1,nvarMT2+1):
                isig = sigvar.GetBinContent(xbin,ybin)
                ibkg = bkgvar.GetBinContent(xbin,1)
                if(isig+ibkg>0): signif = isig/np.sqrt(ibkg+isig)
                else: signif=0
                #print rangPtm[xbin-1], rangMT2[ybin-1],isig#, signif
                signifvar.SetBinContent(xbin, ybin,signif)


        #exit()
        c1 = TCanvas( 'c1', 'Dynamic Filling Example', 200,10, 1600, 900 )
        #bkg2D.Draw('colz')
        #c1.SaveAs('ww.png')

        gStyle.SetPaintTextFormat("4.3f");
        gStyle.SetOptStat(0);

        signifvar.Draw('colz text')
        signifvar.GetYaxis().SetRangeUser(0,400);
        signifvar.GetXaxis().SetRange(0,800);

        #signal2D.Draw('colz')
        #bkgprofY.Draw()#('colz')
        #bkg2D.Draw('colz')
        #signif2D.Draw('colz')
        #signif2D.GetYaxis().SetRangeUser(0,400);
        #signif2D.GetXaxis().SetRangeUser(0,800);
        c1.SaveAs(foldm+"/PtmissvsMT2"+reg+'-'+dm+'.png')
        #os.system("display "+foldm+"/mt2vsPtmiss"+reg+'-'+dm+".png")
        #exit()
        os.system("cp "+optim+'/index.php '+foldm)
cpweb= 'cp -r '+folder+" "+ optim
#os.system(cpweb)
print cpweb    




'''
        #old way to loop over the entries to fill the asymmetric histo
        for ibinf,i in enumerate(range(0,len(binsX))):
            #if(i<binsX[0]):continue
            for jbinf,j in enumerate(range(0,len(binsY))):
               #if(j<binsY[0]):continue
               if jbinf<1 or ibinf<1: continue
               if j>4 or i>4: continue
               for in_i in range(binsX[i-1], binsX[i]):
                   for in_j in range(binsY[j-1],binsY[j]):
                       print i,"in x", in_i,"\t", j,"in y", in_j
                #print ibinf,binsX[i],jbinf, binsY[j]
               #, bkg2D.GetBinContent(i,j)

'''


