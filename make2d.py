import os, sys
from ROOT import TH1D, TH2D, TFile, TTree, TCanvas, gROOT, gStyle
from array import array
import numpy as np
import time
wloc = os.environ['WWW']
start_time= time.time()
gROOT.SetBatch(True)
gStyle.SetPaintTextFormat("4.3f");
gStyle.SetOptStat(0);

#Fill signal and bakground for variable binwidth
def fillvarbins(bkg2D,sig2D,bkgvar,sigvar):
    for ibin in range(1, nPtm):
        for jbin in range(1, nMT2):
            bkg = bkg2D.GetBinContent(ibin, jbin)
            sig = sig2D.GetBinContent(ibin, jbin)
            locMT2 = (jbin-0.5)*(hmaxMT2-hminMT2)/nMT2
            locPtm = (ibin-0.5)*(hmaxPtm-hminPtm)/nPtm
            bkgvar.Fill(locPtm, locMT2, bkg)
            sigvar.Fill(locPtm, locMT2, sig)
            #if sig>0: print locPtm, locMT2, sig


#Fill significance plot for variable binning
def fillsignif(bkgvar,sigvar,signifvar):
    for xbin in range (1,nvarPtm+1): #One is due to empty underflow bin
        for ybin in range(1,nvarMT2+1):
            isig = sigvar.GetBinContent(xbin,ybin)
            ibkg = bkgvar.GetBinContent(xbin,ybin)
            if(isig+ibkg>0): signif = isig/np.sqrt(ibkg+isig)
            else: signif = 0
            signifvar.SetBinContent(xbin, ybin, signif)



            
#Get all histograms/trees within a TFile
def getall(d, basepath="/"):
    for key in d.GetListOfKeys():
        kname = key.GetName()
        if key.IsFolder():
            for i in getall(d.Get(kname), basepath+kname+"/"):
                yield i
        else:
            yield basepath+kname, d.Get(kname)

#Add one histogram
def addhisto(sig, bkg, isSig,histnm, nsig, count=False):
    histo=inpfile.Get(histnm)
    #print histnm, isSig
    if(isSig is True):
        if count is True: nsig+=1
        sig.Add(histo)
    else: bkg.Add(histo)
    return nsig

#Add all histograms given input file
def addhistos(inpfile, dmass,dm, sigunrol,bkgunrol, sigMT2,bkgMT2,sigPtm,bkgPtm):
    allhistos = getall(inpfile)
    nsig      = 0
    print "DM", dm
    #For a given region and dm, add all histograms
    for ihis, histloc in enumerate(allhistos):
        histnm = histloc[0]
        hsplit = histnm.split('mS')
        isSig  = False
        if(reg not in histnm and "comb" not in reg): continue
        if len(hsplit) >1:
            isSig = True
            if "AN" in dm:
                if dmass[dm] not in histnm: continue
            else:
                mass  = hsplit[1].replace('-','_').split('_')
                mS    = int(mass[1])
                mX    = int(mass[3])
                dm_i  = mS-mX
                if(mS<mSmin    or mS>mSmax  ): continue
                if(dm_i<=dmmin or dm_i>dmmax): continue
        if 'ptmissmt2' in histnm:
            #print "nsig", nsig, histnm
            nsig = addhisto(sigunrol,bkgunrol,isSig,histnm, nsig,True)
        elif 'mt2'     in histnm: nsig = addhisto(sigMT2,bkgMT2,isSig,histnm, nsig)
        elif 'ptmiss'  in histnm: nsig = addhisto(sigPtm,bkgPtm,isSig, histnm, nsig)
    print "number of Mass Points:", nsig
    return nsig

#recover the 2D histogram from the unrolled 
def make2D(sigunrol,bkgunrol, sig2D, bkg2D):
    nbin = sigunrol.GetNbinsX()
    x = 0
    y = 0
    n = 0
    
    for i in range(1,nbin):
        if (i % nPtm is 1):
            y+=1
            x =0
        x+=1
        isig = sigunrol.GetBinContent(i)/nsig
        ibkg = bkgunrol.GetBinContent(i)
        if(isig>0): sig2D.SetBinContent(x,y,isig)
        if(ibkg>=0):
            bkg2D.SetBinContent(x,y, ibkg)
            #if(x<10 and y<5):print x, y, ibkg
        if(ibkg+isig>0 and isig>=0): signif2D.SetBinContent(x,y,isig/np.sqrt(ibkg+isig))

#Main area
yearfol='2016'
yearfol='2016-2017-2018'
if '-' in yearfol:
    year=''
else:
    year=yearfol+'_'
yearfol=yearfol+'/'
inpfile = TFile ('../rootfiles/plots_HighPtMissOptimisationRegion_'+year+'SM-T2tt_mS-400to700.root', "READ")
optim   = wloc+"/susy/optimisation/"
folder  ="../Histograms/""significance/"+yearfol
hminMT2 =    0
hmaxMT2 = 1000
hminPtm =    0
hmaxPtm = 2000
nMT2 = 100
nPtm = 100
wPtm = (hmaxPtm-hminPtm)/nPtm
regions = ["combined", "VR1_Tag_sf", "VR1_Tag_em", "VR1_Veto_em", "VR1_Veto_sf"]
dmass   = {"all": [1,700],"ANMP":"mS-450_mX-325","dm_1to125": [1,125], "dm_125to200" : [125,200] , "dm_200to700": [200,700]}

for dm in dmass:
    print "MASS:\t", dm
    for reg in regions:
        dmfol = dm
        if "AN" in dm:
            dmfol = dmass[dm]
            print "Get AN masspoints", dmass[dm]
        else
            dmmin = dmass[dm][0]
            dmmax = dmass[dm][1]
        
        foldm = folder + dmfol + '/' + reg
        os.system("mkdir -p "+foldm)

        #Define necessary histograms
        print "REGION:", reg
        sigunrol = TH1D("sigunrol"+reg+dm, reg, 10000, 0, 10000)
        bkgunrol = TH1D("bkgunrol"+reg+dm, reg, 10000, 0, 10000)
        sigPtm = TH1D("sigPtmiss" +reg+dm, reg, nPtm, 0, hmaxPtm)
        bkgPtm = TH1D("bkgPtmiss" +reg+dm, reg, nPtm, 0, hmaxPtm)
        sigMT2 = TH1D("sigMT2"    +reg+dm, reg, nMT2, 0, hmaxMT2)
        bkgMT2 = TH1D("bkgMT2"    +reg+dm, reg, nMT2, 0, hmaxMT2)

        sig2D    = TH2D("sig2D"   +reg+dm, reg, nPtm, hminPtm, hmaxPtm,nMT2, hminMT2,hmaxMT2)
        bkg2D    = TH2D("bkgs2D"  +reg+dm, reg, nPtm, hminPtm, hmaxPtm,nMT2, hminMT2,hmaxMT2)
        signif2D = TH2D("signif2D"+reg+dm, reg, nPtm, hminPtm, hmaxPtm,nMT2, hminMT2,hmaxMT2)
        mSmin = 400
        mSmax = 700
        signif2D.SetTitle( "m_{T2ll} vs p_{T}^{miss} ("+reg+'-'+dm+")")
        signif2D.SetYTitle("m_{T2ll} [GeV]")
        signif2D.SetXTitle("p_{T}^{miss} [GeV]")

        

        nsig = addhistos(inpfile, dmass,dm, sigunrol,bkgunrol, sigMT2,bkgMT2,sigPtm,bkgPtm)
        make2D(sigunrol,bkgunrol, sig2D, bkg2D)


        #Make profiles (currently unused)
        bkgprofX = bkg2D.ProfileX()
        bkgprofY = bkg2D.ProfileY()
        sigprofY = sig2D.ProfileY()

        #Variable binning if want to reproduce AN binning
        binAN  = True
        varbin = 'varbin'
        if(binAN is True):
            varbin  = 'ANbin'
            rangMT2 = [   0,  20,  40,  60,  80, 100, 120]
            rangPtm = [ 100, 140, 200, 300]

        #Otherwise pick binning from 1D distrib weights
        else:
            totWPtm = bkgPtm.GetSumOfWeights()
            totWMT2 = bkgMT2.GetSumOfWeights()
            wsumMT2 =   0
            wYmax   = 500
            wsumPtm =   0
            binsMT2 = [0]
            rangMT2 = [0]
            binsPtm = [5]
            rangPtm = [100]
            wmaxPtm = 0.05*totWPtm
            wmaxMT2 = 0.005*totWMT2#0.01
        
            print "tot2\t", int(totWPtm), '\t',int(totWMT2)
            for i in range(1,99):
                wbkgMT2 = bkgMT2.GetBinContent(i)
                wbkgPtm = bkgPtm.GetBinContent(i)

                if(wbkgPtm>0): wsumPtm+=wbkgPtm
                if(wbkgMT2>0): wsumMT2+=wbkgMT2
                if wsumPtm>wmaxPtm:
                    wsumPtm = 0
                    locPtm  = i*(hmaxPtm-hminPtm)/nPtm
                    binsPtm.append(i)
                    rangPtm.append(locPtm)
                    print "Y:",locPtm, "\tnew bin", i
                if wsumMT2>wmaxMT2:
                    wsumMT2 = 0
                    locMT2  = i*(hmaxMT2-hminMT2)/nMT2
                    binsMT2.append(i)
                    rangMT2.append(locMT2)
                    print "X:",locMT2, "\tnew bin"
                    
        print "Ptm:\t",rangPtm, "\nMT2:\t",rangMT2
        nvarMT2=len(rangMT2)
        nvarPtm=len(rangPtm)


        #Define 2D histograms with the variable binwidth
        vartitle    = "("+reg+"-"+dm+"); p_{T}^{miss}; m_{T2}^{ll}"
        signiftitle = "s/#sqrt{s+b} "+vartitle
        signifsqtit = "significance squared"+vartitle
        nvarbinx  = len(rangPtm)-1
        nvarbiny  = len(rangMT2)-1
        varbinx   = array('d',rangPtm)
        varbiny   = array('d',rangMT2)
        bkgvar    = TH2D("bkgvar"+reg+dm,"bkg "   +vartitle, nvarbinx, varbinx, nvarbiny, varbiny)
        sigvar    = TH2D("sigvar"+reg+dm,"signal "+vartitle, nvarbinx, varbinx, nvarbiny, varbiny)
        signifvar = TH2D("signifvar"+reg+dm, signiftitle, nvarbinx, varbinx, nvarbiny, varbiny)
        
        fillvarbins(bkg2D,sig2D,bkgvar,sigvar)
        
        fillsignif(bkgvar,sigvar,signifvar)
        
        statt=array('d',7*[0.])
        e= signifvar.Integral()
        print e, signifvar.Sumw2()
        signifvarsq = TH2D("signifvarsq"+reg+dm, signiftitle, nvarbinx, varbinx, nvarbiny, varbiny)
        signifvarsq.Multiply(signifvar,signifvar)#,signifvar)
        print signifvarsq.Integral()



        
        c1 = TCanvas( 'c1', 'Dynamic Filling Example', 200,10, 1200, 900 )


        #exit()
        
        sigvar.Draw('colz text')
        sigvar.GetYaxis().SetRange(0,400);
        sigvar.GetXaxis().SetRange(0,800);
        plot_nm=varbin+"-PtmissvsMT2"#+reg+'-'+dm+'.png'
        c1.SaveAs(foldm+"/"+plot_nm+'_signal.png')


        bkgvar.Draw('colz text')
        bkgvar.GetYaxis().SetRange(0,400);
        bkgvar.GetXaxis().SetRange(0,800);
        c1.SaveAs(foldm+"/"+plot_nm+'_bkg.png')


        signifvar.Draw('colz text')
        signifvar.GetYaxis().SetRange(0,400);
        signifvar.GetXaxis().SetRange(0,800);
        c1.SaveAs(foldm+'/'+plot_nm+'_signif.png')


        signifvarsq.Draw('colz text')
        signifvarsq.GetYaxis().SetRange(0,400);
        signifvarsq.GetXaxis().SetRange(0,800);
        c1.SaveAs(foldm+'/'+plot_nm+'_signifsq.png')

        #exit()
        os.system("cp "+wloc+'/index.php '+foldm)
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


