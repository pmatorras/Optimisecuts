import os, sys,optparse
from ROOT import TH1D, TH2D, TFile, TTree, TCanvas, gROOT, gStyle
from array import array
import numpy as np
import time
from getparams import *
print inpfnm, m1min
start_time= time.time()
gROOT.SetBatch(True)
gStyle.SetPaintTextFormat("4.3f");
gStyle.SetOptStat(0);

dorelerr=opt.relerr
relstr='/'
rangMT2 = [   0,  20,  40,  60,  80, 100, 120]
rangPtm = [ 100, 140, 200, 300]    
if opt.AN is True or True not in [opt.Ptm, opt.MT2]:
    varbin='ANbin'
else:
    varbin='varbin'
    if opt.MT2 is True:
        varbin+='_MT2'
        rangMT2 = [   0,  20,  40,  60,  80, 100, 160]
    if opt.Ptm is True:
        varbin+='_Ptm'
        rangPtm = [ 100, 160, 220,280, 380]
    
print opt.AN


#Main area
regions = ["combined", "VR1_Tag_sf", "VR1_Tag_em", "VR1_Veto_em", "VR1_Veto_sf"]
if dorelerr is True:
    relstr='_errors/'
    if 'T2tt' in sig_nm:
        dmass ={"MP-500_413": "mS-500_mX-413","MP-500_400": "mS-500_mX-400","MP-500_375": "mS-500_mX-375","MP-500_375": "mS-500_mX-375","MP-500_350": "mS-500_mX-350","MP-500_325": "mS-500_mX-325", "dm_1to200": [1,200]} #till 325
    elif 'SlepSnu' in sig_nm:
        dmass = {"all" : [1,1250]}
else:
    if 'T2tt' in sig_nm:
        dmass   = {"all": [1,700],"ANMP":"mS-450_mX-325","dm_1to125": [1,125], "dm_125to200" : [125,200] , "dm_200to700": [200,700]}
    elif 'SlepSnu' in sig_nm:
        dmass  = {#"all":[1,1250],
                  "mC800to1200": [1,1250]}

#Fill signal and bakground for variable binwidth
def fillvarbins(bkg2D,sig2D,bkgvar,sigvar, dorelerr=False, h_sigerrsq=None):
    for ibin in range(1, nPtm):
        for jbin in range(1, nMT2):
            bkg = bkg2D.GetBinContent(ibin, jbin)
            sig = sig2D.GetBinContent(ibin, jbin)
            locMT2 = (jbin-0.5)*(hmaxMT2-hminMT2)/nMT2
            locPtm = (ibin-0.5)*(hmaxPtm-hminPtm)/nPtm
            bkgvar.Fill(locPtm, locMT2, bkg)
            sigvar.Fill(locPtm, locMT2, sig)
            if dorelerr is True:
                sigerr= sig2D.GetBinError(ibin,jbin)
                if(sig>0):
                    if sigerr/sig >1:
                        print locMT2, locPtm, "\t",sig,"\t",sigerr
                    h_sigerrsq.Fill(locPtm,locMT2,sigerr*sigerr)
            #if sig>0: print locPtm, locMT2, sig


#Fill significance plot for variable binning
def fillsignif(bkgvar,sigvar,signifvar, nvarPtm, nvarMT2, dorelerr=False,sigerr=None, sigerrsq=None):
    #print "var", nvarPtm, nvarMT2
    for xbin in range (0,nvarPtm+2): #0 is underflow, 2 overflow + range
        #print "xbin", xbin
        for ybin in range(0,nvarMT2+2):
            isig = sigvar.GetBinContent(xbin,ybin)
            ibkg = bkgvar.GetBinContent(xbin,ybin)
            if(isig+ibkg>0):
                signif = isig/np.sqrt(ibkg+isig)
            
                #print xbin, ybin, round(isig,3) #round(signif,3)
                
            else: signif = 0
            signifvar.SetBinContent(xbin, ybin, signif)
            
            if dorelerr is True:
                ierrsq = sigerrsq.GetBinContent(xbin,ybin)
                if ierrsq>=0:
                    ierr = np.sqrt(ierrsq)
                else:
                    print "ierr negative"
                    exit()
                sigerr.SetBinContent(xbin,ybin,ierr)


            
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
    #print "HISTO-->", histnm, isSig
    if(isSig is True):
        if count is True:
            nsig+=1
            #print histnm, nsig
        sig.Add(histo)
    else: bkg.Add(histo)
    return nsig

#Add all histograms given input file
def addhistos(inpfile, dmass,dm, dmmin,dmmax, reg, sigunrol,bkgunrol, sigMT2,bkgMT2,sigPtm,bkgPtm):
    allhistos = getall(inpfile)
    nsig      = 0
    global m1min, m1max
    
    dmsplit=dm.split(splitnm)
    if len(dmsplit)>1:
        m1min=int(dmsplit[1].split("to")[0])
        m1max=int(dmsplit[1].split("to")[1])
        print "MASSPOINTS IN RANGE", splitnm+":\t[",m1min,",",m1max,"]"
    #For a given region and dm, add all histograms
    for ihis, histloc in enumerate(allhistos):
        histnm = histloc[0]
        hsplit = histnm.split(splitnm)
        isSig  = False
        if(reg not in histnm and "comb" not in reg): continue
        if len(hsplit) >1:
            isSig = True
            if "MP" in dm:
                if dmass[dm] not in histnm: continue
            else:
                #print "else", histnm
                mass  = hsplit[1].replace('-','_').split('_')
                m1    = int(mass[1])
                m2    = int(mass[3])
                dm_i  = m1-m2
                if(m1<m1min    or m1>m1max  ): continue
                if(dm_i<=dmmin or dm_i>dmmax):
                    print "omittin", histnm
                    continue
                
        if 'ptmissmt2' in histnm:
            #print "nsig", nsig, histnm
            nsig = addhisto(sigunrol,bkgunrol,isSig,histnm, nsig,True)
        elif 'mt2'     in histnm: nsig = addhisto(sigMT2,bkgMT2,isSig,histnm, nsig)
        elif 'ptmiss'  in histnm: nsig = addhisto(sigPtm,bkgPtm,isSig, histnm, nsig)
    print "number of Mass Points:", nsig
    return nsig

#Define necessary histograms
def define_histos(reg,dm):
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
    signif2D.SetTitle( "m_{T2ll} vs p_{T}^{miss} ("+reg+'-'+dm+")")
    signif2D.SetYTitle("m_{T2ll} [GeV]")
    signif2D.SetXTitle("p_{T}^{miss} [GeV]")

    return sigunrol,bkgunrol, sigPtm, bkgPtm, sigMT2, bkgMT2, sig2D, bkg2D, signif2D
        

#recover the 2D histogram from the unrolled 
def make2D(sigunrol,bkgunrol, sig2D, bkg2D,signif2D,nsig):
    nbin = sigunrol.GetNbinsX()
    x = 0
    y = 0
    n = 0
    
    for i in range(1,nbin):
        if (i % nPtm is 1):
            y+=1
            x =0
        x+=1
        iserr= sigunrol.GetBinError(i)/nsig
        isig = sigunrol.GetBinContent(i)/nsig
        ibkg = bkgunrol.GetBinContent(i)
        if(isig>=0):
            sig2D.SetBinContent(x,y,isig)
            sig2D.SetBinError(x,y,iserr)
        if(ibkg>=0):
            bkg2D.SetBinContent(x,y, ibkg)
            #if(x<10 and y<5):print x, y, ibkg
        if(ibkg+isig>0 and isig>=0): signif2D.SetBinContent(x,y,isig/np.sqrt(ibkg+isig))

#Variable binning if want to reproduce AN binning
def varbinfromprof():
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
    return rangPtm, rangMT2


#Draw histograms
def draw_histos(sigvar,bkgvar,signifvar, signifvarsq, foldm, varbin, drawerr=False, sigerr=None, sigerrsq=None, sigrelerr=None):
    c1 = TCanvas( 'c1', 'Dynamic Filling Example', 200,10, 1200, 900 )
    os.system('mkdir -p '+foldm)
    plot_nm=varbin+"-PtmissvsMT2"#+reg+'-'+dm+'.png'

    sigvar.Draw('colz text')
    sigvar.GetYaxis().SetRange(0,400);
    sigvar.GetXaxis().SetRange(0,800);
    c1.SaveAs(foldm+"/"+plot_nm+'_signal.png')

    signifvar.Draw('colz text')
    signifvar.GetYaxis().SetRange(0,400);
    signifvar.GetXaxis().SetRange(0,800);
    c1.SaveAs(foldm+'/'+plot_nm+'_signif.png')

    if drawerr is False:
        

        bkgvar.Draw('colz text')
        bkgvar.GetYaxis().SetRange(0,400);
        bkgvar.GetXaxis().SetRange(0,800);
        c1.SaveAs(foldm+"/"+plot_nm+'_bkg.png')

        signifvarsq.Draw('colz text')
        signifvarsq.GetYaxis().SetRange(0,400);
        signifvarsq.GetXaxis().SetRange(0,800);
        c1.SaveAs(foldm+'/'+plot_nm+'_signifsq.png')
    elif drawerr is True:
        sigerr.Draw('colz text')
        sigerr.GetYaxis().SetRange(0,400);
        sigerr.GetXaxis().SetRange(0,800);
        c1.SaveAs(foldm+'/'+plot_nm+'_sigerr.png')

        paintsq=False
        if paintsq is True:
            gStyle.SetPaintTextFormat("2.3g");
            sigerrsq.Draw('colz text')
            sigerrsq.GetYaxis().SetRange(0,400);
            sigerrsq.GetXaxis().SetRange(0,800);
            c1.SaveAs(foldm+'/'+plot_nm+'_sigerrsq.png')

        sigrelerr.Draw('colz text')
        sigrelerr.GetYaxis().SetRange(0,400);
        sigrelerr.GetXaxis().SetRange(0,800);
        c1.SaveAs(foldm+'/'+plot_nm+'_sigrelerr.png')
    
    os.system("cp "+wloc+'/index.php '+foldm)
    
def main():
    for dm in dmass:
        print "MASS:\t", dm
        for reg in regions:
            dmfol = dm
            if "AN" in dm:
                dmfol = dmass[dm]
                print "Get AN masspoints", dmass[dm]
            else:
                dmmin = dmass[dm][0]
                dmmax = dmass[dm][1]

            foldm = folder +relstr + dmfol + '/' + reg
            os.system("mkdir -p "+foldm)


            sigunrol,bkgunrol, sigPtm, bkgPtm, sigMT2, bkgMT2, sig2D, bkg2D, signif2D=define_histos(reg,dm)
            nsig = addhistos(inpfile, dmass,dm,dmmin,dmmax,reg, sigunrol,bkgunrol, sigMT2,bkgMT2,sigPtm,bkgPtm)
            make2D(sigunrol,bkgunrol, sig2D, bkg2D,signif2D, nsig)

            #rangPtm, rangMT2=varbinfromprof()
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
            sigerrsq  = TH2D("sigerr"+reg+dm,"sigerrsq "+vartitle, nvarbinx, varbinx, nvarbiny, varbiny)
            sigerr    = TH2D("sigerr"+reg+dm,"sigerr "+vartitle, nvarbinx, varbinx, nvarbiny, varbiny)
            sigrelerr = TH2D("sigrelerr"+reg+dm,"sigrelerr "+vartitle, nvarbinx, varbinx, nvarbiny, varbiny)
            signifvar = TH2D("signifvar"+reg+dm, signiftitle, nvarbinx, varbinx, nvarbiny, varbiny)
            signifvarsq = TH2D("signifvarsq"+reg+dm, signiftitle, nvarbinx, varbinx, nvarbiny, varbiny)


            fillvarbins(bkg2D,sig2D,bkgvar,sigvar, dorelerr,sigerrsq)

            fillsignif(bkgvar,sigvar,signifvar, nvarPtm, nvarMT2,dorelerr, sigerr, sigerrsq)

            signifvarsq.Multiply(signifvar,signifvar)#,signifvar)
            sigrelerr.Divide(sigerr,sigvar)
            draw_histos(sigvar,bkgvar,signifvar, signifvarsq, foldm, varbin, dorelerr,sigerr, sigerrsq,sigrelerr)
            exit()

    cpweb= 'cp -r '+folder+" "+ optim
    #os.system(cpweb)
    print cpweb    
if __name__ == "__main__":
    main()



