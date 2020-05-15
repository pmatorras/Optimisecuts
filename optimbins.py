import os, sys, time,csv, optparse
import numpy as np
from ROOT      import TH1D, TH2D, TFile, TTree, TCanvas, gROOT, gStyle
from itertools import combinations
from array     import array
from make2d    import *
from getparams import *
start_time = time.time()


regions = ["combined"]#, "VR1_Tag_sf", "VR1_Tag_em", "VR1_Veto_em", "VR1_Veto_sf"]
dmass   = {"dm_1to200": [1,200]}
#{"all": [1,700],"ANMP":"mS-450_mX-325","dm_1to125": [1,125], "dm_125to200" : [125,200] , "dm_200to700": [200,700]}
binOriPtm = [160.0, 220.0, 280.0, 380.0]
binOriMT2 = [80,100,120]
binANPtm = [140, 200, 300]
binANMT2 = [ 80, 100, 120]

print "--------------------------------------------------"
print "VARIABLE TO OPTIMISE:\t", varOptim
print "MC FOR YEAR:\t\t", year
str_AN = ''
if opt.AN is True:
    str_AN = '_ANbin'
    print "USING AN BINNING"
    binOriPtm=binANPtm
    binOriMT2=binANMT2
print "--------------------------------------------------"


nvarMT2 = len(binOriMT2)
nvarPtm = len(binOriPtm)
maxvarPtm = 400 #1.6*rangPtm[-1]
maxvarMT2 = 200
nbinPtm  = len(binOriPtm)-1
nbinMT2  = len(binOriMT2)-1
vbinPtm  = array('d',binOriPtm)
vbinMT2  = array('d',binOriMT2)
print "MT2",nbinMT2

#Draw test distributions
def testdistrib(idx,isMT2, isPtm, binning, nbinMT2, nbinPtm):
    if idx is 0 :
        if isMT2:
            binning=array('d',binANMT2)
            nbinMT2=len(binANMT2)-1
        elif isPtm:
            binning=array('d',binANPtm)
            nbinPtm=len(binANPtm)-1
    if idx is 1 :
        if isMT2:
            binning = array('d',[90.,100.,160.])
            nbinMT2 = 2
        elif isPtm:
            binning=array('d',[100.0,160.0,220.0,280.0, 380.0])
            nbinPtm=4
    if idx is 2:
        if isMT2:
            emptySlot=True
        if isPtm:
            nbinPtm=7
            binning=array('d',[100.0,180.0,240.0,320.0,360.,400.,440., 460.0])
    return nbinMT2,nbinPtm,binning

#Optimise binning for Ptm
def optimBins(varOptim, vartitle):
    global nbinPtm, nbinMT2, vbinPtm, vbinMT2
    nbinMax = 6
    maxPtm  = False
    draw    = True
    isPtm   = False
    isMT2   = False
    nhistos = 0
    hnm     = varOptim
    if "Ptm" in varOptim:
        nbinMax = 6
        maxvar  = maxvarPtm
        binOri  = binOriPtm
        wOri    = wPtm
        isPtm   = True
    elif "MT2" in varOptim:
        nbinMax = 6
        maxvar  = maxvarMT2
        binOri  = binOriMT2
        wOri    = wMT2
        isMT2   = True
    vartitle+=hnm
    print "max possible bin: ", maxvar, wOri
    nran=(maxvar-binOri[0])/wOri +1
    vbins=np.linspace(binOri[0],maxvar, nran)
    print '\n'#nranPtm,rangPtm[0],ptmbins

    significances={}
    #exit()
    for nbin in range(1,nbinMax):
        if  (isPtm): nbinPtm=nbin
        elif(isMT2): nbinMT2=nbin
        if Test is True and nbin is not 4: continue
        allbincombPtm = (combinations(vbins,nbin+1))
        maxsignif=0
        postbinning=[]
        for idx,binning in enumerate(allbincombPtm):
            #if binning[0]>20: continue
                
            if Test is True:
                if idx>1:break
                nbinMT2,nbinPtm,binning =testdistrib(idx,isMT2, isPtm, binning, nbinMT2, nbinPtm)

            varbini= array('d',binning)
            if  (isPtm): vbinPtm=varbini
            elif(isMT2): vbinMT2=varbini
            bkgvari      = TH2D("bkgvar"+hnm+reg+dm,"bkg "   +vartitle, nbinPtm, vbinPtm, nbinMT2, vbinMT2)
            sigvari      = TH2D("sigvar"+hnm+reg+dm,"signal "+vartitle, nbinPtm, vbinPtm, nbinMT2, vbinMT2)
            signifvari   = TH2D("signifvar"+hnm+reg+dm, signiftitle   , nbinPtm, vbinPtm, nbinMT2, vbinMT2)
            signifvarsqi = TH2D("signifvarsq"+hnm+reg+dm, signiftitle , nbinPtm, vbinPtm, nbinMT2, vbinMT2)
            fillvarbins(bkg2D,sig2D,bkgvari,sigvari)
            fillsignif(bkgvari,sigvari,signifvari,nbinPtm,nvarMT2)
            signifvarsqi.Multiply(signifvari,signifvari)#,signifvar)
            
            signifsq = signifvarsqi.Integral(0,nbinPtm+1,0,nbinMT2+1)
            signif   = np.sqrt(signifsq)
            if signif > 1.001*maxsignif:  postbinning=[]
            if signif >  maxsignif:
                maxsignif=signif
                print "NEW MAX"
            if signif > 0.999*maxsignif:
                postbinning.append(binning)
            if idx<3 or signif>0.999*maxsignif :
                sigstr = str(round(signif,5))
                printoutput = str(nbin)+' '+str(idx)+ "\t"+ sigstr+'\t '+str(binning)
                os.system('echo "'+printoutput+ '">>output'+varOptim+".log")
                print printoutput
            foldm  = "test"
            varbin = 'test'+str(idx)
            if Test is True and draw is True: draw_histos(sigvari,bkgvari,signifvari, signifvarsqi, foldm, varbin)
            del bkgvari, sigvari, signifvari, signifvarsqi
        nhistos=idx
        if isPtm: 
            significances[str(nbin)+'_bins']={'signif':round(maxsignif,4), 'possible_Ptmiss_bins':postbinning}
        elif isMT2:
            significances[str(nbin)+'_bins']={'signif':round(maxsignif,4), 'possible_MT2_bins':postbinning}
            
        print("--- %s seconds ---" % (time.time() - start_time))
    return nhistos, significances

def write_bestbin(suffix=''):
    csv_fol="binning/"
    csv_nm="bestbinning"+suffix
    if idx<25: csv_nm+="_test"
    binline=''
    str_AN=''
    if "AN"  in suffix: str_AN='\t(AN binning)'
    if "Ptm" in suffix: binline="MT2 binning:\t"+str(binOriMT2)
    elif "MT2" in suffix: binline="Ptmiss binning:\t"+str(binOriPtm)
    with open(csv_nm+'.csv', mode='w') as employee_file:
        sig_txt = csv.writer(employee_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        sig_txt.writerow(["----------------------------------------------"])
        sig_txt.writerow(['BEST BINNING to maximise significance'])
        sig_txt.writerow([binline+str_AN])
        sig_txt.writerow(['n bins', 'signif.', 'possible binning (nbins not including overflow)'])
        sig_txt.writerow(["----------------------------------------------"])

        for nbin in sorted(significances.iterkeys()):
            sig_txt.writerow([nbin, significances[nbin]])
            sig_txt.writerow(["--------------------------------"])
            print nbin,significances[nbin],"\n-----------------------------"

    os.system("mkdir -p "+csv_fol)
    os.system("mv "+csv_nm+".csv "+csv_fol)
    print "Best binning saved in", csv_fol+csv_nm+".csv"
        




    

#Loop over binning
for dm in dmass:
    print "MASS:\t", dm
    for reg in regions:
        vartitle    = "("+reg+"-"+dm+"); p_{T}^{miss}; m_{T2}^{ll}"
        signiftitle = "s/#sqrt{s+b} "+vartitle
        signifsqtit = "significance squared"+vartitle

        dmfol = dm
        if "AN" in dm:
            dmfol = dmass[dm]
            print "Get AN masspoints", dmass[dm]
        else:
            dmmin = dmass[dm][0]
            dmmax = dmass[dm][1]

        foldm = folder + '/'+ dmfol + '/' + reg
        os.system("mkdir -p "+foldm)

        #Define necessary histograms                                                                                               
        print "REGION:", reg
        sigunrol = TH1D("sigunrol"+reg+dm, reg, 10000, 0, 10000)
        bkgunrol = TH1D("bkgunrol"+reg+dm, reg, 10000, 0, 10000)
        sigPtm = TH1D("sigPtmiss" +reg+dm, reg, nPtm, 0, hmaxPtm)
        bkgPtm = TH1D("bkgPtmiss" +reg+dm, reg, nPtm, 0, hmaxPtm)
        sigMT2 = TH1D("sigMT2"    +reg+dm, reg, nMT2, 0, hmaxMT2)
        bkgMT2 = TH1D("bkgMT2"    +reg+dm, reg, nMT2, 0, hmaxMT2)

        print "2D", nMT2, hminMT2, hmaxMT2
        sig2D    = TH2D("sig2D"   +reg+dm, reg, nPtm, hminPtm, hmaxPtm,nMT2, hminMT2,hmaxMT2)
        bkg2D    = TH2D("bkgs2D"  +reg+dm, reg, nPtm, hminPtm, hmaxPtm,nMT2, hminMT2,hmaxMT2)
        signif2D = TH2D("signif2D"+reg+dm, reg, nPtm, hminPtm, hmaxPtm,nMT2, hminMT2,hmaxMT2)
        mSmin = 400
        mSmax = 700
        signif2D.SetTitle( "m_{T2ll} vs p_{T}^{miss} ("+reg+'-'+dm+")")
	signif2D.SetYTitle("m_{T2ll} [GeV]")
        signif2D.SetXTitle("p_{T}^{miss} [GeV]")



        nsig = addhistos(inpfile, dmass,dm,dmmin,dmmax,reg, sigunrol,bkgunrol, sigMT2,bkgMT2,sigPtm,bkgPtm)
        make2D(sigunrol,bkgunrol, sig2D, bkg2D,signif2D, nsig)
        idx, significances= optimBins(varOptim, vartitle)
        write_bestbin(varOptim+str_AN)
        exit()
