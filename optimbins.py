import os, sys, time,csv
import numpy as np
from ROOT      import TH1D, TH2D, TFile, TTree, TCanvas, gROOT, gStyle
from itertools import combinations
from array     import array
from make2d    import *
        
wloc = os.environ['WWW']
start_time= time.time()
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
wMT2 = (hmaxMT2-hminMT2)/nMT2
regions = ["combined"]#, "VR1_Tag_sf", "VR1_Tag_em", "VR1_Veto_em", "VR1_Veto_sf"]
dmass   = {"dm_1to200": [1,200]}
#{"all": [1,700],"ANMP":"mS-450_mX-325","dm_1to125": [1,125], "dm_125to200" : [125,200] , "dm_200to700": [200,700]}


binOriPtm = [140, 200, 300] #temporary way to assimilate to make2d.py
binOriMT2 = [   0,  20,  40,  60,  80, 100, 120]
nvarMT2 = len(binOriMT2)
nvarPtm = len(binOriPtm)
maxvarPtm = 400 #1.6*rangPtm[-1]
maxvarMT2 = 200
nbinPtm  = len(binOriPtm)-1
nbinMT2  = len(binOriMT2)-1
vbinPtm  = array('d',binOriPtm)
vbinMT2  = array('d',binOriMT2)
print "MT2",nbinMT2

def optimBins(varOptim):
    global nbinPtm, nbinMT2, vbinPtm, vbinMT2
    #Optimise binning for Ptm
    nbinMax = 6
    maxPtm  = False
    draw    = True
    isPtm   = False
    isMT2   = False
    if "Ptm" in varOptim:
        nbinMax = 6
        maxvar  = maxvarPtm
        binOri  = binOriPtm
        wOri    = wPtm
        isPtm   = True
    elif "MT2" in varOptim:
        nbinMax = 10
        maxvar  = maxvarMT2
        binOri  = binOriMT2
        wOri    = wMT2
        isMT2   = True
    print "max Ptmiss", maxvar, wOri
    nran=(maxvar-binOri[0])/wOri +1
    vbins=np.linspace(binOri[0],maxvar, nran)
    print '\n'#nranPtm,rangPtm[0],ptmbins

    significances={}
    #exit()
    for nbin in range(1,nbinMax):
        if  (isPtm): nbinPtm=nbin
        elif(isMT2): nbinMT2=nbin
        if nbinPtm is not 4: continue
        allbincombPtm = (combinations(vbins,nbin+1))
        maxsignif=0
        postbinning=[]
        for idx,binning in enumerate(allbincombPtm):
            print "idx",idx
            if idx>15: break
            #if idx is not 2: continue
            '''
            if idx is 0:
                #nbinPtm=
                nbinPtm=2
                binning=array('d',[140.,200.0,300.0])
            if idx is 1: binning=array('d',[100.0,160.0,220.0,280.0, 380.0])
            if idx is 2:
                nbinPtm=7
                binning=array('d',[100.0,180.0,240.0,320.0,360.,400.,440., 460.0])
            '''
            varbini= array('d',binning)
            if  (isPtm): vbinPtm=varbini
            elif(isMT2): vbinMT2=varbini
            print "inside", nbinPtm, nbinMT2
            bkgvari      = TH2D("bkgvari"+reg+dm,"bkg "   +vartitle, nbinPtm, vbinPtm, nbinMT2, vbinMT2)
            sigvari      = TH2D("sigvari"+reg+dm,"signal "+vartitle, nbinPtm, vbinPtm, nbinMT2, vbinMT2)
            signifvari   = TH2D("signifvari"+reg+dm, signiftitle   , nbinPtm, vbinPtm, nbinMT2, vbinMT2)
            signifvarsqi = TH2D("signifvarsqi"+reg+dm, signiftitle , nbinPtm, vbinPtm, nbinMT2, vbinMT2)
            fillvarbins(bkg2D,sig2D,bkgvari,sigvari)
            fillsignif(bkgvari,sigvari,signifvari,nbinPtm,nvarMT2)
            signifvarsqi.Multiply(signifvari,signifvari)#,signifvar)                                                           
            signifsq=signifvarsqi.Integral(0,nbinPtm+1,0,nbinMT2+1)
            signif=np.sqrt(signifsq)

            if signif >  maxsignif:
                maxsignif=signif
                postbinning=[]
            if signif > 0.999*maxsignif:
                postbinning.append(binning)
            if idx<3 or signif>0.999*maxsignif :print nbinPtm,idx, "\t", signif, varbini

            foldm="test"
            varbin='test'+str(idx)
            if idx<4 and draw is True: draw_histos(sigvari,bkgvari,signifvari, signifvarsqi, foldm, varbin)
            del bkgvari, sigvari, signifvari, signifvarsqi

        significances[str(nbinPtm)+'_bins']={'signif':maxsignif, 'possible_bins':postbinning}
        print("--- %s seconds ---" % (time.time() - start_time))
    return idx, significances

def write_bestbin(suffix=''):
    csv_fol="binning/"
    csv_nm="bestbinningPtm"+suffix
    if idx<25: csv_nm+="_test"
    with open(csv_nm+'.csv', mode='w') as employee_file:
        sig_txt = csv.writer(employee_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        sig_txt.writerow(['Best binning'])
        sig_txt.writerow(['signif', 'possible binning', '(nbins not including overflow)'])

        for nbin in sorted(significances.iterkeys()):
            sig_txt.writerow([nbin, significances[nbin]])
            sig_txt.writerow(["--------------------------------"])
            print nbin,significances[nbin],"\n-----------------------------"

    os.system("mkdir -p "+csv_fol)
    os.system("mv "+csv_nm+".csv "+csv_fol)
    print "Best binning saved in", csv_fol+csv_nm+".csv"
        


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
        varoptim  = 'Ptm'
        print "before function", nbinPtm, nbinMT2
        idx, significances= optimBins(varoptim)
        write_bestbin()
        exit()
