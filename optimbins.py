from make2d import make2D, addhistos, addhisto,fillsignif, fillvarbins
import os, sys
from ROOT import TH1D, TH2D, TFile, TTree, TCanvas, gROOT, gStyle
from array import array
import numpy as np
import time
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
regions = ["combined", "VR1_Tag_sf", "VR1_Tag_em", "VR1_Veto_em", "VR1_Veto_sf"]
dmass   = {"all": [1,700],"ANMP":"mS-450_mX-325","dm_1to125": [1,125], "dm_125to200" : [125,200] , "dm_200to700": [200,700]}
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
        rangPtm = [ 100, 140, 200, 300] #temporary way to assimilate to make2d.py
        rangMT2 = [   0,  20,  40,  60,  80, 100, 120]
        nvarMT2=len(rangMT2)
        nvarPtm=len(rangPtm)
        ptmaxvar=1.6*rangPtm[-1]                                                                                                   
        print ptmaxvar                                                                                                             
        from itertools import combinations


        #Optimise binning for Ptm                                                                                                  
        nvarbinx  = len(rangPtm)-1
        nvarbiny  = len(rangMT2)-1
        varbinx   = array('d',rangPtm)
        varbiny   = array('d',rangMT2)


        nranPtm=(ptmaxvar-rangPtm[0])/wPtm +1
        ptmbins=np.linspace(rangPtm[0],ptmaxvar, nranPtm)
        print '\n'#nranPtm,rangPtm[0],ptmbins                                                                                      
        nbinPtmMax=6
        maxbinsPtm=np.zeros(nbinPtmMax)
        significances={}
        print maxbinsPtm
        #exit()                                                                                                                    
        for nbinPtm in range(1,nbinPtmMax):
            print "i",nbinPtm
            allbincombPtm=(combinations(ptmbins,nbinPtm+1))
            #print "length", len(list(allbincombPtm))                                                                              
            maxsignifsq=0
            postbinning=[]
            for idx,binning in enumerate(allbincombPtm):
                # print "binning",binning, nbinPtm, varbinx, nvarbinx                                                              
                if idx>20: break
                if idx is 0: continue
                varbinxi= array('d',binning)
                bkgvari    = TH2D("bkgvari"+reg+dm,"bkg "   +vartitle, nbinPtm, varbinxi, nvarbiny, varbiny)

                sigvari    = TH2D("sigvari"+reg+dm,"signal "+vartitle, nbinPtm, varbinxi, nvarbiny, varbiny)
                signifvari = TH2D("signifvari"+reg+dm, signiftitle, nbinPtm, varbinxi, nvarbiny, varbiny)
                signifvarsqi = TH2D("signifvarsqi"+reg+dm, signiftitle, nbinPtm, varbinxi, nvarbiny, varbiny)

                fillvarbins(bkg2D,sig2D,bkgvari,sigvari)
                fillsignif(bkgvari,sigvari,signifvari,nvarPtm,nvarMT2)
                signifvarsqi.Multiply(signifvari,signifvari)#,signifvar)                                                           
                signifsq=signifvarsqi.Integral()
                print idx, "\t", signifsq, varbinxi
                if signifsq >  maxsignifsq:
                    maxsignifsq=signifsq
                    postbinning=[]
                if signifsq == maxsignifsq:
                    postbinning.append(binning)
                del bkgvari, sigvari, signifvari, signifvarsqi

            #ee                                                                                                                    
            significances[str(nbinPtm)+'_bins']={'signif':maxsignifsq, 'possible_bins':postbinning}
                #maxbinsPtm[nbinPtm]=maxsignifsq                                                                                   
        print("--- %s seconds ---" % (time.time() - start_time))                                                                   
        import csv
        csv_fol="binning/"
        csv_nm="bestbinningPtm"
        if idx<25: csv_nm+="_test"
        with open(csv_nm+'.csv', mode='w') as employee_file:
            sig_txt = csv.writer(employee_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)

            sig_txt.writerow(['Best binning'])
            sig_txt.writerow(['signif', 'possible binning'])

            for nbin in significances:
                sig_txt.writerow([nbin, significances[nbin]])
                sig_txt.writerow(["--------------------------------"])
                print nbin,significances[nbin],"\n-----------------------------"

        os.system("mkdir -p "+csv_fol)
        os.system("mv "+csv_nm+".csv "+csv_fol)
        print "Best binning saved in", csv_fol+csv_nm+".csv"
        exit()
