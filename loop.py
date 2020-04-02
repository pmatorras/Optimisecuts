import sys
import optparse
from array import array
import ROOT,os
import numpy as np
ROOT.gROOT.SetBatch(0)
exec(open("branches.py").read())
exec(open("cuts.py").read())
#calculate delta r lepton-jet
def dphi_lepjet(lepphi, jetphi1, jetphi2):
    dphi1=abs(lepphi-jetphi1)
    dphi2=abs(lepphi-jetphi2)
    if(jetphi1>-999):
        if(jetphi2>-999):
            dphimin=min(dphi1,dphi2)
            dphimax=max(dphi1,dphi2)
        else:
            dphimin=dphi1
            dphimax=dphi1
    else:
        dphimin=-9
        dphimax=-9
    if(dphimax>np.pi):
        dphimax=2*np.pi-dphimax
    if(dphimin>np.pi): dphimin=2*np.pi-dphimin

    return dphimin, dphimax



def loopentries(sample,tree, samplenm, idxmax):
    counter=0
    print "looping over", samplenm, "\nfor",idxmax, "events"
    evid, maxSV, nSV, njets, nbjets,nCleanjets, nbCleanjets, Dnjetstot, Dnbjetstot, isSF,btagW, bvetoW,\
    PV_x,PV_y,PV_z,PV_npvs,PV_chi2,SV_eta,SV_phi,SV_pt,SV_mass,SV_x,SV_y,SV_z,SV_chi2,\
    MET_sumEt,MET_pt,mll,mt2ll,nLepton,ptmiss,susyMstop,susyMLSP,ISRcut,detall,dRll,detajj,dRjj,\
    dphill,dphijj, dphil1jmin,dphil1jmax,dphil1bmin,dphil1bmax, dphil2jmin,dphil2jmax,dphil2bmin,dphil2bmax,\
    lep1_pt,lep1_eta,lep1_phi,lep2_pt,lep2_eta,lep2_phi, jet1_pt,jet1_eta,jet1_phi,jet2_pt,jet2_eta,jet2_phi,\
    bjet1_pt,bjet1_eta,bjet1_phi,bjet2_pt,bjet2_eta,bjet2_phi= defBranches(tree, samplenm)


    def nBjets(pt, eta, phi,btag, njets, cleanidx):
        nbjets = 0
        samelen= bool(len(pt)== len(btag))
        ijets=range(0,njets)
        if samelen is False: ijets= cleanidx
        for i in range(0,njets):#ijets:
            if(pt[i]>20 and abs(eta[i])<2.4 and btag[ijets[i]]>=0.6321):
                nbjets+=1
                if samelen is True: continue
                if i is 0:
                    #bjet1_pt[0]  = pt[0] #Check leading pttagged option
                    bjet1_eta[0] = eta[0]
                    bjet1_phi[0] = phi[0]
                elif i is 1:
                    #bjet2_pt[0]  = pt[1] #Check leading pttagged option
                    bjet2_eta[0] = eta[1]
                    bjet2_phi[0] = phi[1]
          
        return nbjets
    

    
    for idx, entry in enumerate(sample):
        if(idx>idxmax):
            print "reached max number of events"
            break
        #apply normal cuts
        nCleanJet = entry.nCleanJet
        cleanidx  = entry.CleanJet_jetIdx
        if nCleanJet<1: continue
        
        
        sameflavour = False
        df = False
        btag  = False
        bveto = False
        passedcut,isrcut=defaultcuts(entry,samplenm)
        if passedcut is False:
            continue
        ISRcut[0]=int(isrcut)
        
        CleanJet_pt = entry.CleanJet_pt
        nbCleanjets[0] = nBjets(entry.CleanJet_pt, entry.CleanJet_eta, entry.CleanJet_phi,entry.Jet_btagDeepB, nCleanJet, cleanidx)
        nbjets[0]      = nBjets(entry.Jet_pt, entry.Jet_eta, entry.Jet_phi, entry.Jet_btagDeepB, entry.nJet, cleanidx)
        Dnbjetstot[0]  = nbjets[0]-nbCleanjets[0]
        njets[0]       = entry.nJet
        nCleanjets[0]  = nCleanJet
        Dnjetstot[0]   = entry.nJet-nCleanJet

        #print "dnjets", Dnbjetstot[0]
        
        #Bveto
        btagW[0]=-1
        bvetoW[0]=-1
        isSF[0], btagW[0],bvetoW[0]=flavour_tag(entry, samplenm)
        if isSF[0] ==-999 : continue

        counter+=1
        print "event", idx, "total done", counter 
        evid[0]=idx
        MET_sumEt[0] = entry.MET_sumEt
        MET_pt[0]    = entry.MET_pt
        
        
        PV_x[0]    = entry.PV_x
        PV_y[0]    = entry.PV_y
        PV_z[0]    = entry.PV_z
        PV_npvs[0] = entry.PV_npvs
        PV_chi2[0] = entry.PV_chi2
        
        mll[0]     = entry.mll
        mt2ll[0]   = entry.mt2ll
        nLepton[0] = entry.nLepton
        ptmiss[0]  = entry.ptmiss

        lep1_pt[0]  = entry.Lepton_pt[0]
        lep1_eta[0] = entry.Lepton_eta[0]
        lep1_phi[0] = entry.Lepton_phi[0]
        lep2_pt[0]  = entry.Lepton_pt[1]
        lep2_eta[0] = entry.Lepton_eta[1]
        lep2_phi[0] = entry.Lepton_phi[1]

        dphill[0]  = abs(entry.Lepton_phi[0]-entry.Lepton_phi[1])
        detall[0]  = abs(entry.Lepton_eta[0]-entry.Lepton_eta[1])
        if(dphill[0]>np.pi): dphill[0]=2*np.pi-dphill[0]

        jet1_pt[0]  = entry.CleanJet_pt[0]
        jet1_eta[0] = entry.CleanJet_eta[0]
        jet1_phi[0] = entry.CleanJet_phi[0]
        
        if(nCleanJet>1):
            dphijj[0]  = abs(entry.CleanJet_phi[0]-entry.CleanJet_phi[1])
            detajj[0]  = abs(entry.CleanJet_eta[0]-entry.CleanJet_eta[1])
            jet2_pt[0]  = entry.CleanJet_pt[1]
            jet2_eta[0] = entry.CleanJet_eta[1]
            jet2_phi[0] = entry.CleanJet_phi[1]

        if(dphijj[0]>np.pi): dphijj[0]=2*np.pi-dphijj[0]
        dRll[0] = np.sqrt(dphill[0]*dphill[0]+detall[0]*detall[0])
        dRjj[0] = np.sqrt(dphijj[0]*dphijj[0]+detajj[0]*detajj[0])
        

        signalopts=["T2tt", "TChipmWW","signal"]
        if(samplenm in signalopts): 
            susyMstop[0] = entry.susyMstop
            susyMLSP[0]  = entry.susyMLSP
            bjet1_pt[0]   = entry.leadingPtTagged
            bjet2_pt[0]   = entry.trailingPtTagged
        else:
            bjet1_pt[0]   = entry.leadingPtTagged_btagDeepBM_1c
            bjet2_pt[0]   = entry.trailingPtTagged_btagDeepBM_1c

        #dRll,dRjj
        #dRjj=(

        nSV[0] = len(entry.SV_eta)
        SVentries(entry,nSV[0], SV_eta, SV_phi, SV_pt, SV_mass, SV_x, SV_y, SV_z, SV_chi2)
        
        dphil1jmin[0],dphil1jmax[0]=dphi_lepjet(lep1_phi[0], jet1_phi[0],jet2_phi[0])
        dphil2jmin[0],dphil2jmax[0]=dphi_lepjet(lep2_phi[0], jet1_phi[0],jet2_phi[0])
        dphil1bmin[0],dphil1bmax[0]=dphi_lepjet(lep1_phi[0], bjet1_phi[0],bjet2_phi[0])
        dphil2bmin[0],dphil2bmax[0]=dphi_lepjet(lep2_phi[0], bjet1_phi[0],bjet2_phi[0])
        
        tree.Fill()
        
    print "#############################\n"
    del nSV, SV_eta,SV_pt, SV_mass,SV_x,SV_y, SV_z, SV_chi2
    print tree
    hfile.cd()
    tree.Write()
    return 



if __name__ == '__main__':

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)

    parser.add_option('--nT2tt' , dest='nT2tt' , help='# T2tt entries' , default=100000)
    parser.add_option('--nttbar', dest='nttbar', help='# ttbar entries', default=1000000)
    parser.add_option('--nWW' , dest='nWW' , help='# ww entries' , default=1000000)
    parser.add_option('--nttZ', dest='nttZ', help='# ttZ entries', default=1000000)
    parser.add_option('--nZZ' , dest='nZZ' , help='# ZZ entries' , default=1000000)

    parser.add_option('--signal', dest='signal', help='signal file', default=inputFileT2tt)
    parser.add_option('--ttbar' , dest='ttbar' , help='ttbar file' , default=inputFilettbar)
    parser.add_option('--ttZ', dest='ttZ', help='ttZ file', default=inputFilettZ)
    parser.add_option('--WW' , dest='WW' , help='WW file' , default=inputFileWW)
    parser.add_option('--ZZ' , dest='ZZ' , help='ZZ file' , default=inputFileZZ)
    
    parser.add_option('--output' , dest='output' , help='name of output rootfile', default=cutfilenm)
    parser.add_option('--only' , dest='only' , help='if only one sample wants to be looked', default='no')
    parser.add_option('--except' , dest='ignore' , help='not run one of the samples', default='no')
    

    parser.add_option('--sample'    , dest='sample'    , help='sample to run'    , default='no')
    parser.add_option('--samplefile', dest='samplefile', help='file sample'      , default=inputFileT2tt)
    parser.add_option('--nsample'   , dest='nsample'   , help='# sample entries' , default=1000000)

    parser.add_option('--test' , dest='test' , help='run test on subset', default=False, action='store_true')
    (opt, args) = parser.parse_args()
    print cutfolder, "CUUT"
    os.system("mkdir -p "+ cutfolder)
    folder = './rootfiles/'
    nSam   = int(opt.nsample)
    nttbar = int(opt.nttbar)
    nT2tt  = int(opt.nT2tt)
    nWW    = int(opt.nWW)
    nZZ    = int(opt.nZZ)
    nttZ   = int(opt.nttZ)
    if(opt.test):
        nSam   =  1000
        nT2tt  =  1000
        nttbar = 10000
        nWW    = 10000
        nZZ    =  2000
        nttZ   =  2000
    if('/' in opt.signal): inputFileT2tt=opt.signal
    else:
        inputFileT2tt=folder+opt.signal
    if('/' in opt.ttbar): inputFilettbar=opt.ttbar
    else: inputFilettbar=folder+opt.ttbar
    if('/' in opt.WW): inputFileWW=opt.WW
    else: inputFileWW=folder+opt.WW
    
    if('/' in opt.output): cutfilenm=opt.output
    else: 
        outfol='./Output/'
        cutfilenm=outfol+opt.output
        os.System("mkdir -p"+outfol)

    print "signal from", inputFileT2tt




    if opt.sample is not "no":
        sample=opt.sample
        print "sample is", sample
        tree_sam = ROOT.TTree(sample,"tree with"+sample+"bkg")
        fileloc  = inputFileT2tt#opt.samplefile
        file_sam = ROOT.TFile.Open(fileloc ,"READ")        
        Evs_sam  = file_sam.Get('Events')
        if(opt.test is True): outfol=cutfolder 
        else: outfol="/afs/cern.ch/work/p/pmatorra/private/CMSSW_10_2_14/src/Optimisecuts/Output/"
        if(os.path.exists(outfol) is False): os.system("mkdir -p "+outfol)
        cutfilenm= outfol+"cuts_"+sample+".root"
        hfile = ROOT.TFile(cutfilenm,"RECREATE","Example");

        loopentries(Evs_sam,tree_sam,opt.sample,nSam)

    else:
        hfile = ROOT.TFile(cutfilenm,"RECREATE","Example");

        signal_t = ROOT.TTree("T2tt","tree with signal");
        ttbar_t  = ROOT.TTree("ttbar","tree with ttbar bkg");
        WW_t     = ROOT.TTree("WW","tree with WW bkg")
        ZZ_t     = ROOT.TTree("ZZ","tree with WW bkg")
        ttZ_t    = ROOT.TTree("ttZ","tree with WW bkg")
        
        ttbar_f  = ROOT.TFile.Open(inputFilettbar ,"READ")
        WW_f     = ROOT.TFile.Open(inputFileWW ,"READ")
        ZZ_f     = ROOT.TFile.Open(inputFileZZ ,"READ")
        ttZ_f    = ROOT.TFile.Open(inputFilettZ,"READ")
        signal_f = ROOT.TFile.Open(inputFileT2tt ,"READ")
        
        WW_evs     = WW_f.Get('Events')
        ZZ_evs     = ZZ_f.Get('Events')
        ttZ_evs    = ttZ_f.Get('Events')
        ttbar_evs  = ttbar_f.Get('Events')
        signal_evs = signal_f.Get('Events')

        
        samples= {'ZZ'   :[inputFileZZ, ZZ_t, nZZ],
                  'ttZ'  :[inputFilettZ, ttZ_t, nttZ],
                  'ttbar': [inputFilettbar,ttbar_t, nttbar],
                  'WW'   : [inputFileWW, WW_t, nWW],
                  'T2tt':[inputFileT2tt, signal_t,nT2tt]}
        #print inputFilettbar

        rmsamples=opt.ignore.split('_')
        dosamples=opt.only.split('_')
        print dosamples, rmsamples
        for sample in samples:
            if sample in rmsamples:
                print "Sample", sample, "in --except option"
                continue
            if sample not in dosamples and "no" not in dosamples[0]:
                print "Sample", sample, "\tnot in --only option"
                continue
            if opt.sample not in "no" and sample not in opt.sample: 
                print "Sample", sample, "not in ", opt.sample
                continue

            print samples[sample][0] 
            tree   = samples[sample][1]
            tfile  = ROOT.TFile.Open(samples[sample][0],"READ")
            events = tfile.Get('Events')
            loopentries(events,tree, sample, samples[sample][2])



    print "Writing tree in", cutfilenm
    
