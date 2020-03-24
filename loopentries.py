import sys
import optparse
from array import array
import ROOT,os
import numpy as np
ROOT.gROOT.SetBatch(0)
exec(open("definebranches.py").read())
exec(open("applycuts.py").read())

def loopentries(sample,tree, samplenm, idxmax):
    counter=0
    print "looping over", samplenm, "\nfor",idxmax, "events"
    evid, maxSV, nSV, Dnjetstot,isSF,btagW, bvetoW,MET_sumEt,MET_pt,\
    PV_x,PV_y,PV_z,PV_npvs,PV_chi2,SV_eta,SV_phi,SV_pt,SV_mass,\
    SV_x,SV_y,SV_z,SV_chi2, mll,mt2ll,nLepton,ptmiss,susyMstop,susyMLSP,\
    lep1_pt,lep2_pt,jet1_pt,jet2_pt,dphill,detall,dRll,dphijj,detajj,dRjj= defBranches(tree, samplenm)
    for idx, entry in enumerate(sample):
        if(idx>idxmax):
            print "reached max number of events"
            break
        #apply normal cuts
        if entry.nCleanJet<1: continue
        
        
        sameflavour = False
        df = False
        btag  = False
        bveto = False
        passedcut=defaultcuts(entry,samplenm) 
        if passedcut is False:
            continue
    
        CleanJet_pt = entry.CleanJet_pt
        ncleanJet   = entry.nCleanJet
        Dnjets=entry.nJet-ncleanJet
        Dnjetstot[0]= Dnjets
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

        lep1_pt[0] = entry.Lepton_pt[0]
        lep2_pt[0] = entry.Lepton_pt[1]

        dphill[0]  = abs(entry.Lepton_phi[0]-entry.Lepton_phi[1])
        detall[0]  = abs(entry.Lepton_eta[0]-entry.Lepton_eta[1])
        dphijj[0]  = abs(entry.CleanJet_phi[0]-entry.CleanJet_phi[1])
        detajj[0]  = abs(entry.CleanJet_eta[0]-entry.CleanJet_eta[1])

        if(dphill[0]>np.pi): dphill[0]=2*np.pi-dphill[0]
        if(dphijj[0]>np.pi): dphijj[0]=2*np.pi-dphijj[0]
        dRll[0] = np.sqrt(dphill[0]*dphill[0]+detall[0]*detall[0])
        dRjj[0] = np.sqrt(dphijj[0]*dphijj[0]+detajj[0]*detajj[0])
        

        signalopts=["T2tt", "TChipmWW","signal"]
        if(samplenm in signalopts): 
            susyMstop[0] = entry.susyMstop
            susyMLSP[0]  = entry.susyMLSP
            jet1_pt[0]   = entry.leadingPtTagged
            jet2_pt[0]   = entry.trailingPtTagged
        else:
            jet1_pt[0]   = entry.leadingPtTagged_btagDeepBM_1c
            jet2_pt[0]   = entry.trailingPtTagged_btagDeepBM_1c

        #dRll,dRjj
        #dRjj=(
        
        lastSV=min(len(entry.SV_eta),maxSV)
        nSV[0] = len(entry.SV_eta)
        SVentries(entry,nSV[0], SV_eta, SV_phi, SV_pt, SV_mass, SV_x, SV_y, SV_z, SV_chi2)

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

    parser.add_option('--nttbar' , dest='nttbar' , help='# ttbar entries', default=1000000)
    parser.add_option('--nWW' , dest='nWW' , help='# ww entries', default=1000000)
    parser.add_option('--nT2tt' , dest='nT2tt' , help='# T2tt entries', default=100000)

    parser.add_option('--signal' , dest='signal' , help='signal file', default=inputFilesignal)
    parser.add_option('--ttbar' , dest='ttbar' , help='background file', default=inputFilettbar)
    parser.add_option('--WW' , dest='WW' , help='background file', default=inputFileWW)
    
    parser.add_option('--output' , dest='output' , help='name of output rootfile', default=hfilenm)
    parser.add_option('--only' , dest='only' , help='if only one sample wants to be looked', default='no')
    parser.add_option('--except' , dest='ignore' , help='not run one of the samples', default='no')
    
    parser.add_option('--test' , dest='test' , help='run test on subset', default=False, action='store_true')
    (opt, args) = parser.parse_args()

    folder = '../rootfiles/'
    nttbar = int(opt.nttbar)
    nT2tt  = int(opt.nT2tt)
    nWW    = int(opt.nWW)
    if(opt.test):
        nttbar=10000
        nT2tt=1000
        nWW=10000
    if('/' in opt.signal): inputFileT2tt=opt.signal
    else:
        inputFileT2tt=folder+opt.signal
    if('/' in opt.ttbar): inputFilettbar=opt.ttbar
    else: inputFilettbar=folder+opt.ttbar
    if('/' in opt.WW): inputFileWW=opt.WW
    else: inputFileWW=folder+opt.WW
    
    if('/' in opt.output): hfilenm=opt.output
    else: hfilenm='../Output/'+opt.output

    print "signal from", inputFileT2tt
    hfile = ROOT.TFile(hfilenm,"RECREATE","Example");
    
    signal_t = ROOT.TTree("T2tt","tree with signal");
    ttbar_t  = ROOT.TTree("ttbar","tree with ttbar bkg");
    WW_t     = ROOT.TTree("WW","tree with WW bkg")

    ttbar_f  = ROOT.TFile.Open(inputFilettbar ,"READ")
    WW_f     = ROOT.TFile.Open(inputFileWW ,"READ")
    signal_f = ROOT.TFile.Open(inputFileT2tt ,"READ")
    
    WW_evs     = WW_f.Get('Events')
    ttbar_evs  = ttbar_f.Get('Events')
    signal_evs = signal_f.Get('Events')

    samples= {'ttbar': [inputFilettbar,ttbar_t, nttbar],
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
        print samples[sample][0] 
        tree   = samples[sample][1]
        tfile  = ROOT.TFile.Open(samples[sample][0],"READ")
        events = tfile.Get('Events')
        loopentries(events,tree, sample, samples[sample][2])
    
    print "Writing tree in", hfilenm
    
