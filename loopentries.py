import sys
import optparse
from array import array
import ROOT,os
import numpy as np
ROOT.gROOT.SetBatch(0)
exec(open("definebranches.py").read())
exec(open("applycuts.py").read())

def loopentries(sample,tree, samplenm, idxmax):
    evid, maxSV, nSV, Dnjetstot,isSF,btagW, bvetoW,MET_sumEt,MET_pt,\
    PV_x,PV_y,PV_z,PV_npvs,PV_chi2,SV_eta,SV_phi,SV_pt,SV_mass,\
    SV_x,SV_y,SV_z,SV_chi2, mll,mt2ll,nLepton,ptmiss,susyMstop,susyMLSP= defBranches(tree, samplenm)
    counter=0
    print "looping over", samplenm, "\nfor",idxmax, "events"
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
        #print "crashes here"
        #print type(entry.SV_eta)
        
        mll[0]     = entry.mll
        mt2ll[0]   = entry.mt2ll
        nLepton[0] = entry.nLepton
        ptmiss[0]  = entry.ptmiss
        
        signalopts=["T2tt", "TChipmWW","signal"]
        if(samplenm in signalopts): 
            susyMstop[0] = entry.susyMstop
            susyMLSP[0]  = entry.susyMLSP
        lastSV=min(len(entry.SV_eta),maxSV)
        nSV[0] = len(entry.SV_eta)
        SVentries(entry,nSV[0], SV_eta, SV_phi, SV_pt, SV_mass, SV_x, SV_y, SV_z, SV_chi2)

        tree.Fill()
        
    print "#############################\n"
    del nSV, SV_eta,SV_pt, SV_mass,SV_x,SV_y, SV_z, SV_chi2
    #print nSV
    hfile.cd()
    tree.Write()
    return 



if __name__ == '__main__':

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)

    parser.add_option('--nttbar' , dest='nttbar' , help='# ttbar entries', default=1000000)
    parser.add_option('--nT2tt' , dest='nT2tt' , help='# ttbar\
    entries', default=100000)
    parser.add_option('--signal' , dest='signal' , help='signal file', default=inputFilesignal)
    parser.add_option('--bkg' , dest='bkg' , help='background file', default=inputFilebkg)
    parser.add_option('--output' , dest='output' , help='name of output rootfile', default=hfilenm)
    (opt, args) = parser.parse_args()
    folder='../rootfiles/'
    nttbar=int(opt.nttbar)
    nT2tt=int(opt.nT2tt)
    if('/' in opt.signal): inputFileT2tt=opt.signal
    else: inputFileT2tt=folder+opt.signal
    if('/' in opt.bkg): inputFilettbar=opt.bkg
    else: inputFilettbar=folder+opt.bkg
    if('/' in opt.output): hfilenm=opt.output
    else: hfilenm='../Output/'+opt.output

    print "signal from", inputFilesignal
    hfile = ROOT.TFile(hfilenm,"RECREATE","Example");
    
    signal_t = ROOT.TTree("T2tt","A ROOT tree");
    ttbar_t   = ROOT.TTree("ttbar","A ROOT tree");


    ttbar_f = ROOT.TFile.Open(inputFilettbar ,"READ")
    ttbar_evs= ttbar_f.Get('Events')
    signal_f = ROOT.TFile.Open(inputFilesignal ,"READ")
    signal_evs= signal_f.Get('Events')

    print "nttbar", nttbar
    loopentries(ttbar_evs,ttbar_t, "ttbar",nttbar)

    loopentries(signal_evs, signal_t, "T2tt", nT2tt)
    
    print "Writing tree in", hfilenm

    c1 = ROOT.TCanvas( 'c1', 'Dynamic Filling Example', 200, 10, 1600, 900 )

'''
def makehistos(bkg,signal,name,display=False):
    ymax=1.2*max(bkg.GetMaximum(),signal.GetMaximum())+1
    print "YYYYY MAX", ymax
    bkg.GetYaxis().SetRangeUser(0.,ymax)
    bkg.Draw('e')
    
    signal.SetMarkerColor(2)
    signal.SetLineColor(2)

    signal.Draw('same')
    c1.SaveAs(name+".png")
    #if(display is True): os.system('display '+ name+'.png &')

makehistos(hdnjets_sf_tag_ttbar,hdnjets_sf_tag_T2tt, "hdnjets_sf_tag")#, display=True)
makehistos(hdnjets_df_tag_ttbar,hdnjets_df_tag_T2tt, "hdnjets_df_tag")
makehistos(hdnjets_sf_veto_ttbar,hdnjets_sf_veto_T2tt, "hdnjets_sf_veto")
makehistos(hdnjets_df_veto_ttbar,hdnjets_df_veto_T2tt, "hdnjets_df_veto")
'''
