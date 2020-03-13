import gc
from array import array
import ROOT,os
import numpy as np
ROOT.gROOT.SetBatch(0)
exec(open("definebranches.py").read())
exec(open("applycuts.py").read())

inputFilettbar  = '../rootfiles/nanoLatino_TTTo2L2Nu__part10.root'
inputFilesignal = '../rootfiles/nanoLatino_T2tt__mStop-400to1200__part0.root'

hdnjets_sf_veto_ttbar = ROOT.TH1F( 'hdnjets_sf_veto_ttbar', '#Delta(njet-ncleanjets)', 12, -1, 10 )
hdnjets_df_veto_ttbar = ROOT.TH1F( 'hdnjets_df_veto_ttbar', '#Delta(njet-ncleanjets)', 12, -1, 10 )
hdnjets_sf_tag_ttbar  = ROOT.TH1F( 'hdnjets_sf_tag_ttbar' , '#Delta(njet-ncleanjets)', 12, -1, 10 )
hdnjets_df_tag_ttbar  = ROOT.TH1F( 'hdnjets_df_tag_ttbar' , '#Delta(njet-ncleanjets)', 12, -1, 10 )
hdnjets_sf_veto_T2tt = ROOT.TH1F( 'hdnjets_sf_veto_T2tt', '#Delta(njet-ncleanjets)', 12, -1, 10 )
hdnjets_df_veto_T2tt = ROOT.TH1F( 'hdnjets_df_veto_T2tt', '#Delta(njet-ncleanjets)', 12, -1, 10 )
hdnjets_sf_tag_T2tt  = ROOT.TH1F( 'hdnjets_sf_tag_T2tt' , '#Delta(njet-ncleanjets)', 12, -1, 10 )
hdnjets_df_tag_T2tt  = ROOT.TH1F( 'hdnjets_df_tag_T2tt' , '#Delta(njet-ncleanjets)', 12, -1, 10 )

hfilenm="postcuts.root"
hfile = ROOT.TFile(hfilenm,"RECREATE","Example");
#tree  = ROOT.TTree("tree", "for both")
signal_t = ROOT.TTree("signal","A ROOT tree");
ttbar_t   = ROOT.TTree("bkgs","A ROOT tree");

def loopentries(sample,hdnjets_sf_veto, hdnjets_df_veto, hdnjets_sf_tag, hdnjets_df_tag,tree, samplenm, idxmax=1000):
    evid, maxSV, nSV, Dnjetstot,isSF,btagW, bvetoW,MET_sumEt,MET_pt,\
        PV_x,PV_y,PV_z,PV_npvs,PV_chi2,SV_eta,SV_phi,SV_pt,SV_mass,\
        SV_x,SV_y,SV_z,SV_chi2, mll,mt2ll,nLepton,ptmiss,susyMstop,susyMLSP= defBranches(tree, samplenm)
    counter=0
    print "looping over", samplenm, "\nfor",idxmax, "events"
    for idx, entry in enumerate(sample):
        if(idx>idxmax): break
        #apply normal cuts
        if entry.nCleanJet<1: continue
        
        
        sameflavour = False
        df = False
        btag  = False
        bveto = False
        passedcut=defaultcuts(entry) 
        if passedcut is False:
            continue
    
        CleanJet_pt = entry.CleanJet_pt
        ncleanJet   = entry.nCleanJet
        Dnjets=entry.nJet-ncleanJet
        Dnjetstot[0]= Dnjets
        #Bveto
        btagW[0]=-1
        bvetoW[0]=-1
        sameflavour, bvetoW[0]=flavour_tag(entry, False)
        if sameflavour is None : continue
        if sameflavour is True :
            hdnjets_sf_veto.Fill(Dnjets, bvetoW[0])
            print "sf"
            isSF[0]=1
        if sameflavour is False:
            hdnjets_df_veto.Fill(Dnjets, bvetoW[0])
            isSF[0]=-1
            #Btag
        sameflavour, btagW[0]=flavour_tag(entry, True)
        if sameflavour is None : continue
        if sameflavour is True : hdnjets_sf_tag.Fill(Dnjets, btagW[0])
        if sameflavour is False: hdnjets_df_tag.Fill(Dnjets, btagW[0])
        #print "dnjets", Dnjetstot+
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
            print "SIIIIIIIIIIIIGNAL", susyMstop[0]
        
        lastSV=min(len(entry.SV_eta),maxSV)
        nSV[0] = len(entry.SV_eta)
        SVentries(entry,nSV[0], SV_eta, SV_phi, SV_pt, SV_mass, SV_x, SV_y, SV_z, SV_chi2)

        tree.Fill()
        
    print "#############################\n", nSV, SV_eta
    del nSV, SV_eta,SV_pt, SV_mass,SV_x,SV_y, SV_z, SV_chi2
    #print nSV
    hfile.cd()
    tree.Write()

    return hdnjets_sf_veto, hdnjets_df_veto, hdnjets_sf_tag, hdnjets_df_tag

ttbar = ROOT.TFile.Open(inputFilettbar ,"READ")
ttbar_evs= ttbar.Get('Events')
signal = ROOT.TFile.Open(inputFilesignal ,"READ")
signal_evs= signal.Get('Events')

hdnjets_sf_veto_ttbar, hdnjets_df_veto_ttbar, hdnjets_sf_tag_ttbar, hdnjets_df_tag_ttbar= loopentries(ttbar_evs, hdnjets_sf_veto_ttbar, hdnjets_df_veto_ttbar, hdnjets_sf_tag_ttbar, hdnjets_df_tag_ttbar,ttbar_t, "ttbar",10000)
#gc.collect()
#hfile.cd()
#tree.Write()

hdnjets_sf_veto_T2tt, hdnjets_df_veto_T2tt, hdnjets_sf_tag_T2tt, hdnjets_df_tag_T2tt= loopentries(signal_evs,hdnjets_sf_veto_T2tt, hdnjets_df_veto_T2tt, hdnjets_sf_tag_T2tt, hdnjets_df_tag_T2tt, signal_t, "signal", 1000)

print "Writing tree in", hfilenm

c1 = ROOT.TCanvas( 'c1', 'Dynamic Filling Example', 200, 10, 1600, 900 )

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
'''
makehistos(hdnjets_sf_tag_ttbar,hdnjets_sf_tag_T2tt, "hdnjets_sf_tag")#, display=True)
makehistos(hdnjets_df_tag_ttbar,hdnjets_df_tag_T2tt, "hdnjets_df_tag")
makehistos(hdnjets_sf_veto_ttbar,hdnjets_sf_veto_T2tt, "hdnjets_sf_veto")
makehistos(hdnjets_df_veto_ttbar,hdnjets_df_veto_T2tt, "hdnjets_df_veto")
'''
