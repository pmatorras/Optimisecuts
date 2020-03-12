import gc
from array import array
import ROOT,os
import numpy as np
ROOT.gROOT.SetBatch(0)
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



#Apply current cuts
def applycuts(entry):
    passedcut=False
    if(entry.MET_pt >300 #METCUT
       and (entry.Lepton_isTightElectron_cutBasedMediumPOG[0] + entry.Lepton_isTightMuon_mediumRelIsoTight[0]+entry.Lepton_isTightElectron_cutBasedMediumPOG[1]+entry.Lepton_isTightMuon_mediumRelIsoTight[1])==2 #LeptonID
       and entry.mll>20 and entry.Lepton_pt[0]>25 and entry.Lepton_pt[1]>20  and entry.Lepton_pdgId[0]*entry.Lepton_pdgId[1]<0 #OC
       and entry.CleanJet_pt[0]>150. and entry.CleanJet_pt[0]!=entry.leadingPtTagged and np.arccos(np.cos(entry.MET_phi-entry.CleanJet_phi[0]))>2.5 #ISR
       and entry.nCleanJet>1 ): #ncleanjets
        passedcut=True
#        print entry.MET_pt, passedcut, entry.Lepton_pdgId[0]*entry.Lepton_pdgId[1], entry.Lepton_pdgId[0],entry.Lepton_pdgId[1]
    return passedcut

#Function to cut separately
def flavour_tag(entry, btag):
    df=False
    sf=False
    sameflavour=None
    massZ = 91.1876
    vetoZ = abs(entry.mll-massZ)<15.
    weight= 1
    btagW  = entry.btagWeight_1tag
    bvetoW = 1-entry.btagWeight_1tag
    dfcut = abs(entry.Lepton_pdgId[0])!=abs(entry.Lepton_pdgId[1])
    #print "btag", abs(entry.Lepton_pdgId[0]), abs(entry.Lepton_pdgId[1]), abs(entry.mll-massZ)
    if dfcut is True:
        #print 'df'
        sameflavour=False
        if btag is True:
            weight=btagW
        else:
            weight=bvetoW
                
    else:
        if vetoZ is True:
            #print '--->sf'
            sameflavour=True
            if btag is True:
                weight=btagW
            else:
                weight=bvetoW
#    if True in [sf,df]: print "theres a flavour"
#    print df, sf, "---", sameflavour
    return sameflavour,weight




#print ttbar_evs.ls()#Draw("CleanJet_pt")
#print ttbar_evs.GetEntries()#, CleanJet_pt

def SVentries(entry,maxSV, SV_eta, SV_phi, SV_pt, SV_mass, SV_x, SV_y, SV_z, SV_chi2):
    lastSV=min(len(entry.SV_eta),maxSV)
    '''for iSV in range(0,lastSV):
        SV_eta[iSV]  = entry.SV_eta[iSV]
        SV_phi[iSV]  = entry.SV_phi[iSV]
        SV_pt[iSV]   = entry.SV_pt[iSV]
        SV_mass[iSV] = entry.SV_mass[iSV]
        SV_x[iSV]    = entry.SV_x[iSV]
        SV_y[iSV]    = entry.SV_y[iSV]
        SV_z[iSV]    = entry.SV_z[iSV]
        SV_chi2[iSV] = entry.SV_chi2[iSV]

        print "SV pt[isv]", SV_pt[iSV]
#    if(lastSV<maxSV):
#       for idx,entry in enumerate(entry.SV_eta):
#           if(entry=-999) del entry.SV_eta[idx]
    print len(entry.SV_eta),maxSV
    '''
    return SV_eta, SV_phi, SV_pt, SV_mass, SV_x, SV_y, SV_z, SV_chi2


def loopentries(sample,hdnjets_sf_veto, hdnjets_df_veto, hdnjets_sf_tag, hdnjets_df_tag,tree, samplenm, idxmax=1000):
    #print "TYYPEEES", type(rootfile), type(tree)
    #tree=ROOT.TTree(samplenm)

#    tree = ROOT.TTree(samplenm,"A ROOT tree");
    #print ROOT.hfile.ls()
    maxSV=20
    evid = array('i',[-1])
    nSV  = array('i',[-1])
    Dnjetstot = array('d',[-999])
    isSF   = array('i',[-999])
    btagW  = array('d',[-999])
    bvetoW = array('d',[-999])

    MET_sumEt = array('d',[-999])
    MET_pt    = array('d',[-999])

    PV_x    = array('f',[-999])
    PV_y    = array('f',[-999])
    PV_z    = array('f',[-999])
    PV_npvs = array('f',[-999])
    PV_chi2 = array('f',[-999])

    SV_eta  = array('f',maxSV*[-999])
    

    SV_phi  = array('d',maxSV*[-999])
    SV_pt   = array('d',maxSV*[-999])
    SV_mass = array('d',maxSV*[-999])
    SV_x    = array('d',maxSV*[-999])
    SV_y    = array('d',maxSV*[-999])
    SV_z    = array('d',maxSV*[-999])
    SV_chi2 = array('d',maxSV*[-999])
    
    mll     = array('d',[-999])
    mt2ll   = array('d',[-999])
    nLepton = array('i',[-999])
    ptmiss  = array('d',[-999])

    susyMstop = array('f',[-999])
    #susyMLSP  = array('f',[-999])
    susyMLSP  = array('f',[-999])

    
    tree.Branch ("evid_"     +samplenm,evid     , "evid/I");

    tree.Branch ("Dnjetstot_"+samplenm,Dnjetstot, "Dnjetstot/D");
    tree.Branch ("btagW_"    +samplenm,btagW    , "btagW/D");
    tree.Branch ("bvetoW_"   +samplenm,bvetoW   , "bvetoW/D");
    tree.Branch ("isSF_"     +samplenm,isSF     , "isSF/I");

    tree.Branch ("MET_sumEt_"+samplenm,MET_sumEt, "MET_sumEt/D");
    tree.Branch ("MET_pt_"   +samplenm,MET_pt   , "MET_pt/D");

    tree.Branch ("PV_x_"   +samplenm,PV_x   , "PV_x/F");
    tree.Branch ("PV_y_"   +samplenm,PV_y   , "PV_y/F");
    tree.Branch ("PV_z_"   +samplenm,PV_z   , "PV_z/F");
    tree.Branch ("PV_npvs_"+samplenm,PV_npvs, "PV_npvs/F");
    tree.Branch ("PV_chi2_"+samplenm,PV_chi2, "PV_chi2/F");

    nSV_b=tree.Branch ("nSV_" +samplenm,nSV , "nSV/I");

    
    tree.Branch ("SV_eta_" +samplenm,SV_eta , "SV_eta[nSV]/F");
    
    tree.Branch ("SV_phi_" +samplenm,SV_phi , "SV_phi[nSV]/D");
    
    tree.Branch ("SV_pt_"  +samplenm,SV_pt  , "SV_pt[nSV]/D");
    tree.Branch ("SV_mass_"+samplenm,SV_mass, "SV_mass[nSV]/D");
    tree.Branch ("SV_x_"   +samplenm,SV_x   , "SV_x[nSV]/D");
    tree.Branch ("SV_y_"   +samplenm,SV_y   , "SV_y[nSV]/D");
    tree.Branch ("SV_z_"   +samplenm,SV_z   , "SV_z[nSV]/D");
    tree.Branch ("SV_chi2_"+samplenm,SV_chi2, "SV_chi2[nSV]/D");
    
    tree.Branch ("mll_"       +samplenm,mll    , "mll/D");
    tree.Branch ("mt2ll_"     +samplenm,mt2ll  , "mt2ll/D");
    tree.Branch ("nLepton_"   +samplenm,nLepton, "nLepton/I");
    tree.Branch ("ptmiss_"    +samplenm,ptmiss , "ptmiss/D");

    tree.Branch ("susyMLSP__"   +samplenm,susyMLSP , "susyMLSP/F");
    tree.Branch ("susyMstop_"   +samplenm,susyMstop, "susyMstop/F");
    #tree.Branch ("susyMLSP_ "   +samplenm,susyMLSP , "susyMLSP/F");

    #tree.Branch ("hdnjets_df_tag", hdnjets_df_tag, "Dnjets2/F");
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
        passedcut=applycuts(entry) 
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
        
        #thiseta=array('d',[0])
        #if(len(entry.SV_eta)>(maxSV)): print "WARNING: more Entries in SV than in the vector to be stored in " 
        #        SV_eta, SV_phi, SV_pt, SV_mass, SV_x, SV_y, SV_z, SV_chi2=SVentries(entry,maxSV, SV_eta, SV_phi, SV_pt, SV_mass, SV_x, SV_y, SV_z, SV_chi2)
        lastSV=min(len(entry.SV_eta),maxSV)
        #print lastSV, len(entry.SV_eta), maxSV
        #print "i get to nsv"
        #if(len(entry.SV_eta)>0):
        nSV[0] = len(entry.SV_eta)
        #print nSV[0], counter
        
        for iSV in range(0,2):#lastSV):
            #continue
            print iSV, nSV[0]
            if(iSV<nSV[0]):
                print len(entry.SV_eta)
                print "ii", entry.SV_eta[iSV]
                SV_eta[iSV]  = entry.SV_eta[iSV]
                SV_phi[iSV]  = entry.SV_phi[iSV]
                SV_pt[iSV]   = entry.SV_pt[iSV]
	        SV_mass[iSV] = entry.SV_mass[iSV]
                SV_x[iSV]    = entry.SV_x[iSV]
	        SV_y[iSV]    = entry.SV_y[iSV]
	        SV_z[iSV]    = entry.SV_z[iSV]
	        SV_chi2[iSV] = entry.SV_chi2[iSV]
            
        




        #print SV_eta
        #nSV.Fill()
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

hdnjets_sf_veto_ttbar, hdnjets_df_veto_ttbar, hdnjets_sf_tag_ttbar, hdnjets_df_tag_ttbar= loopentries(ttbar_evs, hdnjets_sf_veto_ttbar, hdnjets_df_veto_ttbar, hdnjets_sf_tag_ttbar, hdnjets_df_tag_ttbar,ttbar_t, "ttbar",10000000)
#gc.collect()
#hfile.cd()
#tree.Write()

hdnjets_sf_veto_T2tt, hdnjets_df_veto_T2tt, hdnjets_sf_tag_T2tt, hdnjets_df_tag_T2tt= loopentries(signal_evs,hdnjets_sf_veto_T2tt, hdnjets_df_veto_T2tt, hdnjets_sf_tag_T2tt, hdnjets_df_tag_T2tt, signal_t, "signal", 1000000)

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
