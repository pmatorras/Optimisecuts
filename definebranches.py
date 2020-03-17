from array import array
import ROOT
hfilenm="output/postcuts.root"
inputFilettbar  = '../rootfiles/nanoLatino_TTTo2L2Nu__part10.root'
inputFilesignal = '../rootfiles/nanoLatino_T2tt__mStop-400to1200__part0.root'
nttbar=100000000
nT2tt =10000000
def defVectors():
    evid = array('i',[-1])
    maxSV=20
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
    susyMLSP  = array('f',[-999])
    
    return evid, maxSV, nSV, Dnjetstot,isSF,btagW, bvetoW,MET_sumEt,MET_pt,\
        PV_x,PV_y,PV_z,PV_npvs,PV_chi2,SV_eta,SV_phi,SV_pt,SV_mass,\
        SV_x,SV_y,SV_z,SV_chi2, mll,mt2ll,nLepton,ptmiss,susyMstop,susyMLSP

def defBranches(tree,samplenm):
    evid, maxSV,nSV, Dnjetstot,isSF,btagW, bvetoW,MET_sumEt,MET_pt,\
        PV_x,PV_y,PV_z,PV_npvs,PV_chi2,SV_eta,SV_phi,SV_pt,SV_mass,\
        SV_x,SV_y,SV_z,SV_chi2,mll,mt2ll,nLepton,ptmiss,susyMstop,susyMLSP= defVectors()
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

    tree.Branch ("nSV_" +samplenm,nSV , "nSV/I");

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
    
    
    return evid, maxSV, nSV, Dnjetstot,isSF,btagW, bvetoW,MET_sumEt,MET_pt,\
        PV_x,PV_y,PV_z,PV_npvs,PV_chi2,SV_eta,SV_phi,SV_pt,SV_mass,\
        SV_x,SV_y,SV_z,SV_chi2,mll,mt2ll,nLepton,ptmiss,susyMstop,susyMLSP
