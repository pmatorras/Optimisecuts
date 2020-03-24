from array import array
import ROOT
hfilenm="output/postcuts.root"
inputFilettbar  = '../rootfiles/nanoLatino_TTTo2L2Nu__part10.root'
inputFileWW     = '../rootfiles/nanoLatino_WWTo2L2Nu__part0.root'
inputFilesignal = '../rootfiles/nanoLatino_T2tt__mStop-400to1200__part1.root'

#nttbar=10000
#nT2tt =1000
def defBranches(tree,samplenm):
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

    nPhoton = array('i',[-999])
    visHTv3 = array('d',[-999])
    dphill  = array('d',[-999])
    lep1_pt = array('d',[-999])
    lep2_pt = array('d',[-999])
    dRll    = array('d',[-999])
    detall  = array('d',[-999])
    dphijj  = array('d',[-999])
    dRjj    = array('d',[-999])
    detajj  = array('d',[-999])
    jet1_pt = array('d',[-999])
    jet2_pt = array('d',[-999])
    #Define tree branches
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

    tree.Branch ("susyMLSP_"   +samplenm,susyMLSP , "susyMLSP/F");
    tree.Branch ("susyMstop_"   +samplenm,susyMstop, "susyMstop/F");
    
    tree.Branch ("nPhoton_" +samplenm,nPhoton , "nPhoton/I");
    tree.Branch ("visHTv3_" +samplenm,visHTv3 , "visHTv3/D");

    tree.Branch ("lep1_pt_"+samplenm,lep1_pt, "lep1_pt/D");
    tree.Branch ("lep2_pt_"+samplenm,lep2_pt, "lep2_pt/D");
    tree.Branch ("jet1_pt_"+samplenm,jet1_pt, "jet1_pt/D");
    tree.Branch ("jet2_pt_"+samplenm,jet2_pt, "jet_pt/D");

    tree.Branch ("dphill_"+samplenm,dphill , "dphill/D");
    tree.Branch ("detall_"+samplenm,detall, "detall/D");
    tree.Branch ("dRll_"  +samplenm,dRll  , "dRll/D");

    tree.Branch ("dphijj_"+samplenm,dphijj , "dphijj/D");
    tree.Branch ("detajj_"+samplenm,detajj, "detajj/D");
    tree.Branch ("dRjj_" +samplenm,dRjj  , "dRjj/D");
    
    
    return evid, maxSV, nSV, Dnjetstot,isSF,btagW, bvetoW,MET_sumEt,MET_pt,\
        PV_x,PV_y,PV_z,PV_npvs,PV_chi2,SV_eta,SV_phi,SV_pt,SV_mass,\
        SV_x,SV_y,SV_z,SV_chi2,mll,mt2ll,nLepton,ptmiss,susyMstop,susyMLSP,\
        lep1_pt,lep2_pt,jet1_pt,jet2_pt,dphill,detall,dRll,dphijj,detajj,dRjj
