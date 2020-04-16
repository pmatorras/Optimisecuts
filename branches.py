from array import array
import ROOT
cutfolder = "./Test/"
cutfilenm = cutfolder+"postcuts.root"
inputFilettbar = '/eos/cms/store/user/scodella/SUSY/Nano/Summer16_102X_nAODv6_Full2016v6/MCSusy2016v6__MCCorr2016Susyv6__susyMT2/nanoLatino_TTTo2L2Nu__part10.root' 
inputFileWW    = '/eos/cms/store/user/scodella/SUSY/Nano/Summer16_102X_nAODv6_Full2016v6/MCSusy2016v6__MCCorr2016Susyv6__susyMT2/nanoLatino_WWTo2L2Nu__part1.root'
inputFileZZ    = '/eos/cms/store/user/scodella/SUSY/Nano/Summer16_102X_nAODv6_Full2016v6/MCSusy2016v6__MCCorr2016Susyv6__susyMT2/nanoLatino_ZZTo2L2Nu_ext1__part0.root'
inputFilettZ   = '/eos/cms/store/user/scodella/SUSY/Nano/Summer16_102X_nAODv6_Full2016v6/MCSusy2016v6__MCCorr2016Susyv6__susyMT2/nanoLatino_TTZToLLNuNu_M-10_ext3__part0.root'

inputFileT2tt  = '/eos/user/s/scodella/SUSY/Nano/Summer16FS_102X_nAODv4_Full2016v4/susyGen__susyW__MCSusy2016FS__MCCorr2016SusyFS__susyMT2FS/nanoLatino_T2tt__mStop-400to1200__part1.root'

#nttbar=10000
#nT2tt =1000
def defBranches(tree,samplenm):
    evid = array('i',[-1])
    maxSV=20
    nSV  = array('i',[-1])
    njets =array('i',[-999])
    nCleanjets=array('i',[-999])
    nbjets =array('i',[-999])
    nbCleanjets=array('i',[-999])
    Dnjetstot = array('i',[-999])
    Dnbjetstot =array('i',[-999])
    
    isSF   = array('i',[-999])
    btagW  = array('d',[-999])
    bvetoW = array('d',[-999])

    MET_sumEt = array('d',[-999])
    MET_pt    = array('d',[-999])
    MET_phi   = array('d',[-999])
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
    ISRcut    = array('i',[-999])
    
    nPhoton = array('i',[-999])
    visHTv3 = array('d',[-999])
    dphill  = array('d',[-999])
    
    lep1_pt = array('d',[-999])
    lep2_pt = array('d',[-999])
    lep1_eta = array('d',[-999])
    lep2_eta = array('d',[-999])
    lep1_phi = array('d',[-999])
    lep2_phi = array('d',[-999])

    jet1_pt = array('d',[-999])
    jet2_pt = array('d',[-999])
    jet1_eta = array('d',[-999])
    jet2_eta = array('d',[-999])
    jet1_phi = array('d',[-999])
    jet2_phi = array('d',[-999])

    bjet1_pt = array('d',[-999])
    bjet2_pt = array('d',[-999])
    bjet1_eta = array('d',[-999])
    bjet2_eta = array('d',[-999])
    bjet1_phi = array('d',[-999])
    bjet2_phi = array('d',[-999])

    dphil1MET = array('d',[-999])
    dphil1jmin= array('d',[-999])
    dphil1jmax= array('d',[-999])
    dphil1bmin= array('d',[-999])
    dphil1bmax= array('d',[-999])

    dphil2MET = array('d',[-999])
    dphil2jmin= array('d',[-999])
    dphil2jmax= array('d',[-999])
    dphil2bmin= array('d',[-999])
    dphil2bmax= array('d',[-999])

    dphilepsumMET = array('d',[-999])
    ptlepsum      = array('d',[-999])
    dRll    = array('d',[-999])
    detall  = array('d',[-999])
    dphijj  = array('d',[-999])
    
    dRjj    = array('d',[-999])
    detajj  = array('d',[-999])






    #Define tree branches
    tree.Branch ("evid_"     +samplenm,evid     , "evid/I");
    
    tree.Branch ("Dnjetstot_"+samplenm,Dnjetstot, "Dnjetstot/i");
    tree.Branch ("njets_"+samplenm,njets, "njets/i");
    tree.Branch ("nbjets_"+samplenm,nbjets, "nbjets/i");
    tree.Branch ("nCleanjets_"+samplenm,nCleanjets, "nCleanjets/i");
    tree.Branch ("nbCleanjets_"+samplenm,nbCleanjets, "nbCleanjets/i");
    tree.Branch ("Dnbjetstot_"+samplenm,Dnbjetstot, "Dnbjetstot/i");
    
    tree.Branch ("btagW_"    +samplenm,btagW    , "btagW/D");
    tree.Branch ("bvetoW_"   +samplenm,bvetoW   , "bvetoW/D");
    tree.Branch ("isSF_"     +samplenm,isSF     , "isSF/I");

    tree.Branch ("MET_sumEt_"+samplenm,MET_sumEt, "MET_sumEt/D");
    tree.Branch ("MET_pt_"   +samplenm,MET_pt   , "MET_pt/D");
    tree.Branch ("MET_phi_"  +samplenm,MET_phi  , "MET_phi/D");

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
    tree.Branch ("ISRcut_"   +samplenm,ISRcut, "ISRcut/I");
    
    tree.Branch ("nPhoton_" +samplenm,nPhoton , "nPhoton/I");
    tree.Branch ("visHTv3_" +samplenm,visHTv3 , "visHTv3/D");

    tree.Branch ("lep1_pt_"+samplenm,lep1_pt, "lep1_pt/D");
    tree.Branch ("lep2_pt_"+samplenm,lep2_pt, "lep2_pt/D");
    tree.Branch ("lep1_phi_"+samplenm,lep1_phi, "lep1_phi/D");
    tree.Branch ("lep2_phi_"+samplenm,lep2_phi, "lep2_phi/D");
    tree.Branch ("lep1_eta_"+samplenm,lep1_eta, "lep1_eta/D");
    tree.Branch ("lep2_eta_"+samplenm,lep2_eta, "lep2_eta/D");

    tree.Branch ("jet1_pt_"+samplenm,jet1_pt, "jet1_pt/D");
    tree.Branch ("jet2_pt_"+samplenm,jet2_pt, "jet2_pt/D");
    tree.Branch ("jet1_phi_"+samplenm,jet1_phi, "jet1_phi/D");
    tree.Branch ("jet2_phi_"+samplenm,jet2_phi, "jet2_phi/D");
    tree.Branch ("jet1_eta_"+samplenm,jet1_eta, "jet1_eta/D");
    tree.Branch ("jet2_eta_"+samplenm,jet2_eta, "jet2_eta/D");

    tree.Branch ("bjet1_pt_"+samplenm,bjet1_pt, "bjet1_pt/D");
    tree.Branch ("bjet2_pt_"+samplenm,bjet2_pt, "bjet2_pt/D");
    tree.Branch ("bjet1_phi_"+samplenm,bjet1_phi, "bjet1_phi/D");
    tree.Branch ("bjet2_phi_"+samplenm,bjet2_phi, "bjet2_phi/D");
    tree.Branch ("bjet1_eta_"+samplenm,bjet1_eta, "bjet1_eta/D");
    tree.Branch ("bjet2_eta_"+samplenm,bjet2_eta, "bjet2_eta/D");

    tree.Branch ("dphill_"+samplenm,dphill , "dphill/D");
    tree.Branch ("detall_"+samplenm,detall, "detall/D");
    tree.Branch ("dRll_"  +samplenm,dRll  , "dRll/D");

    tree.Branch ("dphijj_"+samplenm,dphijj , "dphijj/D");
    tree.Branch ("detajj_"+samplenm,detajj, "detajj/D");
    tree.Branch ("dRjj_" +samplenm,dRjj  , "dRjj/D");

    tree.Branch ("dphil1MET_" +samplenm,dphil1MET , "dphil1MET/D");
    tree.Branch ("dphil1jmin_"+samplenm,dphil1jmin , "dphil1jmin/D");
    tree.Branch ("dphil1jmax_"+samplenm,dphil1jmax , "dphil1jmax/D");
    tree.Branch ("dphil1bmin_"+samplenm,dphil1bmin , "dphil1bmin/D");
    tree.Branch ("dphil1bmax_"+samplenm,dphil1bmax , "dphil1bmax/D");

    tree.Branch ("dphil2MET_" +samplenm,dphil2MET  , "dphil2MET/D");
    tree.Branch ("dphil2jmin_"+samplenm,dphil2jmin , "dphil2jmin/D");
    tree.Branch ("dphil2jmax_"+samplenm,dphil2jmax , "dphil2jmax/D");
    tree.Branch ("dphil2bmin_"+samplenm,dphil2bmin , "dphil2bmin/D");
    tree.Branch ("dphil2bmax_"+samplenm,dphil2bmax , "dphil2bmax/D");

    tree.Branch ("dphilepsumMET_"+samplenm,dphilepsumMET , "dphilepsumMET/D");
    tree.Branch ("ptlepsum_"     +samplenm,ptlepsum      , "ptlepsum/D");

    return evid, maxSV, nSV, njets, nbjets,nCleanjets, nbCleanjets, Dnjetstot, Dnbjetstot, isSF,btagW, bvetoW,\
        PV_x,PV_y,PV_z,PV_npvs,PV_chi2,SV_eta,SV_phi,SV_pt,SV_mass,SV_x,SV_y,SV_z,SV_chi2,\
        MET_sumEt,MET_pt,MET_phi,mll,mt2ll,nLepton,ptmiss,susyMstop,susyMLSP,ISRcut,detall,dRll,detajj,dRjj,\
        dphill,dphijj, dphil1jmin,dphil1jmax,dphil1bmin,dphil1bmax, dphil2jmin,dphil2jmax,dphil2bmin,dphil2bmax,\
        lep1_pt, lep1_eta,lep1_phi,lep2_pt,lep2_eta,lep2_phi, jet1_pt,jet1_eta,jet1_phi,jet2_pt,jet2_eta,jet2_phi,\
        bjet1_pt,bjet1_eta,bjet1_phi,bjet2_pt,bjet2_eta,bjet2_phi, dphil1MET,dphil2MET, dphilepsumMET, ptlepsum
