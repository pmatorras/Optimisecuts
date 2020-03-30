#Apply default cuts
bkgs=['ttbar', 'WW']
def defaultcuts(entry,samplenm):
    if(samplenm in bkgs):
        pttag=bool(entry.CleanJet_pt[0]!=entry.leadingPtTagged_btagDeepBM_1c)
        
    else:
        pttag=bool(entry.CleanJet_pt[0]!=entry.leadingPtTagged)
    
    passedcut=False
    metcut=bool(entry.ptmiss >300)
    lepidcut=bool((entry.Lepton_isTightElectron_cutBasedMediumPOG[0] + entry.Lepton_isTightMuon_mediumRelIsoTight[0]+entry.Lepton_isTightElectron_cutBasedMediumPOG[1]+entry.Lepton_isTightMuon_mediumRelIsoTight[1])==2)
    occut=bool(entry.mll>20 and entry.Lepton_pt[0]>25 and entry.Lepton_pt[1]>20  and entry.Lepton_pdgId[0]*entry.Lepton_pdgId[1]<0)
    isrcut=bool(entry.CleanJet_pt[0]>150. and pttag and np.arccos(np.cos(entry.MET_phi-entry.CleanJet_phi[0]))>2.5)#REMEMBER TO CHANGE THIS
    #ncleanjets=bool(entry.nCleanJet>1)
    conds=[metcut,lepidcut,occut] #,ncleanjets, isrcut 
    if all(conds) is True :
        passedcut=True
    return passedcut, isrcut


#Function to cut separately                                                                                 
def flavour_tag(entry, samplenm):
    df=False
    sf=False
    sameflavour=-999
    massZ = 91.1876
    vetoZ = abs(entry.mll-massZ)<15.
    if(samplenm in bkgs):
        btagW = entry.btagWeight_1tag_btagDeepBM_1c
    else:
        btagW = entry.btagWeight_1tag
    bvetoW = 1-btagW
    dfcut = abs(entry.Lepton_pdgId[0])!=abs(entry.Lepton_pdgId[1])
    if dfcut is True:
        #print 'df'                         
        sameflavour=-1
    else:
        if vetoZ is True:
            sameflavour=1
            
    return sameflavour,btagW,bvetoW
def SVentries(entry,nSV, SV_eta, SV_phi, SV_pt, SV_mass, SV_x, SV_y, SV_z, SV_chi2):
    for iSV in range(0,nSV):#lastSV):            
        SV_eta[iSV]  = entry.SV_eta[iSV]
	SV_phi[iSV]  = entry.SV_phi[iSV]
	SV_pt[iSV]   = entry.SV_pt[iSV]
        SV_mass[iSV] = entry.SV_mass[iSV]
	SV_x[iSV]    = entry.SV_x[iSV]
	SV_y[iSV]    = entry.SV_y[iSV]
        SV_z[iSV]    = entry.SV_z[iSV]
        SV_chi2[iSV] = entry.SV_chi2[iSV]

    return
