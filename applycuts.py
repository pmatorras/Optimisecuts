#Apply default cuts
def defaultcuts(entry):
    passedcut=False
    metcut=bool(entry.MET_pt >300)
    lepidcut=bool((entry.Lepton_isTightElectron_cutBasedMediumPOG[0] + entry.Lepton_isTightMuon_mediumRelIsoTight[0]+entry.Lepton_isTightElectron_cutBasedMediumPOG[1]+entry.Lepton_isTightMuon_mediumRelIsoTight[1])==2)
    occut=bool(entry.mll>20 and entry.Lepton_pt[0]>25 and entry.Lepton_pt[1]>20  and entry.Lepton_pdgId[0]*entry.Lepton_pdgId[1]<0)
    isrcut=bool(entry.CleanJet_pt[0]>150. and entry.CleanJet_pt[0]!=entry.leadingPtTagged and np.arccos(np.cos(entry.MET_phi-entry.CleanJet_phi[0]))>2.5)
    ncleanjets=bool(entry.nCleanJet>1)
    conds=[metcut,lepidcut,occut,isrcut,ncleanjets]
    if all(conds) is True :
        passedcut=True
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
