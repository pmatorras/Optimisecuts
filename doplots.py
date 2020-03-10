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

hfile = ROOT.TFile("postcuts.root","RECREATE","Example");
tree  = ROOT.TTree("tree", "for both")
signal = ROOT.TTree("signal","A ROOT tree");
bkgs   = ROOT.TTree("bkgs","A ROOT tree");




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


def loopentries(sample,hdnjets_sf_veto, hdnjets_df_veto, hdnjets_sf_tag, hdnjets_df_tag, samplenm, idxmax=1000):
    #print "TYYPEEES", type(rootfile), type(tree)
    #tree=ROOT.TTree(samplenm)

#    tree = ROOT.TTree(samplenm,"A ROOT tree");
    #print ROOT.hfile.ls()
    Dnjetstot = array('d',[0])
    isSF   = array('i',[0])
    btagW  = array('d',[0])
    bvetoW = array('d',[0])
    
    tree.Branch ("Dnjetstot_"+samplenm,Dnjetstot, "Dnjetstot/D");
    tree.Branch ("btagW_" +samplenm,btagW , "btagW/D");
    tree.Branch ("bvetoW_"+samplenm,bvetoW, "bvetoW/D");

    tree.Branch ("isSF_"+samplenm,isSF, "isSF/I");

    #tree.Branch ("hdnjets_df_tag", hdnjets_df_tag, "Dnjets2/F");

    print "looping over", sample, "\nfor",idxmax, "events"
    for idx, entry in enumerate(sample):
        if(idx>idxmax): break
        #apply normal cuts
        if entry.nCleanJet<1: continue
        if(idx==7532): print "####", idx, entry.CleanJet_pt[0], #entry.leadingPtTagged, entry.MET_phi , entry.CleanJet_phi[0]
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
        #print "dnjets", Dnjetstot
        tree.Fill()

        #print "output", sameflavour, bweight,"\n------------>IDX:", idx

        #print "ENTRIES--------\n",tree.Scan("Dnjetstot")
    print "fill tfile"
 
    print "#############################\n"#,tree.Scan("Dnjetstot")
    #hfile.Close()

    return hdnjets_sf_veto, hdnjets_df_veto, hdnjets_sf_tag, hdnjets_df_tag
#print (eval(ttbar_evs))

#tree.Scan("Dnjetstot")
ttbar = ROOT.TFile.Open(inputFilettbar ,"READ")
ttbar_evs= ttbar.Get('Events')
signal = ROOT.TFile.Open(inputFilesignal ,"READ")
signal_evs= signal.Get('Events')

#hdnjets_sf_veto_ttbar, hdnjets_df_veto_ttbar, hdnjets_sf_tag_ttbar, hdnjets_df_tag_ttbar= loopentries(ttbar_evs, hdnjets_sf_veto_ttbar, hdnjets_df_veto_ttbar, hdnjets_sf_tag_ttbar, hdnjets_df_tag_ttbar, "ttbar",10000)
hdnjets_sf_veto_T2tt, hdnjets_df_veto_T2tt, hdnjets_sf_tag_T2tt, hdnjets_df_tag_T2tt= loopentries(signal_evs,hdnjets_sf_veto_T2tt, hdnjets_df_veto_T2tt, hdnjets_sf_tag_T2tt, hdnjets_df_tag_T2tt, "signal")
#tree.Scan("Dnjetstot")
print "here?"
hfile.cd()
tree.Write()

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



'''
hdnjets_sf_tag_ttbar.Draw('hist')
hdnjets_sf_tag_T2tt.SetMarkerColor(2)
hdnjets_sf_tag_T2tt.SetLineColor(2)
hdnjets_sf_tag_T2tt.Draw('same')
c1.SaveAs("hdnjets_sf_tag.png")

hdnjets_df_tag_ttbar.Draw('hist')
hdnjets_df_tag_T2tt.SetMarkerColor(2)
hdnjets_df_tag_T2tt.SetLineColor(2)
hdnjets_df_tag_T2tt.Draw('same')
c1.SaveAs("hdnjets_df_tag.png")
'''
print type(ttbar_evs)
'''
hdnjets_df_tag.Draw('hist')
c1.SaveAs("hdnjets_df_tag.png")
hdnjets_sf_veto.Draw('hist')
c1.SaveAs("hdnjets_sf_veto.png")
hdnjets_sf_tag.Draw('hist')
c1.SaveAs("hdnjets_df_veto.png")
'''

#os.system('display hdnjets_sf_tag.png')





#c1.SetLogy()                                                                                                                        
#c1.Modified()
#c1.Update()
#c1.SaveAs("try.png")
#os.system('display try.png')



'''
for idx,entry  in enumerate(ttbar_evs):
    if(idx>10000):
        break
    #cuts for sr3
    metcut=False
    applycuts(entry,metcut)
    print metcut
    if(entry.MET_pt <300): continue
    if(entry.mll<20 or entry.Lepton_pt[0]<25 or entry.Lepton_pt[1]<20 or entry.Lepton_pdgId[0]*entry.Lepton_pdgId[1]>0): continue
    #cuts for isr
    CleanJet_pt=entry.CleanJet_pt
    ncleanJet=len(CleanJet_pt)
    if(ncleanJet<1): continue
    if(entry.CleanJet_pt[0]<150 and entry.CleanJet_pt[0]==entry.leadingPtTagged): continue# and np.arccos(np.cos(entry.MET_phi-entry.CleanJet_phi[0]))<2.5]):continue
#    if(any(isr)): continue
    
    #print idx
    if(abs(entry.Lepton_pdgId[0])==abs(entry.Lepton_pdgId[1])):
        print idx,"sf", entry.nGenJet, entry.nJet, ncleanJet
        hdnjets_sf_tag.Fill(entry.nJet-ncleanJet, entry.btagWeight_1tag)
    if(abs(entry.Lepton_pdgId[0])!=abs(entry.Lepton_pdgId[1])):
        continue
        print "df"
        #print idx, "after"

    njet=entry.nJet
    #print ncleanJet-njet
'''
