#include "BSM_Analysis.h"
#include "interface/selection.h"

int main (int argc, char *argv[])
{

  TApplication app("App",&argc, argv);
  TFile * MassHisto = new TFile("MassHisto.root", "RECREATE");  
  int nDir = 10;
  TDirectory *theDirectory[nDir];
  theDirectory[0] = MassHisto->mkdir("AfterMuonChargeProduct");
  theDirectory[1] = MassHisto->mkdir("AfterTau_decayModeFinding");
  theDirectory[2] = MassHisto->mkdir("AfterTau_IsoLoosDeltaBetaCorr3Hits");
  theDirectory[3] = MassHisto->mkdir("AfterTau_IsoMdiumDeltaBetaCorr3Hits");
  theDirectory[4] = MassHisto->mkdir("AfterTau_IsoTightDeltaBetaCorr3Hits");
  theDirectory[5] = MassHisto->mkdir("AfterTau_AgainstMuonLoose2");
  theDirectory[6] = MassHisto->mkdir("AfterTau_AgainstMuonLoose3");
  theDirectory[7] = MassHisto->mkdir("AfterTau_AgainstMuonTight2");
  theDirectory[8] = MassHisto->mkdir("AfterTau_AgainstMuonTight3");
  theDirectory[9] = MassHisto->mkdir("AfterTau_byTightIsolationMVA3newDMwLT");
  BSM_Analysis BSM_Analysis_(MassHisto, theDirectory, nDir);
}

BSM_Analysis::BSM_Analysis(TFile* theFile, TDirectory *cdDir[], int nDir) 
{
  crateHistoMasps(nDir);
  //load PU weights
  TFile file_PUdata("PUdata.root","read");
  TH1F *PUweights = (TH1F*)file_PUdata.Get("analyzeHiMassTau/NVertices_0");
  PUweights->Scale(1/PUweights->Integral());
  TFile file_PUsim("PUsim.root","read");
  TH1F *PUsim = (TH1F*)file_PUsim.Get("analyzeHiMassTau/NVertices_0");
  PUsim->Scale(1/PUsim->Integral());

  PUweights->Divide(PUsim);

  //configure input file
  //TFile *f = new TFile ("/eos/uscms/store/user/florez/Ntuples_Wjets/OutTree_9_1_p3o.root");
  TFile *f = new TFile ("OutTree.root");
  f->cd("TNT");
  TTree* BOOM = (TTree*)f->Get("TNT/BOOM");

  int nentries = (int) BOOM->GetEntries();
  setBranchAddress(BOOM);

  //for (int i = 0; i < 10000; ++i)
  for (int i = 0; i < nentries; ++i)
    {
      BOOM->GetEntry(i);

       //define global event weight
       double weight =1.;
       weight=PUweights->GetBinContent(PUweights->FindBin(trueInteractions));
       // TLorentz vector to calculate the di-jet invariant mass
       TLorentzVector TagMuon_TL_vec(0., 0., 0., 0.); 
       TLorentzVector ProbeTagMuon_TL_vec(0., 0., 0., 0.);
       double charge_lead = 0.;
       double charge_slead = 0.;
       int lmuon_counter = 0;
       int smuon_counter = 0;
       int pass_tau_id[nDir] = {0};

       // For Trigger
       for (int t = 0 ; t < 10; t++){
         cout <<Trigger_names->at(t)<<endl;
       }       

       for (int j = 0; j < Muon_pt->size(); j++)
         {
            // select Z -> mumu events
            double lead_muon_pt = Muon_pt->at(0);
            if ( (Muon_pt->size() == 2) && (abs(Muon_eta->at(j)) < 2.4) && (Muon_pt->at(j) > 10.0) ){
              if ((Muon_tight->at(j) == 1) && (Muon_isoSum->at(j) < 5.0)){
                  TagMuon_TL_vec.SetPtEtaPhiE(Muon_pt->at(j), Muon_eta->at(j), Muon_phi->at(j), Muon_energy->at(j)); 
                  charge_lead = Muon_charge->at(j);
                  lmuon_counter++;        
              } else {
                  ProbeTagMuon_TL_vec.SetPtEtaPhiE(Muon_pt->at(j), Muon_eta->at(j), Muon_phi->at(j), Muon_energy->at(j));
                  charge_slead = Muon_charge->at(j);
                  smuon_counter++;
                  for (int t = 0; t < Tau_pt->size(); t++)
                    {
                      if( Tau_decayModeFinding->at(t) == 1){pass_tau_id[1] = 1;}
                      if((Tau_decayModeFinding->at(t) == 1) && 
                         (Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(t) == 1)){pass_tau_id[2] = 1;}
                      if((Tau_decayModeFinding->at(t) == 1) &&
                         (Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(t) == 1)){pass_tau_id[3] = 1;}
                      if((Tau_decayModeFinding->at(t) == 1) &&
                         (Tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(t) == 1)){pass_tau_id[4] = 1;}
                      if((Tau_decayModeFinding->at(t) == 1) && (Tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(t) == 1) &&
                         (Tau_againstMuonLoose2->at(t) == 1)){pass_tau_id[5] = 1;} 
                      if((Tau_decayModeFinding->at(t) == 1) && (Tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(t) == 1) &&
                         (Tau_againstMuonLoose3->at(t) == 1)){pass_tau_id[6] = 1;}
                      if((Tau_decayModeFinding->at(t) == 1) && (Tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(t) == 1) &&
                         (Tau_againstMuonTight2->at(t) == 1)){pass_tau_id[7] = 1;}  
                      if((Tau_decayModeFinding->at(t) == 1) && (Tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(t) == 1) &&
                         (Tau_againstMuonTight3->at(t) == 1)){pass_tau_id[8] = 1;}
                      if((Tau_decayModeFinding->at(t) == 1) && (Tau_byTightIsolationMVA3newDMwLT->at(t) == 1) &&
                         (Tau_againstMuonTight3->at(t) == 1)){pass_tau_id[9] = 1;}
                   }
              }               
                
            }
         }

        double charge_product = charge_lead*charge_slead;
        bool passMuonPt = false;
        if ((ProbeTagMuon_TL_vec.Pt() > 25.) || (TagMuon_TL_vec.Pt() > 25.)){passMuonPt = true;}
        if ((charge_product <= -1) && passMuonPt){
          _hmap_diMuon_mass[0]->Fill((TagMuon_TL_vec + ProbeTagMuon_TL_vec).M());
          _hmap_muon_pT[0]->Fill(ProbeTagMuon_TL_vec.Pt());
          _hmap_muon_eta[0]->Fill(ProbeTagMuon_TL_vec.Eta());
        }
        for (int i = 1; i < nDir; i++){
          if ((charge_product == -1) && (passMuonPt) && (pass_tau_id[i] == 1)){
            _hmap_diMuon_mass[i]->Fill((TagMuon_TL_vec + ProbeTagMuon_TL_vec).M());
            _hmap_muon_pT[i]->Fill(ProbeTagMuon_TL_vec.Pt());
            _hmap_muon_eta[i]->Fill(ProbeTagMuon_TL_vec.Eta());
          }
        }
     }
   theFile->cd();
   for (int d = 0; d < nDir; d++)
     {
       cdDir[d]->cd();
       _hmap_diMuon_mass[d]->Write();
       _hmap_muon_pT[d]->Write();
       _hmap_muon_eta[d]->Write();
     } 
   theFile->Close();
}

BSM_Analysis::~BSM_Analysis()
{
  // do anything here that needs to be done at desctruction time
}

void BSM_Analysis::crateHistoMasps (int directories)
{
   for (int i = 0; i < directories; i++)
     {
       _hmap_diMuon_mass[i] = new TH1F("diMuonMass", "m_{#mu, #mu}", 300., 0., 300.); 
       _hmap_muon_pT[i]     = new TH1F("muon_pT", "#mu p_{T}", 300, 0., 300.);
       _hmap_muon_eta[i]    = new TH1F("muon_eta", "#mu #eta", 50, -2.5, 2.5);
     }
}

void BSM_Analysis::setBranchAddress(TTree* BOOM)
{

   // Set object pointer
   
   Muon_pt = 0;
   Muon_eta = 0;
   Muon_phi = 0;
   Muon_p = 0;
   Muon_energy = 0;
   Muon_charge = 0;
   Muon_tight = 0;
   Muon_soft = 0;
   Muon_pf = 0;
   Muon_isoCharged = 0;
   Muon_isoSum = 0;
   Muon_isoCharParPt = 0;
   Muon_chi2 = 0;
   Muon_validHits = 0;
   Muon_validHitsInner = 0;
   Muon_matchedStat = 0;
   Muon_dxy = 0;
   Muon_TLayers = 0;
   Muon_dz = 0;
   Muon_isoNeutralHadron = 0;
   Muon_isoPhoton = 0;
   Muon_isoPU = 0;
   Tau_eta = 0;
   Tau_phi = 0;
   Tau_pt = 0;
   Tau_energy = 0;
   Tau_charge = 0;
   Tau_decayModeFinding = 0;
   Tau_decayModeFindingNewDMs = 0;
   Tau_chargedIsoPtSum = 0;
   Tau_neutralIsoPtSum = 0;
   Tau_againstMuonTight3 = 0;
   Tau_againstElectronMVATightMVA5 = 0;
   Tau_nProngs = 0;
   Tau_puCorrPtSum = 0;
   Tau_byLooseCombinedIsolationDeltaBetaCorr = 0;
   Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
   Tau_byMediumCombinedIsolationDeltaBetaCorr = 0;
   Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
   Tau_byTightCombinedIsolationDeltaBetaCorr = 0;
   Tau_byTightCombinedIsolationDeltaBetaCorr3Hits = 0;
   Tau_byLooseIsolationMVA3newDMwLT = 0;
   Tau_byLooseIsolationMVA3newDMwoLT = 0;
   Tau_byLooseIsolationMva3oldDMwLT = 0;
   Tau_byLooseIsolationMVA3oldDMwoLT = 0;
   Tau_byMediumIsolationMVA3newDMwLT = 0;
   Tau_byMediumIsolationMVA3newDMwoLT = 0;
   Tau_byMediumIsolationMva3oldDMwLT = 0;
   Tau_byMediumIsolationMVA3oldDMwoLT = 0;
   Tau_byTightIsolationMVA3newDMwLT = 0;
   Tau_byTightIsolationMVA3newDMwoLT = 0;
   Tau_byTightIsolationMva3oldDMwLT = 0;
   Tau_byTightIsolationMVA3oldDMwoLT = 0;
   Tau_againstMuonLoose2 = 0;
   Tau_againstMuonLoose3 = 0;
   Tau_againstMuonTight2 = 0;
   Tau_againstElectronMVALooseMVA5 = 0;
   Tau_againstElectronMVAMediumMVA5 = 0;
   Tau_byVLooseCombinedIsolationDeltaBetaCorr = 0;
   Tau_leadChargedCandPt = 0;
   Tau_leadChargedCandCharge = 0;
   Tau_leadChargedCandEta = 0;
   Tau_leadChargedCandPhi = 0;
   Jet_pt = 0;
   Jet_eta = 0;
   Jet_phi = 0;
   Jet_energy = 0;
   Jet_bDiscriminator = 0;
   Jet_mass = 0;
   Jet_neutralHadEnergy = 0;
   Jet_neutralEmEmEnergy = 0;
   Jet_chargedHadronEnergy = 0;
   Jet_chargedEmEnergy = 0;
   Jet_muonEnergy = 0;
   Jet_electronEnergy = 0;
   Jet_photonEnergy = 0;
   UncorrJet_pt = 0;
   Trigger_names = 0;
   // Set branch addresses and branch pointers
   if(!BOOM) return;
   BOOM->SetBranchAddress("Trigger_names", &Trigger_names, &b_Trigger_names);
   BOOM->SetBranchAddress("Muon_pt", &Muon_pt, &b_Muon_pt);
   BOOM->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
   BOOM->SetBranchAddress("Muon_phi", &Muon_phi, &b_Muon_phi);
   BOOM->SetBranchAddress("Muon_p", &Muon_p, &b_Muon_p);
   BOOM->SetBranchAddress("Muon_energy", &Muon_energy, &b_Muon_energy);
   BOOM->SetBranchAddress("Muon_charge", &Muon_charge, &b_Muon_charge);
   BOOM->SetBranchAddress("Muon_tight", &Muon_tight, &b_Muon_tight);
   BOOM->SetBranchAddress("Muon_soft", &Muon_soft, &b_Muon_soft);
   BOOM->SetBranchAddress("Muon_pf", &Muon_pf, &b_Muon_pf);
   BOOM->SetBranchAddress("Muon_isoCharged", &Muon_isoCharged, &b_Muon_isoCharged);
   BOOM->SetBranchAddress("Muon_isoSum", &Muon_isoSum, &b_Muon_isoSum);
   BOOM->SetBranchAddress("Muon_isoCharParPt", &Muon_isoCharParPt, &b_Muon_isoCharParPt);
   BOOM->SetBranchAddress("Muon_chi2", &Muon_chi2, &b_Muon_chi2);
   BOOM->SetBranchAddress("Muon_validHits", &Muon_validHits, &b_Muon_validHits);
   BOOM->SetBranchAddress("Muon_validHitsInner", &Muon_validHitsInner, &b_Muon_validHitsInner);
   BOOM->SetBranchAddress("Muon_matchedStat", &Muon_matchedStat, &b_Muon_matchedStat);
   BOOM->SetBranchAddress("Muon_dxy", &Muon_dxy, &b_Muon_dxy);
   BOOM->SetBranchAddress("Muon_TLayers", &Muon_TLayers, &b_Muon_TLayers);
   BOOM->SetBranchAddress("Muon_dz", &Muon_dz, &b_Muon_dz);
   BOOM->SetBranchAddress("Muon_isoNeutralHadron", &Muon_isoNeutralHadron, &b_Muon_isoNeutralHadron);
   BOOM->SetBranchAddress("Muon_isoPhoton", &Muon_isoPhoton, &b_Muon_isoPhoton);
   BOOM->SetBranchAddress("Muon_isoPU", &Muon_isoPU, &b_Muon_isoPU);
   BOOM->SetBranchAddress("Tau_eta", &Tau_eta, &b_Tau_eta);
   BOOM->SetBranchAddress("Tau_phi", &Tau_phi, &b_Tau_phi);
   BOOM->SetBranchAddress("Tau_pt", &Tau_pt, &b_Tau_pt);
   BOOM->SetBranchAddress("Tau_energy", &Tau_energy, &b_Tau_energy);
   BOOM->SetBranchAddress("Tau_charge", &Tau_charge, &b_Tau_charge);
   BOOM->SetBranchAddress("Tau_decayModeFinding", &Tau_decayModeFinding, &b_Tau_decayModeFinding);
   BOOM->SetBranchAddress("Tau_decayModeFindingNewDMs", &Tau_decayModeFindingNewDMs, &b_Tau_decayModeFindingNewDMs);
   BOOM->SetBranchAddress("Tau_chargedIsoPtSum", &Tau_chargedIsoPtSum, &b_Tau_chargedIsoPtSum);
   BOOM->SetBranchAddress("Tau_neutralIsoPtSum", &Tau_neutralIsoPtSum, &b_Tau_neutralIsoPtSum);
   BOOM->SetBranchAddress("Tau_againstMuonTight3", &Tau_againstMuonTight3, &b_Tau_againstMuonTight3);
   BOOM->SetBranchAddress("Tau_againstElectronMVATightMVA5", &Tau_againstElectronMVATightMVA5, &b_Tau_againstElectronMVATightMVA5);
   BOOM->SetBranchAddress("Tau_nProngs", &Tau_nProngs, &b_Tau_nProngs);
   BOOM->SetBranchAddress("Tau_puCorrPtSum", &Tau_puCorrPtSum, &b_Tau_puCorrPtSum);
   BOOM->SetBranchAddress("Tau_byLooseCombinedIsolationDeltaBetaCorr", &Tau_byLooseCombinedIsolationDeltaBetaCorr, &b_Tau_byLooseCombinedIsolationDeltaBetaCorr);
   BOOM->SetBranchAddress("Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", &Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
   BOOM->SetBranchAddress("Tau_byMediumCombinedIsolationDeltaBetaCorr", &Tau_byMediumCombinedIsolationDeltaBetaCorr, &b_Tau_byMediumCombinedIsolationDeltaBetaCorr);
   BOOM->SetBranchAddress("Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
   BOOM->SetBranchAddress("Tau_byTightCombinedIsolationDeltaBetaCorr", &Tau_byTightCombinedIsolationDeltaBetaCorr, &b_Tau_byTightCombinedIsolationDeltaBetaCorr);
   BOOM->SetBranchAddress("Tau_byTightCombinedIsolationDeltaBetaCorr3Hits", &Tau_byTightCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
   BOOM->SetBranchAddress("Tau_byLooseIsolationMVA3newDMwLT", &Tau_byLooseIsolationMVA3newDMwLT, &b_Tau_byLooseIsolationMVA3newDMwLT);
   BOOM->SetBranchAddress("Tau_byLooseIsolationMVA3newDMwoLT", &Tau_byLooseIsolationMVA3newDMwoLT, &b_Tau_byLooseIsolationMVA3newDMwoLT);
   BOOM->SetBranchAddress("Tau_byLooseIsolationMva3oldDMwLT", &Tau_byLooseIsolationMva3oldDMwLT, &b_Tau_byLooseIsolationMva3oldDMwLT);
   BOOM->SetBranchAddress("Tau_byLooseIsolationMVA3oldDMwoLT", &Tau_byLooseIsolationMVA3oldDMwoLT, &b_Tau_byLooseIsolationMVA3oldDMwoLT);
   BOOM->SetBranchAddress("Tau_byMediumIsolationMVA3newDMwLT", &Tau_byMediumIsolationMVA3newDMwLT, &b_Tau_byMediumIsolationMVA3newDMwLT);
   BOOM->SetBranchAddress("Tau_byMediumIsolationMVA3newDMwoLT", &Tau_byMediumIsolationMVA3newDMwoLT, &b_Tau_byMediumIsolationMVA3newDMwoLT);
   BOOM->SetBranchAddress("Tau_byMediumIsolationMva3oldDMwLT", &Tau_byMediumIsolationMva3oldDMwLT, &b_Tau_byMediumIsolationMva3oldDMwLT);
   BOOM->SetBranchAddress("Tau_byMediumIsolationMVA3oldDMwoLT", &Tau_byMediumIsolationMVA3oldDMwoLT, &b_Tau_byMediumIsolationMVA3oldDMwoLT);
   BOOM->SetBranchAddress("Tau_byTightIsolationMVA3newDMwLT", &Tau_byTightIsolationMVA3newDMwLT, &b_Tau_byTightIsolationMVA3newDMwLT);
   BOOM->SetBranchAddress("Tau_byTightIsolationMVA3newDMwoLT", &Tau_byTightIsolationMVA3newDMwoLT, &b_Tau_byTightIsolationMVA3newDMwoLT);
   BOOM->SetBranchAddress("Tau_byTightIsolationMva3oldDMwLT", &Tau_byTightIsolationMva3oldDMwLT, &b_Tau_byTightIsolationMva3oldDMwLT);
   BOOM->SetBranchAddress("Tau_byTightIsolationMVA3oldDMwoLT", &Tau_byTightIsolationMVA3oldDMwoLT, &b_Tau_byTightIsolationMVA3oldDMwoLT);
   BOOM->SetBranchAddress("Tau_againstMuonLoose2", &Tau_againstMuonLoose2, &b_Tau_againstMuonLoose2);
   BOOM->SetBranchAddress("Tau_againstMuonLoose3", &Tau_againstMuonLoose3, &b_Tau_againstMuonLoose3);
   BOOM->SetBranchAddress("Tau_againstMuonTight2", &Tau_againstMuonTight2, &b_Tau_againstMuonTight2);
   BOOM->SetBranchAddress("Tau_againstElectronMVALooseMVA5", &Tau_againstElectronMVALooseMVA5, &b_Tau_againstElectronMVALooseMVA5);
   BOOM->SetBranchAddress("Tau_againstElectronMVAMediumMVA5", &Tau_againstElectronMVAMediumMVA5, &b_Tau_againstElectronMVAMediumMVA5);
   BOOM->SetBranchAddress("Tau_byVLooseCombinedIsolationDeltaBetaCorr", &Tau_byVLooseCombinedIsolationDeltaBetaCorr, &b_Tau_byVLooseCombinedIsolationDeltaBetaCorr);
   BOOM->SetBranchAddress("Tau_leadChargedCandPt", &Tau_leadChargedCandPt, &b_Tau_leadChargedCandPt);
   BOOM->SetBranchAddress("Tau_leadChargedCandCharge", &Tau_leadChargedCandCharge, &b_Tau_leadChargedCandCharge);
   BOOM->SetBranchAddress("Tau_leadChargedCandEta", &Tau_leadChargedCandEta, &b_Tau_leadChargedCandEta);
   BOOM->SetBranchAddress("Tau_leadChargedCandPhi", &Tau_leadChargedCandPhi, &b_Tau_leadChargedCandPhi);
   BOOM->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
   BOOM->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
   BOOM->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
   BOOM->SetBranchAddress("Jet_energy", &Jet_energy, &b_Jet_energy);
   BOOM->SetBranchAddress("Jet_bDiscriminator", &Jet_bDiscriminator, &b_Jet_bDiscriminator);
   BOOM->SetBranchAddress("Jet_mass", &Jet_mass, &b_Jet_mass);
   BOOM->SetBranchAddress("Jet_neutralHadEnergy", &Jet_neutralHadEnergy, &b_Jet_neutralHadEnergy);
   BOOM->SetBranchAddress("Jet_neutralEmEmEnergy", &Jet_neutralEmEmEnergy, &b_Jet_neutralEmEmEnergy);
   BOOM->SetBranchAddress("Jet_chargedHadronEnergy", &Jet_chargedHadronEnergy, &b_Jet_chargedHadronEnergy);
   BOOM->SetBranchAddress("Jet_chargedEmEnergy", &Jet_chargedEmEnergy, &b_Jet_chargedEmEnergy);
   BOOM->SetBranchAddress("Jet_muonEnergy", &Jet_muonEnergy, &b_Jet_muonEnergy);
   BOOM->SetBranchAddress("Jet_electronEnergy", &Jet_electronEnergy, &b_Jet_electronEnergy);
   BOOM->SetBranchAddress("Jet_photonEnergy", &Jet_photonEnergy, &b_Jet_photonEnergy);
   BOOM->SetBranchAddress("UncorrJet_pt", &UncorrJet_pt, &b_UncorrJet_pt);
   BOOM->SetBranchAddress("npuVertices", &npuVertices, &b_npuVertices);
   BOOM->SetBranchAddress("trueInteractions", &trueInteractions, &b_trueInteractions);
   BOOM->SetBranchAddress("ootnpuVertices", &ootnpuVertices, &b_ootnpuVertices);
   BOOM->SetBranchAddress("npuVerticesp1", &npuVerticesp1, &b_npuVerticesp1);
   BOOM->SetBranchAddress("bestVertices", &bestVertices, &b_bestVertices);
   BOOM->SetBranchAddress("Met_pt", &Met_pt, &b_Met_pt);
   BOOM->SetBranchAddress("Met_sumEt", &Met_sumEt, &b_Met_sumEt);
   BOOM->SetBranchAddress("Met_phi", &Met_phi, &b_Met_phi);
   BOOM->SetBranchAddress("Met_px", &Met_px, &b_Met_px);
   BOOM->SetBranchAddress("Met_py", &Met_py, &b_Met_py);
   BOOM->SetBranchAddress("Met_pz", &Met_pz, &b_Met_pz);
   BOOM->SetBranchAddress("Gen_Met", &Gen_Met, &b_Gen_Met);
   BOOM->SetBranchAddress("Met_shiftedPtUp", &Met_shiftedPtUp, &b_Met_shiftedPtUp);
   BOOM->SetBranchAddress("Met_shiftedPtDown", &Met_shiftedPtDown, &b_Met_shiftedPtDown);
};
