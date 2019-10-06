#define SUSYTauTauAnalysis_cxx
#include "SUSYTauTauAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

void SUSYTauTauAnalysis::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   // create the root file 
   TFile *file=new TFile("histograms_B.root","RECREATE");
 
   Int_t nTauEvents = 0;
   Int_t nTotTaus = 0;
   Int_t nSMTaus = 0;
   Int_t nSUSYTaus = 0;
   Int_t nElectrons = 0;
   Int_t nSMElectrons = 0;
   Int_t TauDecayType = 0;
   vector<Int_t> TotTauIndexVector;
   vector<Int_t> SMTauIndexVector;
   vector<Int_t> TauIndexVector;
   vector<Int_t> SusyIndexVector;
   vector<Int_t> MuonIndexVector;
   vector<Int_t> ElectronIndexVector;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   // for (Long64_t jentry=0; jentry<100.; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (Cut(ientry) < 0) continue;

      // Loop on the true particles
      for (Int_t i=0; i<nMc; i++) {

	Int_t m = mothMc[i];

       	if ( fabs(idMc[i])==15 )
	  TotTauIndexVector.push_back(i);

       	if ( fabs(idMc[i])==15 && idMc[m]<100)
	  SMTauIndexVector.push_back(i);

       	if ( fabs(idMc[i])==15 && ((idMc[m]>1000001 && idMc[m]<1000039) || (idMc[m]>2000001 && idMc[m]<2000015)) )
	  TauIndexVector.push_back(i);

	if ( fabs(idMc[i])==11 && fabs(idMc[m])==15 )
	  MuonIndexVector.push_back(i);

	if ( fabs(idMc[i])==13 && fabs(idMc[m])==15 )
	  ElectronIndexVector.push_back(i);

      } // end of loop on true particles
      // std::cout << endl;

      if (TauIndexVector.size() > 1) {
      	nTauEvents++;
	// Loop over list of Taus
	for (Int_t ii=0; ii<TauIndexVector.size(); ii++) {
	  Int_t TauIndex = TauIndexVector[ii];
	  Int_t MothIndex = mothMc[TauIndex];
	  Int_t GrandMothIndex = mothMc[MothIndex];

	  for (Int_t jj=ii+1; jj<TauIndexVector.size(); jj++) {
	    Int_t SecondTauIndex = TauIndexVector[jj];
	    Int_t SecMothIndex = mothMc[SecondTauIndex];
	    Int_t SecGrandMothIndex = mothMc[SecMothIndex];
	    // cout << MothIndex << "\t" << GrandMothIndex << "\t" << SecMothIndex << "\t" << SecGrandMotherIndex << endl;
	    // check that they come from same SUSY cascade
	    if (MothIndex == SecMothIndex) continue;
	    else if (GrandMothIndex == SecGrandMothIndex) continue;
	    else if (MothIndex == SecGrandMothIndex || GrandMothIndex == MothIndex) {
	      SusyIndexVector.push_back(TauIndex);
	      SusyIndexVector.push_back(SecondTauIndex);
	      continue;
	    }
	    else continue;
	  }
	
	} // end loop over list of Taus
      }

      if (SusyIndexVector.size()>0) {
	// cout << "found " << SusyIndexVector.size() << " SUSY Taus in the event" << endl;
	nSUSYTaus++;
      }

      if (TotTauIndexVector.size()>0)
	nTotTaus++;

      if (SMTauIndexVector.size()>0)
	nSMTaus++;

      if (ElectronIndexVector.size()>0)
	nElectrons++;

      // Loop over Muons
       if (MuonIndexVector.size()>0) {
       	for (Int_t mm=0; mm<MuonIndexVector.size(); mm++) {
       	  Int_t MuonIndex = MuonIndexVector[mm];
       	  Int_t MuonMothIndex = mothMc[MuonIndex];
       	  for (Int_t tt=0; tt<SusyIndexVector.size(); tt++) {
       	    Int_t SusyTauIndex = SusyIndexVector[tt];
	    if (MuonMothIndex==SusyTauIndex)
	      TauDecayType = 1;
       	  }
       	}
       }

      // Loop over Electrons
       if (ElectronIndexVector.size()>0) {
       	for (Int_t ee=0; ee<ElectronIndexVector.size(); ee++) {
       	  Int_t ElectronIndex = ElectronIndexVector[ee];
       	  Int_t ElectronMothIndex = mothMc[ElectronIndex];
       	  for (Int_t tt=0; tt<SusyIndexVector.size(); tt++) {
       	    Int_t SusyTauIndex = SusyIndexVector[tt];
	    // cout << ElectronMothIndex << " and\t" << SusyTauIndex << "\t";
	    if (ElectronMothIndex==SusyTauIndex)
	      TauDecayType = 2;
       	  }
       	}
       }

       // Loop over SM Electrons
       if (ElectronIndexVector.size()>0) {
	 for (Int_t ee=0; ee<ElectronIndexVector.size(); ee++) {
	   Int_t ElectronIndex = ElectronIndexVector[ee];
	   Int_t ElectronMothIndex = mothMc[ElectronIndex];
	   for (Int_t tt=0; tt<SMTauIndexVector.size(); tt++) {
	     Int_t SMTauIndex = SMTauIndexVector[tt];
	     // cout << ElectronMothIndex << " and\t" << SusyTauIndex << "\t";
	     if (ElectronMothIndex==SMTauIndex)
	       nSMElectrons++;
	     // TauDecayType = 2;
	   }
	 }
       }

       if (TauDecayType>0)
	 std::cout << "Tau Decay Type " << TauDecayType << "\t";
       
       TotTauIndexVector.clear();
       SMTauIndexVector.clear();
       TauIndexVector.clear();
       SusyIndexVector.clear();
       MuonIndexVector.clear();
       ElectronIndexVector.clear();

   } // end of loop over events



//    std::cout << "There are " << nSMTaus << " events with SM Taus," << std::endl;   
//    std::cout << "There are " << nTotTaus << " events with Taus," << std::endl;   
   std::cout << "There are " << nSMElectrons << " events with Electrons," << std::endl;   
   std::cout << "There are " << nElectrons << " events with Electrons," << std::endl;   
   std::cout << "There are " << nTauEvents << " events with SUSY Taus, but only " << nSUSYTaus << " are of the same cascade" << std::endl; 

   
   file->Write();
   delete file;
  
} // end of Loop

void SUSYTauTauAnalysis::dumpNewTree(TTree *tree)
{

  cout << " entering bookHist " << endl;
  
  
  cout << "Booking events tree" << endl;
  
  fTree = new TTree("events", "events"); 
  fTree->Branch("run",    &fRunnumber, "run/I");

   fChain->SetBranchAddress("nl1Technical", &nl1Technical, &b_nl1Technical);
   fChain->SetBranchAddress("l1Technical", l1Technical, &b_l1Technical);
   fChain->SetBranchAddress("nl1Global", &nl1Global, &b_nl1Global);
   fChain->SetBranchAddress("l1Global", l1Global, &b_l1Global);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("orbitNumber", &orbitNumber, &b_orbitNumber);
   fChain->SetBranchAddress("nMc", &nMc, &b_nMc);
   fChain->SetBranchAddress("pMc", pMc, &b_pMc);
   fChain->SetBranchAddress("massMc", massMc, &b_massMc);
   fChain->SetBranchAddress("thetaMc", thetaMc, &b_thetaMc);
   fChain->SetBranchAddress("etaMc", etaMc, &b_etaMc);
   fChain->SetBranchAddress("phiMc", phiMc, &b_phiMc);
   fChain->SetBranchAddress("energyMc", energyMc, &b_energyMc);
   fChain->SetBranchAddress("idMc", idMc, &b_idMc);
   fChain->SetBranchAddress("mothMc", mothMc, &b_mothMc);
   fChain->SetBranchAddress("nDauMc", nDauMc, &b_nDauMc);
   fChain->SetBranchAddress("statusMc", statusMc, &b_statusMc);
   fChain->SetBranchAddress("xMc", xMc, &b_xMc);
   fChain->SetBranchAddress("yMc", yMc, &b_yMc);
   fChain->SetBranchAddress("zMc", zMc, &b_zMc);
   fChain->SetBranchAddress("nTrg", &nTrg, &b_nTrg);
   fChain->SetBranchAddress("firedTrg", firedTrg, &b_firedTrg);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("chargeEle", chargeEle, &b_chargeEle);
   fChain->SetBranchAddress("energyEle", energyEle, &b_energyEle);
   fChain->SetBranchAddress("etEle", etEle, &b_etEle);
   fChain->SetBranchAddress("momentumEle", momentumEle, &b_momentumEle);
   fChain->SetBranchAddress("thetaEle", thetaEle, &b_thetaEle);
   fChain->SetBranchAddress("etaEle", etaEle, &b_etaEle);
   fChain->SetBranchAddress("phiEle", phiEle, &b_phiEle);
   fChain->SetBranchAddress("pxEle", pxEle, &b_pxEle);
   fChain->SetBranchAddress("pyEle", pyEle, &b_pyEle);
   fChain->SetBranchAddress("pzEle", pzEle, &b_pzEle);
   fChain->SetBranchAddress("vertexXEle", vertexXEle, &b_vertexXEle);
   fChain->SetBranchAddress("vertexYEle", vertexYEle, &b_vertexYEle);
   fChain->SetBranchAddress("vertexZEle", vertexZEle, &b_vertexZEle);
   fChain->SetBranchAddress("massEle", massEle, &b_massEle);
   fChain->SetBranchAddress("fiducialFlagsEle", fiducialFlagsEle, &b_fiducialFlagsEle);
   fChain->SetBranchAddress("recoFlagsEle", recoFlagsEle, &b_recoFlagsEle);
   fChain->SetBranchAddress("energyCorrectionsEle", energyCorrectionsEle, &b_energyCorrectionsEle);
   fChain->SetBranchAddress("esEnergyEle", esEnergyEle, &b_esEnergyEle);
   fChain->SetBranchAddress("superClusterIndexEle", superClusterIndexEle, &b_superClusterIndexEle);
   fChain->SetBranchAddress("PFsuperClusterIndexEle", PFsuperClusterIndexEle, &b_PFsuperClusterIndexEle);
   fChain->SetBranchAddress("trackIndexEle", trackIndexEle, &b_trackIndexEle);
   fChain->SetBranchAddress("gsfTrackIndexEle", gsfTrackIndexEle, &b_gsfTrackIndexEle);
   fChain->SetBranchAddress("classificationEle", classificationEle, &b_classificationEle);
   fChain->SetBranchAddress("standardClassificationEle", standardClassificationEle, &b_standardClassificationEle);
   fChain->SetBranchAddress("hOverEEle", hOverEEle, &b_hOverEEle);
   fChain->SetBranchAddress("eSuperClusterOverPEle", eSuperClusterOverPEle, &b_eSuperClusterOverPEle);
   fChain->SetBranchAddress("eSeedOverPoutEle", eSeedOverPoutEle, &b_eSeedOverPoutEle);
   fChain->SetBranchAddress("deltaEtaAtVtxEle", deltaEtaAtVtxEle, &b_deltaEtaAtVtxEle);
   fChain->SetBranchAddress("deltaPhiAtVtxEle", deltaPhiAtVtxEle, &b_deltaPhiAtVtxEle);
   fChain->SetBranchAddress("deltaEtaAtCaloEle", deltaEtaAtCaloEle, &b_deltaEtaAtCaloEle);
   fChain->SetBranchAddress("deltaPhiAtCaloEle", deltaPhiAtCaloEle, &b_deltaPhiAtCaloEle);
   fChain->SetBranchAddress("tipEle", tipEle, &b_tipEle);
   fChain->SetBranchAddress("dr03TkSumPtEle", dr03TkSumPtEle, &b_dr03TkSumPtEle);
   fChain->SetBranchAddress("dr03EcalRecHitSumEtEle", dr03EcalRecHitSumEtEle, &b_dr03EcalRecHitSumEtEle);
   fChain->SetBranchAddress("dr03HcalTowerSumEtEle", dr03HcalTowerSumEtEle, &b_dr03HcalTowerSumEtEle);
   fChain->SetBranchAddress("dr04TkSumPtEle", dr04TkSumPtEle, &b_dr04TkSumPtEle);
   fChain->SetBranchAddress("dr04EcalRecHitSumEtEle", dr04EcalRecHitSumEtEle, &b_dr04EcalRecHitSumEtEle);
   fChain->SetBranchAddress("dr04HcalTowerSumEtEle", dr04HcalTowerSumEtEle, &b_dr04HcalTowerSumEtEle);
   fChain->SetBranchAddress("scBasedEcalSum03Ele", scBasedEcalSum03Ele, &b_scBasedEcalSum03Ele);
   fChain->SetBranchAddress("scBasedEcalSum04Ele", scBasedEcalSum04Ele, &b_scBasedEcalSum04Ele);
   fChain->SetBranchAddress("scHaloBasedEcalSum03Ele", scHaloBasedEcalSum03Ele, &b_scHaloBasedEcalSum03Ele);
   fChain->SetBranchAddress("scHaloBasedEcalSum04Ele", scHaloBasedEcalSum04Ele, &b_scHaloBasedEcalSum04Ele);
   fChain->SetBranchAddress("eleIdCutsEle", eleIdCutsEle, &b_eleIdCutsEle);
   fChain->SetBranchAddress("eleIdLikelihoodEle", eleIdLikelihoodEle, &b_eleIdLikelihoodEle);
   fChain->SetBranchAddress("pflowMVAEle", pflowMVAEle, &b_pflowMVAEle);
   fChain->SetBranchAddress("nPFEle", &nPFEle, &b_nPFEle);
   fChain->SetBranchAddress("chargePFEle", chargePFEle, &b_chargePFEle);
   fChain->SetBranchAddress("energyPFEle", energyPFEle, &b_energyPFEle);
   fChain->SetBranchAddress("etPFEle", etPFEle, &b_etPFEle);
   fChain->SetBranchAddress("momentumPFEle", momentumPFEle, &b_momentumPFEle);
   fChain->SetBranchAddress("thetaPFEle", thetaPFEle, &b_thetaPFEle);
   fChain->SetBranchAddress("etaPFEle", etaPFEle, &b_etaPFEle);
   fChain->SetBranchAddress("phiPFEle", phiPFEle, &b_phiPFEle);
   fChain->SetBranchAddress("pxPFEle", pxPFEle, &b_pxPFEle);
   fChain->SetBranchAddress("pyPFEle", pyPFEle, &b_pyPFEle);
   fChain->SetBranchAddress("pzPFEle", pzPFEle, &b_pzPFEle);
   fChain->SetBranchAddress("vertexXPFEle", vertexXPFEle, &b_vertexXPFEle);
   fChain->SetBranchAddress("vertexYPFEle", vertexYPFEle, &b_vertexYPFEle);
   fChain->SetBranchAddress("vertexZPFEle", vertexZPFEle, &b_vertexZPFEle);
   fChain->SetBranchAddress("massPFEle", massPFEle, &b_massPFEle);
   fChain->SetBranchAddress("MvaOutputPFEle", MvaOutputPFEle, &b_MvaOutputPFEle);
   fChain->SetBranchAddress("PS1EnergyPFEle", PS1EnergyPFEle, &b_PS1EnergyPFEle);
   fChain->SetBranchAddress("PS2EnergyPFEle", PS2EnergyPFEle, &b_PS2EnergyPFEle);
   fChain->SetBranchAddress("EcalEnergyPFEle", EcalEnergyPFEle, &b_EcalEnergyPFEle);
   fChain->SetBranchAddress("EcalElectronEnergyPFEle", EcalElectronEnergyPFEle, &b_EcalElectronEnergyPFEle);
   fChain->SetBranchAddress("gsfTrackIndexPFEle", gsfTrackIndexPFEle, &b_gsfTrackIndexPFEle);
   fChain->SetBranchAddress("trackIndexPFEle", trackIndexPFEle, &b_trackIndexPFEle);
   fChain->SetBranchAddress("nSC", &nSC, &b_nSC);
   fChain->SetBranchAddress("nBCSC", nBCSC, &b_nBCSC);
   fChain->SetBranchAddress("nCrystalsSC", nCrystalsSC, &b_nCrystalsSC);
   fChain->SetBranchAddress("iAlgoSC", iAlgoSC, &b_iAlgoSC);
   fChain->SetBranchAddress("rawEnergySC", rawEnergySC, &b_rawEnergySC);
   fChain->SetBranchAddress("energySC", energySC, &b_energySC);
   fChain->SetBranchAddress("etaSC", etaSC, &b_etaSC);
   fChain->SetBranchAddress("thetaSC", thetaSC, &b_thetaSC);
   fChain->SetBranchAddress("phiSC", phiSC, &b_phiSC);
   fChain->SetBranchAddress("e3x3SC", e3x3SC, &b_e3x3SC);
   fChain->SetBranchAddress("e5x5SC", e5x5SC, &b_e5x5SC);
   fChain->SetBranchAddress("eMaxSC", eMaxSC, &b_eMaxSC);
   fChain->SetBranchAddress("e2x2SC", e2x2SC, &b_e2x2SC);
   fChain->SetBranchAddress("e2ndSC", e2ndSC, &b_e2ndSC);
   fChain->SetBranchAddress("covIEtaIEtaSC", covIEtaIEtaSC, &b_covIEtaIEtaSC);
   fChain->SetBranchAddress("covIEtaIPhiSC", covIEtaIPhiSC, &b_covIEtaIPhiSC);
   fChain->SetBranchAddress("covIPhiIPhiSC", covIPhiIPhiSC, &b_covIPhiIPhiSC);
   fChain->SetBranchAddress("hOverESC", hOverESC, &b_hOverESC);
   fChain->SetBranchAddress("recoFlagSC", recoFlagSC, &b_recoFlagSC);
   fChain->SetBranchAddress("channelStatusSC", channelStatusSC, &b_channelStatusSC);
   fChain->SetBranchAddress("timeSC", timeSC, &b_timeSC);
   fChain->SetBranchAddress("chi2ProbSC", chi2ProbSC, &b_chi2ProbSC);
   fChain->SetBranchAddress("seedEnergySC", seedEnergySC, &b_seedEnergySC);
   fChain->SetBranchAddress("idClosProblSC", idClosProblSC, &b_idClosProblSC);
   fChain->SetBranchAddress("sevClosProblSC", sevClosProblSC, &b_sevClosProblSC);
   fChain->SetBranchAddress("fracClosProblSC", fracClosProblSC, &b_fracClosProblSC);
   fChain->SetBranchAddress("trackIndexSC", trackIndexSC, &b_trackIndexSC);
   fChain->SetBranchAddress("trackDeltaRSC", trackDeltaRSC, &b_trackDeltaRSC);
   fChain->SetBranchAddress("trackDeltaPhiSC", trackDeltaPhiSC, &b_trackDeltaPhiSC);
   fChain->SetBranchAddress("trackDeltaEtaSC", trackDeltaEtaSC, &b_trackDeltaEtaSC);
   fChain->SetBranchAddress("gsfTrackIndexSC", gsfTrackIndexSC, &b_gsfTrackIndexSC);
   fChain->SetBranchAddress("gsfTrackDeltaRSC", gsfTrackDeltaRSC, &b_gsfTrackDeltaRSC);
   fChain->SetBranchAddress("gsfTrackDeltaPhiSC", gsfTrackDeltaPhiSC, &b_gsfTrackDeltaPhiSC);
   fChain->SetBranchAddress("gsfTrackDeltaEtaSC", gsfTrackDeltaEtaSC, &b_gsfTrackDeltaEtaSC);
   fChain->SetBranchAddress("pxVtxPropagatedNegChargeSC", pxVtxPropagatedNegChargeSC, &b_pxVtxPropagatedNegChargeSC);
   fChain->SetBranchAddress("pyVtxPropagatedNegChargeSC", pyVtxPropagatedNegChargeSC, &b_pyVtxPropagatedNegChargeSC);
   fChain->SetBranchAddress("pzVtxPropagatedNegChargeSC", pzVtxPropagatedNegChargeSC, &b_pzVtxPropagatedNegChargeSC);
   fChain->SetBranchAddress("pxVtxPropagatedPosChargeSC", pxVtxPropagatedPosChargeSC, &b_pxVtxPropagatedPosChargeSC);
   fChain->SetBranchAddress("pyVtxPropagatedPosChargeSC", pyVtxPropagatedPosChargeSC, &b_pyVtxPropagatedPosChargeSC);
   fChain->SetBranchAddress("pzVtxPropagatedPosChargeSC", pzVtxPropagatedPosChargeSC, &b_pzVtxPropagatedPosChargeSC);
   fChain->SetBranchAddress("nPFSC", &nPFSC, &b_nPFSC);
   fChain->SetBranchAddress("nBCPFSC", nBCPFSC, &b_nBCPFSC);
   fChain->SetBranchAddress("nCrystalsPFSC", nCrystalsPFSC, &b_nCrystalsPFSC);
   fChain->SetBranchAddress("iAlgoPFSC", iAlgoPFSC, &b_iAlgoPFSC);
   fChain->SetBranchAddress("rawEnergyPFSC", rawEnergyPFSC, &b_rawEnergyPFSC);
   fChain->SetBranchAddress("energyPFSC", energyPFSC, &b_energyPFSC);
   fChain->SetBranchAddress("etaPFSC", etaPFSC, &b_etaPFSC);
   fChain->SetBranchAddress("thetaPFSC", thetaPFSC, &b_thetaPFSC);
   fChain->SetBranchAddress("phiPFSC", phiPFSC, &b_phiPFSC);
   fChain->SetBranchAddress("e3x3PFSC", e3x3PFSC, &b_e3x3PFSC);
   fChain->SetBranchAddress("e5x5PFSC", e5x5PFSC, &b_e5x5PFSC);
   fChain->SetBranchAddress("eMaxPFSC", eMaxPFSC, &b_eMaxPFSC);
   fChain->SetBranchAddress("e2x2PFSC", e2x2PFSC, &b_e2x2PFSC);
   fChain->SetBranchAddress("e2ndPFSC", e2ndPFSC, &b_e2ndPFSC);
   fChain->SetBranchAddress("covIEtaIEtaPFSC", covIEtaIEtaPFSC, &b_covIEtaIEtaPFSC);
   fChain->SetBranchAddress("covIEtaIPhiPFSC", covIEtaIPhiPFSC, &b_covIEtaIPhiPFSC);
   fChain->SetBranchAddress("covIPhiIPhiPFSC", covIPhiIPhiPFSC, &b_covIPhiIPhiPFSC);
   fChain->SetBranchAddress("hOverEPFSC", hOverEPFSC, &b_hOverEPFSC);
   fChain->SetBranchAddress("recoFlagPFSC", recoFlagPFSC, &b_recoFlagPFSC);
   fChain->SetBranchAddress("channelStatusPFSC", channelStatusPFSC, &b_channelStatusPFSC);
   fChain->SetBranchAddress("timePFSC", timePFSC, &b_timePFSC);
   fChain->SetBranchAddress("chi2ProbPFSC", chi2ProbPFSC, &b_chi2ProbPFSC);
   fChain->SetBranchAddress("seedEnergyPFSC", seedEnergyPFSC, &b_seedEnergyPFSC);
   fChain->SetBranchAddress("idClosProblPFSC", idClosProblPFSC, &b_idClosProblPFSC);
   fChain->SetBranchAddress("sevClosProblPFSC", sevClosProblPFSC, &b_sevClosProblPFSC);
   fChain->SetBranchAddress("fracClosProblPFSC", fracClosProblPFSC, &b_fracClosProblPFSC);
   fChain->SetBranchAddress("pxVtxPropagatedNegChargePFSC", pxVtxPropagatedNegChargePFSC, &b_pxVtxPropagatedNegChargePFSC);
   fChain->SetBranchAddress("pyVtxPropagatedNegChargePFSC", pyVtxPropagatedNegChargePFSC, &b_pyVtxPropagatedNegChargePFSC);
   fChain->SetBranchAddress("pzVtxPropagatedNegChargePFSC", pzVtxPropagatedNegChargePFSC, &b_pzVtxPropagatedNegChargePFSC);
   fChain->SetBranchAddress("pxVtxPropagatedPosChargePFSC", pxVtxPropagatedPosChargePFSC, &b_pxVtxPropagatedPosChargePFSC);
   fChain->SetBranchAddress("pyVtxPropagatedPosChargePFSC", pyVtxPropagatedPosChargePFSC, &b_pyVtxPropagatedPosChargePFSC);
   fChain->SetBranchAddress("pzVtxPropagatedPosChargePFSC", pzVtxPropagatedPosChargePFSC, &b_pzVtxPropagatedPosChargePFSC);
   fChain->SetBranchAddress("nTrack", &nTrack, &b_nTrack);
   fChain->SetBranchAddress("pxTrack", pxTrack, &b_pxTrack);
   fChain->SetBranchAddress("pyTrack", pyTrack, &b_pyTrack);
   fChain->SetBranchAddress("pzTrack", pzTrack, &b_pzTrack);
   fChain->SetBranchAddress("vtxIndexTrack", vtxIndexTrack, &b_vtxIndexTrack);
   fChain->SetBranchAddress("vtxWeightTrack", vtxWeightTrack, &b_vtxWeightTrack);
   fChain->SetBranchAddress("chargeTrack", chargeTrack, &b_chargeTrack);
   fChain->SetBranchAddress("ptErrorTrack", ptErrorTrack, &b_ptErrorTrack);
   fChain->SetBranchAddress("trackValidHitsTrack", trackValidHitsTrack, &b_trackValidHitsTrack);
   fChain->SetBranchAddress("trackLostHitsTrack", trackLostHitsTrack, &b_trackLostHitsTrack);
   fChain->SetBranchAddress("trackNormalizedChi2Track", trackNormalizedChi2Track, &b_trackNormalizedChi2Track);
   fChain->SetBranchAddress("qualityMaskTrack", qualityMaskTrack, &b_qualityMaskTrack);
   fChain->SetBranchAddress("trackDxyTrack", trackDxyTrack, &b_trackDxyTrack);
   fChain->SetBranchAddress("trackD0Track", trackD0Track, &b_trackD0Track);
   fChain->SetBranchAddress("trackDszTrack", trackDszTrack, &b_trackDszTrack);
   fChain->SetBranchAddress("trackDzTrack", trackDzTrack, &b_trackDzTrack);
   fChain->SetBranchAddress("trackDxyErrorTrack", trackDxyErrorTrack, &b_trackDxyErrorTrack);
   fChain->SetBranchAddress("trackD0ErrorTrack", trackD0ErrorTrack, &b_trackD0ErrorTrack);
   fChain->SetBranchAddress("trackDszErrorTrack", trackDszErrorTrack, &b_trackDszErrorTrack);
   fChain->SetBranchAddress("trackDzErrorTrack", trackDzErrorTrack, &b_trackDzErrorTrack);
   fChain->SetBranchAddress("trackDxyPVTrack", trackDxyPVTrack, &b_trackDxyPVTrack);
   fChain->SetBranchAddress("trackDszPVTrack", trackDszPVTrack, &b_trackDszPVTrack);
   fChain->SetBranchAddress("trackDzPVTrack", trackDzPVTrack, &b_trackDzPVTrack);
   fChain->SetBranchAddress("trackVxTrack", trackVxTrack, &b_trackVxTrack);
   fChain->SetBranchAddress("trackVyTrack", trackVyTrack, &b_trackVyTrack);
   fChain->SetBranchAddress("trackVzTrack", trackVzTrack, &b_trackVzTrack);
   fChain->SetBranchAddress("pxAtOuterTrack", pxAtOuterTrack, &b_pxAtOuterTrack);
   fChain->SetBranchAddress("pyAtOuterTrack", pyAtOuterTrack, &b_pyAtOuterTrack);
   fChain->SetBranchAddress("pzAtOuterTrack", pzAtOuterTrack, &b_pzAtOuterTrack);
   fChain->SetBranchAddress("xAtOuterTrack", xAtOuterTrack, &b_xAtOuterTrack);
   fChain->SetBranchAddress("yAtOuterTrack", yAtOuterTrack, &b_yAtOuterTrack);
   fChain->SetBranchAddress("zAtOuterTrack", zAtOuterTrack, &b_zAtOuterTrack);
   fChain->SetBranchAddress("pxAtInnerTrack", pxAtInnerTrack, &b_pxAtInnerTrack);
   fChain->SetBranchAddress("pyAtInnerTrack", pyAtInnerTrack, &b_pyAtInnerTrack);
   fChain->SetBranchAddress("pzAtInnerTrack", pzAtInnerTrack, &b_pzAtInnerTrack);
   fChain->SetBranchAddress("xAtInnerTrack", xAtInnerTrack, &b_xAtInnerTrack);
   fChain->SetBranchAddress("yAtInnerTrack", yAtInnerTrack, &b_yAtInnerTrack);
   fChain->SetBranchAddress("zAtInnerTrack", zAtInnerTrack, &b_zAtInnerTrack);
   fChain->SetBranchAddress("recHitsSizeTrack", recHitsSizeTrack, &b_recHitsSizeTrack);
   fChain->SetBranchAddress("isPixB1Track", isPixB1Track, &b_isPixB1Track);
   fChain->SetBranchAddress("isPixB2Track", isPixB2Track, &b_isPixB2Track);
   fChain->SetBranchAddress("isPixE1Track", isPixE1Track, &b_isPixE1Track);
   fChain->SetBranchAddress("isPixE2Track", isPixE2Track, &b_isPixE2Track);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsTrack", numberOfValidPixelBarrelHitsTrack, &b_numberOfValidPixelBarrelHitsTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsTrack", numberOfValidPixelEndcapHitsTrack, &b_numberOfValidPixelEndcapHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsTrack", numberOfValidStripTIBHitsTrack, &b_numberOfValidStripTIBHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsTrack", numberOfValidStripTIDHitsTrack, &b_numberOfValidStripTIDHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsTrack", numberOfValidStripTOBHitsTrack, &b_numberOfValidStripTOBHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsTrack", numberOfValidStripTECHitsTrack, &b_numberOfValidStripTECHitsTrack);
   fChain->SetBranchAddress("truncatedDeDxTrack", truncatedDeDxTrack, &b_truncatedDeDxTrack);
   fChain->SetBranchAddress("truncatedDeDxErrorTrack", truncatedDeDxErrorTrack, &b_truncatedDeDxErrorTrack);
   fChain->SetBranchAddress("truncatedDeDxNoMTrack", truncatedDeDxNoMTrack, &b_truncatedDeDxNoMTrack);
   fChain->SetBranchAddress("medianDeDxTrack", medianDeDxTrack, &b_medianDeDxTrack);
   fChain->SetBranchAddress("medianDeDxErrorTrack", medianDeDxErrorTrack, &b_medianDeDxErrorTrack);
   fChain->SetBranchAddress("medianDeDxNoMTrack", medianDeDxNoMTrack, &b_medianDeDxNoMTrack);
   fChain->SetBranchAddress("harmonic2DeDxTrack", harmonic2DeDxTrack, &b_harmonic2DeDxTrack);
   fChain->SetBranchAddress("harmonic2DeDxErrorTrack", harmonic2DeDxErrorTrack, &b_harmonic2DeDxErrorTrack);
   fChain->SetBranchAddress("harmonic2DeDxNoMTrack", harmonic2DeDxNoMTrack, &b_harmonic2DeDxNoMTrack);
   fChain->SetBranchAddress("nGsfTrack", &nGsfTrack, &b_nGsfTrack);
   fChain->SetBranchAddress("pxGsfTrack", pxGsfTrack, &b_pxGsfTrack);
   fChain->SetBranchAddress("pyGsfTrack", pyGsfTrack, &b_pyGsfTrack);
   fChain->SetBranchAddress("pzGsfTrack", pzGsfTrack, &b_pzGsfTrack);
   fChain->SetBranchAddress("chargeGsfTrack", chargeGsfTrack, &b_chargeGsfTrack);
   fChain->SetBranchAddress("ptErrorGsfTrack", ptErrorGsfTrack, &b_ptErrorGsfTrack);
   fChain->SetBranchAddress("trackValidHitsGsfTrack", trackValidHitsGsfTrack, &b_trackValidHitsGsfTrack);
   fChain->SetBranchAddress("trackLostHitsGsfTrack", trackLostHitsGsfTrack, &b_trackLostHitsGsfTrack);
   fChain->SetBranchAddress("trackNormalizedChi2GsfTrack", trackNormalizedChi2GsfTrack, &b_trackNormalizedChi2GsfTrack);
   fChain->SetBranchAddress("qualityMaskGsfTrack", qualityMaskGsfTrack, &b_qualityMaskGsfTrack);
   fChain->SetBranchAddress("trackDxyGsfTrack", trackDxyGsfTrack, &b_trackDxyGsfTrack);
   fChain->SetBranchAddress("trackD0GsfTrack", trackD0GsfTrack, &b_trackD0GsfTrack);
   fChain->SetBranchAddress("trackDszGsfTrack", trackDszGsfTrack, &b_trackDszGsfTrack);
   fChain->SetBranchAddress("trackDzGsfTrack", trackDzGsfTrack, &b_trackDzGsfTrack);
   fChain->SetBranchAddress("trackDxyErrorGsfTrack", trackDxyErrorGsfTrack, &b_trackDxyErrorGsfTrack);
   fChain->SetBranchAddress("trackD0ErrorGsfTrack", trackD0ErrorGsfTrack, &b_trackD0ErrorGsfTrack);
   fChain->SetBranchAddress("trackDszErrorGsfTrack", trackDszErrorGsfTrack, &b_trackDszErrorGsfTrack);
   fChain->SetBranchAddress("trackDzErrorGsfTrack", trackDzErrorGsfTrack, &b_trackDzErrorGsfTrack);
   fChain->SetBranchAddress("trackDxyPVGsfTrack", trackDxyPVGsfTrack, &b_trackDxyPVGsfTrack);
   fChain->SetBranchAddress("trackDszPVGsfTrack", trackDszPVGsfTrack, &b_trackDszPVGsfTrack);
   fChain->SetBranchAddress("trackDzPVGsfTrack", trackDzPVGsfTrack, &b_trackDzPVGsfTrack);
   fChain->SetBranchAddress("trackVxGsfTrack", trackVxGsfTrack, &b_trackVxGsfTrack);
   fChain->SetBranchAddress("trackVyGsfTrack", trackVyGsfTrack, &b_trackVyGsfTrack);
   fChain->SetBranchAddress("trackVzGsfTrack", trackVzGsfTrack, &b_trackVzGsfTrack);
   fChain->SetBranchAddress("pxAtOuterGsfTrack", pxAtOuterGsfTrack, &b_pxAtOuterGsfTrack);
   fChain->SetBranchAddress("pyAtOuterGsfTrack", pyAtOuterGsfTrack, &b_pyAtOuterGsfTrack);
   fChain->SetBranchAddress("pzAtOuterGsfTrack", pzAtOuterGsfTrack, &b_pzAtOuterGsfTrack);
   fChain->SetBranchAddress("xAtOuterGsfTrack", xAtOuterGsfTrack, &b_xAtOuterGsfTrack);
   fChain->SetBranchAddress("yAtOuterGsfTrack", yAtOuterGsfTrack, &b_yAtOuterGsfTrack);
   fChain->SetBranchAddress("zAtOuterGsfTrack", zAtOuterGsfTrack, &b_zAtOuterGsfTrack);
   fChain->SetBranchAddress("pxAtInnerGsfTrack", pxAtInnerGsfTrack, &b_pxAtInnerGsfTrack);
   fChain->SetBranchAddress("pyAtInnerGsfTrack", pyAtInnerGsfTrack, &b_pyAtInnerGsfTrack);
   fChain->SetBranchAddress("pzAtInnerGsfTrack", pzAtInnerGsfTrack, &b_pzAtInnerGsfTrack);
   fChain->SetBranchAddress("xAtInnerGsfTrack", xAtInnerGsfTrack, &b_xAtInnerGsfTrack);
   fChain->SetBranchAddress("yAtInnerGsfTrack", yAtInnerGsfTrack, &b_yAtInnerGsfTrack);
   fChain->SetBranchAddress("zAtInnerGsfTrack", zAtInnerGsfTrack, &b_zAtInnerGsfTrack);
   fChain->SetBranchAddress("recHitsSizeGsfTrack", recHitsSizeGsfTrack, &b_recHitsSizeGsfTrack);
   fChain->SetBranchAddress("isPixB1GsfTrack", isPixB1GsfTrack, &b_isPixB1GsfTrack);
   fChain->SetBranchAddress("isPixB2GsfTrack", isPixB2GsfTrack, &b_isPixB2GsfTrack);
   fChain->SetBranchAddress("isPixE1GsfTrack", isPixE1GsfTrack, &b_isPixE1GsfTrack);
   fChain->SetBranchAddress("isPixE2GsfTrack", isPixE2GsfTrack, &b_isPixE2GsfTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsGsfTrack", numberOfValidPixelBarrelHitsGsfTrack, &b_numberOfValidPixelBarrelHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsGsfTrack", numberOfValidPixelEndcapHitsGsfTrack, &b_numberOfValidPixelEndcapHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsGsfTrack", numberOfValidStripTIBHitsGsfTrack, &b_numberOfValidStripTIBHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsGsfTrack", numberOfValidStripTIDHitsGsfTrack, &b_numberOfValidStripTIDHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsGsfTrack", numberOfValidStripTOBHitsGsfTrack, &b_numberOfValidStripTOBHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsGsfTrack", numberOfValidStripTECHitsGsfTrack, &b_numberOfValidStripTECHitsGsfTrack);
   fChain->SetBranchAddress("chargeModeGsfTrack", chargeModeGsfTrack, &b_chargeModeGsfTrack);
   fChain->SetBranchAddress("pxModeGsfTrack", pxModeGsfTrack, &b_pxModeGsfTrack);
   fChain->SetBranchAddress("pyModeGsfTrack", pyModeGsfTrack, &b_pyModeGsfTrack);
   fChain->SetBranchAddress("pzModeGsfTrack", pzModeGsfTrack, &b_pzModeGsfTrack);
   fChain->SetBranchAddress("recoFlagsGsfTrack", recoFlagsGsfTrack, &b_recoFlagsGsfTrack);
   fChain->SetBranchAddress("nGlobalMuonTrack", &nGlobalMuonTrack, &b_nGlobalMuonTrack);
   fChain->SetBranchAddress("pxGlobalMuonTrack", pxGlobalMuonTrack, &b_pxGlobalMuonTrack);
   fChain->SetBranchAddress("pyGlobalMuonTrack", pyGlobalMuonTrack, &b_pyGlobalMuonTrack);
   fChain->SetBranchAddress("pzGlobalMuonTrack", pzGlobalMuonTrack, &b_pzGlobalMuonTrack);
   fChain->SetBranchAddress("chargeGlobalMuonTrack", chargeGlobalMuonTrack, &b_chargeGlobalMuonTrack);
   fChain->SetBranchAddress("ptErrorGlobalMuonTrack", ptErrorGlobalMuonTrack, &b_ptErrorGlobalMuonTrack);
   fChain->SetBranchAddress("trackValidHitsGlobalMuonTrack", trackValidHitsGlobalMuonTrack, &b_trackValidHitsGlobalMuonTrack);
   fChain->SetBranchAddress("trackLostHitsGlobalMuonTrack", trackLostHitsGlobalMuonTrack, &b_trackLostHitsGlobalMuonTrack);
   fChain->SetBranchAddress("trackNormalizedChi2GlobalMuonTrack", trackNormalizedChi2GlobalMuonTrack, &b_trackNormalizedChi2GlobalMuonTrack);
   fChain->SetBranchAddress("qualityMaskGlobalMuonTrack", qualityMaskGlobalMuonTrack, &b_qualityMaskGlobalMuonTrack);
   fChain->SetBranchAddress("trackDxyGlobalMuonTrack", trackDxyGlobalMuonTrack, &b_trackDxyGlobalMuonTrack);
   fChain->SetBranchAddress("trackD0GlobalMuonTrack", trackD0GlobalMuonTrack, &b_trackD0GlobalMuonTrack);
   fChain->SetBranchAddress("trackDszGlobalMuonTrack", trackDszGlobalMuonTrack, &b_trackDszGlobalMuonTrack);
   fChain->SetBranchAddress("trackDzGlobalMuonTrack", trackDzGlobalMuonTrack, &b_trackDzGlobalMuonTrack);
   fChain->SetBranchAddress("trackDxyErrorGlobalMuonTrack", trackDxyErrorGlobalMuonTrack, &b_trackDxyErrorGlobalMuonTrack);
   fChain->SetBranchAddress("trackD0ErrorGlobalMuonTrack", trackD0ErrorGlobalMuonTrack, &b_trackD0ErrorGlobalMuonTrack);
   fChain->SetBranchAddress("trackDszErrorGlobalMuonTrack", trackDszErrorGlobalMuonTrack, &b_trackDszErrorGlobalMuonTrack);
   fChain->SetBranchAddress("trackDzErrorGlobalMuonTrack", trackDzErrorGlobalMuonTrack, &b_trackDzErrorGlobalMuonTrack);
   fChain->SetBranchAddress("trackDxyPVGlobalMuonTrack", trackDxyPVGlobalMuonTrack, &b_trackDxyPVGlobalMuonTrack);
   fChain->SetBranchAddress("trackDszPVGlobalMuonTrack", trackDszPVGlobalMuonTrack, &b_trackDszPVGlobalMuonTrack);
   fChain->SetBranchAddress("trackDzPVGlobalMuonTrack", trackDzPVGlobalMuonTrack, &b_trackDzPVGlobalMuonTrack);
   fChain->SetBranchAddress("trackVxGlobalMuonTrack", trackVxGlobalMuonTrack, &b_trackVxGlobalMuonTrack);
   fChain->SetBranchAddress("trackVyGlobalMuonTrack", trackVyGlobalMuonTrack, &b_trackVyGlobalMuonTrack);
   fChain->SetBranchAddress("trackVzGlobalMuonTrack", trackVzGlobalMuonTrack, &b_trackVzGlobalMuonTrack);
   fChain->SetBranchAddress("pxAtOuterGlobalMuonTrack", pxAtOuterGlobalMuonTrack, &b_pxAtOuterGlobalMuonTrack);
   fChain->SetBranchAddress("pyAtOuterGlobalMuonTrack", pyAtOuterGlobalMuonTrack, &b_pyAtOuterGlobalMuonTrack);
   fChain->SetBranchAddress("pzAtOuterGlobalMuonTrack", pzAtOuterGlobalMuonTrack, &b_pzAtOuterGlobalMuonTrack);
   fChain->SetBranchAddress("xAtOuterGlobalMuonTrack", xAtOuterGlobalMuonTrack, &b_xAtOuterGlobalMuonTrack);
   fChain->SetBranchAddress("yAtOuterGlobalMuonTrack", yAtOuterGlobalMuonTrack, &b_yAtOuterGlobalMuonTrack);
   fChain->SetBranchAddress("zAtOuterGlobalMuonTrack", zAtOuterGlobalMuonTrack, &b_zAtOuterGlobalMuonTrack);
   fChain->SetBranchAddress("pxAtInnerGlobalMuonTrack", pxAtInnerGlobalMuonTrack, &b_pxAtInnerGlobalMuonTrack);
   fChain->SetBranchAddress("pyAtInnerGlobalMuonTrack", pyAtInnerGlobalMuonTrack, &b_pyAtInnerGlobalMuonTrack);
   fChain->SetBranchAddress("pzAtInnerGlobalMuonTrack", pzAtInnerGlobalMuonTrack, &b_pzAtInnerGlobalMuonTrack);
   fChain->SetBranchAddress("xAtInnerGlobalMuonTrack", xAtInnerGlobalMuonTrack, &b_xAtInnerGlobalMuonTrack);
   fChain->SetBranchAddress("yAtInnerGlobalMuonTrack", yAtInnerGlobalMuonTrack, &b_yAtInnerGlobalMuonTrack);
   fChain->SetBranchAddress("zAtInnerGlobalMuonTrack", zAtInnerGlobalMuonTrack, &b_zAtInnerGlobalMuonTrack);
   fChain->SetBranchAddress("recHitsSizeGlobalMuonTrack", recHitsSizeGlobalMuonTrack, &b_recHitsSizeGlobalMuonTrack);
   fChain->SetBranchAddress("isPixB1GlobalMuonTrack", isPixB1GlobalMuonTrack, &b_isPixB1GlobalMuonTrack);
   fChain->SetBranchAddress("isPixB2GlobalMuonTrack", isPixB2GlobalMuonTrack, &b_isPixB2GlobalMuonTrack);
   fChain->SetBranchAddress("isPixE1GlobalMuonTrack", isPixE1GlobalMuonTrack, &b_isPixE1GlobalMuonTrack);
   fChain->SetBranchAddress("isPixE2GlobalMuonTrack", isPixE2GlobalMuonTrack, &b_isPixE2GlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsGlobalMuonTrack", numberOfValidPixelBarrelHitsGlobalMuonTrack, &b_numberOfValidPixelBarrelHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsGlobalMuonTrack", numberOfValidPixelEndcapHitsGlobalMuonTrack, &b_numberOfValidPixelEndcapHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsGlobalMuonTrack", numberOfValidStripTIBHitsGlobalMuonTrack, &b_numberOfValidStripTIBHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsGlobalMuonTrack", numberOfValidStripTIDHitsGlobalMuonTrack, &b_numberOfValidStripTIDHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsGlobalMuonTrack", numberOfValidStripTOBHitsGlobalMuonTrack, &b_numberOfValidStripTOBHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsGlobalMuonTrack", numberOfValidStripTECHitsGlobalMuonTrack, &b_numberOfValidStripTECHitsGlobalMuonTrack);
   fChain->SetBranchAddress("nSTAMuonTrack", &nSTAMuonTrack, &b_nSTAMuonTrack);
   fChain->SetBranchAddress("pxSTAMuonTrack", pxSTAMuonTrack, &b_pxSTAMuonTrack);
   fChain->SetBranchAddress("pySTAMuonTrack", pySTAMuonTrack, &b_pySTAMuonTrack);
   fChain->SetBranchAddress("pzSTAMuonTrack", pzSTAMuonTrack, &b_pzSTAMuonTrack);
   fChain->SetBranchAddress("chargeSTAMuonTrack", chargeSTAMuonTrack, &b_chargeSTAMuonTrack);
   fChain->SetBranchAddress("ptErrorSTAMuonTrack", ptErrorSTAMuonTrack, &b_ptErrorSTAMuonTrack);
   fChain->SetBranchAddress("trackValidHitsSTAMuonTrack", trackValidHitsSTAMuonTrack, &b_trackValidHitsSTAMuonTrack);
   fChain->SetBranchAddress("trackLostHitsSTAMuonTrack", trackLostHitsSTAMuonTrack, &b_trackLostHitsSTAMuonTrack);
   fChain->SetBranchAddress("trackNormalizedChi2STAMuonTrack", trackNormalizedChi2STAMuonTrack, &b_trackNormalizedChi2STAMuonTrack);
   fChain->SetBranchAddress("qualityMaskSTAMuonTrack", qualityMaskSTAMuonTrack, &b_qualityMaskSTAMuonTrack);
   fChain->SetBranchAddress("trackDxySTAMuonTrack", trackDxySTAMuonTrack, &b_trackDxySTAMuonTrack);
   fChain->SetBranchAddress("trackD0STAMuonTrack", trackD0STAMuonTrack, &b_trackD0STAMuonTrack);
   fChain->SetBranchAddress("trackDszSTAMuonTrack", trackDszSTAMuonTrack, &b_trackDszSTAMuonTrack);
   fChain->SetBranchAddress("trackDzSTAMuonTrack", trackDzSTAMuonTrack, &b_trackDzSTAMuonTrack);
   fChain->SetBranchAddress("trackDxyErrorSTAMuonTrack", trackDxyErrorSTAMuonTrack, &b_trackDxyErrorSTAMuonTrack);
   fChain->SetBranchAddress("trackD0ErrorSTAMuonTrack", trackD0ErrorSTAMuonTrack, &b_trackD0ErrorSTAMuonTrack);
   fChain->SetBranchAddress("trackDszErrorSTAMuonTrack", trackDszErrorSTAMuonTrack, &b_trackDszErrorSTAMuonTrack);
   fChain->SetBranchAddress("trackDzErrorSTAMuonTrack", trackDzErrorSTAMuonTrack, &b_trackDzErrorSTAMuonTrack);
   fChain->SetBranchAddress("trackDxyPVSTAMuonTrack", trackDxyPVSTAMuonTrack, &b_trackDxyPVSTAMuonTrack);
   fChain->SetBranchAddress("trackDszPVSTAMuonTrack", trackDszPVSTAMuonTrack, &b_trackDszPVSTAMuonTrack);
   fChain->SetBranchAddress("trackDzPVSTAMuonTrack", trackDzPVSTAMuonTrack, &b_trackDzPVSTAMuonTrack);
   fChain->SetBranchAddress("trackVxSTAMuonTrack", trackVxSTAMuonTrack, &b_trackVxSTAMuonTrack);
   fChain->SetBranchAddress("trackVySTAMuonTrack", trackVySTAMuonTrack, &b_trackVySTAMuonTrack);
   fChain->SetBranchAddress("trackVzSTAMuonTrack", trackVzSTAMuonTrack, &b_trackVzSTAMuonTrack);
   fChain->SetBranchAddress("pxAtOuterSTAMuonTrack", pxAtOuterSTAMuonTrack, &b_pxAtOuterSTAMuonTrack);
   fChain->SetBranchAddress("pyAtOuterSTAMuonTrack", pyAtOuterSTAMuonTrack, &b_pyAtOuterSTAMuonTrack);
   fChain->SetBranchAddress("pzAtOuterSTAMuonTrack", pzAtOuterSTAMuonTrack, &b_pzAtOuterSTAMuonTrack);
   fChain->SetBranchAddress("xAtOuterSTAMuonTrack", xAtOuterSTAMuonTrack, &b_xAtOuterSTAMuonTrack);
   fChain->SetBranchAddress("yAtOuterSTAMuonTrack", yAtOuterSTAMuonTrack, &b_yAtOuterSTAMuonTrack);
   fChain->SetBranchAddress("zAtOuterSTAMuonTrack", zAtOuterSTAMuonTrack, &b_zAtOuterSTAMuonTrack);
   fChain->SetBranchAddress("pxAtInnerSTAMuonTrack", pxAtInnerSTAMuonTrack, &b_pxAtInnerSTAMuonTrack);
   fChain->SetBranchAddress("pyAtInnerSTAMuonTrack", pyAtInnerSTAMuonTrack, &b_pyAtInnerSTAMuonTrack);
   fChain->SetBranchAddress("pzAtInnerSTAMuonTrack", pzAtInnerSTAMuonTrack, &b_pzAtInnerSTAMuonTrack);
   fChain->SetBranchAddress("xAtInnerSTAMuonTrack", xAtInnerSTAMuonTrack, &b_xAtInnerSTAMuonTrack);
   fChain->SetBranchAddress("yAtInnerSTAMuonTrack", yAtInnerSTAMuonTrack, &b_yAtInnerSTAMuonTrack);
   fChain->SetBranchAddress("zAtInnerSTAMuonTrack", zAtInnerSTAMuonTrack, &b_zAtInnerSTAMuonTrack);
   fChain->SetBranchAddress("recHitsSizeSTAMuonTrack", recHitsSizeSTAMuonTrack, &b_recHitsSizeSTAMuonTrack);
   fChain->SetBranchAddress("isPixB1STAMuonTrack", isPixB1STAMuonTrack, &b_isPixB1STAMuonTrack);
   fChain->SetBranchAddress("isPixB2STAMuonTrack", isPixB2STAMuonTrack, &b_isPixB2STAMuonTrack);
   fChain->SetBranchAddress("isPixE1STAMuonTrack", isPixE1STAMuonTrack, &b_isPixE1STAMuonTrack);
   fChain->SetBranchAddress("isPixE2STAMuonTrack", isPixE2STAMuonTrack, &b_isPixE2STAMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsSTAMuonTrack", numberOfValidPixelBarrelHitsSTAMuonTrack, &b_numberOfValidPixelBarrelHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsSTAMuonTrack", numberOfValidPixelEndcapHitsSTAMuonTrack, &b_numberOfValidPixelEndcapHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsSTAMuonTrack", numberOfValidStripTIBHitsSTAMuonTrack, &b_numberOfValidStripTIBHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsSTAMuonTrack", numberOfValidStripTIDHitsSTAMuonTrack, &b_numberOfValidStripTIDHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsSTAMuonTrack", numberOfValidStripTOBHitsSTAMuonTrack, &b_numberOfValidStripTOBHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsSTAMuonTrack", numberOfValidStripTECHitsSTAMuonTrack, &b_numberOfValidStripTECHitsSTAMuonTrack);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("PVxPV", PVxPV, &b_PVxPV);
   fChain->SetBranchAddress("PVyPV", PVyPV, &b_PVyPV);
   fChain->SetBranchAddress("PVzPV", PVzPV, &b_PVzPV);
   fChain->SetBranchAddress("PVErrxPV", PVErrxPV, &b_PVErrxPV);
   fChain->SetBranchAddress("PVErryPV", PVErryPV, &b_PVErryPV);
   fChain->SetBranchAddress("PVErrzPV", PVErrzPV, &b_PVErrzPV);
   fChain->SetBranchAddress("SumPtPV", SumPtPV, &b_SumPtPV);
   fChain->SetBranchAddress("ndofPV", ndofPV, &b_ndofPV);
   fChain->SetBranchAddress("chi2PV", chi2PV, &b_chi2PV);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("chargeMuon", chargeMuon, &b_chargeMuon);
   fChain->SetBranchAddress("energyMuon", energyMuon, &b_energyMuon);
   fChain->SetBranchAddress("etMuon", etMuon, &b_etMuon);
   fChain->SetBranchAddress("momentumMuon", momentumMuon, &b_momentumMuon);
   fChain->SetBranchAddress("thetaMuon", thetaMuon, &b_thetaMuon);
   fChain->SetBranchAddress("etaMuon", etaMuon, &b_etaMuon);
   fChain->SetBranchAddress("phiMuon", phiMuon, &b_phiMuon);
   fChain->SetBranchAddress("pxMuon", pxMuon, &b_pxMuon);
   fChain->SetBranchAddress("pyMuon", pyMuon, &b_pyMuon);
   fChain->SetBranchAddress("pzMuon", pzMuon, &b_pzMuon);
   fChain->SetBranchAddress("vertexXMuon", vertexXMuon, &b_vertexXMuon);
   fChain->SetBranchAddress("vertexYMuon", vertexYMuon, &b_vertexYMuon);
   fChain->SetBranchAddress("vertexZMuon", vertexZMuon, &b_vertexZMuon);
   fChain->SetBranchAddress("massMuon", massMuon, &b_massMuon);
   fChain->SetBranchAddress("trackIndexMuon", trackIndexMuon, &b_trackIndexMuon);
   fChain->SetBranchAddress("standAloneTrackIndexMuon", standAloneTrackIndexMuon, &b_standAloneTrackIndexMuon);
   fChain->SetBranchAddress("combinedTrackIndexMuon", combinedTrackIndexMuon, &b_combinedTrackIndexMuon);
   fChain->SetBranchAddress("muonIdMuon", muonIdMuon, &b_muonIdMuon);
   fChain->SetBranchAddress("sumPt03Muon", sumPt03Muon, &b_sumPt03Muon);
   fChain->SetBranchAddress("emEt03Muon", emEt03Muon, &b_emEt03Muon);
   fChain->SetBranchAddress("hadEt03Muon", hadEt03Muon, &b_hadEt03Muon);
   fChain->SetBranchAddress("hoEt03Muon", hoEt03Muon, &b_hoEt03Muon);
   fChain->SetBranchAddress("nTrk03Muon", nTrk03Muon, &b_nTrk03Muon);
   fChain->SetBranchAddress("nJets03Muon", nJets03Muon, &b_nJets03Muon);
   fChain->SetBranchAddress("sumPt05Muon", sumPt05Muon, &b_sumPt05Muon);
   fChain->SetBranchAddress("emEt05Muon", emEt05Muon, &b_emEt05Muon);
   fChain->SetBranchAddress("hadEt05Muon", hadEt05Muon, &b_hadEt05Muon);
   fChain->SetBranchAddress("hoEt05Muon", hoEt05Muon, &b_hoEt05Muon);
   fChain->SetBranchAddress("nTrk05Muon", nTrk05Muon, &b_nTrk05Muon);
   fChain->SetBranchAddress("nJets05Muon", nJets05Muon, &b_nJets05Muon);
   fChain->SetBranchAddress("EcalExpDepoMuon", EcalExpDepoMuon, &b_EcalExpDepoMuon);
   fChain->SetBranchAddress("HcalExpDepoMuon", HcalExpDepoMuon, &b_HcalExpDepoMuon);
   fChain->SetBranchAddress("HoExpDepoMuon", HoExpDepoMuon, &b_HoExpDepoMuon);
   fChain->SetBranchAddress("emS9Muon", emS9Muon, &b_emS9Muon);
   fChain->SetBranchAddress("hadS9Muon", hadS9Muon, &b_hadS9Muon);
   fChain->SetBranchAddress("hoS9Muon", hoS9Muon, &b_hoS9Muon);
   fChain->SetBranchAddress("CaloCompMuon", CaloCompMuon, &b_CaloCompMuon);
   fChain->SetBranchAddress("nCaloTowers", &nCaloTowers, &b_nCaloTowers);
   fChain->SetBranchAddress("energyCaloTowers", energyCaloTowers, &b_energyCaloTowers);
   fChain->SetBranchAddress("xCaloTowers", xCaloTowers, &b_xCaloTowers);
   fChain->SetBranchAddress("yCaloTowers", yCaloTowers, &b_yCaloTowers);
   fChain->SetBranchAddress("zCaloTowers", zCaloTowers, &b_zCaloTowers);
   fChain->SetBranchAddress("CALOCaloTowers", CALOCaloTowers, &b_CALOCaloTowers);
   fChain->SetBranchAddress("CaloIndexCaloTowers", CaloIndexCaloTowers, &b_CaloIndexCaloTowers);
   fChain->SetBranchAddress("nMet", &nMet, &b_nMet);
   fChain->SetBranchAddress("chargeMet", chargeMet, &b_chargeMet);
   fChain->SetBranchAddress("energyMet", energyMet, &b_energyMet);
   fChain->SetBranchAddress("etMet", etMet, &b_etMet);
   fChain->SetBranchAddress("momentumMet", momentumMet, &b_momentumMet);
   fChain->SetBranchAddress("thetaMet", thetaMet, &b_thetaMet);
   fChain->SetBranchAddress("etaMet", etaMet, &b_etaMet);
   fChain->SetBranchAddress("phiMet", phiMet, &b_phiMet);
   fChain->SetBranchAddress("pxMet", pxMet, &b_pxMet);
   fChain->SetBranchAddress("pyMet", pyMet, &b_pyMet);
   fChain->SetBranchAddress("pzMet", pzMet, &b_pzMet);
   fChain->SetBranchAddress("vertexXMet", vertexXMet, &b_vertexXMet);
   fChain->SetBranchAddress("vertexYMet", vertexYMet, &b_vertexYMet);
   fChain->SetBranchAddress("vertexZMet", vertexZMet, &b_vertexZMet);
   fChain->SetBranchAddress("massMet", massMet, &b_massMet);
   fChain->SetBranchAddress("nTCMet", &nTCMet, &b_nTCMet);
   fChain->SetBranchAddress("chargeTCMet", chargeTCMet, &b_chargeTCMet);
   fChain->SetBranchAddress("energyTCMet", energyTCMet, &b_energyTCMet);
   fChain->SetBranchAddress("etTCMet", etTCMet, &b_etTCMet);
   fChain->SetBranchAddress("momentumTCMet", momentumTCMet, &b_momentumTCMet);
   fChain->SetBranchAddress("thetaTCMet", thetaTCMet, &b_thetaTCMet);
   fChain->SetBranchAddress("etaTCMet", etaTCMet, &b_etaTCMet);
   fChain->SetBranchAddress("phiTCMet", phiTCMet, &b_phiTCMet);
   fChain->SetBranchAddress("pxTCMet", pxTCMet, &b_pxTCMet);
   fChain->SetBranchAddress("pyTCMet", pyTCMet, &b_pyTCMet);
   fChain->SetBranchAddress("pzTCMet", pzTCMet, &b_pzTCMet);
   fChain->SetBranchAddress("vertexXTCMet", vertexXTCMet, &b_vertexXTCMet);
   fChain->SetBranchAddress("vertexYTCMet", vertexYTCMet, &b_vertexYTCMet);
   fChain->SetBranchAddress("vertexZTCMet", vertexZTCMet, &b_vertexZTCMet);
   fChain->SetBranchAddress("massTCMet", massTCMet, &b_massTCMet);
   fChain->SetBranchAddress("nPFMet", &nPFMet, &b_nPFMet);
   fChain->SetBranchAddress("chargePFMet", chargePFMet, &b_chargePFMet);
   fChain->SetBranchAddress("energyPFMet", energyPFMet, &b_energyPFMet);
   fChain->SetBranchAddress("etPFMet", etPFMet, &b_etPFMet);
   fChain->SetBranchAddress("momentumPFMet", momentumPFMet, &b_momentumPFMet);
   fChain->SetBranchAddress("thetaPFMet", thetaPFMet, &b_thetaPFMet);
   fChain->SetBranchAddress("etaPFMet", etaPFMet, &b_etaPFMet);
   fChain->SetBranchAddress("phiPFMet", phiPFMet, &b_phiPFMet);
   fChain->SetBranchAddress("pxPFMet", pxPFMet, &b_pxPFMet);
   fChain->SetBranchAddress("pyPFMet", pyPFMet, &b_pyPFMet);
   fChain->SetBranchAddress("pzPFMet", pzPFMet, &b_pzPFMet);
   fChain->SetBranchAddress("vertexXPFMet", vertexXPFMet, &b_vertexXPFMet);
   fChain->SetBranchAddress("vertexYPFMet", vertexYPFMet, &b_vertexYPFMet);
   fChain->SetBranchAddress("vertexZPFMet", vertexZPFMet, &b_vertexZPFMet);
   fChain->SetBranchAddress("massPFMet", massPFMet, &b_massPFMet);
   fChain->SetBranchAddress("nGenMet", &nGenMet, &b_nGenMet);
   fChain->SetBranchAddress("chargeGenMet", chargeGenMet, &b_chargeGenMet);
   fChain->SetBranchAddress("energyGenMet", energyGenMet, &b_energyGenMet);
   fChain->SetBranchAddress("etGenMet", etGenMet, &b_etGenMet);
   fChain->SetBranchAddress("momentumGenMet", momentumGenMet, &b_momentumGenMet);
   fChain->SetBranchAddress("thetaGenMet", thetaGenMet, &b_thetaGenMet);
   fChain->SetBranchAddress("etaGenMet", etaGenMet, &b_etaGenMet);
   fChain->SetBranchAddress("phiGenMet", phiGenMet, &b_phiGenMet);
   fChain->SetBranchAddress("pxGenMet", pxGenMet, &b_pxGenMet);
   fChain->SetBranchAddress("pyGenMet", pyGenMet, &b_pyGenMet);
   fChain->SetBranchAddress("pzGenMet", pzGenMet, &b_pzGenMet);
   fChain->SetBranchAddress("vertexXGenMet", vertexXGenMet, &b_vertexXGenMet);
   fChain->SetBranchAddress("vertexYGenMet", vertexYGenMet, &b_vertexYGenMet);
   fChain->SetBranchAddress("vertexZGenMet", vertexZGenMet, &b_vertexZGenMet);
   fChain->SetBranchAddress("massGenMet", massGenMet, &b_massGenMet);
   fChain->SetBranchAddress("nAK5CorrJet", &nAK5CorrJet, &b_nAK5CorrJet);
   fChain->SetBranchAddress("chargeAK5CorrJet", chargeAK5CorrJet, &b_chargeAK5CorrJet);
   fChain->SetBranchAddress("energyAK5CorrJet", energyAK5CorrJet, &b_energyAK5CorrJet);
   fChain->SetBranchAddress("etAK5CorrJet", etAK5CorrJet, &b_etAK5CorrJet);
   fChain->SetBranchAddress("momentumAK5CorrJet", momentumAK5CorrJet, &b_momentumAK5CorrJet);
   fChain->SetBranchAddress("thetaAK5CorrJet", thetaAK5CorrJet, &b_thetaAK5CorrJet);
   fChain->SetBranchAddress("etaAK5CorrJet", etaAK5CorrJet, &b_etaAK5CorrJet);
   fChain->SetBranchAddress("phiAK5CorrJet", phiAK5CorrJet, &b_phiAK5CorrJet);
   fChain->SetBranchAddress("pxAK5CorrJet", pxAK5CorrJet, &b_pxAK5CorrJet);
   fChain->SetBranchAddress("pyAK5CorrJet", pyAK5CorrJet, &b_pyAK5CorrJet);
   fChain->SetBranchAddress("pzAK5CorrJet", pzAK5CorrJet, &b_pzAK5CorrJet);
   fChain->SetBranchAddress("vertexXAK5CorrJet", vertexXAK5CorrJet, &b_vertexXAK5CorrJet);
   fChain->SetBranchAddress("vertexYAK5CorrJet", vertexYAK5CorrJet, &b_vertexYAK5CorrJet);
   fChain->SetBranchAddress("vertexZAK5CorrJet", vertexZAK5CorrJet, &b_vertexZAK5CorrJet);
   fChain->SetBranchAddress("massAK5CorrJet", massAK5CorrJet, &b_massAK5CorrJet);
   fChain->SetBranchAddress("alphaAK5CorrJet", alphaAK5CorrJet, &b_alphaAK5CorrJet);
   fChain->SetBranchAddress("emFracAK5CorrJet", emFracAK5CorrJet, &b_emFracAK5CorrJet);
   fChain->SetBranchAddress("hadFracAK5CorrJet", hadFracAK5CorrJet, &b_hadFracAK5CorrJet);
   fChain->SetBranchAddress("combinedSecondaryVertexBJetTagsAK5CorrJet", combinedSecondaryVertexBJetTagsAK5CorrJet, &b_combinedSecondaryVertexBJetTagsAK5CorrJet);
   fChain->SetBranchAddress("combinedSecondaryVertexMVABJetTagsAK5CorrJet", combinedSecondaryVertexMVABJetTagsAK5CorrJet, &b_combinedSecondaryVertexMVABJetTagsAK5CorrJet);
   fChain->SetBranchAddress("jetBProbabilityBJetTagsAK5CorrJet", jetBProbabilityBJetTagsAK5CorrJet, &b_jetBProbabilityBJetTagsAK5CorrJet);
   fChain->SetBranchAddress("jetProbabilityBJetTagsAK5CorrJet", jetProbabilityBJetTagsAK5CorrJet, &b_jetProbabilityBJetTagsAK5CorrJet);
   fChain->SetBranchAddress("simpleSecondaryVertexBJetTagsAK5CorrJet", simpleSecondaryVertexBJetTagsAK5CorrJet, &b_simpleSecondaryVertexBJetTagsAK5CorrJet);
   fChain->SetBranchAddress("softMuonBJetTagsAK5CorrJet", softMuonBJetTagsAK5CorrJet, &b_softMuonBJetTagsAK5CorrJet);
   fChain->SetBranchAddress("trackCountingHighPurBJetTagsAK5CorrJet", trackCountingHighPurBJetTagsAK5CorrJet, &b_trackCountingHighPurBJetTagsAK5CorrJet);
   fChain->SetBranchAddress("trackCountingHighEffBJetTagsAK5CorrJet", trackCountingHighEffBJetTagsAK5CorrJet, &b_trackCountingHighEffBJetTagsAK5CorrJet);
   fChain->SetBranchAddress("nAK5Jet", &nAK5Jet, &b_nAK5Jet);
   fChain->SetBranchAddress("chargeAK5Jet", chargeAK5Jet, &b_chargeAK5Jet);
   fChain->SetBranchAddress("energyAK5Jet", energyAK5Jet, &b_energyAK5Jet);
   fChain->SetBranchAddress("etAK5Jet", etAK5Jet, &b_etAK5Jet);
   fChain->SetBranchAddress("momentumAK5Jet", momentumAK5Jet, &b_momentumAK5Jet);
   fChain->SetBranchAddress("thetaAK5Jet", thetaAK5Jet, &b_thetaAK5Jet);
   fChain->SetBranchAddress("etaAK5Jet", etaAK5Jet, &b_etaAK5Jet);
   fChain->SetBranchAddress("phiAK5Jet", phiAK5Jet, &b_phiAK5Jet);
   fChain->SetBranchAddress("pxAK5Jet", pxAK5Jet, &b_pxAK5Jet);
   fChain->SetBranchAddress("pyAK5Jet", pyAK5Jet, &b_pyAK5Jet);
   fChain->SetBranchAddress("pzAK5Jet", pzAK5Jet, &b_pzAK5Jet);
   fChain->SetBranchAddress("vertexXAK5Jet", vertexXAK5Jet, &b_vertexXAK5Jet);
   fChain->SetBranchAddress("vertexYAK5Jet", vertexYAK5Jet, &b_vertexYAK5Jet);
   fChain->SetBranchAddress("vertexZAK5Jet", vertexZAK5Jet, &b_vertexZAK5Jet);
   fChain->SetBranchAddress("massAK5Jet", massAK5Jet, &b_massAK5Jet);
   fChain->SetBranchAddress("alphaAK5Jet", alphaAK5Jet, &b_alphaAK5Jet);
   fChain->SetBranchAddress("emFracAK5Jet", emFracAK5Jet, &b_emFracAK5Jet);
   fChain->SetBranchAddress("hadFracAK5Jet", hadFracAK5Jet, &b_hadFracAK5Jet);
   fChain->SetBranchAddress("nAK5PFCorrJet", &nAK5PFCorrJet, &b_nAK5PFCorrJet);
   fChain->SetBranchAddress("chargeAK5PFCorrJet", chargeAK5PFCorrJet, &b_chargeAK5PFCorrJet);
   fChain->SetBranchAddress("energyAK5PFCorrJet", energyAK5PFCorrJet, &b_energyAK5PFCorrJet);
   fChain->SetBranchAddress("etAK5PFCorrJet", etAK5PFCorrJet, &b_etAK5PFCorrJet);
   fChain->SetBranchAddress("momentumAK5PFCorrJet", momentumAK5PFCorrJet, &b_momentumAK5PFCorrJet);
   fChain->SetBranchAddress("thetaAK5PFCorrJet", thetaAK5PFCorrJet, &b_thetaAK5PFCorrJet);
   fChain->SetBranchAddress("etaAK5PFCorrJet", etaAK5PFCorrJet, &b_etaAK5PFCorrJet);
   fChain->SetBranchAddress("phiAK5PFCorrJet", phiAK5PFCorrJet, &b_phiAK5PFCorrJet);
   fChain->SetBranchAddress("pxAK5PFCorrJet", pxAK5PFCorrJet, &b_pxAK5PFCorrJet);
   fChain->SetBranchAddress("pyAK5PFCorrJet", pyAK5PFCorrJet, &b_pyAK5PFCorrJet);
   fChain->SetBranchAddress("pzAK5PFCorrJet", pzAK5PFCorrJet, &b_pzAK5PFCorrJet);
   fChain->SetBranchAddress("vertexXAK5PFCorrJet", vertexXAK5PFCorrJet, &b_vertexXAK5PFCorrJet);
   fChain->SetBranchAddress("vertexYAK5PFCorrJet", vertexYAK5PFCorrJet, &b_vertexYAK5PFCorrJet);
   fChain->SetBranchAddress("vertexZAK5PFCorrJet", vertexZAK5PFCorrJet, &b_vertexZAK5PFCorrJet);
   fChain->SetBranchAddress("massAK5PFCorrJet", massAK5PFCorrJet, &b_massAK5PFCorrJet);
   fChain->SetBranchAddress("chargedHadronEnergyAK5PFCorrJet", chargedHadronEnergyAK5PFCorrJet, &b_chargedHadronEnergyAK5PFCorrJet);
   fChain->SetBranchAddress("neutralHadronEnergyAK5PFCorrJet", neutralHadronEnergyAK5PFCorrJet, &b_neutralHadronEnergyAK5PFCorrJet);
   fChain->SetBranchAddress("chargedEmEnergyAK5PFCorrJet", chargedEmEnergyAK5PFCorrJet, &b_chargedEmEnergyAK5PFCorrJet);
   fChain->SetBranchAddress("neutralEmEnergyAK5PFCorrJet", neutralEmEnergyAK5PFCorrJet, &b_neutralEmEnergyAK5PFCorrJet);
   fChain->SetBranchAddress("neutralMultiplicityAK5PFCorrJet", neutralMultiplicityAK5PFCorrJet, &b_neutralMultiplicityAK5PFCorrJet);
   fChain->SetBranchAddress("chargedMultiplicityAK5PFCorrJet", chargedMultiplicityAK5PFCorrJet, &b_chargedMultiplicityAK5PFCorrJet);
   fChain->SetBranchAddress("muonMultiplicityAK5PFCorrJet", muonMultiplicityAK5PFCorrJet, &b_muonMultiplicityAK5PFCorrJet);
   fChain->SetBranchAddress("nAK5PFJet", &nAK5PFJet, &b_nAK5PFJet);
   fChain->SetBranchAddress("chargeAK5PFJet", chargeAK5PFJet, &b_chargeAK5PFJet);
   fChain->SetBranchAddress("energyAK5PFJet", energyAK5PFJet, &b_energyAK5PFJet);
   fChain->SetBranchAddress("etAK5PFJet", etAK5PFJet, &b_etAK5PFJet);
   fChain->SetBranchAddress("momentumAK5PFJet", momentumAK5PFJet, &b_momentumAK5PFJet);
   fChain->SetBranchAddress("thetaAK5PFJet", thetaAK5PFJet, &b_thetaAK5PFJet);
   fChain->SetBranchAddress("etaAK5PFJet", etaAK5PFJet, &b_etaAK5PFJet);
   fChain->SetBranchAddress("phiAK5PFJet", phiAK5PFJet, &b_phiAK5PFJet);
   fChain->SetBranchAddress("pxAK5PFJet", pxAK5PFJet, &b_pxAK5PFJet);
   fChain->SetBranchAddress("pyAK5PFJet", pyAK5PFJet, &b_pyAK5PFJet);
   fChain->SetBranchAddress("pzAK5PFJet", pzAK5PFJet, &b_pzAK5PFJet);
   fChain->SetBranchAddress("vertexXAK5PFJet", vertexXAK5PFJet, &b_vertexXAK5PFJet);
   fChain->SetBranchAddress("vertexYAK5PFJet", vertexYAK5PFJet, &b_vertexYAK5PFJet);
   fChain->SetBranchAddress("vertexZAK5PFJet", vertexZAK5PFJet, &b_vertexZAK5PFJet);
   fChain->SetBranchAddress("massAK5PFJet", massAK5PFJet, &b_massAK5PFJet);
   fChain->SetBranchAddress("chargedHadronEnergyAK5PFJet", chargedHadronEnergyAK5PFJet, &b_chargedHadronEnergyAK5PFJet);
   fChain->SetBranchAddress("neutralHadronEnergyAK5PFJet", neutralHadronEnergyAK5PFJet, &b_neutralHadronEnergyAK5PFJet);
   fChain->SetBranchAddress("chargedEmEnergyAK5PFJet", chargedEmEnergyAK5PFJet, &b_chargedEmEnergyAK5PFJet);
   fChain->SetBranchAddress("neutralEmEnergyAK5PFJet", neutralEmEnergyAK5PFJet, &b_neutralEmEnergyAK5PFJet);
   fChain->SetBranchAddress("neutralMultiplicityAK5PFJet", neutralMultiplicityAK5PFJet, &b_neutralMultiplicityAK5PFJet);
   fChain->SetBranchAddress("chargedMultiplicityAK5PFJet", chargedMultiplicityAK5PFJet, &b_chargedMultiplicityAK5PFJet);
   fChain->SetBranchAddress("muonMultiplicityAK5PFJet", muonMultiplicityAK5PFJet, &b_muonMultiplicityAK5PFJet);
   fChain->SetBranchAddress("nAK5GenJet", &nAK5GenJet, &b_nAK5GenJet);
   fChain->SetBranchAddress("chargeAK5GenJet", chargeAK5GenJet, &b_chargeAK5GenJet);
   fChain->SetBranchAddress("energyAK5GenJet", energyAK5GenJet, &b_energyAK5GenJet);
   fChain->SetBranchAddress("etAK5GenJet", etAK5GenJet, &b_etAK5GenJet);
   fChain->SetBranchAddress("momentumAK5GenJet", momentumAK5GenJet, &b_momentumAK5GenJet);
   fChain->SetBranchAddress("thetaAK5GenJet", thetaAK5GenJet, &b_thetaAK5GenJet);
   fChain->SetBranchAddress("etaAK5GenJet", etaAK5GenJet, &b_etaAK5GenJet);
   fChain->SetBranchAddress("phiAK5GenJet", phiAK5GenJet, &b_phiAK5GenJet);
   fChain->SetBranchAddress("pxAK5GenJet", pxAK5GenJet, &b_pxAK5GenJet);
   fChain->SetBranchAddress("pyAK5GenJet", pyAK5GenJet, &b_pyAK5GenJet);
   fChain->SetBranchAddress("pzAK5GenJet", pzAK5GenJet, &b_pzAK5GenJet);
   fChain->SetBranchAddress("vertexXAK5GenJet", vertexXAK5GenJet, &b_vertexXAK5GenJet);
   fChain->SetBranchAddress("vertexYAK5GenJet", vertexYAK5GenJet, &b_vertexYAK5GenJet);
   fChain->SetBranchAddress("vertexZAK5GenJet", vertexZAK5GenJet, &b_vertexZAK5GenJet);
   fChain->SetBranchAddress("massAK5GenJet", massAK5GenJet, &b_massAK5GenJet);
   fChain->SetBranchAddress("genPtHat", &genPtHat, &b_genPtHat);
   fChain->SetBranchAddress("genProcessId", &genProcessId, &b_genProcessId);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genAlphaQCD", &genAlphaQCD, &b_genAlphaQCD);
   fChain->SetBranchAddress("genAlphaQED", &genAlphaQED, &b_genAlphaQED);
 
  
  
}
