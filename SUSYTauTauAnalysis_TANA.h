//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar  2 15:06:52 2010 by ROOT version 5.22/00d
// from TTree ntp1/ntp1
// found on file: rfio:///castor/cern.ch/user/e/emanuele/CMST3/Vecbos2010/MC/7TeV/SUSY/LM1/default_MC_16.root
//////////////////////////////////////////////////////////

#ifndef SUSYTauTauAnalysis_h
#define SUSYTauTauAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class SUSYTauTauAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nl1Technical;
   Int_t           l1Technical[64];   //[nl1Technical]
   Int_t           nl1Global;
   Int_t           l1Global[128];   //[nl1Global]
   Int_t           runNumber;
   Int_t           eventNumber;
   Int_t           lumiBlock;
   Int_t           bunchCrossing;
   Int_t           orbitNumber;
   Int_t           nMc;
   Float_t         pMc[1001];   //[nMc]
   Float_t         massMc[1001];   //[nMc]
   Float_t         thetaMc[1001];   //[nMc]
   Float_t         etaMc[1001];   //[nMc]
   Float_t         phiMc[1001];   //[nMc]
   Float_t         energyMc[1001];   //[nMc]
   Int_t           idMc[1001];   //[nMc]
   Int_t           mothMc[1001];   //[nMc]
   Int_t           nDauMc[1001];   //[nMc]
   Int_t           statusMc[1001];   //[nMc]
   Float_t         xMc[1001];   //[nMc]
   Float_t         yMc[1001];   //[nMc]
   Float_t         zMc[1001];   //[nMc]
   Int_t           nTrg;
   UChar_t         firedTrg[100];   //[nTrg]
   Int_t           nEle;
   Int_t           chargeEle[12];   //[nEle]
   Float_t         energyEle[12];   //[nEle]
   Float_t         etEle[12];   //[nEle]
   Float_t         momentumEle[12];   //[nEle]
   Float_t         thetaEle[12];   //[nEle]
   Float_t         etaEle[12];   //[nEle]
   Float_t         phiEle[12];   //[nEle]
   Float_t         pxEle[12];   //[nEle]
   Float_t         pyEle[12];   //[nEle]
   Float_t         pzEle[12];   //[nEle]
   Float_t         vertexXEle[12];   //[nEle]
   Float_t         vertexYEle[12];   //[nEle]
   Float_t         vertexZEle[12];   //[nEle]
   Float_t         massEle[12];   //[nEle]
   Int_t           fiducialFlagsEle[12];   //[nEle]
   Int_t           recoFlagsEle[12];   //[nEle]
   Int_t           energyCorrectionsEle[12];   //[nEle]
   Float_t         esEnergyEle[12];   //[nEle]
   Int_t           superClusterIndexEle[12];   //[nEle]
   Int_t           PFsuperClusterIndexEle[12];   //[nEle]
   Int_t           trackIndexEle[12];   //[nEle]
   Int_t           gsfTrackIndexEle[12];   //[nEle]
   Int_t           classificationEle[12];   //[nEle]
   Int_t           standardClassificationEle[12];   //[nEle]
   Float_t         hOverEEle[12];   //[nEle]
   Float_t         eSuperClusterOverPEle[12];   //[nEle]
   Float_t         eSeedOverPoutEle[12];   //[nEle]
   Float_t         deltaEtaAtVtxEle[12];   //[nEle]
   Float_t         deltaPhiAtVtxEle[12];   //[nEle]
   Float_t         deltaEtaAtCaloEle[12];   //[nEle]
   Float_t         deltaPhiAtCaloEle[12];   //[nEle]
   Float_t         tipEle[12];   //[nEle]
   Float_t         dr03TkSumPtEle[12];   //[nEle]
   Float_t         dr03EcalRecHitSumEtEle[12];   //[nEle]
   Float_t         dr03HcalTowerSumEtEle[12];   //[nEle]
   Float_t         dr04TkSumPtEle[12];   //[nEle]
   Float_t         dr04EcalRecHitSumEtEle[12];   //[nEle]
   Float_t         dr04HcalTowerSumEtEle[12];   //[nEle]
   Float_t         scBasedEcalSum03Ele[12];   //[nEle]
   Float_t         scBasedEcalSum04Ele[12];   //[nEle]
   Float_t         scHaloBasedEcalSum03Ele[12];   //[nEle]
   Float_t         scHaloBasedEcalSum04Ele[12];   //[nEle]
   Int_t           eleIdCutsEle[12];   //[nEle]
   Float_t         eleIdLikelihoodEle[12];   //[nEle]
   UChar_t         pflowMVAEle[12];   //[nEle]
   Int_t           nPFEle;
   Int_t           chargePFEle[24];   //[nPFEle]
   Float_t         energyPFEle[24];   //[nPFEle]
   Float_t         etPFEle[24];   //[nPFEle]
   Float_t         momentumPFEle[24];   //[nPFEle]
   Float_t         thetaPFEle[24];   //[nPFEle]
   Float_t         etaPFEle[24];   //[nPFEle]
   Float_t         phiPFEle[24];   //[nPFEle]
   Float_t         pxPFEle[24];   //[nPFEle]
   Float_t         pyPFEle[24];   //[nPFEle]
   Float_t         pzPFEle[24];   //[nPFEle]
   Float_t         vertexXPFEle[24];   //[nPFEle]
   Float_t         vertexYPFEle[24];   //[nPFEle]
   Float_t         vertexZPFEle[24];   //[nPFEle]
   Float_t         massPFEle[24];   //[nPFEle]
   Float_t         MvaOutputPFEle[24];   //[nPFEle]
   Float_t         PS1EnergyPFEle[24];   //[nPFEle]
   Float_t         PS2EnergyPFEle[24];   //[nPFEle]
   Float_t         EcalEnergyPFEle[24];   //[nPFEle]
   Float_t         EcalElectronEnergyPFEle[24];   //[nPFEle]
   Int_t           gsfTrackIndexPFEle[24];   //[nPFEle]
   Int_t           trackIndexPFEle[24];   //[nPFEle]
   Int_t           nSC;
   Int_t           nBCSC[42];   //[nSC]
   Int_t           nCrystalsSC[42];   //[nSC]
   Int_t           iAlgoSC[42];   //[nSC]
   Float_t         rawEnergySC[42];   //[nSC]
   Float_t         energySC[42];   //[nSC]
   Float_t         etaSC[42];   //[nSC]
   Float_t         thetaSC[42];   //[nSC]
   Float_t         phiSC[42];   //[nSC]
   Float_t         e3x3SC[42];   //[nSC]
   Float_t         e5x5SC[42];   //[nSC]
   Float_t         eMaxSC[42];   //[nSC]
   Float_t         e2x2SC[42];   //[nSC]
   Float_t         e2ndSC[42];   //[nSC]
   Float_t         covIEtaIEtaSC[42];   //[nSC]
   Float_t         covIEtaIPhiSC[42];   //[nSC]
   Float_t         covIPhiIPhiSC[42];   //[nSC]
   Float_t         hOverESC[42];   //[nSC]
   Int_t           recoFlagSC[42];   //[nSC]
   Int_t           channelStatusSC[42];   //[nSC]
   Float_t         timeSC[42];   //[nSC]
   Float_t         chi2ProbSC[42];   //[nSC]
   Float_t         seedEnergySC[42];   //[nSC]
   Int_t           idClosProblSC[42];   //[nSC]
   Int_t           sevClosProblSC[42];   //[nSC]
   Float_t         fracClosProblSC[42];   //[nSC]
   Int_t           trackIndexSC[42];   //[nSC]
   Float_t         trackDeltaRSC[42];   //[nSC]
   Float_t         trackDeltaPhiSC[42];   //[nSC]
   Float_t         trackDeltaEtaSC[42];   //[nSC]
   Int_t           gsfTrackIndexSC[42];   //[nSC]
   Float_t         gsfTrackDeltaRSC[42];   //[nSC]
   Float_t         gsfTrackDeltaPhiSC[42];   //[nSC]
   Float_t         gsfTrackDeltaEtaSC[42];   //[nSC]
   Float_t         pxVtxPropagatedNegChargeSC[42];   //[nSC]
   Float_t         pyVtxPropagatedNegChargeSC[42];   //[nSC]
   Float_t         pzVtxPropagatedNegChargeSC[42];   //[nSC]
   Float_t         pxVtxPropagatedPosChargeSC[42];   //[nSC]
   Float_t         pyVtxPropagatedPosChargeSC[42];   //[nSC]
   Float_t         pzVtxPropagatedPosChargeSC[42];   //[nSC]
   Int_t           nPFSC;
   Int_t           nBCPFSC[8];   //[nPFSC]
   Int_t           nCrystalsPFSC[8];   //[nPFSC]
   Int_t           iAlgoPFSC[8];   //[nPFSC]
   Float_t         rawEnergyPFSC[8];   //[nPFSC]
   Float_t         energyPFSC[8];   //[nPFSC]
   Float_t         etaPFSC[8];   //[nPFSC]
   Float_t         thetaPFSC[8];   //[nPFSC]
   Float_t         phiPFSC[8];   //[nPFSC]
   Float_t         e3x3PFSC[8];   //[nPFSC]
   Float_t         e5x5PFSC[8];   //[nPFSC]
   Float_t         eMaxPFSC[8];   //[nPFSC]
   Float_t         e2x2PFSC[8];   //[nPFSC]
   Float_t         e2ndPFSC[8];   //[nPFSC]
   Float_t         covIEtaIEtaPFSC[8];   //[nPFSC]
   Float_t         covIEtaIPhiPFSC[8];   //[nPFSC]
   Float_t         covIPhiIPhiPFSC[8];   //[nPFSC]
   Float_t         hOverEPFSC[8];   //[nPFSC]
   Int_t           recoFlagPFSC[8];   //[nPFSC]
   Int_t           channelStatusPFSC[8];   //[nPFSC]
   Float_t         timePFSC[8];   //[nPFSC]
   Float_t         chi2ProbPFSC[8];   //[nPFSC]
   Float_t         seedEnergyPFSC[8];   //[nPFSC]
   Int_t           idClosProblPFSC[8];   //[nPFSC]
   Int_t           sevClosProblPFSC[8];   //[nPFSC]
   Float_t         fracClosProblPFSC[8];   //[nPFSC]
   Float_t         pxVtxPropagatedNegChargePFSC[8];   //[nPFSC]
   Float_t         pyVtxPropagatedNegChargePFSC[8];   //[nPFSC]
   Float_t         pzVtxPropagatedNegChargePFSC[8];   //[nPFSC]
   Float_t         pxVtxPropagatedPosChargePFSC[8];   //[nPFSC]
   Float_t         pyVtxPropagatedPosChargePFSC[8];   //[nPFSC]
   Float_t         pzVtxPropagatedPosChargePFSC[8];   //[nPFSC]
   Int_t           nTrack;
   Float_t         pxTrack[413];   //[nTrack]
   Float_t         pyTrack[413];   //[nTrack]
   Float_t         pzTrack[413];   //[nTrack]
   Int_t           vtxIndexTrack[413];   //[nTrack]
   Float_t         vtxWeightTrack[413];   //[nTrack]
   Float_t         chargeTrack[413];   //[nTrack]
   Float_t         ptErrorTrack[413];   //[nTrack]
   Float_t         trackValidHitsTrack[413];   //[nTrack]
   Float_t         trackLostHitsTrack[413];   //[nTrack]
   Float_t         trackNormalizedChi2Track[413];   //[nTrack]
   Int_t           qualityMaskTrack[413];   //[nTrack]
   Float_t         trackDxyTrack[413];   //[nTrack]
   Float_t         trackD0Track[413];   //[nTrack]
   Float_t         trackDszTrack[413];   //[nTrack]
   Float_t         trackDzTrack[413];   //[nTrack]
   Float_t         trackDxyErrorTrack[413];   //[nTrack]
   Float_t         trackD0ErrorTrack[413];   //[nTrack]
   Float_t         trackDszErrorTrack[413];   //[nTrack]
   Float_t         trackDzErrorTrack[413];   //[nTrack]
   Float_t         trackDxyPVTrack[413];   //[nTrack]
   Float_t         trackDszPVTrack[413];   //[nTrack]
   Float_t         trackDzPVTrack[413];   //[nTrack]
   Float_t         trackVxTrack[413];   //[nTrack]
   Float_t         trackVyTrack[413];   //[nTrack]
   Float_t         trackVzTrack[413];   //[nTrack]
   Float_t         pxAtOuterTrack[413];   //[nTrack]
   Float_t         pyAtOuterTrack[413];   //[nTrack]
   Float_t         pzAtOuterTrack[413];   //[nTrack]
   Float_t         xAtOuterTrack[413];   //[nTrack]
   Float_t         yAtOuterTrack[413];   //[nTrack]
   Float_t         zAtOuterTrack[413];   //[nTrack]
   Float_t         pxAtInnerTrack[413];   //[nTrack]
   Float_t         pyAtInnerTrack[413];   //[nTrack]
   Float_t         pzAtInnerTrack[413];   //[nTrack]
   Float_t         xAtInnerTrack[413];   //[nTrack]
   Float_t         yAtInnerTrack[413];   //[nTrack]
   Float_t         zAtInnerTrack[413];   //[nTrack]
   Float_t         recHitsSizeTrack[413];   //[nTrack]
   UChar_t         isPixB1Track[413];   //[nTrack]
   UChar_t         isPixB2Track[413];   //[nTrack]
   UChar_t         isPixE1Track[413];   //[nTrack]
   UChar_t         isPixE2Track[413];   //[nTrack]
   Int_t           numberOfValidPixelBarrelHitsTrack[413];   //[nTrack]
   Int_t           numberOfValidPixelEndcapHitsTrack[413];   //[nTrack]
   Int_t           numberOfValidStripTIBHitsTrack[413];   //[nTrack]
   Int_t           numberOfValidStripTIDHitsTrack[413];   //[nTrack]
   Int_t           numberOfValidStripTOBHitsTrack[413];   //[nTrack]
   Int_t           numberOfValidStripTECHitsTrack[413];   //[nTrack]
   Float_t         truncatedDeDxTrack[413];   //[nTrack]
   Float_t         truncatedDeDxErrorTrack[413];   //[nTrack]
   Float_t         truncatedDeDxNoMTrack[413];   //[nTrack]
   Float_t         medianDeDxTrack[413];   //[nTrack]
   Float_t         medianDeDxErrorTrack[413];   //[nTrack]
   Float_t         medianDeDxNoMTrack[413];   //[nTrack]
   Float_t         harmonic2DeDxTrack[413];   //[nTrack]
   Float_t         harmonic2DeDxErrorTrack[413];   //[nTrack]
   Float_t         harmonic2DeDxNoMTrack[413];   //[nTrack]
   Int_t           nGsfTrack;
   Float_t         pxGsfTrack[36];   //[nGsfTrack]
   Float_t         pyGsfTrack[36];   //[nGsfTrack]
   Float_t         pzGsfTrack[36];   //[nGsfTrack]
   Float_t         chargeGsfTrack[36];   //[nGsfTrack]
   Float_t         ptErrorGsfTrack[36];   //[nGsfTrack]
   Float_t         trackValidHitsGsfTrack[36];   //[nGsfTrack]
   Float_t         trackLostHitsGsfTrack[36];   //[nGsfTrack]
   Float_t         trackNormalizedChi2GsfTrack[36];   //[nGsfTrack]
   Int_t           qualityMaskGsfTrack[36];   //[nGsfTrack]
   Float_t         trackDxyGsfTrack[36];   //[nGsfTrack]
   Float_t         trackD0GsfTrack[36];   //[nGsfTrack]
   Float_t         trackDszGsfTrack[36];   //[nGsfTrack]
   Float_t         trackDzGsfTrack[36];   //[nGsfTrack]
   Float_t         trackDxyErrorGsfTrack[36];   //[nGsfTrack]
   Float_t         trackD0ErrorGsfTrack[36];   //[nGsfTrack]
   Float_t         trackDszErrorGsfTrack[36];   //[nGsfTrack]
   Float_t         trackDzErrorGsfTrack[36];   //[nGsfTrack]
   Float_t         trackDxyPVGsfTrack[36];   //[nGsfTrack]
   Float_t         trackDszPVGsfTrack[36];   //[nGsfTrack]
   Float_t         trackDzPVGsfTrack[36];   //[nGsfTrack]
   Float_t         trackVxGsfTrack[36];   //[nGsfTrack]
   Float_t         trackVyGsfTrack[36];   //[nGsfTrack]
   Float_t         trackVzGsfTrack[36];   //[nGsfTrack]
   Float_t         pxAtOuterGsfTrack[36];   //[nGsfTrack]
   Float_t         pyAtOuterGsfTrack[36];   //[nGsfTrack]
   Float_t         pzAtOuterGsfTrack[36];   //[nGsfTrack]
   Float_t         xAtOuterGsfTrack[36];   //[nGsfTrack]
   Float_t         yAtOuterGsfTrack[36];   //[nGsfTrack]
   Float_t         zAtOuterGsfTrack[36];   //[nGsfTrack]
   Float_t         pxAtInnerGsfTrack[36];   //[nGsfTrack]
   Float_t         pyAtInnerGsfTrack[36];   //[nGsfTrack]
   Float_t         pzAtInnerGsfTrack[36];   //[nGsfTrack]
   Float_t         xAtInnerGsfTrack[36];   //[nGsfTrack]
   Float_t         yAtInnerGsfTrack[36];   //[nGsfTrack]
   Float_t         zAtInnerGsfTrack[36];   //[nGsfTrack]
   Float_t         recHitsSizeGsfTrack[36];   //[nGsfTrack]
   UChar_t         isPixB1GsfTrack[36];   //[nGsfTrack]
   UChar_t         isPixB2GsfTrack[36];   //[nGsfTrack]
   UChar_t         isPixE1GsfTrack[36];   //[nGsfTrack]
   UChar_t         isPixE2GsfTrack[36];   //[nGsfTrack]
   Int_t           numberOfValidPixelBarrelHitsGsfTrack[36];   //[nGsfTrack]
   Int_t           numberOfValidPixelEndcapHitsGsfTrack[36];   //[nGsfTrack]
   Int_t           numberOfValidStripTIBHitsGsfTrack[36];   //[nGsfTrack]
   Int_t           numberOfValidStripTIDHitsGsfTrack[36];   //[nGsfTrack]
   Int_t           numberOfValidStripTOBHitsGsfTrack[36];   //[nGsfTrack]
   Int_t           numberOfValidStripTECHitsGsfTrack[36];   //[nGsfTrack]
   Int_t           chargeModeGsfTrack[36];   //[nGsfTrack]
   Float_t         pxModeGsfTrack[36];   //[nGsfTrack]
   Float_t         pyModeGsfTrack[36];   //[nGsfTrack]
   Float_t         pzModeGsfTrack[36];   //[nGsfTrack]
   Int_t           recoFlagsGsfTrack[36];   //[nGsfTrack]
   Int_t           nGlobalMuonTrack;
   Float_t         pxGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         pyGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         pzGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         chargeGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         ptErrorGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackValidHitsGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackLostHitsGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackNormalizedChi2GlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Int_t           qualityMaskGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackDxyGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackD0GlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackDszGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackDzGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackDxyErrorGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackD0ErrorGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackDszErrorGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackDzErrorGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackDxyPVGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackDszPVGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackDzPVGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackVxGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackVyGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         trackVzGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         pxAtOuterGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         pyAtOuterGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         pzAtOuterGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         xAtOuterGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         yAtOuterGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         zAtOuterGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         pxAtInnerGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         pyAtInnerGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         pzAtInnerGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         xAtInnerGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         yAtInnerGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         zAtInnerGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Float_t         recHitsSizeGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   UChar_t         isPixB1GlobalMuonTrack[7];   //[nGlobalMuonTrack]
   UChar_t         isPixB2GlobalMuonTrack[7];   //[nGlobalMuonTrack]
   UChar_t         isPixE1GlobalMuonTrack[7];   //[nGlobalMuonTrack]
   UChar_t         isPixE2GlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Int_t           numberOfValidPixelBarrelHitsGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Int_t           numberOfValidPixelEndcapHitsGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Int_t           numberOfValidStripTIBHitsGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Int_t           numberOfValidStripTIDHitsGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Int_t           numberOfValidStripTOBHitsGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Int_t           numberOfValidStripTECHitsGlobalMuonTrack[7];   //[nGlobalMuonTrack]
   Int_t           nSTAMuonTrack;
   Float_t         pxSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         pySTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         pzSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         chargeSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         ptErrorSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackValidHitsSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackLostHitsSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackNormalizedChi2STAMuonTrack[13];   //[nSTAMuonTrack]
   Int_t           qualityMaskSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackDxySTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackD0STAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackDszSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackDzSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackDxyErrorSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackD0ErrorSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackDszErrorSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackDzErrorSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackDxyPVSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackDszPVSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackDzPVSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackVxSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackVySTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         trackVzSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         pxAtOuterSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         pyAtOuterSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         pzAtOuterSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         xAtOuterSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         yAtOuterSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         zAtOuterSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         pxAtInnerSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         pyAtInnerSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         pzAtInnerSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         xAtInnerSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         yAtInnerSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         zAtInnerSTAMuonTrack[13];   //[nSTAMuonTrack]
   Float_t         recHitsSizeSTAMuonTrack[13];   //[nSTAMuonTrack]
   UChar_t         isPixB1STAMuonTrack[13];   //[nSTAMuonTrack]
   UChar_t         isPixB2STAMuonTrack[13];   //[nSTAMuonTrack]
   UChar_t         isPixE1STAMuonTrack[13];   //[nSTAMuonTrack]
   UChar_t         isPixE2STAMuonTrack[13];   //[nSTAMuonTrack]
   Int_t           numberOfValidPixelBarrelHitsSTAMuonTrack[13];   //[nSTAMuonTrack]
   Int_t           numberOfValidPixelEndcapHitsSTAMuonTrack[13];   //[nSTAMuonTrack]
   Int_t           numberOfValidStripTIBHitsSTAMuonTrack[13];   //[nSTAMuonTrack]
   Int_t           numberOfValidStripTIDHitsSTAMuonTrack[13];   //[nSTAMuonTrack]
   Int_t           numberOfValidStripTOBHitsSTAMuonTrack[13];   //[nSTAMuonTrack]
   Int_t           numberOfValidStripTECHitsSTAMuonTrack[13];   //[nSTAMuonTrack]
   Int_t           nPV;
   Float_t         PVxPV[3];   //[nPV]
   Float_t         PVyPV[3];   //[nPV]
   Float_t         PVzPV[3];   //[nPV]
   Float_t         PVErrxPV[3];   //[nPV]
   Float_t         PVErryPV[3];   //[nPV]
   Float_t         PVErrzPV[3];   //[nPV]
   Float_t         SumPtPV[3];   //[nPV]
   Float_t         ndofPV[3];   //[nPV]
   Float_t         chi2PV[3];   //[nPV]
   Int_t           nMuon;
   Int_t           chargeMuon[57];   //[nMuon]
   Float_t         energyMuon[57];   //[nMuon]
   Float_t         etMuon[57];   //[nMuon]
   Float_t         momentumMuon[57];   //[nMuon]
   Float_t         thetaMuon[57];   //[nMuon]
   Float_t         etaMuon[57];   //[nMuon]
   Float_t         phiMuon[57];   //[nMuon]
   Float_t         pxMuon[57];   //[nMuon]
   Float_t         pyMuon[57];   //[nMuon]
   Float_t         pzMuon[57];   //[nMuon]
   Float_t         vertexXMuon[57];   //[nMuon]
   Float_t         vertexYMuon[57];   //[nMuon]
   Float_t         vertexZMuon[57];   //[nMuon]
   Float_t         massMuon[57];   //[nMuon]
   Int_t           trackIndexMuon[57];   //[nMuon]
   Int_t           standAloneTrackIndexMuon[57];   //[nMuon]
   Int_t           combinedTrackIndexMuon[57];   //[nMuon]
   Int_t           muonIdMuon[57];   //[nMuon]
   Float_t         sumPt03Muon[57];   //[nMuon]
   Float_t         emEt03Muon[57];   //[nMuon]
   Float_t         hadEt03Muon[57];   //[nMuon]
   Float_t         hoEt03Muon[57];   //[nMuon]
   Float_t         nTrk03Muon[57];   //[nMuon]
   Float_t         nJets03Muon[57];   //[nMuon]
   Float_t         sumPt05Muon[57];   //[nMuon]
   Float_t         emEt05Muon[57];   //[nMuon]
   Float_t         hadEt05Muon[57];   //[nMuon]
   Float_t         hoEt05Muon[57];   //[nMuon]
   Float_t         nTrk05Muon[57];   //[nMuon]
   Float_t         nJets05Muon[57];   //[nMuon]
   Float_t         EcalExpDepoMuon[57];   //[nMuon]
   Float_t         HcalExpDepoMuon[57];   //[nMuon]
   Float_t         HoExpDepoMuon[57];   //[nMuon]
   Float_t         emS9Muon[57];   //[nMuon]
   Float_t         hadS9Muon[57];   //[nMuon]
   Float_t         hoS9Muon[57];   //[nMuon]
   Float_t         CaloCompMuon[57];   //[nMuon]
   Int_t           nCaloTowers;
   Float_t         energyCaloTowers[8455];   //[nCaloTowers]
   Float_t         xCaloTowers[8455];   //[nCaloTowers]
   Float_t         yCaloTowers[8455];   //[nCaloTowers]
   Float_t         zCaloTowers[8455];   //[nCaloTowers]
   Int_t           CALOCaloTowers[8455];   //[nCaloTowers]
   Int_t           CaloIndexCaloTowers[8455];   //[nCaloTowers]
   Int_t           nMet;
   Int_t           chargeMet[1];   //[nMet]
   Float_t         energyMet[1];   //[nMet]
   Float_t         etMet[1];   //[nMet]
   Float_t         momentumMet[1];   //[nMet]
   Float_t         thetaMet[1];   //[nMet]
   Float_t         etaMet[1];   //[nMet]
   Float_t         phiMet[1];   //[nMet]
   Float_t         pxMet[1];   //[nMet]
   Float_t         pyMet[1];   //[nMet]
   Float_t         pzMet[1];   //[nMet]
   Float_t         vertexXMet[1];   //[nMet]
   Float_t         vertexYMet[1];   //[nMet]
   Float_t         vertexZMet[1];   //[nMet]
   Float_t         massMet[1];   //[nMet]
   Int_t           nTCMet;
   Int_t           chargeTCMet[1];   //[nTCMet]
   Float_t         energyTCMet[1];   //[nTCMet]
   Float_t         etTCMet[1];   //[nTCMet]
   Float_t         momentumTCMet[1];   //[nTCMet]
   Float_t         thetaTCMet[1];   //[nTCMet]
   Float_t         etaTCMet[1];   //[nTCMet]
   Float_t         phiTCMet[1];   //[nTCMet]
   Float_t         pxTCMet[1];   //[nTCMet]
   Float_t         pyTCMet[1];   //[nTCMet]
   Float_t         pzTCMet[1];   //[nTCMet]
   Float_t         vertexXTCMet[1];   //[nTCMet]
   Float_t         vertexYTCMet[1];   //[nTCMet]
   Float_t         vertexZTCMet[1];   //[nTCMet]
   Float_t         massTCMet[1];   //[nTCMet]
   Int_t           nPFMet;
   Int_t           chargePFMet[1];   //[nPFMet]
   Float_t         energyPFMet[1];   //[nPFMet]
   Float_t         etPFMet[1];   //[nPFMet]
   Float_t         momentumPFMet[1];   //[nPFMet]
   Float_t         thetaPFMet[1];   //[nPFMet]
   Float_t         etaPFMet[1];   //[nPFMet]
   Float_t         phiPFMet[1];   //[nPFMet]
   Float_t         pxPFMet[1];   //[nPFMet]
   Float_t         pyPFMet[1];   //[nPFMet]
   Float_t         pzPFMet[1];   //[nPFMet]
   Float_t         vertexXPFMet[1];   //[nPFMet]
   Float_t         vertexYPFMet[1];   //[nPFMet]
   Float_t         vertexZPFMet[1];   //[nPFMet]
   Float_t         massPFMet[1];   //[nPFMet]
   Int_t           nGenMet;
   Int_t           chargeGenMet[1];   //[nGenMet]
   Float_t         energyGenMet[1];   //[nGenMet]
   Float_t         etGenMet[1];   //[nGenMet]
   Float_t         momentumGenMet[1];   //[nGenMet]
   Float_t         thetaGenMet[1];   //[nGenMet]
   Float_t         etaGenMet[1];   //[nGenMet]
   Float_t         phiGenMet[1];   //[nGenMet]
   Float_t         pxGenMet[1];   //[nGenMet]
   Float_t         pyGenMet[1];   //[nGenMet]
   Float_t         pzGenMet[1];   //[nGenMet]
   Float_t         vertexXGenMet[1];   //[nGenMet]
   Float_t         vertexYGenMet[1];   //[nGenMet]
   Float_t         vertexZGenMet[1];   //[nGenMet]
   Float_t         massGenMet[1];   //[nGenMet]
   Int_t           nAK5CorrJet;
   Int_t           chargeAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         energyAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         etAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         momentumAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         thetaAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         etaAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         phiAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         pxAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         pyAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         pzAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         vertexXAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         vertexYAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         vertexZAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         massAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         alphaAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         emFracAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         hadFracAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         combinedSecondaryVertexBJetTagsAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         combinedSecondaryVertexMVABJetTagsAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         jetBProbabilityBJetTagsAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         jetProbabilityBJetTagsAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         simpleSecondaryVertexBJetTagsAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         softMuonBJetTagsAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         trackCountingHighPurBJetTagsAK5CorrJet[37];   //[nAK5CorrJet]
   Float_t         trackCountingHighEffBJetTagsAK5CorrJet[37];   //[nAK5CorrJet]
   Int_t           nAK5Jet;
   Int_t           chargeAK5Jet[37];   //[nAK5Jet]
   Float_t         energyAK5Jet[37];   //[nAK5Jet]
   Float_t         etAK5Jet[37];   //[nAK5Jet]
   Float_t         momentumAK5Jet[37];   //[nAK5Jet]
   Float_t         thetaAK5Jet[37];   //[nAK5Jet]
   Float_t         etaAK5Jet[37];   //[nAK5Jet]
   Float_t         phiAK5Jet[37];   //[nAK5Jet]
   Float_t         pxAK5Jet[37];   //[nAK5Jet]
   Float_t         pyAK5Jet[37];   //[nAK5Jet]
   Float_t         pzAK5Jet[37];   //[nAK5Jet]
   Float_t         vertexXAK5Jet[37];   //[nAK5Jet]
   Float_t         vertexYAK5Jet[37];   //[nAK5Jet]
   Float_t         vertexZAK5Jet[37];   //[nAK5Jet]
   Float_t         massAK5Jet[37];   //[nAK5Jet]
   Float_t         alphaAK5Jet[37];   //[nAK5Jet]
   Float_t         emFracAK5Jet[37];   //[nAK5Jet]
   Float_t         hadFracAK5Jet[37];   //[nAK5Jet]
   Int_t           nAK5PFCorrJet;
   Int_t           chargeAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         energyAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         etAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         momentumAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         thetaAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         etaAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         phiAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         pxAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         pyAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         pzAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         vertexXAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         vertexYAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         vertexZAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         massAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         chargedHadronEnergyAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         neutralHadronEnergyAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         chargedEmEnergyAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         neutralEmEnergyAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         neutralMultiplicityAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         chargedMultiplicityAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Float_t         muonMultiplicityAK5PFCorrJet[64];   //[nAK5PFCorrJet]
   Int_t           nAK5PFJet;
   Int_t           chargeAK5PFJet[64];   //[nAK5PFJet]
   Float_t         energyAK5PFJet[64];   //[nAK5PFJet]
   Float_t         etAK5PFJet[64];   //[nAK5PFJet]
   Float_t         momentumAK5PFJet[64];   //[nAK5PFJet]
   Float_t         thetaAK5PFJet[64];   //[nAK5PFJet]
   Float_t         etaAK5PFJet[64];   //[nAK5PFJet]
   Float_t         phiAK5PFJet[64];   //[nAK5PFJet]
   Float_t         pxAK5PFJet[64];   //[nAK5PFJet]
   Float_t         pyAK5PFJet[64];   //[nAK5PFJet]
   Float_t         pzAK5PFJet[64];   //[nAK5PFJet]
   Float_t         vertexXAK5PFJet[64];   //[nAK5PFJet]
   Float_t         vertexYAK5PFJet[64];   //[nAK5PFJet]
   Float_t         vertexZAK5PFJet[64];   //[nAK5PFJet]
   Float_t         massAK5PFJet[64];   //[nAK5PFJet]
   Float_t         chargedHadronEnergyAK5PFJet[64];   //[nAK5PFJet]
   Float_t         neutralHadronEnergyAK5PFJet[64];   //[nAK5PFJet]
   Float_t         chargedEmEnergyAK5PFJet[64];   //[nAK5PFJet]
   Float_t         neutralEmEnergyAK5PFJet[64];   //[nAK5PFJet]
   Float_t         neutralMultiplicityAK5PFJet[64];   //[nAK5PFJet]
   Float_t         chargedMultiplicityAK5PFJet[64];   //[nAK5PFJet]
   Float_t         muonMultiplicityAK5PFJet[64];   //[nAK5PFJet]
   Int_t           nAK5GenJet;
   Int_t           chargeAK5GenJet[29];   //[nAK5GenJet]
   Float_t         energyAK5GenJet[29];   //[nAK5GenJet]
   Float_t         etAK5GenJet[29];   //[nAK5GenJet]
   Float_t         momentumAK5GenJet[29];   //[nAK5GenJet]
   Float_t         thetaAK5GenJet[29];   //[nAK5GenJet]
   Float_t         etaAK5GenJet[29];   //[nAK5GenJet]
   Float_t         phiAK5GenJet[29];   //[nAK5GenJet]
   Float_t         pxAK5GenJet[29];   //[nAK5GenJet]
   Float_t         pyAK5GenJet[29];   //[nAK5GenJet]
   Float_t         pzAK5GenJet[29];   //[nAK5GenJet]
   Float_t         vertexXAK5GenJet[29];   //[nAK5GenJet]
   Float_t         vertexYAK5GenJet[29];   //[nAK5GenJet]
   Float_t         vertexZAK5GenJet[29];   //[nAK5GenJet]
   Float_t         massAK5GenJet[29];   //[nAK5GenJet]
   Double_t        genPtHat;
   Double_t        genProcessId;
   Double_t        genWeight;
   Double_t        genAlphaQCD;
   Double_t        genAlphaQED;

   // List of branches
   TBranch        *b_nl1Technical;   //!
   TBranch        *b_l1Technical;   //!
   TBranch        *b_nl1Global;   //!
   TBranch        *b_l1Global;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_orbitNumber;   //!
   TBranch        *b_nMc;   //!
   TBranch        *b_pMc;   //!
   TBranch        *b_massMc;   //!
   TBranch        *b_thetaMc;   //!
   TBranch        *b_etaMc;   //!
   TBranch        *b_phiMc;   //!
   TBranch        *b_energyMc;   //!
   TBranch        *b_idMc;   //!
   TBranch        *b_mothMc;   //!
   TBranch        *b_nDauMc;   //!
   TBranch        *b_statusMc;   //!
   TBranch        *b_xMc;   //!
   TBranch        *b_yMc;   //!
   TBranch        *b_zMc;   //!
   TBranch        *b_nTrg;   //!
   TBranch        *b_firedTrg;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_chargeEle;   //!
   TBranch        *b_energyEle;   //!
   TBranch        *b_etEle;   //!
   TBranch        *b_momentumEle;   //!
   TBranch        *b_thetaEle;   //!
   TBranch        *b_etaEle;   //!
   TBranch        *b_phiEle;   //!
   TBranch        *b_pxEle;   //!
   TBranch        *b_pyEle;   //!
   TBranch        *b_pzEle;   //!
   TBranch        *b_vertexXEle;   //!
   TBranch        *b_vertexYEle;   //!
   TBranch        *b_vertexZEle;   //!
   TBranch        *b_massEle;   //!
   TBranch        *b_fiducialFlagsEle;   //!
   TBranch        *b_recoFlagsEle;   //!
   TBranch        *b_energyCorrectionsEle;   //!
   TBranch        *b_esEnergyEle;   //!
   TBranch        *b_superClusterIndexEle;   //!
   TBranch        *b_PFsuperClusterIndexEle;   //!
   TBranch        *b_trackIndexEle;   //!
   TBranch        *b_gsfTrackIndexEle;   //!
   TBranch        *b_classificationEle;   //!
   TBranch        *b_standardClassificationEle;   //!
   TBranch        *b_hOverEEle;   //!
   TBranch        *b_eSuperClusterOverPEle;   //!
   TBranch        *b_eSeedOverPoutEle;   //!
   TBranch        *b_deltaEtaAtVtxEle;   //!
   TBranch        *b_deltaPhiAtVtxEle;   //!
   TBranch        *b_deltaEtaAtCaloEle;   //!
   TBranch        *b_deltaPhiAtCaloEle;   //!
   TBranch        *b_tipEle;   //!
   TBranch        *b_dr03TkSumPtEle;   //!
   TBranch        *b_dr03EcalRecHitSumEtEle;   //!
   TBranch        *b_dr03HcalTowerSumEtEle;   //!
   TBranch        *b_dr04TkSumPtEle;   //!
   TBranch        *b_dr04EcalRecHitSumEtEle;   //!
   TBranch        *b_dr04HcalTowerSumEtEle;   //!
   TBranch        *b_scBasedEcalSum03Ele;   //!
   TBranch        *b_scBasedEcalSum04Ele;   //!
   TBranch        *b_scHaloBasedEcalSum03Ele;   //!
   TBranch        *b_scHaloBasedEcalSum04Ele;   //!
   TBranch        *b_eleIdCutsEle;   //!
   TBranch        *b_eleIdLikelihoodEle;   //!
   TBranch        *b_pflowMVAEle;   //!
   TBranch        *b_nPFEle;   //!
   TBranch        *b_chargePFEle;   //!
   TBranch        *b_energyPFEle;   //!
   TBranch        *b_etPFEle;   //!
   TBranch        *b_momentumPFEle;   //!
   TBranch        *b_thetaPFEle;   //!
   TBranch        *b_etaPFEle;   //!
   TBranch        *b_phiPFEle;   //!
   TBranch        *b_pxPFEle;   //!
   TBranch        *b_pyPFEle;   //!
   TBranch        *b_pzPFEle;   //!
   TBranch        *b_vertexXPFEle;   //!
   TBranch        *b_vertexYPFEle;   //!
   TBranch        *b_vertexZPFEle;   //!
   TBranch        *b_massPFEle;   //!
   TBranch        *b_MvaOutputPFEle;   //!
   TBranch        *b_PS1EnergyPFEle;   //!
   TBranch        *b_PS2EnergyPFEle;   //!
   TBranch        *b_EcalEnergyPFEle;   //!
   TBranch        *b_EcalElectronEnergyPFEle;   //!
   TBranch        *b_gsfTrackIndexPFEle;   //!
   TBranch        *b_trackIndexPFEle;   //!
   TBranch        *b_nSC;   //!
   TBranch        *b_nBCSC;   //!
   TBranch        *b_nCrystalsSC;   //!
   TBranch        *b_iAlgoSC;   //!
   TBranch        *b_rawEnergySC;   //!
   TBranch        *b_energySC;   //!
   TBranch        *b_etaSC;   //!
   TBranch        *b_thetaSC;   //!
   TBranch        *b_phiSC;   //!
   TBranch        *b_e3x3SC;   //!
   TBranch        *b_e5x5SC;   //!
   TBranch        *b_eMaxSC;   //!
   TBranch        *b_e2x2SC;   //!
   TBranch        *b_e2ndSC;   //!
   TBranch        *b_covIEtaIEtaSC;   //!
   TBranch        *b_covIEtaIPhiSC;   //!
   TBranch        *b_covIPhiIPhiSC;   //!
   TBranch        *b_hOverESC;   //!
   TBranch        *b_recoFlagSC;   //!
   TBranch        *b_channelStatusSC;   //!
   TBranch        *b_timeSC;   //!
   TBranch        *b_chi2ProbSC;   //!
   TBranch        *b_seedEnergySC;   //!
   TBranch        *b_idClosProblSC;   //!
   TBranch        *b_sevClosProblSC;   //!
   TBranch        *b_fracClosProblSC;   //!
   TBranch        *b_trackIndexSC;   //!
   TBranch        *b_trackDeltaRSC;   //!
   TBranch        *b_trackDeltaPhiSC;   //!
   TBranch        *b_trackDeltaEtaSC;   //!
   TBranch        *b_gsfTrackIndexSC;   //!
   TBranch        *b_gsfTrackDeltaRSC;   //!
   TBranch        *b_gsfTrackDeltaPhiSC;   //!
   TBranch        *b_gsfTrackDeltaEtaSC;   //!
   TBranch        *b_pxVtxPropagatedNegChargeSC;   //!
   TBranch        *b_pyVtxPropagatedNegChargeSC;   //!
   TBranch        *b_pzVtxPropagatedNegChargeSC;   //!
   TBranch        *b_pxVtxPropagatedPosChargeSC;   //!
   TBranch        *b_pyVtxPropagatedPosChargeSC;   //!
   TBranch        *b_pzVtxPropagatedPosChargeSC;   //!
   TBranch        *b_nPFSC;   //!
   TBranch        *b_nBCPFSC;   //!
   TBranch        *b_nCrystalsPFSC;   //!
   TBranch        *b_iAlgoPFSC;   //!
   TBranch        *b_rawEnergyPFSC;   //!
   TBranch        *b_energyPFSC;   //!
   TBranch        *b_etaPFSC;   //!
   TBranch        *b_thetaPFSC;   //!
   TBranch        *b_phiPFSC;   //!
   TBranch        *b_e3x3PFSC;   //!
   TBranch        *b_e5x5PFSC;   //!
   TBranch        *b_eMaxPFSC;   //!
   TBranch        *b_e2x2PFSC;   //!
   TBranch        *b_e2ndPFSC;   //!
   TBranch        *b_covIEtaIEtaPFSC;   //!
   TBranch        *b_covIEtaIPhiPFSC;   //!
   TBranch        *b_covIPhiIPhiPFSC;   //!
   TBranch        *b_hOverEPFSC;   //!
   TBranch        *b_recoFlagPFSC;   //!
   TBranch        *b_channelStatusPFSC;   //!
   TBranch        *b_timePFSC;   //!
   TBranch        *b_chi2ProbPFSC;   //!
   TBranch        *b_seedEnergyPFSC;   //!
   TBranch        *b_idClosProblPFSC;   //!
   TBranch        *b_sevClosProblPFSC;   //!
   TBranch        *b_fracClosProblPFSC;   //!
   TBranch        *b_pxVtxPropagatedNegChargePFSC;   //!
   TBranch        *b_pyVtxPropagatedNegChargePFSC;   //!
   TBranch        *b_pzVtxPropagatedNegChargePFSC;   //!
   TBranch        *b_pxVtxPropagatedPosChargePFSC;   //!
   TBranch        *b_pyVtxPropagatedPosChargePFSC;   //!
   TBranch        *b_pzVtxPropagatedPosChargePFSC;   //!
   TBranch        *b_nTrack;   //!
   TBranch        *b_pxTrack;   //!
   TBranch        *b_pyTrack;   //!
   TBranch        *b_pzTrack;   //!
   TBranch        *b_vtxIndexTrack;   //!
   TBranch        *b_vtxWeightTrack;   //!
   TBranch        *b_chargeTrack;   //!
   TBranch        *b_ptErrorTrack;   //!
   TBranch        *b_trackValidHitsTrack;   //!
   TBranch        *b_trackLostHitsTrack;   //!
   TBranch        *b_trackNormalizedChi2Track;   //!
   TBranch        *b_qualityMaskTrack;   //!
   TBranch        *b_trackDxyTrack;   //!
   TBranch        *b_trackD0Track;   //!
   TBranch        *b_trackDszTrack;   //!
   TBranch        *b_trackDzTrack;   //!
   TBranch        *b_trackDxyErrorTrack;   //!
   TBranch        *b_trackD0ErrorTrack;   //!
   TBranch        *b_trackDszErrorTrack;   //!
   TBranch        *b_trackDzErrorTrack;   //!
   TBranch        *b_trackDxyPVTrack;   //!
   TBranch        *b_trackDszPVTrack;   //!
   TBranch        *b_trackDzPVTrack;   //!
   TBranch        *b_trackVxTrack;   //!
   TBranch        *b_trackVyTrack;   //!
   TBranch        *b_trackVzTrack;   //!
   TBranch        *b_pxAtOuterTrack;   //!
   TBranch        *b_pyAtOuterTrack;   //!
   TBranch        *b_pzAtOuterTrack;   //!
   TBranch        *b_xAtOuterTrack;   //!
   TBranch        *b_yAtOuterTrack;   //!
   TBranch        *b_zAtOuterTrack;   //!
   TBranch        *b_pxAtInnerTrack;   //!
   TBranch        *b_pyAtInnerTrack;   //!
   TBranch        *b_pzAtInnerTrack;   //!
   TBranch        *b_xAtInnerTrack;   //!
   TBranch        *b_yAtInnerTrack;   //!
   TBranch        *b_zAtInnerTrack;   //!
   TBranch        *b_recHitsSizeTrack;   //!
   TBranch        *b_isPixB1Track;   //!
   TBranch        *b_isPixB2Track;   //!
   TBranch        *b_isPixE1Track;   //!
   TBranch        *b_isPixE2Track;   //!
   TBranch        *b_numberOfValidPixelBarrelHitsTrack;   //!
   TBranch        *b_numberOfValidPixelEndcapHitsTrack;   //!
   TBranch        *b_numberOfValidStripTIBHitsTrack;   //!
   TBranch        *b_numberOfValidStripTIDHitsTrack;   //!
   TBranch        *b_numberOfValidStripTOBHitsTrack;   //!
   TBranch        *b_numberOfValidStripTECHitsTrack;   //!
   TBranch        *b_truncatedDeDxTrack;   //!
   TBranch        *b_truncatedDeDxErrorTrack;   //!
   TBranch        *b_truncatedDeDxNoMTrack;   //!
   TBranch        *b_medianDeDxTrack;   //!
   TBranch        *b_medianDeDxErrorTrack;   //!
   TBranch        *b_medianDeDxNoMTrack;   //!
   TBranch        *b_harmonic2DeDxTrack;   //!
   TBranch        *b_harmonic2DeDxErrorTrack;   //!
   TBranch        *b_harmonic2DeDxNoMTrack;   //!
   TBranch        *b_nGsfTrack;   //!
   TBranch        *b_pxGsfTrack;   //!
   TBranch        *b_pyGsfTrack;   //!
   TBranch        *b_pzGsfTrack;   //!
   TBranch        *b_chargeGsfTrack;   //!
   TBranch        *b_ptErrorGsfTrack;   //!
   TBranch        *b_trackValidHitsGsfTrack;   //!
   TBranch        *b_trackLostHitsGsfTrack;   //!
   TBranch        *b_trackNormalizedChi2GsfTrack;   //!
   TBranch        *b_qualityMaskGsfTrack;   //!
   TBranch        *b_trackDxyGsfTrack;   //!
   TBranch        *b_trackD0GsfTrack;   //!
   TBranch        *b_trackDszGsfTrack;   //!
   TBranch        *b_trackDzGsfTrack;   //!
   TBranch        *b_trackDxyErrorGsfTrack;   //!
   TBranch        *b_trackD0ErrorGsfTrack;   //!
   TBranch        *b_trackDszErrorGsfTrack;   //!
   TBranch        *b_trackDzErrorGsfTrack;   //!
   TBranch        *b_trackDxyPVGsfTrack;   //!
   TBranch        *b_trackDszPVGsfTrack;   //!
   TBranch        *b_trackDzPVGsfTrack;   //!
   TBranch        *b_trackVxGsfTrack;   //!
   TBranch        *b_trackVyGsfTrack;   //!
   TBranch        *b_trackVzGsfTrack;   //!
   TBranch        *b_pxAtOuterGsfTrack;   //!
   TBranch        *b_pyAtOuterGsfTrack;   //!
   TBranch        *b_pzAtOuterGsfTrack;   //!
   TBranch        *b_xAtOuterGsfTrack;   //!
   TBranch        *b_yAtOuterGsfTrack;   //!
   TBranch        *b_zAtOuterGsfTrack;   //!
   TBranch        *b_pxAtInnerGsfTrack;   //!
   TBranch        *b_pyAtInnerGsfTrack;   //!
   TBranch        *b_pzAtInnerGsfTrack;   //!
   TBranch        *b_xAtInnerGsfTrack;   //!
   TBranch        *b_yAtInnerGsfTrack;   //!
   TBranch        *b_zAtInnerGsfTrack;   //!
   TBranch        *b_recHitsSizeGsfTrack;   //!
   TBranch        *b_isPixB1GsfTrack;   //!
   TBranch        *b_isPixB2GsfTrack;   //!
   TBranch        *b_isPixE1GsfTrack;   //!
   TBranch        *b_isPixE2GsfTrack;   //!
   TBranch        *b_numberOfValidPixelBarrelHitsGsfTrack;   //!
   TBranch        *b_numberOfValidPixelEndcapHitsGsfTrack;   //!
   TBranch        *b_numberOfValidStripTIBHitsGsfTrack;   //!
   TBranch        *b_numberOfValidStripTIDHitsGsfTrack;   //!
   TBranch        *b_numberOfValidStripTOBHitsGsfTrack;   //!
   TBranch        *b_numberOfValidStripTECHitsGsfTrack;   //!
   TBranch        *b_chargeModeGsfTrack;   //!
   TBranch        *b_pxModeGsfTrack;   //!
   TBranch        *b_pyModeGsfTrack;   //!
   TBranch        *b_pzModeGsfTrack;   //!
   TBranch        *b_recoFlagsGsfTrack;   //!
   TBranch        *b_nGlobalMuonTrack;   //!
   TBranch        *b_pxGlobalMuonTrack;   //!
   TBranch        *b_pyGlobalMuonTrack;   //!
   TBranch        *b_pzGlobalMuonTrack;   //!
   TBranch        *b_chargeGlobalMuonTrack;   //!
   TBranch        *b_ptErrorGlobalMuonTrack;   //!
   TBranch        *b_trackValidHitsGlobalMuonTrack;   //!
   TBranch        *b_trackLostHitsGlobalMuonTrack;   //!
   TBranch        *b_trackNormalizedChi2GlobalMuonTrack;   //!
   TBranch        *b_qualityMaskGlobalMuonTrack;   //!
   TBranch        *b_trackDxyGlobalMuonTrack;   //!
   TBranch        *b_trackD0GlobalMuonTrack;   //!
   TBranch        *b_trackDszGlobalMuonTrack;   //!
   TBranch        *b_trackDzGlobalMuonTrack;   //!
   TBranch        *b_trackDxyErrorGlobalMuonTrack;   //!
   TBranch        *b_trackD0ErrorGlobalMuonTrack;   //!
   TBranch        *b_trackDszErrorGlobalMuonTrack;   //!
   TBranch        *b_trackDzErrorGlobalMuonTrack;   //!
   TBranch        *b_trackDxyPVGlobalMuonTrack;   //!
   TBranch        *b_trackDszPVGlobalMuonTrack;   //!
   TBranch        *b_trackDzPVGlobalMuonTrack;   //!
   TBranch        *b_trackVxGlobalMuonTrack;   //!
   TBranch        *b_trackVyGlobalMuonTrack;   //!
   TBranch        *b_trackVzGlobalMuonTrack;   //!
   TBranch        *b_pxAtOuterGlobalMuonTrack;   //!
   TBranch        *b_pyAtOuterGlobalMuonTrack;   //!
   TBranch        *b_pzAtOuterGlobalMuonTrack;   //!
   TBranch        *b_xAtOuterGlobalMuonTrack;   //!
   TBranch        *b_yAtOuterGlobalMuonTrack;   //!
   TBranch        *b_zAtOuterGlobalMuonTrack;   //!
   TBranch        *b_pxAtInnerGlobalMuonTrack;   //!
   TBranch        *b_pyAtInnerGlobalMuonTrack;   //!
   TBranch        *b_pzAtInnerGlobalMuonTrack;   //!
   TBranch        *b_xAtInnerGlobalMuonTrack;   //!
   TBranch        *b_yAtInnerGlobalMuonTrack;   //!
   TBranch        *b_zAtInnerGlobalMuonTrack;   //!
   TBranch        *b_recHitsSizeGlobalMuonTrack;   //!
   TBranch        *b_isPixB1GlobalMuonTrack;   //!
   TBranch        *b_isPixB2GlobalMuonTrack;   //!
   TBranch        *b_isPixE1GlobalMuonTrack;   //!
   TBranch        *b_isPixE2GlobalMuonTrack;   //!
   TBranch        *b_numberOfValidPixelBarrelHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidPixelEndcapHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidStripTIBHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidStripTIDHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidStripTOBHitsGlobalMuonTrack;   //!
   TBranch        *b_numberOfValidStripTECHitsGlobalMuonTrack;   //!
   TBranch        *b_nSTAMuonTrack;   //!
   TBranch        *b_pxSTAMuonTrack;   //!
   TBranch        *b_pySTAMuonTrack;   //!
   TBranch        *b_pzSTAMuonTrack;   //!
   TBranch        *b_chargeSTAMuonTrack;   //!
   TBranch        *b_ptErrorSTAMuonTrack;   //!
   TBranch        *b_trackValidHitsSTAMuonTrack;   //!
   TBranch        *b_trackLostHitsSTAMuonTrack;   //!
   TBranch        *b_trackNormalizedChi2STAMuonTrack;   //!
   TBranch        *b_qualityMaskSTAMuonTrack;   //!
   TBranch        *b_trackDxySTAMuonTrack;   //!
   TBranch        *b_trackD0STAMuonTrack;   //!
   TBranch        *b_trackDszSTAMuonTrack;   //!
   TBranch        *b_trackDzSTAMuonTrack;   //!
   TBranch        *b_trackDxyErrorSTAMuonTrack;   //!
   TBranch        *b_trackD0ErrorSTAMuonTrack;   //!
   TBranch        *b_trackDszErrorSTAMuonTrack;   //!
   TBranch        *b_trackDzErrorSTAMuonTrack;   //!
   TBranch        *b_trackDxyPVSTAMuonTrack;   //!
   TBranch        *b_trackDszPVSTAMuonTrack;   //!
   TBranch        *b_trackDzPVSTAMuonTrack;   //!
   TBranch        *b_trackVxSTAMuonTrack;   //!
   TBranch        *b_trackVySTAMuonTrack;   //!
   TBranch        *b_trackVzSTAMuonTrack;   //!
   TBranch        *b_pxAtOuterSTAMuonTrack;   //!
   TBranch        *b_pyAtOuterSTAMuonTrack;   //!
   TBranch        *b_pzAtOuterSTAMuonTrack;   //!
   TBranch        *b_xAtOuterSTAMuonTrack;   //!
   TBranch        *b_yAtOuterSTAMuonTrack;   //!
   TBranch        *b_zAtOuterSTAMuonTrack;   //!
   TBranch        *b_pxAtInnerSTAMuonTrack;   //!
   TBranch        *b_pyAtInnerSTAMuonTrack;   //!
   TBranch        *b_pzAtInnerSTAMuonTrack;   //!
   TBranch        *b_xAtInnerSTAMuonTrack;   //!
   TBranch        *b_yAtInnerSTAMuonTrack;   //!
   TBranch        *b_zAtInnerSTAMuonTrack;   //!
   TBranch        *b_recHitsSizeSTAMuonTrack;   //!
   TBranch        *b_isPixB1STAMuonTrack;   //!
   TBranch        *b_isPixB2STAMuonTrack;   //!
   TBranch        *b_isPixE1STAMuonTrack;   //!
   TBranch        *b_isPixE2STAMuonTrack;   //!
   TBranch        *b_numberOfValidPixelBarrelHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidPixelEndcapHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidStripTIBHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidStripTIDHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidStripTOBHitsSTAMuonTrack;   //!
   TBranch        *b_numberOfValidStripTECHitsSTAMuonTrack;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_PVxPV;   //!
   TBranch        *b_PVyPV;   //!
   TBranch        *b_PVzPV;   //!
   TBranch        *b_PVErrxPV;   //!
   TBranch        *b_PVErryPV;   //!
   TBranch        *b_PVErrzPV;   //!
   TBranch        *b_SumPtPV;   //!
   TBranch        *b_ndofPV;   //!
   TBranch        *b_chi2PV;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_chargeMuon;   //!
   TBranch        *b_energyMuon;   //!
   TBranch        *b_etMuon;   //!
   TBranch        *b_momentumMuon;   //!
   TBranch        *b_thetaMuon;   //!
   TBranch        *b_etaMuon;   //!
   TBranch        *b_phiMuon;   //!
   TBranch        *b_pxMuon;   //!
   TBranch        *b_pyMuon;   //!
   TBranch        *b_pzMuon;   //!
   TBranch        *b_vertexXMuon;   //!
   TBranch        *b_vertexYMuon;   //!
   TBranch        *b_vertexZMuon;   //!
   TBranch        *b_massMuon;   //!
   TBranch        *b_trackIndexMuon;   //!
   TBranch        *b_standAloneTrackIndexMuon;   //!
   TBranch        *b_combinedTrackIndexMuon;   //!
   TBranch        *b_muonIdMuon;   //!
   TBranch        *b_sumPt03Muon;   //!
   TBranch        *b_emEt03Muon;   //!
   TBranch        *b_hadEt03Muon;   //!
   TBranch        *b_hoEt03Muon;   //!
   TBranch        *b_nTrk03Muon;   //!
   TBranch        *b_nJets03Muon;   //!
   TBranch        *b_sumPt05Muon;   //!
   TBranch        *b_emEt05Muon;   //!
   TBranch        *b_hadEt05Muon;   //!
   TBranch        *b_hoEt05Muon;   //!
   TBranch        *b_nTrk05Muon;   //!
   TBranch        *b_nJets05Muon;   //!
   TBranch        *b_EcalExpDepoMuon;   //!
   TBranch        *b_HcalExpDepoMuon;   //!
   TBranch        *b_HoExpDepoMuon;   //!
   TBranch        *b_emS9Muon;   //!
   TBranch        *b_hadS9Muon;   //!
   TBranch        *b_hoS9Muon;   //!
   TBranch        *b_CaloCompMuon;   //!
   TBranch        *b_nCaloTowers;   //!
   TBranch        *b_energyCaloTowers;   //!
   TBranch        *b_xCaloTowers;   //!
   TBranch        *b_yCaloTowers;   //!
   TBranch        *b_zCaloTowers;   //!
   TBranch        *b_CALOCaloTowers;   //!
   TBranch        *b_CaloIndexCaloTowers;   //!
   TBranch        *b_nMet;   //!
   TBranch        *b_chargeMet;   //!
   TBranch        *b_energyMet;   //!
   TBranch        *b_etMet;   //!
   TBranch        *b_momentumMet;   //!
   TBranch        *b_thetaMet;   //!
   TBranch        *b_etaMet;   //!
   TBranch        *b_phiMet;   //!
   TBranch        *b_pxMet;   //!
   TBranch        *b_pyMet;   //!
   TBranch        *b_pzMet;   //!
   TBranch        *b_vertexXMet;   //!
   TBranch        *b_vertexYMet;   //!
   TBranch        *b_vertexZMet;   //!
   TBranch        *b_massMet;   //!
   TBranch        *b_nTCMet;   //!
   TBranch        *b_chargeTCMet;   //!
   TBranch        *b_energyTCMet;   //!
   TBranch        *b_etTCMet;   //!
   TBranch        *b_momentumTCMet;   //!
   TBranch        *b_thetaTCMet;   //!
   TBranch        *b_etaTCMet;   //!
   TBranch        *b_phiTCMet;   //!
   TBranch        *b_pxTCMet;   //!
   TBranch        *b_pyTCMet;   //!
   TBranch        *b_pzTCMet;   //!
   TBranch        *b_vertexXTCMet;   //!
   TBranch        *b_vertexYTCMet;   //!
   TBranch        *b_vertexZTCMet;   //!
   TBranch        *b_massTCMet;   //!
   TBranch        *b_nPFMet;   //!
   TBranch        *b_chargePFMet;   //!
   TBranch        *b_energyPFMet;   //!
   TBranch        *b_etPFMet;   //!
   TBranch        *b_momentumPFMet;   //!
   TBranch        *b_thetaPFMet;   //!
   TBranch        *b_etaPFMet;   //!
   TBranch        *b_phiPFMet;   //!
   TBranch        *b_pxPFMet;   //!
   TBranch        *b_pyPFMet;   //!
   TBranch        *b_pzPFMet;   //!
   TBranch        *b_vertexXPFMet;   //!
   TBranch        *b_vertexYPFMet;   //!
   TBranch        *b_vertexZPFMet;   //!
   TBranch        *b_massPFMet;   //!
   TBranch        *b_nGenMet;   //!
   TBranch        *b_chargeGenMet;   //!
   TBranch        *b_energyGenMet;   //!
   TBranch        *b_etGenMet;   //!
   TBranch        *b_momentumGenMet;   //!
   TBranch        *b_thetaGenMet;   //!
   TBranch        *b_etaGenMet;   //!
   TBranch        *b_phiGenMet;   //!
   TBranch        *b_pxGenMet;   //!
   TBranch        *b_pyGenMet;   //!
   TBranch        *b_pzGenMet;   //!
   TBranch        *b_vertexXGenMet;   //!
   TBranch        *b_vertexYGenMet;   //!
   TBranch        *b_vertexZGenMet;   //!
   TBranch        *b_massGenMet;   //!
   TBranch        *b_nAK5CorrJet;   //!
   TBranch        *b_chargeAK5CorrJet;   //!
   TBranch        *b_energyAK5CorrJet;   //!
   TBranch        *b_etAK5CorrJet;   //!
   TBranch        *b_momentumAK5CorrJet;   //!
   TBranch        *b_thetaAK5CorrJet;   //!
   TBranch        *b_etaAK5CorrJet;   //!
   TBranch        *b_phiAK5CorrJet;   //!
   TBranch        *b_pxAK5CorrJet;   //!
   TBranch        *b_pyAK5CorrJet;   //!
   TBranch        *b_pzAK5CorrJet;   //!
   TBranch        *b_vertexXAK5CorrJet;   //!
   TBranch        *b_vertexYAK5CorrJet;   //!
   TBranch        *b_vertexZAK5CorrJet;   //!
   TBranch        *b_massAK5CorrJet;   //!
   TBranch        *b_alphaAK5CorrJet;   //!
   TBranch        *b_emFracAK5CorrJet;   //!
   TBranch        *b_hadFracAK5CorrJet;   //!
   TBranch        *b_combinedSecondaryVertexBJetTagsAK5CorrJet;   //!
   TBranch        *b_combinedSecondaryVertexMVABJetTagsAK5CorrJet;   //!
   TBranch        *b_jetBProbabilityBJetTagsAK5CorrJet;   //!
   TBranch        *b_jetProbabilityBJetTagsAK5CorrJet;   //!
   TBranch        *b_simpleSecondaryVertexBJetTagsAK5CorrJet;   //!
   TBranch        *b_softMuonBJetTagsAK5CorrJet;   //!
   TBranch        *b_trackCountingHighPurBJetTagsAK5CorrJet;   //!
   TBranch        *b_trackCountingHighEffBJetTagsAK5CorrJet;   //!
   TBranch        *b_nAK5Jet;   //!
   TBranch        *b_chargeAK5Jet;   //!
   TBranch        *b_energyAK5Jet;   //!
   TBranch        *b_etAK5Jet;   //!
   TBranch        *b_momentumAK5Jet;   //!
   TBranch        *b_thetaAK5Jet;   //!
   TBranch        *b_etaAK5Jet;   //!
   TBranch        *b_phiAK5Jet;   //!
   TBranch        *b_pxAK5Jet;   //!
   TBranch        *b_pyAK5Jet;   //!
   TBranch        *b_pzAK5Jet;   //!
   TBranch        *b_vertexXAK5Jet;   //!
   TBranch        *b_vertexYAK5Jet;   //!
   TBranch        *b_vertexZAK5Jet;   //!
   TBranch        *b_massAK5Jet;   //!
   TBranch        *b_alphaAK5Jet;   //!
   TBranch        *b_emFracAK5Jet;   //!
   TBranch        *b_hadFracAK5Jet;   //!
   TBranch        *b_nAK5PFCorrJet;   //!
   TBranch        *b_chargeAK5PFCorrJet;   //!
   TBranch        *b_energyAK5PFCorrJet;   //!
   TBranch        *b_etAK5PFCorrJet;   //!
   TBranch        *b_momentumAK5PFCorrJet;   //!
   TBranch        *b_thetaAK5PFCorrJet;   //!
   TBranch        *b_etaAK5PFCorrJet;   //!
   TBranch        *b_phiAK5PFCorrJet;   //!
   TBranch        *b_pxAK5PFCorrJet;   //!
   TBranch        *b_pyAK5PFCorrJet;   //!
   TBranch        *b_pzAK5PFCorrJet;   //!
   TBranch        *b_vertexXAK5PFCorrJet;   //!
   TBranch        *b_vertexYAK5PFCorrJet;   //!
   TBranch        *b_vertexZAK5PFCorrJet;   //!
   TBranch        *b_massAK5PFCorrJet;   //!
   TBranch        *b_chargedHadronEnergyAK5PFCorrJet;   //!
   TBranch        *b_neutralHadronEnergyAK5PFCorrJet;   //!
   TBranch        *b_chargedEmEnergyAK5PFCorrJet;   //!
   TBranch        *b_neutralEmEnergyAK5PFCorrJet;   //!
   TBranch        *b_neutralMultiplicityAK5PFCorrJet;   //!
   TBranch        *b_chargedMultiplicityAK5PFCorrJet;   //!
   TBranch        *b_muonMultiplicityAK5PFCorrJet;   //!
   TBranch        *b_nAK5PFJet;   //!
   TBranch        *b_chargeAK5PFJet;   //!
   TBranch        *b_energyAK5PFJet;   //!
   TBranch        *b_etAK5PFJet;   //!
   TBranch        *b_momentumAK5PFJet;   //!
   TBranch        *b_thetaAK5PFJet;   //!
   TBranch        *b_etaAK5PFJet;   //!
   TBranch        *b_phiAK5PFJet;   //!
   TBranch        *b_pxAK5PFJet;   //!
   TBranch        *b_pyAK5PFJet;   //!
   TBranch        *b_pzAK5PFJet;   //!
   TBranch        *b_vertexXAK5PFJet;   //!
   TBranch        *b_vertexYAK5PFJet;   //!
   TBranch        *b_vertexZAK5PFJet;   //!
   TBranch        *b_massAK5PFJet;   //!
   TBranch        *b_chargedHadronEnergyAK5PFJet;   //!
   TBranch        *b_neutralHadronEnergyAK5PFJet;   //!
   TBranch        *b_chargedEmEnergyAK5PFJet;   //!
   TBranch        *b_neutralEmEnergyAK5PFJet;   //!
   TBranch        *b_neutralMultiplicityAK5PFJet;   //!
   TBranch        *b_chargedMultiplicityAK5PFJet;   //!
   TBranch        *b_muonMultiplicityAK5PFJet;   //!
   TBranch        *b_nAK5GenJet;   //!
   TBranch        *b_chargeAK5GenJet;   //!
   TBranch        *b_energyAK5GenJet;   //!
   TBranch        *b_etAK5GenJet;   //!
   TBranch        *b_momentumAK5GenJet;   //!
   TBranch        *b_thetaAK5GenJet;   //!
   TBranch        *b_etaAK5GenJet;   //!
   TBranch        *b_phiAK5GenJet;   //!
   TBranch        *b_pxAK5GenJet;   //!
   TBranch        *b_pyAK5GenJet;   //!
   TBranch        *b_pzAK5GenJet;   //!
   TBranch        *b_vertexXAK5GenJet;   //!
   TBranch        *b_vertexYAK5GenJet;   //!
   TBranch        *b_vertexZAK5GenJet;   //!
   TBranch        *b_massAK5GenJet;   //!
   TBranch        *b_genPtHat;   //!
   TBranch        *b_genProcessId;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_genAlphaQCD;   //!
   TBranch        *b_genAlphaQED;   //!

   SUSYTauTauAnalysis(TTree *tree=0);
   virtual ~SUSYTauTauAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //   virtual void     Loop();
   virtual void     Loop( Int_t isDetailedHad = 0 );
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   Int_t TauHadTyp;

};

#endif

#ifdef SUSYTauTauAnalysis_cxx
SUSYTauTauAnalysis::SUSYTauTauAnalysis(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
    //      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("rfio:///castor/cern.ch/user/e/emanuele/CMST3/Vecbos2010/MC/7TeV/SUSY/LM1/default_MC_16.root");     
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TauTau/default_MC_16.root"); 
    if (!f) {
      /* f = new TFile("rfio:///castor/cern.ch/user/e/emanuele/CMST3/Vecbos2010/MC/7TeV/SUSY/LM1/default_MC_16.root"); */
      //f = new TFile("~/scratch0/default_MC_16.root");	
      
      //f = new TFile("rfio:/castor/cern.ch/user/e/emanuele/CMST3/Vecbos2010/MC/7TeV/SUSY/LM1/default_MC_16.root");
      //f = new TFile("/tmp/azzolini/TauTau/default_MC_16.root");   // lxplus255
      f = new TFile("~/work/ANALYSIS/TauTau/default_MC_16.root");   // CAT 
    }
    tree = (TTree*)gDirectory->Get("ntp1");
    
    //// vir --- questo funziona: crea root file e tree, ma le variabili sono vuote   
    //// In the following, 'events' is a tree that is copied from 'b0.root' to the new file 'new.root'
    //TFile g("/tmp/azzolini/new.root", "RECREATE");
    //tree->Write();
    //g.Close();
    //// vir --- questo funziona: crea root file e tree, ma le variabili sono vuote   
  
    //babar dice: 
    // In the following, 'events' is a tree that is copied from 'b0.root' to the new file 'new.root'
    // TFile f("b0.root")
    //TFile g("new.root", "RECREATE")
    //((TTree*)f.Get("events"))->Write()
    //g.Close();
    
    // prova 2 - $ROOTSYS/tutorials/tree/copytree2.C
    //Create a new file + a clone of old tree header. Do not copy events
    //TFile *newfile = new TFile("small.root","recreate");
    //TTree *newtree = tree->CloneTree(0);
    
    ////Divert branch fH to a separate file and copy all events
    //newtree->GetBranch("fH")->SetFile("small_fH.root");
    //newtree->CopyEntries(tree);
    //newtree->Print();
    //newfile->Write();
    //delete f;
    //  delete oldfile;
    //delete newfile;
    
    
  }
  
  //  ---- VIR -------
  // -- The file with the big tree
  //TFile *oldfile = new TFile("~/work/ANALYSIS/TauTau/default_MC_16.root");   // CAT 
  //TTree *oldtree = (TTree*)oldfile->Get("ntp1");
  //  Event *event = new Event();
  
  
  Init(tree);
  
  //Create a new file + a clone of old tree header. Do not copy events
  TFile *newfile = new TFile("small.root","recreate");
  TTree *newtree = tree->CloneTree(0);
  
   //Divert branch fH to a separate file and copy all events
  //   newtree->GetBranch("nMc")->SetFile("small_nMc.root");
   newtree->CopyEntries(tree);

   newtree->Print();
   newfile->Write();
   //delete f;
   //delete newfile;
}

SUSYTauTauAnalysis::~SUSYTauTauAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SUSYTauTauAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SUSYTauTauAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void SUSYTauTauAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

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
   Notify();
}

Bool_t SUSYTauTauAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SUSYTauTauAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SUSYTauTauAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SUSYTauTauAnalysis_cxx
