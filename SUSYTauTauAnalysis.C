#define SUSYTauTauAnalysis_cxx
#include "SUSYTauTauAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

void SUSYTauTauAnalysis::Loop(Int_t isDetailedHad)
{
  //   In a ROOT session, you can do:
  //      Root > .L SUSYTauTauAnalysis.C
  //      Root > > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  vector<Int_t> TauIndexVector;
  vector<Int_t> SusyIndexVector;
  vector<Int_t> MuonIndexVector;
  vector<Int_t> ElectronIndexVector;
  vector<Int_t> SusyMuonVector;
  vector<Int_t> SusyElectronVector;
  short TauCounter[4] = {0, 0, 0, 0};
  short TauHadCounter[7] = {0, 0, 0, 0, 0, 0, 0};
  Int_t nCascade =0;
  Int_t nTauEvents = 0;
  Int_t nSUSYTaus = 0;



  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    // for (Long64_t jentry=0; jentry<500.; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (Cut(ientry) < 0) continue;

    // Loop on the true particles and fill the initial vectors of Taus, Muons and Electrons
    for (Int_t i=0; i<nMc; i++) {
      Int_t m = mothMc[i];

      if ( fabs(idMc[i])==15 && ( fabs(idMc[m]) == 15 || (idMc[m]>1000001 && idMc[m]<1000039) || (idMc[m]>2000001 && idMc[m]<2000015)) )
	TauIndexVector.push_back(i);
      else if ( fabs(idMc[i])==13 && (fabs(idMc[m])==15 || fabs(idMc[m])==13) )
	MuonIndexVector.push_back(i);
      else if ( fabs(idMc[i])==11  && (fabs(idMc[m])==15 || fabs(idMc[m])==11) )
	ElectronIndexVector.push_back(i);
      else continue;
    } // end of loop on true particles

    if (TauIndexVector.size() > 1) {
      nTauEvents++;
      // Loop over list of Taus
      for (Int_t ii=0; ii<TauIndexVector.size(); ii++) {
	Int_t TauIndex = TauIndexVector[ii];
	Int_t MothIndex = mothMc[TauIndex];
	Int_t GrandMothIndex = mothMc[MothIndex];
	Int_t GreatGrandMothIndex = mothMc[GrandMothIndex];

	for (Int_t jj=ii+1; jj<TauIndexVector.size(); jj++) {
	  Int_t SecondTauIndex = TauIndexVector[jj];
	  Int_t SecMothIndex = mothMc[SecondTauIndex];
	  Int_t SecGrandMothIndex = mothMc[SecMothIndex];
	  Int_t SecGreatGrandMothIndex = mothMc[SecGrandMothIndex];

	  // check that they come from same SUSY cascade
	  if ( fabs(idMc[MothIndex])!=15 && fabs(idMc[SecMothIndex])!=15 ) {
	    if (MothIndex == SecMothIndex) continue;
	    else if (GrandMothIndex == SecGrandMothIndex) continue;
	    else if (MothIndex == SecGrandMothIndex || GrandMothIndex == MothIndex) {
	      SusyIndexVector.push_back(TauIndex);
	      SusyIndexVector.push_back(SecondTauIndex);
	    }	else continue;

	  } else if ( fabs(idMc[MothIndex])!=15 && fabs(idMc[SecMothIndex])==15 ) {
	    if (MothIndex == SecGrandMothIndex) continue;
	    else if (GrandMothIndex == SecGreatGrandMothIndex) continue;
	    else if (MothIndex == SecGreatGrandMothIndex || GrandMothIndex == SecGrandMothIndex) {
	      SusyIndexVector.push_back(TauIndex);
	      SusyIndexVector.push_back(SecondTauIndex);
	    } else continue;

	  } else if ( fabs(idMc[MothIndex])==15 && fabs(idMc[SecMothIndex])!=15 ) {
	    if (GrandMothIndex == SecMothIndex) continue;
	    else if (GreatGrandMothIndex == SecGrandMothIndex) continue;
	    else if (GrandMothIndex == SecGrandMothIndex || GreatGrandMothIndex == MothIndex) {
	      SusyIndexVector.push_back(TauIndex);
	      SusyIndexVector.push_back(SecondTauIndex);
	    } else continue;

	  } else if ( fabs(idMc[MothIndex])==15 && fabs(idMc[SecMothIndex])==15 ) {
	    if (GrandMothIndex == SecGrandMothIndex) continue;
	    else if (GreatGrandMothIndex == SecGreatGrandMothIndex) continue;
	    else if (GrandMothIndex == SecGreatGrandMothIndex || GreatGrandMothIndex == GrandMothIndex) {
	      SusyIndexVector.push_back(TauIndex);
	      SusyIndexVector.push_back(SecondTauIndex);
	    }
	  } else std::cout << "what am I doing ? " << endl;

	}  // end loop over list of Taus
	  
      } // end loop over list of Taus
    }

    //VIR
    TauHadTyp = -1;
    if (isDetailedHad==1){
      //std::cout << "did you want the HAD breakdown? ? " << endl;  
      //      for (Int_t imc = 0; imc < nMc; ++imc) {
      for (Int_t imc = 0; imc < SusyIndexVector.size(); ++imc) {

	Int_t tauIndex=SusyIndexVector[imc];

	//if (TMath::Abs(idMc[imc]) ==15 ) {
	//std::cout<<" yeahhhhh a tau "<< std::endl;

	// tau -> 2-body decays 
	if (nDauMc[imc]==2) {
	  //std::cout<<" 2 figlie ---------- "<< std::endl;

	  for (Int_t p = 0; p < nMc; ++p) {
	    // tau -> a1 nu.
	    Int_t pmothIndex=mothMc[p];
                                                                              
	    if ( TMath::Abs(idMc[p]) == 20213 && fabs(idMc[pmothIndex]==15)) {

	      Int_t a1motherIndex = mothMc[p];
	      // 		  std::cout<<" dentro a a1   "<< std::endl;
	      for(Int_t r = 0; r < nMc; ++r) {
		if( ((mothMc[r]) == a1motherIndex) && (TMath::Abs(idMc[r]) == 113 ) ) {
		  // 		      std::cout<<" dentro a a1 rho0  "<< std::endl;
		  TauHadTyp=11;
		}
		else if ( ((mothMc[r]) ==a1motherIndex) &&  (TMath::Abs(idMc[r]) == 213)) {
		  // 		      std::cout<<" dentro a a1 rho+ "<< std::endl;
		  TauHadTyp = 12; // 1 prong ...
		}
	      } // i am a rho0 or rho+
	    }// I am a a1

	    // tau -> pi nu.
	    else if (TMath::Abs(idMc[p]) == 211 && fabs(idMc[pmothIndex])==15 && tauIndex==mothMc[p]) {
	      // 		  std::cout<<" dentro a pi nu "<< std::endl;
	      TauHadTyp=13;
	    }
	    // tau -> rho+- nu.
	    else if ( TMath::Abs(idMc[p]) == 213  && fabs(idMc[pmothIndex])==15 && tauIndex==mothMc[p]){
	      // 		  std::cout<<" dentro a rho nu "<< std::endl;
	      TauHadTyp=14;
	    }	    
	  } // p
	} ////////   tau -> 2-body
	     
	// tau -> 3-body decays 
	else if (nDauMc[imc]==3) {
	  for (Int_t lep = 0; lep < nMc; ++lep) {
	    Int_t lepmothIndex=mothMc[lep];

	    // tau -> e nu nu
	    if ( TMath::Abs(idMc[lep]) == 11 &&  fabs(idMc[lepmothIndex])==15 && tauIndex==mothMc[lep]){
	      // 		if ( (TMath::Abs(idMc[lep]) == 11 ||  TMath::Abs(idMc[lep]) == 12 || TMath::Abs(idMc[lep]) == 16) && fabs(idMc[lepmothIndex])==15 && tauIndex==mothMc[lep]){
	      // 		  std::cout<<" dentro a e nu nu "<< std::endl;
	      TauHadTyp=15;
	    }
	    // tau -> mu nu nu
	    if ( TMath::Abs(idMc[lep]) == 13 && fabs(idMc[lepmothIndex])==15 && tauIndex==mothMc[lep]) {
	      // 		if ( (TMath::Abs(idMc[lep]) == 13 || TMath::Abs(idMc[lep]) == 14 || TMath::Abs(idMc[lep]) == 16) && fabs(idMc[lepmothIndex])==15 && tauIndex==mothMc[lep]) {
	      // 		  std::cout<<" dentro a mu nu nu "<< std::endl;
	      TauHadTyp=16;
	    }
	  } //lep
	} ////////   tau -> 3-body

	// tau -> 5-body decays 
	else if( nDauMc[imc]==5 ){
	  // tau -> pi pi pi  ... il quinto figlio qualunque (pi0 e' il piu' freq.)
	  Int_t primarypion = 0;
	  // for (Int_t p = 0; p < 50; ++p) {
	  for (Int_t p = 0; p < nMc; ++p) {

	    Int_t pmothIndex=mothMc[p];

	    if (TMath::Abs(idMc[p]) == 211 && fabs(idMc[pmothIndex])==15 && tauIndex==mothMc[p]){
	      primarypion++;
	    }
	    if ( primarypion == 3) {
	      // 		  std::cout<<" dentro a pi pi pi pi0 nu  "<< std::endl;
	      TauHadTyp = 17;
	    }
	    //	    std::cout<<"TauHadTyp=  "<< TauHadTyp << std::endl;
	  }//p
	}  ////////   tau -> 5-body
	//}//tau=15
      } // loop imc< nMc
	  
      //}// isDetailedHad
     
      if (TauHadTyp>-1) {
	//nCascade++;
	//2 bodies

	if (TauHadTyp==15)   TauHadCounter[0]++;
	if (TauHadTyp==16)   TauHadCounter[1]++;
	if (TauHadTyp==13)   TauHadCounter[2]++;
	if (TauHadTyp==11)   TauHadCounter[3]++;
	if (TauHadTyp==12)   TauHadCounter[4]++;
	if (TauHadTyp==14)   TauHadCounter[5]++;
	if (TauHadTyp==17)   TauHadCounter[6]++;

      }
       

      }// isDetailedHad
     

    if (SusyIndexVector.size()>0) {
      nSUSYTaus++;

      // Loop over Muons
      if (MuonIndexVector.size()>0) {
	for (Int_t mm=0; mm<MuonIndexVector.size(); mm++) {
	  Int_t MuonIndex = MuonIndexVector[mm];
	  Int_t MuonMothIndex = mothMc[MuonIndex];
	  Int_t MuonGrandMothIndex = mothMc[MuonMothIndex];
	  for (Int_t tt=0; tt<SusyIndexVector.size(); tt++) {
	    Int_t SusyTauIndex = SusyIndexVector[tt];
	    if ( MuonMothIndex==SusyTauIndex || (fabs(idMc[MuonMothIndex])==13 && MuonGrandMothIndex==SusyTauIndex) )
	      SusyMuonVector.push_back(mm);
	  }
	}
      }

      // Loop over Electrons
      if (ElectronIndexVector.size()>0) {
	for (Int_t ee=0; ee<ElectronIndexVector.size(); ee++) {
	  Int_t ElectronIndex = ElectronIndexVector[ee];
	  Int_t ElectronMothIndex = mothMc[ElectronIndex];
	  Int_t ElectronGrandMothIndex = mothMc[ElectronMothIndex];
	  for (Int_t tt=0; tt<SusyIndexVector.size(); tt++) {
	    Int_t SusyTauIndex = SusyIndexVector[tt];
	    if ( ElectronMothIndex==SusyTauIndex || (fabs(idMc[ElectronMothIndex])==11 && ElectronGrandMothIndex==SusyTauIndex) )
	      SusyElectronVector.push_back(ee);
	  }
	}
      }

      // 	  if (ElectronIndexVector.size()>0) {
      // 	    for (Int_t tt=0; tt<SusyIndexVector.size(); tt++) {
      // 	      Int_t SusyTauIndex = SusyIndexVector[tt];
      // 	      for (Int_t ee=0; ee<ElectronIndexVector.size(); ee++) {
      // 		Int_t ElectronIndex = ElectronIndexVector[ee];
      // 		Int_t ElectronMothIndex = mothMc[ElectronIndex];
      // 		Int_t ElectronGrandMothIndex = mothMc[ElectronMothIndex];	 
      // 		  if ( ElectronMothIndex==SusyTauIndex || (fabs(idMc[ElectronMothIndex])==11 && ElectronGrandMothIndex==SusyTauIndex) ) {
      // 		  SusyElectronVector.push_back(ee);
		  
      // 		  if (ElectronIndexVector.size()>1) {

      // 		    for (Int_t ee2=ee+1; ee2<ElectronIndexVector.size(); ee2++) {
      // 		      Int_t ElectronIndex2 = ElectronIndexVector[ee2];
      // 		      Int_t ElectronMothIndex2 = mothMc[ElectronIndex2];
      // 		      Int_t ElectronGrandMothIndex2 = mothMc[ElectronMothIndex2];	 	   
		      
      // 		      std::cout << "'elettone1 " << ElectronIndex << " , elettrone2 " << ElectronIndex2 << " ,Madre1 " <<  ElectronMothIndex << " ,Madre2 " <<  ElectronMothIndex2 << endl;

      // 		      if ( ElectronMothIndex2==SusyTauIndex || (fabs(idMc[ElectronMothIndex2])==11 && ElectronGrandMothIndex2==SusyTauIndex) ) {
      // 			SusyElectronVector.clear();
      // 			std::cout << "Pulizie di Primavera " << endl;
      // 		      }
      // 		    }
      // 		    if (SusyElectronVector.size()==0) break;
      // 		  }
      // 		}
      // 	      }
      // 	    }
      // 	  }
	
	
      // 	  if (SusyElectronVector.size()>2 )
      // 	    std::cout << "We have " <<  SusyElectronVector.size() << "Electrons in our final state" << endl;
	  
      // Now classify the sub-cases for Muons and Electrons
      Int_t TauDecayType = 0;
      if (SusyMuonVector.size()==0 && SusyElectronVector.size()==0)
	TauDecayType = 1;
      else if (SusyMuonVector.size()>0 && SusyElectronVector.size()==0)
	TauDecayType = 2;
      else if (SusyMuonVector.size()==0 && SusyElectronVector.size()>0)
	TauDecayType = 3;
      else if (SusyMuonVector.size()>0 && SusyElectronVector.size()>0) {
	for (Int_t mmm=0; mmm<SusyMuonVector.size(); mmm++) {
	  Int_t MuonIndex = MuonIndexVector[mmm];
	  Int_t MuonMothIndex = mothMc[MuonIndex];
	  Int_t MuonGrandMothIndex = mothMc[MuonMothIndex];
	  for (Int_t eee=0; eee<SusyElectronVector.size(); eee++) {
	    Int_t ElectronIndex = ElectronIndexVector[eee];
	    Int_t ElectronMothIndex = mothMc[ElectronIndex];
	    Int_t ElectronGrandMothIndex = mothMc[ElectronMothIndex];	     
	    if (MuonMothIndex==ElectronMothIndex || MuonMothIndex==ElectronGrandMothIndex || MuonGrandMothIndex==ElectronMothIndex || MuonGrandMothIndex==ElectronGrandMothIndex) {
	      continue;
	    } else TauDecayType = 4;
	  }	   
	}
      }



      SusyMuonVector.clear();
      SusyElectronVector.clear();

      // std::cout << "Event number " << jentry << "\t and Tau Decay Type " << TauDecayType << "\t" << endl;

      if (TauDecayType>0) {
	nCascade++;
	if (TauDecayType==1) TauCounter[0]++;
	if (TauDecayType==2) TauCounter[1]++;
	if (TauDecayType==3) TauCounter[2]++;
	if (TauDecayType==4) TauCounter[3]++;
      }

    } // end condition on SusyIndexVector size

       
    TauIndexVector.clear();
    SusyIndexVector.clear();
    MuonIndexVector.clear();
    ElectronIndexVector.clear();


  } // end of loop over events

  
  if (isDetailedHad==1){
    
    std::cout << "-------------------------Hadronic breakdown-----------------------------------" << endl;
    std::cout << "There are " << TauHadCounter[0] << " events with tau -> e nu nu" << endl;
    std::cout << "There are " << TauHadCounter[1] << " events with t tau -> mu nu nu" << endl;
    std::cout << "There are " << TauHadCounter[2] << " events with tau -> pi pi pi .. pi0" << endl;
    std::cout << "There are " << TauHadCounter[3] << " events with tau -> a1 rho" << endl;
    std::cout << "There are " << TauHadCounter[4] << " events with tau -> a1 rho" << endl;
    std::cout << "There are " << TauHadCounter[5] << " events with tau -> pi nu" << endl;
    std::cout << "There are " << TauHadCounter[6] << " events with tau -> rho+- nu" << endl;
    std::cout << "-----------------------------------------------------------------------------------------" << endl;
    
  } // isDetailedHad
  
  std::cout << "There are " << nTauEvents << " events with SUSY Taus, but only " << nSUSYTaus << " are of the same cascade" << std::endl;   
  std::cout << "There are " << TauCounter[0] << " events with Hadronic Tau decays" << endl;
  std::cout << "There are " << TauCounter[1] << " events with 1 Muon" << endl;
  std::cout << "There are " << TauCounter[2] << " events with 1 Electron" << endl;
  std::cout << "There are " << TauCounter[3] << " events with 1 Muon and 1 Electron" << endl;
  std::cout << "There are " << nCascade << " events with a complete SUSY cascade" << endl;
  //TauCounter->Fill();
  //TauHadCounter->Fill();
  
}
