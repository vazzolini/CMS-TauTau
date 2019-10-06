#define SUSYTauTauAnalysis_cxx
#include "SUSYTauTauAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void SUSYTauTauAnalysis::Loop()
{
  //   In a ROOT session, you can do:
  //      Root > .L SUSYTauTauAnalysis.C
  //      Root > > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();


  TFile *file=new TFile("pippo.root","RECREATE"); 
  
  TH1F *TauType=new TH1F("TauType", "Type of the Tau", 200, 0 , 20);
   
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;


    /////
    //    Int_t ipB[] = {-1, -1};
    Int_t ipB[] = {1, 1};
    Int_t typB[] = {7, 7};
    Int_t cnt(-1);

    fB1Index = fB2Index  = -99;


    // -- Find indices of X20's
    for (Int_t imc = 0; imc < nMc; ++imc) {    // Loop on Mc-truth particles

      if ((TMath::Abs(idMc[imc]) == 1000023) || (TMath::Abs(idMc[imc]) == 1000015) || (TMath::Abs(idMc[imc]) == 2000015)) {
	cnt++;
	ipB[cnt] = imc;
	std::cout<<"cnt= "<< cnt << "ipB[cnt]="<< ipB[cnt]<< std::endl;
	if (fX1Index < 0) {
	  fX1Index = imc; 
	} else {
	  fX2Index = imc; 
	}
      }
      if (cnt == 1) break;
    }   // loop nMc
    
    Int_t ib;
    fTauTyp = -1;

    // Meaning of fTauTyp:
    // -- Determine special modes of sl (or leptonic taunu) B decay: 
    //    1 UNUSED
    //    2 UNUSED
    //    3 UNUSED
    //    4 UNUSED
    //    5 TAU NU, TAU->E
    //    6 TAU NU, TAU->MU
    //    7 TAU NU, TAU->PI
    //    8 TAU NU, TAU->A1 (RHO0)
    //    9 TAU NU, TAU->A1 (RHO+-)
    //   10 TAU NU, TAU->RHO+-
    //   11 TAU NU, TAU->3 PI (+1 more particle)
    //  100 TAU NU, TAU->ALTRO
  
  
    for (ib = 0; ib < 2; ib++) {
      std::cout<<"ib "<< ib << "----ipB[ib]="<< ipB[ib]<< "----nDauMc[ipB[ib]]== "<< nDauMc[ipB[ib]]<< std::endl;

      if (nDauMc[ipB[ib]]==2) {
	for (Int_t imc = 0; imc < nMc; ++imc) {
	  if (isTruTau(imc)) {
	    if ((mothMc[imc]-1) == ipB[ib]) {
	      // tau -> 2-body decays 
	      if (nDauMc[imc]==2) {
		for (Int_t p = 0; p < nMc; ++p) {
		  // tau -> a1 nu.
		  if ( ((mothMc[p]-1) == imc) &&  (TMath::Abs(idMc[p]) == 20213)) {
		    for(Int_t r = 0; r < nMc; ++r) {
		      if( ((mothMc[r]-1) == p) && (TMath::Abs(idMc[r]) == 113 ) ) {
			fTauTyp=8;
		      }
		      else if ( ((mothMc[r]-1) == p) &&  (TMath::Abs(idMc[r]) == 213)) {
			fTauTyp = 9; // 1 prong ...
		      }
		    }
		  }
		  // tau -> pi nu.
		  else if ( ((mothMc[p]-1) == imc) &&  (TMath::Abs(idMc[p]) == 211)) {
		    fTauTyp=7;
		  }
		  // tau -> rho+- nu.
		  else if ( ((mothMc[p]-1) == imc) &&  (TMath::Abs(idMc[p]) == 213))
		    fTauTyp=10;
		}
	      }
	      // tau -> 3-body decays 
	      else if (nDauMc[imc]==3) {
		for (Int_t lep = 0; lep < nMc; ++lep) {
		  // tau -> e nu nu
		  if ( ((mothMc[lep]-1) == imc) &&  (TMath::Abs(idMc[lep]) == 11)) {
		    fTauTyp=5;
		  }
		  // tau -> mu nu nu
		  if ( ((mothMc[lep]-1) == imc) &&  (TMath::Abs(idMc[lep]) == 13)) {
		    fTauTyp=6;
		  }
		} 
	      }
	      // tau -> 5-body decays 
	      else if( (nDauMc[imc]==5) ){
		// tau -> pi pi pi  ... il quinto figlio qualunque (pi0 e' il piu' freq.)
		Int_t primarypion = 0;
		for (Int_t p = 0; p < nMc; ++p) {
		  if ( ((mothMc[p]-1) == imc) && (TMath::Abs(idMc[p]) == 211) )
		    primarypion++;
		}
		if ( primarypion == 3) {
		  fTauTyp = 11;
		}
		std::cout<<"fTauTyp = "<< fTauTyp << std::endl;
	      }
	    } // loop moth
	  }
	} // loop nMc
      } // ndau=2 
    } // ib
    //  }

    //std::cout<<"i -particle ="<<i<<"     nMc ="<<nMc<< std::endl;
    //std::cout<<"abs(idMc[i]) =" << abs(idMc[i])<<std::endl;
      
     
  } // jentry
  file->Write();
  delete.file;

} //Loop



// ----------------------------------------------------------------------
void  SUSYTauTauAnalysis::bookHist(int dump) {

  cout << " entering recoilTaunu bookHist " << endl;

  char name[100], title[100];

  if (!fHistFile) {
    cout << "Call recoilTaunu::openHistFile(...) before booking histograms" << endl;
  } else {
    fHistFile->cd();
  }    

  if (dump > 0) {
    cout << "Booking events tree" << endl;
    fDump = dump;
    fTree = new TTree("events", "events"); 

   // MC info
    fTree->Branch("MTauTyp",    &fTauTyp, "MTauTyp/I");
  }
  
  fHistFile->cd();
}


