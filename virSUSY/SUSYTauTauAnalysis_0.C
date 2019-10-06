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

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    Int_t* bdau=new Int_t[30];  //B dau's




    for (Int_t i=0;i<nMc;i++){      // Loop on Mc-truth particles

      //std::cout<<"i -particle ="<<i<<"     nMc ="<<nMc<< std::endl;
      //std::cout<<"abs(idMc[i]) =" << abs(idMc[i])<<std::endl;
      
      //p=2212


      //loop over daughters
      bool tauhadtauhad = 0;
      bool taueXtauhad = 0;
      bool taumuXtauhad = 0;


      /////////////////////////////////////////////////////////////////////////////////////// a la btd 
	//	for (Int_t k=0;k<30;k++) bdau[k]=-9;
	////	if(idMc[i]==1000023){
	//	if(idMc[i]==15){

	//  std::cout<<"I am a  ="<< abs(idMc[i]) << std::endl;
	//  //	std::cout<<"num of daugthers ="<< nDauMc[i] << std::endl;
	//  //	std::cout<<" mother  ="<< mothMC[i] << std::endl;      
	
	  ////X2 0 daughters
	//  giveMeDau(i,bdau);
	//	  //  if (X20dau[0]==0) continue;
	//  if (bdau[0]==0) continue;


	//  std::cout<<"bdau[0]=" << bdau[0] << std::endl; 
	//  //std::cout<<"bdau[0]=" << bdau[0] << "  idMc[bdau[1]]=" << idMc[bdau[1]]<< "  idMc[bdau[2]]=="<< abs(idMc[bdau[2]]) <<  "   idMc[bdau[3]])==" << abs(idMc[bdau[3]])<< std::endl;

	//  if (bdau[0]==3 &&  abs(idMc[bdau[1]])==13 && abs(idMc[bdau[2]])==12 && abs(idMc[bdau[3]])==18) { 
	//    std::cout<<"bdau..."<< std::endl; 
	//    printInfo(bdau);}

	//	//if (bdau[0]==2 &&  abs(idMc[bdau[1]])==421 && abs(idMc[bdau[2]])==211) {cout<<"bdau..."<<endl; printInfo(bdau);}
	//	} //loop idMc   
	////////////////////////////////////////////////////////////////////////////////////// /. a la btd 

	//	if(idMc[i]==1000023){
	if(idMc[i]==15){
	  std::cout<<"I'm a tau ="<< abs(idMc[i]) << "  with num of daugthers ="<< nDauMc[i] << std::endl;
	}
      
	int Xid = abs(idMc[i]);
	if(Xid!=1000023)continue;
      
	if (nDauMc[i] == 3) {
	  Int_t id[3]={0,0};
	  Int_t counter=0;
	  for (Int_t j=0; j<nMc; j++) {
	    Int_t imother = mothMc[j]-1;
	    if (imother == i) {
	      id[counter]=idMc[j];
	      counter++;
	    }
	    if (counter == 3) break;
	  }
	  bool X2stau=0;
	  bool X2stautau=0;

	  if (TMath::Abs(id[0])==100023){ 
	    //std::cout<<"id[0])== "<<TMath::Abs(id[0]) << " << id[1]== " <<TMath::Abs(id[1])<< std::endl;
	  }

	  if ((TMath::Abs(id[0])==100023) 
	      && (TMath::Abs(id[1])==1000015 || TMath::Abs(id[1])==2000015)) 
	    X2stau=1;
	  //std::cout<<" X2stau=" << X2stau << "  Abs(id[1])==" <<TMath::Abs(id[1])<< std::endl;
	  if ((TMath::Abs(id[0])==100023) 
	      && (TMath::Abs(id[1])==1000015 || TMath::Abs(id[1])==2000015)
	      && (TMath::Abs(id[1])==15)) 
	    X2stautau=1;
	  //std::cout<<"X2statau "<< X2stautau << std::endl;
	}  // loop over ndau
 

	//      for(int j=0; j < nDauMc[i]; j++){
	//	int Dau_id = idMc[i][j];
	//	int absDau_id = abs(Dau_id);
	//if(idMc[i]==1000023){std::cout<< Xid << "  &&  " << Dau_id << std::endl; }	
	//	if( Xid == 1000023 && absDau_id == 443 ) tauhadtauhad  =1;
	//if( Xid == 1000023 && absDau_id == 443 ) taueXtauhad   =1;
	//if( Xid == 1000023 && absDau_id == 443 ) taumuXtauhad  =1;
	//      } // loop over ndau
      
      



    
 
      
      
      
    } // loop nMc

  } // jentry
} //Loop

void SUSYTauTauAnalysis::giveMeDau(Int_t index,Int_t* dauarray){
  Int_t ndau=0;
  //std::cout<<"index ="<< index<< "     dauarray = " <<  dauarray  << std::endl;
  std::cout<<"in giveMeDau"  << std::endl;
  if (index>=0) {
    for (Int_t j=0;j<nMc;j++){
      if ((mothMc[j]-1)==index){ 
     	dauarray[++ndau]=j; 
      }
    }
  }
  dauarray[0]=ndau;
  std::cout<<"ndauarray[0] =  "   << ndau<<std::endl;
  //  // ordering from higher to lower abs(lundId)
  for (Int_t j=1;j<=ndau;j++) {
    Int_t tmp=dauarray[j];
    Int_t j_max=j;
    std::cout<<"in giveMeDau  2  "   << "tmp= " << dauarray[j] < < "  j_max=  " << j<<  std::endl;

    for (Int_t k=j+1;k<=ndau;k++) if (abs(idMc[dauarray[j_max]])<abs(idMc[dauarray[k]])) j_max=k;
    dauarray[j]=dauarray[j_max];
    dauarray[j_max]=tmp;
  }
  return;
} //end giveMeDau 

void SUSYTauTauAnalysis::printInfo(Int_t* array){
  std::cout<<" in printINFO !!!!!!!!!!!" << std::endl;
  std::cout<<"array[0] = "<<array[0]<<std::endl;
  for (Int_t i=1;i<=array[0];i++) {
    std::cout<<"array["<<i<<"] = "<<array[i]<<" idMc[array["<<i<<"]] = "<<idMc[array[i]]<<std::endl;
  }
}//end printInfo
