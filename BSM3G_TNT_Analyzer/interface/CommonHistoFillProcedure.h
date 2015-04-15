#ifndef COMMONHISTOFILLPROCEDURE_H
#define COMMONHISTOFILLPROCEDURE_H

#include "foundation.h"

void fillHistoCollection(MyHistoCollection &inputHistoCollection, MyEventCollection inputEventCollection, double weight){

  //fill histograms
  bool verbose =false;
  
  if(verbose)std::cout<<"start filling histograms"<<std::endl;
  
  //Muon histograms
  if(inputEventCollection.muon.size()>0){
    inputHistoCollection.h_muon1pT->Fill(inputEventCollection.muon[0]->pt,weight);
    inputHistoCollection.h_muon1eta->Fill(inputEventCollection.muon[0]->eta,weight);
    if(inputEventCollection.muon.size()>1){
      inputHistoCollection.h_muon2pT->Fill(inputEventCollection.muon[1]->pt,weight);
      inputHistoCollection.h_muon2eta->Fill(inputEventCollection.muon[1]->eta,weight);
    }
  }
  
  if(verbose)std::cout<<"filling histograms successful"<<std::endl;
}

void fillProfileCollection(MyProfileCollection &inputProfileCollection, MyEventCollection inputEventCollection, double pass){

  //fill histograms
  bool verbose =false;
  
  if(verbose)std::cout<<"start filling profiles"<<std::endl;
  
  //Muon histograms
  if(inputEventCollection.muon.size()>0){
    inputProfileCollection.p_muon1pT->Fill(inputEventCollection.muon[0]->pt,pass);
    inputProfileCollection.p_muon1eta->Fill(inputEventCollection.muon[0]->eta,pass);
    if(inputEventCollection.muon.size()>1){
      inputProfileCollection.p_muon2pT->Fill(inputEventCollection.muon[1]->pt,pass);
      inputProfileCollection.p_muon2eta->Fill(inputEventCollection.muon[1]->eta,pass);
    }
  }
  
  if(verbose)std::cout<<"filling profiles successful"<<std::endl;
}

#endif
