#ifndef SELECTION_H
#define SELECTION_H

#include "Objects.h"
#include "TString.h"
#include "CommonHistoFillProcedure.h"

struct Selection{
	//global input
	std::string label;
	MyEventCollection *InputCollection;
	MyHistoCollection *OutputCollection;
	MyProfileCollection *ProfileCollection;
	double weight;
	  	
	//Muon requirements
	float MuonEtaMax;
	
	//inversion arguments
	bool invertMuonRequirements;
	
	Selection(const std::string &inputlabel){
	  label = inputlabel;
	  MuonEtaMax=1.3;
	  invertMuonRequirements = false;
	  weight=0.;
	}
	
	void select(bool verbose){
	  //initalize the cutflow
	  if(verbose)std::cout<<"__________initialize the cutflow histogram__________"<<endl;
	  (*OutputCollection).h_count->Fill("No Cuts",0); 
	  if(invertMuonRequirements)(*OutputCollection).h_count->Fill(TString::Format("#eta^{#mu}>%.1f",MuonEtaMax),0);
	  else (*OutputCollection).h_count->Fill(TString::Format("#eta^{#mu}<=%.1f",MuonEtaMax),0);
	  
	  //Muon Requirements
	  if(verbose)std::cout<<"__________starting muon selection__________"<<endl;
	  if((*InputCollection).muon.size()>1){
	    if(fabs((*InputCollection).muon[0]->eta)>MuonEtaMax || fabs((*InputCollection).muon[1]->eta)>MuonEtaMax){
	    	if(invertMuonRequirements){
			(*OutputCollection).h_count->Fill(TString::Format("#eta^{#mu}>%.1f",MuonEtaMax),weight);
			fillHistoCollection((*OutputCollection), (*InputCollection), weight);
			fillProfileCollection((*ProfileCollection), (*InputCollection), 1.);
			return;
		}
		else{
			fillProfileCollection((*ProfileCollection), (*InputCollection), 0.);
		}
	    }
	  }
	  
	  //everything done, keep passed events for selection
	  if(verbose)std::cout<<"__________all selection requirements passed__________"<<endl;
	  
	  if(!invertMuonRequirements){//check for any of the inversions to be specified
	  	fillHistoCollection((*OutputCollection), (*InputCollection), weight);
		fillProfileCollection((*ProfileCollection), (*InputCollection), 1.);	
	  }
	  else	fillProfileCollection((*ProfileCollection), (*InputCollection), 0.);
	}
};

#endif
