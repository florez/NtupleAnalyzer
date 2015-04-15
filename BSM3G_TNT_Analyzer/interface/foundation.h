#ifndef FOUNDATION_H
#define FOUNDATION_H

//ROOT-based includes
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TProfile.h"

//C++-based includes
#include <vector>
#include <string>

//include physics objects
#include "Objects.h"

//collection of event information
struct MyEventCollection{
	std::string label;
	bool goodVertex;
	bool passedTrigger;
	double eventweight;
	
	//selected objects collections
	vector <muon_s*> muon;
	
	//initialisation
	MyEventCollection(const std::string &inputlabel){
		label 			= inputlabel;
		goodVertex 		= false;
		passedTrigger 		= false;
		double eventweight	= 0.;
	}
	
	//destruction
	void clear(){
		goodVertex	=false;
		passedTrigger	=false;
		eventweight	=0.;
		muon.clear();
	}
};

//collection of histograms
struct MyHistoCollection{
	std::string label;
	TFile *output;
	
	//cutflow
	TH1F *h_count;
	
	//muon histograms
	TH1F *h_muon1pT;
	TH1F *h_muon2pT;
	TH1F *h_muon1eta;
	TH1F *h_muon2eta;
	
	//initialisation
	MyHistoCollection(TFile *_f, const std::string &inputlabel){
		label	=inputlabel;
		output=_f;
		
		output->mkdir(inputlabel.c_str());
		output->cd(inputlabel.c_str());
		
		//cutflow
		h_count	= new TH1F("counts","Cutflow",1,0,1);
		h_count->SetBit(TH1::kCanRebin);
		h_count->SetStats(0);
		
		//muon histograms
		h_muon1pT = new TH1F("h_muon1pT","Leading #mu p_{T}", 99, 10, 1000);
		h_muon1pT->GetXaxis()->SetTitle("p_{T}^{#mu1} [GeV]");
		h_muon1pT->GetYaxis()->SetTitle("events");
		
		h_muon2pT = new TH1F("h_muon2pT","Second to leading #mu p_{T}", 99, 10, 1000);
		h_muon2pT->GetXaxis()->SetTitle("p_{T}^{#mu2} [GeV]");
		h_muon2pT->GetYaxis()->SetTitle("events");
		
		h_muon1eta = new TH1F("h_muon1eta","Leading #mu #eta", 60, -3, 3);
		h_muon1eta->GetXaxis()->SetTitle("#eta^{#mu1} [GeV]");
		h_muon1eta->GetYaxis()->SetTitle("events");
		
		h_muon2eta = new TH1F("h_muon2eta","Second to leading #mu #eta", 60, -3, 3);
		h_muon2eta->GetXaxis()->SetTitle("#eta^{#mu2} [GeV]");
		h_muon2eta->GetYaxis()->SetTitle("events");				
	}
	
	//generate a filesave routine
	void Save(){
		output->cd(label.c_str());
		h_count->Write();
		h_muon1pT->Write();
		h_muon2pT->Write();
		h_muon1eta->Write();
		h_muon2eta->Write();
	}
};

//collection of profiles for uncertainty calculation
struct MyProfileCollection{
	std::string label;
	TFile *output;
	
	//cutflow
	TProfile *p_count;
	
	//muon histograms
	TProfile *p_muon1pT;
	TProfile *p_muon2pT;
	TProfile *p_muon1eta;
	TProfile *p_muon2eta;
	
	//initialisation
	MyProfileCollection(TFile *_f, const std::string &inputlabel){
		label	=inputlabel;
		output=_f;
		
		output->mkdir(inputlabel.c_str());
		output->cd(inputlabel.c_str());
		
		//cutflow
		p_count	= new TProfile("counts","Cutflow",1,0,1);
		p_count->SetBit(TH1::kCanRebin);
		p_count->SetStats(0);
		
		//muon histograms
		p_muon1pT = new TProfile("p_muon1pT","Leading #mu p_{T}", 99, 10, 1000);
		p_muon1pT->GetXaxis()->SetTitle("p_{T}^{#mu1} [GeV]");
		p_muon1pT->GetYaxis()->SetTitle("events");
		
		p_muon2pT = new TProfile("p_muon2pT","Second to leading #mu p_{T}", 99, 10, 1000);
		p_muon2pT->GetXaxis()->SetTitle("p_{T}^{#mu2} [GeV]");
		p_muon2pT->GetYaxis()->SetTitle("events");
		
		p_muon1eta = new TProfile("p_muon1eta","Leading #mu #eta", 60, -3, 3);
		p_muon1eta->GetXaxis()->SetTitle("#eta^{#mu1} [GeV]");
		p_muon1eta->GetYaxis()->SetTitle("events");
		
		p_muon2eta = new TProfile("p_muon2eta","Second to leading #mu #eta", 60, -3, 3);
		p_muon2eta->GetXaxis()->SetTitle("#eta^{#mu2} [GeV]");
		p_muon2eta->GetYaxis()->SetTitle("events");				
	}
	
	//generate a filesave routine
	void Save(){
		output->cd(label.c_str());
		p_count->Write();
		p_muon1pT->Write();
		p_muon2pT->Write();
		p_muon1eta->Write();
		p_muon2eta->Write();
	}	
};

#endif
