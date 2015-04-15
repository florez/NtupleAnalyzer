#ifndef OBJECTS_H
#define OBJECTS_H

#include <iostream>

using namespace std;

struct muon_s{
	double pt;
	double eta;
	
	muon_s(double pT_, double eta_){
		pt	=pT_;
		eta	=eta_;
	}
};

#endif
