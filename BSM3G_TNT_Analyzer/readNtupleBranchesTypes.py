##############################################################
# Author: Andres Florez, Universidad de los Andes, Colombia. #
# Author: Denis Rathjens, University of Hamburg, Germany     #
##############################################################

from optparse import OptionParser

#Parser for command line arguments
parser = OptionParser()
parser.add_option("--inFile", dest="inFile", default="OutTree.root", help="filepath for input ntuple file, default=%default")
parser.add_option("--AnaName",dest="AnaName",default="BSM3G_TNT_Analyzer", help="name of the automatically generated analyzer, default=%default")
parser.add_option("--dir", dest="dir", default="TNT", help="ntuple path to the tree, default=%default")
parser.add_option("--tree", dest="tree", default="BOOM", help="ntuple tree name, default=%default")
(options, args) = parser.parse_args()

print "Input file: " + options.inFile

print "producing list of ntuple branches"

import os, sys

#os.system("python getBrancheslist.py --inFile "+options.inFile+" --dir "+options.dir+" --tree "+options.tree)
from ROOT import gROOT, TF1, TFile, TTree
import sys

gROOT.Reset()
f1 = TFile(options.inFile)
f1.cd(options.dir)
tree = f1.Get(options.dir+"/"+options.tree)
tree.MakeClass("intermediate")

log_h = open(options.AnaName+".h", "w")

print >> log_h, "////////////////////////////////////////////////////////////////////"
print >> log_h, "// Developer: Andres Florez, Universidad de los Andes, Colombia. //"
print >> log_h, "// Developer: Denis Rathjens, University of Hamburg, Germany     //"
print >> log_h, "//////////////////////////////////////////////////////////////////"
print >> log_h, ""
print >> log_h, ""
print >> log_h, "#ifndef "+options.AnaName+"_h" 
print >> log_h, "#define "+options.AnaName+"_h"
print >> log_h, ""
print >> log_h, "#include <iostream>"
print >> log_h, "#include <fstream>"
print >> log_h, "#include <iomanip>"
print >> log_h, "#include <string>"
print >> log_h, "#include <math.h>"
print >> log_h, "#include <stdio.h>"
print >> log_h, "#include <stdlib.h>"
print >> log_h, "#include <vector>"
print >> log_h, ""
print >> log_h, "#include <TROOT.h>"
print >> log_h, "#include <TFile.h>"
print >> log_h, "#include <TBranch.h>"
print >> log_h, "#include <TApplication.h>"
print >> log_h, "#include <TTree.h>"
print >> log_h, "#include <TChain.h>"
print >> log_h, "#include <TH1F.h>"
print >> log_h, ""
print >> log_h, "using namespace std;"
print >> log_h, ""
print >> log_h, "class "+options.AnaName+" {"
print >> log_h, "public :"
print >> log_h, "   "+options.AnaName+"();"
print >> log_h, "   ~"+options.AnaName+"();"
print >> log_h, ""
print >> log_h, "   void setBranchAddress(TTree* "+options.tree+");"

ReadLeaves=0

with open("intermediate.h") as f:
    for line in f:	
    	if ReadLeaves==1:
		if line.find("intermediate(TTree *tree=0);")!=-1:
			print >> log_h, "};"
			print >> log_h, "#endif"
			break
		print >> log_h, line.rstrip("\n")
    	if line.find("// Declaration of leaf types")!=-1:
		ReadLeaves=1
    
print "writing "+options.AnaName+".cc"

log_cc = open(options.AnaName+".cc", "w")

print >> log_cc, "#include \""+options.AnaName+".h\""
print >> log_cc, "#include \"interface/selection.h\""
print >> log_cc, ""
print >> log_cc, "int main (int argc, char *argv[])"
print >> log_cc, "{"
print >> log_cc, ""
print >> log_cc, "  TApplication app(\"App\",&argc, argv);"
print >> log_cc, "  "+options.AnaName+" "+options.AnaName+"_;"
print >> log_cc, ""
print >> log_cc, "}"
print >> log_cc, ""+options.AnaName+"::"+options.AnaName+"()"
print >> log_cc, "{"
print >> log_cc, ""
print >> log_cc, "  //define output file"
print >> log_cc, "  TFile *out=new TFile(\"Output.root\",\"RECREATE\");"
print >> log_cc, ""
print >> log_cc, "  //define entry count histogram"
print >> log_cc, "  TH1F *count=new TH1F(\"count\",\"\",1,0,1);"
print >> log_cc, ""
print >> log_cc, "  //define output collections"
print >> log_cc, "  MyHistoCollection   myHistoColl_TEST	(out, \"TEST\");"
print >> log_cc, "  MyProfileCollection myProfileColl_TEST(out, \"p_TEST\");"
print >> log_cc, ""
print >> log_cc, "  //define collections of selected objects per event"
print >> log_cc, "  MyEventCollection hard(\"hard\");"
print >> log_cc, "  MyEventCollection soft(\"soft\");"
print >> log_cc, ""
print >> log_cc, "  //load PU weights"
print >> log_cc, "  TFile file_PUdata(\"PUdata.root\",\"read\");"
print >> log_cc, "  TH1F *PUweights = (TH1F*)file_PUdata.Get(\"analyzeHiMassTau/NVertices_0\");"
print >> log_cc, "  PUweights->Scale(1/PUweights->Integral());"
print >> log_cc, "  TFile file_PUsim(\"PUsim.root\",\"read\");"
print >> log_cc, "  TH1F *PUsim = (TH1F*)file_PUsim.Get(\"analyzeHiMassTau/NVertices_0\");"
print >> log_cc, "  PUsim->Scale(1/PUsim->Integral());"
print >> log_cc, ""
print >> log_cc, "  PUweights->Divide(PUsim);"
print >> log_cc, ""
print >> log_cc, "  //configure input file"
print >> log_cc, "  TFile *f = new TFile (\"OutTree.root\");"
print >> log_cc, "  f->cd(\"TNT\");"
print >> log_cc, "  TTree* "+options.tree+" = (TTree*)f->Get(\""+options.dir+"/"+options.tree+"\");"
print >> log_cc, ""
print >> log_cc, "  int nentries = (int) "+options.tree+"->GetEntries();"
print >> log_cc, "  setBranchAddress("+options.tree+");"
print >> log_cc, ""  
print >> log_cc, "  for (int i = 0; i < nentries; ++i)"  
print >> log_cc, "    {"
print >> log_cc, "      "+options.tree+"->GetEntry(i);"
print >> log_cc, "	count->Fill(0.5);" 
print >> log_cc, "" 
print >> log_cc, "	//define global event weight"
print >> log_cc, "	double weight =1.;"
print >> log_cc, "	weight=PUweights->GetBinContent(PUweights->FindBin(trueInteractions));" 
print >> log_cc, ""
print >> log_cc, "       for (int j = 0; j < Muon_pt->size(); j++)"
print >> log_cc, "         {"
print >> log_cc, "	      muon_s element(Muon_pt->at(j),Muon_eta->at(j));"
print >> log_cc, "            if(Muon_pt->at(j)>30) hard.muon.push_back(&element);"
print >> log_cc, "            else soft.muon.push_back(&element);"
print >> log_cc, "         }"
print >> log_cc, "	//do some event selection"
print >> log_cc, "  	Selection test(\"test\");"
print >> log_cc, "  	test.InputCollection	= &hard;"
print >> log_cc, "  	test.OutputCollection	= &myHistoColl_TEST;"
print >> log_cc, "  	test.ProfileCollection	= &myProfileColl_TEST;"
print >> log_cc, "  	test.weight=weight;"
print >> log_cc, "  	test.MuonEtaMax=2.1;"
print >> log_cc, ""
print >> log_cc, "  	//the actual selection in a nutshell"
print >> log_cc, "  	test.select(false);"
print >> log_cc, ""
print >> log_cc, "  	//clear used input collections"
print >> log_cc, "  	hard.clear();"
print >> log_cc, "  	soft.clear();"
print >> log_cc, "     }"
print >> log_cc, "  //write histograms&profiles"
print >> log_cc, "  out->cd();"
print >> log_cc, "  count->Write();"
print >> log_cc, "  myHistoColl_TEST.Save();"
print >> log_cc, "  myProfileColl_TEST.Save();"
print >> log_cc, ""
print >> log_cc, "  //close files"
print >> log_cc, "  f->Close();"
print >> log_cc, "  out->Close();"
print >> log_cc, ""
print >> log_cc, "}"
print >> log_cc, ""
print >> log_cc, ""+options.AnaName+"::~"+options.AnaName+"()"
print >> log_cc, "{"
print >> log_cc, ""
print >> log_cc, "  // do anything here that needs to be done at desctruction time"
print >> log_cc, ""
print >> log_cc, "}"
print >> log_cc, ""
print >> log_cc, "void "+options.AnaName+"::setBranchAddress(TTree* "+options.tree+")"
print >> log_cc, "{"
print >> log_cc, ""

ReadBranches=0

with open("intermediate.h") as f:
	for line in f:	
    		if ReadBranches==1:
			if line.find("fChain = tree;")!=-1:
				#print >> log_cc, "   fChain = "+options.tree
				continue
			if line.find("fCurrent = -1")!=-1:
				continue
			if line.find("fChain->SetMakeClass(1);")!=-1:
				continue		
			if line.find("if (!tree) return;")!=-1:
				print >> log_cc, "   if(!"+options.tree+") return;"
				continue				
			if line.find("Notify();")!=-1:
				print >> log_cc, "};"
				break
			print >> log_cc, (line.replace("fChain", options.tree)).rstrip("\n")				
    		if line.find("// (once per file to be processed).")!=-1:
			ReadBranches=1
		
os.system("rm intermediate.C intermediate.h")		
		
print "All done! Run: \"make "+options.AnaName+" to compile the analyzer skeleton"		

