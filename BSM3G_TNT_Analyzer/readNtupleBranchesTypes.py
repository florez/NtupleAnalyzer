##############################################################
# Author: Andres Florez, Universidad de los Andes, Colombia. #
##############################################################

from collections import deque

log_h = open("BSM3G_TNT_Analyzer.h", "w")

print >> log_h, "////////////////////////////////////////////////////////////////////"
print >> log_h, "// Developer: Andres Florez, Universidad de los Andes, Colombia. //"
print >> log_h, "//////////////////////////////////////////////////////////////////"
print >> log_h, ""
print >> log_h, ""
print >> log_h, "#ifndef BSM3G_TNT_Analyzer_h" 
print >> log_h, "#define BSM3G_TNT_Analyzer_h"
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
print >> log_h, ""
print >> log_h, "using namespace std;"
print >> log_h, ""
print >> log_h, "class BSM3G_TNT_Analyzer {"
print >> log_h, "public :"
print >> log_h, "   BSM3G_TNT_Analyzer();"
print >> log_h, "   ~BSM3G_TNT_Analyzer();"
print >> log_h, ""
print >> log_h, "   void setBranchAddress(TTree* Tree);"

queue   = deque()

with open("outputBranchFile") as f:
    for line in f:
        for word in line.split():
           word = word.replace(":", "")
           queue.append(word)


index = 3

for x in range(len(queue)):
    if x == index:
        getVector =  queue[index+2]

        if getVector.find('vector') != -1:
             print >> log_h, "  ", queue[index+2],"*", queue[index],";"
        if getVector.find('/D') != -1:
             print >> log_h, "   double ", queue[index],";" 
        if getVector.find('/F') != -1:
             print >> log_h, "   float ", queue[index],";"
        if getVector.find('/I') != -1:
             print >> log_h, "   int ", queue[index],";"
        
        index+=31

print >> log_h, ""

index_b = 3

print >> log_h, "};"
print >> log_h, "#endif"

log_cc = open("BSM3G_TNT_Analyzer.cc", "w")

print >> log_cc, "#include \"BSM3G_TNT_Analyzer.h\""
print >> log_cc, ""
print >> log_cc, "int main (int argc, char *argv[])"
print >> log_cc, "{"
print >> log_cc, ""
print >> log_cc, "  TApplication app(\"App\",&argc, argv);"
print >> log_cc, "  gROOT->ProcessLine(\"#include <vector>\");"
print >> log_cc, "  BSM3G_TNT_Analyzer BSM3G_TNT_analyzer;"
print >> log_cc, ""
print >> log_cc, "}"
print >> log_cc, "BSM3G_TNT_Analyzer::BSM3G_TNT_Analyzer()"
print >> log_cc, "{"
print >> log_cc, ""
print >> log_cc, "  TFile *f = new TFile (\"OutTree.root\");"
print >> log_cc, "  f->cd(\"TNT\");"
print >> log_cc, "  TTree* tree = (TTree*)f->Get(\"TNT/BOOM\");"
print >> log_cc, ""
print >> log_cc, "  int nentries = (int) tree->GetEntries();"
print >> log_cc, "  setBranchAddress(tree);"
print >> log_cc, ""  
print >> log_cc, "  for (int i = 0; i < nentries; ++i)"  
print >> log_cc, "    {"
print >> log_cc, "       tree->GetEntry(i);"
print >> log_cc, "       for (int j = 0; j < Muon_pt->size(); j++)"
print >> log_cc, "         {"
print >> log_cc, "            cout <<Muon_pt->at(j)<<endl;"
print >> log_cc, "         }"
print >> log_cc, "     }"
print >> log_cc, ""
print >> log_cc, "}"
print >> log_cc, ""
print >> log_cc, "BSM3G_TNT_Analyzer::~BSM3G_TNT_Analyzer()"
print >> log_cc, "{"
print >> log_cc, ""
print >> log_cc, "  // do anything here that needs to be done at desctruction time"
print >> log_cc, ""
print >> log_cc, "}"
print >> log_cc, ""
print >> log_cc, "void BSM3G_TNT_Analyzer::setBranchAddress(TTree* Tree)"
print >> log_cc, "{"
print >> log_cc, ""

index = 3

for x in range(len(queue)):
    if x == index:
        setbranches = "  Tree->SetBranchAddress(\""+queue[index]+"\", &"+queue[index]+");"
        print >> log_cc, setbranches 
        index+=31

print >> log_cc, "}"

