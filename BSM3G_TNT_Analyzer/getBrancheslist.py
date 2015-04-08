##############################################################
# Author: Andres Florez, Universidad de los Andes, Colombia. #
##############################################################

from ROOT import gROOT, TCanvas, TF1, TFile, TTree
import sys

gROOT.Reset()
f1 = TFile("OutTree.root")
f1.cd("TNT")
tree = f1.Get("TNT/BOOM")
tree.Print()
