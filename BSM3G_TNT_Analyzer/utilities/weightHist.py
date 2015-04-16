#!/usr/bin/env python

#original script from Lukas Vanelderen at https://github.com/lveldere/VBF-LS-tau-tools/tree/master/drawHist
#modified by Denis Rathjens

import mytools,shutil
import sys,os,glob,math
import ROOT as rt
from optparse import OptionParser

# parse command line arguments
parser = OptionParser()
parser.add_option("--xsecfile",dest="xsecfile",default="xsections.txt",help="path to file listing MC sample cross sections and pre selection efficiencies, default=%default")
parser.add_option("--lumifile",dest="lumifile",default="lumis.txt",help="path to file listing integrated luminosities for data samples in pb, default=%default")
parser.add_option("--lumi",dest="lumi",help="integrated luminosity in pb, overrides option --lumifile")
parser.add_option("--idir",dest="idir",default="input",help="path to input directory with root files, default=%default")
parser.add_option("--odir",dest="odir",default="",help="path to output directory with root files, default=IDIR_weighted")
parser.add_option("--counthistpath",dest="counthistpath",default="count",help="total number of processed MC events is take from first bin of histogram COUNTHISTPATH, default=%default")
parser.add_option("--force",dest="force",action="store_true",default=False,help="override ODIR it it already exists, default=%default")
(options, args) = parser.parse_args()

options.idir = options.idir.rstrip("/")
if len(options.odir)== 0:
    options.odir = options.idir + "_weighted"
sys.argv.append("-b")

# check options
def checkLocationOption(location,title,option):
    if not os.path.exists(location):
        print "ERROR: given {0}, {1}, does not exist.".format(title,location)
        print "       Specify {0} with option {1}.".format(title,option)
        sys.exit()

checkLocationOption(options.xsecfile,"cross section file","--xsecfile")
if options.lumi is None:
    checkLocationOption(options.lumifile,"luminosity file","--lumiFile. Alternatively, specifiy luminosity with option --lumi")
checkLocationOption(options.idir,"input directory","--idir")

# read integrated luminosity
lumi =0
if options.lumi is not None:
    print "reading luminosity from command line"
    lumi = float(options.lumi)
else:
    print "reading luminosity from " + options.lumifile
    FILE = open(options.lumifile)
    lines = FILE.read().strip().split("\n")
    FILE.close()
    for line in lines:
        if line.find("#")==0:
            continue
        elements = line.split()
        lumi += float(elements[-1])
print "total lumi: ", lumi

# read xsections and efficiencies
print "reading cross sections from " + options.xsecfile
FILE = open(options.xsecfile)
lines = FILE.read().strip().split("\n")
FILE.close()

effxsec = dict()
for line in lines:
    if line.find("#")==0:
        continue
    elements = line.split()
    sample = elements[0]
    xsection = float(elements[1])
    efficiency = float(elements[2])
    effxsec.update([[sample,xsection*efficiency]])

# list files
print "INPUT DIR: " + options.idir
_file = glob.glob(options.idir + "/*.root")
file = dict()
for __file in _file:
    _sample = os.path.split(__file)[1].replace(".root","")
    file.update([[_sample,__file]])
for _sample in sorted(effxsec.keys()):
    if not _sample in file:
        print "WARNING: no histogram file available for sample", _sample

# create output directory
if os.path.exists(options.odir):
    if options.force:
        shutil.rmtree(options.odir)
    else:
        print "ERROR: ODIR already exists, use option --force to override"
        sys.exit()
print "OUTPUT DIR: " + options.odir
os.mkdir(options.odir)

# apply weights
for _sample in sorted(file.keys()):
    _file = file[_sample]
    weight = None
    # list file content
    tfile = rt.TFile.Open(_file)
    path  = mytools.listRootFile(tfile)
    # check availability of effxsec
    if not _sample in effxsec:
        print "WARNING: no xsection available for sample",_sample,"applying no weight"
    else:
        _effxsec = effxsec[_sample]
        # get total number of events processed
        if not options.counthistpath in path:
            print "ERROR: no such object '" + options.counthistpath + "' in file " + _file
            print "       specify correct option --counthistpath"
            sys.exit()
        countHist = tfile.Get(options.counthistpath)
        if not countHist.InheritsFrom("TH1"):
            print "ERROR: object '" + options.counthistpath + "' in file " + _file + " is no histogram"
            print "       specify correct option --counthistpath"
            sys.exit()
        N = countHist.GetBinContent(1)
        # calculate weight
        if N > 0:
            weight = lumi*_effxsec/N
        else:
            weight = 0
        #print "sample {0}: N={1},effxsec={2} => w={3}".format(_sample,N,_effxsec,weight)
    
    # create output root file
    otfile = rt.TFile.Open(options.odir + "/" + _sample + ".root","RECREATE")

    for _path in path:
        o = tfile.Get(_path)

        # create directories
        if o.InheritsFrom("TDirectory"):
            otfile.mkdir(_path)
            continue

        # clone, weight and write histograms
        if o.InheritsFrom("TH1"):
	  	if o.InheritsFrom("TProfile"):#catch TProfiles
	    		_o = o.Clone()
			print _o
	    		if weight is not None:
				for i in range(0,_o.GetNbinsX()):
					#print "bin "+str(i+1)+": "+str(_o.GetBinEntries(i+1))+" with "+str(_o.GetBinContent(i+1))+ " times "+str(weight)
					_o.SetBinContent(i+1, _o.GetBinContent(i+1)*weight)
		  			_o.SetBinEntries(i+1, _o.GetBinEntries(i+1)*weight)
	    		otfile.cd(os.path.split(_path)[0])
	    		_o.Write()
	  	else:
            		_o = o.Clone()
			if _path.find("/")!=-1 and _path.find("h_")!=-1:
				_pathProfile="p_"+_path
				_pathProfile=_pathProfile.replace("/h_","/p_")
				p = tfile.Get(_pathProfile)
				#print "original path: "+_path
				#print "profile path: "+_pathProfile
				for i in range(0, _o.GetNbinsX()):
					eff = p.GetBinContent(i+1)
					N   = p.GetBinEntries(i+1)
					if N>0:
						_o.SetBinError(i+1, math.sqrt(eff*(1-eff)/N))
						#print "N: "+str(N)+", eff: "+str(eff)+" = "+str(math.sqrt(eff*(1-eff)/N))
					else:
						_o.SetBinError(i+1, 0)	
            		if not _o.GetSumw2():
                		_o.Sumw2()
            		if weight is not None:
                		_o.Scale(weight)
            		otfile.cd(os.path.split(_path)[0])
            		_o.Write()
	    
	if o.InheritsFrom("TGraphAsymmErrors"):
            _o = o.Clone()
            if weight is not None:
                for x in range(0, _o.GetXaxis().GetNbins()):
		    xcoord=rt.Double(0.)
		    ycoord=rt.Double(0.)
		    _o.GetPoint(x, xcoord, ycoord)
		    _o.SetPoint(x, xcoord, ycoord*weight)
		    _o.SetPointEYlow(x, _o.GetErrorYlow(x)*weight)
		    _o.SetPointEYhigh(x, _o.GetErrorYhigh(x)*weight)
            otfile.cd(os.path.split(_path)[0])
            _o.Write()

    # closing ...
    tfile.Close()
    otfile.Close()
            
        
        
