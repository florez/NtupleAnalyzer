The BSM3G_TNT_Analyzer is supposed to be a platform providing generic analysis scripts on BSM3G_TNTuples.
In order to utilize this platform, a BSM3G_TNTuple .root-file is needed.

The python script "readNtupleBranches.py" generates a script capable of reading in all branches in the [input].root file.
To do so, run in a shell:
python readNtupleBranches.py --inFile [input].root --AnaName [YourAnalyzerName]
Please replace [input] and [YourAnalyzerName] with valid choices.

Afterwards, compile the new barebone analyzer script with make [YourAnalyzerName] (do not add a ".cc")

The barebone script provides a few general structures that can be used for analysis in the "interface" folder and a few utilities.
The interface contains a barebone class "MyEventCollection" in "foundations.h" that is supposed to contain all event properties of interest. 
General event properties are directly initialized in that class, different objects are defined in "objects.h". 
Here, only the example of minimal properties of a muon is given. This can easily be extended as needed. 
In the analyzer script, each event's properties after object selection are stored in the MyEventCollection container. 
They are then subjected to a configurable "Selection" as defined in "selection.h".
This Selection performs the configurable cuts and stores cut efficiencies in a "MyProfileCollection" container and histograms after all cuts in a "MyHistoCollection" container. 
Both of these containers are specified in "foundation.h" and each contains the same name histograms with a "p_" for profile and an "h_" for histograms. 
It is important to maintain this nomenclature, in order to correctly use the utilities provided.
The way the profiles and histograms are filled in case of a pass or fail argument from the selector is provided in "CommonHistoFillProcedure.h" with a few safeguards.

After running the barebone script, an output file is generated and can then be further manipulated with a few utilities. 
If you have several samples with different weights and a certain amount of luminosity in terms of data to weigh it to, the "weight.py" script copies a folder of root files 
and automatically adds weights specified in "xsection.txt". Be sure to have the name minus ".root" of each dataset-file in this file. 
Then, the first entry after a space is the cross section, the second entry after a space is the ntuplization efficiency (if everything runs well, this should always be 1). 
The luminosity all datasets are weighted to are listed in "lumis.txt" and simply summed. Datasets not specified are not weighted and just copied. 
Moreover, this script changes all uncertainties per histogram bin from the poissonian error to the efficiency error before applying the reweighting.
