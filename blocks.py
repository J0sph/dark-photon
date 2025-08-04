#! /usr/bin/python                                                            

################# blocks.py #########################################         
#                                                                   #                  
#    Run via: lb-conda default python blocks.py                     #         
#                                                                   #        
#####################################################################         
## If not on lxplus, need to run :
## - kinit -f fevolle@CERN.CH                                
## - source /cvmfs/lhcb.cern.ch/lib/LbEnv

################# IMPORT ############################################       
import os, sys
import argparse
from ROOT import  RDataFrame, RDF, TFile, TChain, RooRealVar, RooArgList, RooArgSet, RooGaussian, RooCBShape, RooLinkedList, RooArgusBG, RooAddPdf, RooFit, RooDataSet, RooGenericPdf, gROOT, EnableImplicitMT
import yaml
from datetime import datetime
import uproot
import numpy as np
import argparse

gROOT.SetBatch(True)
EnableImplicitMT()


root_file_range = ""

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--range", help="File block range (a, b, c, ...)", required=True)
args = parser.parse_args()

import string
root_file_range = args.range.lower()  # 'a', 'b', ...
range_index = string.ascii_lowercase.index(root_file_range)
block_size = 100
start_file = range_index * block_size
end_file = start_file + block_size
#start_file = 0
#end_file = 1000

################# SETTING #########################################


dataBlock_dict = {#4: ["data_24c2a_magdown_qee", "magdown", "2024", "FillNumber >= 9808 && FillNumber <= 9910" ],
                  #3 : ["data_24c2a_magup_qee", "magup", "2024", "FillNumber >= 9911 && FillNumber <= 9943" ],
                  #2 : ["data_24c2a_magup_qee", "magup", "2024", "FillNumber >= 9945 && FillNumber <= 9978" ],
                  1 : ["data_24c2a_magup_qee", "magup", "2024", "FillNumber >= 9982 && FillNumber <= 10056" ],
                  #5 : ["data_24c3a_magup_qee", "magup", "2024", "FillNumber >= 10059 && FillNumber <= 10102"],
                  #6: ["data_24c3a_magdown_qee", "magdown", "2024", "FillNumber >= 10104 && FillNumber <= 10190" ],
                  #7: ["data_24c4a_magdown_qee", "magdown", "2024", "FillNumber >= 10197 && FillNumber <= 10213" ],
                  #8 : ["data_24c4a_magup_qee", "magup", "2024", "FillNumber >= 10214 && FillNumber <= 10232" ],
                  #9: ["data_25c1_magdown_qee", "magdown", "2025", "FillNumber >= 10489 && FillNumber <= 10732" ]
}       



track_type   = "Prompt" # Prompt or Displaced
samesign = False
fd_tag = track_type
if samesign == True:
    fd_tag += "SS"
    
mcfile = "Dalitz" # 3-body decay of pi0, as in SM. 
selstring = "hlt1_BremOverlap_Ecal_PID_Kin_Rho_Topo"


################# SELECTION #########################################       


def get_selection(sel, fd_tag, string, fillnbr = None, data_year = None):
    print("Get selection for", fd_tag, " : ", sel)
    selBremOverlap = "!(em_HASBREMADDED==1 && gamma_CaloNeutralID==em_BREMHYPOID ) && !(ep_HASBREMADDED==1 && gamma_CaloNeutralID==ep_BREMHYPOID)"
    selEcal = "!(gamma_CaloNeutralCol > 22 && gamma_CaloNeutralCol < 41 && gamma_CaloNeutralRow > 24 && gamma_CaloNeutralRow < 39)"
    selPID = "ep_PID_E>2 && em_PID_E>2 && em_PROBNN_GHOST<0.4 && ep_PROBNN_GHOST<0.4 &&em_PROBNN_PI<0.2 && ep_PROBNN_PI<0.2"
    selKin = "KS0_PT>1000 && KS0_MAX_DOCA<0.1 && KS0_OpeningAngle>0.0001" #&& KS0_OpeningAngle<0.02"                       
    selRho = "rho_BPV_IP < 0.1 && rho_OpeningAngle>0.001" # && rho_OpeningAngle < 0.015"                                    
    selTopo = "em_QOVERP>-0.3 && ep_QOVERP<0.3 && ep_CHI2_DOF<2 && em_CHI2_DOF<2"

    selList = []
    if fillnbr != None:
        selList.append(fillnbr)
    
    if "truM" in sel:
        if "minBias" not in string:
            # pi0 peak in minBias visible for 10, 50 and 60                     
            truM = "(rho_BKGCAT <= 10 || rho_BKGCAT == 50 || rho_BKGCAT == 60)"
        else:
            # Do not exclude 60 here since no peak in minBias                   
            truM ="(rho_BKGCAT != 10 && rho_BKGCAT != 50)"
        selList.append(truM)
    if "hlt1" in sel:

        dec = "displacedDecision"
        if fd_tag == "Prompt" or fd_tag == "PromptSS":
            dec = "promptDecision"
            if data_year == "2025":
                dec = "NoIPDecision"
        ss = ""
        if "SS" in fd_tag:
            ss = "_SS"

        if data_year != "2025":
            lines = []
            for i in range(1, 5):
                lines.append(f"KS0_Hlt1DiElectronLowMass{ss}_massSlice{i}_{dec}_TOS")
        else:
            lines = [f"KS0_Hlt1DiElectronLowMass{ss}_{dec}_TOS"]

        hlt1sel = ' || '.join(lines)
        selList.append(f"( {hlt1sel} )")
        
    if "BremOverlap" in sel:
        selList.append(selBremOverlap)
    if "Ecal" in sel:
        selList.append(selEcal)
    if "PID" in sel:
        selList.append(selPID)
    if "Kin" in sel:
        selList.append(selKin)
    if "Rho" in sel:
        selList.append(selRho)
    if "Topo" in sel:
        selList.append(selTopo)

    selString = ' && '.join(selList)
    print(f"selString: {selString}")
    return selString

    
################# FUNCTIONS #########################################       

def get_pfns(data_type, data_name=None, data_polarity=None, data_year=None):
    from apd import AnalysisData

    datasets = AnalysisData("qee", "pi02eegamma_r3")

    if data_type == "Data":
        return datasets(polarity=data_polarity, eventtype="94000000",
                        datatype=data_year, filetype="qee_funtuple.root",
                        name=data_name)
    
def setup_files(data_type, track_type, sel="", samesign=False,
data_name = None, data_polarity = None, data_year = None, output_dir = None,start_index=0, end_index=10):

    #def get_min_max(tree, branch_name):
    #    values = []
    #    for entry in tree:
    #        val = getattr(entry, branch_name, None)
    #        if val is not None:
    #            values.append(val)
    #    if values:
    #        return min(values), max(values)
    #    else:
    #        return None, None
    
    ofilename = f"{data_type}_{track_type}"+"_SS"*samesign+f"_{root_file_range}"+".root"

    if output_dir:
        ofilename = os.path.join(output_dir, ofilename)

    if os.path.exists(ofilename):
        return ofilename


    print(f"Merging files into {ofilename}")
    ifilenames = get_pfns(data_type=data_type, data_name= data_name, data_polarity=data_polarity, data_year=data_year)
    size = len(ifilenames)
    print(f"Number of files: {size}")

    print(f"The last file was: {ifilenames[start_index-1]}")

    if end_index < size:
        print(f"\nThe next file to be processed in the next block is:\n{ifilenames[end_index]}")
    else:
        print("\nThere are no more files after this block.")

    print(f"From file: {start_index}: {ifilenames[start_index]}")

    ifilenames = ifilenames[start_index:end_index]

    if end_index < size:
        print(f"To file {end_index-1}: {ifilenames[-1]}")
    else:
        print(f"To file {size}: {ifilenames[-1]}")

    print(f"Number of files: {len(ifilenames)} for sub_block_{root_file_range}")
    
    if len(ifilenames) == 0:
        print(f"No files in range {start_index}:{end_index}")
        return None

    if data_year == "2025":
        treename = "Rho_Tuple_NoIP"
    else:
        treename = f"Rho_Tuple_{track_type}"
    if samesign: treename += "SS"
    treename += "/DecayTree"
    
    chain = TChain(treename)
    i = 0

    for ifilename in ifilenames:
        chain.Add(ifilename)
        print(ifilename)
    
        i += 1

    dirpath = os.path.dirname(ofilename)
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
        
    ofile = TFile(ofilename, "recreate")

    #if chain.GetEntries() > 0:
    #    chain.SetBranchStatus("*", 0)
    #    chain.SetBranchStatus("FillNumber", 1)
    #    min_before, max_before = get_min_max(chain, "FillNumber")
    #    print(f"Before selection: fillnumber[{min_before}, {max_before}]")
    #    chain.SetBranchStatus("*", 1)

    otree = chain.CopyTree(sel)

    #if otree.GetEntries() > 0:
    #    otree.SetBranchStatus("*", 0)
    #    otree.SetBranchStatus("FillNumber", 1)
    #    min_after, max_after = get_min_max(otree, "FillNumber")
    #    print(f"After selection: fillnumber[{min_after}, {max_after}]")
    #    otree.SetBranchStatus("*", 1)
    
    otree.SetName("DecayTree")
    print(f"Merger done, output {data_type} TTree has {otree.GetEntries()} entries")
    otree.Write()
    ofile.Close()

    return ofilename

################# MAIN ###############################################       
    

if __name__ == "__main__":

    print(">> START")
    import time

    start_time = time.time()


    for dataBlock in dataBlock_dict.keys():
        print("\n >> Start with data block", dataBlock)
        data_string, data_polarity, data_year, fillnbr_selection  = dataBlock_dict[dataBlock]
        print(f"name: {data_string}, polarity: {data_polarity}, year: {data_year}, fill number: {fillnbr_selection}")

       
        block_dir_data = f"/eos/lhcb/user/c/castillj/DarkPhoton-samples/block_{dataBlock}"


        expected_file = os.path.join(block_dir_data, f"Data_{fd_tag}"+"_SS"*samesign+f"_{root_file_range}"+".root")  
        #print(expected_file)

        if not os.path.exists(expected_file):
            
            print("\n>> Get selection for data file")
            datasel = get_selection(selstring, fd_tag, "", fillnbr= fillnbr_selection, data_year=data_year)

            print("Combine files in one data sample")
            datafilename = setup_files(
                data_type="Data",
                track_type=track_type,
                sel=datasel,
                samesign=samesign,
                data_name=data_string,
                data_polarity=data_polarity,
                data_year=data_year,
                output_dir=block_dir_data,
                start_index=start_file,
                end_index=end_file
            )
            print(datafilename)
            

            print(f">> BLOCK {dataBlock} DONE")

        else: 
            print(">> BLOCK ALREADY EXISTS")


    print(">> DONE")

    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"Time: {elapsed_time:.2f} s")

## EOF
