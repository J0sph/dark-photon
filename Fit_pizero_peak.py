#! /usr/bin/python                                                            

################# Fit_pizero_peak.py ################################         
#                                                                   #         
#    Author: Felicia Volle                                          #         
#    Date:   06/02/2025                                             #         
#    Description: fit eegamma invariant mass distribution           #         
#                                                                   #         
#    Run via: lb-conda default python Fit_pizero_peak.py            #         
#                                                                   #        
#####################################################################         
## If not on lxplus, need to run :
## - kinit -f fevolle@CERN.CH                                
## - source /cvmfs/lhcb.cern.ch/lib/LbEnv

################# IMPORT ############################################       
import os, sys
import argparse
from ROOT import RDataFrame, RDF, TLegend, TCanvas, TFile, TChain, RooRealVar, RooArgList, RooArgSet, RooGaussian, RooCBShape, RooLinkedList, RooArgusBG, RooAddPdf, RooFit, RooDataSet, RooCrystalBall, RooGenericPdf, gROOT, EnableImplicitMT, kRed, kGreen, kViolet, kBlue, kDashed, TPaveText, gStyle
import pickle
import yaml

gROOT.SetBatch(True)
EnableImplicitMT()
#PLOTTING
gROOT.ProcessLine(".x lhcbStyle.C")
print("--> LHCbStyle loaded")


################# SETTING #########################################


dataBlock_dict = {4: ["data_24c2a_magdown_qee", "magdown", "2024", "FillNumber >= 9808 && FillNumber <= 9910" ],
                  3 : ["data_24c2a_magup_qee", "magup", "2024", "FillNumber >= 9911 && FillNumber <= 9943" ],
                  2 : ["data_24c2a_magup_qee", "magup", "2024", "FillNumber >= 9945 && FillNumber <= 9978" ],
                  1 : ["data_24c2a_magup_qee", "magup", "2024", "FillNumber >= 9982 && FillNumber <= 10056" ],
                  5 : ["data_24c3a_magup_qee", "magup", "2024", "FillNumber >= 10059 && FillNumber <= 10102"],
                  6: ["data_24c3a_magdown_qee", "magdown", "2024", "FillNumber >= 10104 && FillNumber <= 10190" ],
                  7: ["data_24c4a_magdown_qee", "magdown", "2024", "FillNumber >= 10197 && FillNumber <= 10213" ],
                  8 : ["data_24c4a_magup_qee", "magup", "2024", "FillNumber >= 10214 && FillNumber <= 10232" ],
                  9: ["data_25c1_magdown_qee", "magdown", "2025", "FillNumber >= 10489 && FillNumber <= 10732" ]
}       



track_type   = "Prompt" # Prompt or Displaced
samesign = False
fd_tag = track_type
if samesign == True:
    fd_tag += "SS"
    
mcfile = "Dalitz" # 3-body decay of pi0, as in SM. 
selstring = "hlt1_BremOverlap_Ecal_PID_Kin_Rho_Topo"
filenbr = -1 

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
    
    return datasets(polarity="magup", eventtype="39122948",
                    version="v1r2683",
                    datatype="2024", filetype="funtuple_qee.root",
                    name="mcblock12_magup_pi0_dalitz_turbo_ftuple")

def setup_files(data_type, track_type, sel="", samesign=False, filenbr=1,
data_name = None, data_polarity = None, data_year = None, output_dir = None):


    
    ofilename = f"{data_type}_{track_type}"+"_SS"*samesign+".root"

    if output_dir and data_type == "Data" :
        ofilename = os.path.join(output_dir, ofilename)

    if os.path.exists(ofilename):
        return ofilename


    print(f"Merging files into {ofilename}")
    ifilenames = get_pfns(data_type=data_type, data_name= data_name, data_polarity=data_polarity, data_year=data_year)
    if data_type == "Data" and not samesign:
        ifilenames = ifilenames[:2]
    if data_year == "2025":
        treename = "Rho_Tuple_NoIP"
    else:
        treename = f"Rho_Tuple_{track_type}"
    if samesign: treename += "SS"
    treename += "/DecayTree"
    
    chain = TChain(treename)
    i = 0
    for ifilename in ifilenames:
        if filenbr==-1:
            chain.Add(ifilename)
        elif i<filenbr:
            chain.Add(ifilename)
        else:
            break
        i += 1
        
    ofile = TFile(ofilename, "recreate")
    otree = chain.CopyTree(sel)
    otree.SetName("DecayTree")
    print(f"Merger done, output {data_type} TTree has {otree.GetEntries()} entries")
    otree.Write()
    ofile.Close()

    return ofilename


def fit_data(mcfilename, datafilename, track_type, samesign=False, output_dir=None):
    mctreename = "DecayTree"

    mcfile = TFile(mcfilename)
    mctree = mcfile.Get(mctreename)

    massvar = RooRealVar("rho_M", "rho_M", 50, 400)#50, 250)
    listvars = RooArgSet(massvar)

    mcdata = RooDataSet("mcdata", "mcdata", listvars, RooFit.Import(mctree))
    print("MC data has", mcdata.sumEntries(), "entries")
    
    ### Create the fit model
    a1 = RooRealVar("a1", "a1", 1.18, 0.1, 100)
    a2 = RooRealVar("a2", "a2", 1.46, 0.1, 100)
    n1 = RooRealVar("n1", "n1", 1.2, 0.1, 25)
    n2 = RooRealVar("n2", "n2", 1.22, 0.1, 100)

    mean = RooRealVar("mean", "mean", 135, 10, 155)

    sigma = RooRealVar("sigma", "sigma", 14, 1, 100)

    sigmodel = RooCrystalBall("sigmodel", "sigmodel", massvar, mean, sigma, a1, n1, a2, n2)
    #sigmodel = RooCrystalBall("sigmodel", "sigmodel", massvar, mean, sigma, a1, n1, doubleSided=False)


    # Fit options we can think about later
    #mean_constr = RooGaussian("mean_constr", "mean_constr", mean,
    #                             RooFit.RooConst(135), RooFit.RooConst(5))

    ### Create the fit of the model to the mc
    print("\n> Fit MC sample")
    mcresult = sigmodel.fitTo(mcdata, RooFit.Save(True), RooFit.Verbose(False))
    #, RooFit.Minos(1),
    #RooFit.ExternalConstraints(rt.RooArgSet(mean_constr)))

    mcresult.Print()

    ### Save result from above in a dictionnary
    d = {}
    d["a1"] = a1.getVal()
    d["a1_unc"] = a1.getError()
    d["a2"] = a2.getVal()
    d["a2_unc"] = a2.getError()
    d["n1"] = n1.getVal()
    d["n1_unc"] = n1.getError()
    d["n2"] = n2.getVal()
    d["n2_unc"] = n2.getError()
    d["MCfit_status"] = mcresult.status()

    
    ### Plotting
    #gROOT.ProcessLine(".x lhcbStyle.C")
    #print("--> LHCbStyle loaded")
    

    c1 = TCanvas("c1", "c1", 800, 600)
    frame = massvar.frame()
    mcdata.plotOn(frame, RooFit.Binning(100))
    sigmodel.plotOn(frame)

    maxi=frame.GetMaximum()
    frame.SetMaximum(maxi*1.2) 

    frame.SetXTitle("m(e^{+}e^{-}#gamma) [MeV/c^{2}]")
    bin_width = (massvar.getMax() - massvar.getMin()) / 100
    frame.SetYTitle(f"Candidates / ({bin_width:.1f} MeV/c^{{2}})")
    
    frame.Draw()

    lhcbName = TPaveText(gStyle.GetPadLeftMargin() + 0.03,
                         0.87 - gStyle.GetPadTopMargin(),
                         gStyle.GetPadLeftMargin() + 0.38,
                         0.95 - gStyle.GetPadTopMargin(),
                         "BRNDC")
    
    lhcbName.AddText("LHCb Simulation")
    
    lhcbName.SetFillColor(0)
    lhcbName.SetFillStyle(0)
    lhcbName.SetTextAlign(12)
    lhcbName.SetTextSize(0.06)
    lhcbName.SetBorderSize(0)
    lhcbName.Draw()

    c1.Update()
    plotofmcfit = f"MCfit_{fd_tag}.pdf"
    c1.SaveAs(plotofmcfit)
    c1.Close()
    print("Saved figure : ", plotofmcfit)

    mcfile.Close()
    
    ### Get the data
    datatreename = "DecayTree"

    print("Open ", datafilename)
    datafile = TFile(datafilename)
    datatree = datafile.Get(datatreename)

    data = RooDataSet("data", "data", listvars, RooFit.Import(datatree))
    print("Data sample has ", data.sumEntries(), "entries")
    
    ### Create total model of signal + background
    # Set shape parameter of the signal shape to the MC parameter from the previsous fit
    for var in [a1, n1, a2, n2]:
        var.setConstant(True)
        var.Print()
        
    #m0 = RooRealVar("m0", "m0", 1000, 581, 10000)
    #c  = RooRealVar("c", "c", -20, -100, -1)

    #bkgmodel = RooArgusBG("bkgmodel", "bkgmodel", massvar, m0, c)
    #q2_expr = "TMath:Sqrt( (rho_M*rho_M - (mA + mB)*(mA + mB)) * (rho_M*rho_M - (mA - mB)*(mA - mB)) ) / (2*rho_M)"
    q2_expr = "(rho_M-50)**A * exp(-B*(rho_M-50))"
    mA = RooRealVar("A", "A", 0.5, 1e-3, 75)
    mB = RooRealVar("B", "B", 1e-6, 1e-2)
    bkgmodel = RooGenericPdf("bkgmodel", "bkgmodel",
                             q2_expr, RooArgList(massvar, mA, mB))

    nsig = RooRealVar("nsig", "nsig", 1, data.numEntries())
    nbkg = RooRealVar("nbkg", "nbkg", 1, data.numEntries())

    model = RooAddPdf("model", "model", RooArgSet(sigmodel, bkgmodel),
                      RooArgSet(nsig, nbkg))

    ### Fit the model to the data
    print("\n> Fit data sample")
    result = model.fitTo(data, RooFit.Save(True), RooFit.Verbose(False),
                         RooFit.Extended(True))
    #RooFit.ExternalConstraints(rt.RooArgSet(mean_constr)))

    result.Print()
    
    ### Plot the result
    c2 = TCanvas("c2", "c2", 800, 600)
    frame2 = massvar.frame(RooFit.Name("Full fit"), RooFit.Title(""),
                           RooFit.Bins(100))

    data.plotOn(frame2, RooFit.Binning(100))

    frame2.Clone("frame 2")
    #model.getComponents().Print("v")

    model.plotOn(frame2,
                 RooFit.Name("model"),
                 RooFit.LineColor(kViolet+4))


    
    model.plotOn(frame2, 
                 RooFit.Components("sigmodel"),
                 RooFit.Name("sigmodel"),
                 RooFit.LineColor(kRed+2),
                 RooFit.LineStyle(kDashed))

    model.plotOn(frame2,
                 RooFit.Components("bkgmodel"),
                 RooFit.Name("bkgmodel"),
                 RooFit.LineColor(kBlue+2),
                 RooFit.LineStyle(kDashed))

    frame2.Draw()

    ### Add a legend
    leg = TLegend(0.71, 0.35, 0.91, 0.53)
    leg.SetTextSize(0.055)
    leg.AddEntry(frame2.findObject("model"), "Total", "l")
    leg.AddEntry(frame2.findObject("sigmodel"), "#pi^{0}", "l")
    leg.AddEntry(frame2.findObject("bkgmodel"), "Background", "l")
    leg.Draw()

    maxi2=frame2.GetMaximum()
    frame2.SetMaximum(maxi2*1.2) 

    frame2.SetXTitle("m(e^{+}e^{-}#gamma) [MeV/c^{2}]")
    frame2.SetYTitle(f"Candidates / ({bin_width:.1f} MeV/c^{{2}})")

    lhcbName2 = TPaveText(gStyle.GetPadLeftMargin() + 0.03,
                         0.87 - gStyle.GetPadTopMargin(),
                         gStyle.GetPadLeftMargin() + 0.38,
                         0.95 - gStyle.GetPadTopMargin(),
                         "BRNDC")
    
    lhcbName2.AddText("LHCb Preliminary")
    
    lhcbName2.SetFillColor(0)
    lhcbName2.SetFillStyle(0)
    lhcbName2.SetTextAlign(12)
    lhcbName2.SetTextSize(0.06)
    lhcbName2.SetBorderSize(0)
    lhcbName2.Draw()
    
    c2.Update()
    figurename = f"DataFit_{fd_tag}.pdf"
    if output_dir: 
        figurename = os.path.join(output_dir, figurename)
    c2.SaveAs(figurename)
    print("Save figure : ", figurename)
    c2.Close()
    
    ### Print the result
    print("\n --- Fit summary ---")
    print(" Nsig = ", nsig.getVal(), "+-", nsig.getError())
    print(" Nbkg = ", nbkg.getVal(), "+-", nbkg.getError())
    print(" mean = ", mean.getVal(), "+-", mean.getError())
    print(" sigma = ", sigma.getVal(), "+-", sigma.getError())

    ### Save result from above in a dictionnary
    d["Nsig"] = nsig.getVal()
    d["Nsig_unc"] = nsig.getError()
    # ToDo: add other fit variable values
    d["mean"] =mean.getVal()
    d["mean_unc"] = mean.getError()
    d["sigma"] =sigma.getVal()
    d["sigma_unc"] = sigma.getError()
    d["mA"] = mA.getVal()
    d["mA_unc"] = mA.getError()
    d["mB"] = mB.getVal()
    d["mB_unc"] = mB.getError()
    d["Nbkg"] = nbkg.getVal()
    d["Nbkg_unc"] = nbkg.getError()
    d["DataFit_status"] = result.status()

    
    # Save result dictionnary into a yaml file
    yfile = "FitResults.yml"
    if output_dir: 
        yfile = os.path.join(output_dir, yfile)
    with open(yfile, 'w') as outfile:
        yaml.dump(d, outfile, default_flow_style=False)
        print("\nSaved yml file in: ", yfile)
    
    return True


################# MAIN ###############################################       
    

if __name__ == "__main__":
    print(">> START")

    print("\n>> Get selection for MC file")
    mcsel = get_selection("truM_"+selstring, fd_tag, mcfile)

    print("Combine files in one MC sample")
    mcfilename = setup_files(data_type="MC", track_type=track_type,
                             sel=mcsel, samesign=samesign, filenbr=1)

    print(mcfilename)


    for dataBlock in dataBlock_dict.keys():
        print("\n >> Start with data block", dataBlock)
        data_string, data_polarity, data_year, fillnbr_selection  = dataBlock_dict[dataBlock]
        print(f"name: {data_string}, polarity: {data_polarity}, year: {data_year}, fill number: {fillnbr_selection}")

        block_dir = f"block_{dataBlock}"

        expected_files = [f"{block_dir}/DataFit_{fd_tag}.pdf",
                        f"{block_dir}/FitResults.yml"]

        if os.path.exists(block_dir) and all(os.path.exists(f) for f in expected_files):
            print(f">> Skipping block {dataBlock}: already processed.")
            continue
        else:
            if not os.path.exists(block_dir):
                os.mkdir(block_dir)

        print("\n>> Get selection for data file")
        datasel = get_selection(selstring, fd_tag, "", fillnbr= fillnbr_selection, data_year=data_year)

        print("Combine files in one data sample")
        datafilename = setup_files(data_type="Data", track_type=track_type,
                                sel=datasel, samesign=samesign, filenbr=filenbr,
                                data_name = data_string, data_polarity = data_polarity, data_year = data_year, output_dir=block_dir)
        print(datafilename)
        
        print("\n>> Start fitting")
        fit_ok = fit_data(mcfilename=mcfilename, datafilename=datafilename,
                    track_type=track_type, samesign=samesign,output_dir=block_dir)

        print(f">> BLOCK {dataBlock} DONE")


    print(">> DONE")

## EOF
