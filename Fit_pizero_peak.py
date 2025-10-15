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
from ROOT import TLatex, TLine, RooBernstein, RooFormulaVar, RooBukinPdf, RooNovosibirsk, RooChebychev, RDataFrame, RDF, TLegend, TCanvas, TFile, TChain, RooRealVar, RooArgList, RooArgSet, RooGaussian, RooCBShape, RooLinkedList, RooArgusBG, RooAddPdf, RooFit, RooDataSet, RooCrystalBall, RooGenericPdf, gROOT, EnableImplicitMT,kOrange, kRed, kGreen, kViolet, kBlue, kDashed, TPaveText, gStyle
import pickle
import yaml
from datetime import datetime
import time
import uproot
import numpy as np
import argparse

gROOT.SetBatch(True)
EnableImplicitMT()
#PLOTTING
gROOT.ProcessLine(".x lhcbStyle.C")
print("--> LHCbStyle loaded")

#root_file_range = ""

################# SETTING #########################################


dataBlock_dict = {4: ["data_24c2a_magdown_qee", "magdown", "2024", "FillNumber >= 9808 && FillNumber <= 9910", 0.7849178894141272, 0.04709507336484763 ],
                  3 : ["data_24c2a_magup_qee", "magup", "2024", "FillNumber >= 9911 && FillNumber <= 9943", 0.6526399361492664, 0.03915839616895598],
                  2 : ["data_24c2a_magup_qee", "magup", "2024", "FillNumber >= 9945 && FillNumber <= 9978", 0.5797862448573195, 0.03478717469143917],
                  1 : ["data_24c2a_magup_qee", "magup", "2024", "FillNumber >= 9982 && FillNumber <= 10056", 1.113538792778741, 0.06681232756672445],
                  5 : ["data_24c3a_magup_qee", "magup", "2024", "FillNumber >= 10059 && FillNumber <= 10102", 1.0930350319410267, 0.0655821019164616],
                  6: ["data_24c3a_magdown_qee", "magdown", "2024", "FillNumber >= 10104 && FillNumber <= 10190", 0.9296609005470799, 0.05577965403282479],
                  7: ["data_24c4a_magdown_qee", "magdown", "2024", "FillNumber >= 10197 && FillNumber <= 10213", 0.6962158423024445, 0.04177295053814667],
                  8 : ["data_24c4a_magup_qee", "magup", "2024", "FillNumber >= 10214 && FillNumber <= 10232", 0.41624918761675584, 0.02497495125700535],
                  9: ["data_25c1_magdown_qee", "magdown", "2025", "FillNumber >= 10489 && FillNumber <= 10732", 0.41534220, 0 ]
}       



track_type   = "Prompt" # Prompt or Displaced
samesign = False
fd_tag = track_type
if samesign == True:
    fd_tag += "SS"
    
mcfile = "Dalitz" # 3-body decay of pi0, as in SM. 
selstring = "hlt1_BremOverlap_Ecal_PID_Kin_Rho_Topo"
filenumber = -1

################# SELECTION #########################################       


def get_selection(sel, fd_tag, string, fillnbr = None, data_year = None, brew_category=None):
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

    if brew_category:
        selList.append(f"(em_HASBREMADDED + ep_HASBREMADDED) == {brew_category}")

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
    
    if data_type == "MC_pi0":
        return datasets(polarity="magup", eventtype="39122948",
                        version="v1r2683",
                        datatype="2024", filetype="funtuple_qee.root",
                        name="mcblock12_magup_pi0_dalitz_turbo_ftuple")

    if data_type == "MC_eta":
        return datasets(polarity="magup", eventtype="39122432",
                        version="v1r2683",
                        datatype="2024", filetype="funtuple_qee.root",
                        name="mcblock12_magup_eta_dalitz_turbo_ftuple")

def setup_files(data_type, track_type, sel="", samesign=False, filenbr=1,
data_name = None, data_polarity = None, data_year = None, output_dir = None):
    
    ofilename = f"{data_type}_{track_type}"+"_SS"*samesign+".root"

    if output_dir:
        ofilename = os.path.join(output_dir, ofilename)

    if os.path.exists(ofilename):
        return ofilename


    print(f"Merging files into {ofilename}")
    ifilenames = get_pfns(data_type=data_type, data_name= data_name, data_polarity=data_polarity, data_year=data_year)
    print(f"Number of files: {len(ifilenames)}")
    #if data_type == "Data" and not samesign and root_file_range == "a":
    #    ifilenames = ifilenames[:2]

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
            print(ifilename)
        elif i<filenbr:
            chain.Add(ifilename)
            print(ifilename)
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


def fit_data(mcfilenames, datafilename, track_type, samesign=False, output_dir=None, luminosity=1, luminosity_uncertainty=0):
    
    mctreename = "DecayTree"

    mcfile_pi0 = TFile(mcfilenames[0])
    mctree_pi0 = mcfile_pi0.Get(mctreename)

    mcfile_eta = TFile(mcfilenames[1])
    mctree_eta = mcfile_eta.Get(mctreename)

    massvar = RooRealVar("rho_M", "rho_M", 50, 600) 


    listvars = RooArgSet(massvar)
    

    mcdata_pi0 = RooDataSet("mcdata_pi0", "mcdata_pi0", listvars, RooFit.Import(mctree_pi0))
    print("MC data has", mcdata_pi0.sumEntries(), "entries")

    mcdata_eta = RooDataSet("mcdata_eta", "mcdata_eta", listvars, RooFit.Import(mctree_eta))
    print("MC data has", mcdata_eta.sumEntries(), "entries") 

    a1_pi0 = RooRealVar("a1_pi0", "a1_pi0", 1.5, 0.5, 5.0)
    n1_pi0 = RooRealVar("n1_pi0", "n1_pi0", 29, 0, 30)
    a2_pi0 = RooRealVar("a2_pi0", "a2_pi0", 1.5, 0.5, 5.0)
    n2_pi0 = RooRealVar("n2_pi0", "n2_pi0", 2.0, 0.01, 30.0)
    sigma_pi0 = RooRealVar("sigma_pi0", "sigma_pi0", 18.22, 2, 20)
    mean_pi0  = RooRealVar("mean_pi0", "mean_pi0", 130, 120, 140) 

    sigmodel_pi0 = RooCrystalBall("sigmodel_pi0", "sigmodel_pi0", massvar, mean_pi0, sigma_pi0, a1_pi0, n1_pi0, a2_pi0, n2_pi0)

    mean_eta  = RooRealVar("mean_eta", "mean_eta", 550, 530, 600)  
    sigma_eta = RooRealVar("sigma_eta", "sigma_eta", 12, 5, 50)    
    a1_eta    = RooRealVar("a1_eta", "a1_eta", 1.5, 0.1, 5)        
    n1_eta    = RooRealVar("n1_eta", "n1_eta", 2.0, 0.5, 100)       
    a2_eta    = RooRealVar("a2_eta", "a2_eta", 2.0, 0.5, 5)        
    n2_eta    = RooRealVar("n2_eta", "n2_eta", 5.0, 0.5, 20)   
   
    sigmodel_eta = RooCrystalBall("sigmodel_eta", "sigmodel_eta", massvar, mean_eta, sigma_eta, a1_eta, n1_eta, a2_eta, n2_eta )
 
    # Fit options we can think about later
    #mean_constr = RooGaussian("mean_constr", "mean_constr", mean,
    #                             RooFit.RooConst(135), RooFit.RooConst(5))


    #mean = RooRealVar("mean", "mean", 135, 120, 150)    
    #sigma = RooRealVar("sigma", "sigma", 10, 0.1, 100)             
    #tau = RooRealVar("tau", "tau", 1, -10.0, 10.0) 

    #sigmodel = RooNovosibirsk("novosibirsk", "Novosibirsk Function", massvar, mean, sigma, tau)  


    ### Create the fit of the model to the mc
    print("\n> Fit MC sample")
    mcresult_pi0 = sigmodel_pi0.fitTo(mcdata_pi0, RooFit.Save(True), RooFit.Verbose(False))
    #, RooFit.Minos(1),
    #RooFit.ExternalConstraints(rt.RooArgSet(mean_constr)))
    mcresult_pi0.Print()

    mcresult_eta = sigmodel_eta.fitTo(mcdata_eta, RooFit.Save(True), RooFit.Verbose(False))
    #, RooFit.Minos(1),
    #RooFit.ExternalConstraints(rt.RooArgSet(mean_constr)))
    mcresult_eta.Print()

    #chi2_MC = mcresult.chi2()
    #print(f"Chi2 from MC fit: {chi2_MC:.2f}")
    

    ### Save result from above in a dictionnary
    d = {}
    d["a1_MC_pi0"] = a1_pi0.getVal()
    d["a1_MC_pi0_unc"] = a1_pi0.getError()
    d["a2_MC_pi0"] = a2_pi0.getVal()
    d["a2_MC_pi0_unc"] = a2_pi0.getError()
    d["n1_MC_pi0"] = n1_pi0.getVal()
    d["n1_MC_pi0_unc"] = n1_pi0.getError()
    d["n2_MC_pi0"] = n2_pi0.getVal()
    d["n2_MC_pi0_unc"] = n2_pi0.getError()
    d["MCfit_pi0_status"] = mcresult_pi0.status()
    d["a1_MC_eta"]=a1_eta.getVal()
    d["a1_MC_eta_unc"] = a1_eta.getError()
    d["n1_MC_eta"]=n1_eta.getVal()
    d["n1_MC_eta_unc"] = n1_eta.getError()
    d["a2_MC_eta"]=a2_eta.getVal()
    d["a2_MC_eta_unc"] = a2_eta.getError()
    d["n2_MC_eta"]=n2_eta.getVal()
    d["n2_MC_eta_unc"] = n2_eta.getError()
    d["MCfit_eta_status"] = mcresult_eta.status()
    d["mean_MC_pi0"] = mean_pi0.getVal()
    d["mean_MC_pi0_unc"] = mean_pi0.getError()
    d["sigma_MC_pi0"] = sigma_pi0.getVal()
    d["sigma_MC_pi0_unc"] = sigma_pi0.getError()
    d["mean_MC_eta"] = mean_eta.getVal()
    d["mean_MC_eta_unc"] = mean_eta.getError()
    d["sigma_MC_eta"] = sigma_eta.getVal()
    d["sigma_MC_eta_unc"] = sigma_eta.getError()
 
    ### Plotting
    #gROOT.ProcessLine(".x lhcbStyle.C")
    #print("--> LHCbStyle loaded")
    
    num_bins = 100

    c0 = TCanvas("c0", "c0", 800, 600)
    frame = massvar.frame(RooFit.Bins(num_bins))
    mcdata_pi0.plotOn(frame, RooFit.Binning(num_bins))
    sigmodel_pi0.plotOn(frame)

    paramMC_pi0 = sigmodel_pi0.getParameters(mcdata_pi0)
    for param in paramMC_pi0:
        print(f"{param.GetName()} (constant: {param.isConstant()})")
    numParamMC_pi0 = sum(1 for p in paramMC_pi0 if not p.isConstant())
    print(f"Number of free parameters: {numParamMC_pi0}")
    chi2_mc_pi0 = frame.chiSquare(numParamMC_pi0)
    print(f"Chi2/NDF (MC) = {chi2_mc_pi0:.2f}")
    d["chi2_ndf_MC_pi0"] = chi2_mc_pi0



    maxi=frame.GetMaximum()
    frame.SetMaximum(maxi*1.2) 

    frame.SetXTitle("m(e^{+}e^{-}#gamma) [MeV/c^{2}]")
    bin_width = (massvar.getMax() - massvar.getMin()) / num_bins
    frame.SetYTitle(f"Candidates / ({bin_width:.1f} MeV/c^{{2}})")
    
    frame.Draw()

    
    Pi0_mass = 134.9768  # MeV/c^2
    line = TLine(Pi0_mass, 0, Pi0_mass, maxi * 1.2)
    line.SetLineColor(kViolet-2)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw()

    leg = TLegend(0.45, 0.35, 0.70, 0.48) 
    leg.SetFillColor(0)   
    leg.SetFillStyle(0)
    leg.SetTextSize(0.05)
    leg.AddEntry(line, f"PDG m(#pi^{{0}}) = {Pi0_mass:.2f} MeV/c^{{2}}", "l")
    leg.Draw()

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

    c0.Update()
    plotofmcfit_pi0 = f"MCfit_pi0_{fd_tag}.pdf"
    c0.SaveAs(plotofmcfit_pi0)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    c0.SaveAs(f"tests/MCFit_pi0_{fd_tag}_{timestamp}.pdf")

    c0.Close()
    print("Saved figure : ", plotofmcfit_pi0)

    mcfile_pi0.Close()

    ################################################################################
    
    c1 = TCanvas("c1", "c1", 800, 650)
    frame1 = massvar.frame(RooFit.Bins(num_bins))
    mcdata_eta.plotOn(frame1, RooFit.Binning(num_bins))
    sigmodel_eta.plotOn(frame1)

    paramMC_eta = sigmodel_eta.getParameters(mcdata_eta)
    for param in paramMC_eta:
        print(f"{param.GetName()} (constant: {param.isConstant()})")
    numParamMC_eta = sum(1 for p in paramMC_eta if not p.isConstant())
    print(f"Number of free parameters: {numParamMC_eta}")
    chi2_mc_eta = frame1.chiSquare(numParamMC_eta)
    print(f"Chi2/NDF (MC) = {chi2_mc_eta:.2f}")
    d["chi2_ndf_MC_eta"] = chi2_mc_eta



    maxi=frame1.GetMaximum()
    frame1.SetMaximum(maxi*1.2) 

    frame1.SetXTitle("m(e^{+}e^{-}#gamma) [MeV/c^{2}]")
    bin_width = (massvar.getMax() - massvar.getMin()) / num_bins
    frame1.SetYTitle(f"Candidates / ({bin_width:.1f} MeV/c^{{2}})")
    
    frame1.Draw()

    
    eta_mass = 547.862  # MeV/c^2
    line = TLine(eta_mass, 0, eta_mass, maxi * 1.2)
    line.SetLineColor(kViolet)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw()

    leg = TLegend(0.15, 0.35, 0.40, 0.48) 
    leg.SetFillColor(0)   
    leg.SetFillStyle(0)
    leg.SetTextSize(0.05)
    leg.AddEntry(line, f"PDG m(#eta) = {eta_mass:.2f} MeV/c^{{2}}", "l")
    leg.Draw()

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
    plotofmcfit_eta = f"MCfit_eta_{fd_tag}.pdf"
    c1.SaveAs(plotofmcfit_eta)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    c1.SaveAs(f"tests/MCFit_eta_{fd_tag}_{timestamp}.pdf")

    c1.Close()
    print("Saved figure : ", plotofmcfit_eta)

    mcfile_eta.Close()
    

    ############################################################

    ### Get the data
    datatreename = "DecayTree"

    print("Open ", datafilename)
    datafile = TFile(datafilename)
    datatree = datafile.Get(datatreename)

    data = RooDataSet("data", "data", listvars, RooFit.Import(datatree))
    print("Data sample has ", data.sumEntries(), "entries")
    
    ### Create total model of signal + background
    # Set shape parameter of the signal shape to the MC parameter from the previsous fit

    #for var in [a1, n1, a2, n2]:
    for var in [a1_pi0, n1_pi0, a2_pi0, n2_pi0, a1_eta, n1_eta, a2_eta, n2_eta]:
    #for var in [xi,rhoL,rhoR]:
        var.setConstant(False)
        var.Print()
   

     
    #m0 = RooRealVar("m0", "m0", 1000, 581, 10000)
    #c  = RooRealVar("c", "c", -20, -100, -1)
    #bkgmodel = RooArgusBG("bkgmodel", "bkgmodel", massvar, m0, c)


    #q2_expr = "TMath:Sqrt( (rho_M*rho_M - (mA + mB)*(mA + mB)) * (rho_M*rho_M - (f- mB)*(mA - mB)) ) / (2*rho_M)"
    #q2_expr = "(rho_M-50)**A * exp(-B*(rho_M-50))"
    #mA = RooRealVar("A", "A", 0.5, 1e-3, 75)
    #mB = RooRealVar("B", "B", 1e-6, 1e-2)
    #bkgmodel = RooGenericPdf("bkgmodel", "bkgmodel", q2_expr, RooArgList(massvar, mA, mB))

   
    #ch0 = RooRealVar("ch0", "ch0", 0.751191741710717, 0.5, 0.9)
    #ch1 = RooRealVar("ch1", "ch1", -0.35691311920190966, -0.5, -0.1)
    #ch2 = RooRealVar("ch2", "ch2", 0.0049958121811287904, 0.001, 0.005)
    #ch3 = RooRealVar("ch3", "ch3", 0.060527392457907926, 0.04, 0.09)
    #ch4 = RooRealVar("ch4", "ch4", -0.05251852686173577, -0.08, -0.02)

    #ch0 = RooRealVar("ch0", "ch0", 0.7525440472321306, 0.75, 0.753)
    #ch1 = RooRealVar("ch1", "ch1", -0.3558356978417028, -0.36, -0.35)
    #ch2 = RooRealVar("ch2", "ch2", 0.004296718425100122, 0.0040, 0.005)
    #ch3 = RooRealVar("ch3", "ch3", 0.06000170413087846, 0.0600, 0.0610)
    #ch4 = RooRealVar("ch4", "ch4", -0.05273207149794215, -0.053, -0.052)


    #ch0 = RooRealVar("ch0", "ch0",  0.7711694296866772)
    #ch1 = RooRealVar("ch1", "ch1", -0.33352209944666444)
    #ch2 = RooRealVar("ch2", "ch2", 0.01111889040916377)
    #ch3 = RooRealVar("ch3", "ch3", 0.0474602837772222)
    #ch4 = RooRealVar("ch4", "ch4", -0.07318425045367961)
    #ch5 = RooRealVar("ch5", "ch5", -0.002736062704096411)




    ch0 = RooRealVar("ch0","ch0", 0.7434911687971489,  0.74, 0.75)
    ch1 = RooRealVar("ch1","ch1", -0.3205845207842097, -0.33, -0.32)
    ch2 = RooRealVar("ch2","ch2", 0.015866619377069718, 0.015, 0.016)
    ch3 = RooRealVar("ch3","ch3", 0.037254298471043404, 0.037, 0.038)
    ch4 = RooRealVar("ch4","ch4",-0.04589722568446994, -0.046, -0.045)

    ch0.setConstant(False)
    ch1.setConstant(False)
    ch2.setConstant(False)
    ch3.setConstant(False)
    ch4.setConstant(False)


    bkgmodel = RooChebychev("bkgmodel", "bkgmodel", massvar, RooArgList(ch0,ch1,ch2,ch3,ch4))

    bkg_positive = RooGenericPdf("bkg_positive", "pow(@0, 2)", RooArgList(bkgmodel))

    #bkgmodel = RooBernstein("bkgmodel", "Bernstein background", massvar, RooArgList(ch0, ch1, ch2, ch3))

    nsig_pi0 = RooRealVar("nsig_pi0", "nsig_pi0", 1, data.numEntries())
    nsig_eta = RooRealVar("nsig_eta", "nsig_eta", 1, data.numEntries())
    nbkg = RooRealVar("nbkg", "nbkg", 1, data.numEntries())

    model = RooAddPdf("model", "model", RooArgSet(sigmodel_pi0, sigmodel_eta, bkgmodel),
                      RooArgSet(nsig_pi0, nsig_eta, nbkg))

    ### Fit the model to the data

    massvar.setBins(100)  
    binnedData = data.binnedClone()

    print("\n> Fit data sample")

    #massvar.setRange("sideband_low", 50, 490)
    

    #bkgmodel.fitTo(data, RooFit.Range("sideband_low"))

    mean_pi0.setRange(130, 131)
    mean_pi0.setVal(130.650240664456)
    mean_pi0.setConstant(False)

    sigma_pi0.setRange(15, 16)
    sigma_pi0.setVal(15.365004872083228)
    sigma_pi0.setConstant(False)

    mean_eta.setRange(554, 555)
    mean_eta.setVal(554.1517079855686)
    mean_eta.setConstant(False)

    sigma_eta.setRange(17, 18)
    sigma_eta.setVal(17.17650651317862)
    sigma_eta.setConstant(False)

    nsig_pi0.setRange(2380000, 2390000)
    nsig_pi0.setVal(2386659.681494587)
    nsig_pi0.setConstant(False)


    nsig_eta.setRange(478000, 479000)
    nsig_eta.setVal(478760.6042866103)
    nsig_eta.setConstant(False)

    nbkg.setRange(21800000, 21900000)
    nbkg.setVal(21810309.161668904)
    nbkg.setConstant(False)

    a1_pi0.setRange(1.0, 2.0)
    a1_pi0.setVal(1.2607008259193968)
    a1_pi0.setConstant(True)

    a2_pi0.setRange(1.0, 2.0)
    a2_pi0.setVal(1.4348632550705949)
    a2_pi0.setConstant(True)

    n1_pi0.setRange(140, 150)
    n1_pi0.setVal(148.17958343873667)
    n1_pi0.setConstant(True)

    n2_pi0.setRange(0.1, 1.0)
    n2_pi0.setVal(0.9369535265475946)
    n2_pi0.setConstant(True)


    a1_eta.setRange(0.1, 1.0)
    a1_eta.setVal(0.49567005180653423)
    a1_eta.setConstant(True)

    a2_eta.setRange(0.5, 1.0)
    a2_eta.setVal(0.6357319533251446)
    a2_eta.setConstant(True)

    n1_eta.setRange(90.0, 150.0)
    n1_eta.setVal(124.01896103245652)
    n1_eta.setConstant(True)

    n2_eta.setRange(1.0, 20.0)
    n2_eta.setVal(14.141213236032844)
    n2_eta.setConstant(True)

   
    print("***mean_pi0 value before InitialDataFit plot:", mean_pi0.getVal())

    ### Plot the result
    c3 = TCanvas("c3", "c3", 800, 600)
    frame3 = massvar.frame(RooFit.Name("Full fit"), RooFit.Title(""),
                           RooFit.Bins(num_bins))

    data.plotOn(frame3, RooFit.Binning(num_bins))

    #frame2.Clone("frame 2")
    model.getComponents().Print("v")

    print("**mean_pi0 value before InitialDataFit plot:", mean_pi0.getVal())


    model.plotOn(frame3,
                 RooFit.Name("model"),
                 RooFit.LineColor(kViolet+4))

    
    model.plotOn(frame3, 
                 RooFit.Components("sigmodel_pi0"),
                 RooFit.Name("sigmodel_pi0"),
                 RooFit.LineColor(kRed+2),
                 RooFit.LineStyle(kDashed))

    model.plotOn(frame3, 
                 RooFit.Components("sigmodel_eta"),
                 RooFit.Name("sigmodel_eta"),
                 RooFit.LineColor(kOrange+7),
                 RooFit.LineStyle(kDashed))

    model.plotOn(frame3,
                 RooFit.Components("bkgmodel"),
                 RooFit.Name("bkgmodel"),
                 RooFit.LineColor(kBlue+2),
                 RooFit.LineStyle(kDashed))

    paramData = model.getParameters(data)
    for param in paramData:
        print(f"{param.GetName()} (constant: {param.isConstant()})")
    numParamData = sum(1 for p in paramData if not p.isConstant())
    print(f"Number of free parameters: {numParamData}")
    chi2_data = frame3.chiSquare(numParamData)
    print(f"Chi2/NDF (Data) = {chi2_data:.2f}")
    d["chi2_ndf_Data"] = chi2_data

    yaxis = frame3.GetYaxis()
    yaxis.SetNoExponent(False)
    yaxis.SetMaxDigits(3)
    #yaxis.SetTitleOffset(1.4)   

    frame3.Draw()

    maxi=frame3.GetMaximum()
    
    #line_pi0 = TLine(Pi0_mass, 0, Pi0_mass, maxi * 1.2)
    #line_pi0.SetLineColor(kGreen + 2)
    #line_pi0.SetLineWidth(2)
    #line_pi0.SetLineStyle(2)
    #line_pi0.Draw()

    #line_eta = TLine(eta_mass, 0, eta_mass, maxi * 1.2)
    #line_eta.SetLineColor(kGreen + 2)
    #line_eta.SetLineWidth(2)
    #line_eta.SetLineStyle(2)
    #line_eta.Draw()

    #leg = TLegend(0.65, 0.75, 0.90, 0.88)
    #leg.SetTextSize(0.05)
    #leg.AddEntry(line_pi0, "PDG m(#pi^{0})", "l")
    #leg.AddEntry(line_eta, "PDG m(#eta)", "l")
    #leg.Draw()



    ### Add a legend
    leg = TLegend(0.68, 0.35, 0.88, 0.63)
    leg.SetFillColor(0)   
    leg.SetFillStyle(0) 
    leg.SetTextSize(0.055)
    leg.AddEntry(frame3.findObject("model"), "Total", "l")
    leg.AddEntry(frame3.findObject("sigmodel_pi0"), "#pi^{0}", "l")
    leg.AddEntry(frame3.findObject("sigmodel_eta"), "#eta", "l")
    leg.AddEntry(frame3.findObject("bkgmodel"), "Background", "l")
    lum = luminosity
    lum_unc = luminosity_uncertainty

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.06) 
    #latex.DrawLatex(0.68, 0.30, f"#scale[0.8]{{\\mathcal{L} = {lum:.1f}\\,fb^{{-1}}}}")
    latex.DrawLatex(0.68, 0.30, f"#scale[0.8]{{#font[300]{{L}} = {lum:.1f} fb^{{-1}}}}")

    leg.Draw()

    maxi3=frame3.GetMaximum()
    frame3.SetMaximum(maxi3*1.2) 

    frame3.SetXTitle("m(e^{+}e^{-}#gamma) [MeV/c^{2}]")
    frame3.SetYTitle(f"Candidates / ({bin_width:.1f} MeV/c^{{2}})")
    

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
    
    c3.Update()
    figurename = f"InitialDataFit_{fd_tag}.pdf"
    if output_dir: 
        figurename = os.path.join(output_dir, figurename)
    c3.SaveAs(figurename)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    c3.SaveAs(f"tests/InitialDataFit_{fd_tag}_{timestamp}.pdf")
   

    print("Save figure : ", figurename)
    c3.Close()



    result = model.fitTo(binnedData, RooFit.Save(True), RooFit.Verbose(False),
                         RooFit.Extended(True))
    #RooFit.ExternalConstraints(rt.RooArgSet(mean_constr)))

    result.Print()

    #chi2_Data = result.chi2()
    #print(f"Chi2 from Data fit: {chi2_Data:.2f}")
    #d["chi2_Data"] = chi2_Data
    
    ### Plot the result
    c2 = TCanvas("c2", "c2", 800, 600)
    frame2 = massvar.frame(RooFit.Name("Full fit"), RooFit.Title(""),
                           RooFit.Bins(num_bins))

    data.plotOn(frame2, RooFit.Binning(num_bins))

    #frame2.Clone("frame 2")
    model.getComponents().Print("v")

    model.plotOn(frame2,
                 RooFit.Name("model"),
                 RooFit.LineColor(kViolet+4))

    
    model.plotOn(frame2, 
                 RooFit.Components("sigmodel_pi0"),
                 RooFit.Name("sigmodel_pi0"),
                 RooFit.LineColor(kRed+2),
                 RooFit.LineStyle(kDashed))

    model.plotOn(frame2, 
                 RooFit.Components("sigmodel_eta"),
                 RooFit.Name("sigmodel_eta"),
                 RooFit.LineColor(kOrange+7),
                 RooFit.LineStyle(kDashed))

    model.plotOn(frame2,
                 RooFit.Components("bkgmodel"),
                 RooFit.Name("bkgmodel"),
                 RooFit.LineColor(kBlue+2),
                 RooFit.LineStyle(kDashed))

    paramData = model.getParameters(data)
    for param in paramData:
        print(f"{param.GetName()} (constant: {param.isConstant()})")
    numParamData = sum(1 for p in paramData if not p.isConstant())
    print(f"Number of free parameters: {numParamData}")
    chi2_data = frame2.chiSquare(numParamData)
    print(f"Chi2/NDF (Data) = {chi2_data:.2f}")
    d["chi2_ndf_Data"] = chi2_data

    yaxis = frame2.GetYaxis()
    yaxis.SetNoExponent(False)
    yaxis.SetMaxDigits(3)
    #yaxis.SetTitleOffset(1.4)   

    frame2.Draw()

    maxi=frame2.GetMaximum()
    
    #line_pi0 = TLine(Pi0_mass, 0, Pi0_mass, maxi * 1.2)
    #line_pi0.SetLineColor(kGreen + 2)
    #line_pi0.SetLineWidth(2)
    #line_pi0.SetLineStyle(2)
    #line_pi0.Draw()

    #line_eta = TLine(eta_mass, 0, eta_mass, maxi * 1.2)
    #line_eta.SetLineColor(kGreen + 2)
    #line_eta.SetLineWidth(2)
    #line_eta.SetLineStyle(2)
    #line_eta.Draw()

    #leg = TLegend(0.65, 0.75, 0.90, 0.88)
    #leg.SetTextSize(0.05)
    #leg.AddEntry(line_pi0, "PDG m(#pi^{0})", "l")
    #leg.AddEntry(line_eta, "PDG m(#eta)", "l")
    #leg.Draw()



    ### Add a legend
    leg = TLegend(0.68, 0.35, 0.88, 0.63)
    leg.SetFillColor(0)   
    leg.SetFillStyle(0) 
    leg.SetTextSize(0.055)
    leg.AddEntry(frame2.findObject("model"), "Total", "l")
    leg.AddEntry(frame2.findObject("sigmodel_pi0"), "#pi^{0}", "l")
    leg.AddEntry(frame2.findObject("sigmodel_eta"), "#eta", "l")
    leg.AddEntry(frame2.findObject("bkgmodel"), "Background", "l")
    lum = luminosity
    lum_unc = luminosity_uncertainty

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.06) 
    #latex.DrawLatex(0.68, 0.30, f"#scale[0.8]{{\\mathcal{L} = {lum:.1f}\\,fb^{{-1}}}}")
    latex.DrawLatex(0.68, 0.30, f"#scale[0.8]{{#font[300]{{L}} = {lum:.1f} fb^{{-1}}}}")

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

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    c2.SaveAs(f"tests/DataFit_{fd_tag}_{timestamp}.pdf")
   

    print("Save figure : ", figurename)
    c2.Close()

    ### Print the result
    print("\n --- Fit summary ---")
    print(" Nsig_pi0 = ", nsig_pi0.getVal(), "+-", nsig_pi0.getError())
    print(" Nsig_eta = ", nsig_eta.getVal(), "+-", nsig_eta.getError())
    print(" Nbkg = ", nbkg.getVal(), "+-", nbkg.getError())
    print(" mean_pi0 = ", mean_pi0.getVal(), "+-", mean_pi0.getError())
    print(" mean_eta = ", mean_eta.getVal(), "+-", mean_eta.getError())
    print(" sigma_pi0 = ", sigma_pi0.getVal(), "+-", sigma_pi0.getError())
    print(" sigma_eta = ", sigma_eta.getVal(), "+-", sigma_eta.getError())

    ### Save result from above in a dictionnary

    d["sig/bkg_pi0"] = nsig_pi0.getVal()/nbkg.getVal()
    d["Nsig_pi0"] = nsig_pi0.getVal()
    d["Nsig_pi0_unc"] = nsig_pi0.getError()
    d["Nsig_eta"] = nsig_eta.getVal()
    d["Nsig_eta_unc"] = nsig_eta.getError()
    d["Luminosity"] = lum
    d["Luminosity_unc"] = lum_unc
    d["Nsig_pi0/L"] = nsig_pi0.getVal()/lum
    d["Nsig_pi0/L_unc"] = float((nsig_pi0.getVal()/lum)*np.sqrt((nsig_pi0.getError()/nsig_pi0.getVal())**2 + (lum_unc/lum)**2))
    #d["Nsig_eta/L"] = nsig_eta.getVal()/lum
    #d["Nsig_eta/L_unc"] = float((nsig_eta.getVal()/lum)*np.sqrt((nsig_eta.getError()/nsig_eta.getVal())**2 + (lum_unc/lum)**2))
    # ToDo: add other fit variable values
    d["mean_pi0"] =mean_pi0.getVal()
    d["mean_pi0_unc"] = mean_pi0.getError()
    d["sigma_pi0"] =sigma_pi0.getVal()
    d["sigma_pi0_unc"] = sigma_pi0.getError()
    d["mean_eta"] =mean_eta.getVal()
    d["mean_eta_unc"] = mean_eta.getError()
    d["sigma_eta"] =sigma_eta.getVal()
    d["sigma_eta_unc"] = sigma_eta.getError()
    d["a1_pi0"] = a1_pi0.getVal()
    d["a1_pi0_unc"] = a1_pi0.getError()
    d["a2_pi0"] = a2_pi0.getVal()
    d["a2_pi0_unc"] = a2_pi0.getError()
    d["n1_pi0"] = n1_pi0.getVal()
    d["n1_pi0_unc"] = n1_pi0.getError()
    d["n2_pi0"] = n2_pi0.getVal()
    d["n2_pi0_unc"] = n2_pi0.getError()
    d["a1_eta"]=a1_eta.getVal()
    d["a1_eta_unc"] = a1_eta.getError()
    d["n1_eta"]=n1_eta.getVal()
    d["n1_eta_unc"] = n1_eta.getError()
    d["a2_eta"]=a2_eta.getVal()
    d["a2_eta_unc"] = a2_eta.getError()
    d["n2_eta"]=n2_eta.getVal()
    d["n2_eta_unc"] = n2_eta.getError()
    #d["mA"] = mA.getVal()
    #d["mA_unc"] = mA.getError()
    #d["mB"] = mB.getVal()
    #d["mB_unc"] = mB.getError()
    d["Nbkg"] = nbkg.getVal()
    d["Nbkg_unc"] = nbkg.getError()
    d["Nbkg/L"] = nbkg.getVal()/lum
    d["Nbkg/L_unc"] = float((nbkg.getVal()/lum)*np.sqrt((nbkg.getError()/nbkg.getVal())**2 + (lum_unc/lum)**2))
    d["DataFit_status"] = result.status()
    d["ch0"] = ch0.getVal()
    d["ch0_unc"] = ch0.getError()
    d["ch1"] = ch1.getVal()
    d["ch1_unc"] = ch1.getError()
    d["ch2"] = ch2.getVal()
    d["ch2_unc"] = ch2.getError()
    d["ch3"] = ch3.getVal()
    d["ch3_unc"] = ch3.getError()
    d["ch4"] = ch4.getVal()
    d["ch4_unc"] = ch4.getError()
    #d["ch5"] = ch5.getVal()
    #d["ch5_unc"] = ch5.getError()
    #d["ch6"] = ch6.getVal()
    #d["ch6_unc"] = ch6.getError()
    
    # Save result dictionnary into a yaml file
    yfile = "FitResults.yml"
    if output_dir: 
        yfile = os.path.join(output_dir, yfile)
    with open(yfile, 'w') as outfile:
        yaml.dump(d, outfile, default_flow_style=False)
        print("\nSaved yml file in: ", yfile)


    if not os.path.exists("tests"):
        os.makedirs("tests")
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    yfile_test = os.path.join("tests", f"FitResults_{timestamp}.yml")
    with open(yfile_test, 'w') as test_outfile:
        yaml.dump(d, test_outfile, default_flow_style=False)
        print("Also saved yml file in: ", yfile_test)
    
    return True


################# MAIN ###############################################       
    

if __name__ == "__main__":

    dir_MC = "/eos/lhcb/user/c/castillj/DarkPhoton-samples/MC"
    if not os.path.exists(dir_MC):
        os.makedirs(dir_MC)

    print(">> START")

    start_time = time.time()

    print("\n>> Get selection for MC file")
    mcsel = get_selection("truM_"+selstring, fd_tag, mcfile)

    print("Combine files in one MC sample")
    mcfilename_pi0 = setup_files(data_type="MC_pi0", track_type=track_type,
                             sel=mcsel, samesign=samesign, filenbr=filenumber, output_dir=dir_MC)
    mcfilename_eta = setup_files(data_type="MC_eta", track_type=track_type,
                             sel=mcsel, samesign=samesign, filenbr=filenumber, output_dir=dir_MC)

    print(mcfilename_pi0)
    print(mcfilename_eta)


    for dataBlock in dataBlock_dict.keys():
        print("\n >> Start with data block", dataBlock)
        data_string, data_polarity, data_year, fillnbr_selection, lumi, lumi_unc  = dataBlock_dict[dataBlock]
        print(f"name: {data_string}, polarity: {data_polarity}, year: {data_year}, fill number: {fillnbr_selection}")

        if filenumber == -1:
            block_dir_data = f"/eos/lhcb/user/c/castillj/DarkPhoton-samples/block_{dataBlock}"
        else:
            block_dir_data = f"block_{dataBlock}"

        block_dir_fit = f"block_{dataBlock}"

        expected_files = [f"{block_dir_fit}/DataFit_{fd_tag}.pdf",
                        f"{block_dir_fit}/FitResults.yml"]

        if os.path.exists(block_dir_fit) and all(os.path.exists(f) for f in expected_files):
            print(f">> Skipping block {dataBlock}: already processed.")
            continue
        else:
            if not os.path.exists(block_dir_fit):
                os.mkdir(block_dir_fit)

        #if filenumber == -1:
        #    expected_file = os.path.join(block_dir_data, f"Data_{fd_tag}"+"_SS"*samesign+".root")
        #else:
        #    expected_file = os.path.join(block_dir_data, f"Data_{fd_tag}"+"_SS"*samesign+".root") 

        expected_file = os.path.join(block_dir_data, f"Data_{fd_tag}"+"_SS"*samesign+".root") 

        if not os.path.exists(expected_file):
            
            print("\n>> Get selection for data file")
            datasel = get_selection(selstring, fd_tag, "", fillnbr= fillnbr_selection, data_year=data_year)

            print("Combine files in one data sample")
            datafilename = setup_files(data_type="Data", track_type=track_type,
                                    sel=datasel, samesign=samesign, filenbr=filenumber,
                                    data_name = data_string, data_polarity = data_polarity, data_year = data_year, output_dir=block_dir_data)
            print(datafilename)
            
            print("\n>> Start fitting")
            fit_ok = fit_data(mcfilenames=[mcfilename_pi0,mcfilename_eta], datafilename=datafilename,
                        track_type=track_type, samesign=samesign,output_dir=block_dir_fit, luminosity=lumi, luminosity_uncertainty=lumi_unc)

            #fit_ok = fit_data(mcfilenames=[mcfilename_pi0], datafilename=datafilename,
            #            track_type=track_type, samesign=samesign,output_dir=block_dir_fit, luminosity=lumi, luminosity_uncertainty=lumi_unc)

            print(f">> BLOCK {dataBlock} DONE")

        else:
            datafilename = expected_file
            print(datafilename)
            print("\n>> Start fitting")
            fit_ok = fit_data(mcfilenames=[mcfilename_pi0,mcfilename_eta], datafilename=datafilename,
                        track_type=track_type, samesign=samesign,output_dir=block_dir_fit, luminosity=lumi, luminosity_uncertainty=lumi_unc)

            #fit_ok = fit_data(mcfilenames=[mcfilename_pi0], datafilename=datafilename,
            #            track_type=track_type, samesign=samesign,output_dir=block_dir_fit, luminosity=lumi, luminosity_uncertainty=lumi_unc)

            print(f">> BLOCK {dataBlock} DONE")



    print(">> DONE")

    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"Time: {elapsed_time:.2f} s")

## EOF
