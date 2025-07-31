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
from ROOT import TLatex, TLine, RooBernstein, RooBukinPdf, RooNovosibirsk, RooChebychev, RDataFrame, RDF, TLegend, TCanvas, TFile, TChain, RooRealVar, RooArgList, RooArgSet, RooGaussian, RooCBShape, RooLinkedList, RooArgusBG, RooAddPdf, RooFit, RooDataSet, RooCrystalBall, RooGenericPdf, gROOT, EnableImplicitMT, kRed, kGreen, kViolet, kBlue, kDashed, TPaveText, gStyle
import pickle
import yaml
from datetime import datetime
import uproot
import numpy as np

gROOT.SetBatch(True)
EnableImplicitMT()
#PLOTTING
gROOT.ProcessLine(".x lhcbStyle.C")
print("--> LHCbStyle loaded")


################# SETTING #########################################


dataBlock_dict = {#4: ["data_24c2a_magdown_qee", "magdown", "2024", "FillNumber >= 9808 && FillNumber <= 9910" ],
                  #3 : ["data_24c2a_magup_qee", "magup", "2024", "FillNumber >= 9911 && FillNumber <= 9943" ],
                  #2 : ["data_24c2a_magup_qee", "magup", "2024", "FillNumber >= 9945 && FillNumber <= 9978" ],
                  #1 : ["data_24c2a_magup_qee", "magup", "2024", "FillNumber >= 9982 && FillNumber <= 10056" ],
                  #5 : ["data_24c3a_magup_qee", "magup", "2024", "FillNumber >= 10059 && FillNumber <= 10102"],
                  #6: ["data_24c3a_magdown_qee", "magdown", "2024", "FillNumber >= 10104 && FillNumber <= 10190" ],
                  #7: ["data_24c4a_magdown_qee", "magdown", "2024", "FillNumber >= 10197 && FillNumber <= 10213" ],
                  8 : ["data_24c4a_magup_qee", "magup", "2024", "FillNumber >= 10214 && FillNumber <= 10232" ],
                  #9: ["data_25c1_magdown_qee", "magdown", "2025", "FillNumber >= 10489 && FillNumber <= 10732" ]
}       



track_type   = "Prompt" # Prompt or Displaced
samesign = False
fd_tag = track_type
if samesign == True:
    fd_tag += "SS"
    
mcfile = "Dalitz" # 3-body decay of pi0, as in SM. 
selstring = "hlt1_BremOverlap_Ecal_PID_Kin_Rho_Topo"
filenumber = 1

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
    if fillnbr != None:s
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
    
    ofilename = f"{data_type}_{track_type}"+"_SS"*samesign+".root"

    if output_dir:
        ofilename = os.path.join(output_dir, ofilename)

    if os.path.exists(ofilename):
        return ofilename


    print(f"Merging files into {ofilename}")
    ifilenames = get_pfns(data_type=data_type, data_name= data_name, data_polarity=data_polarity, data_year=data_year)
    #if data_type == "Data" and not samesign:
        #ifilenames = ifilenames[:2]
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


def fit_data(mcfilename, datafilename, track_type, samesign=False, output_dir=None, luminosity=1, luminosity_uncertainty=0):
    mctreename = "DecayTree"

    mcfile = TFile(mcfilename)
    mctree = mcfile.Get(mctreename)

    massvar = RooRealVar("rho_M", "rho_M", 50, 600) #50, 400 #70,400/500
    listvars = RooArgSet(massvar)

    mcdata = RooDataSet("mcdata", "mcdata", listvars, RooFit.Import(mctree))
    print("MC data has", mcdata.sumEntries(), "entries")
    
    ### Create the fit model
    #a1 = RooRealVar("a1", "a1", 1.1, 0.1, 5)
    #a2 = RooRealVar("a2", "a2", 1.45, 0.1, 5)
    #n1 = RooRealVar("n1", "n1", 144, 100, 200)
    #n2 = RooRealVar("n2", "n2", 1.22, 0.1, 10)

    #mean = RooRealVar("mean", "mean", 130, 120, 155)

    #sigma = RooRealVar("sigma", "sigma", 15.38, 2, 20)

    #sigmodel = RooCrystalBall("sigmodel", "sigmodel", massvar, mean, sigma, a1, n1, a2, n2)

    a = RooRealVar("alpha", "alpha", -1.60, -2, 2)
    n  = RooRealVar("n", "n", 1.12, 0.1, 2) 
    sigma = RooRealVar("sigma", "sigma", 18.22, 2, 20)
    mean  = RooRealVar("mean", "mean", 130, 120, 140)
   
    sigmodel = RooCrystalBall("sigmodel", "sigmodel", massvar, mean, sigma, a, n, False)


    #mean2  = RooRealVar("mean2", "mean2", 540, 500, 600)  
    #sigma2 = RooRealVar("sigma2", "sigma2", 15.0, 5, 30)
    #a2 = RooRealVar("alpha2", "alpha2", -1.0, -5, 0)
    #n2     = RooRealVar("n2", "n2", 1.5, 0.1, 5)

    #sigmodel2 = RooCrystalBall("sigmodel2", "sigmodel2", massvar, mean2, sigma2, a2, n2)
    
    # Fit options we can think about later
    #mean_constr = RooGaussian("mean_constr", "mean_constr", mean,
    #                             RooFit.RooConst(135), RooFit.RooConst(5))



   #mean = RooRealVar("mean", "mean", 135, 120, 150)    
   #sigma = RooRealVar("sigma", "sigma", 10, 0.1, 100)             
   #tau = RooRealVar("tau", "tau", 1, -10.0, 10.0) 

   #sigmodel = RooNovosibirsk("novosibirsk", "Novosibirsk Function", massvar, mean, sigma, tau)  


    # Parámetros del Bukin PDF
    #mean = RooRealVar("mean", "mean", 130, 120, 155)         
    #sigma = RooRealVar("sigma", "sigma", 15.38, 2, 20)                
    #xi = RooRealVar("xi", "xi", -0.1, -2, 0)                
    #rhoL = RooRealVar("rhoL", "rhoL", 0.1, 0.001, 1.0)          
    #rhoR = RooRealVar("rhoR", "rhoR", 0.1, 0.001, 1.0)          

    # Definir modelo de señal con Bukin
    #sigmodel = RooBukinPdf("sigmodel", "sigmodel", massvar, mean, sigma, xi, rhoL, rhoR)

    ### Create the fit of the model to the mc
    print("\n> Fit MC sample")
    mcresult = sigmodel.fitTo(mcdata, RooFit.Save(True), RooFit.Verbose(False))
    #, RooFit.Minos(1),
    #RooFit.ExternalConstraints(rt.RooArgSet(mean_constr)))

    mcresult.Print()

    #chi2_MC = mcresult.chi2()
    #print(f"Chi2 from MC fit: {chi2_MC:.2f}")
    

    ### Save result from above in a dictionnary
    d = {}
    #d["a1"] = a1.getVal()
    #d["a1_unc"] = a1.getError()
    #d["a2"] = a2.getVal()
    #d["a2_unc"] = a2.getError()
    #d["n1"] = n1.getVal()
    #d["n1_unc"] = n1.getError()
    #d["n2"] = n2.getVal()
    #d["n2_unc"] = n2.getError()
    d["a"]=a.getVal()
    d["a_unc"] = a.getError()
    d["n"]=n.getVal()
    d["n_unc"] = n.getError()
    d["MCfit_status"] = mcresult.status()
    #d["chi2_MC"] = chi2_MC
    #d["Tau"] = tau.getVal()
    #d["Tau_unc"] = tau.getError()
    #d["xi"] = xi.getVal()
    #d["xi_unc"] = xi.getError()
    #d["rhoL"] = rhoL.getVal()
    #d["rhoL_unc"] = rhoL.getError()
    #d["rhoR"] = rhoR.getVal()
    #d["rhoR_unc"] = rhoR.getError()

    
    ### Plotting
    #gROOT.ProcessLine(".x lhcbStyle.C")
    #print("--> LHCbStyle loaded")
    
    num_bins = 100

    c1 = TCanvas("c1", "c1", 800, 600)
    frame = massvar.frame(RooFit.Binning(num_bins))
    mcdata.plotOn(frame, RooFit.Binning(num_bins))
    sigmodel.plotOn(frame)

    paramMC = sigmodel.getParameters(mcdata)
    for param in paramMC:
        print(f"{param.GetName()} (constant: {param.isConstant()})")
    numParamMC = sum(1 for p in paramMC if not p.isConstant())
    print(f"Number of free parameters: {numParamMC}")
    chi2_mc = frame.chiSquare(numParamMC)
    print(f"Chi2/NDF (MC) = {chi2_mc:.2f}")
    d["chi2_ndf_MC"] = chi2_mc



    maxi=frame.GetMaximum()
    frame.SetMaximum(maxi*1.2) 

    frame.SetXTitle("m(e^{+}e^{-}#gamma) [MeV/c^{2}]")
    bin_width = (massvar.getMax() - massvar.getMin()) / num_bins
    frame.SetYTitle(f"Candidates / ({bin_width:.1f} MeV/c^{{2}})")
    
    frame.Draw()

    # Línea vertical del valor teórico
    Pi0_mass = 134.9768  # MeV/c^2
    line = TLine(Pi0_mass, 0, Pi0_mass, maxi * 1.2)
    line.SetLineColor(kGreen + 2)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw()

    # Leyenda
    #leg = TLegend(0.65, 0.75, 0.90, 0.88)  # Ajusta según espacio disponible
    #leg.SetTextSize(0.05)
    #leg.AddEntry(line, "PDG m(#pi^{{0}})", "l")
    #leg.Draw()

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

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    c1.SaveAs(f"tests/MCFit_{fd_tag}_{timestamp}.pdf")

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

    #for var in [a1, n1, a2, n2]:
    for var in [a, n]:
    #for var in [xi,rhoL,rhoR]:
        var.setConstant(True)
        var.Print()

    
        
    #m0 = RooRealVar("m0", "m0", 1000, 581, 10000)
    #c  = RooRealVar("c", "c", -20, -100, -1)

    #bkgmodel = RooArgusBG("bkgmodel", "bkgmodel", massvar, m0, c)

    #q2_expr = "TMath:Sqrt( (rho_M*rho_M - (mA + mB)*(mA + mB)) * (rho_M*rho_M - (f- mB)*(mA - mB)) ) / (2*rho_M)"


    #q2_expr = "(rho_M-50)**A * exp(-B*(rho_M-50))"
    #mA = RooRealVar("A", "A", 0.5, 1e-3, 75)
    #mB = RooRealVar("B", "B", 1e-6, 1e-2)
    #bkgmodel = RooGenericPdf("bkgmodel", "bkgmodel", q2_expr, RooArgList(massvar, mA, mB))

    ch0 = RooRealVar("ch0", "ch0", 0.73)
    ch1 = RooRealVar("ch1", "ch1", -0.29)
    ch2 = RooRealVar("ch2", "ch2", 0.02)
    ch3 = RooRealVar("ch3", "ch3", -0.01)
    ch4 = RooRealVar("ch4", "ch4", -0.06)
    ch5 = RooRealVar("ch5", "ch5", -0.03)
    ch6 = RooRealVar("ch6", "ch6", -0.04)

    ch0.setConstant(False)
    ch1.setConstant(False)
    ch2.setConstant(False)
    ch3.setConstant(False)
    ch4.setConstant(False)
    ch5.setConstant(False)
    ch6.setConstant(False)
    

    bkgmodel = RooChebychev("bkgmodel", "bkgmodel", massvar, RooArgList(ch0,ch1,ch2,ch3,ch4,ch5,ch6))

    #bkgmodel = RooBernstein("bkgmodel", "bkgmodel", massvar, RooArgList(ch0, ch1, ch2))

    nsig = RooRealVar("nsig", "nsig", 1, data.numEntries())
    nbkg = RooRealVar("nbkg", "nbkg", 1, data.numEntries())

    model = RooAddPdf("model", "model", RooArgSet(sigmodel, bkgmodel),
                      RooArgSet(nsig, nbkg))

    ### Fit the model to the data

    
    massvar.setBins(100)  
    binnedData = data.binnedClone()

    print("\n> Fit data sample")



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
                           RooFit.Binning(num_bins))

    data.plotOn(frame2, RooFit.Binning(num_bins))

    #frame2.Clone("frame 2")
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
    

    line = TLine(Pi0_mass, 0, Pi0_mass, maxi * 1.2)
    line.SetLineColor(kGreen + 2)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw()

    # Leyenda
    #leg = TLegend(0.65, 0.75, 0.90, 0.88)  # Ajusta según espacio disponible
    #leg.SetTextSize(0.05)
    #leg.AddEntry(line, "PDG m(#pi^0)", "l")
    #leg.Draw()


    ### Add a legend
    leg = TLegend(0.68, 0.35, 0.88, 0.53)
    leg.SetTextSize(0.055)
    leg.AddEntry(frame2.findObject("model"), "Total", "l")
    leg.AddEntry(frame2.findObject("sigmodel"), "#pi^{0}", "l")
    leg.AddEntry(frame2.findObject("bkgmodel"), "Background", "l")
    lum = luminosity
    lum_unc = luminosity_uncertainty

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.06) 
    latex.DrawLatex(0.68, 0.30, f"#scale[0.8]{{#it{{L}} = {lum:.1f}\\,fb^{{-1}}}}")

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
    print(" Nsig = ", nsig.getVal(), "+-", nsig.getError())
    print(" Nbkg = ", nbkg.getVal(), "+-", nbkg.getError())
    print(" mean = ", mean.getVal(), "+-", mean.getError())
    #print(" sigma = ", sigma.getVal(), "+-", sigma.getError())

    ### Save result from above in a dictionnary

    d["Nsig"] = nsig.getVal()
    d["Nsig_unc"] = nsig.getError()
    d["Luminosity"] = lum
    d["Luminosity_unc"] = lum_unc
    d["Nsig/L"] = nsig.getVal()/lum
    d["Nsig/L_unc"] = float((nsig.getVal()/lum)*np.sqrt((nsig.getError()/nsig.getVal())**2 + (lum_unc/lum)**2))
    # ToDo: add other fit variable values
    d["mean"] =mean.getVal()
    d["mean_unc"] = mean.getError()
    d["sigma"] =sigma.getVal()
    d["sigma_unc"] = sigma.getError()
    #d["mA"] = mA.getVal()
    #d["mA_unc"] = mA.getError()
    #d["mB"] = mB.getVal()
    #d["mB_unc"] = mB.getError()
    d["Nbkg"] = nbkg.getVal()
    d["Nbkg_unc"] = nbkg.getError()
    d["Nbkg/L"] = nsig.getVal()/lum
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
    d["ch5"] = ch5.getVal()
    d["ch5_unc"] = ch5.getError()
    d["ch6"] = ch6.getVal()
    d["ch6_unc"] = ch6.getError()
    
    
    
    
    

    
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

    print("\n>> Get selection for MC file")
    mcsel = get_selection("truM_"+selstring, fd_tag, mcfile)

    print("Combine files in one MC sample")
    mcfilename = setup_files(data_type="MC_pi0", track_type=track_type,
                             sel=mcsel, samesign=samesign, filenbr=filenumber, output_dir=dir_MC)
    mcfilename2 = setup_files(data_type="MC_eta", track_type=track_type,
                             sel=mcsel, samesign=samesign, filenbr=filenumber, output_dir=dir_MC)

    print(mcfilename)
    print(mcfilename2)


    for dataBlock in dataBlock_dict.keys():
        print("\n >> Start with data block", dataBlock)
        data_string, data_polarity, data_year, fillnbr_selection  = dataBlock_dict[dataBlock]
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
            fit_ok = fit_data(mcfilename=mcfilename, datafilename=datafilename,
                        track_type=track_type, samesign=samesign,output_dir=block_dir_fit)

            print(f">> BLOCK {dataBlock} DONE")

        else:
            datafilename = expected_file
            print(datafilename)
            print("\n>> Start fitting")
            fit_ok = fit_data(mcfilename=mcfilename, datafilename=datafilename,
                        track_type=track_type, samesign=samesign,output_dir=block_dir_fit)

            print(f">> BLOCK {dataBlock} DONE")



    print(">> DONE")

## EOF
