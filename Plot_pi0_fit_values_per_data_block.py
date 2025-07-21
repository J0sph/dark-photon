#!/usr/bin/env python3

################# Plot_pi0_fit_values_per_data_block.py #####################        
#                                                                           #         
#    Author: Felicia Volle for Joseph                                       #         
#    Date:   27/06/2025                                                     #         
#    Description: plot fit values of pi0 mass fit per data block            #         
#                                                                           #         
#    Run via: lb-conda default python Plot_pi0_fit_values_per_data_block.py #         
#                                                                           #        
#############################################################################

################# IMPORT ############################################

import os, sys, yaml
import numpy as np
from ROOT import TGraphErrors, TH1F, TCanvas, TLegend, gROOT, gStyle, kBlack, kRed, kBlue, kPink, kCyan, kGreen, kWhite, kAzure, kOrange, kViolet, kSpring, kMagenta, TLine, TPaveText, TBox
from array import array 
from Fit_pizero_peak import dataBlock_dict
################# SETTING ###########################################                               

data_periods = list(dataBlock_dict.keys())

################# FUNCTIONS #########################################                                

def get_dict(ymlpath):
    with open(ymlpath, "r") as outfile:
        d = yaml.load(outfile, Loader=yaml.SafeLoader)
    print("FitResult dictionnary opened from : "+ymlpath)
    return d

def Plotting(observable, values, uncertainties, maxi, mini, dataBlock_list, theory_line=None, mean_line=None, mean_uncertainty=None, LHCblabel=""):
    nbr = len(values)

    c = TCanvas('c', 'c', 600, 450)
    #c.SetLeftMargin(0.15)   
    c.SetBottomMargin(0.18) 

   
    x    = array('f', [m+0.5 for m in range(nbr)])
    xerr = array('f', [0]*(nbr))
    y    = array('f', values)
    yerr = array('f', uncertainties)

    name = observable
        
    h = TH1F(name, name, nbr, 0, nbr)
    for iBin in range(h.GetNbinsX()):
        label = str(dataBlock_list[iBin])
        h.GetXaxis().SetBinLabel(iBin+1, label)


    h.SetMaximum( maxi + (maxi-mini)*2 )
    h.SetMinimum( mini - (maxi-mini)*1 )

    h.GetXaxis().SetTitle("Data blocks")
    h.GetXaxis().SetTitleOffset(0.7)
    h.GetXaxis().SetLabelSize(1.5*0.04)
    #h.GetXaxis().SetLabelOffset(0.010)
    
    if observable == "mean" or observable == "sigma":
        y_title = f"{observable} (MeV/c^{{2}})"
    else:
        y_title = observable  

    h.GetYaxis().SetTitle(y_title)
    h.GetYaxis().SetTitleOffset(1)
    h.GetYaxis().SetNoExponent(False)

    h.Draw()

    #TLegend(x1, y1, x2, y2)
    legend = TLegend(0.40, 0.20, 0.80, 0.40) 
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.055) 

    if mean_uncertainty is not None:
        band = TBox(0, mean_line - mean_uncertainty, nbr, mean_line + mean_uncertainty)
        band.SetFillColorAlpha(kGreen -9, 1)  
        band.SetFillStyle(1001)           
        band.Draw("same")
        legend.AddEntry(band, f"\\pm {mean_uncertainty:.2f} MeV/c^{{2}}", "f")

    if mean_line is not None:
        line_mean = TLine(0, mean_line, nbr, mean_line)
        line_mean.SetLineColor(kGreen + 2)
        line_mean.SetLineStyle(1)  
        line_mean.SetLineWidth(2)
        line_mean.Draw("same")
        legend.AddEntry(line_mean, f"Mean = {mean_line:.2f} MeV/c^{{2}}", "l")
     
    if theory_line is not None:
        line = TLine(0, theory_line, nbr, theory_line)
        line.SetLineColor(kRed)
        line.SetLineStyle(2)  
        line.SetLineWidth(2)
        line.Draw("same")     

    graph = TGraphErrors(nbr, x, y, xerr, yerr)
    graph.SetMarkerColor(kAzure+1)
    graph.SetLineColor(kAzure+1)

    graph.Draw("P")
 

    if observable == "mean" or observable == "sigma": 
        legend.AddEntry(graph, "Data", "lep")
        if theory_line is not None:
            legend.AddEntry(line, f"PDG m(#pi^{{0}}) = {theory_line:.2f} MeV/c^{{2}}", "l")
        legend.Draw()


    lhcbName = TPaveText(gStyle.GetPadLeftMargin() + 0.03,
                         0.87 - gStyle.GetPadTopMargin(),
                         gStyle.GetPadLeftMargin() + 0.38,
                         0.95 - gStyle.GetPadTopMargin(),
                         "BRNDC")
    if LHCblabel == "LHCb_data":
        lhcbName.AddText("LHCb Preliminary")
    elif LHCblabel == "LHCb_MC":
        lhcbName.AddText("LHCb Simulation")
    else:
        print("Invalid LHCb label key word : ", LHCblabel)
    lhcbName.SetFillColor(0)
    lhcbName.SetFillStyle(0)
    lhcbName.SetTextAlign(12)
    lhcbName.SetTextSize(0.055)
    lhcbName.SetBorderSize(0)
    lhcbName.Draw()

    c.Update()

    if not os.path.exists("DataBlockSummary"):
        os.mkdir("DataBlockSummary")
    c.SaveAs(f"DataBlockSummary/{observable}.pdf")
    c.Close()
    return 0
    

################# MAIN ###############################################                                 

if __name__ == '__main__':
    gROOT.ProcessLine(".x lhcbStyle.C")
    print("--> LHCbStyle loaded")

    result_val_dict = {}
    maxi, mini = 0, 0


    # True, generated observable values from Anja
    ymlfiles, datablocks = [], []
    for data in data_periods:
        ymlfile   = get_dict(f"block_{data}/FitResults.yml")
        datablock = data
        ymlfiles.append(ymlfile)
        datablocks.append(datablock)

    Pi0_mass = 134.9768  # MeV/c^2


    for observable in ["mean", "sigma", "Nsig"]:
        print("\nGet values of observable:", observable)
        values, uncertainties, events = [], [], []
        for i in range(len(ymlfiles)):
            fitvaldict = ymlfiles[i]
            datablock = datablocks[i]

            fitval = fitvaldict[f"{observable}"]
            fitunc = fitvaldict[f"{observable}_unc"]
            event =  fitvaldict["Nsig"]
            
            values.append(fitval)
            uncertainties.append(fitunc)
            events.append(event)
            

            print(f"For data block {datablock} the {observable} is {fitval} +- {fitunc}")
        maxi = max(values)
        mini = min(values)
        
        if observable == "mean":
            total_events = np.sum(events)
            #mean_mass = np.sum((events / total_events) * values)
            #std_error = np.sqrt(np.sum(events * (values - mean_mass)**2)/ total_events)
            weights = 1/(np.array(uncertainties)**2)
            mean_mass = np.sum(weights*values)/np.sum(weights)
            error_mass = 1/np.sqrt(np.sum(weights))
            Plotting(observable, values, uncertainties, maxi, mini, data_periods, theory_line=Pi0_mass, mean_line=mean_mass, mean_uncertainty=error_mass, LHCblabel="LHCb_data")

        elif observable == "sigma":
            total_events = np.sum(events)
            #mean_sigma = np.sum((events / total_events) * values)
            #std_error_sigma = np.sqrt(np.sum(events * (values - mean_sigma)**2)/ total_events)
            weights = 1/(np.array(uncertainties)**2)
            mean_sigma = np.sum(weights*values)/np.sum(weights)
            error_sigma = 1/np.sqrt(np.sum(weights))
            Plotting(observable, values, uncertainties, maxi, mini, data_periods, mean_line=mean_sigma, mean_uncertainty=error_sigma, LHCblabel="LHCb_data")
            
        else:
            Plotting(observable, values, uncertainties, maxi, mini, data_periods, LHCblabel="LHCb_data")

    print("\nDONE!")



#EOF