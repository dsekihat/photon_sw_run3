import os, sys, shutil
import math
import argparse
import numpy as np
import ctypes
import yaml
import datetime
import ROOT
from ROOT import TFile, THashList, TH1F, TF1, TCanvas, TPad, TLegend, TPaveText
from ROOT import kWhite, kBlack
from ROOT import gStyle, gROOT, gSystem
from histo_manager import slice_histogram, rebin_histogram
from nm_fitter import NMFitter

# python analyze_material_budget.py -d AnalysisResults_HL_147811.root -m AnalysisResults_HL_147812.root -c configs/config_pp_13.6TeV_LHC22f_material.yml;
# python analyze_material_budget.py -d AnalysisResults_HL_152286_LHC22f.root -m AnalysisResults_HL_152154.root -c configs/config_pp_13.6TeV_LHC22f_material.yml;

parser = argparse.ArgumentParser('script to analyze mateial budget for PCM');
parser.add_argument("-d", "--input_data" , default="AnalysisResults_data.root", type=str, help="path to the root file of data which you want to analyze", required=True)
parser.add_argument("-m", "--input_mc" , default="AnalysisResults_mc.root", type=str, help="path to the root file of MC which you want to analyze", required=True)
parser.add_argument("-c", "--config", default="config.yml", type=str, help="path to the *.yml configuration file", required=True)
parser.add_argument("-s", "--suffix"  , default="" , type=str, help="suffix for output file name", required=False)
args = parser.parse_args();

#________________________________________________
def analyze(filename_data, filename_mc, config, suffix):
    rootfile_data = TFile.Open(filename_data, "READ");
    rootfile_mc   = TFile.Open(filename_mc  , "READ");
    rootdire_data = rootfile_data.Get("material-budget");
    rootdire_mc   = rootfile_mc  .Get("material-budget-mc");
    #rootdire_data.ls();
    #rootdire_mc  .ls();

    outname = "material_budget_correction_{0}_{1}TeV_{2}{3}.root".format(config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix);
    outfile = TFile(outname, "RECREATE");

    ssname = config["data"]["subsystems"][0]['name'];
    print("subsystem = ", ssname);
    list_pair_data = rootdire_data.Get("Pair").FindObject(ssname);
    list_pair_mc   = rootdire_mc  .Get("Pair").FindObject(ssname);
    #list_pair_data.ls();
    #list_pair_mc  .ls();

    list_event_data = rootdire_data.Get("Event").FindObject(ssname);
    list_event_mc   = rootdire_mc  .Get("Event").FindObject(ssname);
    #list_event_data.ls();
    #list_event_mc  .ls();
    h1mult_data = list_event_data.FindObject("hMultNTracksPV");
    h1mult_mc   = list_event_mc  .FindObject("hMultNTracksPV");
    nev_data = h1mult_data.GetEntries();
    nev_mc   = h1mult_mc  .GetEntries();
    print(nev_data, nev_mc);
    nch_data = h1mult_data.GetMean();
    nch_mc   = h1mult_mc  .GetMean();
    print(nch_data, nch_mc);

    list_v0_data = rootdire_data.Get("V0");
    list_v0_mc   = rootdire_mc  .Get("V0");
    #list_v0_data.ls();
    #list_v0_mc  .ls();

    arr_pt_probe = np.array(config["common"]["pt_probe"], dtype=np.float64);
    print(arr_pt_probe);
    tagnames = config["data"]["subsystems"][0]['tagnames']
    print(tagnames);
    probenames = config["data"]["subsystems"][0]['probenames']
    print(probenames);

    for it in range(0, len(tagnames)):
        tagname = tagnames[it];
        for ip in range(0, len(probenames)):
            probename_denominator = probenames[ip][0];
            probename_numerator = probenames[ip][1];
            cutname_denominator = tagname + "_" + probename_denominator;
            cutname_numerator   = tagname + "_" + probename_numerator;
            outlist = THashList();
            outlist.SetOwner(True);
            outlist.SetName(tagname + "__" + probename_denominator + "_" + probename_numerator);

            #cutname_mc_denominator = tagname + "_lowB_" + probename_denominator; #only for temporary solution
            #cutname_mc_numerator   = tagname + "_lowB_" + probename_numerator  ; #only for temporary solution

            list_pair_cut_denominator_data = list_pair_data.FindObject(cutname_denominator).FindObject("nocut");
            list_pair_cut_numerator_data   = list_pair_data.FindObject(cutname_numerator).FindObject("nocut");
            list_pair_cut_denominator_mc = list_pair_mc.FindObject(cutname_denominator).FindObject("nocut");
            list_pair_cut_numerator_mc   = list_pair_mc.FindObject(cutname_numerator).FindObject("nocut");

            hs_same_denominator_data = list_pair_cut_denominator_data.FindObject("hs_conv_point_same").Clone("hs_conv_point_same_{0}_data".format(probename_denominator));
            hs_mix_denominator_data = list_pair_cut_denominator_data.FindObject("hs_conv_point_mix").Clone("hs_conv_point_mix_{0}_data".format(probename_denominator));
            hs_same_numerator_data    = list_pair_cut_numerator_data.FindObject("hs_conv_point_same").Clone("hs_conv_point_same_{0}_data".format(probename_numerator));
            hs_mix_numerator_data    = list_pair_cut_numerator_data.FindObject("hs_conv_point_mix").Clone("hs_conv_point_mix_{0}_data".format(probename_numerator));
            hs_same_denominator_mc = list_pair_cut_denominator_mc.FindObject("hs_conv_point_same").Clone("hs_conv_point_same_{0}_mc".format(probename_denominator)); # MC truth
            hs_same_numerator_mc = list_pair_cut_numerator_mc.FindObject("hs_conv_point_same").Clone("hs_conv_point_same_{0}_mc".format(probename_numerator)); # MC truth

            #outlist.Add(hs_same_denominator_data);
            #outlist.Add(hs_same_numerator_data);
            #outlist.Add(hs_mix_denominator_data);
            #outlist.Add(hs_mix_numerator_data);
            #outlist.Add(hs_same_denominator_mc);
            #outlist.Add(hs_same_numerator_mc);

            h2same_denominator_data = hs_same_denominator_data.Projection(2,0);
            h2same_numerator_data   = hs_same_numerator_data  .Projection(2,0);
            h2mix_denominator_data  = hs_mix_denominator_data .Projection(2,0);
            h2mix_numerator_data    = hs_mix_numerator_data   .Projection(2,0);
            h2same_denominator_mc   = hs_same_denominator_mc  .Projection(2,0);
            h2same_numerator_mc     = hs_same_numerator_mc    .Projection(2,0);

            h2same_denominator_data.SetName(hs_same_denominator_data.GetName().replace("hs_conv_point_", "h2"));
            h2same_numerator_data  .SetName(hs_same_numerator_data  .GetName().replace("hs_conv_point_", "h2"));
            h2mix_denominator_data .SetName(hs_mix_denominator_data .GetName().replace("hs_conv_point_", "h2"));
            h2mix_numerator_data   .SetName(hs_mix_numerator_data   .GetName().replace("hs_conv_point_", "h2"));
            h2same_denominator_mc  .SetName(hs_same_denominator_mc  .GetName().replace("hs_conv_point_", "h2"));
            h2same_numerator_mc    .SetName(hs_same_numerator_mc    .GetName().replace("hs_conv_point_", "h2"));
            h2same_denominator_data.RebinX(2);
            h2same_numerator_data  .RebinX(2);
            h2mix_denominator_data .RebinX(2);
            h2mix_numerator_data   .RebinX(2);
            h2same_denominator_mc  .RebinX(2);
            h2same_numerator_mc    .RebinX(2);

            outlist.Add(h2same_denominator_data);
            outlist.Add(h2same_numerator_data);
            outlist.Add(h2mix_denominator_data);
            outlist.Add(h2mix_numerator_data);
            outlist.Add(h2same_denominator_mc);
            outlist.Add(h2same_numerator_mc);

            #2D analysis in meeg vs. pTg
            h1yield_denominator_data = TH1F("h1yield_{0}_data".format(probename_denominator), "raw yield;p_{T,#gamma} (GeV/c);#frac{1}{N_{ev}} #frac{dN}{dp_{T,#gamma}} (GeV/c)^{-1}", len(arr_pt_probe)-1, arr_pt_probe);
            h1mean_denominator_data  = TH1F("h1mean_{0}_data".format(probename_denominator), "mean;p_{T,#gamma} (GeV/c);peak position (GeV/c^{2})", len(arr_pt_probe)-1, arr_pt_probe);
            h1width_denominator_data = TH1F("h1width_{0}_data".format(probename_denominator), "width;p_{T,#gamma} (GeV/c);peak width (GeV/c^{2})", len(arr_pt_probe)-1, arr_pt_probe);
            h1yield_numerator_data = TH1F("h1yield_{0}_data".format(probename_numerator), "raw yield;p_{T,#gamma} (GeV/c);#frac{1}{N_{ev}} #frac{dN}{dp_{T,#gamma}} (GeV/c)^{-1}", len(arr_pt_probe)-1, arr_pt_probe);
            h1mean_numerator_data  = TH1F("h1mean_{0}_data".format(probename_numerator), "mean;p_{T,#gamma} (GeV/c);peak position (GeV/c^{2})", len(arr_pt_probe)-1, arr_pt_probe);
            h1width_numerator_data = TH1F("h1width_{0}_data".format(probename_numerator), "width;p_{T,#gamma} (GeV/c);peak width (GeV/c^{2})", len(arr_pt_probe)-1, arr_pt_probe);

            h1yield_denominator_mc = TH1F("h1yield_{0}_mc".format(probename_denominator), "raw yield;p_{T,#gamma} (GeV/c);#frac{1}{N_{ev}} #frac{dN}{dp_{T,#gamma}} (GeV/c)^{-1}", len(arr_pt_probe)-1, arr_pt_probe);
            h1mean_denominator_mc  = TH1F("h1mean_{0}_mc".format(probename_denominator), "mean;p_{T,#gamma} (GeV/c);peak position (GeV/c^{2})", len(arr_pt_probe)-1, arr_pt_probe);
            h1width_denominator_mc = TH1F("h1width_{0}_mc".format(probename_denominator), "width;p_{T,#gamma} (GeV/c);peak width (GeV/c^{2})", len(arr_pt_probe)-1, arr_pt_probe);
            h1yield_numerator_mc = TH1F("h1yield_{0}_mc".format(probename_numerator), "raw yield;p_{T,#gamma} (GeV/c);#frac{1}{N_{ev}} #frac{dN}{dp_{T,#gamma}} (GeV/c)^{-1}", len(arr_pt_probe)-1, arr_pt_probe);
            h1mean_numerator_mc  = TH1F("h1mean_{0}_mc".format(probename_numerator), "mean;p_{T,#gamma} (GeV/c);peak position (GeV/c^{2})", len(arr_pt_probe)-1, arr_pt_probe);
            h1width_numerator_mc = TH1F("h1width_{0}_mc".format(probename_numerator), "width;p_{T,#gamma} (GeV/c);peak width (GeV/c^{2})", len(arr_pt_probe)-1, arr_pt_probe);

            h1yield_denominator_data.Sumw2();
            h1mean_denominator_data .Sumw2();
            h1width_denominator_data.Sumw2();
            h1yield_numerator_data  .Sumw2();
            h1mean_numerator_data   .Sumw2();
            h1width_numerator_data  .Sumw2();
            h1yield_denominator_mc.Sumw2();
            h1mean_denominator_mc .Sumw2();
            h1width_denominator_mc.Sumw2();
            h1yield_numerator_mc  .Sumw2();
            h1mean_numerator_mc   .Sumw2();
            h1width_numerator_mc  .Sumw2();

            for ipt in range(0, len(arr_pt_probe)-1):
                pt1 = arr_pt_probe[ipt];
                pt2 = arr_pt_probe[ipt+1];

                #for data
                h1same_denominator_data = slice_histogram(h2same_denominator_data, pt1, pt2, "X", False);
                h1same_denominator_data.SetName("{0}_pt{1:d}".format(h2same_denominator_data.GetName().replace("h2", "h1"),ipt));
                h1mix_denominator_data = slice_histogram(h2mix_denominator_data, pt1, pt2, "X", False);
                h1mix_denominator_data.SetName("{0}_pt{1:d}".format(h2mix_denominator_data.GetName().replace("h2", "h1"),ipt));
                nmf = NMFitter(h1same_denominator_data, h1mix_denominator_data, "cb", "pol1");
                nmf.set_parameters(0.135, 0.005, 0.6, 12, True);
                [fitresult, h1sig, h1bkg, h1ratio, f1sig, f1bkg, f1total] = nmf.fit("SME", "", 0.04, 0.24);
                h1ratio.SetName("{0}".format(h1same_denominator_data.GetName().replace("same", "ratio")));
                h1sig.SetName("{0}".format(h1same_denominator_data.GetName().replace("same", "sig")));
                h1bkg.SetName("{0}".format(h1same_denominator_data.GetName().replace("same", "bkg")));
                f1sig.SetName("{0}".format(h1same_denominator_data.GetName().replace("h1same", "f1sig")));
                f1bkg.SetName("{0}".format(h1same_denominator_data.GetName().replace("h1same", "f1bkg")));
                f1total.SetName("{0}".format(h1same_denominator_data.GetName().replace("h1same", "f1total")));
                outlist.Add(h1same_denominator_data);
                outlist.Add(h1mix_denominator_data);
                outlist.Add(h1ratio);
                outlist.Add(h1sig);
                outlist.Add(h1bkg);
                outlist.Add(f1sig);
                outlist.Add(f1bkg);
                outlist.Add(f1total);

                mean = f1sig.GetParameter(1);
                mean_err = f1sig.GetParError(1);
                sigma = f1sig.GetParameter(2);
                sigma_err = f1sig.GetParError(2);
                bin1 = h1sig.FindBin(mean - 3 * sigma);
                bin2 = h1sig.FindBin(mean + 3 * sigma);
                n_err = ctypes.c_double(0);
                n = h1sig.IntegralAndError(bin1, bin2, n_err, "");
                h1yield_denominator_data.SetBinContent(ipt+1, n);
                h1yield_denominator_data.SetBinError(ipt+1, n_err);
                h1mean_denominator_data.SetBinContent(ipt+1, mean);
                h1mean_denominator_data.SetBinError(ipt+1, mean_err);
                h1width_denominator_data.SetBinContent(ipt+1, sigma);
                h1width_denominator_data.SetBinError(ipt+1, sigma_err);

                h1same_numerator_data = slice_histogram(h2same_numerator_data, pt1, pt2, "X", False);
                h1same_numerator_data.SetName("{0}_pt{1:d}".format(h2same_numerator_data.GetName().replace("h2", "h1"),ipt));
                h1mix_numerator_data = slice_histogram(h2mix_numerator_data, pt1, pt2, "X", False);
                h1mix_numerator_data.SetName("{0}_pt{1:d}".format(h2mix_numerator_data.GetName().replace("h2", "h1"),ipt));
                nmf = NMFitter(h1same_numerator_data, h1mix_numerator_data, "cb", "pol1");
                nmf.set_parameters(0.135, 0.005, 0.6, 12, True);
                [fitresult, h1sig, h1bkg, h1ratio, f1sig, f1bkg, f1total] = nmf.fit("SME", "", 0.04, 0.24);
                h1ratio.SetName("{0}".format(h1same_numerator_data.GetName().replace("same", "ratio")));
                h1sig.SetName("{0}".format(h1same_numerator_data.GetName().replace("same", "sig")));
                h1bkg.SetName("{0}".format(h1same_numerator_data.GetName().replace("same", "bkg")));
                f1sig.SetName("{0}".format(h1same_numerator_data.GetName().replace("h1same", "f1sig")));
                f1bkg.SetName("{0}".format(h1same_numerator_data.GetName().replace("h1same", "f1bkg")));
                f1total.SetName("{0}".format(h1same_numerator_data.GetName().replace("h1same", "f1total")));
                outlist.Add(h1same_numerator_data);
                outlist.Add(h1mix_numerator_data);
                outlist.Add(h1ratio);
                outlist.Add(h1sig);
                outlist.Add(h1bkg);
                outlist.Add(f1sig);
                outlist.Add(f1bkg);
                outlist.Add(f1total);
                mean = f1sig.GetParameter(1);
                mean_err = f1sig.GetParError(1);
                sigma = f1sig.GetParameter(2);
                sigma_err = f1sig.GetParError(2);
                bin1 = h1sig.FindBin(mean - 3 * sigma);
                bin2 = h1sig.FindBin(mean + 3 * sigma);
                n_err = ctypes.c_double(0);
                n = h1sig.IntegralAndError(bin1, bin2, n_err, "");
                h1yield_numerator_data.SetBinContent(ipt+1, n);
                h1yield_numerator_data.SetBinError(ipt+1, n_err);
                h1mean_numerator_data.SetBinContent(ipt+1, mean);
                h1mean_numerator_data.SetBinError(ipt+1, mean_err);
                h1width_numerator_data.SetBinContent(ipt+1, sigma);
                h1width_numerator_data.SetBinError(ipt+1, sigma_err);

                #for MC truth
                h1same_denominator_mc = slice_histogram(h2same_denominator_mc, pt1, pt2, "X", False);
                h1same_denominator_mc.SetName("{0}_pt{1:d}".format(h2same_denominator_mc.GetName().replace("h2", "h1"),ipt));
                nmf = NMFitter(h1same_denominator_mc, None, "cb", "none");
                nmf.set_parameters(0.135, 0.005, 0.6, 12, True);
                [fitresult, h1sig, _, _, f1sig, _, _] = nmf.fit("SME", "", 0.04, 0.24);
                h1sig.SetName("{0}".format(h1same_denominator_mc.GetName().replace("same", "sig")));
                f1sig.SetName("{0}".format(h1same_denominator_mc.GetName().replace("h1same", "f1sig")));
                outlist.Add(h1same_denominator_mc);
                outlist.Add(h1sig);
                outlist.Add(f1sig);
                mean = f1sig.GetParameter(1);
                mean_err = f1sig.GetParError(1);
                sigma = f1sig.GetParameter(2);
                sigma_err = f1sig.GetParError(2);
                bin1 = h1sig.FindBin(mean - 3 * sigma);
                bin2 = h1sig.FindBin(mean + 3 * sigma);
                n_err = ctypes.c_double(0);
                n = h1sig.IntegralAndError(bin1, bin2, n_err, "");
                h1yield_denominator_mc.SetBinContent(ipt+1, n);
                h1yield_denominator_mc.SetBinError(ipt+1, n_err);
                h1mean_denominator_mc.SetBinContent(ipt+1, mean);
                h1mean_denominator_mc.SetBinError(ipt+1, mean_err);
                h1width_denominator_mc.SetBinContent(ipt+1, sigma);
                h1width_denominator_mc.SetBinError(ipt+1, sigma_err);

                h1same_numerator_mc = slice_histogram(h2same_numerator_mc, pt1, pt2, "X", False);
                h1same_numerator_mc.SetName("{0}_pt{1:d}".format(h2same_numerator_mc.GetName().replace("h2", "h1"),ipt));
                nmf = NMFitter(h1same_numerator_mc, None, "cb", "none");
                nmf.set_parameters(0.135, 0.005, 0.6, 12, True);
                [fitresult, h1sig, _, _, f1sig, _, _] = nmf.fit("SME", "", 0.04, 0.24);
                h1sig.SetName("{0}".format(h1same_numerator_mc.GetName().replace("same", "sig")));
                f1sig.SetName("{0}".format(h1same_numerator_mc.GetName().replace("h1same", "f1sig")));
                outlist.Add(h1same_numerator_mc);
                outlist.Add(h1sig);
                outlist.Add(f1sig);
                mean = f1sig.GetParameter(1);
                mean_err = f1sig.GetParError(1);
                sigma = f1sig.GetParameter(2);
                sigma_err = f1sig.GetParError(2);
                bin1 = h1sig.FindBin(mean - 3 * sigma);
                bin2 = h1sig.FindBin(mean + 3 * sigma);
                n_err = ctypes.c_double(0);
                n = h1sig.IntegralAndError(bin1, bin2, n_err, "");
                h1yield_numerator_mc.SetBinContent(ipt+1, n);
                h1yield_numerator_mc.SetBinError(ipt+1, n_err);
                h1mean_numerator_mc.SetBinContent(ipt+1, mean);
                h1mean_numerator_mc.SetBinError(ipt+1, mean_err);
                h1width_numerator_mc.SetBinContent(ipt+1, sigma);
                h1width_numerator_mc.SetBinError(ipt+1, sigma_err);

            h1yield_denominator_data.Scale(1./nev_data);
            h1yield_numerator_data.Scale(1./nev_data);
            h1yield_denominator_mc.Scale(1./nev_mc);
            h1yield_numerator_mc.Scale(1./nev_mc);
            h1yield_denominator_data.Scale(1., "width");
            h1yield_numerator_data  .Scale(1., "width");
            h1yield_denominator_mc  .Scale(1., "width");
            h1yield_numerator_mc    .Scale(1., "width");

            outlist.Add(h1yield_denominator_data);
            outlist.Add(h1mean_denominator_data);
            outlist.Add(h1width_denominator_data);
            outlist.Add(h1yield_numerator_data);
            outlist.Add(h1mean_numerator_data);
            outlist.Add(h1width_numerator_data);
            outlist.Add(h1yield_denominator_mc);
            outlist.Add(h1mean_denominator_mc);
            outlist.Add(h1width_denominator_mc);
            outlist.Add(h1yield_numerator_mc);
            outlist.Add(h1mean_numerator_mc);
            outlist.Add(h1width_numerator_mc);

            h1frac_data = h1yield_numerator_data.Clone("h1frac_data");
            h1frac_data.Reset();
            h1frac_data.SetTitle("N_{#gamma #rightarrow ee}^{total}/N_{#gamma #rightarrow ee}^{W wires IB} in Data");
            h1frac_data.SetYTitle("N_{#gamma #rightarrow ee}^{total}/N_{#gamma #rightarrow ee}^{W wires IB}");
            h1frac_data.Divide(h1yield_numerator_data, h1yield_denominator_data, 1., 1., "B");
            outlist.Add(h1frac_data);

            h1frac_mc = h1yield_numerator_mc.Clone("h1frac_mc");
            h1frac_mc.Reset();
            h1frac_mc.SetTitle("N_{#gamma #rightarrow ee}^{total}/N_{#gamma #rightarrow ee}^{W wires IB} in M.C.");
            h1frac_mc.SetYTitle("N_{#gamma #rightarrow ee}^{total}/N_{#gamma #rightarrow ee}^{W wires IB}");
            h1frac_mc.Divide(h1yield_numerator_mc, h1yield_denominator_mc, 1., 1., "B");
            outlist.Add(h1frac_mc);

            h1dr = h1frac_data.Clone("h1double_ratio");
            h1dr.Reset();
            h1dr.SetTitle("double ratio");
            h1dr.SetYTitle("#frac{N_{#gamma #rightarrow ee}^{total}/N_{#gamma #rightarrow ee}^{W wires IB} in Data}{N_{#gamma #rightarrow ee}^{total}/N_{#gamma #rightarrow ee}^{W wires IB} in M.C.}");
            h1dr.Divide(h1frac_data, h1frac_mc, 1., 1., "G");
            outlist.Add(h1dr);

            #compare the number of photon conversions on W wires between data and MC.
            list_v0_probe_data = list_v0_data.FindObject(probename_denominator);
            list_v0_probe_mc   = list_v0_mc  .FindObject(probename_denominator);
            hs_data = list_v0_probe_data.FindObject("hs_conv_point").Clone("hs_conv_point_data");
            hs_mc   = list_v0_probe_mc  .FindObject("hs_conv_point").Clone("hs_conv_point_mc");

            h1pt_data_org = hs_data.Projection(0);
            h1pt_data = rebin_histogram(h1pt_data_org, arr_pt_probe, False, False);
            h1pt_data.SetName("h1pt_{0}_data".format(probename_denominator));
            h1pt_data.Scale(1., "width");
            h1pt_data.Scale(1./nev_data/nch_data);
            outlist.Add(h1pt_data);

            h1pt_mc_org = hs_mc.Projection(0);
            h1pt_mc = rebin_histogram(h1pt_mc_org, arr_pt_probe, False, False);
            h1pt_mc.SetName("h1pt_{0}_mc".format(probename_denominator));
            h1pt_mc.Scale(1., "width");
            h1pt_mc.Scale(1./nev_mc/nch_mc);
            outlist.Add(h1pt_mc);

            h1ratio_wire = h1pt_data.Clone("h1ratio_{0}".format(probename_denominator));
            h1ratio_wire.SetTitle("ratio of photon conversions on wires");
            h1ratio_wire.SetYTitle("#frac{N_{#gamma #rightarrow ee}^{W wires IB} in Data}{N_{#gamma #rightarrow ee}^{W wires IB} in M.C.}");
            h1ratio_wire.Reset();
            h1ratio_wire.Divide(h1pt_data, h1pt_mc, 1., 1., "G");
            outlist.Add(h1ratio_wire);

            h1corr = h1dr.Clone("h1corr");
            h1corr.Reset();
            h1corr.SetTitle("material budget correction");
            h1corr.SetYTitle("correction factor");
            h1corr.Multiply(h1dr, h1ratio_wire, 1., 1., "B");
            outlist.Add(h1corr); #multiply this factor to MC in your analysis.

            outfile.WriteTObject(outlist);
            outlist.Clear();

    rootfile_data.Close();
    rootfile_mc  .Close();
    outfile.Close();
#________________________________________________
if __name__ == "__main__":
    filename_data = args.input_data;
    filename_mc   = args.input_mc;

    config = "";
    with open(args.config, "r", encoding="utf-8") as config_yml:
        config = yaml.safe_load(config_yml)
    analyze(filename_data, filename_mc, config, args.suffix);

#________________________________________________
