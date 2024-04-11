import os, sys, shutil
import math
import argparse
import numpy as np
import ctypes
import yaml
import ROOT
ROOT.gROOT.SetBatch(True);
from ROOT import TFile, THashList, TF1
#from analyze_pair import analyze_ptspectrum
from pair_analyzer import PairAnalyzer

parser = argparse.ArgumentParser('Example program');
parser.add_argument("-i", "--input" , default="AnalysisResults.root", type=str, help="path to the root file you want to analyze", required=True)
parser.add_argument("-c", "--config", default="config.yml", type=str, help="path to the *.yml configuration file", required=True)
parser.add_argument("-t", "--type"  , default="data" , type=str, help="run type [data or mc]", required=True)
parser.add_argument("-s", "--suffix"  , default="" , type=str, help="suffix for output file name", required=False)
args = parser.parse_args();

filename = args.input;
with open(args.config, "r", encoding="utf-8") as config_yml:
    config = yaml.safe_load(config_yml)
#_________________________________________________________________________________________
def run(filename, config, ismc, suffix=""):
    print(sys._getframe().f_code.co_name);
    arr_pt = np.array(config["common"]["pt_bin"],dtype=float);
    print("pT binning = ",arr_pt);
    print("ismc = ",ismc);
    print("input = ",filename);
    rootfile = TFile.Open(filename,"READ");
    particle = config["common"]["particle"];
    #print(config);

    list_fit_sig = config["common"]["fit_sig"];
    list_fit_bkg = config["common"]["fit_bkg"];

    list_fit_min = config["common"]["fit_min"];
    list_fit_max = config["common"]["fit_max"];
    if len(list_fit_min) != len(list_fit_max):
        return;

    list_integral_min = config["common"]["integral_sigma_min"];
    list_integral_max = config["common"]["integral_sigma_max"];
    if len(list_integral_min) != len(list_integral_max):
        return;

    nsys = len(config[args.type]['subsystems']);
    #print(nsys); 

    particle = config["common"]["particle"];

    if config["common"]["do_ptspectrum"] == True:
        outname = "{0}_{1}_ptspectrum_{2}_{3}TeV_{4}{5}.root".format(particle, args.type, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix);
        print("output file name = ",outname);
        outfile = TFile(outname,"RECREATE");
        if ismc:
            ana_pi0 = PairAnalyzer(particle, filename, "pi0eta-to-gammagamma-mc", ismc);
            ana_pi0.set_arr_pt(arr_pt);
            for isys in range(0,nsys):
                ssname = config[args.type]['subsystems'][isys]['name']; #subsystem name
                ana_pi0.set_subsystem(ssname);
                ana_pi0.set_xtitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
                ana_pi0.set_ytitle("#it{p}_{T,#gamma#gamma} (GeV/#it{c}^{2})");
                print("analyze subsystem", ssname);
                cutnames = config[args.type]["subsystems"][isys]['cutnames']
                print("cutnames", cutnames); 
                nc = len(cutnames);
                outlist_ss = THashList();
                outlist_ss.SetName(ssname);
                outlist_ss.SetOwner(True);
                for ic in range(0,nc):
                    cutname = cutnames[ic];
                    ana_pi0.set_cutname(cutname);
                    outlist_cut = THashList();
                    outlist_cut.SetName(cutname);
                    outlist_cut.SetOwner(True);
                    outlist_ss.Add(outlist_cut);
                    for isig in list_fit_sig:
                        for ibkg in list_fit_bkg:
                            ana_pi0.set_fit_function(isig, ibkg);
                            ana_pi0.fix_cb_n(12.0);
                            ana_pi0.fix_cb_alpha(0.6);
                            outlist_func = THashList();
                            outlist_func.SetName(isig + "_" + ibkg);
                            outlist_func.SetOwner(True);
                            outlist_cut.Add(outlist_func);
                            for ir in range(0, len(list_fit_min)):
                                fit_min = list_fit_min[ir];
                                fit_max = list_fit_max[ir];
                                ana_pi0.set_fit_range(fit_min, fit_max);
                                outlist_fit_range = THashList();
                                outlist_fit_range.SetName("fit_{0:3.2f}_{1:3.2f}_GeVc2".format(fit_min, fit_max));
                                outlist_fit_range.SetOwner(True);
                                outlist_func.Add(outlist_fit_range);
                                for iint in range(0, len(list_integral_min)):
                                    integral_sigma_min = list_integral_min[iint];
                                    integral_sigma_max = list_integral_max[iint];
                                    ana_pi0.set_integral_range(integral_sigma_min, integral_sigma_max);
                                    outlist_int_range = ana_pi0.analyze_ptspectrum_efficiency();
                                    outlist_int_range.SetName("integral_{0:2.1f}_{1:2.1f}_sigma".format(integral_sigma_min, integral_sigma_max));
                                    #outlist_int_range.SetOwner(True);
                                    outlist_fit_range.Add(outlist_int_range);
                outfile.WriteTObject(outlist_ss);
                outlist_ss.Clear();
            del ana_pi0;
        else:
            ana_pi0 = PairAnalyzer(particle, filename, "pi0eta-to-gammagamma", ismc);
            ana_pi0.set_arr_pt(arr_pt);
            for isys in range(0,nsys):
                ssname = config[args.type]['subsystems'][isys]['name']; #subsystem name
                ana_pi0.set_subsystem(ssname);
                ana_pi0.set_xtitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
                ana_pi0.set_ytitle("#it{p}_{T,#gamma#gamma} (GeV/#it{c}^{2})");
                print("analyze subsystem", ssname);
                cutnames = config[args.type]["subsystems"][isys]['cutnames']
                print("cutnames", cutnames); 
                nc = len(cutnames);
                outlist_ss = THashList();
                outlist_ss.SetName(ssname);
                outlist_ss.SetOwner(True);
                for ic in range(0,nc):
                    cutname = cutnames[ic];
                    ana_pi0.set_cutname(cutname);
                    outlist_cut = THashList();
                    outlist_cut.SetName(cutname);
                    outlist_cut.SetOwner(True);
                    outlist_ss.Add(outlist_cut);
                    for isig in list_fit_sig:
                        for ibkg in list_fit_bkg:
                            ana_pi0.set_fit_function(isig, ibkg);
                            ana_pi0.fix_cb_n(12.0);
                            ana_pi0.fix_cb_alpha(0.6);
                            outlist_func = THashList();
                            outlist_func.SetName(isig + "_" + ibkg);
                            outlist_func.SetOwner(True);
                            outlist_cut.Add(outlist_func);
                            for ir in range(0, len(list_fit_min)):
                                fit_min = list_fit_min[ir];
                                fit_max = list_fit_max[ir];
                                ana_pi0.set_fit_range(fit_min, fit_max);
                                outlist_fit_range = THashList();
                                outlist_fit_range.SetName("fit_{0:3.2f}_{1:3.2f}_GeVc2".format(fit_min, fit_max));
                                outlist_fit_range.SetOwner(True);
                                outlist_func.Add(outlist_fit_range);
                                for iint in range(0, len(list_integral_min)):
                                    integral_sigma_min = list_integral_min[iint];
                                    integral_sigma_max = list_integral_max[iint];
                                    ana_pi0.set_integral_range(integral_sigma_min, integral_sigma_max);
                                    outlist_int_range = ana_pi0.analyze_ptspectrum();
                                    outlist_int_range.SetName("integral_{0:2.1f}_{1:2.1f}_sigma".format(integral_sigma_min, integral_sigma_max));
                                    #outlist_int_range.SetOwner(True);
                                    outlist_fit_range.Add(outlist_int_range);
                outfile.WriteTObject(outlist_ss);
                outlist_ss.Clear();
            del ana_pi0;
        outfile.Close();
    else:
        print("please check what to do in",args.config);

    rootfile.Close();
#_________________________________________________________________________________________
ismc = False;
if args.type == "data":
    ismc = False;
elif args.type == "mc":
    ismc = True;
else:
    print("unknown type.sys.exit()");
    sys.exit();
run(filename,config,ismc,args.suffix);
