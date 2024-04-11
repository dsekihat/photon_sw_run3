import os, sys, shutil
import math
import argparse
import numpy as np
import ctypes
import yaml
import ROOT
ROOT.gROOT.SetBatch(True);
from ROOT import TFile, THashList, TF1
from single_photon_analyzer import SinglePhotonAnalyzer

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
    nsys = len(config[args.type]['subsystems']);
    #print(nsys); 
    particle = config["common"]["particle"];

    if config["common"]["do_ptspectrum"] == True:
        outname = "{0}_{1}_ptspectrum_{2}_{3}TeV_{4}{5}.root".format(particle, args.type, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix);
        print("output file name = ",outname);
        outfile = TFile(outname,"RECREATE");
        if ismc:
            #ana_pi0 = SinglePhotonAnalyzer(particle, filename, "single-photon-mc", ismc);
            ana_pi0 = SinglePhotonAnalyzer(particle, filename, "pcm-qc-mc", ismc);
            ana_pi0.set_arr_pt(arr_pt);
            for isys in range(0,nsys):
                ssname = config[args.type]['subsystems'][isys]['name']; #subsystem name
                ana_pi0.set_subsystem(ssname);
                ana_pi0.set_xtitle("#it{p}_{T,#gamma} (GeV/#it{c})");
                ana_pi0.set_ytitle("#it{p}_{T,#gamma} (GeV/#it{c})");
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
                    outlist_cut = ana_pi0.analyze_ptspectrum_efficiency();
                    outlist_cut.SetName(cutname);
                    outlist_cut.SetOwner(True);
                    outlist_ss.Add(outlist_cut);
                outfile.WriteTObject(outlist_ss);
                outlist_ss.Clear();
            del ana_pi0;
        else:
            ana_pi0 = SinglePhotonAnalyzer(particle, filename, "single-photon", ismc);
            ana_pi0.set_arr_pt(arr_pt);
            for isys in range(0,nsys):
                ssname = config[args.type]['subsystems'][isys]['name']; #subsystem name
                ana_pi0.set_subsystem(ssname);
                ana_pi0.set_xtitle("#it{p}_{T,#gamma} (GeV/#it{c})");
                ana_pi0.set_ytitle("#it{p}_{T,#gamma} (GeV/#it{c})");
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
                    outlist_cut = ana_pi0.analyze_ptspectrum();
                    outlist_cut.SetName(cutname);
                    outlist_cut.SetOwner(True);
                    outlist_ss.Add(outlist_cut);
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
