import os, sys, shutil
import numpy as np
import math
import ctypes
import ROOT
from ROOT import TFile, TDirectory, THashList, TF1, TH1F
from histo_manager import slice_histogram, rebin_histogram, get_bkg_subtracted, get_ratio
#import file_manager

class SinglePhotonAnalyzer:
    def __init__(self):
        print("default constructor is called");
    def __init__(self, particle, filename, dirname, ismc):
        print("target particle = {0} , filename = {1} , dirname = {2}".format(particle, filename, dirname));
        self.particle = particle;
        self.rootfile = TFile.Open(filename, "READ");
        self.rootdir = self.rootfile.Get(dirname);
        self.list_ev = self.rootdir.Get("Event");
        #self.list_v0 = self.rootdir.Get("Photon");
        self.list_v0 = self.rootdir.Get("V0");
        self.arr_pt = np.array([0,1,2,3,4,5], dtype=float);
        self.ismc = ismc;
        self.list_gen = None;
        if self.ismc:
            self.list_gen = self.rootdir.Get("Generated");

        self.xtitle = "#it{p}_{T,#gamma} (GeV/#it{c})";
        self.ytitle = "#it{p}_{T,#gamma} (GeV/#it{c})";

        #self.list_v0.Print();

    def __del__(self):
        if self.rootfile.IsOpen():
            print("close input root file.");
            self.rootfile.Close();

    def set_arr_pt(self, arr_pt):
        print("pT array = ", arr_pt);
        self.arr_pt = arr_pt;

    def set_subsystem(self, ssname):
        self.ssname = ssname;
        #self.list_ev_ss = self.list_ev.FindObject(ssname);
        self.list_ev_ss = self.list_ev;
        #self.list_v0_ss = self.list_v0.FindObject(ssname);
        self.list_v0_ss = self.list_v0;

    def set_cutname(self, cutname):
        self.cutname = cutname;
        if self.list_ev_ss is None or self.list_v0_ss is None:
            print("Please define subsystem name first!");
            return None;
        self.list_v0_ss_cut = self.list_v0_ss.FindObject(cutname);

    def set_xtitle(self, title):
        self.xtitle = title;

    def set_ytitle(self, title):
        self.ytitle = title;

    def analyze_ptspectrum(self): #this is main function
        print(sys._getframe().f_code.co_name);
        outlist = THashList();
        #outlist.SetOwner(True);
        outlist.SetName("outlist");

        h1ev = self.list_ev_ss.FindObject("hCollisionCounter").Clone("h1ev");
        h1pt_org = self.list_v0_ss_cut.FindObject("hPt").Clone("hPt");
        h1pt_org.Sumw2();
        h1pt_org.SetDirectory(0);
        nev = h1ev.GetBinContent(4);
        print("nev = {0:e}".format(nev));
   
        h1yield = rebin_histogram(h1pt_org, self.arr_pt, True, False);
        h1yield.SetName("h1yield");
        h1yield.SetTitle("raw yield of photon candidates");
        h1yield.SetXTitle(self.xtitle);
        h1yield.SetYTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T,#gamma}} (GeV/#it{c})^{-1}");
        h1yield.Scale(1/nev); 
        outlist.Add(h1yield);
        return outlist;

    def analyze_ptspectrum_efficiency(self): #this is main function
        print(sys._getframe().f_code.co_name);
        outlist = THashList();
        #outlist.SetOwner(True);
        outlist.SetName("outlist");

        self.list_v0_ss_cut.Print();
        h1ev   = self.list_ev_ss.FindObject("hCollisionCounter").Clone("h1ev");
        h1pt_primary_org = self.list_v0_ss_cut.FindObject("hPt_Photon_Primary").Clone("hPt_Photon_Primary");
        h1pt_wd_org      = self.list_v0_ss_cut.FindObject("hPt_Photon_FromWD").Clone("hPt_Photon_FromWD");
        h1pt_primary_org.Sumw2();
        h1pt_primary_org.SetDirectory(0);
        h1pt_wd_org.Sumw2();
        h1pt_wd_org.SetDirectory(0);
        nev = h1ev.GetBinContent(4);
        print("nev = {0:e}".format(nev));
        outlist.Add(h1ev);
   
        h1dndpt_primary = rebin_histogram(h1pt_primary_org, self.arr_pt, True, False);
        h1dndpt_primary.SetName("h1dndpt_primary");
        h1dndpt_primary.Scale(1/nev); 
        h1dndpt_primary.SetXTitle(self.xtitle);
        h1dndpt_primary.SetYTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T,#gamma}} (GeV/#it{c})^{-1}");
        outlist.Add(h1dndpt_primary);

        h1dndpt_wd = rebin_histogram(h1pt_wd_org, self.arr_pt, True, False);
        h1dndpt_wd.SetName("h1dndpt_wd");
        h1dndpt_wd.Scale(1/nev); 
        h1dndpt_wd.SetXTitle(self.xtitle);
        h1dndpt_wd.SetYTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T,#gamma}} (GeV/#it{c})^{-1}");
        outlist.Add(h1dndpt_wd);

        h1dndpt_all = h1dndpt_primary.Clone("h1dndpt_all");
        h1dndpt_all.Add(h1dndpt_wd, 1.);
        outlist.Add(h1dndpt_all);
        h1fd = h1dndpt_wd.Clone("h1fd");
        h1fd.SetTitle("feed down");
        h1fd.SetYTitle("#frac{#gamma from K^{0}_{S} or #Lambda}{all #gamma}");
        h1fd.Sumw2();
        h1fd.Reset();
        h1fd.Divide(h1dndpt_wd, h1dndpt_all, 1., 1., "B");

        #f1fd = TF1("f1fd","[0] + [1]*exp(-x/[2])",0,10);
        f1fd = TF1("f1fd","[0] + [1]/x",0.1,10);
        f1fd.SetNpx(1000);
        f1fd.SetParameters(0.02, 0.02);
        h1fd.Fit(f1fd,"SME","",0,10);
        outlist.Add(h1fd);
        outlist.Add(f1fd);

        #purity
        h1pt_org = self.list_v0_ss_cut.FindObject("hPt").Clone("hPt");
        h1pt_org.Sumw2();
        h1pt_org.SetDirectory(0);
        h1dndpt_candidate = rebin_histogram(h1pt_org, self.arr_pt, True, False);
        h1dndpt_candidate.SetName("h1dndpt_candidate");
        h1dndpt_candidate.SetName("photon candidate");
        h1dndpt_candidate.Scale(1/nev); 
        h1dndpt_candidate.SetXTitle(self.xtitle);
        h1dndpt_candidate.SetYTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T,#gamma}} (GeV/#it{c})^{-1}");
        outlist.Add(h1dndpt_candidate);

        h1purity = h1dndpt_primary.Clone("h1purity");
        h1purity.Reset();
        h1purity.Sumw2();
        h1purity.Divide(h1dndpt_primary, h1dndpt_candidate, 1., 1., "B");
        h1purity.SetTitle("photon purity");
        h1purity.SetYTitle("purity");
        outlist.Add(h1purity);

 
        #Next, generated information
        h1pt_org  = self.list_gen.FindObject("hPt_Photon").Clone("hPt_Photon");
        h1y_org   = self.list_gen.FindObject("hY_Photon").Clone("hY_Photon");
        h1phi_org = self.list_gen.FindObject("hPhi_Photon").Clone("hPhi_Photon");
        h1pt_org .Sumw2();
        h1y_org  .Sumw2();
        h1phi_org.Sumw2();
        h1pt_org.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1pt_org.SetYTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}");
        outlist.Add(h1pt_org);
        outlist.Add(h1y_org);
        outlist.Add(h1phi_org);
        h1dndpt_gen = rebin_histogram(h1pt_org, self.arr_pt, True, False);
        h1dndpt_gen.Scale(1/nev);
        h1dndpt_gen.SetName("h1dndpt_gen");
        h1dndpt_gen.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1dndpt_gen.SetYTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}");
        outlist.Add(h1dndpt_gen);

        h1eff = h1dndpt_primary.Clone("h1eff");
        h1eff.Sumw2();
        h1eff.SetTitle("efficiency");
        h1eff.SetYTitle("acc. #times rec. efficiency");
        h1eff.Reset();
        h1eff.Divide(h1dndpt_primary, h1dndpt_gen, 1., 1., "B");
        outlist.Add(h1eff);
        return outlist;

#___________________________________________________________________
