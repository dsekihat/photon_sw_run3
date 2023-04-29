import os, sys, shutil
import numpy as np
import math
import ctypes
import ROOT
from ROOT import TFile, TDirectory, THashList, TF1, TH1F
from histo_manager import slice_histogram, rebin_histogram, get_bkg_subtracted, get_ratio
#import file_manager

class PairAnalyzer:
    def __init__(self):
        print("default constructor is called");
    def __init__(self, meson, filename, dirname):
        print("target meson = {0} , filename = {1} , dirname = {2}".format(meson, filename, dirname));
        self.meson = meson;
        self.rootfile = TFile.Open(filename, "READ");
        self.rootdir = self.rootfile.Get(dirname);
        self.list_ev = self.rootdir.Get("Event");
        self.list_pair = self.rootdir.Get("Pair");
        self.arr_pt = np.array([0,1,2,3,4,5], dtype=float);

        self.f1sig = TF1("f1sig","crystalball(0)", 0, 1);
        self.f1bkg = TF1("f1bkg","pol1(0)",0,1);
        self.f1total = TF1("f1total","f1sig + f1bkg", 0, 1);
        self.f1total.SetNpx(1000);
        self.fit_min = 0.04;
        self.fit_max = 0.24;
        self.integral_min_sigma = -3.0;
        self.integral_max_sigma = +3.0;
        self.xtitle = "#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})";
        self.ytitle = "#it{p}_{T,#gamma#gamma} (GeV/#it{c}^{2})";

    def __del__(self):
        if self.rootfile.IsOpen():
            print("close input root file.");
            self.rootfile.Close();

    def set_arr_pt(self, arr_pt):
        print("pT array = ", arr_pt);
        self.arr_pt = arr_pt;

    def set_subsystem(self, ssname):
        self.ssname = ssname;
        self.list_ev_ss   = self.list_ev.FindObject(ssname);
        self.list_pair_ss = self.list_pair.FindObject(ssname);

    def set_cutname(self, cutname):
        self.cutname = cutname;
        if self.list_ev_ss is None or self.list_pair_ss is None:
            print("Please define subsystem name first!");
            return None;
        #self.list_ev_ss_cut   = self.list_ev_ss.FindObject(cutname);
        self.list_pair_ss_cut = self.list_pair_ss.FindObject(cutname);

    def set_fit_range(self, fit_min, fit_max):
        self.fit_min = fit_min;
        self.fit_max = fit_max;

    def set_integral_range(self, integral_min, integral_max):
        self.integral_min_sigma = integral_min;
        self.integral_max_sigma = integral_max;

    def set_fit_function(self, sig, bkg):
        if sig == "cb":
            self.f1sig = TF1("f1sig","crystalball(0)", 0,1);
            if bkg == "pol1":
                self.f1bkg   = TF1("f1bkg","pol1(0)", 0,1);
                self.f1total = TF1("f1total","crystalball(0) + pol1(5)", 0,1);
            elif bkg == "pol2":
                self.f1bkg   = TF1("f1bkg","pol2(0)", 0,1);
                self.f1total = TF1("f1total","crystalball(0) + pol2(5)", 0,1);

        self.f1sig.SetNpx(1000);
        self.f1bkg.SetNpx(1000);
        self.f1total.SetNpx(1000);
        print("initially, ", self.f1total.GetExpFormula(""));

    def set_xtitle(self, title):
        self.xtitle = title;

    def set_ytitle(self, title):
        self.ytitle = title;

    def analyze_ptspectrum(self): #this is main function
        #print(sys._getframe().f_code.co_name);
        outlist = THashList();
        #outlist.SetOwner(True);
        outlist.SetName("outlist");

        h1ev   = self.list_ev_ss.FindObject("hCollisionCounter").Clone("h1ev");
        h2same = self.list_pair_ss_cut.FindObject("hMggPt_Same").Clone("h2same");
        h2mix  = self.list_pair_ss_cut.FindObject("hMggPt_Mixed").Clone("h2mix");
        h2same.Sumw2();
        h2mix .Sumw2();
        h2same.SetDirectory(0);
        h2mix .SetDirectory(0);
        if self.meson == "pi0":
            h2same.RebinX(2);
            h2mix .RebinX(2);
        elif self.meson == "eta":
            h2same.RebinX(4);
            h2mix .RebinX(4);
    
        nev = h1ev.GetBinContent(4);
        print("nev = {0:e}".format(nev));
    
        outlist.Add(h1ev);
        outlist.Add(h2same);
        outlist.Add(h2mix);

        npt = len(self.arr_pt);
    
        h1yield = TH1F("h1yield","raw yield"  ,npt-1, self.arr_pt);
        h1mean  = TH1F("h1mean" ,"peak mean"  ,npt-1, self.arr_pt);
        h1sigma = TH1F("h1sigma","peak sigma" ,npt-1, self.arr_pt);
        h1alpha = TH1F("h1alpha","alpha of CB",npt-1, self.arr_pt);
        h1n     = TH1F("h1n"    ,"n of CB"    ,npt-1, self.arr_pt);
    
        h1yield.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1yield.SetYTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}");
        h1mean.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1mean.SetYTitle("peak mean (GeV/#it{c}^{2})");
        h1sigma.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1sigma.SetYTitle("peak sigma (GeV/#it{c}^{2})");
        h1alpha.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1alpha.SetYTitle("#alpha of CB");
        h1n.SetXTitle("#it{p}_{T} (GeV/#it{c})");
        h1n.SetYTitle("n of CB");
    
        for i in range(0, npt-1):
            pt1 = self.arr_pt[i];
            pt2 = self.arr_pt[i+1];
    
            h1same = slice_histogram(h2same, pt1, pt2, "X", False);
            h1same.SetName("h1mgg_same_pt{0}".format(i));
            h1same.SetTitle("m_{{#gamma#gamma}}^{{same}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
            h1mix  = slice_histogram(h2mix , pt1, pt2, "X", False);
            h1mix.SetTitle("m_{{#gamma#gamma}}^{{mix}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
            h1mix.SetName("h1mgg_mix_pt{0}".format(i));
            h1same.SetDirectory(0);
            h1mix .SetDirectory(0);
    
            npair_same = h1same.GetEntries();
            npair_mix  = h1mix.GetEntries();
            if npair_mix < 1e-6:
                continue;
    
            h1mix.Scale(npair_same/npair_mix);
    
            h1ratio = get_ratio(h1same, h1mix);
            h1ratio.SetName("h1mgg_ratio_pt{0}".format(i));
            h1ratio.SetTitle("m_{{#gamma#gamma}}^{{ratio}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
            h1ratio .SetDirectory(0);
   
            height = 1.0;
            mean_init = 0.130;
            sigma_init = 0.008;
            if "pi0" in self.meson:
                mean_init = 0.135;
                sigma = 0.008;
            elif "eta" in self.meson:
                mean_init = 0.548;
                sigma_init = 0.012;
            bin_mean = h1ratio.FindBin(mean_init);
            height = h1ratio.GetBinContent(bin_mean) - 1.0;
    
            f1total = self.f1total.Clone("f1ratio_pt{0}".format(i));
            f1total.SetParameter(0,height);
            f1total.SetParameter(1,mean_init);
            f1total.SetParameter(2,sigma_init);
            f1total.SetParameter(3,0.6);
            f1total.SetParameter(4,0.1);
            f1total.SetParameter(5,0.1);
            f1total.SetParameter(6,-0.1);
            f1total.SetParLimits(0,1e-3,100);
            f1total.SetParLimits(1, mean_init - 3 * sigma_init, mean_init + 3 * sigma_init);
            f1total.SetParLimits(2, 0.5 * sigma_init, 2 * sigma_init);
            f1total.FixParameter(3,0.6);
            #f1total.FixParameter(4,1e+10);
            #f1total.SetParLimits(3,0, 100);
            f1total.SetParLimits(4,0, 100);
            h1ratio.Fit(f1total,"SME","",self.fit_min, self.fit_max);
    
            h1bkg = h1mix.Clone("h1bkg");
            f1bkg = self.f1bkg.Clone("f1bkg_pt{0}".format(i));
            npar_sig = self.f1sig.GetNpar();
            npar_bkg = self.f1bkg.GetNpar();
            for ip in range(npar_bkg):
                f1bkg.FixParameter(ip, f1total.GetParameter(ip + npar_sig));
            h1bkg.Multiply(f1bkg);
            h1bkg.SetName("h1mgg_bkg_pt{0}".format(i));
            h1bkg.SetTitle("m_{{#gamma#gamma}}^{{bkg}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
            h1bkg .SetDirectory(0);
    
            h1sig = get_bkg_subtracted(h1same, h1bkg);
            h1sig.SetName("h1mgg_sig_pt{0}".format(i));
            h1sig.SetTitle("m_{{#gamma#gamma}}^{{sig}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
            h1sig .SetDirectory(0);
    
            f1sig = self.f1sig.Clone("f1sig_pt{0}".format(i));
            f1sig.SetNpx(1000);
            for ip in range(npar_sig):
                f1sig.SetParameter(ip, f1total.GetParameter(ip));
    
            height = h1sig.GetBinContent(bin_mean);
            f1sig.SetParameter(0,height);
            f1sig.SetParameter(1,mean_init);
            f1sig.SetParameter(2,sigma_init);
            f1sig.SetParameter(3,0.6);
            f1sig.SetParameter(4,1);
            f1sig.SetParLimits(0,1,1e+6);
            f1sig.SetParLimits(1, mean_init - 3 * sigma_init, mean_init + 3 * sigma_init);
            f1sig.SetParLimits(2, 0.5 * sigma_init, 2 * sigma_init);
            f1sig.FixParameter(3,0.6);
            f1sig.SetParLimits(4,0, 100);
            h1sig.Fit(f1sig,"SME","",self.fit_min, self.fit_max);
    
            mean      = f1sig.GetParameter(1);
            mean_err  = f1sig.GetParError(1);
            sigma     = f1sig.GetParameter(2);
            sigma_err = f1sig.GetParError(2);
            alpha     = f1sig.GetParameter(3);
            alpha_err = f1sig.GetParError(3);
            n         = f1sig.GetParameter(4);
            n_err     = f1sig.GetParError(4);
    
            h1mean.SetBinContent(i+1, mean);
            h1mean.SetBinError(i+1, mean_err);
            h1sigma.SetBinContent(i+1, sigma);
            h1sigma.SetBinError(i+1, sigma_err);
            h1alpha.SetBinContent(i+1, alpha);
            h1alpha.SetBinError(i+1, alpha_err);
            h1n.SetBinContent(i+1, n);
            h1n.SetBinError(i+1, n_err);
    
            bin1 = h1sig.FindBin(mean + self.integral_min_sigma * sigma);
            bin2 = h1sig.FindBin(mean + self.integral_max_sigma * sigma);
            ry_err = ctypes.c_double(0);
            ry = h1sig.IntegralAndError(bin1,bin2,ry_err,"");
            h1yield.SetBinContent(i+1, ry);
            h1yield.SetBinError(i+1, ry_err);
    
            outlist.Add(h1same);
            outlist.Add(h1mix);
            outlist.Add(h1ratio);
            outlist.Add(f1total);
            outlist.Add(h1bkg);
            outlist.Add(h1sig);
            outlist.Add(f1sig);
            outlist.Add(f1bkg);
           
        h1yield.Scale(1/nev); 
        h1yield.Scale(1., "width"); 
        outlist.Add(h1yield);
        outlist.Add(h1mean );
        outlist.Add(h1sigma);
        outlist.Add(h1alpha);
        outlist.Add(h1n    );
    
        return outlist;

#___________________________________________________________________
if __name__ == "__main__":
    arr_pt = np.array([0.2,0.3,0.4,0.5], dtype=float);
    #ana = PairAnalyzer("pi0", "AnalysisResults_HL_75289.root", "pi0eta-to-gammagamma", "PCMPCM", "qc_qc", arr_pt);
    ana = PairAnalyzer("pi0", "AnalysisResults_HL_75289.root", "pi0eta-to-gammagamma");
    del ana;

    f1sig = TF1("f1sig","crystalball(0)",0,1);
    f1bkg = TF1("f1bkg","[0] + [1]*x",0,1);
    f1total = TF1("f1total","f1sig+f1bkg",0,1);
    print("initially, ", f1total.GetExpFormula("p"));

    f1total.SetParameters(1, 0.135, 0.005, 0.6, 1, 1 , 1);
    #f1total.SetParameters(1, 1,1, 0.6,  0.135, 1, 0.005);
    f1total.Draw();
    print("later, ", f1total.GetExpFormula("p"));
