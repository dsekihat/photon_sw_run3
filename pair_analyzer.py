import os, sys, shutil
import numpy as np
import math
import ctypes
import ROOT
from ROOT import TFile, TDirectory, THashList, TF1, TH1F
from histo_manager import slice_histogram, rebin_histogram, get_bkg_subtracted, get_ratio
from nm_fitter import NMFitter

class PairAnalyzer:
    def __init__(self):
        print("default constructor is called");

    def __init__(self, particle, filename, dirname, ismc):
        print("target particle = {0} , filename = {1} , dirname = {2}".format(particle, filename, dirname));
        self.particle = particle;
        self.rootfile = TFile.Open(filename, "READ");
        self.rootdir = self.rootfile.Get(dirname);
        self.list_ev = self.rootdir.Get("Event");
        self.list_pair = self.rootdir.Get("Pair");
        self.arr_pt = np.array([0,1,2,3,4,5], dtype=float);
        self.ismc = ismc;
        self.list_gen = None;
        if self.ismc:
            self.list_gen = self.rootdir.Get("Generated");

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
        self.is_cb_n_fixed = False;
        self.is_cb_alpha_fixed = False;
        self.cb_alpha = 0.0;
        self.cb_n = 0.0;

    def __del__(self):
        if self.rootfile.IsOpen():
            print("close input root file.");
            self.rootfile.Close();

    def set_arr_pt(self, arr_pt):
        self.arr_pt = arr_pt;

    def set_subsystem(self, ssname):
        self.ssname = ssname;
        self.list_ev_ss   = self.list_ev.FindObject(ssname).FindObject("after");
        self.list_pair_ss = self.list_pair.FindObject(ssname);

    def set_cutname(self, cutname):
        self.cutname = cutname;
        if self.list_ev_ss is None or self.list_pair_ss is None:
            print("Please define subsystem name first!");
            return None;
        #self.list_pair_ss_cut = self.list_pair_ss.FindObject(cutname);
        self.list_pair_ss_cut = self.list_pair_ss.FindObject(cutname).FindObject("nocut");

    def set_fit_range(self, fit_min, fit_max):
        self.fit_min = fit_min;
        self.fit_max = fit_max;

    def set_integral_range(self, integral_min, integral_max):
        self.integral_min_sigma = integral_min;
        self.integral_max_sigma = integral_max;

    def set_fit_function(self, sig, bkg):
        self.sig_name = sig;
        self.bkg_name = bkg;
        if sig == "cb":
            self.f1sig = TF1("f1sig","crystalball(0)", 0,1);
            if bkg == "pol1":
                self.f1bkg   = TF1("f1bkg","pol1(0)", 0,1);
                self.f1total = TF1("f1total","crystalball(0) + pol1(5)", 0,1);
            elif bkg == "pol2":
                self.f1bkg   = TF1("f1bkg","pol2(0)", 0,1);
                self.f1total = TF1("f1total","crystalball(0) + pol2(5)", 0,1);
            self.f1total.SetParameter(3,0.5);
            self.f1total.SetParameter(4,10.0);
            self.f1sig.SetParameter(3,0.5);
            self.f1sig.SetParameter(4,10.0);
            self.f1total.SetParLimits(3,0, 100);
            self.f1total.SetParLimits(4,0, 100);
            self.f1sig.SetParLimits(3,0, 100);
            self.f1sig.SetParLimits(4,0, 100);
        self.f1sig.SetNpx(1000);
        self.f1bkg.SetNpx(1000);
        self.f1total.SetNpx(1000);
        print("initially, ", self.f1total.GetExpFormula(""));

    def fix_cb_alpha(self, alpha):
        self.cb_alpha = alpha;
        self.is_cb_alpha_fixed = True;
        self.f1sig.FixParameter(3, alpha);
        self.f1total.FixParameter(3, alpha);
    def fix_cb_n(self, n):
        self.is_cb_n_fixed = True;
        self.cb_n = n;
        self.f1sig.FixParameter(4, n);
        self.f1total.FixParameter(4, n);

    def set_xtitle(self, title):
        self.xtitle = title;

    def set_ytitle(self, title):
        self.ytitle = title;

    def analyze_ptspectrum(self): #this is main function
        print(sys._getframe().f_code.co_name);
        outlist = THashList();
        #outlist.SetOwner(True);
        outlist.SetName("outlist");

        h1ev   = self.list_ev_ss.FindObject("hCollisionCounter").Clone("h1ev");
        h2same = self.list_pair_ss_cut.FindObject("hMggPt_Same").Clone("h2same");
        h2mix  = self.list_pair_ss_cut.FindObject("hMggPt_Mixed").Clone("h2mix");
        if self.ssname == "EMCEMC":
            h2mix  = self.list_pair_ss_cut.FindObject("hMggPt_Same_RotatedBkg").Clone("h2mix");

        h2same.Sumw2();
        h2mix .Sumw2();
        h2same.SetDirectory(0);
        h2mix .SetDirectory(0);
        if self.particle == "pi0":
            h2same.RebinX(2);
            h2mix .RebinX(2);
        elif self.particle == "eta":
            h2same.RebinX(2);
            h2mix .RebinX(2);
  
        nev = h1ev.GetBinContent(7);
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
            bw = h1same.GetBinWidth(1);
            h1same.SetXTitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
            h1same.SetYTitle("counts / {0:d} MeV/#it{{c}}^{{2}}".format(int(bw*1e+3)));
            h1mix  = slice_histogram(h2mix , pt1, pt2, "X", False);
            h1mix.SetName("h1mgg_mix_pt{0}".format(i));
            h1mix.SetTitle("m_{{#gamma#gamma}}^{{mix}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
            h1mix.SetXTitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
            h1mix.SetYTitle("counts / {0:d} MeV/#it{{c}}^{{2}}".format(int(bw*1e+3)));
            h1same.SetDirectory(0);
            h1mix .SetDirectory(0);
 
            npair_same = h1same.GetEntries();
            npair_mix  = h1mix.GetEntries();
            h1mix_scaled = h1mix.Clone("h1mix_scaled_pt{0}".format(i));
            h1mix_scaled.Scale(npair_same/npair_mix);

            height = 1.0;
            mean_init = 0.130;
            sigma_init = 0.008;
            if "pi0" in self.particle:
                mean_init = 0.132;
                sigma_init = 0.005;
            elif "eta" in self.particle:
                mean_init = 0.545;
                sigma_init = 0.008;

            nmf = NMFitter(h1same, h1mix_scaled, self.sig_name, self.bkg_name); 
            nmf.set_parameters(mean_init, sigma_init, self.cb_alpha, self.cb_n, self.is_cb_alpha_fixed, self.is_cb_n_fixed);
            fit_result = nmf.fit("SME", "SME", "", self.fit_min, self.fit_max);
            h1sig = fit_result[1];
            h1bkg = fit_result[2];
            h1ratio = fit_result[3];
            f1sig = fit_result[4];
            f1bkg = fit_result[5];
            f1total = fit_result[6];

            h1bkg.SetName("h1mgg_bkg_pt{0}".format(i));
            h1bkg.SetTitle("m_{{#gamma#gamma}}^{{bkg}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
            h1bkg.SetXTitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
            h1bkg.SetYTitle("counts / {0:d} MeV/#it{{c}}^{{2}}".format(int(bw*1e+3)));
            h1bkg .SetDirectory(0);

            h1sig.SetName("h1mgg_sig_pt{0}".format(i));
            h1sig.SetTitle("m_{{#gamma#gamma}}^{{sig}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
            h1sig.SetXTitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
            h1sig.SetYTitle("counts / {0:d} MeV/#it{{c}}^{{2}}".format(int(bw*1e+3)));
            h1sig.SetDirectory(0);

            h1ratio.SetName("h1mgg_ratio_pt{0}".format(i));
            h1ratio.SetTitle("m_{{#gamma#gamma}}^{{ratio}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
            h1ratio.SetXTitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
            h1ratio.SetYTitle("counts / {0:d} MeV/#it{{c}}^{{2}}".format(int(bw*1e+3)));
            h1ratio.SetDirectory(0);

            f1sig.SetName("f1sig_pt{0}".format(i));
            f1sig.SetName("f1sig_pt{0}".format(i));
            f1total.SetName("f1ratio_pt{0}".format(i));
 
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

    def analyze_ptspectrum_efficiency(self, isTaggingPi0 = False): #this is main function
        print(sys._getframe().f_code.co_name);
        outlist = THashList();
        #outlist.SetOwner(True);
        outlist.SetName("outlist");

        h1ev   = self.list_ev_ss.FindObject("hCollisionCounter").Clone("h1ev");
        h2same = self.list_pair_ss_cut.FindObject("hMggPt_Pi0_Primary").Clone("hMggPt_Pi0_Primary");
        h2same_fd = self.list_pair_ss_cut.FindObject("hMggPt_Pi0_FromWD").Clone("hMggPt_Pi0_FromWD");
        if "eta" in self.particle.lower():
            h2same = self.list_pair_ss_cut.FindObject("hMggPt_Eta_Primary").Clone("hMggPt_Eta_Primary");

        h2same.Sumw2();
        h2same.SetDirectory(0);
        if self.particle == "pi0":
            h2same.RebinX(2);
            h2same_fd.RebinX(2);
        elif self.particle == "eta":
            h2same.RebinX(2);
            h2same_fd.RebinX(2);
    
        nev = h1ev.GetBinContent(7);
        print("nev = {0:e}".format(nev));
    
        outlist.Add(h1ev);
        outlist.Add(h2same);
        if self.particle == "pi0":
            outlist.Add(h2same_fd);

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
        h1yield_fd = h1yield.Clone("h1yield_fd");
    
        for i in range(0, npt-1):
            pt1 = self.arr_pt[i];
            pt2 = self.arr_pt[i+1];
    
            h1same = slice_histogram(h2same, pt1, pt2, "X", False);
            h1same.SetName("h1mgg_same_pt{0}".format(i));
            h1same.SetTitle("m_{{#gamma#gamma}}^{{same}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
            bw = h1same.GetBinWidth(1);
            h1same.SetXTitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
            h1same.SetYTitle("counts / {0:d} MeV/#it{{c}}^{{2}}".format(int(bw*1e+3)));
            h1same.SetDirectory(0);

            height = 1.0;
            mean_init = 0.130;
            sigma_init = 0.008;
            if "pi0" in self.particle:
                mean_init = 0.132;
                sigma_init = 0.005;
            elif "eta" in self.particle:
                mean_init = 0.545;
                sigma_init = 0.008;
            bin_mean = h1same.FindBin(mean_init);
            height   = h1same.GetBinContent(bin_mean);
            h1sig = h1same.Clone("h1mgg_pi0_pt{0}".format(i));
            h1sig.SetName("h1mgg_sig_pt{0}".format(i));
            h1sig.SetTitle("m_{{#gamma#gamma}}^{{sig}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
            h1sig.SetXTitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
            h1sig.SetYTitle("counts / {0:d} MeV/#it{{c}}^{{2}}".format(int(bw*1e+3)));
            h1sig .SetDirectory(0);
    
            f1sig = self.f1sig.Clone("f1sig_pt{0}".format(i));
            f1sig.SetNpx(1000);
    
            height = h1sig.GetBinContent(bin_mean);
            f1sig.SetParameter(0,height);
            f1sig.SetParameter(1,mean_init);
            f1sig.SetParameter(2,sigma_init);
            #f1sig.SetParameter(3,0.6);
            #f1sig.SetParameter(4,1);
            f1sig.SetParLimits(0,1,1e+6);
            f1sig.SetParLimits(1, mean_init - 3 * sigma_init, mean_init + 3 * sigma_init);
            f1sig.SetParLimits(2, 0.5 * sigma_init, 2 * sigma_init);
            #f1sig.FixParameter(3,0.5);
            #f1sig.FixParameter(4,10.0);
            #f1sig.SetParLimits(3,0, 100);
            #f1sig.SetParLimits(4,0, 100);
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
            outlist.Add(h1sig);
            outlist.Add(f1sig);

            if self.particle == "pi0":
                h1same_fd = slice_histogram(h2same_fd, pt1, pt2, "X", False);
                h1same_fd.SetName("h1mgg_same_fd_pt{0}".format(i));
                h1same_fd.SetTitle("m_{{#gamma#gamma}}^{{same}} weak decay, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
                h1same_fd.SetXTitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
                h1same_fd.SetYTitle("counts / {0:d} MeV/#it{{c}}^{{2}}".format(int(bw*1e+3)));
                h1same_fd.SetDirectory(0);
                outlist.Add(h1same_fd);
                bin1 = h1same_fd.FindBin(mean + self.integral_min_sigma * sigma);
                bin2 = h1same_fd.FindBin(mean + self.integral_max_sigma * sigma);
                ry_err = ctypes.c_double(0);
                ry = h1same_fd.IntegralAndError(bin1,bin2,ry_err,"");
                h1yield_fd.SetBinContent(i+1, ry);
                h1yield_fd.SetBinError(i+1, ry_err);

        h1yield.Scale(1/nev); 
        h1yield.Scale(1., "width"); 
        outlist.Add(h1yield);
        outlist.Add(h1mean );
        outlist.Add(h1sigma);
        outlist.Add(h1alpha);
        outlist.Add(h1n    );

        if self.particle == "pi0":
            h1yield_fd.Scale(1/nev); 
            h1yield_fd.Scale(1., "width"); 
            outlist.Add(h1yield_fd);
            h1sum = h1yield_fd.Clone("h1sum");
            h1sum.Sumw2();
            h1sum.Add(h1yield, 1.);
            h1fd = h1sum.Clone("h1fd");
            h1fd.Reset();
            h1fd.Sumw2();
            h1fd.Divide(h1yield_fd, h1sum, 1., 1., "B");
            h1fd.SetTitle("feed down correction for #pi^{0}");
            h1fd.SetYTitle("#frac{rec. #pi^{0} from weak decay}{all rec. #pi^{0}}");
            outlist.Add(h1fd);

        #Next, generated information
        if isTaggingPi0:
            list_pcm_cut = self.rootdir.Get("PCM").FindObject("qc");
            h1pt_org = list_pcm_cut.FindObject("hPt_v0photon_Pi0_Primary").Clone("hPt_v0photon_Pi0_Primary");
            h1pt_org .Sumw2();
            outlist.Add(h1pt_org);
            #arr_tmp = np.array(h1pt_org.GetXaxis().GetXbins(), dtype=float);
            #print(arr_tmp);
            h1dndpt = rebin_histogram(h1pt_org, self.arr_pt, True, False);
            h1dndpt.Scale(1/nev);
            h1dndpt.SetName("h1dndpt_gen");
            h1dndpt.SetXTitle("#it{p}_{T,#gamma} (GeV/#it{c})");
            h1dndpt.SetYTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T,#gamma}} (GeV/#it{c})^{-1}");
            outlist.Add(h1dndpt);
            h1eff = h1yield.Clone("h1eff");
            h1eff.Sumw2();
            h1eff.Reset();
            h1eff.SetTitle("conditional acceptance");
            h1eff.SetXTitle("#it{p}_{T,#gamma} (GeV/#it{c})");
            h1eff.SetYTitle("conditional acceptance <#varepsilon_{#gamma}f>");
            h1eff.Reset();
            h1eff.Divide(h1yield, h1dndpt, 1., 1., "B");
            outlist.Add(h1eff);
        else:    
            h1pt_org = None;
            h1y_org = None;
            h1phi_org = None;
            if "pi0" in self.particle.lower():
                h1pt_org  = self.list_gen.FindObject(self.ssname).FindObject("hPt_Pi0").Clone("hPt_Pi0");
                h1y_org   = self.list_gen.FindObject(self.ssname).FindObject("hY_Pi0").Clone("hY_Pi0");
                h1phi_org = self.list_gen.FindObject(self.ssname).FindObject("hPhi_Pi0").Clone("hPhi_Pi0");
                h1pt_acc_org  = self.list_gen.FindObject(self.ssname).FindObject("hPt_Pi0_Acc").Clone("hPt_Pi0_Acc");
                h1y_acc_org   = self.list_gen.FindObject(self.ssname).FindObject("hY_Pi0_Acc").Clone("hY_Pi0_Acc");
                h1phi_acc_org = self.list_gen.FindObject(self.ssname).FindObject("hPhi_Pi0_Acc").Clone("hPhi_Pi0_Acc");
            elif "eta" in self.particle.lower():
                h1pt_org  = self.list_gen.FindObject(self.ssname).FindObject("hPt_Eta").Clone("hPt_Eta");
                h1y_org   = self.list_gen.FindObject(self.ssname).FindObject("hY_Eta").Clone("hY_Eta");
                h1phi_org = self.list_gen.FindObject(self.ssname).FindObject("hPhi_Eta").Clone("hPhi_Eta");
                h1pt_acc_org  = self.list_gen.FindObject(self.ssname).FindObject("hPt_Eta_Acc").Clone("hPt_Eta_Acc");
                h1y_acc_org   = self.list_gen.FindObject(self.ssname).FindObject("hY_Eta_Acc").Clone("hY_Eta_Acc");
                h1phi_acc_org = self.list_gen.FindObject(self.ssname).FindObject("hPhi_Eta_Acc").Clone("hPhi_Eta_Acc");
            else:
                print("Please choose pi0 or eta.");
                return outlist;
            h1pt_org .Sumw2();
            h1y_org  .Sumw2();
            h1phi_org.Sumw2();
            h1pt_acc_org .Sumw2();
            h1y_acc_org  .Sumw2();
            h1phi_acc_org.Sumw2();
            outlist.Add(h1pt_org);
            outlist.Add(h1y_org);
            outlist.Add(h1phi_org);
            outlist.Add(h1pt_acc_org);
            outlist.Add(h1y_acc_org);
            outlist.Add(h1phi_acc_org);

            h1dndpt = rebin_histogram(h1pt_org, self.arr_pt, True, False);
            h1dndpt.Scale(1/nev);
            h1dndpt.SetName("h1dndpt_gen");
            h1dndpt.SetXTitle("#it{p}_{T} (GeV/#it{c})");
            h1dndpt.SetYTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}");
            outlist.Add(h1dndpt);

            h1dndpt_acc = rebin_histogram(h1pt_acc_org, self.arr_pt, True, False);
            h1dndpt_acc.Scale(1/nev);
            h1dndpt_acc.SetName("h1dndpt_gen_acc");
            h1dndpt_acc.SetXTitle("#it{p}_{T} (GeV/#it{c})");
            h1dndpt_acc.SetYTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}");
            outlist.Add(h1dndpt_acc);

            h1acc = h1dndpt_acc.Clone("h1acc");
            h1acc.Sumw2();
            h1acc.Reset();
            h1acc.SetTitle("br. #times acceptance");
            h1acc.SetYTitle("br. #times acceptance");
            h1acc.Reset();
            h1acc.Divide(h1dndpt_acc, h1dndpt, 1., 1., "B");
            outlist.Add(h1acc);

            h1rec = h1yield.Clone("hrec");
            h1rec.Sumw2();
            h1rec.Reset();
            h1rec.SetTitle("rec. efficiency");
            h1rec.SetYTitle("rec. efficiency");
            h1rec.Reset();
            h1rec.Divide(h1yield, h1dndpt_acc, 1., 1., "B");
            outlist.Add(h1rec);

            h1eff = h1yield.Clone("h1eff");
            h1eff.Sumw2();
            h1eff.Reset();
            h1eff.SetTitle("efficiency");
            h1eff.SetYTitle("br. #times acc. #times rec. efficiency");
            h1eff.Reset();
            h1eff.Divide(h1yield, h1dndpt, 1., 1., "B");
            outlist.Add(h1eff);
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
