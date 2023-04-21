import os, sys, shutil
import math
import numpy as np
import ctypes
import ROOT
from ROOT import TH1D, TH2D, TH3D, TList, TF1
from ROOT import gROOT, gSystem, gStyle, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange
from histo_manager import slice_histogram, rebin_histogram, get_bkg_subtracted, get_ratio

#_________________________________________________________________________________________
def analyze_ptspectrum(rootfile, ssname, cutname, arr_pt):
    print(sys._getframe().f_code.co_name);
    outlist = TList();
    outlist.SetName(cutname);
    print(arr_pt);
    npt = len(arr_pt);

    dir_main = rootfile.Get("pi0eta-to-gammagamma");
    list_pair = dir_main.Get("Pair");
    list_ev   = dir_main.Get("Event");

    list_pair.Print();
    list_ev  .Print();
    
    list_pair_subsystem = list_pair.FindObject(ssname);
    list_ev_subsystem   = list_ev.FindObject(ssname);
    list_pair_subsystem.Print();
    list_ev_subsystem.Print();

    list_cut = list_pair_subsystem.FindObject(cutname);
    list_cut.Print();

    h1ev   = list_ev_subsystem.FindObject("hCollisionCounter");
    h2same = list_cut.FindObject("hMggPt_Same");
    h2mix  = list_cut.FindObject("hMggPt_Mixed");
    h2same.Sumw2();
    h2mix .Sumw2();
    h2same.SetDirectory(0);
    h2mix .SetDirectory(0);
    h2same.RebinX(2);
    h2mix .RebinX(2);

    nev = h1ev.GetBinContent(4);
    print("nev = {0:e}".format(nev));

    outlist.Add(h1ev);
    outlist.Add(h2same);
    outlist.Add(h2mix);

    h1yield = TH1D("h1yield","raw yield",npt-1, arr_pt);
    h1mean  = TH1D("h1mean" ,"peak mean",npt-1, arr_pt);
    h1sigma = TH1D("h1sigma","peak sigma",npt-1, arr_pt);
    h1alpha = TH1D("h1alpha","alpha of CB",npt-1, arr_pt);
    h1n     = TH1D("h1n"    ,"n of CB",npt-1, arr_pt);

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
        pt1 = arr_pt[i];
        pt2 = arr_pt[i+1];

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

        f1ratio = TF1("f1ratio_pt{0}".format(i),"crystalball(0) + pol1(5)",0,1); #height, mu, sigma, alpha, n
        f1ratio.SetNpx(1000);
        #print(f1ratio.GetExpFormula(""));

        bin135 = h1ratio.FindBin(0.135);
        bin200 = h1ratio.FindBin(0.2);
        height = h1ratio.GetBinContent(bin135) - h1ratio.GetBinContent(bin200);

        f1ratio.SetParameter(0,height);
        f1ratio.SetParameter(1,0.13);
        f1ratio.SetParameter(2,0.01);
        f1ratio.SetParameter(3,0.5);
        f1ratio.SetParameter(4,0.1);
        f1ratio.SetParameter(5,0.1);
        f1ratio.SetParameter(6,-0.1);
        f1ratio.SetParLimits(0,1e-3,10);
        f1ratio.SetParLimits(1,0.12, 0.16);
        f1ratio.SetParLimits(2,0.003, 0.02);
        f1ratio.FixParameter(3,0.6);
        #f1ratio.FixParameter(4,1e+10);
        #f1ratio.SetParLimits(3,0, 100);
        f1ratio.SetParLimits(4,0, 100);
        h1ratio.Fit(f1ratio,"SME","",0.04, 0.24);

        h1bkg = h1mix.Clone("h1bkg");
        f1bkg = TF1("f1bkg_pt{0}".format(i),"pol1(0)",0,1);
        for ip in range(f1bkg.GetNpar()):
            f1bkg.FixParameter(ip, f1ratio.GetParameter(ip + 5));
        h1bkg.Multiply(f1bkg);
        h1bkg.SetName("h1mgg_bkg_pt{0}".format(i));
        h1bkg.SetTitle("m_{{#gamma#gamma}}^{{bkg}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
        h1bkg .SetDirectory(0);

        h1sig = get_bkg_subtracted(h1same, h1bkg);
        h1sig.SetName("h1mgg_sig_pt{0}".format(i));
        h1sig.SetTitle("m_{{#gamma#gamma}}^{{sig}}, {0:2.1f} < #it{{p}}_{{T,#gamma#gamma}} < {1:2.1f} GeV/#it{{c}}".format(pt1, pt2));
        h1sig .SetDirectory(0);

        f1sig = TF1("f1sig_pt{0}".format(i),"crystalball(0)",0,1);
        f1sig.SetNpx(1000);
        for ip in range(f1sig.GetNpar()):
            f1sig.SetParameter(ip, f1ratio.GetParameter(ip));

        bin135 = h1sig.FindBin(0.135);
        height = h1sig.GetBinContent(bin135);
        f1sig.SetParameter(0,height);
        f1sig.SetParameter(1,0.13);
        f1sig.SetParameter(2,0.01);
        f1sig.SetParameter(3,0.5);
        f1sig.SetParameter(4,1);
        f1sig.SetParLimits(0,1,1e+7);
        f1sig.SetParLimits(1,0.12, 0.16);
        f1sig.SetParLimits(2,0.003, 0.02);
        f1sig.FixParameter(3,0.6);
        f1sig.SetParLimits(4,0, 100);
        h1sig.Fit(f1sig,"SME","",0.04, 0.24);

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

        bin1 = h1sig.FindBin(mean - 3*sigma);
        bin2 = h1sig.FindBin(mean + 3*sigma);
        ry_err = ctypes.c_double(0);
        ry = h1sig.IntegralAndError(bin1,bin2,ry_err,"");
        h1yield.SetBinContent(i+1, ry);
        h1yield.SetBinError(i+1, ry_err);

        outlist.Add(h1same);
        outlist.Add(h1mix);
        outlist.Add(h1ratio);
        outlist.Add(f1ratio);
        outlist.Add(h1bkg);
        outlist.Add(h1sig);
        outlist.Add(f1sig);
       
    h1yield.Scale(1/nev); 
    h1yield.Scale(1., "width"); 
    outlist.Add(h1yield);
    outlist.Add(h1mean );
    outlist.Add(h1sigma);
    outlist.Add(h1alpha);
    outlist.Add(h1n    );

    return outlist;
#_________________________________________________________________________________________
def analyze_mee_ptee_efficiency(rootfile,cutname,arr_mee, arr_ptee):
    print(sys._getframe().f_code.co_name);
    outlist = TList();
    outlist.SetName(cutname);
    return outlist;
#_________________________________________________________________________________________
