import os, sys, shutil
import math
import argparse
import numpy as np
import ctypes
import yaml
import ROOT
#ROOT.gROOT.SetBatch(True);
from ROOT import TFile, THashList, TF1, TGraphAsymmErrors, TMultiGraph
from pair_analyzer import PairAnalyzer
from histo_manager import slice_histogram, rebin_histogram, get_bkg_subtracted, get_ratio, get_R_factor, get_corrected_bkg, get_corrected_bkg_simple, get_significance
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue

#parser = argparse.ArgumentParser('Example program');
#parser.add_argument("-i", "--input" , default="AnalysisResults.root", type=str, help="path to the root file you want to analyze", required=True)
#parser.add_argument("-c", "--config", default="config.yml", type=str, help="path to the *.yml configuration file", required=True)
#parser.add_argument("-t", "--type"  , default="data" , type=str, help="run type [data or mc]", required=True)
#parser.add_argument("-s", "--suffix"  , default="" , type=str, help="suffix for output file name", required=False)
#args = parser.parse_args();
#
#filename = args.input;
#with open(args.config, "r", encoding="utf-8") as config_yml:
#    config = yaml.safe_load(config_yml)

#____________________________________________________________________________________________
def analyze_mee_vs_single_dca3d(filename, cutname, arr_mee):
    print(sys._getframe().f_code.co_name);

    #arr_pair_dca3d_min   = np.array([  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64);
    #arr_pair_dca3d_max   = np.array([999.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0], dtype=np.float64);
    #arr_single_dca3d_min = np.array([  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64);
    #arr_single_dca3d_max = np.array([999.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0], dtype=np.float64);

    arr_pair_dca3d_min   = np.array([  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64);
    arr_pair_dca3d_max   = np.array([999.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0], dtype=np.float64);
    arr_single_dca3d_min = np.array([  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64);
    arr_single_dca3d_max = np.array([999.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0], dtype=np.float64);

    #arr_pair_dca3d_min   = np.array([  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64);
    #arr_pair_dca3d_max   = np.array([999.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0], dtype=np.float64);
    #arr_single_dca3d_min = np.array([  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64);
    #arr_single_dca3d_max = np.array([999.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0], dtype=np.float64);

    arr_mee_min = np.array([0.0,  0.0, 0.04, 0.08, 0.14, 0.14, 0.5, 1.1, 1.1, 2.7], dtype=np.float64);
    arr_mee_max = np.array([3.5, 0.14, 0.14, 0.14, 0.30, 0.50, 1.1, 2.5, 2.7, 3.2], dtype=np.float64);
    id_prompt = 3;
    id_nonprompt = 7;

    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("dalitz-ee-qc");
    rootdire.ls();
    list_ee = rootdire.Get("DalitzEE");
    list_ee.ls();
    list_ee_cut = list_ee.FindObject("Cent_FT0C_0.00_80.00").FindObject(cutname);
    list_ee_cut.ls();

    list_ev = rootdire.Get("Event").FindObject("Cent_FT0C_0.00_80.00");
    list_ev.ls();
    h1mult = list_ev.FindObject("hMultNTracksPV");
    nev = h1mult.GetEntries();
    print("nev = {0:f} B events".format(nev/1e+9));

    outname = "mee_dca3d_{0}_{1}TeV_{2}.root".format("PbPb", 5.36, "LHC23zs");
    print(outname);
    outfile = TFile(outname, "RECREATE");
    outlist = THashList();
    outlist.SetName(cutname);
    outlist.SetOwner(True);

    hs_dilepton_uls_dca_same  = list_ee_cut.FindObject("hs_dilepton_uls_dca_same");
    hs_dilepton_lspp_dca_same = list_ee_cut.FindObject("hs_dilepton_lspp_dca_same");
    hs_dilepton_lsmm_dca_same = list_ee_cut.FindObject("hs_dilepton_lsmm_dca_same");
    hs_dilepton_uls_dca_same  = list_ee_cut.FindObject("hs_dilepton_uls_same");
    hs_dilepton_lspp_dca_same = list_ee_cut.FindObject("hs_dilepton_lspp_same");
    hs_dilepton_lsmm_dca_same = list_ee_cut.FindObject("hs_dilepton_lsmm_same");
    hs_dilepton_uls_dca_same .Sumw2();
    hs_dilepton_lspp_dca_same.Sumw2();
    hs_dilepton_lsmm_dca_same.Sumw2();

    #hs_dilepton_uls_same = list_ee_cut.FindObject("hs_dilepton_uls_same");
    #h2tmp = hs_dilepton_uls_same.Projection(0,3);
    #outlist.Add(h2tmp);

    #hs_dilepton_uls_dca_mix = list_ee_cut.FindObject("hs_dilepton_uls_dca_mix");
    #hs_dilepton_lspp_dca_mix = list_ee_cut.FindObject("hs_dilepton_lspp_dca_mix");
    #hs_dilepton_lsmm_dca_mix = list_ee_cut.FindObject("hs_dilepton_lsmm_dca_mix");

    list_h1m_pair_dca = [];
    list_h1m_single_dca = [];

    n_pair_dca = len(arr_pair_dca3d_min);
    for idca in range(0, n_pair_dca):
        dca_min = arr_pair_dca3d_min[idca];
        dca_max = arr_pair_dca3d_max[idca];
        outlist_mee = get_mee_sliced_by_pair_dca(hs_dilepton_uls_dca_same, hs_dilepton_lspp_dca_same, hs_dilepton_lsmm_dca_same, dca_min, dca_max, nev, arr_mee);
        h1m = outlist_mee.FindObject("h1sig").Clone("h1sig_pair_dca{0:d}".format(idca));
        list_h1m_pair_dca.append(h1m);
        outlist.Add(outlist_mee);

#    n_single_dca = len(arr_single_dca3d_min);
#    for idca in range(0, n_single_dca):
#        dca_min = arr_single_dca3d_min[idca];
#        dca_max = arr_single_dca3d_max[idca];
#        outlist_mee = get_mee_sliced_by_single_dca(hs_dilepton_uls_dca_same, hs_dilepton_lspp_dca_same, hs_dilepton_lsmm_dca_same, dca_min, dca_max, nev, arr_mee);
#        h1m = outlist_mee.FindObject("h1sig").Clone("h1sig_single_dca{0:d}".format(idca));
#        list_h1m_single_dca.append(h1m);
#        outlist.Add(outlist_mee);
#
#    n_mee = len(arr_mee_min);
#    for im in range(0, n_mee):
#        m_min = arr_mee_min[im];
#        m_max = arr_mee_max[im];
#        outlist_pair_dca = get_pair_dca_sliced_by_mee(hs_dilepton_uls_dca_same, hs_dilepton_lspp_dca_same, hs_dilepton_lsmm_dca_same, m_min, m_max, nev);
#        outlist.Add(outlist_pair_dca);
#
#    for im in range(0, n_mee):
#        m_min = arr_mee_min[im];
#        m_max = arr_mee_max[im];
#        outlist_single_dca = get_single_dca_sliced_by_mee(hs_dilepton_uls_dca_same, hs_dilepton_lspp_dca_same, hs_dilepton_lsmm_dca_same, m_min, m_max, nev);
#        outlist.Add(outlist_single_dca);
#
#    mg_pair = get_eff_correlation(list_h1m_pair_dca, nev, kBlue+4, 20);
#    mg_pair.SetName("mg_pair");
#    outlist.Add(mg_pair);
#
#    mg_single = get_eff_correlation(list_h1m_single_dca, nev, kRed+4, 21);
#    mg_single.SetName("mg_single");
#    outlist.Add(mg_single);

    outfile.WriteTObject(outlist);
    outlist.Clear();

    rootfile.Close();
    outfile.Close();
#____________________________________________________________________________________________
def get_mee_sliced_by_pair_dca(hs_uls_same, hs_lspp_same, hs_lsmm_same, dca_min, dca_max, nev, arr_mee):
    outlist = THashList();
    outlist.SetName("mee_in_pair_dca_{0:2.1f}_{1:2.1f}_sigma".format(dca_min, dca_max));
    outlist.SetOwner(True);

    #bin1 = hs_uls_same.GetAxis(3).FindBin(dca_min + 1e-3);
    #bin2 = hs_uls_same.GetAxis(3).FindBin(dca_max - 1e-3);
    bin1 = hs_uls_same.GetAxis(2).FindBin(dca_min + 1e-3);
    bin2 = hs_uls_same.GetAxis(2).FindBin(dca_max - 1e-3);

    #hs_uls_same .GetAxis(3).SetRange(bin1, bin2);
    #hs_lspp_same.GetAxis(3).SetRange(bin1, bin2);
    #hs_lsmm_same.GetAxis(3).SetRange(bin1, bin2);
    hs_uls_same .GetAxis(2).SetRange(bin1, bin2);
    hs_lspp_same.GetAxis(2).SetRange(bin1, bin2);
    hs_lsmm_same.GetAxis(2).SetRange(bin1, bin2);

    #mass spectrum integrated over DCA
    h1m_uls_same_org  = hs_uls_same.Projection(0);
    h1m_lspp_same_org = hs_lspp_same.Projection(0);
    h1m_lsmm_same_org = hs_lsmm_same.Projection(0);
    h1m_uls_same_org .SetName("h1m_uls_same_org");
    h1m_lspp_same_org.SetName("h1m_lspp_same_org");
    h1m_lsmm_same_org.SetName("h1m_lsmm_same_org");
    h1m_uls_same  = rebin_histogram(h1m_uls_same_org, arr_mee, False, False);
    h1m_lspp_same = rebin_histogram(h1m_lspp_same_org, arr_mee, False, False);
    h1m_lsmm_same = rebin_histogram(h1m_lsmm_same_org, arr_mee, False, False);

    h1m_uls_same .SetName("h1m_uls_same");
    h1m_lspp_same.SetName("h1m_lspp_same");
    h1m_lsmm_same.SetName("h1m_lsmm_same");
    outlist.Add(h1m_uls_same );
    outlist.Add(h1m_lspp_same);
    outlist.Add(h1m_lsmm_same);

    #h1m_uls_mix_org = hs_dilepton_uls_dca_mix.Projection(0);
    #h1m_lspp_mix_org = hs_dilepton_lspp_dca_mix.Projection(0);
    #h1m_lsmm_mix_org = hs_dilepton_lsmm_dca_mix.Projection(0);
    #h1m_uls_mix  = rebin_histogram(h1m_uls_mix_org, arr_mee, False, False);
    #h1m_lspp_mix = rebin_histogram(h1m_lspp_mix_org, arr_mee, False, False);
    #h1m_lsmm_mix = rebin_histogram(h1m_lsmm_mix_org, arr_mee, False, False);
    #h1m_uls_mix .SetName("h1m_uls_mix");
    #h1m_lspp_mix.SetName("h1m_lspp_mix");
    #h1m_lsmm_mix.SetName("h1m_lsmm_mix");
    #outlist.Add(h1m_uls_mix );
    #outlist.Add(h1m_lspp_mix);
    #outlist.Add(h1m_lsmm_mix);
    #h1R = get_R_factor(h1m_uls_mix, None, h1m_lspp_mix, h1m_lsmm_mix);
    #h1R.SetName("h1R");
    #outlist.Add(h1R);
    #h1bkg = get_corrected_bkg(h1R, h1m_lspp_same, h1m_lsmm_same);

    h1bkg = get_corrected_bkg_simple(1.0, 0.0, h1m_lspp_same, h1m_lsmm_same);
    h1bkg.SetName("h1bkg");
    h1sig = get_bkg_subtracted(h1m_uls_same, h1bkg);
    h1sig.SetName("h1sig");
    h1bkg_tmp = h1bkg.Clone("h1bkg_tmp"); #only for significance
    h1sig_tmp = h1sig.Clone("h1sig_tmp"); #only for significance


    h1m_uls_same.Scale(1./nev);
    h1m_uls_same.Scale(1., "width");
    h1bkg.Scale(1./nev);
    h1bkg.Scale(1., "width");
    h1sig.Scale(1./nev);
    h1sig.Scale(1., "width");
    h1sig.SetYTitle("#frac{1}{N_{ev}} #frac{dN}{dm_{ee}} (GeV/c^{2})^{-1}");

    f1total = TF1("f1total", "crystalball(0) + crystalball(5) + expo(10)", 0.4, 1.2);
    f1total.SetNpx(1000);
    f1total.SetParameter(0, 1e-6);
    f1total.SetParameter(1, 0.782);
    f1total.SetParameter(2, 0.02);
    f1total.FixParameter(3, 0.3);
    f1total.FixParameter(4, 3);

    f1total.SetParameter(5, 1e-6);
    f1total.SetParameter(6, 1.019);
    f1total.SetParameter(7, 0.02);
    f1total.FixParameter(8, 0.4);
    f1total.FixParameter(9, 4);

    f1total.SetParameter(10, -11);
    f1total.SetParameter(11, -3);
    f1total.SetParLimits(0, 1e-6, 1e-4);
    f1total.SetParLimits(1, 0.7, 0.8);
    f1total.SetParLimits(5, 1e-6, 1e-4);
    f1total.SetParLimits(6, 0.95, 1.05);
    #h1sig.Fit(f1total, "SMEI", "", 0.5, 1.1);

    h1sb = get_ratio(h1sig, h1bkg);
    h1sb.SetName("h1sbratio");
    h1sb.SetYTitle("S/B ratio");

    h1significance = get_significance(h1sig_tmp, h1bkg_tmp);
    h1significance.SetName("h1significance");
    h1significance.SetYTitle("significance S/#sqrt{S + 2B}");

    h1eval = h1sb.Clone("h1eval");
    h1eval.Reset();
    h1eval.Multiply(h1sb, h1significance, 1., 1., "B");
    h1eval.SetYTitle("S/B #times S/#sqrt{S + 2B}");

    outlist.Add(h1bkg);
    outlist.Add(h1sig);
    outlist.Add(h1sb);
    outlist.Add(h1significance);
    outlist.Add(h1eval);
    outlist.Add(f1total);

    hs_uls_same .GetAxis(0).SetRange(0, -1);
    hs_lspp_same.GetAxis(0).SetRange(0, -1);
    hs_lsmm_same.GetAxis(0).SetRange(0, -1);
    hs_uls_same .GetAxis(1).SetRange(0, -1);
    hs_lspp_same.GetAxis(1).SetRange(0, -1);
    hs_lsmm_same.GetAxis(1).SetRange(0, -1);
    hs_uls_same .GetAxis(2).SetRange(0, -1);
    hs_lspp_same.GetAxis(2).SetRange(0, -1);
    hs_lsmm_same.GetAxis(2).SetRange(0, -1);
    hs_uls_same .GetAxis(3).SetRange(0, -1);
    hs_lspp_same.GetAxis(3).SetRange(0, -1);
    hs_lsmm_same.GetAxis(3).SetRange(0, -1);
    return outlist;

#____________________________________________________________________________________________
def get_mee_sliced_by_single_dca(hs_uls_same, hs_lspp_same, hs_lsmm_same, dca_min, dca_max, nev, arr_mee):
    outlist = THashList();
    outlist.SetName("mee_in_single_dca_{0:2.1f}_{1:2.1f}_sigma".format(dca_min, dca_max));
    outlist.SetOwner(True);

    bin1 = hs_uls_same.GetAxis(1).FindBin(dca_min + 1e-3);
    bin2 = hs_uls_same.GetAxis(1).FindBin(dca_max - 1e-3);

    hs_uls_same .GetAxis(1).SetRange(bin1, bin2);
    hs_lspp_same.GetAxis(1).SetRange(bin1, bin2);
    hs_lsmm_same.GetAxis(1).SetRange(bin1, bin2);
    hs_uls_same .GetAxis(2).SetRange(bin1, bin2);
    hs_lspp_same.GetAxis(2).SetRange(bin1, bin2);
    hs_lsmm_same.GetAxis(2).SetRange(bin1, bin2);

    #mass spectrum integrated over DCA
    h1m_uls_same_org  = hs_uls_same.Projection(0);
    h1m_lspp_same_org = hs_lspp_same.Projection(0);
    h1m_lsmm_same_org = hs_lsmm_same.Projection(0);
    h1m_uls_same_org .SetName("h1m_uls_same_org");
    h1m_lspp_same_org.SetName("h1m_lspp_same_org");
    h1m_lsmm_same_org.SetName("h1m_lsmm_same_org");
    h1m_uls_same  = rebin_histogram(h1m_uls_same_org, arr_mee, False, False);
    h1m_lspp_same = rebin_histogram(h1m_lspp_same_org, arr_mee, False, False);
    h1m_lsmm_same = rebin_histogram(h1m_lsmm_same_org, arr_mee, False, False);

    h1m_uls_same .SetName("h1m_uls_same");
    h1m_lspp_same.SetName("h1m_lspp_same");
    h1m_lsmm_same.SetName("h1m_lsmm_same");
    outlist.Add(h1m_uls_same );
    outlist.Add(h1m_lspp_same);
    outlist.Add(h1m_lsmm_same);

    #h1m_uls_mix_org = hs_dilepton_uls_dca_mix.Projection(0);
    #h1m_lspp_mix_org = hs_dilepton_lspp_dca_mix.Projection(0);
    #h1m_lsmm_mix_org = hs_dilepton_lsmm_dca_mix.Projection(0);
    #h1m_uls_mix  = rebin_histogram(h1m_uls_mix_org, arr_mee, False, False);
    #h1m_lspp_mix = rebin_histogram(h1m_lspp_mix_org, arr_mee, False, False);
    #h1m_lsmm_mix = rebin_histogram(h1m_lsmm_mix_org, arr_mee, False, False);
    #h1m_uls_mix .SetName("h1m_uls_mix");
    #h1m_lspp_mix.SetName("h1m_lspp_mix");
    #h1m_lsmm_mix.SetName("h1m_lsmm_mix");
    #outlist.Add(h1m_uls_mix );
    #outlist.Add(h1m_lspp_mix);
    #outlist.Add(h1m_lsmm_mix);
    #h1R = get_R_factor(h1m_uls_mix, None, h1m_lspp_mix, h1m_lsmm_mix);
    #h1R.SetName("h1R");
    #outlist.Add(h1R);
    #h1bkg = get_corrected_bkg(h1R, h1m_lspp_same, h1m_lsmm_same);

    h1bkg = get_corrected_bkg_simple(1.0, 0.0, h1m_lspp_same, h1m_lsmm_same);
    h1bkg.SetName("h1bkg");
    h1bkg_tmp = h1bkg.Clone("h1bkg_tmp"); #only for significance
    h1sig = get_bkg_subtracted(h1m_uls_same, h1bkg);
    h1sig.SetName("h1sig");
    h1sig.SetYTitle("#frac{1}{N_{ev}} #frac{dN}{dm_{ee}} (GeV/c^{2})^{-1}");
    h1sig_tmp = h1sig.Clone("h1sig_tmp"); #only for significance

    h1m_uls_same.Scale(1./nev);
    h1m_uls_same.Scale(1., "width");
    h1bkg.Scale(1./nev);
    h1bkg.Scale(1., "width");
    h1sig.Scale(1./nev);
    h1sig.Scale(1., "width");

    f1total = TF1("f1total", "crystalball(0) + crystalball(5) + expo(10)", 0.4, 1.2);
    f1total.SetNpx(1000);
    f1total.SetParameter(0, 1e-6);
    f1total.SetParameter(1, 0.782);
    f1total.SetParameter(2, 0.02);
    f1total.FixParameter(3, 0.3);
    f1total.FixParameter(4, 3);

    f1total.SetParameter(5, 1e-6);
    f1total.SetParameter(6, 1.019);
    f1total.SetParameter(7, 0.02);
    f1total.FixParameter(8, 0.4);
    f1total.FixParameter(9, 4);

    f1total.SetParameter(10, -11);
    f1total.SetParameter(11, -3);
    f1total.SetParLimits(0, 1e-6, 1e-4);
    f1total.SetParLimits(1, 0.7, 0.8);
    f1total.SetParLimits(5, 1e-6, 1e-4);
    f1total.SetParLimits(6, 0.95, 1.05);
    #h1sig.Fit(f1total, "SMEI", "", 0.5, 1.1);

    h1sb = get_ratio(h1sig, h1bkg);
    h1sb.SetName("h1sbratio");
    h1sb.SetYTitle("S/B ratio");

    h1significance = get_significance(h1sig_tmp, h1bkg_tmp);
    h1significance.SetName("h1significance");
    h1significance.SetYTitle("significance S/#sqrt{S + 2B}");

    h1eval = h1sb.Clone("h1eval");
    h1eval.Reset();
    h1eval.Multiply(h1sb, h1significance, 1., 1., "B");
    h1eval.SetYTitle("S/B #times S/#sqrt{S + 2B}");

    outlist.Add(h1bkg);
    outlist.Add(h1sig);
    outlist.Add(h1sb);
    outlist.Add(h1significance);
    outlist.Add(h1eval);
    outlist.Add(f1total);

    hs_uls_same .GetAxis(0).SetRange(0, -1);
    hs_lspp_same.GetAxis(0).SetRange(0, -1);
    hs_lsmm_same.GetAxis(0).SetRange(0, -1);
    hs_uls_same .GetAxis(1).SetRange(0, -1);
    hs_lspp_same.GetAxis(1).SetRange(0, -1);
    hs_lsmm_same.GetAxis(1).SetRange(0, -1);
    hs_uls_same .GetAxis(2).SetRange(0, -1);
    hs_lspp_same.GetAxis(2).SetRange(0, -1);
    hs_lsmm_same.GetAxis(2).SetRange(0, -1);
    hs_uls_same .GetAxis(3).SetRange(0, -1);
    hs_lspp_same.GetAxis(3).SetRange(0, -1);
    hs_lsmm_same.GetAxis(3).SetRange(0, -1);

    return outlist;
#____________________________________________________________________________________________
def get_pair_dca_sliced_by_mee(hs_uls_same, hs_lspp_same, hs_lsmm_same, m_min, m_max, nev):
    outlist = THashList();
    outlist.SetName("pair_dca_in_mee_{0:3.2f}_{1:3.2f}_GeV".format(m_min, m_max));

    bin1 = hs_uls_same.GetAxis(0).FindBin(m_min + 1e-3);
    bin2 = hs_uls_same.GetAxis(0).FindBin(m_max - 1e-3);

    hs_uls_same .GetAxis(0).SetRange(bin1, bin2);
    hs_lspp_same.GetAxis(0).SetRange(bin1, bin2);
    hs_lsmm_same.GetAxis(0).SetRange(bin1, bin2);

    #mass spectrum integrated over DCA
    h1m_uls_same  = hs_uls_same.Projection(3);
    h1m_lspp_same = hs_lspp_same.Projection(3);
    h1m_lsmm_same = hs_lsmm_same.Projection(3);

    h1m_uls_same .SetName("h1m_uls_same");
    h1m_lspp_same.SetName("h1m_lspp_same");
    h1m_lsmm_same.SetName("h1m_lsmm_same");
    outlist.Add(h1m_uls_same );
    outlist.Add(h1m_lspp_same);
    outlist.Add(h1m_lsmm_same);

    #h1m_uls_mix_org = hs_dilepton_uls_dca_mix.Projection(0);
    #h1m_lspp_mix_org = hs_dilepton_lspp_dca_mix.Projection(0);
    #h1m_lsmm_mix_org = hs_dilepton_lsmm_dca_mix.Projection(0);
    #h1m_uls_mix  = rebin_histogram(h1m_uls_mix_org, arr_mee, False, False);
    #h1m_lspp_mix = rebin_histogram(h1m_lspp_mix_org, arr_mee, False, False);
    #h1m_lsmm_mix = rebin_histogram(h1m_lsmm_mix_org, arr_mee, False, False);
    #h1m_uls_mix .SetName("h1m_uls_mix");
    #h1m_lspp_mix.SetName("h1m_lspp_mix");
    #h1m_lsmm_mix.SetName("h1m_lsmm_mix");
    #outlist.Add(h1m_uls_mix );
    #outlist.Add(h1m_lspp_mix);
    #outlist.Add(h1m_lsmm_mix);
    #h1R = get_R_factor(h1m_uls_mix, None, h1m_lspp_mix, h1m_lsmm_mix);
    #h1R.SetName("h1R");
    #outlist.Add(h1R);
    #h1bkg = get_corrected_bkg(h1R, h1m_lspp_same, h1m_lsmm_same);

    h1bkg = get_corrected_bkg_simple(1.0, 0.0, h1m_lspp_same, h1m_lsmm_same);
    h1bkg.SetName("h1bkg");
    h1sig = get_bkg_subtracted(h1m_uls_same, h1bkg);
    h1sig.SetName("h1sig");

    #h1m_uls_same.Scale(1./nev);
    h1m_uls_same.Scale(1., "width");
    #h1bkg.Scale(1./nev);
    h1bkg.Scale(1., "width");
    #h1sig.Scale(1./nev);
    nsig = h1sig.Integral(1, 999, "");
    h1sig.Scale(1/nsig);
    h1sig.Scale(1., "width");
    h1sig.SetYTitle("#frac{1}{N_{ee}} #frac{dN}{dDCA_{ee}^{3D}} (#sigma)^{-1}");
    #print("integral = ", h1sig.Integral(1, 999, "width"));

    outlist.Add(h1bkg);
    outlist.Add(h1sig);

    hs_uls_same .GetAxis(0).SetRange(0, -1);
    hs_lspp_same.GetAxis(0).SetRange(0, -1);
    hs_lsmm_same.GetAxis(0).SetRange(0, -1);
    hs_uls_same .GetAxis(1).SetRange(0, -1);
    hs_lspp_same.GetAxis(1).SetRange(0, -1);
    hs_lsmm_same.GetAxis(1).SetRange(0, -1);
    hs_uls_same .GetAxis(2).SetRange(0, -1);
    hs_lspp_same.GetAxis(2).SetRange(0, -1);
    hs_lsmm_same.GetAxis(2).SetRange(0, -1);
    hs_uls_same .GetAxis(3).SetRange(0, -1);
    hs_lspp_same.GetAxis(3).SetRange(0, -1);
    hs_lsmm_same.GetAxis(3).SetRange(0, -1);

    return outlist;
#____________________________________________________________________________________________
def get_single_dca_sliced_by_mee(hs_uls_same, hs_lspp_same, hs_lsmm_same, m_min, m_max, nev):
    outlist = THashList();
    outlist.SetName("single_dca_in_mee_{0:3.2f}_{1:3.2f}_GeV".format(m_min, m_max));

    bin1 = hs_uls_same.GetAxis(0).FindBin(m_min + 1e-3);
    bin2 = hs_uls_same.GetAxis(0).FindBin(m_max - 1e-3);

    hs_uls_same .GetAxis(0).SetRange(bin1, bin2);
    hs_lspp_same.GetAxis(0).SetRange(bin1, bin2);
    hs_lsmm_same.GetAxis(0).SetRange(bin1, bin2);

    #mass spectrum integrated over DCA
    h1m_uls_same  = hs_uls_same.Projection(2);#electron
    h1m_lspp_same = hs_lspp_same.Projection(2);#electron
    h1m_lsmm_same = hs_lsmm_same.Projection(2);#electron
    h1m_uls_same .SetName("h1m_uls_same");
    h1m_lspp_same.SetName("h1m_lspp_same");
    h1m_lsmm_same.SetName("h1m_lsmm_same");

    h1m_uls_same_pos  = hs_uls_same.Projection(1);#positron
    h1m_lspp_same_pos = hs_lspp_same.Projection(1);#positron
    h1m_lsmm_same_pos = hs_lsmm_same.Projection(1);#positron
    h1m_uls_same_pos .SetName("h1m_uls_same_pos");
    h1m_lspp_same_pos.SetName("h1m_lspp_same_pos");
    h1m_lsmm_same_pos.SetName("h1m_lsmm_same_pos");

    h1m_uls_same .Add(h1m_uls_same_pos , 1.);
    h1m_lspp_same.Add(h1m_lspp_same_pos, 1.);
    h1m_lsmm_same.Add(h1m_lsmm_same_pos, 1.);

    outlist.Add(h1m_uls_same );
    outlist.Add(h1m_lspp_same);
    outlist.Add(h1m_lsmm_same);

    #h1m_uls_mix_org = hs_dilepton_uls_dca_mix.Projection(0);
    #h1m_lspp_mix_org = hs_dilepton_lspp_dca_mix.Projection(0);
    #h1m_lsmm_mix_org = hs_dilepton_lsmm_dca_mix.Projection(0);
    #h1m_uls_mix  = rebin_histogram(h1m_uls_mix_org, arr_mee, False, False);
    #h1m_lspp_mix = rebin_histogram(h1m_lspp_mix_org, arr_mee, False, False);
    #h1m_lsmm_mix = rebin_histogram(h1m_lsmm_mix_org, arr_mee, False, False);
    #h1m_uls_mix .SetName("h1m_uls_mix");
    #h1m_lspp_mix.SetName("h1m_lspp_mix");
    #h1m_lsmm_mix.SetName("h1m_lsmm_mix");
    #outlist.Add(h1m_uls_mix );
    #outlist.Add(h1m_lspp_mix);
    #outlist.Add(h1m_lsmm_mix);
    #h1R = get_R_factor(h1m_uls_mix, None, h1m_lspp_mix, h1m_lsmm_mix);
    #h1R.SetName("h1R");
    #outlist.Add(h1R);
    #h1bkg = get_corrected_bkg(h1R, h1m_lspp_same, h1m_lsmm_same);

    h1bkg = get_corrected_bkg_simple(1.0, 0.0, h1m_lspp_same, h1m_lsmm_same);
    h1bkg.SetName("h1bkg");
    h1sig = get_bkg_subtracted(h1m_uls_same, h1bkg);
    h1sig.SetName("h1sig");

    #h1m_uls_same.Scale(1./nev);
    h1m_uls_same.Scale(1., "width");
    #h1bkg.Scale(1./nev);
    h1bkg.Scale(1., "width");
    #h1sig.Scale(1./nev);
    nsig = h1sig.Integral(1, 999, "");
    h1sig.Scale(1/nsig);
    h1sig.Scale(1., "width");
    h1sig.SetYTitle("#frac{1}{N_{e}} #frac{dN}{dDCA_{e}^{3D}} (#sigma)^{-1}");

    outlist.Add(h1bkg);
    outlist.Add(h1sig);

    hs_uls_same .GetAxis(0).SetRange(0, -1);
    hs_lspp_same.GetAxis(0).SetRange(0, -1);
    hs_lsmm_same.GetAxis(0).SetRange(0, -1);
    hs_uls_same .GetAxis(1).SetRange(0, -1);
    hs_lspp_same.GetAxis(1).SetRange(0, -1);
    hs_lsmm_same.GetAxis(1).SetRange(0, -1);
    hs_uls_same .GetAxis(2).SetRange(0, -1);
    hs_lspp_same.GetAxis(2).SetRange(0, -1);
    hs_lsmm_same.GetAxis(2).SetRange(0, -1);
    hs_uls_same .GetAxis(3).SetRange(0, -1);
    hs_lspp_same.GetAxis(3).SetRange(0, -1);
    hs_lsmm_same.GetAxis(3).SetRange(0, -1);

    return outlist;

#____________________________________________________________________________________________
def get_eff_correlation(list_h1m, nev, color, marker):
    m_min_prompt = 0.08;
    m_max_prompt = 0.14;
    m_min_nonprompt = 1.1;
    m_max_nonprompt = 2.5;

    mg = TMultiGraph();
    mg.SetName("mg_eff_corr");
    #mg.GetXaxis().SetTitle("1 - efficiency for nonprompt source in {0:3.2f} < m_{{ee}} < {1:3.2f} GeV/c^{{2}}".format(m_min_nonprompt, m_max_nonprompt));
    mg.GetXaxis().SetTitle("rejection power = 1/eff_{{nonprompt}} in {0:3.2f} < m_{{ee}} < {1:3.2f} GeV/c^{{2}}".format(m_min_nonprompt, m_max_nonprompt));
    mg.GetYaxis().SetTitle("eff_{{prompt}} in {0:3.2f} < m_{{ee}} < {1:3.2f} GeV/c^{{2}}".format(m_min_prompt, m_max_prompt));

    bin1_prompt = list_h1m[0].FindBin(m_min_prompt + 1e-3);
    bin2_prompt = list_h1m[0].FindBin(m_max_prompt - 1e-3);
    nall_prompt_err = ctypes.c_double(0);
    nall_prompt = list_h1m[0].IntegralAndError(bin1_prompt, bin2_prompt, nall_prompt_err, "width");

    bin1_nonprompt = list_h1m[0].FindBin(m_min_nonprompt + 1e-3);
    bin2_nonprompt = list_h1m[0].FindBin(m_max_nonprompt - 1e-3);
    nall_nonprompt_err = ctypes.c_double(0);
    nall_nonprompt = list_h1m[0].IntegralAndError(bin1_nonprompt, bin2_nonprompt, nall_nonprompt_err, "width");

    nd = len(list_h1m);
    for i in range(1, nd):
        bin1_prompt = list_h1m[i].FindBin(m_min_prompt + 1e-3);
        bin2_prompt = list_h1m[i].FindBin(m_max_prompt - 1e-3);
        n_prompt_err = ctypes.c_double(0);
        n_prompt = list_h1m[i].IntegralAndError(bin1_prompt, bin2_prompt, n_prompt_err, "width");
        eff_prompt = n_prompt/nall_prompt;
        eff_prompt_err = math.sqrt(eff_prompt * (1 - eff_prompt)/nall_prompt/nev);

        bin1_nonprompt = list_h1m[i].FindBin(m_min_nonprompt + 1e-3);
        bin2_nonprompt = list_h1m[i].FindBin(m_max_nonprompt - 1e-3);
        n_nonprompt_err = ctypes.c_double(0);
        n_nonprompt = list_h1m[i].IntegralAndError(bin1_nonprompt, bin2_nonprompt, n_nonprompt_err, "width");
        eff_nonprompt = n_nonprompt/nall_nonprompt;
        eff_nonprompt_err = math.sqrt(eff_nonprompt * (1 - eff_nonprompt)/nall_nonprompt/nev);
        rp = 1/eff_nonprompt;
        rp_err = rp * eff_nonprompt_err/eff_nonprompt;

        g1 = TGraphAsymmErrors();
        g1.SetName("g1_dca{0:d}".format(i));
        #g1.SetPoint(0, 1 - eff_nonprompt, eff_prompt);
        g1.SetPoint(0, rp, eff_prompt);
        g1.SetPointError(0, rp_err, rp_err, eff_prompt_err, eff_prompt_err);
        g1.SetMarkerStyle(marker);
        g1.SetMarkerSize(1.5);
        g1.SetLineWidth(2);
        g1.SetMarkerColor(color - i);
        g1.SetLineColor(color - i);
        mg.Add(g1);

    return mg;

#____________________________________________________________________________________________
if __name__ == "__main__":
    #filename = "AnalysisResults_HL_151847.root"; # LHC22o
    #filename = "AnalysisResults_HL_153470.root"; # LHC23zs
    #filename = "AnalysisResults_HL_163030.root"; # LHC23zz PbPb at 5.36 TeV
    filename = "AnalysisResults_HL_168431.root"; # LHC22o pass6
    cutname = "mee_all_tpchadrejortofreq";
    arr_mee = np.array([
            0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14,
            0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5,
            0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04, 1.06, 1.08, 1.1,
            1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
            2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
            2.8, 2.9, 2.95, 3.00, 3.05, 3.10, 3.15, 3.2, 3.3, 3.4, 3.5
    ], dtype=np.float64);
    #arr_ptee = np.array([0.0, 10.0], dtype=np.float64);
    analyze_mee_vs_single_dca3d(filename, cutname, arr_mee);


