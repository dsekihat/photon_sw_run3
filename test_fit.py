import os, sys, shutil
import numpy as np
import math
import ctypes
import ROOT
from ROOT import TFile, TDirectory, THashList, TF1, TH1F
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue
from histo_manager import slice_histogram, rebin_histogram
from nm_fitter import NMFitter
from true_shape_fitter import TrueShapeFitter

#_____________________________________________________________________________________
def run(filename_data, filename_mc):
    rootfile_data = TFile.Open(filename_data, "READ");
    rootfile_mc   = TFile.Open(filename_mc  , "READ");

    ptmin = 1.0;
    ptmax = 100.0;

    outfile = TFile("20240209_signal_shape_fit_tmp.root", "RECREATE");

    rootdir_nm_data = rootfile_data.Get("pi0eta-to-gammagamma");
    rootdir_nm_mc   = rootfile_mc  .Get("pi0eta-to-gammagamma-mc");
    
    rootlist_pair_data = rootdir_nm_data.Get("Pair").FindObject("PCMPCM").FindObject("qc_qc").FindObject("nocut");
    rootlist_pair_mc   = rootdir_nm_mc  .Get("Pair").FindObject("PCMPCM").FindObject("qc_qc").FindObject("nocut");
    rootlist_pair_data.ls();
    rootlist_pair_mc  .ls();

    h2same_data = rootlist_pair_data.FindObject("hMggPt_Same");
    h2mix_data  = rootlist_pair_data.FindObject("hMggPt_Mixed");
    h2same_mc_primary = rootlist_pair_mc.FindObject("hMggPt_Pi0_Primary");
    h2same_mc_wd  = rootlist_pair_mc.FindObject("hMggPt_Pi0_FromWD");
    h2same_data.RebinX(2);
    h2mix_data .RebinX(2);
    h2same_mc_primary.RebinX(2);
    h2same_mc_wd .RebinX(2);
    #feed down correction can be done by tempalte fit? i.e. TFractionFit
    outfile.WriteTObject(h2same_data);
    outfile.WriteTObject(h2mix_data );
    outfile.WriteTObject(h2same_mc_primary);
    outfile.WriteTObject(h2same_mc_wd);

    #for data
    h1same_data = slice_histogram(h2same_data, ptmin, ptmax, "X", False);
    h1same_data.SetName("{0}_pt{1:d}".format(h2same_data.GetName().replace("h2", "h1"),0));
    h1mix_data = slice_histogram(h2mix_data, ptmin, ptmax, "X", False);
    h1mix_data.SetName("{0}_pt{1:d}".format(h2mix_data.GetName().replace("h2", "h1"),0));

    nmf = NMFitter(h1same_data, h1mix_data, "cb", "pol1");
    nmf.set_parameters(0.135, 0.005, 0.6, 12, True);
    [fitresult, h1sig, h1bkg, h1ratio, f1sig, f1bkg, f1total] = nmf.fit("SME", "", 0.04, 0.24);

    #for MC
    h1same_mc_primary = slice_histogram(h2same_mc_primary, ptmin, ptmax, "X", False);
    h1same_mc_primary.SetName("{0}_pt{1:d}".format(h2same_mc_primary.GetName().replace("h2", "h1"),0));
    h1same_mc_wd = slice_histogram(h2same_mc_wd, ptmin, ptmax, "X", False);
    h1same_mc_wd.SetName("{0}_pt{1:d}".format(h2same_mc_wd.GetName().replace("h2", "h1"),0));

    list_h1true = [h1same_mc_primary, h1same_mc_wd];
    tsf = TrueShapeFitter(list_h1true);
    f1true = TF1("f1true", tsf, 0, 1, len(list_h1true));
    f1true.SetNpx(1000);
    f1true.SetLineColor(kBlue+1);
    f1true.SetParameters(10, 1);
    h1sig.Fit(f1true, "SMEI", "", 0.04, 0.24);

    outfile.WriteTObject(h1same_data);
    outfile.WriteTObject(h1mix_data );
    outfile.WriteTObject(h1ratio);
    outfile.WriteTObject(h1bkg);
    outfile.WriteTObject(h1sig);
    outfile.WriteTObject(f1sig);
    outfile.WriteTObject(f1bkg);
    outfile.WriteTObject(f1total);
    outfile.WriteTObject(f1true);
    outfile.WriteTObject(h1same_mc_primary);
    outfile.WriteTObject(h1same_mc_wd);

    rootfile_data.Close();
    rootfile_mc  .Close();
    outfile.Close();
#_____________________________________________________________________________________
if __name__ == "__main__":
    filename_data = "AnalysisResults_HL_163030.root";
    filename_mc   = "AnalysisResults_HL_162145.root";
    run(filename_data, filename_mc);
