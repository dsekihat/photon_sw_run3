import os, sys, shutil
import math
import argparse
import numpy as np
import ctypes
import yaml
import ROOT
ROOT.gROOT.SetBatch(True);
from ROOT import TFile, THashList, TF1
from histo_manager import slice_histogram, rebin_histogram, get_bkg_subtracted, get_ratio, get_R_factor, get_corrected_bkg, get_corrected_bkg_simple, get_significance
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue

#_______________________________________________________________________
def analyze_conditional_acceptance(filename):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("tagging-pi0-mc");
    rootdire.ls();

    list_photon = rootdire.Get("PCM");
    list_pair   = rootdire.Get("Pair").FindObject("PCMDalitz");

    list_photon_cut = list_photon.FindObject("qc");
    list_pair_cut = list_pair.FindObject("qc_mee_0_120_tpconly").FindObject("nocut");
    list_photon_cut.ls();
    list_pair_cut.ls();

    outfile = TFile("tmp_conditional_acceptance_pp_13.6TeV.root", "RECREATE");

    arr_pt = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 6.0, 10], dtype=np.float64);

    h1pt_denominator_org = list_photon_cut.FindObject("hPt_v0photon_Pi0_Primary");
    h1pt_denominator = rebin_histogram(h1pt_denominator_org, arr_pt, False, False);
    h1pt_denominator.SetName("h1pt_denominator");

    h2m_pt_org = list_pair_cut.FindObject("hMggPt_Pi0_Primary");
    h1pt_numerator_org = slice_histogram(h2m_pt_org, 0.0, 1.0, "Y", False);
    h1pt_numerator = rebin_histogram(h1pt_numerator_org, arr_pt, False, False);
    h1pt_numerator.SetName("h1pt_numerator");

    outfile.WriteTObject(h1pt_denominator_org);
    outfile.WriteTObject(h1pt_numerator_org);
    outfile.WriteTObject(h1pt_denominator);
    outfile.WriteTObject(h1pt_numerator);

    h1eff = h1pt_numerator.Clone("h1eff");
    h1eff.Sumw2();
    h1eff.Reset();
    h1eff.Divide(h1pt_numerator, h1pt_denominator, 1., 1., "B");
    outfile.WriteTObject(h1eff);

    outfile.Close();
    rootfile.Close();
#_______________________________________________________________________
def analyze_single_photon(filename):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("pcm-qc");
    rootdire.ls();

    outfile = TFile("tmp_pcm_photon_pt_pp_13.6TeV.root", "RECREATE");

    list_ev = rootdire.Get("Event");
    list_ev.ls();
    h1z = list_ev.FindObject("hZvtx_after");
    nev = h1z.GetEntries();
    print(nev);

    list_v0_cut = rootdire.Get("V0").FindObject("qc");
    list_v0_cut.ls();
    h1pt = list_v0_cut.FindObject("hPt");
    h1pt.SetName("h1pt");
    h1pt.RebinX(100);
    h1pt.Scale(1/nev);
    nev_exp = 1e+12 * 60e-3; # 1pb^{-1} x 60 mb

    h1pt.Scale(nev_exp);
    h1pt.SetYTitle("expected raw counts in 1 pb^{-1}");
    #h1pt.SetYTitle("raw counts / N_{ev}");
    h1pt.SetDirectory(0);

    outfile.WriteTObject(h1z);
    outfile.WriteTObject(h1pt);
    outfile.Close();

    rootfile.Close();

#_______________________________________________________________________
if __name__ == "__main__":
    #filename = "AnalysisResults_HL_155759.root";
    #analyze_conditional_acceptance(filename);

    filename = "AnalysisResults_HL_155758.root";
    analyze_single_photon(filename);
