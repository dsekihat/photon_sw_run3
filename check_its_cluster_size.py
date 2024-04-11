import os, sys, shutil
import math
import argparse
import numpy as np
import ctypes
import yaml
import datetime
import ROOT
from ROOT import TFile, THashList, TF1, TCanvas, TPad, TLegend, TPaveText
from ROOT import kWhite, kBlack
from ROOT import gStyle, gROOT, gSystem

#____________________________________________________________________________________________
def check(filename):
    print(sys._getframe().f_code.co_name);
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("pcm-qc");
    rootdire.ls();
    list_track = rootdire.Get("V0Leg");
    #list_track.ls();
    list_track_itstpc = list_track.FindObject("qc_ITSTPC");
    list_track_itsonly = list_track.FindObject("qc_ITSonly");

    h1_itstpc  = list_track_itstpc .FindObject("hMeanClusterSizeITS");
    n_track_itstpc = h1_itstpc.GetEntries();
    h1_itstpc.Scale(1/n_track_itstpc);
    h1_itsonly = list_track_itsonly.FindObject("hMeanClusterSizeITS");
    n_track_itsonly = h1_itsonly.GetEntries();
    h1_itsonly.Scale(1/n_track_itsonly);
    h1_itstpc .RebinX(5);
    h1_itsonly.RebinX(5);

    h1_itstpc.Draw("h");
    h1_itsonly.Draw("h,E1same");
    h1_itstpc .SetDirectory(0);
    h1_itsonly.SetDirectory(0);

    ROOT.SetOwnership(h1_itstpc, False);
    ROOT.SetOwnership(h1_itsonly, False);

    #rootfile.Close();

#____________________________________________________________________________________________
if __name__ == "__main__":
    #filename = "AnalysisResults_HL_163030.root"; # LHC23zz PbPb at 5.36 TeV
    filename = "AnalysisResults_HL_163031.root"; # LHC23zs pass2 pp at 13.6 TeV
    check(filename);
