import os, sys, shutil
import math
import argparse
import numpy as np
import ctypes
import yaml
import ROOT
ROOT.gROOT.SetBatch(False);
from ROOT import TFile, THashList, TF1, TMath

filename = "AnalysisResults_HBT_20230429.root";
rootfile = TFile.Open(filename, "READ");
rootfile.ls();
rootdir = rootfile.Get("photon-hbt");
rootdir.ls();
list_pair = rootdir.Get("Pair");
list_pair.Print();
list_pair_ss = list_pair.FindObject("PCMPCM");
list_pair_ss.Print();
list_pair_ss_cut = list_pair_ss.FindObject("qc_qc");
list_pair_ss_cut.Print();

hs_q_same = list_pair_ss_cut.FindObject("hs_q_same");
hs_q_mix = list_pair_ss_cut.FindObject("hs_q_mix");
hs_q_same.Sumw2();
hs_q_mix .Sumw2();

ktmin = 0.4;
ktmax = 0.6;
bin0 = hs_q_same.GetAxis(4).FindBin(ktmin + 1e-3);
bin1 = hs_q_same.GetAxis(4).FindBin(ktmax - 1e-3);
print(bin0, bin1);

#hs_q_same.GetAxis(1).SetRangeUser(-0.03, +0.03);#qout
#hs_q_mix .GetAxis(1).SetRangeUser(-0.03, +0.03);#qout
##hs_q_same.GetAxis(2).SetRangeUser(-0.03, +0.03);#qlong
##hs_q_mix .GetAxis(2).SetRangeUser(-0.03, +0.03);#qlong
#hs_q_same.GetAxis(3).SetRangeUser(-0.03, +0.03);#qside
#hs_q_mix .GetAxis(3).SetRangeUser(-0.03, +0.03);#qside

hs_q_same.GetAxis(4).SetRange(bin0, bin1);#kT
hs_q_mix .GetAxis(4).SetRange(bin0, bin1);#kT
h1qinv_same = hs_q_same.Projection(0);
h1qinv_mix  = hs_q_mix .Projection(0);
npair_same = h1qinv_same.GetEntries();
npair_mix  = h1qinv_mix .GetEntries();

h1qinv_same.Scale(1/npair_same);
h1qinv_mix .Scale(1/npair_mix);
#h1qinv_same.Draw();

h1cf = h1qinv_same.Clone("h1cf");
h1cf.Reset();
h1cf.Divide(h1qinv_same, h1qinv_mix, 1., 1., "G");
h1cf.Draw();

ROOT.SetOwnership(h1qinv_same, False);
ROOT.SetOwnership(h1qinv_mix, False);
ROOT.SetOwnership(h1cf, False);
h1qinv_same.SetDirectory(0);
h1qinv_mix .SetDirectory(0);
h1cf .SetDirectory(0);

sf = TMath.Hbar() * TMath.C() / TMath.Qe() /1e+9 / 1e-15; #GeV x fm

f1 = TF1("f1","1 + [0] * exp(-[1]*[1]*x*x /0.197/0.197)",0,0.1);
f1.SetNpx(1000);
f1.SetParameters(0.5,1);
#f1.SetParNames("#lambda_{inv}","R_{inv} (fm)");
h1cf.Fit(f1,"SME","",0, 0.1);
ROOT.SetOwnership(f1, False);
rootfile.Close();
