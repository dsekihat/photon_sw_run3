import os, sys, shutil
import math
import argparse
import numpy as np
import ctypes
import datetime
import yaml
import ROOT
ROOT.gROOT.SetBatch(False);
from ROOT import TFile, THashList, TF1, TMath, TCanvas, TPaveText
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

#______________________________________________________________
def hbt_qinv(ssname, cutname, kt1, kt2, period):
    filename = "AnalysisResults_HL_82418.root";
    rootfile = TFile.Open(filename, "READ");
    rootfile.ls();
    rootdir = rootfile.Get("photon-hbt");
    rootdir.ls();
    list_pair = rootdir.Get("Pair");
    list_pair.Print();
    list_pair_ss = list_pair.FindObject(ssname);
    list_pair_ss.Print();
    list_pair_ss_cut = list_pair_ss.FindObject(cutname);
    list_pair_ss_cut.Print();
    
    hs_q_same = list_pair_ss_cut.FindObject("hs_q_same");
    hs_q_mix = list_pair_ss_cut.FindObject("hs_q_mix");
    hs_q_same.Sumw2();
    hs_q_mix .Sumw2();
    
    ktmin = kt1;
    ktmax = kt2;
    bin0 = hs_q_same.GetAxis(4).FindBin(ktmin + 1e-3);
    bin1 = hs_q_same.GetAxis(4).FindBin(ktmax - 1e-3);
    print(bin0, bin1);
    
    #hs_q_same.GetAxis(1).SetRangeUser(-0.03, +0.03);#qlong
    #hs_q_mix .GetAxis(1).SetRangeUser(-0.03, +0.03);#qlong
    #hs_q_same.GetAxis(2).SetRangeUser(-0.03, +0.03);#qout
    #hs_q_mix .GetAxis(2).SetRangeUser(-0.03, +0.03);#qout
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
    h1cf.SetYTitle("C(q)");
    h1cf.Divide(h1qinv_same, h1qinv_mix, 1., 1., "G");
    
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
    f1.SetParNames("#lambda_{inv}","R_{inv} (fm)");
    h1cf.Fit(f1,"SME","",0, 0.06);
    ROOT.SetOwnership(f1, False);

    la = f1.GetParameter(0);
    la_err = f1.GetParError(0);
    R = f1.GetParameter(1);
    R_err = f1.GetParError(1);

    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetTicks(1,1);
    c1.SetMargin(0.12,0.05,0.1,0.05);

    frame1 = c1.DrawFrame(0.0, 0.5, 0.3, 5.0);
    frame1.GetXaxis().SetTitle("q_{inv} (GeV/c)");
    frame1.GetYaxis().SetTitle("C(q)");
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.5);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(True);
    frame1.GetYaxis().SetMaxDigits(3);

    txt = TPaveText(0.15,0.75,0.4,0.9,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0:3.2f} < #it{{k}}_{{T}} < {1:3.2f} GeV/#it{{c}}".format(kt1, kt2));
    txt.AddText(ssname);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    txt_fit = TPaveText(0.5,0.75,0.8,0.9,"NDC");
    txt_fit.SetFillColor(kWhite);
    txt_fit.SetFillStyle(0);
    txt_fit.SetBorderSize(0);
    txt_fit.SetTextAlign(12);#middle,left
    txt_fit.SetTextFont(42);#helvetica
    txt_fit.SetTextSize(0.035);
    txt_fit.AddText("C(q) = 1 + #lambda_{inv}exp(#minus R_{inv}^{2} q_{inv}^{2})");
    txt_fit.AddText("#lambda_{{inv}} = {0:3.2f} #pm {1:3.2f}".format(la, la_err));
    txt_fit.AddText("R_{{inv}} = {0:3.2f} #pm {1:3.2f} (fm)".format(R, R_err));
    txt_fit.Draw();
    ROOT.SetOwnership(txt,False);

    h1cf.SetMarkerStyle(20);
    h1cf.SetMarkerColor(kBlack);
    h1cf.SetLineColor(kBlack);
    h1cf.Draw("E0same");

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    #c1.SaveAs("{0}_pp_13.6TeV_{1}_pass3_{2}_{3}_Cq_kt{4:3.2f}_{5:3.2f}GeV.eps".format(date, period, ssname, cutname, kt1, kt2));
    #c1.SaveAs("{0}_pp_13.6TeV_{1}_pass3_{2}_{3}_Cq_kt{4:3.2f}_{5:3.2f}GeV.pdf".format(date, period, ssname, cutname, kt1, kt2));
    #c1.SaveAs("{0}_pp_13.6TeV_{1}_pass3_{2}_{3}_Cq_kt{4:3.2f}_{5:3.2f}GeV.png".format(date, period, ssname, cutname, kt1, kt2));

    rootfile.Close();
#______________________________________________________________
def hbt_2d_out_side():
    filename = "AnalysisResults_HL_82418.root";
    rootfile = TFile.Open(filename, "READ");
    rootfile.ls();
    rootdir = rootfile.Get("photon-hbt");
    rootdir.ls();
    list_pair = rootdir.Get("Pair");
    list_pair.Print();
    list_pair_ss = list_pair.FindObject("PCMPCM");
    list_pair_ss.Print();
    #list_pair_ss_cut = list_pair_ss.FindObject("qc_qc");
    list_pair_ss_cut = list_pair_ss.FindObject("analysis_wo_mee_analysis_wo_mee");
    #list_pair_ss_cut = list_pair_ss.FindObject("analysis_wo_mee_test03");
    list_pair_ss_cut.Print();
    
    hs_q_same = list_pair_ss_cut.FindObject("hs_q_same");
    hs_q_mix = list_pair_ss_cut.FindObject("hs_q_mix");
    hs_q_same.Sumw2();
    hs_q_mix .Sumw2();
    
    ktmin = 0.4;
    ktmax = 0.8;
    bin0 = hs_q_same.GetAxis(4).FindBin(ktmin + 1e-3);
    bin1 = hs_q_same.GetAxis(4).FindBin(ktmax - 1e-3);

    hs_q_same.GetAxis(1).SetRangeUser(-0.03, +0.03);#qlong
    hs_q_mix .GetAxis(1).SetRangeUser(-0.03, +0.03);#qlong
    #hs_q_same.GetAxis(2).SetRangeUser(-0.03, +0.03);#qout
    #hs_q_mix .GetAxis(2).SetRangeUser(-0.03, +0.03);#qout
    #hs_q_same.GetAxis(3).SetRangeUser(-0.03, +0.03);#qside
    #hs_q_mix .GetAxis(3).SetRangeUser(-0.03, +0.03);#qside
    
    hs_q_same.GetAxis(4).SetRange(bin0, bin1);#kT
    hs_q_mix .GetAxis(4).SetRange(bin0, bin1);#kT

    h2same = hs_q_same.Projection(3,2);
    h2mix  = hs_q_mix .Projection(3,2);

    h2same.SetDirectory(0);
    h2mix.SetDirectory(0);
    ROOT.SetOwnership(h2same, False);
    ROOT.SetOwnership(h2mix, False);

    npair_same = h2same.GetEntries();
    npair_mix  = h2mix .GetEntries();
    
    h2same.Scale(1/npair_same);
    h2mix .Scale(1/npair_mix);
    
    h2cf = h2same.Clone("h2cf");
    h2cf.Reset();
    h2cf.Divide(h2same, h2mix, 1., 1., "G");
    h2cf.Draw("colz");
    h2cf.SetDirectory(0);
    ROOT.SetOwnership(h2cf, False);

#    ROOT.SetOwnership(h1qinv_same, False);
#    ROOT.SetOwnership(h1qinv_mix, False);
#    ROOT.SetOwnership(h1cf, False);
#    h1qinv_same.SetDirectory(0);
#    h1qinv_mix .SetDirectory(0);
#    h1cf .SetDirectory(0);
#    
#    sf = TMath.Hbar() * TMath.C() / TMath.Qe() /1e+9 / 1e-15; #GeV x fm
#    
#    f1 = TF1("f1","1 + [0] * exp(-[1]*[1]*x*x /0.197/0.197)",0,0.1);
#    f1.SetNpx(1000);
#    f1.SetParameters(0.5,1);
#    f1.SetParNames("#lambda_{inv}","R_{inv} (fm)");
#    h1cf.Fit(f1,"SME","",0, 0.06);
#    ROOT.SetOwnership(f1, False);
#    rootfile.Close();


if __name__ == "__main__":
    period = "LHC22q";
    cutname = "analysis_wo_mee_analysis_wo_mee";
    hbt_qinv("PCMPCM", cutname, 0.4,0.6, period);
    #hbt_qinv("PCMPHOS", "analysis_wo_mee_test03", 0.8,1.0, period);
    #hbt_2d_out_side();
