import re
import numpy as np
import datetime
import math
import ROOT
import ctypes
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

#________________________________________________
def make_common_style(g1,marker,size,color,width=1,fill=0):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color);
    g1.SetMarkerSize(size);
    g1.SetLineColor(color);
    g1.SetLineWidth(width);
    g1.SetFillColor(color);
    g1.SetFillStyle(fill);
#________________________________________________
def check_material_phi(filename_data, filename_mc, cutname, etaid, rid, period_data, period_mc, suffix):
    rootfile_data = TFile.Open(filename_data, "READ");
    rootfile_mc   = TFile.Open(filename_mc  , "READ");

    list_data = rootfile_data.Get(cutname);
    list_mc = rootfile_mc.Get(cutname);

    h1data = list_data.FindObject("h1phi_eta{0:d}_r{1:d}".format(etaid, rid));
    h1mc = list_mc.FindObject("h1phi_eta{0:d}_r{1:d}".format(etaid, rid));
    h1data.SetDirectory(0);
    h1mc.SetDirectory(0);

    title = h1data.GetTitle();
    list_range = re.compile(r"[-+]?\d*\.\d+|\d+").findall(title);
    print(list_range);
    eta_min = float(list_range[0]);
    eta_max = float(list_range[1]);
    rxy_min = float(list_range[2]);
    rxy_max = float(list_range[3]);

    make_common_style(h1data, 20, 1.0, kRed+1, 1, 0);
    make_common_style(h1mc  , 20, 1.0, kBlue+1, 1, 0);
    ROOT.SetOwnership(h1data, False);
    ROOT.SetOwnership(h1mc, False);

    ymax = max(h1data.GetMaximum() , h1mc.GetMaximum()) * 1.7;
    ymin = max(h1data.GetMinimum() , h1mc.GetMinimum()) * -0.2;

    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.Divide(1,2,1e-3,1e-3);
    p1 = c1.cd(1);
    p1.SetPad(0,0.35,1,1);
    p1.SetMargin(0.15,0.02,0.0,0.06);
    p1.SetTicks(1,1);

    frame1 = p1.DrawFrame(0, ymin, TMath.TwoPi(), ymax);
    frame1.GetXaxis().SetTitle("conversion point #it{#varphi} (rad.)");
    frame1.GetYaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d^{3}#it{N}_{#gamma}}{d#it{r}_{xy} d#it{#eta} d#it{#varphi}} (cm #upoint rad.)^{#minus1}");
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetMaxDigits(3);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    ROOT.SetOwnership(frame1,False);

    h1mc.Draw("E0Hsame");
    h1data.Draw("E0Hsame");

    txt = TPaveText(0.2,0.72,0.4,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0:2.1f} < #it{{#eta}} < {1:2.1f}".format(eta_min, eta_max));
    txt.AddText("{0:2.1f} < #it{{r}}_{{xy}} < {1:2.1f} cm".format(rxy_min, rxy_max));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    leg = TLegend(0.2,0.6,0.4,0.7);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    #leg.AddEntry(h1data, "Data #gamma candidates (LHC22f pass4 520472)","LP");
    #leg.AddEntry(h1mc  , "M.C. rec. #gamma (LHC23d1f 520472)","LP");
    leg.AddEntry(h1data, "Data #gamma candidates (LHC22f pass4)","LP");
    leg.AddEntry(h1mc  , "M.C. rec. #gamma (LHC23d1f)","LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);


    p2 = c1.cd(2);
    p2.SetPad(0,0,1,0.35);
    p2.SetMargin(0.15,0.02,0.22,0.0);
    p2.SetTicks(1,1);

    frame2 = p2.DrawFrame(0,0.,TMath.TwoPi(),3.2);
    frame2.GetXaxis().SetTitle("#it{#varphi} (rad.)");
    frame2.GetYaxis().SetTitle("#frac{Data}{M.C.}");
    frame2.GetXaxis().SetTitleSize(0.10);
    frame2.GetYaxis().SetTitleSize(0.10);
    frame2.GetXaxis().SetTitleOffset(1.0);
    frame2.GetYaxis().SetTitleOffset(0.6);
    frame2.GetXaxis().SetLabelSize(0.1);
    frame2.GetYaxis().SetLabelSize(0.1);
    frame2.GetYaxis().CenterTitle(True);
    ROOT.SetOwnership(frame2,False);

    h1ratio = h1data.Clone("h1ratio");
    h1ratio.Reset();
    h1ratio.Sumw2();
    h1ratio.Divide(h1data, h1mc, 1., 1., "G");
    h1ratio.Draw("E0same");
    h1ratio.SetDirectory(0);
    ROOT.SetOwnership(h1ratio,False);

    line1 = TLine(0,1,TMath.TwoPi(),1);
    line1.SetLineColor(kBlack);
    line1.SetLineStyle(1);
    line1.SetLineWidth(2);
    line1.Draw("");
    ROOT.SetOwnership(line1,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_{2}_material_budget_vs_phi_eta{3}_r{4}_{5}{6}.eps".format(date, period_data, period_mc, etaid, rid, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_{2}_material_budget_vs_phi_eta{3}_r{4}_{5}{6}.pdf".format(date, period_data, period_mc, etaid, rid, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_{2}_material_budget_vs_phi_eta{3}_r{4}_{5}{6}.png".format(date, period_data, period_mc, etaid, rid, cutname, suffix));

#________________________________________________
if __name__ == "__main__":
    period_data = "LHC22f";
    period_mc   = "LHC23d1k";
    filename_data = "material_budget_data_pp_13.6TeV_{0}.root".format(period_data);
    filename_mc   = "material_budget_mc_pp_13.6TeV_{0}.root".format(period_mc);

    cutname = "wwire_ib";
    cutname = "qc";
    check_material_phi(filename_data, filename_mc, cutname,  0, 0, period_data, period_mc, ""); #3 ITS inner barrel layers should be visible.
    check_material_phi(filename_data, filename_mc, cutname,  0, 1, period_data, period_mc, ""); #Tungsate wire should be visible.
    check_material_phi(filename_data, filename_mc, cutname,  1, 0, period_data, period_mc, ""); #3 ITS inner barrel layers should be visible.
    check_material_phi(filename_data, filename_mc, cutname,  1, 1, period_data, period_mc, ""); #Tungsate wire should be visible.
    #check_material_phi(filename_data, filename_mc, cutname,  20, 12, period_data, period_mc, "");
    #check_material_phi(filename_data, filename_mc, cutname,  20, 14, period_data, period_mc, "");
    #check_material_phi(filename_data, filename_mc, cutname,  20, 15, period_data, period_mc, "");
    #check_material_phi(filename_data, filename_mc, cutname,  20, 8, period_data, period_mc, "");
    #check_material_phi(filename_data, filename_mc, cutname,  20, 5, period_data, period_mc, "");
    #check_material_phi(filename_data, filename_mc, cutname,  20, 3, period_data, period_mc, "");
    #check_material_phi(filename_data, filename_mc, cutname,  10, 8, period_data, period_mc, "");
#________________________________________________
