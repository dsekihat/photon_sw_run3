import numpy as np
import datetime
import ROOT
from  file_manager import FileManager
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan
from histo_manager import slice_histogram, rebin_histogram
from painter import make_common_style
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

#_____________________________________________________________________
def draw_xy(filename, subsystem, period, cutname, suffix):
    rootfile = TFile.Open(filename,"READ");
    rootdir = rootfile.Get("pcm-qc");
    rootdir.ls();
    list_v0  = rootdir.Get("V0");
    list_v0.Print();
    list_cut = list_v0.FindObject(cutname);

    h2 = list_cut.FindObject("hGammaRxy_recalc");
    ROOT.SetOwnership(h2, False);
    h2.SetDirectory(0);
    h2.SetContour(100);
 
    gStyle.SetPalette(55);
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetTicks(1,1);
    c1.SetLogz(1);
    c1.SetMargin(0.12,0.12,0.1,0.05);

    frame1 = c1.DrawFrame(-100,-100,100,100);
    frame1.GetXaxis().SetTitle("conversion point X (cm)");
    frame1.GetYaxis().SetTitle("conversion point Y (cm)");
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

    h2.Draw("colz, same");

    #leg = TLegend(0.6,0.7,0.8,0.9);
    #leg.SetBorderSize(0);
    #leg.SetFillColor(kWhite);
    #leg.SetTextSize(0.035);
    #leg.SetHeader("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem));
    #leg.AddEntry(h1same , "same","LP");
    #leg.AddEntry(h1bkg  , "bkg (mixed event)","LP");
    #leg.AddEntry(h1sig  , "signal","LP");
    #leg.Draw("");
    #ROOT.SetOwnership(leg,False);

    #txt = TPaveText(0.15,0.7,0.4,0.9,"NDC");
    #txt.SetFillColor(kWhite);
    #txt.SetFillStyle(0);
    #txt.SetBorderSize(0);
    #txt.SetTextAlign(12);#middle,left
    #txt.SetTextFont(42);#helvetica
    #txt.SetTextSize(0.035);
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("{0} pass3".format(period));
    ##txt.AddText("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem));
    #txt.AddText(ptstr[25:]);
    #txt.AddText(suffix);
    #txt.Draw();
    #ROOT.SetOwnership(txt,False);


    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_rxy_{2}_{3}.eps".format(date, period, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_rxy_{2}_{3}.pdf".format(date, period, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_rxy_{2}_{3}.png".format(date, period, cutname, suffix));

    rootfile.Close();

#_____________________________________________________________________
def draw_mgg_pi0(filename, subsystem, period, cutname, pt_id=3, suffix=""):
    rootfile = TFile.Open(filename,"READ");
    list_ss  = rootfile.Get(subsystem);
    list_ss_cut = list_ss.FindObject(cutname);
    list_ss_cut_fit           = list_ss_cut  .FindObject("cb_pol1");
    list_ss_cut_fit_range     = list_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_ss_cut_fit_range_int = list_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");

    h1same = list_ss_cut_fit_range_int.FindObject("h1mgg_same_pt{0}".format(pt_id));
    h1bkg  = list_ss_cut_fit_range_int.FindObject("h1mgg_bkg_pt{0}".format(pt_id));
    h1sig  = list_ss_cut_fit_range_int.FindObject("h1mgg_sig_pt{0}".format(pt_id));
    h1same.SetDirectory(0);
    h1bkg.SetDirectory(0);
    h1sig.SetDirectory(0);
    ROOT.SetOwnership(h1same, False);
    ROOT.SetOwnership(h1bkg, False);
    ROOT.SetOwnership(h1sig, False);

    ptstr = h1same.GetTitle();
    #print(ptstr[25:]);
 
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetTicks(1,1);
    c1.SetMargin(0.12,0.05,0.1,0.05);

    h1same.GetXaxis().SetRangeUser(0.02,0.3);
    ymax = h1same.GetMaximum() * 1.4;
    ymin = h1same.GetMaximum() * -0.05;
    mw = h1same.GetBinWidth(1);

    frame1 = c1.DrawFrame(0.04,ymin,0.24,ymax);
    frame1.GetXaxis().SetTitle("#it{m}_{#gamma#gamma} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("Counts / {0:d} MeV/#it{{c}}^{{2}}".format(int(mw * 1e+3)));
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

    make_common_style(h1same  , 20, 1.0, kBlack, 1, 0);
    make_common_style(h1bkg   , 20, 1.0, kBlue+1, 1, 0);
    make_common_style(h1sig   , 20, 1.0, kRed+1, 1, 0);

    h1same.Draw("E0same");
    h1bkg. Draw("E0same");
    h1sig. Draw("E0same");

    leg = TLegend(0.6,0.7,0.8,0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.035);
    leg.SetHeader("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem));
    leg.AddEntry(h1same , "same","LP");
    leg.AddEntry(h1bkg  , "bkg (mixed event)","LP");
    leg.AddEntry(h1sig  , "signal","LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.15,0.7,0.4,0.9,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0} pass4".format(period));
    #txt.AddText("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem));
    txt.AddText(ptstr[25:]);
    txt.AddText(suffix);
    txt.Draw();
    ROOT.SetOwnership(txt,False);


    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_mgg_pi0_{2}_{3}_Pt{4}{5}.eps".format(date, subsystem, period, cutname, pt_id, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_mgg_pi0_{2}_{3}_Pt{4}{5}.pdf".format(date, subsystem, period, cutname, pt_id, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_mgg_pi0_{2}_{3}_Pt{4}{5}.png".format(date, subsystem, period, cutname, pt_id, suffix));

    rootfile.Close();

#_____________________________________________________________________
def draw_raw_yield_pi0(subsystem, period, cutname):
    rootfile_loose   = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withAnyTracks.root".format(period),"READ");
    rootfile_tpconly = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withTPConlyTracks.root".format(period),"READ");
    rootfile_tpconlybar = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withTPConlyTracksBar.root".format(period),"READ");
    rootfile_lkb     = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withLKB.root".format(period),"READ");
    rootfile_run2    = TFile.Open("data_PCMResultsFullCorrection_PP.root","READ");

    list_pcm_pcm_loose   = rootfile_loose  .Get(subsystem);
    list_pcm_pcm_tpconly = rootfile_tpconly.Get(subsystem);
    list_pcm_pcm_tpconlybar = rootfile_tpconlybar.Get(subsystem);
    list_pcm_pcm_lkb     = rootfile_lkb    .Get(subsystem);
    dir_pcm_pcm_run2     = rootfile_run2   .Get("Pi013TeV");

    list_pcm_pcm_loose_cut   = list_pcm_pcm_loose  .FindObject(cutname);
    list_pcm_pcm_tpconly_cut = list_pcm_pcm_tpconly.FindObject(cutname);
    list_pcm_pcm_tpconlybar_cut = list_pcm_pcm_tpconlybar.FindObject(cutname);
    list_pcm_pcm_lkb_cut     = list_pcm_pcm_lkb    .FindObject(cutname);

    list_pcm_pcm_loose_cut_fit   = list_pcm_pcm_loose_cut  .FindObject("cb_pol1");
    list_pcm_pcm_tpconly_cut_fit = list_pcm_pcm_tpconly_cut.FindObject("cb_pol1");
    list_pcm_pcm_tpconlybar_cut_fit = list_pcm_pcm_tpconlybar_cut.FindObject("cb_pol1");
    list_pcm_pcm_lkb_cut_fit     = list_pcm_pcm_lkb_cut    .FindObject("cb_pol1");

    list_pcm_pcm_loose_cut_fit_range   = list_pcm_pcm_loose_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_pcm_pcm_tpconly_cut_fit_range = list_pcm_pcm_tpconly_cut_fit.FindObject("fit_0.04_0.24_GeVc2");
    list_pcm_pcm_tpconlybar_cut_fit_range = list_pcm_pcm_tpconlybar_cut_fit.FindObject("fit_0.04_0.24_GeVc2");
    list_pcm_pcm_lkb_cut_fit_range     = list_pcm_pcm_lkb_cut_fit    .FindObject("fit_0.04_0.24_GeVc2");

    list_pcm_pcm_loose_cut_fit_range_int   = list_pcm_pcm_loose_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    list_pcm_pcm_tpconly_cut_fit_range_int = list_pcm_pcm_tpconly_cut_fit_range.FindObject("integral_-3.0_3.0_sigma");
    list_pcm_pcm_tpconlybar_cut_fit_range_int = list_pcm_pcm_tpconlybar_cut_fit_range.FindObject("integral_-3.0_3.0_sigma");
    list_pcm_pcm_lkb_cut_fit_range_int     = list_pcm_pcm_lkb_cut_fit_range    .FindObject("integral_-3.0_3.0_sigma");

    h1_pcm_pcm_loose   = list_pcm_pcm_loose_cut_fit_range_int.FindObject("h1yield");
    h1_pcm_pcm_tpconly = list_pcm_pcm_tpconly_cut_fit_range_int.FindObject("h1yield");
    h1_pcm_pcm_tpconlybar = list_pcm_pcm_tpconlybar_cut_fit_range_int.FindObject("h1yield");
    h1_pcm_pcm_lkb     = list_pcm_pcm_lkb_cut_fit_range_int.FindObject("h1yield");
    h1_pcm_pcm_run2    = dir_pcm_pcm_run2 .Get("RAWYieldPerEventsPi0_INT7");
    h1_pcm_pcm_loose.SetDirectory(0);
    h1_pcm_pcm_tpconly.SetDirectory(0);
    h1_pcm_pcm_tpconlybar.SetDirectory(0);
    h1_pcm_pcm_lkb .SetDirectory(0);
    h1_pcm_pcm_run2.SetDirectory(0);
    ROOT.SetOwnership(h1_pcm_pcm_loose, False);
    ROOT.SetOwnership(h1_pcm_pcm_tpconly, False);
    ROOT.SetOwnership(h1_pcm_pcm_tpconlybar, False);
    ROOT.SetOwnership(h1_pcm_pcm_lkb , False);
    ROOT.SetOwnership(h1_pcm_pcm_run2, False);
 
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetLogx(1);
    c1.SetLogy(1);
    c1.SetTicks(1,1);
    c1.SetMargin(0.15,0.02,0.1,0.02);
    frame1 = c1.DrawFrame(0.4,1e-8,20,1e-2);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}");
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetTitleOffset(1.3);
    frame1.GetYaxis().SetTitleOffset(2.0);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(True);

    make_common_style(h1_pcm_pcm_run2   , 20, 1.0, kBlack, 1, 0);
    make_common_style(h1_pcm_pcm_lkb    , 20, 1.0, kGreen+2, 1, 0);
    make_common_style(h1_pcm_pcm_tpconly, 20, 1.0, kBlue+1, 1, 0);
    make_common_style(h1_pcm_pcm_tpconlybar, 24, 1.0, kBlue+1, 1, 0);
    make_common_style(h1_pcm_pcm_loose  , 20, 1.0, kRed+1, 1, 0);

    h1_pcm_pcm_run2.Draw("E0same");
    h1_pcm_pcm_lkb .Draw("E0same");
    h1_pcm_pcm_tpconly.Draw("E0same");
    h1_pcm_pcm_tpconlybar.Draw("E0same");
    h1_pcm_pcm_loose.Draw("E0same");

    leg = TLegend(0.2,0.75,0.4,0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.03);
    leg.SetHeader("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem));
    leg.AddEntry(h1_pcm_pcm_run2   , "pp at #sqrt{#it{s}} = 13 TeV","LP");
    leg.AddEntry(h1_pcm_pcm_loose  , "pp at #sqrt{{#it{{s}}}} = 13.6 TeV, {0} pass4, Any tracks".format(period),"LP");
    leg.AddEntry(h1_pcm_pcm_tpconly, "pp at #sqrt{{#it{{s}}}} = 13.6 TeV, {0} pass4, TPConly tracks".format(period),"LP");
    leg.AddEntry(h1_pcm_pcm_tpconlybar, "pp at #sqrt{{#it{{s}}}} = 13.6 TeV, {0} pass4, Anti-TPConly tracks".format(period),"LP");
    leg.AddEntry(h1_pcm_pcm_lkb    , "pp at #sqrt{{#it{{s}}}} = 13.6 TeV, {0} pass4, LK builder".format(period),"LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_raw_yield_{2}_{3}.eps".format(date, subsystem, period, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_raw_yield_{2}_{3}.pdf".format(date, subsystem, period, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_raw_yield_{2}_{3}.png".format(date, subsystem, period, cutname));

    rootfile_loose.Close();
    rootfile_tpconly.Close();
    rootfile_tpconlybar.Close();
    rootfile_lkb .Close();
    rootfile_run2.Close();

#_____________________________________________________________________
def draw_peak_mean_pi0(subsystem, period, cutname):

    rootfile_loose   = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withAnyTracks.root".format(period),"READ");
    rootfile_tpconly = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withTPConlyTracks.root".format(period),"READ");
    rootfile_tpconlybar = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withTPConlyTracksBar.root".format(period),"READ");
    rootfile_lkb     = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withLKB.root".format(period),"READ");
    rootfile_run2    = TFile.Open("data_PCMResultsFullCorrection_PP.root","READ");

    list_pcm_pcm_loose   = rootfile_loose  .Get(subsystem);
    list_pcm_pcm_tpconly = rootfile_tpconly.Get(subsystem);
    list_pcm_pcm_tpconlybar = rootfile_tpconlybar.Get(subsystem);
    list_pcm_pcm_lkb     = rootfile_lkb    .Get(subsystem);
    dir_pcm_pcm_run2     = rootfile_run2   .Get("Pi013TeV");

    list_pcm_pcm_loose_cut   = list_pcm_pcm_loose  .FindObject(cutname);
    list_pcm_pcm_tpconly_cut = list_pcm_pcm_tpconly.FindObject(cutname);
    list_pcm_pcm_tpconlybar_cut = list_pcm_pcm_tpconlybar.FindObject(cutname);
    list_pcm_pcm_lkb_cut     = list_pcm_pcm_lkb    .FindObject(cutname);

    list_pcm_pcm_loose_cut_fit   = list_pcm_pcm_loose_cut  .FindObject("cb_pol1");
    list_pcm_pcm_tpconly_cut_fit = list_pcm_pcm_tpconly_cut.FindObject("cb_pol1");
    list_pcm_pcm_tpconlybar_cut_fit = list_pcm_pcm_tpconlybar_cut.FindObject("cb_pol1");
    list_pcm_pcm_lkb_cut_fit     = list_pcm_pcm_lkb_cut    .FindObject("cb_pol1");

    list_pcm_pcm_loose_cut_fit_range   = list_pcm_pcm_loose_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_pcm_pcm_tpconly_cut_fit_range = list_pcm_pcm_tpconly_cut_fit.FindObject("fit_0.04_0.24_GeVc2");
    list_pcm_pcm_tpconlybar_cut_fit_range = list_pcm_pcm_tpconlybar_cut_fit.FindObject("fit_0.04_0.24_GeVc2");
    list_pcm_pcm_lkb_cut_fit_range     = list_pcm_pcm_lkb_cut_fit    .FindObject("fit_0.04_0.24_GeVc2");

    list_pcm_pcm_loose_cut_fit_range_int   = list_pcm_pcm_loose_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    list_pcm_pcm_tpconly_cut_fit_range_int = list_pcm_pcm_tpconly_cut_fit_range.FindObject("integral_-3.0_3.0_sigma");
    list_pcm_pcm_tpconlybar_cut_fit_range_int = list_pcm_pcm_tpconlybar_cut_fit_range.FindObject("integral_-3.0_3.0_sigma");
    list_pcm_pcm_lkb_cut_fit_range_int     = list_pcm_pcm_lkb_cut_fit_range    .FindObject("integral_-3.0_3.0_sigma");


    h1_pcm_pcm_loose   = list_pcm_pcm_loose_cut_fit_range_int.FindObject("h1mean");
    h1_pcm_pcm_tpconly = list_pcm_pcm_tpconly_cut_fit_range_int.FindObject("h1mean");
    h1_pcm_pcm_tpconlybar = list_pcm_pcm_tpconlybar_cut_fit_range_int.FindObject("h1mean");
    h1_pcm_pcm_lkb     = list_pcm_pcm_lkb_cut_fit_range_int.FindObject("h1mean");
    h1_pcm_pcm_run2    = dir_pcm_pcm_run2 .Get("Pi0_Mass_data_INT7");
    h1_pcm_pcm_loose.SetDirectory(0);
    h1_pcm_pcm_tpconly.SetDirectory(0);
    h1_pcm_pcm_tpconlybar.SetDirectory(0);
    h1_pcm_pcm_lkb .SetDirectory(0);
    h1_pcm_pcm_run2.SetDirectory(0);
    ROOT.SetOwnership(h1_pcm_pcm_loose, False);
    ROOT.SetOwnership(h1_pcm_pcm_tpconly, False);
    ROOT.SetOwnership(h1_pcm_pcm_tpconlybar, False);
    ROOT.SetOwnership(h1_pcm_pcm_lkb , False);
    ROOT.SetOwnership(h1_pcm_pcm_run2, False);
 
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetLogx(1);
    c1.SetTicks(1,1);
    c1.SetMargin(0.15,0.02,0.1,0.02);
    frame1 = c1.DrawFrame(0.4,0.12,20, 0.15);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("peak mean (GeV/#it{c}^{2})");
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetTitleOffset(1.3);
    frame1.GetYaxis().SetTitleOffset(2.0);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(True);

    make_common_style(h1_pcm_pcm_run2   , 20, 1.0, kBlack, 1, 0);
    make_common_style(h1_pcm_pcm_lkb    , 20, 1.0, kGreen+2, 1, 0);
    make_common_style(h1_pcm_pcm_tpconly, 20, 1.0, kBlue+1, 1, 0);
    make_common_style(h1_pcm_pcm_tpconlybar, 24, 1.0, kBlue+1, 1, 0);
    make_common_style(h1_pcm_pcm_loose  , 20, 1.0, kRed+1, 1, 0);

    h1_pcm_pcm_run2.Draw("E0same");
    h1_pcm_pcm_lkb .Draw("E0same");
    h1_pcm_pcm_tpconly.Draw("E0same");
    h1_pcm_pcm_tpconlybar.Draw("E0same");
    h1_pcm_pcm_loose.Draw("E0same");

    leg = TLegend(0.2,0.75,0.4,0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.03);
    leg.SetHeader("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem));
    leg.AddEntry(h1_pcm_pcm_run2   , "pp at #sqrt{#it{s}} = 13 TeV","LP");
    leg.AddEntry(h1_pcm_pcm_loose  , "pp at #sqrt{{#it{{s}}}} = 13.6 TeV, {0} pass4, Any tracks".format(period),"LP");
    leg.AddEntry(h1_pcm_pcm_tpconly, "pp at #sqrt{{#it{{s}}}} = 13.6 TeV, {0} pass4, TPConly tracks".format(period),"LP");
    leg.AddEntry(h1_pcm_pcm_tpconlybar, "pp at #sqrt{{#it{{s}}}} = 13.6 TeV, {0} pass4, Anti-TPConly tracks".format(period),"LP");
    leg.AddEntry(h1_pcm_pcm_lkb    , "pp at #sqrt{{#it{{s}}}} = 13.6 TeV, {0} pass4, LK builder".format(period),"LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_peak_mean_{2}_{3}.eps".format(date, subsystem, period, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_peak_mean_{2}_{3}.pdf".format(date, subsystem, period, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_peak_mean_{2}_{3}.png".format(date, subsystem, period, cutname));

    rootfile_loose.Close();
    rootfile_tpconly.Close();
    rootfile_tpconlybar.Close();
    rootfile_lkb .Close();
    rootfile_run2.Close();

#_____________________________________________________________________
def draw_peak_sigma_pi0(subsystem, period, cutname):
    rootfile_loose   = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withAnyTracks.root".format(period),"READ");
    rootfile_tpconly = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withTPConlyTracks.root".format(period),"READ");
    rootfile_tpconlybar = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withTPConlyTracksBar.root".format(period),"READ");
    rootfile_lkb     = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withLKB.root".format(period),"READ");
    rootfile_run2    = TFile.Open("data_PCMResultsFullCorrection_PP.root","READ");

    list_pcm_pcm_loose   = rootfile_loose  .Get(subsystem);
    list_pcm_pcm_tpconly = rootfile_tpconly.Get(subsystem);
    list_pcm_pcm_tpconlybar = rootfile_tpconlybar.Get(subsystem);
    list_pcm_pcm_lkb     = rootfile_lkb    .Get(subsystem);
    dir_pcm_pcm_run2     = rootfile_run2   .Get("Pi013TeV");

    list_pcm_pcm_loose_cut   = list_pcm_pcm_loose  .FindObject(cutname);
    list_pcm_pcm_tpconly_cut = list_pcm_pcm_tpconly.FindObject(cutname);
    list_pcm_pcm_tpconlybar_cut = list_pcm_pcm_tpconlybar.FindObject(cutname);
    list_pcm_pcm_lkb_cut     = list_pcm_pcm_lkb    .FindObject(cutname);

    list_pcm_pcm_loose_cut_fit   = list_pcm_pcm_loose_cut  .FindObject("cb_pol1");
    list_pcm_pcm_tpconly_cut_fit = list_pcm_pcm_tpconly_cut.FindObject("cb_pol1");
    list_pcm_pcm_tpconlybar_cut_fit = list_pcm_pcm_tpconlybar_cut.FindObject("cb_pol1");
    list_pcm_pcm_lkb_cut_fit     = list_pcm_pcm_lkb_cut    .FindObject("cb_pol1");

    list_pcm_pcm_loose_cut_fit_range   = list_pcm_pcm_loose_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_pcm_pcm_tpconly_cut_fit_range = list_pcm_pcm_tpconly_cut_fit.FindObject("fit_0.04_0.24_GeVc2");
    list_pcm_pcm_tpconlybar_cut_fit_range = list_pcm_pcm_tpconlybar_cut_fit.FindObject("fit_0.04_0.24_GeVc2");
    list_pcm_pcm_lkb_cut_fit_range     = list_pcm_pcm_lkb_cut_fit    .FindObject("fit_0.04_0.24_GeVc2");

    list_pcm_pcm_loose_cut_fit_range_int   = list_pcm_pcm_loose_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    list_pcm_pcm_tpconly_cut_fit_range_int = list_pcm_pcm_tpconly_cut_fit_range.FindObject("integral_-3.0_3.0_sigma");
    list_pcm_pcm_tpconlybar_cut_fit_range_int = list_pcm_pcm_tpconlybar_cut_fit_range.FindObject("integral_-3.0_3.0_sigma");
    list_pcm_pcm_lkb_cut_fit_range_int     = list_pcm_pcm_lkb_cut_fit_range    .FindObject("integral_-3.0_3.0_sigma");


    h1_pcm_pcm_loose   = list_pcm_pcm_loose_cut_fit_range_int.FindObject("h1sigma");
    h1_pcm_pcm_tpconly = list_pcm_pcm_tpconly_cut_fit_range_int.FindObject("h1sigma");
    h1_pcm_pcm_tpconlybar = list_pcm_pcm_tpconlybar_cut_fit_range_int.FindObject("h1sigma");
    h1_pcm_pcm_lkb     = list_pcm_pcm_lkb_cut_fit_range_int.FindObject("h1sigma");
    h1_pcm_pcm_run2    = dir_pcm_pcm_run2 .Get("Pi0_Width_data_INT7");
    h1_pcm_pcm_loose.SetDirectory(0);
    h1_pcm_pcm_tpconly.SetDirectory(0);
    h1_pcm_pcm_tpconlybar.SetDirectory(0);
    h1_pcm_pcm_lkb .SetDirectory(0);
    h1_pcm_pcm_run2.SetDirectory(0);
    ROOT.SetOwnership(h1_pcm_pcm_loose, False);
    ROOT.SetOwnership(h1_pcm_pcm_tpconly, False);
    ROOT.SetOwnership(h1_pcm_pcm_tpconlybar, False);
    ROOT.SetOwnership(h1_pcm_pcm_lkb , False);
    ROOT.SetOwnership(h1_pcm_pcm_run2, False);
 
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetLogx(1);
    c1.SetTicks(1,1);
    c1.SetMargin(0.15,0.02,0.1,0.02);
    frame1 = c1.DrawFrame(0.4,0.,20, 0.02);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("peak sigma (GeV/#it{c}^{2})");
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetTitleOffset(1.3);
    frame1.GetYaxis().SetTitleOffset(2.0);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(True);

    make_common_style(h1_pcm_pcm_run2   , 20, 1.0, kBlack, 1, 0);
    make_common_style(h1_pcm_pcm_lkb    , 20, 1.0, kGreen+2, 1, 0);
    make_common_style(h1_pcm_pcm_tpconly, 20, 1.0, kBlue+1, 1, 0);
    make_common_style(h1_pcm_pcm_tpconlybar, 24, 1.0, kBlue+1, 1, 0);
    make_common_style(h1_pcm_pcm_loose  , 20, 1.0, kRed+1, 1, 0);

    h1_pcm_pcm_run2.Draw("E0same");
    h1_pcm_pcm_lkb .Draw("E0same");
    h1_pcm_pcm_tpconly.Draw("E0same");
    h1_pcm_pcm_tpconlybar.Draw("E0same");
    h1_pcm_pcm_loose.Draw("E0same");

    leg = TLegend(0.2,0.75,0.4,0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.03);
    leg.SetHeader("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem));
    leg.AddEntry(h1_pcm_pcm_run2   , "pp at #sqrt{#it{s}} = 13 TeV","LP");
    leg.AddEntry(h1_pcm_pcm_loose  , "pp at #sqrt{{#it{{s}}}} = 13.6 TeV, {0} pass4, Any tracks".format(period),"LP");
    leg.AddEntry(h1_pcm_pcm_tpconly, "pp at #sqrt{{#it{{s}}}} = 13.6 TeV, {0} pass4, TPConly tracks".format(period),"LP");
    leg.AddEntry(h1_pcm_pcm_tpconlybar, "pp at #sqrt{{#it{{s}}}} = 13.6 TeV, {0} pass4, Anti-TPConly tracks".format(period),"LP");
    leg.AddEntry(h1_pcm_pcm_lkb    , "pp at #sqrt{{#it{{s}}}} = 13.6 TeV, {0} pass4, LK builder".format(period),"LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_peak_sigma_{2}_{3}.eps".format(date, subsystem, period, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_peak_sigma_{2}_{3}.pdf".format(date, subsystem, period, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_peak_sigma_{2}_{3}.png".format(date, subsystem, period, cutname));

    rootfile_loose.Close();
    rootfile_tpconly.Close();
    rootfile_tpconlybar.Close();
    rootfile_lkb .Close();
    rootfile_run2.Close();

#_____________________________________________________________________
#_____________________________________________________________________
if __name__ == "__main__":
    draw_raw_yield_pi0( "PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee");
    draw_peak_mean_pi0( "PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee");
    draw_peak_sigma_pi0("PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee");
    draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withAnyTracks.root", "PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee", 1, "AnyTracks");
    draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withTPConlyTracks.root", "PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee", 1, "TPConlyTracks");
    draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withTPConlyTracksBar.root", "PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee", 1, "TPConlyTracksBar");
    draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withLKB.root", "PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee", 1, "LKB");
    #draw_xy("AnalysisResults_HL_78451.root", "PCMPCM", "LHC22m", "analysis_wo_mee", "AnyTracks");
    #draw_xy("AnalysisResults_HL_78449.root", "PCMPCM", "LHC22m", "analysis_wo_mee", "TPConlyTracks");
    #draw_xy("AnalysisResults_HL_78448.root", "PCMPCM", "LHC22m", "analysis_wo_mee", "TPConlyTracksBar");
    #draw_xy("AnalysisResults_HL_78450.root", "PCMPCM", "LHC22m", "analysis_wo_mee", "LKB");


    #draw_raw_yield_pi0( "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee");
    #draw_peak_mean_pi0( "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee");
    #draw_peak_sigma_pi0("PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee");
    #draw_mgg_pi0("output_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTracks.root", "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 3, "AnyTrack");
    #draw_mgg_pi0("output_data_ptspectrum_pp_13.6TeV_LHC22f_V0withTPConlyTrack.root", "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 3, "TPConlyTrack");
    #draw_mgg_pi0("output_data_ptspectrum_pp_13.6TeV_LHC22f_V0withLKB.root", "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 1, "LKB");

    #draw_xy("AnalysisResults_HL_75289.root", "PCMPCM", "LHC22q", "analysis_wo_mee", "AnyTrack");
    #draw_xy("AnalysisResults_HL_75290.root", "PCMPCM", "LHC22q", "analysis_wo_mee", "LKB");
    #draw_xy("AnalysisResults_HL_75291.root", "PCMPCM", "LHC22q", "analysis_wo_mee", "TPConlyTrack");

    #draw_raw_yield_pi0("PCMPCM", "LHC22q", "analysis_analysis");
    #draw_raw_yield_pi0("PCMPCM", "LHC22q", "qc_qc");

    #draw_raw_yield_pi0("PCMPHOS", "LHC22q", "qc_test03");
