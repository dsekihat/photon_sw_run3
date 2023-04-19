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
def draw_raw_yield_pi0(subsystem, cutname):

    rootfile_loose   = TFile.Open("output_data_ptspectrum_pp_13.6TeV_LHC22f_V0withLooseTrack.root","READ");
    rootfile_tpconly = TFile.Open("output_data_ptspectrum_pp_13.6TeV_LHC22f_V0withTPConlyTrack.root","READ");
    rootfile_lkb     = TFile.Open("output_data_ptspectrum_pp_13.6TeV_LHC22f_V0withLKB.root","READ");
    rootfile_run2    = TFile.Open("data_PCMResultsFullCorrection_PP.root","READ");

    list_pcm_pcm_loose   = rootfile_loose  .Get(subsystem);
    list_pcm_pcm_tpconly = rootfile_tpconly.Get(subsystem);
    list_pcm_pcm_lkb     = rootfile_lkb    .Get(subsystem);
    dir_pcm_pcm_run2     = rootfile_run2   .Get("Pi013TeV");

    list_pcm_pcm_loose_cut   = list_pcm_pcm_loose  .FindObject(cutname);
    list_pcm_pcm_tpconly_cut = list_pcm_pcm_tpconly.FindObject(cutname);
    list_pcm_pcm_lkb_cut     = list_pcm_pcm_lkb    .FindObject(cutname);

    h1_pcm_pcm_loose   = list_pcm_pcm_loose_cut.FindObject("h1yield");
    h1_pcm_pcm_tpconly = list_pcm_pcm_tpconly_cut.FindObject("h1yield");
    h1_pcm_pcm_lkb     = list_pcm_pcm_lkb_cut.FindObject("h1yield");
    h1_pcm_pcm_run2    = dir_pcm_pcm_run2 .Get("RAWYieldPerEventsPi0_INT7");
    h1_pcm_pcm_loose.SetDirectory(0);
    h1_pcm_pcm_tpconly.SetDirectory(0);
    h1_pcm_pcm_lkb .SetDirectory(0);
    h1_pcm_pcm_run2.SetDirectory(0);
    ROOT.SetOwnership(h1_pcm_pcm_loose, False);
    ROOT.SetOwnership(h1_pcm_pcm_tpconly, False);
    ROOT.SetOwnership(h1_pcm_pcm_lkb , False);
    ROOT.SetOwnership(h1_pcm_pcm_run2, False);
 
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetLogx(1);
    c1.SetLogy(1);
    c1.SetTicks(1,1);
    c1.SetMargin(0.15,0.02,0.1,0.02);
    frame1 = c1.DrawFrame(0.4,1e-8,20,1e-1);
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
    make_common_style(h1_pcm_pcm_loose  , 20, 1.0, kRed+1, 1, 0);

    h1_pcm_pcm_run2.Draw("E0same");
    h1_pcm_pcm_lkb .Draw("E0same");
    h1_pcm_pcm_tpconly.Draw("E0same");
    h1_pcm_pcm_loose.Draw("E0same");

    leg = TLegend(0.2,0.8,0.4,0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.03);
    leg.SetHeader("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem));
    leg.AddEntry(h1_pcm_pcm_run2   , "pp at #sqrt{#it{s}} = 13 TeV","LP");
    leg.AddEntry(h1_pcm_pcm_loose  , "pp at #sqrt{#it{s}} = 13.6 TeV, LHC22f pass3, Loose tracks","LP");
    leg.AddEntry(h1_pcm_pcm_tpconly, "pp at #sqrt{#it{s}} = 13.6 TeV, LHC22f pass3, TPConly tracks","LP");
    leg.AddEntry(h1_pcm_pcm_lkb    , "pp at #sqrt{#it{s}} = 13.6 TeV, LHC22f pass3, LK builder","LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_LHC22f_pass3_pi0_raw_yield_{1}_{2}.eps".format(date, subsystem, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_LHC22f_pass3_pi0_raw_yield_{1}_{2}.pdf".format(date, subsystem, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_LHC22f_pass3_pi0_raw_yield_{1}_{2}.png".format(date, subsystem, cutname));

    rootfile_loose.Close();
    rootfile_tpconly.Close();
    rootfile_lkb .Close();
    rootfile_run2.Close();

#_____________________________________________________________________
#_____________________________________________________________________
def draw_peak_mean_pi0():

    rootfile_myv0 = TFile.Open("output_data_ptspectrum_pp_13.6TeV_LHC22q_V0_with_dsekihat.root","READ");
    rootfile_lkb  = TFile.Open("output_data_ptspectrum_pp_13.6TeV_LHC22q_V0_with_LKB.root","READ");
    rootfile_run2  = TFile.Open("data_PCMResultsFullCorrection_PP.root","READ");

    list_pcm_pcm_myv0 = rootfile_myv0.Get("PCMPCM");
    list_pcm_pcm_lkb  = rootfile_lkb .Get("PCMPCM");
    dir_pcm_pcm_run2  = rootfile_run2.Get("Pi013TeV");
    dir_pcm_pcm_run2.ls();

    h1_pcm_pcm_myv0 = list_pcm_pcm_myv0.FindObject("h1mean");
    h1_pcm_pcm_lkb  = list_pcm_pcm_lkb .FindObject("h1mean");
    h1_pcm_pcm_run2 = dir_pcm_pcm_run2 .Get("Pi0_Mass_data_INT7");
    h1_pcm_pcm_myv0.SetDirectory(0);
    h1_pcm_pcm_lkb .SetDirectory(0);
    h1_pcm_pcm_run2.SetDirectory(0);
    ROOT.SetOwnership(h1_pcm_pcm_myv0, False);
    ROOT.SetOwnership(h1_pcm_pcm_lkb , False);
    ROOT.SetOwnership(h1_pcm_pcm_run2, False);
 
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetLogx(1);
    c1.SetTicks(1,1);
    c1.SetMargin(0.15,0.02,0.1,0.02);
    frame1 = c1.DrawFrame(0.4,0.12,20,0.15);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("peak position (GeV/#it{c}^{2})");
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetTitleOffset(1.3);
    frame1.GetYaxis().SetTitleOffset(2.0);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(True);

    make_common_style(h1_pcm_pcm_myv0, 20, 1.0, kRed+1, 1, 0);
    make_common_style(h1_pcm_pcm_lkb , 20, 1.0, kBlue+1, 1, 0);
    make_common_style(h1_pcm_pcm_run2, 20, 1.0, kBlack, 1, 0);

    h1_pcm_pcm_run2.Draw("E0same");
    h1_pcm_pcm_lkb .Draw("E0same");
    h1_pcm_pcm_myv0.Draw("E0same");

    leg = TLegend(0.2,0.8,0.4,0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.03);
    leg.SetHeader("#pi^{0} #rightarrow #gamma#gamma , PCM-PCM");
    leg.AddEntry(h1_pcm_pcm_run2, "pp at #sqrt{#it{s}} = 13 TeV","LP");
    leg.AddEntry(h1_pcm_pcm_myv0, "pp at #sqrt{#it{s}} = 13.6 TeV, LHC22q pass3, Daiki's V0 finder","LP");
    leg.AddEntry(h1_pcm_pcm_lkb, "pp at #sqrt{#it{s}} = 13.6 TeV, LHC22q pass3, LK builder","LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_LHC22q_pass3_pi0_peak_mean.eps".format(date));
    c1.SaveAs("{0}_pp_13.6TeV_LHC22q_pass3_pi0_peak_mean.pdf".format(date));
    c1.SaveAs("{0}_pp_13.6TeV_LHC22q_pass3_pi0_peak_mean.png".format(date));

    rootfile_myv0.Close();
    rootfile_lkb .Close();
    rootfile_run2.Close();

#_____________________________________________________________________
#_____________________________________________________________________
def draw_peak_sigma_pi0():

    rootfile_myv0 = TFile.Open("output_data_ptspectrum_pp_13.6TeV_LHC22q_V0_with_dsekihat.root","READ");
    rootfile_lkb  = TFile.Open("output_data_ptspectrum_pp_13.6TeV_LHC22q_V0_with_LKB.root","READ");
    rootfile_run2  = TFile.Open("data_PCMResultsFullCorrection_PP.root","READ");

    list_pcm_pcm_myv0 = rootfile_myv0.Get("PCMPCM");
    list_pcm_pcm_lkb  = rootfile_lkb .Get("PCMPCM");
    dir_pcm_pcm_run2  = rootfile_run2.Get("Pi013TeV");
    dir_pcm_pcm_run2.ls();

    h1_pcm_pcm_myv0 = list_pcm_pcm_myv0.FindObject("h1sigma");
    h1_pcm_pcm_lkb  = list_pcm_pcm_lkb .FindObject("h1sigma");
    h1_pcm_pcm_run2 = dir_pcm_pcm_run2 .Get("Pi0_Width_data_INT7");
    h1_pcm_pcm_myv0.SetDirectory(0);
    h1_pcm_pcm_lkb .SetDirectory(0);
    h1_pcm_pcm_run2.SetDirectory(0);
    ROOT.SetOwnership(h1_pcm_pcm_myv0, False);
    ROOT.SetOwnership(h1_pcm_pcm_lkb , False);
    ROOT.SetOwnership(h1_pcm_pcm_run2, False);
 
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetLogx(1);
    c1.SetTicks(1,1);
    c1.SetMargin(0.15,0.02,0.1,0.02);
    frame1 = c1.DrawFrame(0.4,0.,20,0.02);
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

    make_common_style(h1_pcm_pcm_myv0, 20, 1.0, kRed+1, 1, 0);
    make_common_style(h1_pcm_pcm_lkb , 20, 1.0, kBlue+1, 1, 0);
    make_common_style(h1_pcm_pcm_run2, 20, 1.0, kBlack, 1, 0);

    h1_pcm_pcm_run2.Draw("E0same");
    h1_pcm_pcm_lkb .Draw("E0same");
    h1_pcm_pcm_myv0.Draw("E0same");

    leg = TLegend(0.2,0.8,0.4,0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.03);
    leg.SetHeader("#pi^{0} #rightarrow #gamma#gamma , PCM-PCM");
    leg.AddEntry(h1_pcm_pcm_run2, "pp at #sqrt{#it{s}} = 13 TeV","LP");
    leg.AddEntry(h1_pcm_pcm_myv0, "pp at #sqrt{#it{s}} = 13.6 TeV, LHC22q pass3, Daiki's V0 finder","LP");
    leg.AddEntry(h1_pcm_pcm_lkb, "pp at #sqrt{#it{s}} = 13.6 TeV, LHC22q pass3, LK builder","LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_LHC22q_pass3_pi0_peak_sigma.eps".format(date));
    c1.SaveAs("{0}_pp_13.6TeV_LHC22q_pass3_pi0_peak_sigma.pdf".format(date));
    c1.SaveAs("{0}_pp_13.6TeV_LHC22q_pass3_pi0_peak_sigma.png".format(date));

    rootfile_myv0.Close();
    rootfile_lkb .Close();
    rootfile_run2.Close();

#_____________________________________________________________________
#_____________________________________________________________________
if __name__ == "__main__":
    #draw_raw_yield_pi0("PCMPCM", "qc_qc");
    draw_raw_yield_pi0("PCMPHOS", "qc_test03");
    #draw_peak_mean_pi0();
    #draw_peak_sigma_pi0();
