import numpy as np
import datetime
import math
import ROOT
import ctypes
from  file_manager import FileManager
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kGray
from histo_manager import slice_histogram, rebin_histogram, convert_dn2n, get_ratio
from painter import make_common_style
from bin_shift_correction import apply_bin_shift_correction_x
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

#_____________________________________________________________________
def draw_eta_to_pi0_ratio(filename_pi0, filename_eta, suffix=""):
    rootfile_pi0 = TFile.Open(filename_pi0, "READ");
    rootfile_eta = TFile.Open(filename_eta, "READ");

    h1xsection_pi0_org = rootfile_pi0.Get("h1iy");
    h1xsection_eta_org = rootfile_eta.Get("h1iy_eta");
    h1xsection_pi0_org.Scale(TMath.TwoPi());
    h1xsection_eta_org.Scale(TMath.TwoPi());
    h1xsection_pi0_org.Scale(1/59.4);
    h1xsection_eta_org.Scale(1/59.4);
    h1xsection_pi0_org.Scale(1/1.8);
    h1xsection_eta_org.Scale(1/1.8);

    h1n_pi0_org = convert_dn2n(h1xsection_pi0_org);
    h1n_eta = convert_dn2n(h1xsection_eta_org);

    arr_pt = np.array(h1n_eta.GetXaxis().GetXbins(), dtype=np.float64);
    print("arr_pt = ", arr_pt);

    h1n_pi0 = rebin_histogram(h1n_pi0_org, arr_pt, False, False);
    h1ratio = h1n_eta.Clone("h1ratio");
    h1ratio = get_ratio(h1n_eta, h1n_pi0, "G");
    make_common_style(h1ratio, 20, 1.0, kRed+1, 1, 0);
    h1ratio.SetDirectory(0);
    ROOT.SetOwnership(h1ratio,False);

    rootfile_run2 = TFile.Open("data_PCMResultsFullCorrection_PP.root","READ");
    dir_eta_run2  = rootfile_run2.Get("Eta13TeV");
    g1stat_run2    = dir_eta_run2.Get("graphEtaToPi0StatError");
    g1syst_run2    = dir_eta_run2.Get("EtaToPi0SystError");
    ROOT.SetOwnership(g1stat_run2, False);
    ROOT.SetOwnership(g1syst_run2, False);
    make_common_style(g1stat_run2, 24, 1.0, kGray+2, 1, 0);
    make_common_style(g1syst_run2, 24, 1.0, kGray+2, 1, 0);

    c1 = TCanvas("c0", "c0",0,0,800,800);
    c1.SetMargin(0.12, 0.03, 0.12, 0.03);
    c1.SetTicks(1,1);
   
    frame1 = c1.DrawFrame(0., 0, 10, 1.0);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#eta/#pi^{0}");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.5);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    g1stat_run2.Draw("PZ");
    g1syst_run2.Draw("P2Z");
    h1ratio.Draw("E0,same");

    leg = TLegend(0.6,0.8,0.8,0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(g1syst_run2, "pp at #sqrt{#it{s}} = 13 TeV", "P");
    leg.AddEntry(h1ratio, "pp at #sqrt{#it{s}} = 13.6 TeV", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    
    rootfile_pi0.Close();
    rootfile_eta.Close();
    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_eta_to_pi0_ratio{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_eta_to_pi0_ratio{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_eta_to_pi0_ratio{1}.png".format(date, suffix));

#_____________________________________________________________________
def draw_photon_conversion_rec_efficiency_rxy(filename, cutname, period, suffix):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("pcm-qc-mc");
    rootlist_gen = rootdire.Get("Generated");
    rootlist_gen.ls();
    rootlist_rec = rootdire.Get("V0").FindObject(cutname);
    rootlist_rec.ls();

    h2rz_rec = rootlist_rec.FindObject("hRZ_Photon_Primary_MC");
    h2rz_gen = rootlist_gen.FindObject("hPhotonRZ");
    h2rz_rec.Sumw2();
    h2rz_gen.Sumw2();
    h2rz_rec.SetDirectory(0);
    h2rz_gen.SetDirectory(0);
    ROOT.SetOwnership(h2rz_rec, False);
    ROOT.SetOwnership(h2rz_gen, False);
    h1z_rec = h2rz_rec.ProjectionY("h1z_rec");
    h1z_gen = h2rz_gen.ProjectionY("h1z_gen");
    h1z_gen.RebinX(2);

    c1 = TCanvas("c0", "c0",0,0,800,800);
    c1.SetMargin(0.12, 0.02, 0.12, 0.02);
    c1.SetTicks(1,1);
   
    frame1 = c1.DrawFrame(0.,0,90,0.2);
    frame1.GetXaxis().SetTitle("R_{xy} (cm)");
    frame1.GetYaxis().SetTitle("reconstruction efficiency");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.6);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1p = h1z_rec.Clone("h1p");
    h1p.Sumw2();
    h1p.Reset();
    h1p.Divide(h1z_rec, h1z_gen, 1., 1., "B");
    h1p.SetDirectory(0);
    ROOT.SetOwnership(h1p, False);
    h1p.Draw("E0same,h");

    make_common_style(h1p, 20, 1.0, kBlack, 1, 0);

    txt = TPaveText(0.6,0.7,0.8,0.95,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    txt.AddText("ALICE simulation");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("#gamma #rightarrow e^{+}e^{#minus}");
    txt.AddText("|#eta_{#gamma}| < 0.9");
    txt.AddText(period);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    rootfile.Close();
    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_photon_conversion_rec_eff_rxy_{2}{3}.eps".format(date, period, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_photon_conversion_rec_eff_rxy_{2}{3}.pdf".format(date, period, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_photon_conversion_rec_eff_rxy_{2}{3}.png".format(date, period, cutname, suffix));

#_____________________________________________________________________
def draw_photon_conversion_rec_efficiency(filename, cutname, period, suffix):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("pcm-qc-mc");
    rootlist_gen = rootdire.Get("Generated");
    #rootlist_gen.ls();
    rootlist_rec = rootdire.Get("V0").FindObject(cutname);
    rootlist_rec.ls();
    arr_pt = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6,7,8,9,10], dtype=float);

    h1photon_conv_rec_org = rootlist_rec.FindObject("hPt_Photon_Primary");
    h1photon_converted_org = rootlist_gen.FindObject("hPt_ConvertedPhoton");
    h1photon_conv_rec_org.Sumw2();
    h1photon_converted_org.Sumw2();
    h1photon_conv_rec_org.SetDirectory(0);
    h1photon_converted_org.SetDirectory(0);
    ROOT.SetOwnership(h1photon_conv_rec_org, False);
    ROOT.SetOwnership(h1photon_converted_org, False);

    h1photon_converted = rebin_histogram(h1photon_converted_org, arr_pt, False, False);
    h1photon_conv_rec = rebin_histogram(h1photon_conv_rec_org, arr_pt, False, False);
    h1photon_conv_rec.SetDirectory(0);
    h1photon_converted.SetDirectory(0);
    ROOT.SetOwnership(h1photon_conv_rec , False);
    ROOT.SetOwnership(h1photon_converted, False);

    c1 = TCanvas("c0", "c0",0,0,800,800);
    c1.SetMargin(0.12, 0.02, 0.12, 0.02);
    c1.SetTicks(1,1);
   
    frame1 = c1.DrawFrame(0.,0,10,0.2);
    frame1.GetXaxis().SetTitle("#it{p}_{T,#gamma} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("reconstruction efficiency");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.6);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1p = h1photon_conv_rec.Clone("h1p");
    h1p.Sumw2();
    h1p.Reset();
    h1p.Divide(h1photon_conv_rec, h1photon_converted, 1., 1., "B");
    h1p.SetDirectory(0);
    ROOT.SetOwnership(h1p, False);
    h1p.Draw("E0same");

    make_common_style(h1p, 20, 1.0, kBlack, 1, 0);

    txt = TPaveText(0.6,0.7,0.8,0.95,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    txt.AddText("ALICE simulation");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("#gamma #rightarrow e^{+}e^{#minus}");
    txt.AddText("|#eta_{#gamma}| < 0.9");
    txt.AddText(period);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    rootfile.Close();
    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_photon_conversion_rec_eff_{2}{3}.eps".format(date, period, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_photon_conversion_rec_eff_{2}{3}.pdf".format(date, period, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_photon_conversion_rec_eff_{2}{3}.png".format(date, period, cutname, suffix));

#_____________________________________________________________________
#_____________________________________________________________________
def draw_photon_conversion_probability(filename, period, suffix):
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("pcm-qc-mc");
    rootlist_gen = rootdire.Get("Generated");
    rootlist_gen.ls();
    arr_pt = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6,7,8,9,10], dtype=float);

    h1photon_all_org = rootlist_gen.FindObject("hPt_Photon");
    h1photon_converted_org = rootlist_gen.FindObject("hPt_ConvertedPhoton");
    h1photon_all_org.Sumw2();
    h1photon_converted_org.Sumw2();
    h1photon_all_org.SetDirectory(0);
    h1photon_converted_org.SetDirectory(0);
    ROOT.SetOwnership(h1photon_all_org, False);
    ROOT.SetOwnership(h1photon_converted_org, False);

    h1photon_converted = rebin_histogram(h1photon_converted_org, arr_pt, False, False);
    h1photon_all       = rebin_histogram(h1photon_all_org, arr_pt, False, False);
    h1photon_all      .SetDirectory(0);
    h1photon_converted.SetDirectory(0);
    ROOT.SetOwnership(h1photon_all      , False);
    ROOT.SetOwnership(h1photon_converted, False);

    c1 = TCanvas("c0", "c0",0,0,800,800);
    c1.SetMargin(0.12, 0.02, 0.12, 0.02);
    c1.SetTicks(1,1);
   
    frame1 = c1.DrawFrame(0.,0,10,0.2);
    frame1.GetXaxis().SetTitle("#it{p}_{T,#gamma} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("photon conversion probability");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.6);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1p = h1photon_all.Clone("h1p");
    h1p.Sumw2();
    h1p.Reset();
    h1p.Divide(h1photon_converted, h1photon_all, 1., 1., "B");
    h1p.SetDirectory(0);
    ROOT.SetOwnership(h1p, False);
    h1p.Draw("E0same");

    make_common_style(h1p, 20, 1.0, kBlack, 1, 0);

    txt = TPaveText(0.6,0.7,0.8,0.95,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    txt.AddText("ALICE simulation");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("#gamma #rightarrow e^{+}e^{#minus}");
    txt.AddText("|#eta_{#gamma}| < 0.9");
    txt.AddText(period);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    rootfile.Close();
    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_photon_conversion_probability{2}.eps".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_photon_conversion_probability{2}.pdf".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_photon_conversion_probability{2}.png".format(date, period, suffix));

#_____________________________________________________________________
def draw_photon_efficiency(filename, period, suffix):
    rootfile = TFile.Open(filename, "READ");
    rootlist_ss = rootfile.Get("PCM");
    cutnames = ["qc", "qc_ITSTPC", "qc_ITSonly"];
    nc = len(cutnames);
    list_h1eff = [];

    c1 = TCanvas("c0", "c0",0,0,800,800);
    c1.SetMargin(0.12, 0.03, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogx(1);
    c1.SetLogy(1);
   
    frame1 = c1.DrawFrame(0.1,1e-5,10,1e-2);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("acc. #times rec. efficiency");
    frame1.GetXaxis().SetTitleOffset(1.3);
    frame1.GetYaxis().SetTitleOffset(1.5);
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    for ic in range(0, nc):
        cutname = cutnames[ic];
        rootlist = rootlist_ss.FindObject(cutname);
        h1eff = rootlist.FindObject("h1eff");
        h1eff.SetDirectory(0);
        ROOT.SetOwnership(h1eff, False);
        h1eff.Draw("E0same");
        list_h1eff.append(h1eff);

    txt = TPaveText(0.15,0.84,0.4,0.94,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("LHC23d1k");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    make_common_style(list_h1eff[0], 20, 1.0, kBlack, 1, 0);
    make_common_style(list_h1eff[1], 20, 1.0, kRed+1, 1, 0);
    make_common_style(list_h1eff[2], 20, 1.0, kBlue+1, 1, 0);

    leg = TLegend(0.3,0.3,0.5,0.4);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(list_h1eff[0],"combined", "LP");
    leg.AddEntry(list_h1eff[1],"V0 with ITS-TPC matched tracks", "LP");
    leg.AddEntry(list_h1eff[2],"V0 with ITSonly tracks", "LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

 
    rootfile.Close();
    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_photon_efficiency{2}.eps".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_photon_efficiency{2}.pdf".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_photon_efficiency{2}.png".format(date, period, suffix));


#_____________________________________________________________________

#_____________________________________________________________________
def compare_nch_pv(filename_nominalB, filename_lowB, suffix):
    rootfile_nominalB = TFile.Open(filename_nominalB, "READ");
    rootdir_nominalB = rootfile_nominalB.Get("pcm-qc");
    rootlist_nominalB = rootdir_nominalB.Get("Event");
    h1nch_nominalB = rootlist_nominalB.FindObject("hMultNTracksPV");
    #rootlist_nominalB.ls();
    h1nch_nominalB.SetDirectory(0);
    ROOT.SetOwnership(h1nch_nominalB, False);
    nev1 = h1nch_nominalB.GetEntries();
    h1nch_nominalB.Scale(1/nev1);

    rootfile_lowB = TFile.Open(filename_lowB, "READ");
    rootdir_lowB = rootfile_lowB.Get("pcm-qc");
    rootlist_lowB = rootdir_lowB.Get("Event");
    h1nch_lowB = rootlist_lowB.FindObject("hMultNTracksPV");
    #rootlist_lowB.ls();
    h1nch_lowB.SetDirectory(0);
    ROOT.SetOwnership(h1nch_lowB, False);
    nev2 = h1nch_lowB.GetEntries();
    h1nch_lowB.Scale(1/nev2);

    c1 = TCanvas("c0", "c0",0,0,800,800);
    c1.SetMargin(0.12, 0.08, 0.1, 0.03);
    c1.SetTicks(1,1);
    c1.SetLogy(1);
   
    frame1 = c1.DrawFrame(-0.5,1e-11,1000.5,1);
    frame1.GetXaxis().SetTitle("N_{ch}^{PV} in |#eta| < 0.8");
    frame1.GetYaxis().SetTitle("probability density");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.6);
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMaxDigits(3);
    ROOT.SetOwnership(frame1,False);

    h1nch_lowB.Draw("E0same");
    h1nch_nominalB.Draw("E0same");
    make_common_style(h1nch_nominalB, 20, 1.0, kBlue+1, 1, 0);
    make_common_style(h1nch_lowB    , 20, 1.0, kRed+1, 1, 0);

    txt = TPaveText(0.2,0.82,0.4,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    leg = TLegend(0.2,0.6,0.4,0.7);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(h1nch_nominalB ,"B = 0.5 T (LHC22o pass4)", "LP");
    leg.AddEntry(h1nch_lowB   ,"B = 0.2T (LHC23zo, zp pass1)", "LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

 
    rootfile_nominalB.Close();
    rootfile_lowB.Close();
    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_nch08{1}.eps".format(date, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_nch08{1}.pdf".format(date, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_nch08{1}.png".format(date, suffix));


#_____________________________________________________________________
def draw_tpc_dedx_after(filename_data, period, suffix):
    rootfile = TFile.Open(filename_data, "READ");

    rootdir1_mumu = rootfile.Get("dalitz-mumu-qc");
    list_event = rootdir1_mumu.Get("Event");
    list_event.ls();
    h1z = list_event.FindObject("hZvtx_after"); 
    nev = h1z.GetEntries();

    list_track = rootdir1_mumu.Get("Track").FindObject("mmumu_0_1100_lowB");
    list_track.ls();
    h2 = list_track.FindObject("hTPCdEdx"); 
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);

    gStyle.SetPalette(55);
    c1 = TCanvas("c0", "c0",0,0,800,800);
    c1.SetMargin(0.1, 0.14, 0.1, 0.03);
    c1.SetTicks(1,1);
    c1.SetLogx(1);
    c1.SetLogz(1);
    
    h2.GetXaxis().SetTitleOffset(1.3);
    h2.GetYaxis().SetTitleOffset(1.3);
    h2.GetYaxis().SetTitle("TPC dE/dx");
    h2.SetContour(100);
    h2.GetXaxis().SetMoreLogLabels(1);
    h2.GetXaxis().SetRangeUser(0.03, 1);
    h2.GetYaxis().SetRangeUser(0., 200);
    #ntrack = h2.GetEntries();
    h2.Scale(1./nev);
    h2.SetMinimum(1e-8);
    h2.SetMaximum(2e-3);
    h2.Draw("colz");

    rootfile.Close();
    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_TPCdEdx{2}_after.eps".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_TPCdEdx{2}_after.pdf".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_TPCdEdx{2}_after.png".format(date, period, suffix));

#_____________________________________________________________________
def draw_tpc_dedx(filename_data, period, suffix):
    rootfile = TFile.Open(filename_data, "READ");
    rootdir1 = rootfile.Get("skimmer-primary-muon");
    rootdir2 = rootdir1.Get("Track");
    #rootdir2.ls();

    rootdir1_mumu = rootfile.Get("dalitz-mumu-qc");
    rootdir2_mumu = rootdir1_mumu.Get("Event");
    rootdir2_mumu.ls();
    h1z = rootdir2_mumu.FindObject("hZvtx_before"); 
    nev = h1z.GetEntries();

    h2 = rootdir2.Get("hTPCdEdx_Pin_before");
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);

    gStyle.SetPalette(55);
    c1 = TCanvas("c0", "c0",0,0,800,800);
    c1.SetMargin(0.1, 0.14, 0.1, 0.03);
    c1.SetTicks(1,1);
    c1.SetLogx(1);
    c1.SetLogz(1);
    
    h2.GetXaxis().SetTitleOffset(1.3);
    h2.GetYaxis().SetTitleOffset(1.3);
    h2.SetContour(100);
    h2.GetXaxis().SetMoreLogLabels(1);
    h2.GetXaxis().SetRangeUser(0.03, 1);
    h2.GetYaxis().SetRangeUser(0., 200);
    #ntrack = h2.GetEntries();
    h2.Scale(1./nev);
    h2.SetMinimum(1e-8);
    h2.SetMaximum(2e-3);
    h2.Draw("colz");

    rootfile.Close();
    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_TPCdEdx{2}.eps".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_TPCdEdx{2}.pdf".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_TPCdEdx{2}.png".format(date, period, suffix));

#_____________________________________________________________________
#_____________________________________________________________________
def draw_tof_beta_after(filename_data, period, suffix):
    rootfile = TFile.Open(filename_data, "READ");
    rootdir1_mumu = rootfile.Get("dalitz-mumu-qc");
    list_event = rootdir1_mumu.Get("Event");
    list_event.ls();
    h1z = list_event.FindObject("hZvtx_after"); 
    nev = h1z.GetEntries();


    list_track = rootdir1_mumu.Get("Track").FindObject("mmumu_0_1100_lowB");
    list_track.ls();

    h2 = list_track.FindObject("hTOFbeta");
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);

    gStyle.SetPalette(55);
    c1 = TCanvas("c0", "c0",0,0,800,800);
    c1.SetMargin(0.1, 0.14, 0.1, 0.03);
    c1.SetTicks(1,1);
    #c1.SetLogx(1);
    c1.SetLogz(1);
    
    h2.GetXaxis().SetTitleOffset(1.3);
    h2.GetYaxis().SetTitleOffset(1.3);
    h2.SetContour(100);
    h2.GetXaxis().SetMoreLogLabels(1);
    h2.GetXaxis().SetRangeUser(0.1, 1);
    h2.GetYaxis().SetRangeUser(0.6, 1.1);
    #ntrack = h2.GetEntries();
    h2.Scale(1./nev);
    h2.SetMinimum(1e-8);
    h2.SetMaximum(2e-3);
    h2.Draw("colz");

    #f1el = TF1("f1el", "sqrt(pow(x,2)/(pow(x,2) + pow([0],2)))", 0.1, 10);
    #f1el.SetNpx(1000);
    #f1el.SetParameter(0, 511e-6);
    #f1el.Draw("same");
    #f1el.SetLineColor(kBlack);
    #ROOT.SetOwnership(f1el, False);

    #f1mu = TF1("f1mu", "sqrt(pow(x,2)/(pow(x,2) + pow([0],2)))", 0.1, 10);
    #f1mu.SetNpx(1000);
    #f1mu.SetParameter(0, 0.106);
    #f1mu.Draw("same");
    #f1mu.SetLineColor(kBlack);
    #ROOT.SetOwnership(f1mu, False);

    #f1pi = TF1("f1pi", "sqrt(pow(x,2)/(pow(x,2) + pow([0],2)))", 0.1, 10);
    #f1pi.SetNpx(1000);
    #f1pi.SetParameter(0, 0.139);
    #f1pi.Draw("same");
    #f1pi.SetLineColor(kBlack);
    #ROOT.SetOwnership(f1pi, False);

    #f1ka = TF1("f1ka", "sqrt(pow(x,2)/(pow(x,2) + pow([0],2)))", 0.1, 10);
    #f1ka.SetNpx(1000);
    #f1ka.SetParameter(0, 0.494);
    #f1ka.Draw("same");
    #f1ka.SetLineColor(kBlack);
    #ROOT.SetOwnership(f1ka, False);

    #f1pr = TF1("f1pr", "sqrt(pow(x,2)/(pow(x,2) + pow([0],2)))", 0.1, 10);
    #f1pr.SetNpx(1000);
    #f1pr.SetParameter(0, 0.938);
    #f1pr.Draw("same");
    #f1pr.SetLineColor(kBlack);
    #ROOT.SetOwnership(f1pr, False);


    rootfile.Close();
    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_TOFbeta{2}_after.eps".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_TOFbeta{2}_after.pdf".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_TOFbeta{2}_after.png".format(date, period, suffix));

#_____________________________________________________________________
def draw_tof_beta(filename_data, period, suffix):
    rootfile = TFile.Open(filename_data, "READ");
    rootdir1 = rootfile.Get("skimmer-primary-muon");
    rootdir2 = rootdir1.Get("Track");
    #rootdir2.ls();

    rootdir1_mumu = rootfile.Get("dalitz-mumu-qc");
    rootdir2_mumu = rootdir1_mumu.Get("Event");
    rootdir2_mumu.ls();
    h1z = rootdir2_mumu.FindObject("hZvtx_before"); 
    nev = h1z.GetEntries();

    h2 = rootdir2.Get("hTOFbeta_Pin_before");
    h2.SetDirectory(0);
    ROOT.SetOwnership(h2, False);

    gStyle.SetPalette(55);
    c1 = TCanvas("c0", "c0",0,0,800,800);
    c1.SetMargin(0.1, 0.14, 0.1, 0.03);
    c1.SetTicks(1,1);
    #c1.SetLogx(1);
    c1.SetLogz(1);
    
    h2.GetXaxis().SetTitleOffset(1.3);
    h2.GetYaxis().SetTitleOffset(1.3);
    h2.SetContour(100);
    h2.GetXaxis().SetMoreLogLabels(1);
    h2.GetXaxis().SetRangeUser(0.1, 1);
    h2.GetYaxis().SetRangeUser(0.6, 1.1);
    #ntrack = h2.GetEntries();
    h2.Scale(1./nev);
    h2.SetMinimum(1e-8);
    h2.SetMaximum(2e-3);
    h2.Draw("colz");

    #f1el = TF1("f1el", "sqrt(pow(x,2)/(pow(x,2) + pow([0],2)))", 0.1, 10);
    #f1el.SetNpx(1000);
    #f1el.SetParameter(0, 511e-6);
    #f1el.Draw("same");
    #f1el.SetLineColor(kBlack);
    #ROOT.SetOwnership(f1el, False);

    #f1mu = TF1("f1mu", "sqrt(pow(x,2)/(pow(x,2) + pow([0],2)))", 0.1, 10);
    #f1mu.SetNpx(1000);
    #f1mu.SetParameter(0, 0.106);
    #f1mu.Draw("same");
    #f1mu.SetLineColor(kBlack);
    #ROOT.SetOwnership(f1mu, False);

    #f1pi = TF1("f1pi", "sqrt(pow(x,2)/(pow(x,2) + pow([0],2)))", 0.1, 10);
    #f1pi.SetNpx(1000);
    #f1pi.SetParameter(0, 0.139);
    #f1pi.Draw("same");
    #f1pi.SetLineColor(kBlack);
    #ROOT.SetOwnership(f1pi, False);

    #f1ka = TF1("f1ka", "sqrt(pow(x,2)/(pow(x,2) + pow([0],2)))", 0.1, 10);
    #f1ka.SetNpx(1000);
    #f1ka.SetParameter(0, 0.494);
    #f1ka.Draw("same");
    #f1ka.SetLineColor(kBlack);
    #ROOT.SetOwnership(f1ka, False);

    #f1pr = TF1("f1pr", "sqrt(pow(x,2)/(pow(x,2) + pow([0],2)))", 0.1, 10);
    #f1pr.SetNpx(1000);
    #f1pr.SetParameter(0, 0.938);
    #f1pr.Draw("same");
    #f1pr.SetLineColor(kBlack);
    #ROOT.SetOwnership(f1pr, False);


    rootfile.Close();
    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_TOFbeta{2}.eps".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_TOFbeta{2}.pdf".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_TOFbeta{2}.png".format(date, period, suffix));

#_____________________________________________________________________
def check_material_phi(filename_data, filename_mc, cutname, period, rmin, rmax, suffix):
    rootfile_data = TFile.Open(filename_data, "READ");
    rootfile_mc   = TFile.Open(filename_mc  , "READ");

    rootdir_mc_rec = rootfile_mc.Get("pcm-qc-mc");
    list_v0_mc_rec = rootdir_mc_rec.Get("V0");
    list_ev_mc_rec = rootdir_mc_rec.Get("Event");
    list_cut_mc_rec = list_v0_mc_rec.FindObject(cutname);
    list_cut_mc_rec.Print();

    h1ev_mc_rec = list_ev_mc_rec.FindObject("hCollisionCounter").Clone("h1ev");
    nev_mc_rec = h1ev_mc_rec.GetBinContent(4);
    print("nev_mc_rec = {0:e}".format(nev_mc_rec));

    h1nch_mc_rec = list_ev_mc_rec.FindObject("hMultNTracksPV");
    nch_mc_rec = h1nch_mc_rec.GetMean();
    print("nch pv = ", nch_mc_rec);
    dr = rmax - rmin;

    h2_mc_rec = list_cut_mc_rec.FindObject("hGammaRPhi_recalc");
    h2_mc_rec.Sumw2();
    h2_mc_rec.RebinX(5);
    h2_mc_rec.SetDirectory(0);
    ROOT.SetOwnership(h2_mc_rec, False);

    bin0 = h2_mc_rec.GetYaxis().FindBin(rmin + 1e-3);
    bin1 = h2_mc_rec.GetYaxis().FindBin(rmax - 1e-3);
    h1_mc_rec = h2_mc_rec.ProjectionX("h1_mc_rec", bin0, bin1, "");
    h1_mc_rec.SetDirectory(0);
    ROOT.SetOwnership(h1_mc_rec, False);
    h1_mc_rec.Scale(1,"width");
    h1_mc_rec.Scale(1/dr);
    h1_mc_rec.Scale(1/nev_mc_rec);
    h1_mc_rec.Scale(1/nch_mc_rec);#nch
    make_common_style(h1_mc_rec, 20, 1.0, kBlue+1, 1, 0);

    rootdir_data = rootfile_data.Get("pcm-qc");
    list_v0_data = rootdir_data.Get("V0");
    list_ev_data = rootdir_data.Get("Event");
    list_cut_data = list_v0_data.FindObject(cutname);
    list_cut_data.Print();

    h1ev_data = list_ev_data.FindObject("hCollisionCounter").Clone("h1ev");
    nev_data = h1ev_data.GetBinContent(4);
    print("nev_data = {0:e}".format(nev_data));

    h1nch_data = list_ev_data.FindObject("hMultNTracksPV");
    nch_data = h1nch_data.GetMean();
    print("nch pv = ", nch_data);

    h2_data = list_cut_data.FindObject("hGammaRPhi_recalc");
    h2_data.Sumw2();
    h2_data.RebinX(5);
    h2_data.SetDirectory(0);
    ROOT.SetOwnership(h2_data, False);
    bin0 = h2_data.GetYaxis().FindBin(rmin + 1e-3);
    bin1 = h2_data.GetYaxis().FindBin(rmax - 1e-3);
    h1_data = h2_data.ProjectionX("h1_data", bin0, bin1, "");
    h1_data.SetDirectory(0);
    ROOT.SetOwnership(h1_data, False);
    h1_data.Scale(1,"width");
    h1_data.Scale(1/dr);
    h1_data.Scale(1/nev_data);
    h1_data.Scale(1/nch_data);
    make_common_style(h1_data, 20, 1.0, kRed+1, 1, 0);

    ymax = max(h1_data.GetMaximum() , h1_mc_rec.GetMaximum()) * 1.6;
    ymin = max(h1_data.GetMinimum() , h1_mc_rec.GetMinimum()) * -0.1;

    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.Divide(1,2,0,0);
    p1 = c1.cd(1);
    p1.SetPad(0,0.35,1,1);
    p1.SetMargin(0.15,0.02,0.0,0.06);
    p1.SetTicks(1,1);
    #p1.SetLogy(1);

    frame1 = p1.DrawFrame(0, ymin, TMath.TwoPi(), ymax);
    frame1.GetXaxis().SetTitle("conversion point #it{#varphi} (rad.)");
    frame1.GetYaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}_{#gamma}}{d#it{#varphi} d#it{R}_{xy}} (rad. #upoint cm)^{-1}");
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetMaxDigits(3);
    #frame1.GetXaxis().SetLabelOffset(0.01);
    #frame1.GetYaxis().SetLabelOffset(0.01);
    ROOT.SetOwnership(frame1,False);

    #h1_mc_gen.Draw("E0h,same");
    h1_mc_rec.Draw("E0h,same");
    h1_data.Draw("E0h,same");

    txt = TPaveText(0.2,0.72,0.4,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("|#it{#eta}_{#gamma}| < 0.9");
    txt.AddText("{0:2.1f} < #it{{R}}_{{xy}} < {1:2.1f} cm".format(rmin, rmax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    leg = TLegend(0.2,0.6,0.4,0.7);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    leg.AddEntry(h1_data   ,"Data #gamma candidates (LHC22f pass4 520472)","LP");
    leg.AddEntry(h1_mc_rec ,"M.C. rec. primary #gamma (LHC23d1f 520472)","LP");
    #leg.AddEntry(h1_mc_gen ,"M.C. gen. (LHC23d1)","LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);


    p2 = c1.cd(2);
    p2.SetPad(0,0,1,0.35);
    p2.SetMargin(0.15,0.02,0.22,0.0);
    p2.SetTicks(1,1);
    #p2.SetLogy(1);

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

    line1 = TLine(0,1,TMath.TwoPi(),1);
    line1.SetLineColor(kBlack);
    line1.SetLineStyle(1);
    line1.SetLineWidth(2);
    line1.Draw("");
    ROOT.SetOwnership(line1,False);

    h1ratio = h1_data.Clone("h1ratio");
    h1ratio.Reset();
    h1ratio.Sumw2();
    h1ratio.Divide(h1_data, h1_mc_rec, 1., 1., "G");
    h1ratio.Draw("E0,same");
    h1ratio.SetDirectory(0);
    ROOT.SetOwnership(h1ratio,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_material_budget_vs_phi_rxy{2:2.1f}_{3:2.1f}cm_{4}{5}.eps".format(date, period, rmin, rmax, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_material_budget_vs_phi_rxy{2:2.1f}_{3:2.1f}cm_{4}{5}.pdf".format(date, period, rmin, rmax, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_material_budget_vs_phi_rxy{2:2.1f}_{3:2.1f}cm_{4}{5}.png".format(date, period, rmin, rmax, cutname, suffix));

    rootfile_data.Close();
    rootfile_mc  .Close();
#_____________________________________________________________________
#_____________________________________________________________________
def check_material_rxy(filename_data, filename_mc, cutname, period, suffix):
    rootfile_data = TFile.Open(filename_data, "READ");
    rootfile_mc   = TFile.Open(filename_mc  , "READ");

    rootdir_mc  = rootfile_mc.Get("pcm-qc-mc");
    list_gen = rootdir_mc.Get("Generated");
    list_gen.Print();

    h1ev_gen = list_gen.FindObject("hCollisionCounter").Clone("h1ev");
    nev_gen = h1ev_gen.GetBinContent(4);
    print("nev_gen = {0:e}".format(nev_gen));


    h2_mc_gen = list_gen.FindObject("hPhotonRZ");
    h2_mc_gen.Sumw2();
    h1_mc_gen = h2_mc_gen.ProjectionY("h1_mc_gen");
    h2_mc_gen.SetDirectory(0);
    h1_mc_gen.SetDirectory(0);
    ROOT.SetOwnership(h2_mc_gen,False);
    ROOT.SetOwnership(h1_mc_gen,False);
    h1_mc_gen.RebinX(5);
    h1_mc_gen.Scale(1,"width");
    h1_mc_gen.Scale(1/nev_gen);
    h1_mc_gen.Scale(1/11.1);#nch
    make_common_style(h1_mc_gen, 20, 1.0, kBlack, 1, 0);

    rootdir_mc_rec = rootfile_mc.Get("pcm-qc-mc");
    list_v0_mc_rec = rootdir_mc_rec.Get("V0");
    list_ev_mc_rec = rootdir_mc_rec.Get("Event");
    list_cut_mc_rec = list_v0_mc_rec.FindObject(cutname);
    list_cut_mc_rec.Print();

    h1ev_mc_rec = list_ev_mc_rec.FindObject("hCollisionCounter").Clone("h1ev");
    nev_mc_rec = h1ev_mc_rec.GetBinContent(4);
    print("nev_mc_rec = {0:e}".format(nev_mc_rec));

    h1nch_mc_rec = list_ev_mc_rec.FindObject("hMultNTracksPV");
    nch_mc_rec = h1nch_mc_rec.GetMean();
    print("nch pv = ", nch_mc_rec);


    h2_mc_rec = list_cut_mc_rec.FindObject("hGammaRPhi_recalc");
    h2_mc_rec.Sumw2();
    h2_mc_rec.RebinX(2);
    h2_mc_rec.SetDirectory(0);
    ROOT.SetOwnership(h2_mc_rec, False);
    h1_mc_rec = h2_mc_rec.ProjectionY("h1_mc_rec");
    h1_mc_rec.SetDirectory(0);
    ROOT.SetOwnership(h1_mc_rec, False);
    h1_mc_rec.Scale(1,"width");
    h1_mc_rec.Scale(1/nev_mc_rec);
    h1_mc_rec.Scale(1/nch_mc_rec);#nch
    make_common_style(h1_mc_rec, 20, 1.0, kBlue+1, 1, 0);

    rootdir_data = rootfile_data.Get("pcm-qc");
    list_v0_data = rootdir_data.Get("V0");
    list_ev_data = rootdir_data.Get("Event");
    list_cut_data = list_v0_data.FindObject(cutname);
    list_cut_data.Print();

    h1ev_data = list_ev_data.FindObject("hCollisionCounter").Clone("h1ev");
    nev_data = h1ev_data.GetBinContent(4);
    print("nev_data = {0:e}".format(nev_data));

    h1nch_data = list_ev_data.FindObject("hMultNTracksPV");
    nch_data = h1nch_data.GetMean();
    print("nch pv = ", nch_data);

    h2_data = list_cut_data.FindObject("hGammaRPhi_recalc");
    h2_data.Sumw2();
    h2_data.RebinX(2);
    h2_data.SetDirectory(0);
    ROOT.SetOwnership(h2_data, False);
    h1_data = h2_data.ProjectionY("h1_data");
    h1_data.SetDirectory(0);
    ROOT.SetOwnership(h1_data, False);
    h1_data.Scale(1,"width");
    h1_data.Scale(1/nev_data);
    h1_data.Scale(1/nch_data);
    make_common_style(h1_data, 20, 1.0, kRed+1, 1, 0);

    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.Divide(1,2,0,0);
    p1 = c1.cd(1);
    p1.SetPad(0,0.35,1,1);
    p1.SetMargin(0.16,0.03,0.0,0.03);
    p1.SetTicks(1,1);
    p1.SetLogy(1);

    frame1 = p1.DrawFrame(0,4e-8,90,1e-2);
    frame1.GetXaxis().SetTitle("conversion radius #it{R}_{xy} (cm)");
    #frame1.GetYaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}> #upoint #it{N}_{ev}} #frac{d#it{N}_{#gamma}}{d#it{R}_{xy}} (cm)^{-1}");
    frame1.GetYaxis().SetTitle("#frac{1}{<#it{N}_{ch}^{PV}>} #frac{1}{#it{N}_{ev}} #frac{d#it{N}_{#gamma}}{d#it{R}_{xy}} (cm)^{-1}");
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    #frame1.GetXaxis().SetLabelOffset(0.01);
    #frame1.GetYaxis().SetLabelOffset(0.01);
    ROOT.SetOwnership(frame1,False);

    #h1_mc_gen.Draw("E0h,same");
    h1_mc_rec.Draw("E0h,same");
    h1_data.Draw("E0h,same");


    txt = TPaveText(0.2,0.8,0.4,0.95,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.05);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("|#it{#eta}_{#gamma}| < 0.9");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    leg = TLegend(0.2,0.7,0.4,0.8);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);
    leg.AddEntry(h1_data   ,"Data #gamma candidates (LHC22f pass4 520472)","LP");
    leg.AddEntry(h1_mc_rec ,"M.C. rec. primary #gamma (LHC23d1f 520472)","LP");
    #leg.AddEntry(h1_mc_gen ,"M.C. gen. (LHC23d1)","LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);


    p2 = c1.cd(2);
    p2.SetPad(0,0,1,0.35);
    p2.SetMargin(0.16,0.03,0.22,0.0);
    p2.SetTicks(1,1);
    #p2.SetLogy(1);

    frame2 = p2.DrawFrame(0,0.,90,3.2);
    frame2.GetXaxis().SetTitle("#it{R}_{xy} (cm)");
    frame2.GetYaxis().SetTitle("#frac{Data}{M.C.}");
    frame2.GetXaxis().SetTitleSize(0.10);
    frame2.GetYaxis().SetTitleSize(0.10);
    frame2.GetXaxis().SetTitleOffset(1.0);
    frame2.GetYaxis().SetTitleOffset(0.7);
    frame2.GetXaxis().SetLabelSize(0.1);
    frame2.GetYaxis().SetLabelSize(0.1);
    frame2.GetYaxis().CenterTitle(True);
    ROOT.SetOwnership(frame2,False);

    line1 = TLine(0,1,90,1);
    line1.SetLineColor(kBlack);
    line1.SetLineStyle(1);
    line1.SetLineWidth(2);
    line1.Draw("");
    ROOT.SetOwnership(line1,False);

    h1ratio = h1_data.Clone("h1ratio");
    h1ratio.Reset();
    h1ratio.Sumw2();
    h1ratio.Divide(h1_data, h1_mc_rec, 1., 1., "G");
    h1ratio.Draw("E0,same");
    h1ratio.SetDirectory(0);
    ROOT.SetOwnership(h1ratio,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_material_budget_vs_rxy_{2}{3}.eps".format(date, period, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_material_budget_vs_rxy_{2}{3}.pdf".format(date, period, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_material_budget_vs_rxy_{2}{3}.png".format(date, period, cutname, suffix));

    rootfile_data.Close();
    rootfile_mc  .Close();
#_____________________________________________________________________
def draw_rz_mc_gen(filename, subsystem, period, suffix=""):
    print("reading...", filename);
    rootfile = TFile.Open(filename,"READ");
    rootdir = rootfile.Get("material-budget-mc");
    rootdir.ls();
    list_gen  = rootdir.Get("Generated");
    list_gen.Print();

    h1ev = list_gen.FindObject("hCollisionCounter").Clone("h1ev");
    nev = h1ev.GetBinContent(4);
    print("nev = {0:e}".format(nev));

    h2 = list_gen.FindObject("hPhotonRZ");
    ROOT.SetOwnership(h2, False);
    h2.SetDirectory(0);
    h2.SetContour(100);

#    for ix in range(0, h2.GetNbinsX()):
#        for iy in range(0, h2.GetNbinsY()):
#            x = h2.GetXaxis().GetBinCenter(ix + 1);
#            y = h2.GetYaxis().GetBinCenter(iy + 1);
#            #dx = h2.GetXaxis().GetBinWidth(ix + 1);
#            #dy = h2.GetYaxis().GetBinWidth(iy + 1);
#            #n = h2.GetBinContent(ix+1, iy+1);
#            #h2.SetBinContent(ix+1, iy+1, n/dx/dy);
#            rxy = math.sqrt(x*x + y*y);
#            if rxy > 90.:
#                h2.SetBinContent(ix+1, iy+1, 0);
 
    gStyle.SetPalette(55);
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetTicks(1,1);
    c1.SetLogz(1);
    c1.SetMargin(0.11,0.12,0.1,0.03);
    #h2.RebinX(4);
    #h2.RebinY(4);

    frame1 = c1.DrawFrame(-100,0,100,100);
    frame1.GetXaxis().SetTitle("conversion point Z (cm)");
    frame1.GetYaxis().SetTitle("conversion point R_{xy} (cm)");
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.5);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    #frame1.GetXaxis().SetMoreLogLabels(True);
    #frame1.GetYaxis().SetMaxDigits(3);

    h2.Scale(1/nev);
    print(h2.GetMaximum());
    h2.SetMinimum(1e-7);
    h2.SetMaximum(3e-3);
    h2.Draw("colz, same");

#    f09 = TF1("f09", "abs(x) * TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-0.9)))", -100, +100);
#    f09.SetNpx(1000);
#    f09.SetLineColor(kBlack);
#    f09.SetLineWidth(2);
#    f09.SetLineStyle(1);
#    f09.Draw("same");
#    ROOT.SetOwnership(f09, False);
#
#    f12 = TF1("f12", "abs(x) * TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-1.2)))", -100, +100);
#    f12.SetNpx(1000);
#    f12.SetLineColor(kBlack);
#    f12.SetLineWidth(2);
#    f12.SetLineStyle(2);
#    f12.Draw("same");
#    ROOT.SetOwnership(f12, False);


    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_MC_{1}_rz{2}.eps".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_MC_{1}_rz{2}.pdf".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_MC_{1}_rz{2}.png".format(date, period, suffix));

    rootfile.Close();

#_____________________________________________________________________
#_____________________________________________________________________
def draw_xy_mc_gen(filename, subsystem, period, suffix=""):
    print("reading...", filename);
    rootfile = TFile.Open(filename,"READ");
    rootdir = rootfile.Get("material-budget-mc");
    rootdir.ls();
    list_gen  = rootdir.Get("Generated");
    list_gen.Print();

    h1ev = list_gen.FindObject("hCollisionCounter").Clone("h1ev");
    nev = h1ev.GetBinContent(4);
    print("nev = {0:e}".format(nev));

    h2 = list_gen.FindObject("hPhotonRxy");
    ROOT.SetOwnership(h2, False);
    h2.SetDirectory(0);
    h2.SetContour(100);

#    for ix in range(0, h2.GetNbinsX()):
#        for iy in range(0, h2.GetNbinsY()):
#            x = h2.GetXaxis().GetBinCenter(ix + 1);
#            y = h2.GetYaxis().GetBinCenter(iy + 1);
#            #dx = h2.GetXaxis().GetBinWidth(ix + 1);
#            #dy = h2.GetYaxis().GetBinWidth(iy + 1);
#            #n = h2.GetBinContent(ix+1, iy+1);
#            #h2.SetBinContent(ix+1, iy+1, n/dx/dy);
#            rxy = math.sqrt(x*x + y*y);
#            if rxy > 90.:
#                h2.SetBinContent(ix+1, iy+1, 0);
 
    gStyle.SetPalette(55);
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetTicks(1,1);
    c1.SetLogz(1);
    c1.SetMargin(0.11,0.12,0.1,0.03);
    h2.RebinX(4);
    h2.RebinY(4);

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
    #frame1.GetXaxis().SetMoreLogLabels(True);
    #frame1.GetYaxis().SetMaxDigits(3);

    h2.Scale(1/nev);
    print(h2.GetMaximum());
    h2.SetMinimum(1e-7);
    h2.SetMaximum(3e-3);
    h2.Draw("colz, same");

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_MC_{1}_rxy{2}.eps".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_MC_{1}_rxy{2}.pdf".format(date, period, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_MC_{1}_rxy{2}.png".format(date, period, suffix));

    rootfile.Close();

#_____________________________________________________________________
def draw_rz_ITSTPC(filename, subsystem, period, cutname, is_recalculated=True, suffix=""):
    print("reading...", filename);
    rootfile = TFile.Open(filename,"READ");
    rootdir = rootfile.Get("pcm-qc");
    #rootdir.ls();
    list_v0  = rootdir.Get("V0");
    list_ev = rootdir.Get("Event");
    #list_v0.Print();
    list_cut = list_v0.FindObject(cutname);

    h1ev = list_ev.FindObject("hCollisionCounter").Clone("h1ev");
    nev = h1ev.GetBinContent(4);
    print("nev = {0:e}".format(nev));

    h1nch = list_ev.FindObject("hMultNTracksPV").Clone("h1nch");
    nch = h1nch.GetMean();

    h2 = list_cut.FindObject("hRadius");
    ROOT.SetOwnership(h2, False);
    h2.SetDirectory(0);
    h2.SetContour(1000);
    h2.GetZaxis().SetTitle("N_{#gamma}/N_{ev}/<N_{ch}>");
    h2.GetZaxis().SetTitleSize(0.04);
    h2.GetZaxis().SetTitleOffset(1.5);
 
    gStyle.SetPalette(55);
    c1 = TCanvas("c0","c0",0,0,850,750);
    c1.SetTicks(1,1);
    c1.SetLogz(1);
    c1.SetMargin(0.11,0.18,0.10,0.02);

    frame1 = c1.DrawFrame(-100,0,100,100);
    frame1.GetXaxis().SetTitle("conversion point Z (cm)");
    frame1.GetYaxis().SetTitle("conversion point R_{xy} (cm)");
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    #frame1.GetXaxis().SetMoreLogLabels(True);
    #frame1.GetYaxis().SetMaxDigits(3);
    h2.Scale(1/nev/nch);
    print(h2.GetMinimum());
    print(h2.GetMaximum());
    h2.SetMinimum(6e-12);
    h2.SetMaximum(4e-6);
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

    txt = TPaveText(0.12, 0.15, 0.4, 0.25, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica bold
    txt.SetTextSize(0.03);
    #txt.AddText("ALICE Performance");
    #txt.AddText("pp, #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("B = 0.2 T");
    #txt.AddText("Pb-Pb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("B = 0.5 T");
    #txt.AddText(cutname);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    #txt_pcm = TPaveText(0.4,0.92,0.7,1.0,"NDC");
    #txt_pcm.SetFillColor(kWhite);
    #txt_pcm.SetFillStyle(0);
    #txt_pcm.SetBorderSize(0);
    #txt_pcm.SetTextAlign(12);#middle,left
    #txt_pcm.SetTextFont(42);#helvetica
    #txt_pcm.SetTextSize(0.035);
    #txt_pcm.AddText("Photon Conversion Method (#gamma #rightarrow e^{+}e^{#minus})");
    #txt_pcm.AddText("|#it{#eta}_{#gamma}| < 0.9");
    #txt_pcm.Draw();
    #ROOT.SetOwnership(txt_pcm,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass2_rxy_{2}_{4}.eps".format(date, period, cutname, is_recalculated, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass2_rxy_{2}_{4}.pdf".format(date, period, cutname, is_recalculated, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass2_rxy_{2}_{4}.png".format(date, period, cutname, is_recalculated, suffix));
    #c1.SaveAs("{0}_PbPb_5.36TeV_{1}_pass2_rxy_{2}_{4}.eps".format(date, period, cutname, is_recalculated, suffix));
    #c1.SaveAs("{0}_PbPb_5.36TeV_{1}_pass2_rxy_{2}_{4}.pdf".format(date, period, cutname, is_recalculated, suffix));
    #c1.SaveAs("{0}_PbPb_5.36TeV_{1}_pass2_rxy_{2}_{4}.png".format(date, period, cutname, is_recalculated, suffix));

    rootfile.Close();

#_____________________________________________________________________
def draw_xy_ITSTPC(filename, subsystem, period, cutname, is_recalculated=True, suffix=""):
    print("reading...", filename);
    rootfile = TFile.Open(filename,"READ");
    rootdir = rootfile.Get("pcm-qc");
    #rootdir.ls();
    list_v0  = rootdir.Get("V0");
    list_ev = rootdir.Get("Event");
    #list_v0.Print();
    list_cut = list_v0.FindObject(cutname);

    h1ev = list_ev.FindObject("hCollisionCounter").Clone("h1ev");
    nev = h1ev.GetBinContent(4);
    print("nev = {0:e}".format(nev));

    h1nch = list_ev.FindObject("hMultNTracksPV").Clone("h1nch");
    nch = h1nch.GetMean();

    h2 = list_cut.FindObject("hGammaRxy");
    if is_recalculated:
        h2 = list_cut.FindObject("hGammaRxy_recalc");
    ROOT.SetOwnership(h2, False);
    h2.SetDirectory(0);
    h2.SetContour(1000);
    h2.GetZaxis().SetTitle("N_{#gamma}/N_{ev}/<N_{ch}>");
    h2.GetZaxis().SetTitleSize(0.04);
    h2.GetZaxis().SetTitleOffset(1.5);
 
    gStyle.SetPalette(55);
    c1 = TCanvas("c0","c0",0,0,850,750);
    c1.SetTicks(1,1);
    c1.SetLogz(1);
    c1.SetMargin(0.11,0.18,0.10,0.02);

    frame1 = c1.DrawFrame(-100,-100,100,100);
    frame1.GetXaxis().SetTitle("conversion point X (cm)");
    frame1.GetYaxis().SetTitle("conversion point Y (cm)");
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    #frame1.GetXaxis().SetMoreLogLabels(True);
    #frame1.GetYaxis().SetMaxDigits(3);
    h2.Scale(1/nev/nch);
    print(h2.GetMinimum());
    print(h2.GetMaximum());
    h2.SetMinimum(6e-12);
    h2.SetMaximum(4e-6);
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

    txt = TPaveText(0.12, 0.88, 0.4, 0.96,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica bold
    txt.SetTextSize(0.03);
    #txt.AddText("ALICE Performance");
    #txt.AddText("pp, #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("B = 0.2 T");
    txt.AddText("Pb-Pb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("B = 0.5 T");
    #txt.AddText(cutname);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    #txt_pcm = TPaveText(0.4,0.92,0.7,1.0,"NDC");
    #txt_pcm.SetFillColor(kWhite);
    #txt_pcm.SetFillStyle(0);
    #txt_pcm.SetBorderSize(0);
    #txt_pcm.SetTextAlign(12);#middle,left
    #txt_pcm.SetTextFont(42);#helvetica
    #txt_pcm.SetTextSize(0.035);
    #txt_pcm.AddText("Photon Conversion Method (#gamma #rightarrow e^{+}e^{#minus})");
    #txt_pcm.AddText("|#it{#eta}_{#gamma}| < 0.9");
    #txt_pcm.Draw();
    #ROOT.SetOwnership(txt_pcm,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    #c1.SaveAs("{0}_pp_13.6TeV_{1}_pass2_rxy_{2}_{4}.eps".format(date, period, cutname, is_recalculated, suffix));
    #c1.SaveAs("{0}_pp_13.6TeV_{1}_pass2_rxy_{2}_{4}.pdf".format(date, period, cutname, is_recalculated, suffix));
    #c1.SaveAs("{0}_pp_13.6TeV_{1}_pass2_rxy_{2}_{4}.png".format(date, period, cutname, is_recalculated, suffix));
    c1.SaveAs("{0}_PbPb_5.36TeV_{1}_pass2_rxy_{2}_{4}.eps".format(date, period, cutname, is_recalculated, suffix));
    c1.SaveAs("{0}_PbPb_5.36TeV_{1}_pass2_rxy_{2}_{4}.pdf".format(date, period, cutname, is_recalculated, suffix));
    c1.SaveAs("{0}_PbPb_5.36TeV_{1}_pass2_rxy_{2}_{4}.png".format(date, period, cutname, is_recalculated, suffix));

    rootfile.Close();

#_____________________________________________________________________
def draw_xy(filename, subsystem, period, cutname, is_recalculated=True, suffix=""):
    print("reading...", filename);
    rootfile = TFile.Open(filename,"READ");
    rootdir = rootfile.Get("pcm-qc");
    #rootdir.ls();
    list_v0  = rootdir.Get("V0");
    list_ev = rootdir.Get("Event");
    #list_v0.Print();
    list_cut = list_v0.FindObject(cutname);

    h1ev = list_ev.FindObject("hCollisionCounter").Clone("h1ev");
    nev = h1ev.GetBinContent(4);
    print("nev = {0:e}".format(nev));

    h2 = list_cut.FindObject("hGammaRxy");
    if is_recalculated:
        h2 = list_cut.FindObject("hGammaRxy_recalc");
    ROOT.SetOwnership(h2, False);
    h2.SetDirectory(0);
    h2.SetContour(1000);
    h2.GetZaxis().SetTitle("N_{#gamma}/N_{ev}");
    h2.GetZaxis().SetTitleSize(0.04);
    h2.GetZaxis().SetTitleOffset(0.7);
 
    gStyle.SetPalette(55);
    c1 = TCanvas("c0","c0",0,0,850,850);
    c1.SetTicks(1,1);
    c1.SetLogz(1);
    c1.SetMargin(0.11,0.12,0.10,0.09);

    frame1 = c1.DrawFrame(-100,-100,100,100);
    frame1.GetXaxis().SetTitle("conversion point X (cm)");
    frame1.GetYaxis().SetTitle("conversion point Y (cm)");
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    #frame1.GetXaxis().SetMoreLogLabels(True);
    #frame1.GetYaxis().SetMaxDigits(3);
    h2.Scale(1/nev);
    print(h2.GetMinimum());
    print(h2.GetMaximum());
    h2.SetMinimum(1e-8);
    h2.SetMaximum(5e-4);
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

    txt = TPaveText(0.10,0.92,0.4,1.0,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica bold
    txt.SetTextSize(0.035);
    txt.AddText("ALICE Performance");
    txt.AddText("pp, #sqrt{#it{s}} = 13.6 TeV");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    txt_pcm = TPaveText(0.4,0.92,0.7,1.0,"NDC");
    txt_pcm.SetFillColor(kWhite);
    txt_pcm.SetFillStyle(0);
    txt_pcm.SetBorderSize(0);
    txt_pcm.SetTextAlign(12);#middle,left
    txt_pcm.SetTextFont(42);#helvetica
    txt_pcm.SetTextSize(0.035);
    txt_pcm.AddText("Photon Conversion Method (#gamma #rightarrow e^{+}e^{#minus})");
    txt_pcm.AddText("|#it{#eta}_{#gamma}| < 0.9");
    txt_pcm.Draw();
    ROOT.SetOwnership(txt_pcm,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_rxy_{2}_recalc{3:d}_{4}.eps".format(date, period, cutname, is_recalculated, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_rxy_{2}_recalc{3:d}_{4}.pdf".format(date, period, cutname, is_recalculated, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_rxy_{2}_recalc{3:d}_{4}.png".format(date, period, cutname, is_recalculated, suffix));

    rootfile.Close();

#_____________________________________________________________________
def draw_raw_yield_inc_gamma_low_high(filename_low, filename_high, period_low, period_high, subsystem, cutname, suffix=""):
    print("reading...", filename_low, filename_high);
    arr_pt = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6,7,8,9,10], dtype=float);

    rootfile_low = TFile.Open(filename_low,"READ");
    rootdir_low = rootfile_low.Get("pcm-qc");
    list_low_v0 = rootdir_low.Get("V0");
    list_low_ev = rootdir_low.Get("Event");

    h1ev_low = list_low_ev.FindObject("hCollisionCounter").Clone("h1ev");
    nev_low = h1ev_low.GetBinContent(4);
    print("nev_low = {0:e}".format(nev_low));

    list_low_v0_cut = list_low_v0.FindObject(cutname);
    h1pt_low_org = list_low_v0_cut.FindObject("hPt").Clone("hPt_{0}_org".format(cutname));
    h1pt_low = rebin_histogram(h1pt_low_org, arr_pt, True, False);
    make_common_style(h1pt_low  , 20, 1.2, kRed+1, 1, 0);
    h1pt_low.Scale(1/nev_low);
    h1pt_low.SetDirectory(0);
    ROOT.SetOwnership(h1pt_low,False);

    rootfile_high = TFile.Open(filename_high,"READ");
    rootdir_high = rootfile_high.Get("pcm-qc");
    list_high_v0 = rootdir_high.Get("V0");
    list_high_ev = rootdir_high.Get("Event");

    h1ev_high = list_high_ev.FindObject("hCollisionCounter").Clone("h1ev");
    nev_high = h1ev_high.GetBinContent(4);
    print("nev_high = {0:e}".format(nev_high));

    list_high_v0_cut = list_high_v0.FindObject(cutname);
    h1pt_high_org = list_high_v0_cut.FindObject("hPt").Clone("hPt_{0}_org".format(cutname));
    h1pt_high = rebin_histogram(h1pt_high_org, arr_pt, True, False);
    make_common_style(h1pt_high  , 20, 1.2, kBlue+1, 1, 0);
    h1pt_high.Scale(1/nev_high);
    h1pt_high.SetDirectory(0);
    ROOT.SetOwnership(h1pt_high,False);

    gStyle.SetPalette(55);
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.Divide(1,2,0,0);

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1 ,1);
    p1.SetMargin(0.14,0.02,0.,0.02);
    p1.SetTicks(1,1);
    p1.SetLogy(1);
    frame1 = p1.DrawFrame(0,6e-8,10,1e+0);
    frame1.GetXaxis().SetTitle("#it{p}_{T,#gamma} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T,#gamma}} (GeV/#it{c})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.25);
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    ROOT.SetOwnership(frame1,False);
    #frame1.GetXaxis().SetMoreLogLabels(True);
    #frame1.GetYaxis().SetMaxDigits(3);

    h1pt_low.Draw("E0,same");
    h1pt_high.Draw("E0,same");

    collrate_low  = "";
    collrate_high = "";

    if "22q" in period_low:
        collrate_low = "15 kHz";
    elif "22f" in period_low:
        collrate_low = "10 kHz";

    if "22o" in period_high or "22m" in period_high :
        collrate_high = "500 kHz";
    elif "22q" in period_high:
        collrate_high = "15 kHz";


    txt = TPaveText(0.65,0.70,0.9,0.94,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica bold
    txt.SetTextSize(0.052);
    txt.AddText("Work in Progress");
    txt.AddText("pp, #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("PCM (#gamma #rightarrow e^{+}e^{#minus})");
    txt.AddText("|#it{#eta}_{#gamma}| < 0.9");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    cutname_tmp = "";
    if "ITSTPC" in cutname:
        cutname_tmp = "V0s with ITS-TPC matched tracks";
    elif "TPConly" in cutname:
        cutname_tmp = "V0s with TPConly tracks";
    elif "ITSonly" in cutname:
        cutname_tmp = "V0s with ITSonly tracks";
    elif "analysis" in cutname:
        cutname_tmp = "V0s with Any tracks";

    leg = TLegend(0.45,0.53,0.65,0.68);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.052);
    leg.SetHeader(cutname_tmp);
    leg.AddEntry(h1pt_low , "{0}, {1}".format(period_low, collrate_low), "P");
    leg.AddEntry(h1pt_high, "{0}, {1}".format(period_high, collrate_high), "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1 ,0.35);
    p2.SetMargin(0.14,0.02,0.25,0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0, 0.1, 10, 2.1);
    frame2.GetXaxis().SetTitle("#it{p}_{T,#gamma} (GeV/#it{c})");
    frame2.GetYaxis().SetTitle("#frac{{{0}}}{{{1}}}".format(period_high, period_low));
    frame2.GetXaxis().SetTitleOffset(1.1);
    frame2.GetYaxis().SetTitleOffset(0.6);
    frame2.GetXaxis().SetTitleSize(0.10);
    frame2.GetYaxis().SetTitleSize(0.10);
    frame2.GetXaxis().SetLabelSize(0.10);
    frame2.GetYaxis().SetLabelSize(0.10);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().CenterTitle(True);
    ROOT.SetOwnership(frame2,False);

    h1ratio = h1pt_high.Clone("h1ratio");
    h1ratio.Reset();
    h1ratio.Sumw2();
    h1ratio.Divide(h1pt_high, h1pt_low, 1., 1., "G");
    h1ratio.Draw("E0same");
    h1ratio.SetDirectory(0);
    ROOT.SetOwnership(h1ratio, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_pass4_inc_gamma_{1}_{2}_{3}{4}.eps".format(date, period_low, period_high, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_pass4_inc_gamma_{1}_{2}_{3}{4}.pdf".format(date, period_low, period_high, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_pass4_inc_gamma_{1}_{2}_{3}{4}.png".format(date, period_low, period_high, cutname, suffix));

    rootfile_low.Close();
    rootfile_high.Close();

#_____________________________________________________________________
def draw_raw_yield_inc_gamma(filename, subsystem, period, cutname_tmp, suffix=""):
    print("reading...", filename);
    rootfile = TFile.Open(filename,"READ");
    rootdir = rootfile.Get("pcm-qc");
    #rootdir.ls();
    list_v0  = rootdir.Get("V0");
    list_ev = rootdir.Get("Event");
    #list_v0.Print();

    h1ev = list_ev.FindObject("hCollisionCounter").Clone("h1ev");
    nev = h1ev.GetBinContent(4);
    print("nev = {0:e}".format(nev));

    cutnames = ["qc_ITSTPC", "qc_TPConly"];
    colors = [kBlack, kRed+1];
 
    gStyle.SetPalette(55);
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.Divide(1,2,0,0);

    p1 = c1.cd(1);
    p1.SetPad(0, 0.35, 1 ,1);
    p1.SetMargin(0.14,0.02,0.,0.02);
    p1.SetTicks(1,1);
    p1.SetLogy(1);
    frame1 = p1.DrawFrame(0,6e-8,10,1e+0);
    frame1.GetXaxis().SetTitle("#it{p}_{T,#gamma} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T,#gamma}} (GeV/#it{c})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.25);
    frame1.GetXaxis().SetTitleSize(0.05);
    frame1.GetYaxis().SetTitleSize(0.05);
    frame1.GetXaxis().SetLabelSize(0.05);
    frame1.GetYaxis().SetLabelSize(0.05);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    ROOT.SetOwnership(frame1,False);
    #frame1.GetXaxis().SetMoreLogLabels(True);
    #frame1.GetYaxis().SetMaxDigits(3);

    arr_pt = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6,7,8,9,10], dtype=float);

    list_h1pt = [];

    for ic in range(0, len(cutnames)):
        cutname = cutnames[ic];
        list_v0_cut = list_v0.FindObject(cutname);
        h1pt_org = list_v0_cut.FindObject("hPt").Clone("hPt_{0}_org".format(cutname));
        h1pt = rebin_histogram(h1pt_org, arr_pt, True, False);
        make_common_style(h1pt  , 20, 1.2, colors[ic], 1, 0);
        h1pt.Scale(1/nev);
        h1pt.Draw("E0,same");
        ROOT.SetOwnership(h1pt_org, False);
        ROOT.SetOwnership(h1pt, False);
        h1pt_org.SetDirectory(0);
        h1pt.SetDirectory(0);
        list_h1pt.append(h1pt);

    txt = TPaveText(0.65,0.70,0.9,0.94,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica bold
    txt.SetTextSize(0.052);
    txt.AddText("Work in Progress");
    txt.AddText("pp, #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("PCM (#gamma #rightarrow e^{+}e^{#minus})");
    txt.AddText("|#it{#eta}_{#gamma}| < 0.9");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    leg = TLegend(0.3,0.58,0.5,0.68);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.052);
    leg.AddEntry(list_h1pt[0], "#gamma candidates with ITS-TPC matched tracks", "P");
    leg.AddEntry(list_h1pt[1], "#gamma candidates with TPConly tracks", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);


    p2 = c1.cd(2);
    p2.SetPad(0, 0, 1 ,0.35);
    p2.SetMargin(0.14,0.02,0.25,0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0,1.,10, 6.5);
    frame2.GetXaxis().SetTitle("#it{p}_{T,#gamma} (GeV/#it{c})");
    frame2.GetYaxis().SetTitle("#frac{TPConly}{ITS-TPC}");
    frame2.GetXaxis().SetTitleOffset(1.1);
    frame2.GetYaxis().SetTitleOffset(0.6);
    frame2.GetXaxis().SetTitleSize(0.10);
    frame2.GetYaxis().SetTitleSize(0.10);
    frame2.GetXaxis().SetLabelSize(0.10);
    frame2.GetYaxis().SetLabelSize(0.10);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().CenterTitle(True);
    ROOT.SetOwnership(frame2,False);

    h1ratio = list_h1pt[1].Clone("h1ratio");
    h1ratio.Reset();
    h1ratio.Sumw2();
    h1ratio.Divide(list_h1pt[1], list_h1pt[0], 1., 1., "G");
    h1ratio.Draw("E0same");
    h1ratio.SetDirectory(0);
    ROOT.SetOwnership(h1ratio, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_inc_gamma_{2}{3}.eps".format(date, period, cutname_tmp, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_inc_gamma_{2}{3}.pdf".format(date, period, cutname_tmp, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_inc_gamma_{2}{3}.png".format(date, period, cutname_tmp, suffix));

    rootfile.Close();

#_____________________________________________________________________
def draw_mgg_eta(filename, subsystem, period, cutname, pt_id=3, suffix=""):
    rootfile = TFile.Open(filename,"READ");
    list_ss  = rootfile.Get(subsystem);
    list_ss_cut = list_ss.FindObject(cutname);
    list_ss_cut_fit           = list_ss_cut  .FindObject("cb_pol1");
    list_ss_cut_fit_range     = list_ss_cut_fit  .FindObject("fit_0.35_0.75_GeVc2");
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
    c1.SetMargin(0.10,0.04,0.1,0.04);

    h1same.GetXaxis().SetRangeUser(0.4,0.7);
    ymax = h1same.GetMaximum() * 1.5;
    #ymax = h1same.GetMaximum() * 2.0;#only for PHOS
    ymin = h1same.GetMaximum() * -0.05;
    mw = h1same.GetBinWidth(1);
    h1same.GetXaxis().SetRangeUser(0,1);

    frame1 = c1.DrawFrame(0.4,ymin,0.7,ymax);
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

    make_common_style(h1same  , 20, 1.2, kBlack, 2, 0);
    make_common_style(h1bkg   , 24, 1.2, kBlue+1, 2, 0);
    make_common_style(h1sig   , 21, 1.2, kRed+1, 2, 0);

    h1same.Draw("E0same");
    h1bkg. Draw("E0same");
    h1sig. Draw("E0same");
    f1sig = h1sig.GetFunction("f1sig_pt{0}".format(pt_id));
    f1sig.SetLineColor(kRed+1);
    f1sig.SetLineWidth(2);
    f1sig.SetNpx(1000);

    subsystem_tmp = "";
    if subsystem == "PCMPCM":
        subsystem_tmp = "PCM-PCM";
    elif subsystem == "PHOSPHOS":
        subsystem_tmp = "PHOS-PHOS";
    elif subsystem == "EMCEMC":
        subsystem_tmp = "EMC-EMC";

    leg = TLegend(0.55,0.72,0.7,0.92);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.035);
    leg.SetHeader("#eta #rightarrow #gamma#gamma , {0}".format(subsystem_tmp));
    leg.AddEntry(h1same , "signal + background","P");
    if subsystem == "EMCEMC":
        leg.AddEntry(h1bkg  , "rotation background","P");
    else:
        leg.AddEntry(h1bkg  , "mixed-event background","P");
    leg.AddEntry(h1sig  , "signal","P");
    leg.AddEntry(f1sig  , "fit to signal","L");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.12,0.80,0.4,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    #txt.AddText("ALICE Performance");
    txt.AddText("pp, #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("{0} pass4".format(period));
    #txt.AddText("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem));
    txt.AddText(ptstr[25:]);
    #txt.AddText(suffix);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_mgg_eta_{2}_{3}_Pt{4}{5}.eps".format(date, period, subsystem, cutname, pt_id, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_mgg_eta_{2}_{3}_Pt{4}{5}.pdf".format(date, period, subsystem, cutname, pt_id, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_mgg_eta_{2}_{3}_Pt{4}{5}.png".format(date, period, subsystem, cutname, pt_id, suffix));

    rootfile.Close();

#_____________________________________________________________________
def draw_mgg_pi0_pbpb(filename, subsystem, period, cutname, pt_id=0, suffix=""):
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
    c1.SetMargin(0.10,0.04,0.1,0.04);

    h1same.GetXaxis().SetRangeUser(0.02,0.3);
    ymax = h1same.GetMaximum() * 1.2;
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

    make_common_style(h1same  , 20, 1.2, kBlack , 2, 0);
    make_common_style(h1bkg   , 24, 1.2, kBlue+1, 2, 0);
    make_common_style(h1sig   , 21, 1.2, kRed+1 , 2, 0);

    sf = 1;
    f1sig = h1sig.GetFunction("f1sig_pt{0}".format(pt_id));
    f1sig.SetLineColor(kRed+1);
    f1sig.SetLineWidth(2);
    height = f1sig.GetParameter(0);
    f1sig.SetParameter(0, height * sf);
    mean  = f1sig.GetParameter(1);
    mean_err  = f1sig.GetParError(1);
    sigma = f1sig.GetParameter(2);
    sigma_err = f1sig.GetParError(2);

    bin0 = h1sig.FindBin(mean - 3 * sigma); 
    bin1 = h1sig.FindBin(mean + 3 * sigma); 
    nsig = 0;
    nsig_err = ctypes.c_double(0);
    nsig = h1sig.IntegralAndError(bin0, bin1, nsig_err, "");

    h1sig.Scale(sf);
    h1same.Draw("E0same");
    h1bkg. Draw("E0same");
    h1sig. Draw("E0same");

    leg = TLegend(0.55,0.72,0.7,0.92);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.035);
    #leg.SetHeader("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem));
    leg.SetHeader("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format("PCM-PCM"));
    leg.AddEntry(h1same , "signal + background","P");
    if subsystem == "EMCEMC":
        leg.AddEntry(h1bkg  , "rotated bkg","P");
    else:
        leg.AddEntry(h1bkg  , "mixed-event background","P");
    leg.AddEntry(h1sig  , "signal","P");
    leg.AddEntry(f1sig  , "fit to signal","L");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.12,0.80,0.4,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("0#minus100%");
    #txt.AddText("{0}".format(period));
    txt.AddText(ptstr[25:]);
    #txt.AddText(suffix);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    #txt_sf = TPaveText(0.4,0.25,0.6,0.35,"NDC");
    #txt_sf.SetFillColor(kWhite);
    #txt_sf.SetFillStyle(0);
    #txt_sf.SetBorderSize(0);
    #txt_sf.SetTextAlign(12);#middle,left
    #txt_sf.SetTextFont(42);#helvetica
    #txt_sf.SetTextSize(0.035);
    #txt_sf.SetTextColor(kRed+1);
    #txt_sf.AddText("signal #times {0:2.1f} for visibility".format(sf));
    #txt_sf.Draw();
    #ROOT.SetOwnership(txt_sf, False);

    txt_nsig = TPaveText(0.15,0.55,0.4,0.7,"NDC");
    txt_nsig.SetFillColor(kWhite);
    txt_nsig.SetFillStyle(0);
    txt_nsig.SetBorderSize(0);
    txt_nsig.SetTextAlign(12);#middle,left
    txt_nsig.SetTextFont(42);#helvetica
    txt_nsig.SetTextSize(0.035);
    txt_nsig.AddText("N (#pm3#sigma) = {0:3.1f} #pm {1:3.1f}".format(nsig, nsig_err.value));
    txt_nsig.AddText("#mu = {0:2.1f} #pm {1:2.1f} MeV/#it{{c}}^{{2}}".format(mean*1e+3, mean_err*1e+3));
    txt_nsig.AddText("#sigma = {0:2.1f} #pm {1:2.1f} MeV/#it{{c}}^{{2}}".format(sigma*1e+3, sigma_err*1e+3));
    #txt_nsig.Draw();
    ROOT.SetOwnership(txt_nsig, False);

    #txt_nsig = TPaveText(0.15,0.55,0.4,0.72,"NDC");
    #txt_nsig.SetFillColor(kWhite);
    #txt_nsig.SetFillStyle(0);
    #txt_nsig.SetBorderSize(0);
    #txt_nsig.SetTextAlign(12);#middle,left
    #txt_nsig.SetTextFont(42);#helvetica
    #txt_nsig.SetTextSize(0.035);
    #txt_nsig.AddText("N_{{#pi^{{0}}}} (#pm3#sigma) = {0:3.1f} #pm {1:3.1f}".format(nsig, nsig_err.value));
    #txt_nsig.AddText("#mu_{{#pi^{{0}}}} = {0:2.1f} #pm {1:2.1f} MeV/#it{{c}}^{{2}}".format(mean*1e+3, mean_err*1e+3));
    #txt_nsig.AddText("#sigma_{{#pi^{{0}}}} = {0:2.1f} #pm {1:2.1f} MeV/#it{{c}}^{{2}}".format(sigma*1e+3, sigma_err*1e+3));
    #txt_nsig.Draw();
    #ROOT.SetOwnership(txt_nsig, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_PbPb_5.36TeV_{1}_mgg_pi0_{2}_{3}_Pt{4}{5}_wN.eps".format(date, period, subsystem, cutname, pt_id, suffix));
    c1.SaveAs("{0}_PbPb_5.36TeV_{1}_mgg_pi0_{2}_{3}_Pt{4}{5}_wN.pdf".format(date, period, subsystem, cutname, pt_id, suffix));
    c1.SaveAs("{0}_PbPb_5.36TeV_{1}_mgg_pi0_{2}_{3}_Pt{4}{5}_wN.png".format(date, period, subsystem, cutname, pt_id, suffix));

    rootfile.Close();

#_____________________________________________________________________
def draw_mgg_pi0(filename, subsystem, period, cutname, pt_id=3, suffix=""):
    rootfile = TFile.Open(filename,"READ");
    list_ss  = rootfile.Get(subsystem);
    list_ss_cut = list_ss.FindObject(cutname);
    list_ss_cut_fit           = list_ss_cut  .FindObject("cb_pol1");
    list_ss_cut_fit_range     = list_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    #list_ss_cut_fit_range     = list_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
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
    c1.SetMargin(0.10,0.04,0.1,0.04);

    h1same.GetXaxis().SetRangeUser(0.04,0.24);
    ymax = h1same.GetMaximum() * 1.5;
    #ymax = h1same.GetMaximum() * 2.0;
    ymin = h1same.GetMaximum() * -0.05;
    mw = h1same.GetBinWidth(1);
    #ymax = 1700;
    #ymin = -50;

    #frame1 = c1.DrawFrame(0.06,ymin,0.22,ymax);
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

    make_common_style(h1same  , 20, 1.2, kBlack, 2, 0);
    make_common_style(h1bkg   , 24, 1.2, kBlue+1, 2, 0);
    make_common_style(h1sig   , 21, 1.2, kRed+1, 2, 0);

    f1sig = h1sig.GetFunction("f1sig_pt{0}".format(pt_id));
    f1sig.SetLineColor(kRed+1);
    f1sig.SetLineWidth(2);
    f1sig.SetNpx(1000);
    #height = f1sig.GetParameter(0);
    #f1sig.SetParameter(0, height * 2);
    #h1sig.Scale(2);
    mean  = f1sig.GetParameter(1);
    mean_err  = f1sig.GetParError(1);
    sigma = f1sig.GetParameter(2);
    sigma_err = f1sig.GetParError(2);

    bin0 = h1sig.FindBin(mean - 3 * sigma); 
    bin1 = h1sig.FindBin(mean + 3 * sigma); 
    nsig = 0;
    nsig_err = ctypes.c_double(0);
    nsig = h1sig.IntegralAndError(bin0, bin1, nsig_err, "");

    nbkg = 0;
    nbkg_err = ctypes.c_double(0);
    nbkg = h1bkg.IntegralAndError(bin0, bin1, nbkg_err, "");

    h1same.Draw("E0same");
    h1bkg. Draw("E0same");
    h1sig. Draw("E0same");

    subsystem_tmp = "";
    if subsystem == "PCMPCM":
        subsystem_tmp = "PCM-PCM";
    elif subsystem == "PHOSPHOS":
        subsystem_tmp = "PHOS-PHOS";
    elif subsystem == "EMCEMC":
        subsystem_tmp = "EMC-EMC";

    leg = TLegend(0.55,0.72,0.7,0.92);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.035);

    if "DalitzEE" in subsystem :
        leg.SetHeader("#pi^{0} #rightarrow ee#gamma");
    else:
        leg.SetHeader("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem_tmp));
    leg.AddEntry(h1same , "signal + background","P");
    if subsystem == "EMCEMC":
        leg.AddEntry(h1bkg  , "rotation background","P");
    else:
        leg.AddEntry(h1bkg  , "mixed-event background","P");
    leg.AddEntry(h1sig  , "signal","P");
    leg.AddEntry(f1sig  , "fit to signal","L");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.12,0.80,0.4,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    #txt.AddText("ALICE Performance");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("{0} pass4".format(period));
    #txt.AddText("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem));
    txt.AddText(ptstr[25:]);
    #txt.AddText(suffix);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    #txt_nsig = TPaveText(0.15,0.6,0.4,0.7,"NDC");
    #txt_nsig.SetFillColor(kWhite);
    #txt_nsig.SetFillStyle(0);
    #txt_nsig.SetBorderSize(0);
    #txt_nsig.SetTextAlign(12);#middle,left
    #txt_nsig.SetTextFont(42);#helvetica
    #txt_nsig.SetTextSize(0.035);
    #txt_nsig.AddText("S (#pm3#sigma) = {0:3.1f} #pm {1:3.1f}".format(nsig, nsig_err.value));
    #txt_nsig.AddText("B (#pm3#sigma) = {0:3.1f} #pm {1:3.1f}".format(nbkg, nbkg_err.value));
    ##txt_nsig.AddText("S/B = {0:4.2f}".format(nsig/nbkg));
    ##txt_nsig.AddText("S/#sqrt{{S+2B}} = {0:4.2f}".format(nsig/math.sqrt(nsig + 2*nbkg)));
    ##txt_nsig.AddText("#mu = {0:2.1f} #pm {1:2.1f} MeV/#it{{c}}^{{2}}".format(mean*1e+3, mean_err*1e+3));
    ##txt_nsig.AddText("#sigma = {0:2.1f} #pm {1:2.1f} MeV/#it{{c}}^{{2}}".format(sigma*1e+3, sigma_err*1e+3));
    #txt_nsig.Draw();
    #ROOT.SetOwnership(txt_nsig, False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_mgg_pi0_{2}_{3}_Pt{4}{5}.eps".format(date, period, subsystem, cutname, pt_id, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_mgg_pi0_{2}_{3}_Pt{4}{5}.pdf".format(date, period, subsystem, cutname, pt_id, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_mgg_pi0_{2}_{3}_Pt{4}{5}.png".format(date, period, subsystem, cutname, pt_id, suffix));

    rootfile.Close();

#_____________________________________________________________________
def draw_raw_yield_run2(subsystem, period, cutname, suffix=""):
    rootfile_pi0_run2 = TFile.Open("data_PCMResultsFullCorrection_PP.root","READ");
    dir_pi0_run2  = rootfile_pi0_run2.Get("Pi013TeV");
    h1_pi0_run2    = dir_pi0_run2 .Get("RAWYieldPerEventsPi0_INT7");
    h1_pi0_run2.SetDirectory(0);
    ROOT.SetOwnership(h1_pi0_run2, False);

    rootfile_eta_run2 = TFile.Open("data_PCMResultsFullCorrection_PP.root","READ");
    dir_eta_run2  = rootfile_eta_run2.Get("Eta13TeV");
    h1_eta_run2    = dir_eta_run2 .Get("RAWYieldPerEventsEta_INT7");
    h1_eta_run2.SetDirectory(0);
    ROOT.SetOwnership(h1_eta_run2, False);


    rootfile_pi0_run3 = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_20240410.root".format(period),"READ");
    list_pi0_run3 = rootfile_pi0_run3.Get(subsystem);
    list_pi0_run3_cut               = list_pi0_run3  .FindObject(cutname);
    list_pi0_run3_cut_fit           = list_pi0_run3_cut  .FindObject("cb_pol1");
    list_pi0_run3_cut_fit_range     = list_pi0_run3_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_pi0_run3_cut_fit_range_int = list_pi0_run3_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1_pi0_run3   = list_pi0_run3_cut_fit_range_int.FindObject("h1yield");
    h1_pi0_run3.SetDirectory(0);
    ROOT.SetOwnership(h1_pi0_run3, False);

    rootfile_eta_run3 = TFile.Open("eta_data_ptspectrum_pp_13.6TeV_{0}_20240410.root".format(period),"READ");
    list_eta_run3 = rootfile_eta_run3.Get(subsystem);
    list_eta_run3_cut               = list_eta_run3  .FindObject(cutname);
    list_eta_run3_cut_fit           = list_eta_run3_cut  .FindObject("cb_pol2");
    list_eta_run3_cut_fit_range     = list_eta_run3_cut_fit  .FindObject("fit_0.40_0.70_GeVc2");
    list_eta_run3_cut_fit_range_int = list_eta_run3_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1_eta_run3   = list_eta_run3_cut_fit_range_int.FindObject("h1yield");
    h1_eta_run3.SetDirectory(0);
    ROOT.SetOwnership(h1_eta_run3, False);

 
    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetLogx(1);
    c1.SetLogy(1);
    c1.SetTicks(1,1);
    c1.SetMargin(0.16,0.03,0.12,0.03);
    frame1 = c1.DrawFrame(0.4,1e-8,10,1e-2);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    #frame1.GetXaxis().SetLabelOffset(0.01);
    #frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(True);

    make_common_style(h1_pi0_run3  , 20, 1.0, kRed+1, 1, 0);
    make_common_style(h1_pi0_run2  , 24, 1.0, kRed+1, 1, 0);
    make_common_style(h1_eta_run3  , 20, 1.0, kBlue+1, 1, 0);
    make_common_style(h1_eta_run2  , 24, 1.0, kBlue+1, 1, 0);

    for i in range(0, h1_pi0_run3.GetNbinsX()):
        if h1_pi0_run3.GetBinCenter(i+1) > 10:
            h1_pi0_run3.SetBinContent(i+1, 0);
            h1_pi0_run3.SetBinError(i+1, 0);

    for i in range(0, h1_eta_run3.GetNbinsX()):
        if h1_eta_run3.GetBinCenter(i+1) > 6:
            h1_eta_run3.SetBinContent(i+1, 0);
            h1_eta_run3.SetBinError(i+1, 0);


    h1_pi0_run2.Draw("E0same");
    h1_pi0_run3.Draw("E0same");
    h1_eta_run2.Draw("E0same");
    h1_eta_run3.Draw("E0same");

    txt = TPaveText(0.2,0.8,0.4,0.95,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("{0} pass6".format(period));
    txt.AddText("#pi^{{0}}/#eta #rightarrow #gamma#gamma , {0}".format(subsystem));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    leg = TLegend(0.6,0.8,1.0,0.95);
    leg.SetBorderSize(0);
    leg.SetNColumns(2);
    leg.SetFillStyle(0);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.04);
    leg.AddEntry(h1_pi0_run3, "#pi^{0}","LP");
    leg.AddEntry(h1_pi0_run2, "#pi^{0} Run2","LP");
    leg.AddEntry(h1_eta_run3, "#eta","LP");
    leg.AddEntry(h1_eta_run2, "#eta Run2","LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass6_nm_raw_yield_{2}_{3}{4}.eps".format(date, subsystem, period, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass6_nm_raw_yield_{2}_{3}{4}.pdf".format(date, subsystem, period, cutname, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass6_nm_raw_yield_{2}_{3}{4}.png".format(date, subsystem, period, cutname, suffix));

    rootfile_pi0_run3.Close();
    rootfile_pi0_run2.Close();

#_____________________________________________________________________
def draw_raw_yield_pi0_high_low(subsystem, period_low, period_high, cutname):
    #rootfile_lkb_low  = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withLKB.root".format(period_low),"READ");
    #rootfile_lkb_high = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withLKB.root".format(period_high),"READ");
    rootfile_lkb_low  = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withAnyTrack_20230828_nsw1.root".format(period_low),"READ");
    rootfile_lkb_high = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withAnyTrack_20230828_nsw1.root".format(period_high),"READ");
    rootfile_any_low  = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withAnyTrack_20230828_nsw1.root".format(period_low),"READ");
    rootfile_any_high = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withAnyTrack_20230828_nsw1.root".format(period_high),"READ");

    list_pcm_pcm_lkb_low  = rootfile_lkb_low .Get(subsystem);
    list_pcm_pcm_lkb_high = rootfile_lkb_high.Get(subsystem);
    list_pcm_pcm_any_low  = rootfile_any_low .Get(subsystem);
    list_pcm_pcm_any_high = rootfile_any_high.Get(subsystem);

    list_pcm_pcm_lkb_low_cut      = list_pcm_pcm_lkb_low .FindObject(cutname);
    list_pcm_pcm_lkb_high_cut    = list_pcm_pcm_lkb_high.FindObject(cutname);
    list_pcm_pcm_any_low_cut = list_pcm_pcm_any_low .FindObject(cutname);
    list_pcm_pcm_any_high_cut        = list_pcm_pcm_any_high.FindObject(cutname);

    list_pcm_pcm_lkb_low_cut_fit = list_pcm_pcm_lkb_low_cut.FindObject("cb_pol1");
    list_pcm_pcm_lkb_high_cut_fit = list_pcm_pcm_lkb_high_cut.FindObject("cb_pol1");
    list_pcm_pcm_any_low_cut_fit = list_pcm_pcm_any_low_cut.FindObject("cb_pol1");
    list_pcm_pcm_any_high_cut_fit = list_pcm_pcm_any_high_cut.FindObject("cb_pol1");

    list_pcm_pcm_lkb_low_cut_fit_range  = list_pcm_pcm_lkb_low_cut_fit.FindObject("fit_0.04_0.24_GeVc2");
    list_pcm_pcm_lkb_high_cut_fit_range = list_pcm_pcm_lkb_high_cut_fit.FindObject("fit_0.04_0.24_GeVc2");
    list_pcm_pcm_any_low_cut_fit_range  = list_pcm_pcm_any_low_cut_fit.FindObject("fit_0.04_0.24_GeVc2");
    list_pcm_pcm_any_high_cut_fit_range = list_pcm_pcm_any_high_cut_fit.FindObject("fit_0.04_0.24_GeVc2");

    list_pcm_pcm_lkb_low_cut_fit_range_int  = list_pcm_pcm_lkb_low_cut_fit_range.FindObject("integral_-3.0_3.0_sigma");
    list_pcm_pcm_lkb_high_cut_fit_range_int = list_pcm_pcm_lkb_high_cut_fit_range.FindObject("integral_-3.0_3.0_sigma");
    list_pcm_pcm_any_low_cut_fit_range_int  = list_pcm_pcm_any_low_cut_fit_range.FindObject("integral_-3.0_3.0_sigma");
    list_pcm_pcm_any_high_cut_fit_range_int = list_pcm_pcm_any_high_cut_fit_range.FindObject("integral_-3.0_3.0_sigma");

    h1_pcm_pcm_lkb_low  = list_pcm_pcm_lkb_low_cut_fit_range_int .FindObject("h1yield");
    h1_pcm_pcm_lkb_high = list_pcm_pcm_lkb_high_cut_fit_range_int.FindObject("h1yield");
    h1_pcm_pcm_any_low  = list_pcm_pcm_any_low_cut_fit_range_int .FindObject("h1yield");
    h1_pcm_pcm_any_high = list_pcm_pcm_any_high_cut_fit_range_int.FindObject("h1yield");
    h1_pcm_pcm_lkb_low .SetDirectory(0);
    h1_pcm_pcm_lkb_high.SetDirectory(0);
    h1_pcm_pcm_any_low .SetDirectory(0);
    h1_pcm_pcm_any_high.SetDirectory(0);
    ROOT.SetOwnership(h1_pcm_pcm_lkb_low , False);
    ROOT.SetOwnership(h1_pcm_pcm_lkb_high, False);
    ROOT.SetOwnership(h1_pcm_pcm_any_low , False);
    ROOT.SetOwnership(h1_pcm_pcm_any_high, False);
 
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

    make_common_style(h1_pcm_pcm_lkb_low , 20, 1.2, kBlue+1, 2, 0);
    make_common_style(h1_pcm_pcm_lkb_high, 24, 1.2, kBlue+1, 2, 0);
    make_common_style(h1_pcm_pcm_any_low , 20, 1.2, kRed+1, 2, 0);
    make_common_style(h1_pcm_pcm_any_high, 24, 1.2, kRed+1, 2, 0);

    h1_pcm_pcm_lkb_low .Draw("E0same");
    h1_pcm_pcm_lkb_high.Draw("E0same");
    h1_pcm_pcm_any_low .Draw("E0same");
    h1_pcm_pcm_any_high.Draw("E0same");

    leg = TLegend(0.4,0.75,0.6,0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.03);
    #leg.SetHeader("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format(subsystem));
    leg.SetHeader("#pi^{{0}} #rightarrow #gamma#gamma , {0}".format("PCM-PCM"));
    #leg.AddEntry(h1_pcm_pcm_lkb_low , "{0} pass4 15 kHz, V0 with LKB"        .format(period_low),"LP");
    #leg.AddEntry(h1_pcm_pcm_lkb_high, "{0} pass4 500 kHz, V0 with LKB"        .format(period_high),"LP");
    leg.AddEntry(h1_pcm_pcm_any_low , "{0} pass4 15 kHz, V0s with any tracks" .format(period_low),"LP");
    leg.AddEntry(h1_pcm_pcm_any_high, "{0} pass4 500 kHz, V0s with any tracks" .format(period_high),"LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.15,0.9,0.4,0.95,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.Draw();
    ROOT.SetOwnership(txt,False);


    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_raw_yield_{2}_{3}_{4}.eps".format(date, period_low, period_high, subsystem, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_raw_yield_{2}_{3}_{4}.pdf".format(date, period_low, period_high, subsystem, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_raw_yield_{2}_{3}_{4}.png".format(date, period_low, period_high, subsystem, cutname));

    rootfile_lkb_low .Close();
    rootfile_lkb_high .Close();
    rootfile_any_low .Close();
    rootfile_any_high .Close();
#_____________________________________________________________________
def draw_raw_yield_pi0(subsystem, period, cutname):
    #rootfile_loose      = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withAnyTracks.root".format(period),"READ");
    rootfile_tpconly    = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withTPConlyTracks.root".format(period),"READ");
    rootfile_tpconlybar = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withTPConlyTracksBar.root".format(period),"READ");
    #rootfile_lkb        = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withLKB.root".format(period),"READ");
    rootfile_run2       = TFile.Open("data_PCMResultsFullCorrection_PP.root","READ");
    rootfile_loose      = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withAnyTrack.root".format("LHC22q"),"READ");
    rootfile_lkb        = TFile.Open("pi0_data_ptspectrum_pp_13.6TeV_{0}_V0withLKB.root".format("LHC22q"),"READ");

    period = "LHC22q";

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
def draw_meeg_pi0(filename, subsystem, period, cutname, pt_id=3, suffix=""):
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

    h1same.GetXaxis().SetRangeUser(0.06,0.24);
    ymax = h1same.GetMaximum() * 1.5;
    ymin = h1same.GetMaximum() * -0.05;
    mw = h1same.GetBinWidth(1);
    h1same.GetXaxis().UnZoom();

    frame1 = c1.DrawFrame(0.04, ymin, 0.24,ymax);
    #frame1.GetXaxis().SetTitle("#it{m}_{ee#gamma} (GeV/#it{c}^{2})");
    frame1.GetXaxis().SetTitle("#it{m}_{#gamma#gamma*} (GeV/#it{c}^{2})");
    #frame1.GetXaxis().SetTitle("#it{m}_{eeee} (GeV/#it{c}^{2})");
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
    ROOT.SetOwnership(frame1,False);

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
    if "Dalitz" in subsystem:
        leg.SetHeader("#pi^{{0}} #rightarrow #gamma#gamma* , {0}".format(subsystem));
    else:
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
    txt.AddText("{0} pass1".format(period));
    txt.AddText(ptstr[25:].replace("#gamma#gamma","#gamma"));
    if "Dalitz" in subsystem:
        txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(0.0, 0.12));
    #txt.AddText(suffix);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_meeg_pi0_{2}_{3}_Pt{4}{5}.eps".format(date, period, subsystem, cutname, pt_id, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_meeg_pi0_{2}_{3}_Pt{4}{5}.pdf".format(date, period, subsystem, cutname, pt_id, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_meeg_pi0_{2}_{3}_Pt{4}{5}.png".format(date, period, subsystem, cutname, pt_id, suffix));

    rootfile.Close();

#_____________________________________________________________________
def draw_tagging(filename_inc, ssname_inc, cutname_inc,
        filename_tag, ssname1_tag, cutname1_tag, ssname2_tag, cutname2_tag, period, suffix=""):
    print("reading...", filename_inc, filename_tag);

    rootfile_inc = TFile.Open(filename_inc, "READ");
    rootlist_inc_ss = rootfile_inc.Get(ssname_inc);
    rootlist_inc_ss_cut = rootlist_inc_ss.FindObject(cutname_inc);
    h1dndpt_inc = rootlist_inc_ss_cut.FindObject("h1yield");
    h1dndpt_inc.SetDirectory(0);
    ROOT.SetOwnership(h1dndpt_inc,False);
    make_common_style(h1dndpt_inc  , 20, 1.0, kRed+1, 1, 0);

    funcname = "cb_pol1";
    fitname = "fit_0.04_0.24_GeVc2";
    intname = "integral_-3.0_3.0_sigma";

    rootfile_tag = TFile.Open(filename_tag, "READ");
    rootlist_tag_ss1 = rootfile_tag.Get(ssname1_tag);
    rootlist_tag_ss1_cut1 = rootlist_tag_ss1.FindObject(cutname1_tag);
    h1dndpt_tag1 = rootlist_tag_ss1_cut1.FindObject(funcname).FindObject(fitname).FindObject(intname).FindObject("h1yield");
    h1dndpt_tag1.SetDirectory(0);
    ROOT.SetOwnership(h1dndpt_tag1,False);
    make_common_style(h1dndpt_tag1, 20, 1.0, kBlue+1, 1, 0);

    rootlist_tag_ss2 = rootfile_tag.Get(ssname2_tag);
    rootlist_tag_ss2_cut2 = rootlist_tag_ss2.FindObject(cutname2_tag);
    h1dndpt_tag2 = rootlist_tag_ss2_cut2.FindObject(funcname).FindObject(fitname).FindObject(intname).FindObject("h1yield");
    h1dndpt_tag2.SetDirectory(0);
    ROOT.SetOwnership(h1dndpt_tag2,False);
    ROOT.SetOwnership(h1dndpt_tag2,False);
    make_common_style(h1dndpt_tag2, 20, 1.0, kGreen+2, 1, 0);

    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetTicks(1,1);
    c1.SetLogy(1);
    c1.SetMargin(0.15,0.02,0.1,0.02);
    frame1 = c1.DrawFrame(0., 1e-7, 6, 1);
    frame1.GetXaxis().SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T,ee}} (GeV/#it{c})^{-1}");
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.9);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(True);
    frame1.GetYaxis().SetMaxDigits(3);

    h1dndpt_inc.Draw("E0same");
    h1dndpt_tag1.Draw("E0same");
    h1dndpt_tag2.Draw("E0same");

    #txt = TPaveText(0.2,0.85,0.4,0.95,"NDC");
    #txt.SetFillColor(kWhite);
    #txt.SetFillStyle(0);
    #txt.SetBorderSize(0);
    #txt.SetTextAlign(12);#middle,left
    #txt.SetTextFont(42);#helvetica
    #txt.SetTextSize(0.035);
    #txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    #txt.AddText("{0} pass4".format(period));
    #txt.Draw();
    #ROOT.SetOwnership(txt,False);

    leg = TLegend(0.45,0.75,0.65,0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetTextSize(0.035);
    leg.SetHeader("pp at #sqrt{#it{s}} = 13.6 TeV");
    leg.AddEntry(h1dndpt_inc , "#gamma^{{inc}} candidates with {0}".format(ssname_inc),"LP");
    leg.AddEntry(h1dndpt_tag1, "#gamma^{{tagged #pi^{{0}}}} with {0}".format(ssname1_tag),"LP");
    leg.AddEntry(h1dndpt_tag2, "#gamma^{{tagged #pi^{{0}}}} with {0}".format(ssname2_tag),"LP");
    #leg.AddEntry(h1dndpt_tag1, "#gamma^{{tagged}} from #pi^{{0}} with {0}".format(ssname1_tag),"LP");
    #leg.AddEntry(h1dndpt_tag2, "#gamma^{{tagged}} from #pi^{{0}} with {0}".format(ssname2_tag),"LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_{2}_{3}{4}.eps".format(date, period, ssname_inc, cutname_inc, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_{2}_{3}{4}.pdf".format(date, period, ssname_inc, cutname_inc, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_{2}_{3}{4}.png".format(date, period, ssname_inc, cutname_inc, suffix));

    rootfile_inc.Close();
    rootfile_tag.Close();
#_____________________________________________________________________
def draw_hbt_1d_qinv(filename, ssname, period, cutname, ktid, suffix):

    rootfile = TFile.Open(filename, "READ");
    list_ss  = rootfile.Get(ssname);
    list_ss_cut = list_ss.FindObject(cutname);
    list_ss_cut_fit     = list_ss_cut.FindObject("fit_0.00_0.08_GeVc2");
    #list_ss_cut_fit.Print();

    h1 = list_ss_cut_fit.FindObject('h1qinv_cf_kt{0}'.format(ktid));
    h1.SetDirectory(0);
    ROOT.SetOwnership(h1,False);
    ktstr = h1.GetTitle();
    print(ktstr);
    make_common_style(h1, 20, 1.0, kBlack, 1, 0);
    ymax = h1.GetMaximum() * 1.5;
    ymin = 1 - h1.GetMaximum() * 0.2;

    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetTicks(1,1);
    c1.SetMargin(0.1,0.03,0.1,0.02);
    frame1 = c1.DrawFrame(0.,ymin, 0.3, ymax);
    frame1.GetXaxis().SetTitle("#it{q}_{inv} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("C(#it{q})");
    frame1.GetXaxis().SetTitleSize(0.035);
    frame1.GetYaxis().SetTitleSize(0.035);
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.3);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    h1.Draw("E0same");

    txt = TPaveText(0.15,0.77,0.4,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText(ktstr[6:]);
    txt.AddText(ssname);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    f1 = list_ss_cut_fit.FindObject("f1_cf_kt{0}".format(ktid));
    la = f1.GetParameter(0);
    la_err = f1.GetParError(0);
    R = f1.GetParameter(1);
    R_err = f1.GetParError(1);

    txt_fit = TPaveText(0.5,0.77,0.8,0.92,"NDC");
    txt_fit.SetFillColor(kWhite);
    txt_fit.SetFillStyle(0);
    txt_fit.SetBorderSize(0);
    txt_fit.SetTextAlign(12);#middle,left
    txt_fit.SetTextFont(42);#helvetica
    txt_fit.SetTextSize(0.035);
    txt_fit.AddText("C(#it{q}) = 1 + #it{#lambda}_{inv} exp(#minus #it{R}_{inv}^{2} #it{q}_{inv}^{2})");
    txt_fit.AddText("#lambda_{{inv}} = {0:3.2f} #pm {1:3.2f}".format(la, la_err));
    txt_fit.AddText("R_{{inv}} = {0:3.2f} #pm {1:3.2f} (fm)".format(R, R_err));
    txt_fit.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_{2}_{3}_hbt_qinv_kT{4}{5}.eps".format(date, period, ssname, cutname, ktid, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_{2}_{3}_hbt_qinv_kT{4}{5}.pdf".format(date, period, ssname, cutname, ktid, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_{2}_{3}_hbt_qinv_kT{4}{5}.png".format(date, period, ssname, cutname, ktid, suffix));

    rootfile.Close();

#_____________________________________________________________________
def draw_hbt_3d_qosl(filename, ssname, period, cutname, ktid, qmin, qmax, suffix):

    rootfile = TFile.Open(filename, "READ");
    list_ss  = rootfile.Get(ssname);
    list_ss_cut = list_ss.FindObject(cutname);
    list_ss_cut_fit     = list_ss_cut.FindObject("fit_0.00_0.08_GeVc2");
    #list_ss_cut_fit.Print();

    h3 = list_ss_cut_fit.FindObject('h3_cf_kt{0}'.format(ktid));
    h3.SetDirectory(0);
    ROOT.SetOwnership(h3,False);

    bin0 = h3.GetXaxis().FindBin(qmin + 1e-6);
    bin1 = h3.GetXaxis().FindBin(qmax - 1e-6);
    nbin = bin1 - bin0 + 1; #note : nbin is the same for x,y,z
    print("nbin = ",nbin);

    h1out = h3.ProjectionY("h1out", bin0, bin1, bin0, bin1, "");
    h1out.SetDirectory(0);
    h1out.Scale(1/nbin/nbin);
    ROOT.SetOwnership(h1out,False);
    print("y = ", h1out.GetXaxis().GetTitle());

    h1side = h3.ProjectionZ("h1side", bin0, bin1, bin0, bin1, "");
    h1side.SetDirectory(0);
    h1side.Scale(1/nbin/nbin);
    ROOT.SetOwnership(h1side,False);
    print("z = ", h1side.GetXaxis().GetTitle());

    h1long = h3.ProjectionX("h1long", bin0, bin1, bin0, bin1, "");
    h1long.SetDirectory(0);
    h1long.Scale(1/nbin/nbin);
    ROOT.SetOwnership(h1long,False);
    print("x = ", h1long.GetXaxis().GetTitle());

    ktstr = h3.GetTitle();
    print(ktstr);
    #make_common_style(h1, 20, 1.0, kBlack, 1, 0);
    #ymax = h1.GetMaximum() * 1.5;
    #ymin = 1 - h1.GetMaximum() * 0.2;
    ymax =3;
    ymin =0.5;
    make_common_style(h1out , 20, 1.0, kBlack, 1, 0);
    make_common_style(h1side, 20, 1.0, kBlack, 1, 0);
    make_common_style(h1long, 20, 1.0, kBlack, 1, 0);

    list_h1 = [h1out, h1side, h1long];
    list_axis = ["#it{q}_{out} (GeV/#it{c})","#it{q}_{side} (GeV/#it{c})", "#it{q}_{long} (GeV/#it{c})"];

    f3 = list_ss_cut_fit.FindObject("f3_cf_kt{0}".format(ktid));
    h3_fit = f3.CreateHistogram();
    print(h3_fit.GetName());
    h3_fit.SetName("h3_fit");

    bin0 = h3_fit.GetXaxis().FindBin(qmin + 1e-6);
    bin1 = h3_fit.GetXaxis().FindBin(qmax - 1e-6);
    nbin = bin1 - bin0 + 1; #note : nbin is the same for x,y,z
    print("nbin_fit = ",nbin);

    print(h3_fit.GetNbinsX(), h3_fit.GetXaxis().GetXmin(), h3_fit.GetXaxis().GetXmax());
    print(h3_fit.GetNbinsY(), h3_fit.GetYaxis().GetXmin(), h3_fit.GetYaxis().GetXmax());
    print(h3_fit.GetNbinsZ(), h3_fit.GetZaxis().GetXmin(), h3_fit.GetZaxis().GetXmax());

    for ix in range(0, h3_fit.GetNbinsX()):
        for iy in range(0, h3_fit.GetNbinsY()):
            for iz in range(0, h3_fit.GetNbinsZ()):
                x = h3_fit.GetXaxis().GetBinCenter(ix+1);
                y = h3_fit.GetYaxis().GetBinCenter(iy+1);
                z = h3_fit.GetZaxis().GetBinCenter(iz+1);
                cf = f3.Eval(x,y,z);
                h3_fit.SetBinContent(ix+1, iy+1, iz+1, cf);
                h3_fit.SetBinError(ix+1, iy+1, iz+1, 0);
    h1out_fit  = h3_fit.ProjectionY("h1out_fit" , bin0, bin1, bin0, bin1, "");
    h1side_fit = h3_fit.ProjectionZ("h1side_fit", bin0, bin1, bin0, bin1, "");
    h1long_fit = h3_fit.ProjectionX("h1long_fit", bin0, bin1, bin0, bin1, "");

    h1out_fit.SetDirectory(0);
    h1out_fit.Scale(1/nbin/nbin);
    ROOT.SetOwnership(h1out_fit, False);

    h1side_fit.SetDirectory(0);
    h1side_fit.Scale(1/nbin/nbin);
    ROOT.SetOwnership(h1side_fit, False);

    h1long_fit.SetDirectory(0);
    h1long_fit.Scale(1/nbin/nbin);
    ROOT.SetOwnership(h1long_fit, False);

    make_common_style(h1out_fit , 20, 1.0, kRed+1, 1, 0);
    make_common_style(h1side_fit, 20, 1.0, kRed+1, 1, 0);
    make_common_style(h1long_fit, 20, 1.0, kRed+1, 1, 0);
    h1out_fit.SetLineWidth(1);
    h1side_fit.SetLineWidth(1);
    h1long_fit.SetLineWidth(1);
    list_h1_fit = [h1out_fit, h1side_fit, h1long_fit];

    c1 = TCanvas("c0","c0",0,0,1800,600);
    c1.SetMargin(0.,0.,0.,0.);
    c1.Divide(3, 1, 1e-6, 1e-6, 0);
    for i in range(0,3):
        p1 = c1.cd(i+1);
        p1.SetTicks(1,1);
        if i == 0:
            p1.SetPad(0,0,0.36,1);
            p1.SetMargin(0.12,0.0,0.12,0.03);
        elif i == 1:
            p1.SetPad(0.36,0,0.68,1);
            p1.SetMargin(0.0,0.0,0.12,0.03);
        elif i == 2:
            p1.SetPad(0.68,0,1,1);
            p1.SetMargin(0.0,0.02,0.12,0.03);
        frame1 = p1.DrawFrame(-0.28,ymin, +0.28, ymax);
        frame1.GetXaxis().SetTitle(list_axis[i]);
        frame1.GetYaxis().SetTitle("C(#it{q})");
        frame1.GetXaxis().SetTitleSize(0.05);
        frame1.GetYaxis().SetTitleSize(0.05);
        frame1.GetXaxis().SetTitleOffset(1.0);
        frame1.GetYaxis().SetTitleOffset(1.2);
        frame1.GetXaxis().SetLabelSize(0.05);
        frame1.GetYaxis().SetLabelSize(0.05);
        ROOT.SetOwnership(frame1,False);
        ROOT.SetOwnership(p1,False);
        list_h1[i].Draw("E0same");
        list_h1_fit[i].Draw("chist,same");

    c1.cd(1);
    txt = TPaveText(0.15,0.75,0.4,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.045);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText(ktstr[6:]);
    txt.AddText(ssname);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    la = f3.GetParameter(0);
    la_err = f3.GetParError(0);

    rout = f3.GetParameter(2);
    rout_err = f3.GetParError(2);
    rside = f3.GetParameter(3);
    rside_err = f3.GetParError(3);
    rlong = f3.GetParameter(1);
    rlong_err = f3.GetParError(1);

    txt_fit = TPaveText(0.15,0.5,0.4,0.75,"NDC");
    txt_fit.SetFillColor(kWhite);
    txt_fit.SetFillStyle(0);
    txt_fit.SetBorderSize(0);
    txt_fit.SetTextAlign(12);#middle,left
    txt_fit.SetTextFont(42);#helvetica
    txt_fit.SetTextSize(0.04);
    txt_fit.AddText("C(#it{q}) = 1 + #it{#lambda} exp(#minus #it{R}_{out}^{2} #it{q}_{out}^{2} #minus #it{R}_{side}^{2} #it{q}_{side}^{2} #minus #it{R}_{long}^{2} #it{q}_{long}^{2})");
    #txt_fit.AddText("C(#it{q}_{out},#it{q}_{side},#it{q}_{long}) = 1 + #it{#lambda}_{inv} exp(#minus #it{R}_{out}^{2} #it{q}_{out}^{2} #minus #it{R}_{side}^{2} #it{q}_{side}^{2} #minus #it{R}_{long}^{2} #it{q}_{long}^{2})");
    txt_fit.AddText("#lambda = {0:3.2f} #pm {1:3.2f}".format(la, la_err));
    txt_fit.AddText("#it{{R}}_{{out}} = {0:3.2f} #pm {1:3.2f} (fm)".format(rout, rout_err));
    txt_fit.AddText("#it{{R}}_{{side}} = {0:3.2f} #pm {1:3.2f} (fm)".format(rside, rside_err));
    txt_fit.AddText("#it{{R}}_{{long}} = {0:3.2f} #pm {1:3.2f} (fm)".format(rlong, rlong_err));
    txt_fit.Draw();
    ROOT.SetOwnership(txt,False);

    c1.cd(1);
    txt_range_out = TPaveText(0.5,0.75,0.7,0.92,"NDC");
    txt_range_out.SetFillColor(kWhite);
    txt_range_out.SetFillStyle(0);
    txt_range_out.SetBorderSize(0);
    txt_range_out.SetTextAlign(12);#middle,left
    txt_range_out.SetTextFont(42);#helvetica
    txt_range_out.SetTextSize(0.045);
    txt_range_out.AddText("{0:3.2f} < #it{{q}}_{{side}} < {1:3.2f} GeV/#it{{c}}".format(qmin, qmax));
    txt_range_out.AddText("{0:3.2f} < #it{{q}}_{{long}} < {1:3.2f} GeV/#it{{c}}".format(qmin, qmax));
    txt_range_out.Draw();
    ROOT.SetOwnership(txt_range_out,False);

    c1.cd(2);
    txt_range_side = TPaveText(0.45,0.75,0.7,0.92,"NDC");
    txt_range_side.SetFillColor(kWhite);
    txt_range_side.SetFillStyle(0);
    txt_range_side.SetBorderSize(0);
    txt_range_side.SetTextAlign(12);#middle,left
    txt_range_side.SetTextFont(42);#helvetica
    txt_range_side.SetTextSize(0.045);
    txt_range_side.AddText("{0:3.2f} < #it{{q}}_{{out}} < {1:3.2f} GeV/#it{{c}}".format(qmin, qmax));
    txt_range_side.AddText("{0:3.2f} < #it{{q}}_{{long}} < {1:3.2f} GeV/#it{{c}}".format(qmin, qmax));
    txt_range_side.Draw();
    ROOT.SetOwnership(txt_range_side,False);

    c1.cd(3);
    txt_range_long = TPaveText(0.45,0.75,0.7,0.92,"NDC");
    txt_range_long.SetFillColor(kWhite);
    txt_range_long.SetFillStyle(0);
    txt_range_long.SetBorderSize(0);
    txt_range_long.SetTextAlign(12);#middle,left
    txt_range_long.SetTextFont(42);#helvetica
    txt_range_long.SetTextSize(0.045);
    txt_range_long.AddText("{0:3.2f} < #it{{q}}_{{out}} < {1:3.2f} GeV/#it{{c}}".format(qmin, qmax));
    txt_range_long.AddText("{0:3.2f} < #it{{q}}_{{side}} < {1:3.2f} GeV/#it{{c}}".format(qmin, qmax));
    txt_range_long.Draw();
    ROOT.SetOwnership(txt_range_long,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_{2}_{3}_hbt_3d_qosl_kT{4}_q{5:3.2f}_{6:3.2f}{7}.eps".format(date, period, ssname, cutname, ktid, qmin, qmax, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_{2}_{3}_hbt_3d_qosl_kT{4}_q{5:3.2f}_{6:3.2f}{7}.pdf".format(date, period, ssname, cutname, ktid, qmin, qmax, suffix));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_{2}_{3}_hbt_3d_qosl_kT{4}_q{5:3.2f}_{6:3.2f}{7}.png".format(date, period, ssname, cutname, ktid, qmin, qmax, suffix));

    rootfile.Close();

#_____________________________________________________________________
def compare_eta_peak(filename_data, filename_mc, period_data, period_mc, ssname, cutname):
    rootfile_run3_data = TFile.Open(filename_data, "READ");
    list_run3_data_ss     = rootfile_run3_data.Get(ssname);
    list_run3_data_ss_cut = list_run3_data_ss.FindObject(cutname);
    list_run3_data_ss_cut_fit           = list_run3_data_ss_cut  .FindObject("cb_pol2");
    list_run3_data_ss_cut_fit_range     = list_run3_data_ss_cut_fit  .FindObject("fit_0.40_0.70_GeVc2");
    list_run3_data_ss_cut_fit_range_int = list_run3_data_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1mean_data = list_run3_data_ss_cut_fit_range_int.FindObject("h1mean");
    h1mean_data.SetDirectory(0);
    ROOT.SetOwnership(h1mean_data,False);
    h1sigma_data = list_run3_data_ss_cut_fit_range_int.FindObject("h1sigma");
    h1sigma_data.SetDirectory(0);
    ROOT.SetOwnership(h1sigma_data,False);

    rootfile_run3_mc   = TFile.Open(filename_mc  , "READ");
    list_run3_mc_ss     = rootfile_run3_mc.Get(ssname);
    list_run3_mc_ss_cut = list_run3_mc_ss.FindObject(cutname);
    list_run3_mc_ss_cut_fit           = list_run3_mc_ss_cut  .FindObject("cb_pol2");
    list_run3_mc_ss_cut_fit_range     = list_run3_mc_ss_cut_fit  .FindObject("fit_0.40_0.70_GeVc2");
    list_run3_mc_ss_cut_fit_range_int = list_run3_mc_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1mean_mc = list_run3_mc_ss_cut_fit_range_int.FindObject("h1mean");
    h1mean_mc.SetDirectory(0);
    ROOT.SetOwnership(h1mean_mc,False);
    h1sigma_mc = list_run3_mc_ss_cut_fit_range_int.FindObject("h1sigma");
    h1sigma_mc.SetDirectory(0);
    ROOT.SetOwnership(h1sigma_mc,False);

    h1mean_mc   .Sumw2();
    h1mean_data .Sumw2();
    h1sigma_mc  .Sumw2();
    h1sigma_data.Sumw2();
    h1mean_mc  .Scale(1e+3);
    h1mean_data.Scale(1e+3);
    h1sigma_mc  .Scale(1e+3);
    h1sigma_data.Scale(1e+3);

    c1 = TCanvas("c1","c1",0,0,800,800);
    c1.SetMargin(0.,0.,0.,0.);
    c1.Divide(1,2,1e-6, 1e-6, 0);

    p1 = c1.cd(1);
    p1.SetPad(0,0.55,1,1);
    p1.SetMargin(0.14,0.03,0.,0.04);
    p1.SetTicks(1,1);
    frame1 = p1.DrawFrame(0.,501, 10, 581);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("peak position (MeV/#it{c}^{2})");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(0.85);
    frame1.GetXaxis().SetTitleSize(0.08);
    frame1.GetYaxis().SetTitleSize(0.08);
    frame1.GetXaxis().SetLabelSize(0.08);
    frame1.GetYaxis().SetLabelSize(0.08);
    frame1.GetXaxis().SetMoreLogLabels(1);
    frame1.GetYaxis().CenterTitle(1);
    frame1.GetYaxis().SetNdivisions(510);
    make_common_style(h1mean_data, 20, 1.0, kBlue+1, 1, 0);
    make_common_style(h1mean_mc  , 24, 1.0, kBlue+1, 1, 0);

    h1mean_mc.Draw("E0same");
    h1mean_data.Draw("E0same");

    p2 = c1.cd(2);
    p2.SetPad(0,0,1,0.55);
    p2.SetMargin(0.14,0.03,0.2,0.0);
    p2.SetTicks(1,1);
    frame2 = p2.DrawFrame(0., 0 , 10, 41);
    frame2.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame2.GetYaxis().SetTitle("peak width (MeV/#it{c}^{2})");
    frame2.GetXaxis().SetTitleOffset(1.2);
    frame2.GetYaxis().SetTitleOffset(1.0);
    frame2.GetXaxis().SetTitleSize(0.065);
    frame2.GetYaxis().SetTitleSize(0.065);
    frame2.GetXaxis().SetLabelSize(0.065);
    frame2.GetYaxis().SetLabelSize(0.065);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(1);
    frame2.GetYaxis().SetNdivisions(510);
    make_common_style(h1sigma_data, 20, 1.0, kBlue+1, 1, 0);
    make_common_style(h1sigma_mc  , 24, 1.0, kBlue+1, 1, 0);

    h1sigma_mc.Draw("E0same");
    h1sigma_data.Draw("E0same");

    ROOT.SetOwnership(frame1,False);
    ROOT.SetOwnership(frame2,False);

    c1.cd(1);
    leg = TLegend(0.55,0.6,0.75,0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.08);
    leg.SetHeader("#eta , PCMPCM");
    leg.AddEntry(h1mean_data ,"Data {0}".format(period_data) ,"LP");
    leg.AddEntry(h1mean_mc   ,"M.C. {0}".format(period_mc)   ,"LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.2,0.8,0.4,0.9,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.08);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_{2}_eta_peak_param_{3}_{4}.eps".format(date, period_data, period_mc, ssname, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_{2}_eta_peak_param_{3}_{4}.pdf".format(date, period_data, period_mc, ssname, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_{2}_eta_peak_param_{3}_{4}.png".format(date, period_data, period_mc, ssname, cutname));

#_____________________________________________________________________
def compare_pi0_peak_rate(filename_low, filename_high, ssname, cutname):
    rootfile_run3_low = TFile.Open(filename_low, "READ");
    list_run3_low_ss     = rootfile_run3_low.Get(ssname);
    list_run3_low_ss_cut = list_run3_low_ss.FindObject(cutname);
    list_run3_low_ss_cut_fit           = list_run3_low_ss_cut  .FindObject("cb_pol1");
    list_run3_low_ss_cut_fit_range     = list_run3_low_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_run3_low_ss_cut_fit_range_int = list_run3_low_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1mean_low = list_run3_low_ss_cut_fit_range_int.FindObject("h1mean");
    h1mean_low.SetDirectory(0);
    ROOT.SetOwnership(h1mean_low,False);
    h1sigma_low = list_run3_low_ss_cut_fit_range_int.FindObject("h1sigma");
    h1sigma_low.SetDirectory(0);
    ROOT.SetOwnership(h1sigma_low,False);

    rootfile_run3_high   = TFile.Open(filename_high  , "READ");
    list_run3_high_ss     = rootfile_run3_high.Get(ssname);
    list_run3_high_ss_cut = list_run3_high_ss.FindObject(cutname);
    list_run3_high_ss_cut_fit           = list_run3_high_ss_cut  .FindObject("cb_pol1");
    list_run3_high_ss_cut_fit_range     = list_run3_high_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_run3_high_ss_cut_fit_range_int = list_run3_high_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1mean_high = list_run3_high_ss_cut_fit_range_int.FindObject("h1mean");
    h1mean_high.SetDirectory(0);
    ROOT.SetOwnership(h1mean_high,False);
    h1sigma_high = list_run3_high_ss_cut_fit_range_int.FindObject("h1sigma");
    h1sigma_high.SetDirectory(0);
    ROOT.SetOwnership(h1sigma_high,False);

    h1mean_high   .Sumw2();
    h1mean_low .Sumw2();
    h1sigma_high  .Sumw2();
    h1sigma_low.Sumw2();
    h1mean_high  .Scale(1e+3);
    h1mean_low.Scale(1e+3);
    h1sigma_high  .Scale(1e+3);
    h1sigma_low.Scale(1e+3);

    c1 = TCanvas("c1","c1",0,0,800,800);
    c1.SetMargin(0.,0.,0.,0.);
    c1.Divide(1,2,1e-6, 1e-6, 0);

    p1 = c1.cd(1);
    p1.SetPad(0,0.55,1,1);
    p1.SetMargin(0.12,0.03,0.,0.04);
    p1.SetTicks(1,1);
    #p1.SetLogx(1);
    frame1 = p1.DrawFrame(0.0,125, 10, 141);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("peak position (MeV/#it{c}^{2})");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(0.68);
    frame1.GetXaxis().SetTitleSize(0.09);
    frame1.GetYaxis().SetTitleSize(0.09);
    frame1.GetXaxis().SetLabelSize(0.09);
    frame1.GetYaxis().SetLabelSize(0.09);
    frame1.GetXaxis().SetMoreLogLabels(1);
    frame1.GetYaxis().CenterTitle(1);
    frame1.GetYaxis().SetNdivisions(510);
    make_common_style(h1mean_low, 20, 1.2, kRed+1, 2, 0);
    make_common_style(h1mean_high  , 21, 1.2, kBlue+1, 2, 0);

    h1mean_high.Draw("E0same");
    h1mean_low.Draw("E0same");

    p2 = c1.cd(2);
    p2.SetPad(0,0,1,0.55);
    p2.SetMargin(0.12,0.03,0.2,0.0);
    p2.SetTicks(1,1);
    #p2.SetLogx(1);
    frame2 = p2.DrawFrame(0.0, 2 , 10, 17);
    frame2.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame2.GetYaxis().SetTitle("peak width (MeV/#it{c}^{2})");
    frame2.GetXaxis().SetTitleOffset(1.2);
    frame2.GetYaxis().SetTitleOffset(0.8);
    frame2.GetXaxis().SetTitleSize(0.075);
    frame2.GetYaxis().SetTitleSize(0.075);
    frame2.GetXaxis().SetLabelSize(0.075);
    frame2.GetYaxis().SetLabelSize(0.075);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(1);
    frame2.GetYaxis().SetNdivisions(510);
    make_common_style(h1sigma_low, 20, 1.2, kRed+1, 2, 0);
    make_common_style(h1sigma_high  , 21, 1.2, kBlue+1, 2, 0);

    h1sigma_high.Draw("E0same");
    h1sigma_low.Draw("E0same");

    ROOT.SetOwnership(frame1,False);
    ROOT.SetOwnership(frame2,False);

    c1.cd(1);
    leg = TLegend(0.55,0.6,0.75,0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.08);
    leg.SetHeader("#pi^{0} , PCM-PCM");
    leg.AddEntry(h1mean_low ,"Data LHC22q, 15 kHz"  ,"LP");
    leg.AddEntry(h1mean_high   ,"Data LHC22o, 500 kHz","LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.2,0.8,0.4,0.9,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.08);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    period = "LHC22q";
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_peak_param_low_high_{2}_{3}.eps".format(date, period, ssname, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_peak_param_low_high_{2}_{3}.pdf".format(date, period, ssname, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pass4_pi0_peak_param_low_high_{2}_{3}.png".format(date, period, ssname, cutname));

#_____________________________________________________________________
def compare_pi0_peak_pbpb(filename_data, filename_mc, ssname, cutname):
    rootfile_run3_data = TFile.Open(filename_data, "READ");
    list_run3_data_ss     = rootfile_run3_data.Get(ssname);
    list_run3_data_ss_cut = list_run3_data_ss.FindObject(cutname);
    list_run3_data_ss_cut_fit           = list_run3_data_ss_cut  .FindObject("cb_pol1");
    list_run3_data_ss_cut_fit_range     = list_run3_data_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_run3_data_ss_cut_fit_range_int = list_run3_data_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1mean_data = list_run3_data_ss_cut_fit_range_int.FindObject("h1mean");
    h1mean_data.SetDirectory(0);
    ROOT.SetOwnership(h1mean_data,False);
    h1sigma_data = list_run3_data_ss_cut_fit_range_int.FindObject("h1sigma");
    h1sigma_data.SetDirectory(0);
    ROOT.SetOwnership(h1sigma_data,False);

    rootfile_run3_mc   = TFile.Open(filename_mc  , "READ");
    list_run3_mc_ss     = rootfile_run3_mc.Get(ssname);
    list_run3_mc_ss_cut = list_run3_mc_ss.FindObject(cutname);
    list_run3_mc_ss_cut_fit           = list_run3_mc_ss_cut  .FindObject("cb_pol1");
    list_run3_mc_ss_cut_fit_range     = list_run3_mc_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_run3_mc_ss_cut_fit_range_int = list_run3_mc_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1mean_mc = list_run3_mc_ss_cut_fit_range_int.FindObject("h1mean");
    h1mean_mc.SetDirectory(0);
    ROOT.SetOwnership(h1mean_mc,False);
    h1sigma_mc = list_run3_mc_ss_cut_fit_range_int.FindObject("h1sigma");
    h1sigma_mc.SetDirectory(0);
    ROOT.SetOwnership(h1sigma_mc,False);

    h1mean_mc   .Sumw2();
    h1mean_data .Sumw2();
    h1sigma_mc  .Sumw2();
    h1sigma_data.Sumw2();
    h1mean_mc  .Scale(1e+3);
    h1mean_data.Scale(1e+3);
    h1sigma_mc  .Scale(1e+3);
    h1sigma_data.Scale(1e+3);

    c1 = TCanvas("c1","c1",0,0,800,800);
    c1.SetMargin(0.,0.,0.,0.);
    c1.Divide(1,2,1e-6, 1e-6, 0);

    p1 = c1.cd(1);
    p1.SetPad(0,0.52,1,1);
    p1.SetMargin(0.12,0.03,0.,0.04);
    p1.SetTicks(1,1);
    #p1.SetLogx(1);
    frame1 = p1.DrawFrame(0.0,127, 10, 141);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("peak position (MeV/#it{c}^{2})");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(0.7);
    frame1.GetXaxis().SetTitleSize(0.08);
    frame1.GetYaxis().SetTitleSize(0.08);
    frame1.GetXaxis().SetLabelSize(0.08);
    frame1.GetYaxis().SetLabelSize(0.08);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    frame1.GetYaxis().CenterTitle(1);
    frame1.GetYaxis().SetNdivisions(510);
    make_common_style(h1mean_data, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1mean_mc  , 24, 1.2, kBlack, 2, 0);

    h1mean_mc.Draw("E0same");
    h1mean_data.Draw("E0same");

    p2 = c1.cd(2);
    p2.SetPad(0,0,1,0.52);
    p2.SetMargin(0.12,0.03,0.2,0.0);
    p2.SetTicks(1,1);
    #p2.SetLogx(1);
    frame2 = p2.DrawFrame(0.0, 2 , 10, 17);
    frame2.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame2.GetYaxis().SetTitle("peak width (MeV/#it{c}^{2})"); #gaussian sigma
    frame2.GetXaxis().SetTitleOffset(1.2);
    frame2.GetYaxis().SetTitleOffset(0.7);
    frame2.GetXaxis().SetTitleSize(0.075);
    frame2.GetYaxis().SetTitleSize(0.075);
    frame2.GetXaxis().SetLabelSize(0.075);
    frame2.GetYaxis().SetLabelSize(0.075);
    frame2.GetXaxis().SetLabelOffset(0.01);
    frame2.GetYaxis().SetLabelOffset(0.01);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(1);
    frame2.GetYaxis().SetNdivisions(510);
    make_common_style(h1sigma_data, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1sigma_mc  , 24, 1.2, kBlack, 2, 0);

    h1sigma_mc.Draw("E0same");
    h1sigma_data.Draw("E0same");

    ROOT.SetOwnership(frame1,False);
    ROOT.SetOwnership(frame2,False);

    c1.cd(1);
    leg = TLegend(0.55,0.67,0.75,0.92);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.07);
    leg.SetHeader("#pi^{0} #rightarrow #gamma#gamma , PCM-PCM");
    leg.AddEntry(h1mean_data ,"Data LHC23zzf_pass2"  ,"LP");
    leg.AddEntry(h1mean_mc   ,"M.C. LHC23k6_TRDfix","LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.15,0.77,0.4,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.07);
    txt.AddText("Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("centrality 0#minus100%");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);    
    period = "LHC23zzf";
    c1.SaveAs("{0}_PbPb_5.36TeV_{1}_pass2_pi0_peak_param_{2}_{3}.eps".format(date, period, ssname, cutname));
    c1.SaveAs("{0}_PbPb_5.36TeV_{1}_pass2_pi0_peak_param_{2}_{3}.pdf".format(date, period, ssname, cutname));
    c1.SaveAs("{0}_PbPb_5.36TeV_{1}_pass2_pi0_peak_param_{2}_{3}.png".format(date, period, ssname, cutname));

#_____________________________________________________________________
def compare_pi0_peak(filename_data, filename_mc, period_data, period_mc, ssname, cutname):
    rootfile_run3_data = TFile.Open(filename_data, "READ");
    list_run3_data_ss     = rootfile_run3_data.Get(ssname);
    list_run3_data_ss_cut = list_run3_data_ss.FindObject(cutname);
    list_run3_data_ss_cut_fit           = list_run3_data_ss_cut  .FindObject("cb_pol1");
    list_run3_data_ss_cut_fit_range     = list_run3_data_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_run3_data_ss_cut_fit_range_int = list_run3_data_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1mean_data = list_run3_data_ss_cut_fit_range_int.FindObject("h1mean");
    h1mean_data.SetDirectory(0);
    ROOT.SetOwnership(h1mean_data,False);
    h1sigma_data = list_run3_data_ss_cut_fit_range_int.FindObject("h1sigma");
    h1sigma_data.SetDirectory(0);
    ROOT.SetOwnership(h1sigma_data,False);

    rootfile_run3_mc   = TFile.Open(filename_mc  , "READ");
    list_run3_mc_ss     = rootfile_run3_mc.Get(ssname);
    list_run3_mc_ss_cut = list_run3_mc_ss.FindObject(cutname);
    list_run3_mc_ss_cut_fit           = list_run3_mc_ss_cut  .FindObject("cb_pol1");
    list_run3_mc_ss_cut_fit_range     = list_run3_mc_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_run3_mc_ss_cut_fit_range_int = list_run3_mc_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1mean_mc = list_run3_mc_ss_cut_fit_range_int.FindObject("h1mean");
    h1mean_mc.SetDirectory(0);
    ROOT.SetOwnership(h1mean_mc,False);
    h1sigma_mc = list_run3_mc_ss_cut_fit_range_int.FindObject("h1sigma");
    h1sigma_mc.SetDirectory(0);
    ROOT.SetOwnership(h1sigma_mc,False);

    h1mean_mc   .Sumw2();
    h1mean_data .Sumw2();
    h1sigma_mc  .Sumw2();
    h1sigma_data.Sumw2();
    h1mean_mc  .Scale(1e+3);
    h1mean_data.Scale(1e+3);
    h1sigma_mc  .Scale(1e+3);
    h1sigma_data.Scale(1e+3);

    c1 = TCanvas("c1","c1",0,0,800,800);
    c1.SetMargin(0.,0.,0.,0.);
    c1.Divide(1,2,1e-6, 1e-6, 0);

    p1 = c1.cd(1);
    p1.SetPad(0,0.55,1,1);
    p1.SetMargin(0.12,0.03,0.,0.04);
    p1.SetTicks(1,1);
    #p1.SetLogx(1);
    frame1 = p1.DrawFrame(0.0,125, 10, 141);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("peak position (MeV/#it{c}^{2})");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(0.68);
    frame1.GetXaxis().SetTitleSize(0.09);
    frame1.GetYaxis().SetTitleSize(0.09);
    frame1.GetXaxis().SetLabelSize(0.09);
    frame1.GetYaxis().SetLabelSize(0.09);
    frame1.GetXaxis().SetMoreLogLabels(1);
    frame1.GetYaxis().CenterTitle(1);
    frame1.GetYaxis().SetNdivisions(510);
    make_common_style(h1mean_data, 20, 1.1, kRed+1, 2, 0);
    make_common_style(h1mean_mc  , 24, 1.1, kRed+1, 2, 0);

    h1mean_mc.Draw("E0same");
    h1mean_data.Draw("E0same");

    p2 = c1.cd(2);
    p2.SetPad(0,0,1,0.55);
    p2.SetMargin(0.12,0.03,0.2,0.0);
    p2.SetTicks(1,1);
    #p2.SetLogx(1);
    frame2 = p2.DrawFrame(0.0, 2 , 10, 17);
    frame2.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame2.GetYaxis().SetTitle("peak width (MeV/#it{c}^{2})"); #gaussian sigma
    frame2.GetXaxis().SetTitleOffset(1.2);
    frame2.GetYaxis().SetTitleOffset(0.8);
    frame2.GetXaxis().SetTitleSize(0.075);
    frame2.GetYaxis().SetTitleSize(0.075);
    frame2.GetXaxis().SetLabelSize(0.075);
    frame2.GetYaxis().SetLabelSize(0.075);
    frame2.GetXaxis().SetMoreLogLabels(1);
    frame2.GetYaxis().CenterTitle(1);
    frame2.GetYaxis().SetNdivisions(510);
    make_common_style(h1sigma_data, 20, 1.2, kRed+1, 2, 0);
    make_common_style(h1sigma_mc  , 24, 1.2, kRed+1, 2, 0);

    h1sigma_mc.Draw("E0same");
    h1sigma_data.Draw("E0same");

    ROOT.SetOwnership(frame1,False);
    ROOT.SetOwnership(frame2,False);

    c1.cd(1);
    leg = TLegend(0.55,0.6,0.75,0.9);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.08);
    leg.SetHeader("#pi^{0} , PCM-PCM");
    leg.AddEntry(h1mean_data ,"Data {0}".format(period_data)  ,"LP");
    leg.AddEntry(h1mean_mc   ,"M.C. {0}".format(period_mc),"LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.2,0.8,0.4,0.9,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.08);
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_{2}_pi0_peak_param_{3}_{4}.eps".format(date, period_data, period_mc, ssname, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_{2}_pi0_peak_param_{3}_{4}.pdf".format(date, period_data, period_mc, ssname, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_{2}_pi0_peak_param_{3}_{4}.png".format(date, period_data, period_mc, ssname, cutname));

#_____________________________________________________________________
def draw_efficiency(filename_mc_pi0, filename_mc_eta, period, ssname, cutname):

    rootfile_pi0 = TFile.Open(filename_mc_pi0, "READ");
    list_pi0_ss     = rootfile_pi0.Get(ssname);
    list_pi0_ss_cut = list_pi0_ss.FindObject(cutname);
    list_pi0_ss_cut_fit           = list_pi0_ss_cut  .FindObject("cb_pol1");
    list_pi0_ss_cut_fit_range     = list_pi0_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_pi0_ss_cut_fit_range_int = list_pi0_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1eff_pi0 = list_pi0_ss_cut_fit_range_int.FindObject("h1eff");
    h1eff_pi0.SetDirectory(0);
    ROOT.SetOwnership(h1eff_pi0,False);

    rootfile_eta = TFile.Open(filename_mc_eta, "READ");
    list_eta_ss     = rootfile_eta.Get(ssname);
    list_eta_ss_cut = list_eta_ss.FindObject(cutname);
    list_eta_ss_cut_fit           = list_eta_ss_cut  .FindObject("cb_pol1");
    list_eta_ss_cut_fit_range     = list_eta_ss_cut_fit  .FindObject("fit_0.40_0.70_GeVc2");
    list_eta_ss_cut_fit_range_int = list_eta_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1eff_eta = list_eta_ss_cut_fit_range_int.FindObject("h1eff");
    h1eff_eta.SetDirectory(0);
    ROOT.SetOwnership(h1eff_eta,False);

    h1eff_pi0.Sumw2();
    h1eff_eta.Sumw2();

    #run2 efficiency
    filename_run2 = "data_PCMResultsFullCorrection_PP.root";
    rootfile_run2 = TFile.Open(filename_run2, "READ");
    dir_pi0_run2 = rootfile_run2.Get("Pi013TeV");
    dir_eta_run2 = rootfile_run2.Get("Eta13TeV");
    h1eff_pi0_run2 = dir_pi0_run2.Get("EffTimesAccPi0_INT7");
    h1eff_eta_run2 = dir_eta_run2.Get("EffTimesAccEta_INT7");
    h1eff_pi0_run2.SetDirectory(0);
    h1eff_eta_run2.SetDirectory(0);
    ROOT.SetOwnership(h1eff_pi0_run2,False);
    ROOT.SetOwnership(h1eff_eta_run2,False);
    h1eff_pi0_run2.Scale(1./TMath.TwoPi());
    h1eff_eta_run2.Scale(1./TMath.TwoPi());
    h1eff_pi0_run2.Scale(1./1.6);
    h1eff_eta_run2.Scale(1./1.6);
    h1eff_pi0_run2.Scale(0.99); #branching ratio
    h1eff_eta_run2.Scale(0.39); #branching ratio

    c1 = TCanvas("c1","c1",0,0,800,800);
    c1.SetMargin(0.11,0.02,0.11,0.02);
    c1.SetTicks(1,1);
    c1.SetLogx(1);
    c1.SetLogy(1);
    frame1 = c1.DrawFrame(0.4,1e-8, 10, 0.1);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("acc. #times rec. efficiency");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.0);
    frame1.GetYaxis().SetLabelOffset(0.0);
    frame1.GetXaxis().SetMoreLogLabels(1);

    make_common_style(h1eff_pi0_run2, 24, 1.0, kRed+1 , 1, 0);
    make_common_style(h1eff_eta_run2, 24, 1.0, kBlue+1, 1, 0);
    h1eff_pi0_run2.Draw("E0same");
    h1eff_eta_run2.Draw("E0same");

    make_common_style(h1eff_pi0, 20, 1.0, kRed+1 , 1, 0);
    make_common_style(h1eff_eta, 20, 1.0, kBlue+1, 1, 0);
    h1eff_pi0.Draw("E0same");
    h1eff_eta.Draw("E0same");


    ROOT.SetOwnership(frame1,False);
    leg = TLegend(0.6,0.82,1.0,0.94);
    leg.SetNColumns(2);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(h1eff_pi0 ,"#pi^{0}" ,"LP");
    leg.AddEntry(h1eff_pi0_run2 ,"#pi^{0} Run2" ,"LP");
    leg.AddEntry(h1eff_eta ,"#eta"    ,"LP");
    leg.AddEntry(h1eff_eta_run2 ,"#eta Run2"    ,"LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.15,0.77,0.4,0.94,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    txt.AddText("ALICE simulation");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText("LHC24b1");
    txt.AddText("anchored to LHC22o apass6");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_efficiency_{2}_{3}.eps".format(date, period, ssname, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_efficiency_{2}_{3}.pdf".format(date, period, ssname, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_efficiency_{2}_{3}.png".format(date, period, ssname, cutname));

#_____________________________________________________________________
def draw_cross_section_eta(filename_data, filename_mc, period, ssname, cutname):

    filename_run2 = "data_PCMResultsFullCorrection_PP.root";
    rootfile_run2 = TFile.Open(filename_run2, "READ");
    dir_pi0_run2 = rootfile_run2.Get("Pi013TeV");
    g1_run2_stat = dir_pi0_run2.Get("graphInvCrossSectionPi0");#in pb GeV-2
    g1_run2_syst = dir_pi0_run2.Get("InvCrossSectionPi0Sys");# in pb GeV-2
    ROOT.SetOwnership(g1_run2_stat,False);
    ROOT.SetOwnership(g1_run2_syst,False);
    np = g1_run2_stat.GetN();
    sf = 1e-12/1e-3;
    for ip in range(0, np):
        y = g1_run2_stat.GetPointY(ip);
        y_stat_high = g1_run2_stat.GetErrorYhigh(ip);
        y_stat_low  = g1_run2_stat.GetErrorYlow(ip);
        y_syst_high = g1_run2_syst.GetErrorYhigh(ip);
        y_syst_low  = g1_run2_syst.GetErrorYlow(ip);

        g1_run2_stat.SetPointY(ip, y * sf);
        g1_run2_stat.SetPointEYhigh(ip, y_stat_high * sf);
        g1_run2_stat.SetPointEYlow(ip, y_stat_low * sf);

        g1_run2_syst.SetPointY(ip, y * sf);
        g1_run2_syst.SetPointEYhigh(ip, y_syst_high * sf);
        g1_run2_syst.SetPointEYlow(ip, y_syst_low * sf);

    #run2 eta
    dir_eta_run2 = rootfile_run2.Get("Eta13TeV");
    g1_eta_run2_stat = dir_eta_run2.Get("graphInvCrossSectionEta");#in pb GeV-2
    g1_eta_run2_syst = dir_eta_run2.Get("InvCrossSectionEtaSys");# in pb GeV-2
    ROOT.SetOwnership(g1_eta_run2_stat,False);
    ROOT.SetOwnership(g1_eta_run2_syst,False);
    np = g1_eta_run2_stat.GetN();
    sf = 1e-12/1e-3;
    for ip in range(0, np):
        y = g1_eta_run2_stat.GetPointY(ip);
        y_stat_high = g1_eta_run2_stat.GetErrorYhigh(ip);
        y_stat_low  = g1_eta_run2_stat.GetErrorYlow(ip);
        y_syst_high = g1_eta_run2_syst.GetErrorYhigh(ip);
        y_syst_low  = g1_eta_run2_syst.GetErrorYlow(ip);

        g1_eta_run2_stat.SetPointY(ip, y * sf);
        g1_eta_run2_stat.SetPointEYhigh(ip, y_stat_high * sf);
        g1_eta_run2_stat.SetPointEYlow(ip, y_stat_low * sf);

        g1_eta_run2_syst.SetPointY(ip, y * sf);
        g1_eta_run2_syst.SetPointEYhigh(ip, y_syst_high * sf);
        g1_eta_run2_syst.SetPointEYlow(ip, y_syst_low * sf);

#    rootfile_run3_data = TFile.Open(filename_data, "READ");
#    list_run3_data_ss     = rootfile_run3_data.Get(ssname);
#    list_run3_data_ss_cut = list_run3_data_ss.FindObject(cutname);
#    list_run3_data_ss_cut_fit           = list_run3_data_ss_cut  .FindObject("cb_pol1");
#    list_run3_data_ss_cut_fit_range     = list_run3_data_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
#    list_run3_data_ss_cut_fit_range_int = list_run3_data_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
#    h1yield = list_run3_data_ss_cut_fit_range_int.FindObject("h1yield");
#    h1yield.SetDirectory(0);
#    ROOT.SetOwnership(h1yield,False);
#
#    rootfile_run3_mc   = TFile.Open(filename_mc  , "READ");
#    list_run3_mc_ss     = rootfile_run3_mc.Get(ssname);
#    list_run3_mc_ss_cut = list_run3_mc_ss.FindObject(cutname);
#    list_run3_mc_ss_cut_fit           = list_run3_mc_ss_cut  .FindObject("cb_pol1");
#    list_run3_mc_ss_cut_fit_range     = list_run3_mc_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
#    list_run3_mc_ss_cut_fit_range_int = list_run3_mc_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
#    h1eff = list_run3_mc_ss_cut_fit_range_int.FindObject("h1eff");
#    h1eff.SetDirectory(0);
#    ROOT.SetOwnership(h1eff,False);
#
#    h1iy = h1yield.Clone("h1iy");
#    h1iy.Reset();
#    h1iy.Sumw2();
#    h1iy.Divide(h1yield, h1eff, 1., 1., "G");
#    for i in range(0, h1iy.GetNbinsX()):
#        n = h1iy.GetBinContent(i+1);
#        n_err = h1iy.GetBinError(i+1);
#        pt = h1iy.GetBinCenter(i+1);
#        h1iy.SetBinContent(i+1, n/pt);
#        h1iy.SetBinError(i+1, n_err/pt);
#
#    h1iy.Scale(1/TMath.TwoPi());
#    h1iy.Scale(1/1.8);
#    h1iy.Scale(59.4); #TVX trigger cross section in mb
#    h1iy.SetDirectory(0);
#    g1iy = apply_bin_shift_correction_x(h1iy);
#    ROOT.SetOwnership(g1iy, False);
#    ROOT.SetOwnership(h1iy, False);
#
#    #for inclusive photon
#    filename_photon_data = "photon_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_xsection.root";
#    filename_photon_mc   = "photon_mc_ptspectrum_pp_13.6TeV_LHC23d10_V0withAnyTrack_xsection.root";
#
#    rootfile_photon_data = TFile.Open(filename_photon_data, "READ");
#    list_photon_data = rootfile_photon_data.Get("PCM");
#    list_photon_data_cut = list_photon_data.FindObject("analysis_wo_mee");
#    h1yield_photon = list_photon_data_cut.FindObject("h1yield");
#    h1yield_photon.SetDirectory(0);
#    ROOT.SetOwnership(h1yield_photon, False);
#
#    rootfile_photon_mc   = TFile.Open(filename_photon_mc  , "READ");
#    list_photon_mc = rootfile_photon_mc.Get("PCM");
#    list_photon_mc_cut = list_photon_mc.FindObject("analysis_wo_mee");
#    h1eff_photon = list_photon_mc_cut.FindObject("h1eff");
#    h1eff_photon.SetDirectory(0);
#    ROOT.SetOwnership(h1eff_photon, False);
#    h1purity_photon = list_photon_mc_cut.FindObject("h1purity");
#
#    h1iy_photon = h1yield_photon.Clone("h1iy_photon");
#    h1iy_photon.Reset();
#    h1iy_photon.Sumw2();
#    h1iy_photon.Divide(h1yield_photon, h1eff_photon, 1., 1., "G");
#    for i in range(0, h1iy.GetNbinsX()):
#        n     = h1iy_photon.GetBinContent(i+1);
#        n_err = h1iy_photon.GetBinError(i+1);
#        pt    = h1iy_photon.GetBinCenter(i+1);
#        h1iy_photon.SetBinContent(i+1, n/pt);
#        h1iy_photon.SetBinError(i+1, n_err/pt);
#    h1iy_photon.Multiply(h1purity_photon);
#    h1iy_photon.Scale(0.97); #feed down
#    h1iy_photon.Scale(1/TMath.TwoPi());
#    h1iy_photon.Scale(1/1.8);
#    h1iy_photon.Scale(59.4); #TVX trigger cross section in mb
#    h1iy_photon.SetDirectory(0);
#    ROOT.SetOwnership(h1iy_photon,False);
#    make_common_style(h1iy_photon, 20, 1.0, kBlue+1, 1, 0);
#
#    outfile = TFile("tmp_output_pi0_gamma.root","RECREATE");
#    outfile.WriteTObject(h1iy);
#    outfile.WriteTObject(h1iy_photon);
#    outfile.Close();

    #for eta
    rootfile_eta_run3_data = TFile.Open(filename_data, "READ");
    list_eta_run3_data_ss     = rootfile_eta_run3_data.Get(ssname);
    list_eta_run3_data_ss_cut = list_eta_run3_data_ss.FindObject(cutname);
    list_eta_run3_data_ss_cut_fit           = list_eta_run3_data_ss_cut  .FindObject("cb_pol1");
    list_eta_run3_data_ss_cut_fit_range     = list_eta_run3_data_ss_cut_fit  .FindObject("fit_0.40_0.70_GeVc2");
    list_eta_run3_data_ss_cut_fit_range_int = list_eta_run3_data_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1yield_eta = list_eta_run3_data_ss_cut_fit_range_int.FindObject("h1yield");
    h1yield_eta.SetDirectory(0);
    ROOT.SetOwnership(h1yield_eta,False);

    rootfile_eta_run3_mc = TFile.Open(filename_mc, "READ");
    list_eta_run3_mc_ss     = rootfile_eta_run3_mc.Get(ssname);
    list_eta_run3_mc_ss_cut = list_eta_run3_mc_ss.FindObject(cutname);
    list_eta_run3_mc_ss_cut_fit           = list_eta_run3_mc_ss_cut  .FindObject("cb_pol1");
    list_eta_run3_mc_ss_cut_fit_range     = list_eta_run3_mc_ss_cut_fit  .FindObject("fit_0.40_0.70_GeVc2");
    list_eta_run3_mc_ss_cut_fit_range_int = list_eta_run3_mc_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1eff_eta = list_eta_run3_mc_ss_cut_fit_range_int.FindObject("h1eff");
    h1eff_eta.SetDirectory(0);
    ROOT.SetOwnership(h1eff_eta,False);

    h1iy_eta = h1yield_eta.Clone("h1iy_eta");
    h1iy_eta.Reset();
    h1iy_eta.Sumw2();
    h1iy_eta.Divide(h1yield_eta, h1eff_eta, 1., 1., "G");
    for i in range(0, h1iy_eta.GetNbinsX()):
        n     = h1iy_eta.GetBinContent(i+1);
        n_err = h1iy_eta.GetBinError(i+1);
        pt    = h1iy_eta.GetBinCenter(i+1);
        h1iy_eta.SetBinContent(i+1, n/pt);
        h1iy_eta.SetBinError(i+1, n_err/pt);

    h1iy_eta.Scale(1/TMath.TwoPi());
    h1iy_eta.Scale(1/1.8);
    h1iy_eta.Scale(59.4); #TVX trigger cross section in mb
    h1iy_eta.SetDirectory(0);
    h1iy_eta.SetTitle("cross section");
    h1iy_eta.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    h1iy_eta.GetYaxis().SetTitle("#it{E} #frac{d^{3}#it{#sigma}}{d#it{p}^{3}} (mb GeV^{#minus2} #it{c}^{3})");
    ROOT.SetOwnership(h1iy_eta, False);

    outfile = TFile("output_eta_cross_section.root","RECREATE");
    outfile.WriteTObject(h1yield_eta);
    outfile.WriteTObject(h1eff_eta);
    outfile.WriteTObject(h1iy_eta);
    outfile.Close();

    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetMargin(0.14,0.02,0.11,0.02);
    c1.SetTicks(1,1);
    c1.SetLogx(1);
    c1.SetLogy(1);
    frame1 = c1.DrawFrame(0.1,1e-6, 30, 1e+3);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#it{E} #frac{d^{3}#it{#sigma}}{d#it{p}^{3}} (mb GeV^{#minus2} #it{c}^{3})");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.5);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    #make_common_style(g1_run2_stat, 24, 1.0, kGray+2, 1, 0);
    #make_common_style(g1_run2_syst, 24, 1.0, kGray+2, 1, 0);
    #g1_run2_stat.Draw("PZ");
    #g1_run2_syst.Draw("P2Z");


    make_common_style(g1_eta_run2_stat, 24, 1.0, kGray+2, 1, 0);
    make_common_style(g1_eta_run2_syst, 24, 1.0, kGray+2, 1, 0);
    g1_eta_run2_stat.Draw("PZ");
    g1_eta_run2_syst.Draw("P2Z");

    #make_common_style(h1iy, 20, 1.0, kRed+1, 1, 0);
    #make_common_style(g1iy, 20, 1.0, kRed+1, 1, 0);
    ##g1iy.Draw("PZ");
    #h1iy.Draw("E1same");


    make_common_style(h1iy_eta, 20, 1.0, kBlue+1, 1, 0);
    h1iy_eta.Draw("E1same");

    #h1iy_photon.Draw("E1same");
    #h1yield_eta.Draw("E1same");
    #h1eff_eta.Draw("E1same");

    #f1 = TF1("f1","[0]*pow(1 + x/[1],-[2])",0,10);#for d2N/dpt/dy
    #f1.SetLineColor(kRed+1);
    #f1.SetLineStyle(7);
    #f1.SetNpx(1000);
    #f1.SetParameters(12000,0.6,6);
    #f1.SetParLimits(0,0,1e+5);
    #f1.SetParLimits(1,0,10);
    #f1.SetParLimits(2,0,10);
    #g1iy.Fit(f1, "SME", "", 0.6, 8);
    #ROOT.SetOwnership(f1,False);

    #leg_pi0 = TLegend(0.55,0.8,0.75,0.92);
    #leg_pi0.SetBorderSize(0);
    #leg_pi0.SetFillColor(kWhite);
    #leg_pi0.SetFillStyle(0);
    #leg_pi0.SetTextSize(0.04);
    #leg_pi0.SetHeader("#pi^{0} , PCM-PCM");
    #leg_pi0.AddEntry(h1iy  ,"pp at #sqrt{#it{s}} = 13.6 TeV","LP");
    #leg_pi0.AddEntry(g1_run2_syst  ,"pp at #sqrt{#it{s}} = 13 TeV","FP");
    #leg_pi0.Draw("");
    #ROOT.SetOwnership(leg_pi0,False);

    leg_eta = TLegend(0.2,0.35,0.4,0.5);
    leg_eta.SetBorderSize(0);
    leg_eta.SetFillColor(kWhite);
    leg_eta.SetFillStyle(0);
    leg_eta.SetTextSize(0.04);
    leg_eta.SetHeader("#eta , PCMPCM");
    leg_eta.AddEntry(h1iy_eta  ,"pp at #sqrt{#it{s}} = 13.6 TeV","LP");
    leg_eta.AddEntry(g1_eta_run2_syst  ,"pp at #sqrt{#it{s}} = 13 TeV","FP");
    leg_eta.Draw("");
    ROOT.SetOwnership(leg_eta,False);

    #A = f1.GetParameter(0);
    #A_err = f1.GetParError(0);
    #T = f1.GetParameter(1);
    #T_err = f1.GetParError(1);
    #n = f1.GetParameter(2);
    #n_err = f1.GetParError(2);

    #txt = TPaveText(0.2,0.2,0.4,0.5,"NDC");
    #txt.SetFillColor(kWhite);
    #txt.SetFillStyle(0);
    #txt.SetBorderSize(0);
    #txt.SetTextAlign(12);#middle,left
    #txt.SetTextFont(42);#helvetica
    #txt.SetTextSize(0.04);
    #txt.AddText("Hagedorn fit");
    #txt.AddText("A #times (1 + #frac{p_{T}}{T})^{-n}");
    #txt.AddText("A = {0:3.2f} #pm {1:3.2f} (mb GeV^{{-2}})".format(A, A_err));
    #txt.AddText("T = {0:3.2f} #pm {1:3.2f} (GeV)".format(T, T_err));
    #txt.AddText("n = {0:3.2f} #pm {1:3.2f}".format(n, n_err));
    #txt.Draw();
    #ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_eta_cross_section_{2}_{3}.eps".format(date, period, ssname, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_eta_cross_section_{2}_{3}.pdf".format(date, period, ssname, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_eta_cross_section_{2}_{3}.png".format(date, period, ssname, cutname));

#_____________________________________________________________________
def draw_cross_section_pi0(filename_data, filename_mc, period, ssname, cutname):

    filename_run2 = "data_PCMResultsFullCorrection_PP.root";
    rootfile_run2 = TFile.Open(filename_run2, "READ");
    dir_pi0_run2 = rootfile_run2.Get("Pi013TeV");
    g1_run2_stat = dir_pi0_run2.Get("graphInvCrossSectionPi0");#in pb GeV-2
    g1_run2_syst = dir_pi0_run2.Get("InvCrossSectionPi0Sys");# in pb GeV-2
    ROOT.SetOwnership(g1_run2_stat,False);
    ROOT.SetOwnership(g1_run2_syst,False);
    np = g1_run2_stat.GetN();
    sf = 1e-12/1e-3;
    for ip in range(0, np):
        y = g1_run2_stat.GetPointY(ip);
        y_stat_high = g1_run2_stat.GetErrorYhigh(ip);
        y_stat_low  = g1_run2_stat.GetErrorYlow(ip);
        y_syst_high = g1_run2_syst.GetErrorYhigh(ip);
        y_syst_low  = g1_run2_syst.GetErrorYlow(ip);

        g1_run2_stat.SetPointY(ip, y * sf);
        g1_run2_stat.SetPointEYhigh(ip, y_stat_high * sf);
        g1_run2_stat.SetPointEYlow(ip, y_stat_low * sf);

        g1_run2_syst.SetPointY(ip, y * sf);
        g1_run2_syst.SetPointEYhigh(ip, y_syst_high * sf);
        g1_run2_syst.SetPointEYlow(ip, y_syst_low * sf);

    #run2 eta
    dir_eta_run2 = rootfile_run2.Get("Eta13TeV");
    g1_eta_run2_stat = dir_eta_run2.Get("graphInvCrossSectionEta");#in pb GeV-2
    g1_eta_run2_syst = dir_eta_run2.Get("InvCrossSectionEtaSys");# in pb GeV-2
    ROOT.SetOwnership(g1_eta_run2_stat,False);
    ROOT.SetOwnership(g1_eta_run2_syst,False);
    np = g1_eta_run2_stat.GetN();
    sf = 1e-12/1e-3 * 0.01;
    for ip in range(0, np):
        y = g1_eta_run2_stat.GetPointY(ip);
        y_stat_high = g1_eta_run2_stat.GetErrorYhigh(ip);
        y_stat_low  = g1_eta_run2_stat.GetErrorYlow(ip);
        y_syst_high = g1_eta_run2_syst.GetErrorYhigh(ip);
        y_syst_low  = g1_eta_run2_syst.GetErrorYlow(ip);

        g1_eta_run2_stat.SetPointY(ip, y * sf);
        g1_eta_run2_stat.SetPointEYhigh(ip, y_stat_high * sf);
        g1_eta_run2_stat.SetPointEYlow(ip, y_stat_low * sf);

        g1_eta_run2_syst.SetPointY(ip, y * sf);
        g1_eta_run2_syst.SetPointEYhigh(ip, y_syst_high * sf);
        g1_eta_run2_syst.SetPointEYlow(ip, y_syst_low * sf);

    rootfile_run3_data = TFile.Open(filename_data, "READ");
    list_run3_data_ss     = rootfile_run3_data.Get(ssname);
    list_run3_data_ss.Print();
    list_run3_data_ss_cut = list_run3_data_ss.FindObject(cutname);
    list_run3_data_ss_cut_fit           = list_run3_data_ss_cut  .FindObject("cb_pol1");
    list_run3_data_ss_cut_fit_range     = list_run3_data_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_run3_data_ss_cut_fit_range_int = list_run3_data_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1yield = list_run3_data_ss_cut_fit_range_int.FindObject("h1yield");
    h1yield.SetDirectory(0);
    ROOT.SetOwnership(h1yield,False);

    rootfile_run3_mc   = TFile.Open(filename_mc  , "READ");
    list_run3_mc_ss     = rootfile_run3_mc.Get(ssname);
    list_run3_mc_ss_cut = list_run3_mc_ss.FindObject(cutname);
    list_run3_mc_ss_cut_fit           = list_run3_mc_ss_cut  .FindObject("cb_pol1");
    list_run3_mc_ss_cut_fit_range     = list_run3_mc_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_run3_mc_ss_cut_fit_range_int = list_run3_mc_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1eff = list_run3_mc_ss_cut_fit_range_int.FindObject("h1eff");
    h1eff.SetDirectory(0);
    ROOT.SetOwnership(h1eff,False);

    h1iy = h1yield.Clone("h1iy");
    h1iy.Reset();
    h1iy.Sumw2();
    h1iy.Divide(h1yield, h1eff, 1., 1., "G");
    for i in range(0, h1iy.GetNbinsX()):
        n = h1iy.GetBinContent(i+1);
        n_err = h1iy.GetBinError(i+1);
        pt = h1iy.GetBinCenter(i+1);
        h1iy.SetBinContent(i+1, n/pt);
        h1iy.SetBinError(i+1, n_err/pt);

    h1iy.Scale(1/TMath.TwoPi());
    h1iy.Scale(1/1.8);
    h1iy.Scale(59.4); #TVX trigger cross section in mb
    h1iy.SetDirectory(0);
    g1iy = apply_bin_shift_correction_x(h1iy);
    ROOT.SetOwnership(g1iy, False);
    ROOT.SetOwnership(h1iy, False);
    h1iy.SetTitle("cross section");
    h1iy.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    h1iy.GetYaxis().SetTitle("#it{E} #frac{d^{3}#it{#sigma}}{d#it{p}^{3}} (mb GeV^{#minus2} #it{c}^{3})");


    #for inclusive photon
    filename_photon_data = "photon_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_xsection.root";
    filename_photon_mc   = "photon_mc_ptspectrum_pp_13.6TeV_LHC23d10_V0withAnyTrack_xsection.root";

    rootfile_photon_data = TFile.Open(filename_photon_data, "READ");
    list_photon_data = rootfile_photon_data.Get("PCM");
    list_photon_data_cut = list_photon_data.FindObject("analysis_wo_mee");
    h1yield_photon = list_photon_data_cut.FindObject("h1yield");
    h1yield_photon.SetDirectory(0);
    ROOT.SetOwnership(h1yield_photon, False);

    rootfile_photon_mc   = TFile.Open(filename_photon_mc  , "READ");
    list_photon_mc = rootfile_photon_mc.Get("PCM");
    list_photon_mc_cut = list_photon_mc.FindObject("analysis_wo_mee");
    h1eff_photon = list_photon_mc_cut.FindObject("h1eff");
    h1eff_photon.SetDirectory(0);
    ROOT.SetOwnership(h1eff_photon, False);
    h1purity_photon = list_photon_mc_cut.FindObject("h1purity");

    h1iy_photon = h1yield_photon.Clone("h1iy_photon");
    h1iy_photon.Reset();
    h1iy_photon.Sumw2();
    h1iy_photon.Divide(h1yield_photon, h1eff_photon, 1., 1., "G");
    for i in range(0, h1iy.GetNbinsX()):
        n     = h1iy_photon.GetBinContent(i+1);
        n_err = h1iy_photon.GetBinError(i+1);
        pt    = h1iy_photon.GetBinCenter(i+1);
        h1iy_photon.SetBinContent(i+1, n/pt);
        h1iy_photon.SetBinError(i+1, n_err/pt);
    h1iy_photon.Multiply(h1purity_photon);
    h1iy_photon.Scale(0.97); #feed down
    h1iy_photon.Scale(1/TMath.TwoPi());
    h1iy_photon.Scale(1/1.8);
    h1iy_photon.Scale(59.4); #TVX trigger cross section in mb
    h1iy_photon.SetDirectory(0);
    ROOT.SetOwnership(h1iy_photon,False);
    make_common_style(h1iy_photon, 20, 1.0, kBlue+1, 1, 0);

    outfile = TFile("output_pi0_cross_section.root","RECREATE");
    outfile.WriteTObject(h1yield);
    outfile.WriteTObject(h1eff);
    outfile.WriteTObject(h1iy);
    outfile.Close();


    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetMargin(0.14,0.02,0.11,0.02);
    c1.SetTicks(1,1);
    c1.SetLogx(1);
    c1.SetLogy(1);
    frame1 = c1.DrawFrame(0.1,1e-6, 30, 1e+3);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#it{E} #frac{d^{3}#it{#sigma}}{d#it{p}^{3}} (mb GeV^{#minus2} #it{c}^{3})");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.5);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    make_common_style(g1_run2_stat, 24, 1.0, kGray+2, 1, 0);
    make_common_style(g1_run2_syst, 24, 1.0, kGray+2, 1, 0);
    g1_run2_stat.Draw("PZ");
    g1_run2_syst.Draw("P2Z");


    make_common_style(h1iy, 20, 1.0, kRed+1, 1, 0);
    make_common_style(g1iy, 20, 1.0, kRed+1, 1, 0);
    #g1iy.Draw("PZ");
    h1iy.Draw("E1same");


    #make_common_style(h1iy_eta, 20, 1.0, kBlue+1, 1, 0);
    #h1iy_eta.Draw("E1same");

    #h1iy_photon.Draw("E1same");
    #h1yield_eta.Draw("E1same");
    #h1eff_eta.Draw("E1same");

    #f1 = TF1("f1","[0]*pow(1 + x/[1],-[2])",0,10);#for d2N/dpt/dy
    #f1.SetLineColor(kRed+1);
    #f1.SetLineStyle(7);
    #f1.SetNpx(1000);
    #f1.SetParameters(12000,0.6,6);
    #f1.SetParLimits(0,0,1e+5);
    #f1.SetParLimits(1,0,10);
    #f1.SetParLimits(2,0,10);
    #g1iy.Fit(f1, "SME", "", 0.6, 8);
    #ROOT.SetOwnership(f1,False);

    leg_pi0 = TLegend(0.55,0.8,0.75,0.92);
    leg_pi0.SetBorderSize(0);
    leg_pi0.SetFillColor(kWhite);
    leg_pi0.SetFillStyle(0);
    leg_pi0.SetTextSize(0.04);
    leg_pi0.SetHeader("#pi^{0} , PCM-PCM");
    leg_pi0.AddEntry(h1iy  ,"pp at #sqrt{#it{s}} = 13.6 TeV","LP");
    leg_pi0.AddEntry(g1_run2_syst  ,"pp at #sqrt{#it{s}} = 13 TeV","FP");
    leg_pi0.Draw("");
    ROOT.SetOwnership(leg_pi0,False);


    #A = f1.GetParameter(0);
    #A_err = f1.GetParError(0);
    #T = f1.GetParameter(1);
    #T_err = f1.GetParError(1);
    #n = f1.GetParameter(2);
    #n_err = f1.GetParError(2);

    #txt = TPaveText(0.2,0.2,0.4,0.5,"NDC");
    #txt.SetFillColor(kWhite);
    #txt.SetFillStyle(0);
    #txt.SetBorderSize(0);
    #txt.SetTextAlign(12);#middle,left
    #txt.SetTextFont(42);#helvetica
    #txt.SetTextSize(0.04);
    #txt.AddText("Hagedorn fit");
    #txt.AddText("A #times (1 + #frac{p_{T}}{T})^{-n}");
    #txt.AddText("A = {0:3.2f} #pm {1:3.2f} (mb GeV^{{-2}})".format(A, A_err));
    #txt.AddText("T = {0:3.2f} #pm {1:3.2f} (GeV)".format(T, T_err));
    #txt.AddText("n = {0:3.2f} #pm {1:3.2f}".format(n, n_err));
    #txt.Draw();
    #ROOT.SetOwnership(txt,False);



    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pi0_cross_section_{2}_{3}.eps".format(date, period, ssname, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pi0_cross_section_{2}_{3}.pdf".format(date, period, ssname, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_pi0_cross_section_{2}_{3}.png".format(date, period, ssname, cutname));

#_____________________________________________________________________
#_____________________________________________________________________
def draw_cross_section_pi0_pbpb(filename_data, filename_mc, ssname, cutname):

    rootfile_run2 = TFile.Open("charged_pion_PbPb5.02TeV_0090.root", "READ")
    h1pi_run2 = rootfile_run2.Get("h1all");
    h1pi_run2.SetDirectory(0);
    ROOT.SetOwnership(h1pi_run2, False);
    make_common_style(h1pi_run2, 24, 1.2, kBlack, 2, 0);
    h1pi_run2.GetXaxis().SetRangeUser(0.6, 10);

    rootfile_run3_data = TFile.Open(filename_data, "READ");
    list_run3_data_ss     = rootfile_run3_data.Get(ssname);
    list_run3_data_ss_cut = list_run3_data_ss.FindObject(cutname);
    list_run3_data_ss_cut_fit           = list_run3_data_ss_cut  .FindObject("cb_pol1");
    list_run3_data_ss_cut_fit_range     = list_run3_data_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_run3_data_ss_cut_fit_range_int = list_run3_data_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1yield = list_run3_data_ss_cut_fit_range_int.FindObject("h1yield");
    h1yield.SetDirectory(0);
    ROOT.SetOwnership(h1yield,False);

    rootfile_run3_mc   = TFile.Open(filename_mc  , "READ");
    list_run3_mc_ss     = rootfile_run3_mc.Get(ssname);
    list_run3_mc_ss_cut = list_run3_mc_ss.FindObject(cutname);
    list_run3_mc_ss_cut_fit           = list_run3_mc_ss_cut  .FindObject("cb_pol1");
    list_run3_mc_ss_cut_fit_range     = list_run3_mc_ss_cut_fit  .FindObject("fit_0.04_0.24_GeVc2");
    list_run3_mc_ss_cut_fit_range_int = list_run3_mc_ss_cut_fit_range  .FindObject("integral_-3.0_3.0_sigma");
    h1eff = list_run3_mc_ss_cut_fit_range_int.FindObject("h1eff");
    h1eff.SetDirectory(0);
    ROOT.SetOwnership(h1eff,False);

    h1iy = h1yield.Clone("h1iy");
    h1iy.Reset();
    h1iy.Sumw2();
    h1iy.Divide(h1yield, h1eff, 1., 1., "G");
    for i in range(0, h1iy.GetNbinsX()):
        n = h1iy.GetBinContent(i+1);
        n_err = h1iy.GetBinError(i+1);
        pt = h1iy.GetBinCenter(i+1);
        h1iy.SetBinContent(i+1, n/pt);
        h1iy.SetBinError(i+1, n_err/pt);

    h1iy.Scale(1/TMath.TwoPi());
    h1iy.Scale(1/1.8); # dy = 1.8
    h1iy.SetDirectory(0);
    g1iy = apply_bin_shift_correction_x(h1iy);
    ROOT.SetOwnership(g1iy, False);
    ROOT.SetOwnership(h1iy, False);

    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetMargin(0.16,0.03,0.11,0.02);
    c1.SetTicks(1,1);
    c1.SetLogx(1);
    c1.SetLogy(1);
    frame1 = c1.DrawFrame(0.6, 1e-5, 10, 1e+3);
    frame1.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#frac{1}{2#pi N_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T} d#it{p}_{T} d#it{y}} (GeV/#it{c})^{#minus2}");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.7);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1pi_run2.Draw("E0same");

    make_common_style(h1iy, 20, 1.2, kRed+1, 2, 0);
    make_common_style(g1iy, 20, 1.2, kRed+1, 2, 0);
    #g1iy.Draw("PZ");
    h1iy.Draw("E0same");

    leg = TLegend(0.18, 0.85,0.35,0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.AddEntry(h1iy, "#pi^{0} : Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.36 TeV, 0#minus100% (This work)", "LP");
    leg.AddEntry(h1pi_run2, "#frac{#pi^{+}+#pi^{#minus}}{2} : Pb#minusPb at #sqrt{#it{s}_{NN}} = 5.02 TeV, 0#minus90%", "LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.68,0.82,0.88,0.85,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("PRC 101, 044907");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    period = "LHC23zzf";
    c1.SaveAs("{0}_PbPb_5.36TeV_{1}_pass2_{2}_{3}.eps".format(date, period, ssname, cutname));
    c1.SaveAs("{0}_PbPb_5.36TeV_{1}_pass2_{2}_{3}.pdf".format(date, period, ssname, cutname));
    c1.SaveAs("{0}_PbPb_5.36TeV_{1}_pass2_{2}_{3}.png".format(date, period, ssname, cutname));

#_____________________________________________________________________
def compare_rxy(filename_pp, filename_pbpb, cutname):
    rootfile_pp   = TFile.Open(filename_pp  , "READ");
    rootfile_pbpb = TFile.Open(filename_pbpb, "READ");

    rootfile_pp  .ls();
    rootfile_pbpb.ls();

    rootdire_pp   = rootfile_pp  .Get("pcm-qc");
    rootdire_pbpb = rootfile_pbpb.Get("pcm-qc");
    rootdire_pp  .ls();
    rootdire_pbpb.ls();

    rootlist_v0_pp   = rootdire_pp  .Get("V0");
    rootlist_v0_pbpb = rootdire_pbpb.Get("V0");

    rootlist_v0_cut_pp   = rootlist_v0_pp  .FindObject(cutname);
    rootlist_v0_cut_pbpb = rootlist_v0_pbpb.FindObject(cutname);
    rootlist_v0_cut_pp  .ls();
    rootlist_v0_cut_pbpb.ls();

    h2_pp   = rootlist_v0_cut_pp  .FindObject("hMassGamma");
    h2_pbpb = rootlist_v0_cut_pbpb.FindObject("hMassGamma");
    h1r_pp   = h2_pp  .ProjectionX("h1r_pp");
    h1r_pbpb = h2_pbpb.ProjectionX("h1r_pbpb");

    h2_pp  .SetDirectory(0);
    h2_pbpb.SetDirectory(0);
    ROOT.SetOwnership(h2_pp  ,False);
    ROOT.SetOwnership(h2_pbpb,False);

    h1r_pp  .SetDirectory(0);
    h1r_pbpb.SetDirectory(0);
    ROOT.SetOwnership(h1r_pp  ,False);
    ROOT.SetOwnership(h1r_pbpb,False);

    rootlist_ev_pp   = rootdire_pp  .Get("Event");
    rootlist_ev_pbpb = rootdire_pbpb.Get("Event");
    rootlist_ev_pp  .ls();
    rootlist_ev_pbpb.ls();

    h1mult_pp   = rootlist_ev_pp  .FindObject("hMultNTracksPV");
    h1mult_pbpb = rootlist_ev_pbpb.FindObject("hMultNTracksPV");
    nev_pp   = h1mult_pp  .GetEntries();
    nev_pbpb = h1mult_pbpb.GetEntries();
    mean_mult_pp   = h1mult_pp  .GetMean();
    mean_mult_pbpb = h1mult_pbpb.GetMean();

    print(nev_pp, nev_pbpb, mean_mult_pp, mean_mult_pbpb);

    h1r_pp  .Scale(1/nev_pp);
    h1r_pbpb.Scale(1/nev_pbpb);
    h1r_pp  .Scale(1/mean_mult_pp);
    h1r_pbpb.Scale(1/mean_mult_pbpb);
    h1r_pp  .Scale(1, "width");
    h1r_pbpb.Scale(1, "width");

    c1 = TCanvas("c0","c0",0,0,800,800);
    c1.SetMargin(0.16,0.03,0.09,0.05);
    c1.SetTicks(1,1);
    #c1.SetLogy(1);
    frame1 = c1.DrawFrame(0, 0, +90, 7e-4);
    frame1.GetXaxis().SetTitle("conversion radius r_{xy} (cm)");
    frame1.GetYaxis().SetTitle("#frac{1}{N_{ev}} #frac{1}{<N_{ch}^{PV}>} #frac{dN_{#gamma}}{dr_{xy}} (cm)^{-1}");
    frame1.GetXaxis().SetTitleOffset(1.0);
    frame1.GetYaxis().SetTitleOffset(1.7);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetMaxDigits(3);
    ROOT.SetOwnership(frame1,False);

    make_common_style(h1r_pp  , 20, 1.0, kRed+1 , 2, 0);
    make_common_style(h1r_pbpb, 21, 1.0, kBlue+1, 2, 0);

    h1r_pp  .Draw("E0same,h");
    h1r_pbpb.Draw("E0same,h");

    leg = TLegend(0.2,0.8,0.4,0.92);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.SetHeader("pp at 13.6 TeV");
    #leg.SetHeader("centFT0C < 80%");
    #leg.AddEntry(h1r_pp,"LHC22o pass4", "LP");
    #leg.AddEntry(h1r_pp,"LHC22o pass4", "LP");
    #leg.AddEntry(h1r_pbpb,"LHC23zzf pass1, cent FT0C < 80%", "LP");
    #leg.AddEntry(h1r_pp,"LHC23zzg pass1 544028 30 kHz", "LP");
    #leg.AddEntry(h1r_pbpb,"LHC23zzf pass1 544013 6 kHz", "LP");
    #leg.AddEntry(h1r_pp,"LHC22o pass4 500 kHz (SVertexer)", "LP");
    leg.AddEntry(h1r_pp,"LHC22f pass4 10 kHz (SVertexer)", "LP");
    leg.AddEntry(h1r_pbpb,"LHC22f pass4 10 kHz (Daiki's pairing)", "LP");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_rxy_comparison_rate_pp_LHC22f_{1}.eps".format(date, cutname));
    c1.SaveAs("{0}_rxy_comparison_rate_pp_LHC22f_{1}.pdf".format(date, cutname));
    c1.SaveAs("{0}_rxy_comparison_rate_pp_LHC22f_{1}.png".format(date, cutname));

#    rootfile_pp  .Close();
#    rootfile_pbpb.Close();

#_____________________________________________________________________
if __name__ == "__main__":
    #draw_raw_yield_pi0( "PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee");
    #draw_raw_yield_pi0( "PCMPCM", "LHC22m", "qc_qc");
    #draw_raw_yield_pi0( "PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee");
    #draw_peak_mean_pi0( "PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee");
    #draw_peak_sigma_pi0("PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTrack_KF_PV_M0.root", "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 0, "KF_AnyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTrack.root"         , "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 0, "AnyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTrack_KF_PV_M0.root", "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 1, "KF_AnyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTrack.root"         , "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 1, "AnyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTrack_KF_PV_M0.root", "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 2, "KF_AnyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTrack.root"         , "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 2, "AnyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTrack_KF_PV_M0.root", "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 3, "KF_AnyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTrack.root"         , "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 3, "AnyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTrack_KF_PV_M0.root", "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 4, "KF_AnyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTrack.root"         , "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 4, "AnyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withAnyTracks.root", "PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee", 1, "AnyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withTPConlyTracks.root", "PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee", 1, "TPConlyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withTPConlyTracksBar.root", "PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee", 1, "TPConlyTracksBar");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withLKB.root", "PCMPCM", "LHC22m", "analysis_wo_mee_analysis_wo_mee", 1, "LKB");

    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTrack_KF_PV_M0_PID3.root"        , "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 0, "KF_AnyTracks_PID3");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTrack_KF_PV_M0_PID3_ITSonly.root", "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 0, "KF_AnyTracks_PID3_ITSonly");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTrack_KF_PV_M0_PID3.root"        , "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 1, "KF_AnyTracks_PID3");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTrack_KF_PV_M0_PID3_ITSonly.root", "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 1, "KF_AnyTracks_PID3_ITSonly");

    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withAnyTrack_nsw1_20230721.root", "PCMPCM", "LHC22m", "analysis_analysis", 11, "nsw1");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withAnyTrack_nsw5_20230721.root", "PCMPCM", "LHC22m", "analysis_analysis", 11, "nsw5");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22o_V0withAnyTrack_nsw1_20230721.root", "PCMPCM", "LHC22o", "analysis_analysis", 11, "nsw1");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22o_V0withAnyTrack_nsw5_20230721.root", "PCMPCM", "LHC22o", "analysis_analysis", 11, "nsw5");
#    draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withAnyTrack_nsw1_20230721.root", "PCMPCM", "LHC22m", "analysis_analysis", 6, "nsw1");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withAnyTrack_nsw5_20230721.root", "PCMPCM", "LHC22m", "analysis_analysis", 6, "nsw5");
#    draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22o_V0withAnyTrack_nsw1_20230721.root", "PCMPCM", "LHC22o", "analysis_analysis", 6, "nsw1");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22o_V0withAnyTrack_nsw5_20230721.root", "PCMPCM", "LHC22o", "analysis_analysis", 6, "nsw5");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22o_V0withAnyTrack_20230828_nsw1.root", "PCMPCM", "LHC22o", "analysis_analysis", 4, "nsw1");

    #draw_xy("AnalysisResults_HL_78451.root", "PCMPCM", "LHC22m", "analysis_wo_mee", "AnyTracks");
    #draw_xy("AnalysisResults_HL_78449.root", "PCMPCM", "LHC22m", "analysis_wo_mee", "TPConlyTracks");
    #draw_xy("AnalysisResults_HL_78448.root", "PCMPCM", "LHC22m", "analysis_wo_mee", "TPConlyTracksBar");
    #draw_xy("AnalysisResults_HL_78450.root", "PCMPCM", "LHC22m", "analysis_wo_mee", "LKB");

    #draw_xy("AnalysisResults_HL_88933.root", "PCMPCM", "LHC22f", "analysis_wo_mee", False, "AnyTracks");
    #draw_xy("AnalysisResults_HL_94252.root", "PCMPCM", "LHC22f", "analysis_wo_mee", False, "KF_AnyTracks");
    #draw_xy("AnalysisResults_HL_94252.root", "PCMPCM", "LHC22f", "analysis_wo_mee", True, "KF_AnyTracks");

    #draw_xy("AnalysisResults_HL_85687.root", "PCMPCM", "LHC22q", "analysis_wo_mee", False, "AnyTracks");
    #draw_xy("AnalysisResults_HL_85687.root", "PCMPCM", "LHC22q", "analysis_wo_mee", True , "AnyTracks");

    #draw_meeg_pi0("tagging_pi0_data_ptspectrum_pp_13.6TeV_LHC23zo_V0withPCB_20231108.root", "PCMDalitz", "LHC23zo", "qc_mee_0_120_tpconly_lowB", 0, "PCB");
    #draw_meeg_pi0("tagging_pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230806_nsw1.root", "PCMPCMibw", "LHC22q", "analysis1690_wwire_ib", 3, "AnyTrack");
    #draw_meeg_pi0("tagging_pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230806_nsw1.root", "PCMPHOS", "LHC22q", "analysis1690_test03", 1, "AnyTrack");
    #draw_meeg_pi0("tagging_pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230725_nsw1.root", "PCMPHOS", "LHC22q", "analysis1690_test03", 0, "AnyTrack");
    #draw_meeg_pi0("tagging_pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230725_nsw1.root", "PCMPHOS", "LHC22q", "analysis1690_test03", 1, "AnyTrack");
    #draw_meeg_pi0("tagging_pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230725_nsw1.root", "PCMPCMibw", "LHC22q", "analysis1690_analysis116", 0, "AnyTrack");
    #draw_meeg_pi0("tagging_pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230725_nsw1.root", "PCMPCMibw", "LHC22q", "analysis1690_analysis116", 1, "AnyTrack");
    #draw_meeg_pi0("tagging_pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230725_nsw1.root", "PCMPCMibw", "LHC22q", "analysis1690_wwire_ib", 0, "AnyTrack");
    #draw_meeg_pi0("tagging_pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230725_nsw1.root", "PCMPCMibw", "LHC22q", "analysis1690_wwire_ib", 1, "AnyTrack");
    #draw_meeg_pi0("tagging_pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack.root", "PCMPHOS", "LHC22q", "analysis_wo_mee_test03", 0, "AnyTracks");
    #draw_meeg_pi0("tagging_pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack.root", "PCMEMC", "LHC22q", "analysis_wo_mee_standard", 0, "AnyTracks");
    #draw_tagging(
    #    "photon_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230725_nsw1.root", "PCM", "analysis1690",
    #    "tagging_pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230725_nsw1.root",  "PCMPCMibw", "analysis1690_analysis116", "PCMPHOS", "analysis1690_test03",
    #    "LHC22q","");
    #draw_tagging(
    #    "photon_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230725_nsw1.root", "PCM", "qc",
    #    "tagging_pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230725_nsw1.root",  "PCMPCMibw", "analysis1690_analysis116", "PCMPHOS", "analysis1690_test03",
    #    "LHC22q","");

    #PCMPHOS
    ssname = "PCMPHOS";
    period = "LHC22m";
    #draw_raw_yield_pi0( ssname, period, "analysis_wo_mee_test03");
    #draw_peak_mean_pi0( ssname, period, "analysis_wo_mee_test03");
    #draw_peak_sigma_pi0(ssname, period, "analysis_wo_mee_test03");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withAnyTracks.root"       , ssname, period, "analysis_wo_mee_test03", 3, "AnyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withTPConlyTracks.root"   , ssname, period, "analysis_wo_mee_test03", 3, "TPConlyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withTPConlyTracksBar.root", ssname, period, "analysis_wo_mee_test03", 3, "TPConlyTracksBar");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withLKB.root"             , ssname, period, "analysis_wo_mee_test03", 3, "LKB");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22m_V0withAnyTracks.root"       , "PCMEMC", period, "analysis_wo_mee_standard", 5, "AnyTracks");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack.root"       , "PCMPCM", "LHC22q", "analysis_wo_mee_analysis_wo_mee", 6, "");
    #draw_mgg_eta("eta_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack.root"       , "PCMPCM", "LHC22q", "analysis_wo_mee_analysis_wo_mee", 0, "");
    #draw_mgg_eta("eta_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack.root"       , "PCMPCM", "LHC22q", "analysis_wo_mee_analysis_wo_mee", 1, "");
    #draw_mgg_eta("eta_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack.root"       , "PCMPCM", "LHC22q", "analysis_wo_mee_analysis_wo_mee", 2, "");
#    draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack.root"       , "EMCEMC", "LHC22q", "standard_standard", 21, ""); #For performance plots QM23
    #draw_mgg_eta("eta_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack.root"       , "EMCEMC", "LHC22q", "standard_standard", 7, "");
#    draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack.root"       , "PHOSPHOS", "LHC22q", "test03_test03", 18, ""); #For performance plots QM23
    #draw_mgg_eta("eta_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack.root"       , "PHOSPHOS", "LHC22q", "test03_test03", 5, "");

    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230711.root"       , "EMCEMC", "LHC22q", "standard_standard", 20, "");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230711.root"       , "PHOSPHOS", "LHC22q", "test03_test03", 20, "");
#    draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230806_nsw1.root"       , "PCMPCM", "LHC22q", "analysis_analysis", 6, ""); #For performance plots QM23
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230828_nsw1.root"       , "PCMPCM", "LHC22q", "analysis_analysis", 0, ""); #For FSP meeting 2023
    #draw_mgg_eta("eta_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230828_nsw1.root"       , "PCMPCM", "LHC22q", "analysis_analysis", 0, ""); #For FSP meeting 2023
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22o_V0withPCB_20231108.root"       , "PCMPCM", "LHC22o", "qc_qc", 0, "");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC23zo_V0withPCB_20231108.root"       , "PCMPCM", "LHC22zo", "qc_qc", 0, "");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22o_V0withPCB_20231108.root"       , "PCMPCM", "LHC22o" , "qc_qc", 7, "");
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC23zo_V0withPCB_20231108.root"      , "PCMPCM", "LHC22zo", "qc_qc", 7, "");
#    draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22o_V0withPCB_20231108_1bigbin.root"       , "PCMPCM", "LHC22o" , "qc_qc", 0, "1bigbin");
#    draw_mgg_eta("eta_data_ptspectrum_pp_13.6TeV_LHC22o_V0withPCB_20231108_1bigbin.root"       , "PCMPCM", "LHC22o" , "qc_qc", 0, "1bigbin");
#    draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC23zo_V0withPCB_20231108_1bigbin.root"       , "PCMPCM", "LHC23zo" , "qc_qc", 0, "1bigbin");
#    draw_mgg_eta("eta_data_ptspectrum_pp_13.6TeV_LHC23zo_V0withPCB_20231108_1bigbin.root"       , "PCMPCM", "LHC23zo" , "qc_qc", 0, "1bigbin");
#    draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC23zo_V0withPCB_20231205.root"       , "PCMPCM", "LHC23zo" , "qc_qc", 1, "");
#    draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC23zs_V0withPCB_20240213.root"       , "PCMPCM", "LHC23zs" , "qc_qc", 0, ""); # for Felix DPG 2024 march

    #draw_hbt_1d_qinv("photon_hbt_1d_data_pp_13.6TeV_LHC22q_V0withAnyTrack.root"       , "PCMPHOS", "LHC22q", "analysis_wo_mee_test03", 0, "");
    #draw_hbt_1d_qinv("photon_hbt_1d_data_pp_13.6TeV_LHC22q_V0withAnyTrack.root"       , "PCMPHOS", "LHC22q", "analysis_wo_mee_test03", 1, "");
    #draw_hbt_1d_qinv("photon_hbt_1d_data_pp_13.6TeV_LHC22q_V0withAnyTrack.root"       , "PCMPHOS", "LHC22q", "analysis_wo_mee_test03", 2, "");
    #draw_hbt_1d_qinv("photon_hbt_1d_data_pp_13.6TeV_LHC22q_V0withAnyTrack.root"       , "PCMPHOS", "LHC22q", "analysis_wo_mee_test03", 3, "");

    #qmin = -0.04;
    #qmax = +0.04;
    #draw_hbt_3d_qosl("photon_hbt_3d_data_pp_13.6TeV_LHC22q_V0withAnyTrack.root" , "PCMPHOS", "LHC22q", "analysis_wo_mee_test03"         , 0, qmin, qmax, "");
    #draw_hbt_3d_qosl("photon_hbt_3d_data_pp_13.6TeV_LHC22q_V0withAnyTrack.root" , "PCMPCM" , "LHC22q", "analysis_wo_mee_analysis_wo_mee", 0, qmin, qmax, "");

    #qmin = -0.02;
    #qmax = +0.02;
    #draw_hbt_3d_qosl("photon_hbt_3d_data_pp_13.6TeV_LHC22q_V0withAnyTrack.root"       , "PCMPHOS", "LHC22q", "analysis_wo_mee_test03", 0, qmin, qmax, "");
    #draw_hbt_3d_qosl("photon_hbt_3d_data_pp_13.6TeV_LHC22q_V0withAnyTrack.root"       , "PCMPCM", "LHC22q", "analysis_wo_mee_analysis_wo_mee", 0, qmin, qmax, "");


    #draw_raw_yield_pi0( "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee");
    #draw_peak_mean_pi0( "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee");
    #draw_peak_sigma_pi0("PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee");
    #draw_mgg_pi0("output_data_ptspectrum_pp_13.6TeV_LHC22f_V0withAnyTracks.root", "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 3, "AnyTrack");
    #draw_mgg_pi0("output_data_ptspectrum_pp_13.6TeV_LHC22f_V0withTPConlyTrack.root", "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 3, "TPConlyTrack");
    #draw_mgg_pi0("output_data_ptspectrum_pp_13.6TeV_LHC22f_V0withLKB.root", "PCMPCM", "LHC22f", "analysis_wo_mee_analysis_wo_mee", 1, "LKB");


    #draw_xy_ITSTPC("AnalysisResults_HL_137977.root", "PCMPCM", "LHC22o"  , "qc", True, "PCB");
    #draw_xy_ITSTPC("AnalysisResults_HL_137977.root", "PCMPCM", "LHC22o"  , "qc_ITSTPC", True, "PCB");
    #draw_xy_ITSTPC("AnalysisResults_HL_137977.root", "PCMPCM", "LHC22o"  , "qc_ITSonly", True, "PCB");
    #draw_xy_ITSTPC("AnalysisResults_HL_133366.root", "PCMPCM", "LHC22o"  , "qc", True, "PCB");
    #draw_xy_ITSTPC("AnalysisResults_HL_133022.root", "PCMPCM", "LHC23zo" , "qc", True, "PCB");
    #draw_xy_ITSTPC("AnalysisResults_HL_132818.root", "PCMPCM", "LHC23zzh", "qc_ITSTPC", True, "PCB_cpass1");
    #draw_xy_ITSTPC("AnalysisResults_HL_134949.root", "PCMPCM", "LHC23zzh", "qc_ITSTPC", True, "PCB_cpass2");
    #draw_xy_ITSTPC("AnalysisResults_HL_134949.root", "PCMPCM", "LHC23zzh", "qc", True, "PCB_cpass2");
    #draw_xy_ITSTPC("AnalysisResults_alien20_new.root", "PCMPCM", "LHC23zzk", "qc", True, "PCB_apass1_relval");
    #draw_xy_ITSTPC("AnalysisResults_HL_148948.root", "PCMPCM", "LHC23zzfghik", "qc", False, "PCB_apass1");
    #draw_xy_ITSTPC("AnalysisResults_HL_166514.root", "PCMPCM", "LHC23zs", "qc", False, "apass2"); #for Felix, 2024 march
    #draw_xy_ITSTPC("AnalysisResults_HL_168029.root", "PCMPCM", "LHC23zs", "qc", False, "apass3"); #for Felix, 2024 march
    #draw_xy_ITSTPC("AnalysisResults_HL_165787.root", "PCMPCM", "LHC23zzh", "qc", False, "apass2"); #for Felix, 2024 march
    #draw_rz_ITSTPC("AnalysisResults_HL_166514.root", "PCMPCM", "LHC23zs", "qc", False, "apass2"); #for Felix, 2024 march
    #draw_rz_ITSTPC("AnalysisResults_HL_166514.root", "PCMPCM", "LHC23zs", "nocut", False, "apass2"); #for Felix, 2024 march

    #draw_xy("AnalysisResults_HL_111543.root", "PCMPCM", "LHC22f", "analysis180", True, "AnyTrack"); #for performance plot QM23
    #draw_xy("AnalysisResults_HL_111543.root", "PCMPCM", "LHC22f", "qc", True, "AnyTrack"); #for performance plot QM23
    #draw_xy("AnalysisResults_HL_111543.root", "PCMPCM", "LHC22f", "qc_ITSTPC", True, "AnyTrack"); #for performance plot QM23
    #draw_xy("AnalysisResults_HL_111543.root", "PCMPCM", "LHC22f", "qc_TPConly", True, "AnyTrack"); #for performance plot QM23
    #draw_xy("AnalysisResults_HL_111543.root", "PCMPCM", "LHC22f", "qc_ITSonly", True, "AnyTrack"); #for performance plot QM23
    #draw_xy("AnalysisResults_HL_104885.root", "PCMPCM", "LHC22f", "analysis", True, "AnyTrack");
    #draw_xy("AnalysisResults_HL_75289.root", "PCMPCM", "LHC22q", "analysis_wo_mee", "AnyTrack");
    #draw_xy("AnalysisResults_HL_75290.root", "PCMPCM", "LHC22q", "analysis_wo_mee", "LKB");
    #draw_xy("AnalysisResults_HL_75291.root", "PCMPCM", "LHC22q", "analysis_wo_mee", "TPConlyTrack");
    #draw_xy_mc("AnalysisResults_HL_84266.root", "PCMPCM", "LHC23d1d", "analysis_wo_mee", "TPConlyTrack");
    #draw_xy_mc("AnalysisResults_LHC23d1f_520472.root", "PCMPCM", "LHC23d1f", "analysis_wo_mee", "AnyTrack");
    #draw_xy_mc_gen("AnalysisResults_HL_125889.root", "PCMPCM", "LHC23d1f", "");
    #draw_rz_mc_gen("AnalysisResults_HL_125889.root", "PCMPCM", "LHC23d1f", "");
    #check_material_rxy("AnalysisResults_HL_85104.root", "AnalysisResults_HL_84266.root", "analysis_wo_mee", "LHC23d1d", "TPConlyTrack");
    #check_material_rxy("AnalysisResults_HL_85104.root", "AnalysisResults_HL_84266.root", "qc", "LHC23d1d", "TPConlyTrack");
    #check_material_rxy("AnalysisResults_HL_88397_520472.root", "AnalysisResults_LHC23d1f_520472.root", "analysis_wo_mee", "LHC23d1f", "AnyTrack");
    #check_material_phi("AnalysisResults_HL_88397_520472.root", "AnalysisResults_LHC23d1f_520472.root", "analysis_wo_mee", "LHC23d1f", 1, 15, "AnyTrack");
    #check_material_phi("AnalysisResults_HL_88397_520472.root", "AnalysisResults_LHC23d1f_520472.root", "analysis_wo_mee", "LHC23d1f", 30, 35, "AnyTrack");
    #check_material_phi("AnalysisResults_HL_88397_520472.root", "AnalysisResults_LHC23d1f_520472.root", "analysis_wo_mee", "LHC23d1f", 35, 40, "AnyTrack");
    #check_material_phi("AnalysisResults_HL_88397_520472.root", "AnalysisResults_LHC23d1f_520472.root", "analysis_wo_mee", "LHC23d1f", 45, 55, "AnyTrack");
    #check_material_phi("AnalysisResults_HL_88397_520472.root", "AnalysisResults_LHC23d1f_520472.root", "analysis_wo_mee", "LHC23d1f", 80, 90, "AnyTrack");

    #draw_raw_yield_inc_gamma("AnalysisResults_HL_117651.root", "PCM", "LHC22q", "analysis", "");
    #draw_raw_yield_inc_gamma_low_high("AnalysisResults_HL_117651.root", "AnalysisResults_HL_117656.root", "LHC22q", "LHC22o", "PCM", "analysis", "");
    #draw_raw_yield_inc_gamma_low_high("AnalysisResults_HL_117651.root", "AnalysisResults_HL_117656.root", "LHC22q", "LHC22o", "PCM", "qc", "");
    #draw_raw_yield_inc_gamma_low_high("AnalysisResults_HL_117651.root", "AnalysisResults_HL_117656.root", "LHC22q", "LHC22o", "PCM", "qc_ITSTPC", "");
    #draw_raw_yield_inc_gamma_low_high("AnalysisResults_HL_117651.root", "AnalysisResults_HL_117656.root", "LHC22q", "LHC22o", "PCM", "qc_TPConly", "");
    #draw_raw_yield_inc_gamma_low_high("AnalysisResults_HL_117651.root", "AnalysisResults_HL_117656.root", "LHC22q", "LHC22o", "PCM", "qc_ITSonly", "");
    #draw_raw_yield_inc_gamma_low_high("AnalysisResults_HL_117568.root", "AnalysisResults_HL_117651.root", "LHC22f", "LHC22q", "PCM", "analysis", "");
    #draw_raw_yield_inc_gamma_low_high("AnalysisResults_HL_117568.root", "AnalysisResults_HL_117651.root", "LHC22f", "LHC22q", "PCM", "qc_ITSTPC", "");

    #draw_raw_yield_pi0("PCMPCM", "LHC22q", "analysis_analysis");
    #draw_raw_yield_pi0("PCMPCM", "LHC22q", "qc_qc");
    #draw_raw_yield_pi0("PCMPHOS", "LHC22q", "qc_test03");
    #draw_raw_yield_run2("PCMPCM", "LHC22o", "qc_qc");
    #draw_raw_yield_pi0_high_low("PCMPCM", "LHC22q", "LHC22o", "analysis_analysis");
    #draw_raw_yield_pi0_high_low("PCMPCM", "LHC22q", "LHC22o", "qc_ITSTPC_qc_ITSTPC");
    #draw_raw_yield_pi0_high_low("PCMPCM", "LHC22q", "LHC22o", "qc_TPConly_qc_TPConly");

    #draw_cross_section_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC23d1f_V0withAnyTrack_20230828_nsw1.root", "pi0_mc_ptspectrum_pp_13.6TeV_LHC23d1f_V0withAnyTrack_20230828_nsw1.root", "PCMPCM", "analysis_analysis");
    #compare_pi0_peak("pi0_data_ptspectrum_pp_13.6TeV_LHC23d1f_V0withAnyTrack_20230828_nsw1.root", "pi0_mc_ptspectrum_pp_13.6TeV_LHC23d1f_V0withAnyTrack_20230828_nsw1.root", "PCMPCM", "analysis_analysis");
    #compare_pi0_peak_rate("pi0_data_ptspectrum_pp_13.6TeV_LHC22q_V0withAnyTrack_20230828_nsw1.root", "pi0_data_ptspectrum_pp_13.6TeV_LHC22o_V0withAnyTrack_20230828_nsw1.root", "PCMPCM", "analysis_analysis");

    filename_data = "pi0_data_ptspectrum_pp_13.6TeV_LHC22o_20240410.root";
    filename_mc   = "pi0_mc_ptspectrum_pp_13.6TeV_LHC24b1_20240410_TPCPIDNN.root";
    ssname = "PCMPCM";
    cutname = "qc_qc";
    period_data = "LHC22o";
    period_mc = "LHC24b1";
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22o_20240401.root"       , "PCMPCM", "LHC22o" , "qc_qc", 0, "");
    compare_pi0_peak(filename_data, filename_mc, period_data, period_mc, ssname, cutname);
    #draw_cross_section_pi0(filename_data, filename_mc, period_data, ssname, cutname);
    ssname = "PCMDalitzEE";
    cutname = "qc_mee0_60_minpt100_maxeta09_tpconly_lowB";
    period_data = "LHC22o";
    period_mc = "LHC24b1";
    #draw_mgg_pi0("pi0_data_ptspectrum_pp_13.6TeV_LHC22o_20240401.root"       , ssname, period_data, cutname, 1, "");
    #compare_pi0_peak(filename_data, filename_mc, period_data, period_mc, ssname, cutname);
    #draw_cross_section_pi0(filename_data, filename_mc, period_data, ssname, cutname);


    #draw_raw_yield_run2("PCMPCM", "LHC22o", "qc_qc");

    filename_data_eta = "eta_data_ptspectrum_pp_13.6TeV_LHC22o_20240410.root";
    filename_mc_eta   = "eta_mc_ptspectrum_pp_13.6TeV_LHC24b1_20240410_TPCPIDNN.root";
    ssname = "PCMPCM";
    cutname = "qc_qc";
    period_data = "LHC22o";
    period_mc = "LHC24b1";
    #draw_mgg_eta("eta_data_ptspectrum_pp_13.6TeV_LHC22o_20240401.root"       , "PCMPCM", "LHC22o" , "qc_qc", 0, "");
    #compare_eta_peak(filename_data_eta, filename_mc_eta, period_data, period_mc, ssname, cutname);
    #draw_cross_section_eta(filename_data_eta, filename_mc_eta, period_data, ssname, cutname);

    filename_pi0_xsection = "output_pi0_cross_section.root";
    filename_eta_xsection = "output_eta_cross_section.root";
    #draw_eta_to_pi0_ratio(filename_pi0_xsection, filename_eta_xsection);

    filename_mc_pi0 = "pi0_mc_ptspectrum_pp_13.6TeV_LHC24b1_20240410_TPCPIDNN.root";
    filename_mc_eta = "eta_mc_ptspectrum_pp_13.6TeV_LHC24b1_20240410_TPCPIDNN.root";
    period = "LHC24b1";
    #draw_efficiency(filename_mc_pi0, filename_mc_eta, period, "PCMPCM", "qc_qc");

    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC22s_V0withAnyTrack.root"         , "PCMPCM", "LHC22s", "analysis_analysis", 0, "AnyTrack");
    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC22s_V0withAnyTrack_KF_PV_M0.root", "PCMPCM", "LHC22s", "analysis_analysis", 0, "KF_AnyTrack");
    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC22s_V0withAnyTrack.root"         , "PCMPCM", "LHC22s", "analysis_wo_mee_analysis_wo_mee", 0, "AnyTrack");
    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC22s_V0withAnyTrack_KF_PV_M0.root", "PCMPCM", "LHC22s", "analysis_wo_mee_analysis_wo_mee", 0, "KF_AnyTrack");
    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC22s_V0withAnyTrack.root"         , "PCMPCM", "LHC22s", "qc_qc", 0, "AnyTrack");
    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC22s_V0withAnyTrack_KF_PV_M0.root", "PCMPCM", "LHC22s", "qc_qc", 0, "KF_AnyTrack");

    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC23zz_V0withAnyTrack_125314.root"         , "PCMPCM", "LHC23zz", "analysis_analysis", 0, "AnyTrack");
    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC23zzh_V0withPCB_20231108.root"         , "PCMPCM", "LHC23zzh", "qc_ITSTPC_qc_ITSTPC", 0, "_PCB");
    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC23zzh_V0withPCB_20231108.root"         , "PCMPCM", "LHC23zzh_cpass1", "qc_qc", 0, "_PCB");
    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC23zzh_V0withPCB_20231110.root"         , "PCMPCM", "LHC23zzh", "qc_ITSTPC_qc_ITSTPC", 0, "_PCB");
    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC23zzh_V0withPCB_20231114.root"         , "PCMPCM", "LHC23zzh", "qc_ITSTPC_qc_ITSTPC", 0, "_PCB");

    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC23zzk_V0withPCB_20240108.root"         , "PCMPCM", "", "qc_qc", 0, "_PCB");
    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC23zzk_V0withPCB_20240108.root"         , "PCMPCM", "", "qc_qc", 1, "_PCB");
    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC23zzk_V0withPCB_20240108.root"         , "PCMPCM", "", "qc_qc", 12, "_PCB");
    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC23zzf_V0withPCB_20240212.root"         , "PCMPCM", "LHC23zzf", "qc_qc", 0, "_PCB"); # for Felix DPG 2024 march
    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC23zzf_V0withPCB_20240212.root"         , "PCMPCM", "LHC23zzf", "qc_qc", 1, "_PCB"); # for Felix DPG 2024 march
    #draw_mgg_pi0_pbpb("pi0_data_ptspectrum_PbPb_5.36TeV_LHC23zzf_V0withPCB_20240212.root"         , "PCMPCM", "LHC23zzf", "qc_qc", 2, "_PCB"); # for Felix DPG 2024 march


    filename_data = "pi0_data_ptspectrum_PbPb_5.36TeV_LHC23zzf_V0withPCB_20240212.root";
    filename_mc = "pi0_mc_ptspectrum_PbPb_5.36TeV_LHC23k6_TRDfix_V0withPCM_20240221.root";
    #draw_cross_section_pi0_pbpb(filename_data, filename_mc, "PCMPCM", "qc_qc");
    #compare_pi0_peak_pbpb(filename_data, filename_mc, "PCMPCM", "qc_qc");

    #draw_photon_conversion_probability("AnalysisResults_HL_195792.root", "LHC24b1", "_TPCPIDNN");
    #draw_photon_conversion_rec_efficiency("AnalysisResults_HL_195792.root", "qc", "LHC24b1", "_TPCPIDNN");
    #draw_photon_conversion_rec_efficiency_rxy("AnalysisResults_HL_195792.root", "qc", "LHC24b1", "_TPCPIDNN");
    #draw_photon_conversion_probability("AnalysisResults_HL_196011.root", "LHC24b1", "_TPCPID");
    #draw_photon_conversion_rec_efficiency("AnalysisResults_HL_196011.root", "qc", "LHC24b1", "_TPCPID");
    #draw_photon_conversion_rec_efficiency_rxy("AnalysisResults_HL_196011.root", "qc", "LHC24b1", "_TPCPID");

    #draw_tof_beta("AnalysisResults_HL_133022.root", "LHC23zo_zp", "");
    #draw_tof_beta("AnalysisResults_HL_137979.root", "LHC23zo_zp", "");
    #draw_tpc_dedx("AnalysisResults_HL_137979.root", "LHC23zo_zp", "");
    #draw_tof_beta_after("AnalysisResults_HL_137979.root", "LHC23zo_zp", "");
    #draw_tpc_dedx_after("AnalysisResults_HL_137979.root", "LHC23zo_zp", "");

    #compare_nch_pv("AnalysisResults_HL_133366.root", "AnalysisResults_HL_133022.root", "");

    #draw_photon_efficiency("photon_mc_ptspectrum_pp_13.6TeV_LHC23d1k_V0withPCB_20231108.root", "LHC23d1k", "");

    #compare_rxy("AnalysisResults_HL_137977.root", "/AnalysisResults_LHC23zzf_apass1_FT0_80.root", "qc");
    #compare_rxy("../tmp_mb/AnalysisResults_LHC22f_pass1_pcb.root", "../tmp_mb/AnalysisResults_LHC22f_pass1.root", "qc");
    #compare_rxy("../tmp_mb/AnalysisResults_LHC22f_pass1_pcb.root", "../tmp_mb/AnalysisResults_LHC22f_pass1_anytrack_withITS.root", "qc");
    #compare_rxy("AnalysisResults_HL_137977.root", "../tmp_mb/AnalysisResults_LHC22f_pass1.root", "qc");
    #compare_rxy("AnalysisResults_LHC23zzg_apass1_FT0_80.root", "AnalysisResults_LHC23zzf_apass1_FT0_80.root", "qc");
    #compare_rxy("AnalysisResults_LHC23zzg_apass1_FT0_80.root", "AnalysisResults_LHC23zzf_apass1_FT0_80.root", "qc");

