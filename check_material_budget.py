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
from histo_manager import slice_histogram, rebin_histogram
from nm_fitter import NMFitter

# python check_material_budget.py -i AnalysisResults_HL_117411.root -t data -c configs/config_pp_13.6TeV_LHC22f_material.yml;
# python check_material_budget.py -i AnalysisResults_HL_117398.root -t mc -c configs/config_pp_13.6TeV_LHC22f_material.yml
# python check_material_budget.py -i AnalysisResults_HL_120928.root -t mc -c configs/config_pp_13.6TeV_LHC22f_material.yml

parser = argparse.ArgumentParser('script to analyze mateial budget for PCM');
parser.add_argument("-i", "--input" , default="AnalysisResults.root", type=str, help="path to the root file you want to analyze", required=True)
parser.add_argument("-c", "--config", default="config.yml", type=str, help="path to the *.yml configuration file", required=True)
parser.add_argument("-t", "--type"  , default="data" , type=str, help="run type [data or mc]", required=True)
parser.add_argument("-s", "--suffix"  , default="" , type=str, help="suffix for output file name", required=False)
args = parser.parse_args();

#________________________________________________
def get_bin_edges(axis):
    list_edges = [];
    for i in range(0, axis.GetNbins()+1):
        list_edges.append(axis.GetBinLowEdge(i+1));
    return list_edges;

#________________________________________________
def analyze_pair(list_ev, list_pair, arr_rxy_probe, arr_eta_probe, arr_pt_tag, arr_pt_probe):
    ROOT.gROOT.SetBatch(True);
    print(sys._getframe().f_code.co_name);
    list_ev.Print();
    list_pair.Print();

    outlist = THashList();
    outlist.SetOwner(True);
    outlist.SetName("outlist");

    h1mult = list_ev.FindObject("hMultNTracksPV").Clone("h1mult");
    nev = h1mult.GetEntries();
    nch = h1mult.GetMean();
    print("nev = {0:e} M events, <nch> = {1}".format(nev/1e+6, nch));
    outlist.Add(h1mult);

    hs_same = list_pair.FindObject("hs_conv_point_same").Clone("hs_same"); #Clone is called to avoid acess deleted object. 0:pT, 1:rxy, 2:phi, 3:eta
    outlist.Add(hs_same);
    h2same = hs_same.Projection(5, 0);
    outlist.Add(h2same);

    hs_mix = list_pair.FindObject("hs_conv_point_mix").Clone("hs_mix"); #Clone is called to avoid acess deleted object. 0:pT, 1:rxy, 2:phi, 3:eta
    outlist.Add(hs_mix);
    h2mix = hs_mix.Projection(5, 0);
    outlist.Add(h2mix);

    npt_probe = len(arr_pt_probe);
    for ipt in range(0, npt_probe-1):
        ptmin = arr_pt_probe[ipt];
        ptmax = arr_pt_probe[ipt+1];

        h1same = slice_histogram(h2same, ptmin, ptmax, "X", False);
        h1same.SetName("h1same_{0:d}".format(ipt));
        h1mix = slice_histogram(h2mix, ptmin, ptmax, "X", False);
        h1mix.SetName("h1mix_{0:d}".format(ipt));
        h1same.RebinX(2);
        h1mix.RebinX(2);
        outlist.Add(h1same);
        outlist.Add(h1mix);

        nmf = NMFitter(h1same, h1mix, "cb", "pol1");
        nmf.set_parameters(0.135, 0.005, 0.6, 12, True);
        [fitresult, h1sig, h1bkg, h1ratio, f1sig, f1bkg, f1total] = nmf.fit("SME", "", 0.04, 0.24);

        #h1ratio.SetName("h1ratio_{0:d}".format(ipt));
        h1sig.SetName("h1sig_{0:d}".format(ipt));
        #h1bkg.SetName("h1bkg_{0:d}".format(ipt));
        f1sig.SetName("f1sig_{0:d}".format(ipt));
        f1bkg.SetName("f1bkg_{0:d}".format(ipt));
        f1total.SetName("f1total_{0:d}".format(ipt));

        #outlist.Add(h1ratio);
        outlist.Add(h1sig);
        #outlist.Add(h1bkg);
        outlist.Add(f1sig);
        outlist.Add(f1bkg);
        outlist.Add(f1total);

    return outlist;
#________________________________________________
def analyze_single_photon(list_ev, list_v0, arr_rxy, arr_eta):
    ROOT.gROOT.SetBatch(True);
    print(sys._getframe().f_code.co_name);
    list_ev.Print();
    list_v0.Print();

    outlist = THashList();
    outlist.SetOwner(True);
    outlist.SetName("outlist");

    h1mult = list_ev.FindObject("hMultNTracksPV").Clone("h1mult");
    nev = h1mult.GetEntries();
    nch = h1mult.GetMean();
    print("nev = {0:e} M events, <nch> = {1}".format(nev/1e+6, nch));
    outlist.Add(h1mult);
    hs = list_v0.FindObject("hs_conv_point").Clone("hs"); #Clone is called to avoid acess deleted object. 0:pT, 1:rxy, 2:phi, 3:eta
    outlist.Add(hs);

    #arr_rxy = np.array(get_bin_edges(hs.GetAxis(0)), dtype=float);
    #arr_phi = np.array(get_bin_edges(hs.GetAxis(1)), dtype=float);
    #arr_eta = np.array(get_bin_edges(hs.GetAxis(2)), dtype=float);
    #print(arr_rxy);
    #print(arr_phi);
    #print(arr_eta);
    hs.Sumw2();
    hs.Scale(1./nev);
    hs.Scale(1./nch);

    #first, eta loop
    for ieta in range(0, len(arr_eta)-1):
        eta1 = arr_eta[ieta];
        eta2 = arr_eta[ieta+1];
        deta = eta2 - eta1;
        bin_eta1 = hs.GetAxis(3).FindBin(eta1 + 1e-3);
        bin_eta2 = hs.GetAxis(3).FindBin(eta2 - 1e-3);

        hs.GetAxis(3).SetRange(bin_eta1, bin_eta2); #don't call SetRangeUser
        h2rphi = hs.Projection(1, 2, "");
        h2rphi.SetName("h2rphi_eta{0:d}".format(ieta));
        h2rphi.SetTitle("conversion point #it{{r}}_{{xy}} vs. #it{{#varphi}} in {0:3.2f} < #it{{#eta}} < {1:3.2f}".format(eta1, eta2));
        h2rphi.SetXTitle("#varphi (rad.)");
        h2rphi.SetYTitle("r_{xy} (cm)");
        h2rphi.SetZTitle("#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} #frac{dN_{#gamma}}{d#eta}");
        h2rphi.Scale(1/deta);
        outlist.Add(h2rphi);

        #second, rxy loop
        for ir in range(0, len(arr_rxy)-1):
            r1 = arr_rxy[ir];
            r2 = arr_rxy[ir+1];
            dr = r2 - r1;
            bin_r1 = h2rphi.GetYaxis().FindBin(r1 + 1e-6);
            bin_r2 = h2rphi.GetYaxis().FindBin(r2 - 1e-6);
            # print("integrate over", r1 , r2, "cm", bin_r1, bin_r2);

            h1phi = h2rphi.ProjectionX("h1phi_eta{0:d}_r{1:d}".format(ieta, ir), bin_r1, bin_r2, "");
            h1phi.SetTitle("conversion point #it{{#varphi}} in {0:3.2f} < #it{{#eta}} < {1:3.2f} , {2:3.2f} < #it{{r}}_{{xy}} < {3:3.2f} cm".format(eta1, eta2, r1, r2));
            h1phi.SetXTitle("#varphi (rad.)");
            h1phi.SetYTitle("#frac{1}{<N_{ch}^{PV}>} #frac{1}{N_{ev}} #frac{d^{3}N_{#gamma}}{dr_{xy} d#eta d#varphi} (cm #upoint rad.)^{-1}");
            h1phi.Scale(dr);
            h1phi.Scale(1, "width");
            outlist.Add(h1phi);

    hs.GetAxis(3).SetRange(0,0);
    return outlist;
#________________________________________________
def run(filename, config, isMC, suffix):
    print(sys._getframe().f_code.co_name);
    print("reading...", filename);

    taskname = "material-budget";
    if ismc:
        taskname = "material-budget-mc";

    rootfile = TFile.Open(filename, "READ");
    #rootfile.ls();
    rootdire = rootfile.Get(taskname);
    #rootdire.ls();
    list_v0 = rootdire.Get("V0");
    list_pair = rootdire.Get("Pair").FindObject("PCMPCM");
    list_ev = rootdire.Get("Event").FindObject("PCMPCM");
    #list_v0.ls();
    #list_pair.ls();

    outname = "material_budget_{0}_{1}_{2}TeV_{3}{4}.root".format(args.type, config["common"]["system"], config["common"]["energy"], config["common"]["period"], suffix);
    if ismc:
        outname = "material_budget_{0}_{1}_{2}TeV_{3}{4}.root".format(args.type, config["common"]["system"], config["common"]["energy"], config["common"]["period_mc"], suffix);
    print("out file name = ", outname);
    outfile = TFile(outname, "RECREATE");

    arr_rxy = config["common"]["rxy_bin"];
    arr_eta = config["common"]["eta_bin"];
    arr_pt_tag = config["common"]["pt_tag"];
    arr_pt_probe = config["common"]["pt_probe"];

#    cutnames = config[args.type]["subsystems"][0]['cutnames']
#    tagnames = config[args.type]["subsystems"][0]['tagnames']
#    print("tagnames", tagnames); 
#    print("cutnames", cutnames); 
#    ntag = len(tagnames);
#    nc = len(cutnames);

    for ic in range(0, 1):
        #cutname = cutnames[ic];
        cutname = "wwire_ib";
        cutname = "qc";
        print("analyzing single photn ...", cutname);
        list_v0_cut = list_v0.FindObject(cutname);
        outlist_v0 = analyze_single_photon(list_ev, list_v0_cut, arr_rxy, arr_eta);
        outlist_v0.SetName(cutname);
        outlist_v0.SetOwner(True);
        outfile.WriteTObject(outlist_v0);
        outlist_v0.Clear();

#    for it in range(0, ntag):
#        tagname = tagnames[it];
#        for ic in range(0, nc):
#            probename = cutnames[ic];
#            print("analyzing pair ...", probename);
#            tpname = tagname + "_" + probename;
#            list_pair_cut = list_pair.FindObject(tpname).FindObject("nocut");
#            outlist_pair = analyze_pair(list_ev, list_pair_cut, arr_rxy, arr_eta, arr_pt_tag, arr_pt_probe);
#            outlist_pair.SetName("pair_" + tpname);
#            outlist_pair.SetOwner(True);
#            outfile.WriteTObject(outlist_pair);
#            outlist_pair.Clear();

    outfile.Close();
    rootfile.Close();
#________________________________________________
def draw_resolution(filename, period, cutname):
    ROOT.gROOT.SetBatch(False);
    rootfile = TFile.Open(filename, "READ");
    rootdire = rootfile.Get("pcm-qc-mc");
    list_v0 = rootdire.Get("V0");
    list_ev = rootdire.Get("Event");
    h1vtx = list_ev.FindObject("hZvtx_after");
    nev = h1vtx.GetEntries();
    print("nev = {0:e} M events".format(nev/1e+6));

    print("analyzing single photn ...", cutname);
    list_v0_cut = list_v0.FindObject(cutname);

    gStyle.SetPalette(55);
    gStyle.SetOptTitle(0);
    gStyle.SetOptStat(0);
    h2reso_pt  = list_v0_cut.FindObject("hPtGen_DeltaPtOverPtGen");
    h2reso_eta = list_v0_cut.FindObject("hPtGen_DeltaEta");
    h2reso_phi = list_v0_cut.FindObject("hPtGen_DeltaPhi");
    h2reso_pt .Sumw2();
    h2reso_eta.Sumw2();
    h2reso_phi.Sumw2();
    h2reso_pt .RebinY(5);
    h2reso_eta.RebinY(5);
    h2reso_phi.RebinY(5);
    nphoton = h2reso_pt .GetEntries();
    nphoton = h2reso_eta.GetEntries();
    nphoton = h2reso_phi.GetEntries();
    h2reso_pt .Scale(1/nphoton);
    h2reso_eta.Scale(1/nphoton);
    h2reso_phi.Scale(1/nphoton);
    #h2reso_pt .Scale(1/nev);
    #h2reso_eta.Scale(1/nev);
    #h2reso_phi.Scale(1/nev);
    h2reso_pt .SetContour(100);
    h2reso_eta.SetContour(100);
    h2reso_phi.SetContour(100);
    h2reso_pt .GetZaxis().SetLabelSize(0.04);
    h2reso_eta.GetZaxis().SetLabelSize(0.04);
    h2reso_phi.GetZaxis().SetLabelSize(0.04);
    h2reso_pt .GetZaxis().SetTitle("probability density");
    h2reso_eta.GetZaxis().SetTitle("probability density");
    h2reso_phi.GetZaxis().SetTitle("probability density");
    h2reso_pt .GetZaxis().SetTitleOffset(1.8);
    h2reso_eta.GetZaxis().SetTitleOffset(1.8);
    h2reso_phi.GetZaxis().SetTitleOffset(1.8);
    h2reso_pt .GetZaxis().SetRangeUser(1e-7, 1e-1);
    h2reso_eta.GetZaxis().SetRangeUser(1e-7, 1e-1);
    h2reso_phi.GetZaxis().SetRangeUser(1e-7, 1e-1);
    h2reso_pt .SetDirectory(0);
    h2reso_eta.SetDirectory(0);
    h2reso_phi.SetDirectory(0);
    ROOT.SetOwnership(h2reso_pt , False);
    ROOT.SetOwnership(h2reso_eta, False);
    ROOT.SetOwnership(h2reso_phi, False);

    c1 = TCanvas("c1", "c1", 0, 0, 1500, 500);
    c1.SetMargin(0, 0, 0, 0);
    c1.Divide(3,1, 0.01, 0.01);

    p1 = c1.cd(1);
    p1.SetMargin(0.13, 0.18, 0.11, 0.05);
    p1.SetTicks(1,1);
    p1.SetLogx(1);
    p1.SetLogz(1);
    frame1 = p1.DrawFrame(0.1, -1, 10, +1);
    frame1.GetXaxis().SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("(#it{p}_{T}^{rec} #minus #it{p}_{T}^{gen})/#it{p}_{T}^{gen}");
    frame1.GetXaxis().SetTitleOffset(1.2);
    frame1.GetYaxis().SetTitleOffset(1.5);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetMoreLogLabels(True);
    h2reso_pt.Draw("colz,same");

    p2 = c1.cd(2);
    p2.SetMargin(0.13, 0.18, 0.11, 0.05);
    p2.SetTicks(1,1);
    p2.SetLogx(1);
    p2.SetLogz(1);
    frame2 = p2.DrawFrame(0.1, -0.5, 10, +0.5);
    frame2.GetXaxis().SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    frame2.GetYaxis().SetTitle("#it{#eta}^{rec} #minus #it{#eta}^{gen}");
    frame2.GetXaxis().SetTitleOffset(1.2);
    frame2.GetYaxis().SetTitleOffset(1.5);
    frame2.GetXaxis().SetTitleSize(0.04);
    frame2.GetYaxis().SetTitleSize(0.04);
    frame2.GetXaxis().SetLabelSize(0.04);
    frame2.GetYaxis().SetLabelSize(0.04);
    frame2.GetXaxis().SetMoreLogLabels(True);
    h2reso_eta.Draw("colz,same");

    track_type_str = "";
    if "ITSTPC" in cutname:
        track_type_str = "#gamma #rightarrow e^{+}e^{#minus} with ITS-TPC matched tracks";
    elif "TPConly" in cutname:
        track_type_str = "#gamma #rightarrow e^{+}e^{#minus} with TPConly tracks";

    txt = TPaveText(0.15,0.72,0.5,0.92,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.04);
    txt.AddText("ALICE simulation");
    txt.AddText("pp at #sqrt{#it{s}} = 13.6 TeV");
    txt.AddText(period);
    txt.AddText("Momentum resolution of #gamma");
    txt.AddText(track_type_str);
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    p3 = c1.cd(3);
    p3.SetMargin(0.13, 0.18, 0.11, 0.05);
    p3.SetTicks(1,1);
    p3.SetLogx(1);
    p3.SetLogz(1);
    frame3 = p3.DrawFrame(0.1, -0.5, 10, +0.5);
    frame3.GetXaxis().SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    frame3.GetYaxis().SetTitle("#it{#varphi}^{rec} #minus #it{#varphi}^{gen} (rad.)");
    frame3.GetXaxis().SetTitleOffset(1.2);
    frame3.GetYaxis().SetTitleOffset(1.5);
    frame3.GetXaxis().SetTitleSize(0.04);
    frame3.GetYaxis().SetTitleSize(0.04);
    frame3.GetXaxis().SetLabelSize(0.04);
    frame3.GetYaxis().SetLabelSize(0.04);
    frame3.GetXaxis().SetMoreLogLabels(True);
    h2reso_phi.Draw("colz,same");

    ROOT.SetOwnership(c1, False);
    c1.Modified();
    c1.Update();
    date = datetime.date.today().strftime("%Y%m%d");
    c1.SaveAs("{0}_pp_13.6TeV_{1}_resolution_{2}.eps".format(date, period, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_resolution_{2}.pdf".format(date, period, cutname));
    c1.SaveAs("{0}_pp_13.6TeV_{1}_resolution_{2}.png".format(date, period, cutname));

    rootfile.Close();
#________________________________________________
if __name__ == "__main__":
    ismc = False;
    if args.type == "data":
        ismc = False;
    elif args.type == "mc":
        ismc = True;
    else:
        print("unknown type.sys.exit()");
        sys.exit();
    filename = args.input;

    config = "";
    with open(args.config, "r", encoding="utf-8") as config_yml:
        config = yaml.safe_load(config_yml)
    run(filename, config, ismc, args.suffix);

    if ismc:
        period = config["common"]["period_mc"]
        #draw_resolution(filename, period, "qc_ITSTPC");
        #draw_resolution(filename, period, "qc_TPConly");
        #draw_resolution(filename, period, "qc_ITSonly");

        #draw_resolution(filename, "LHC23d1k", "qc_ITSTPC");
        #draw_resolution(filename, "LHC23d1k", "qc_TPConly");
        #draw_resolution(filename, "LHC23d1k", "qc_ITSonly");
#________________________________________________
