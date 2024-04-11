import os, sys, shutil
import numpy as np
import datetime
import math
import ctypes
import ROOT
#ROOT.gROOT.SetBatch(True);
from ROOT import TFile, TDirectory, TList, THnSparse, TH1, TH2, TH3, TMath, TCanvas, TPad, TPaveText, TLegend, TPython, TLine, TH1D, TF1, TGraphErrors, TGraph2DErrors, TGraphAsymmErrors, TObjArray, TArrow, TVirtualFitter
from ROOT import gROOT, gSystem, gStyle, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kMagenta, kCyan, kGray, kAzure

#_______________________________________________
def convert_d2ndptdy_to_iy(g1org):
    g1 = TGraphAsymmErrors();
    g1.SetName("g1iy_{0}".format(g1org.GetName()));
    g1.GetXaxis().SetTitle("p_{T} (GeV/c)");
    g1.GetYaxis().SetTitle("#frac{1}{2#pi} #frac{d^{2}N}{p_{T} dp_{T} dy} (GeV/c)^{-2}");

    x = ctypes.c_double(0);
    y = ctypes.c_double(0);
    n = g1org.GetN();
    for i in range(0,n):
        g1org.GetPoint(i,x,y);
        exl = g1org.GetErrorXlow(i);
        exh = g1org.GetErrorXhigh(i);
        eyl = g1org.GetErrorYlow(i);
        eyh = g1org.GetErrorYhigh(i);
        #print(x.value,y.value);
        #print(exl,exh);
        g1.SetPoint(i,x.value,y.value/TMath.TwoPi()/x.value);
        g1.SetPointError(i,exl,exh,eyl/TMath.TwoPi()/x.value, eyh/TMath.TwoPi()/x.value);
    return g1;
#_______________________________________________
def convert_th1_to_graph(h1):
    g1 = TGraphAsymmErrors();
    g1.SetName("g1_{0}".format(h1.GetName()));
    g1.GetXaxis().SetTitle(h1.GetXaxis().GetTitle());
    g1.GetYaxis().SetTitle(h1.GetYaxis().GetTitle());

    n = h1.GetNbinsX();
    for i in range(0,n):
        y = h1.GetBinContent(i+1);
        y_err = h1.GetBinError(i+1);
        x = h1.GetBinCenter(i+1);

        exl = x - h1.GetBinLowEdge(i+1);
        exh = h1.GetBinLowEdge(i+2) - x;
        #print(x,exl,exh);
        g1.SetPoint(i,x,y);
        #g1.SetPointError(i,exl,exh,y_err,y_err);
        g1.SetPointError(i,0,0,y_err,y_err);
    return g1;
#_______________________________________________
def correct_bin_shift_x_graph(h1,f1):
    #print(h1,f1);
    n = 10;
    arr_pT = np.array(h1.GetXaxis().GetXbins());
    g1 = convert_th1_to_graph(h1);

    #print(f1.GetExpFormula());
    f2 = TF1("f2","x * ({0})".format(f1.GetExpFormula()),0,10);
    f3 = TF1("f3","x*x * ({0})".format(f1.GetExpFormula()),0,10);

    list_g1 = [];
    #print(f2.GetExpFormula());
    for i in range(0,n):
        g1.Fit(f1,"QSME","",1,5);
        #print("i = ",i,f1.GetExpFormula("p"));

        for ipar in range(f1.GetNpar()):
            f2.SetParameter(ipar,f1.GetParameter(ipar));
            f2.SetParError(ipar ,f1.GetParError(ipar));
            f3.SetParameter(ipar,f1.GetParameter(ipar));
            f3.SetParError(ipar ,f1.GetParError(ipar));

        #print(f2.GetExpFormula("p"));
        #print("p0 = ", f1.GetParameter(0));
        for j in range(0,len(arr_pT)-1):
            x = ctypes.c_double(0);
            y = ctypes.c_double(0);
            g1.GetPoint(j,x,y);
            pTcenter = x.value;
            pT1      = arr_pT[j];
            pT2      = arr_pT[j+1];
            pTcenter = (pT2 + pT1)/2.;
            dpT = pT2 - pT1;
            eyl = g1.GetErrorYlow(j);
            eyh = g1.GetErrorYhigh(j);
            n_int  = f1.Integral(pT1,pT2)/dpT;
            nx_int = f2.Integral(pT1,pT2)/dpT;
            nx2_int = f3.Integral(pT1,pT2)/dpT;
            fluc = math.sqrt(abs( nx2_int - nx_int*nx_int ));
            #print(n_int,nx_int,nx2_int,fluc);
            n_int_err  = f1.IntegralError(pT1,pT2)/dpT;
            nx_int_err = f2.IntegralError(pT1,pT2)/dpT;
            #print(n_int,nx_int,n_int_err,nx_int_err);
            #print(i, n_int_err/n_int,nx_int_err/nx_int);
            x_new = nx_int/n_int; #average pT
            x_new_err = math.sqrt( math.pow(n_int_err/n_int,2) + math.pow(nx_int_err/nx_int,2) ) * x_new;
            #print(i, x_new, fluc);
            #print(i,x_new,x_new_err);

            exl = x_new - pT1;
            exh = pT2 - x_new;
            #exl = 0.;
            #exh = 0.;

            #print(i,j,pT1,pT2,x.value,x_new,y.value,y_new);
            g1.SetPoint(j,x_new,y.value);
            g1.SetPointError(j,exl,exh,eyl,eyh);
        list_g1.append(g1.Clone("g{0}".format(i)));
    return list_g1;
#_______________________________________________
def apply_bin_shift_correction_x(h1iy):
    h1iy_stat = h1iy.Clone("h1iy_stat");

    f1 = TF1("f1","x * [0]*pow(1 + x*x/[1]/[1],-[2])",0,10);#for d2N/dpt/dy
    f1.SetParameters(500,0.6,3);
    f1.SetParLimits(0,0,10000);
    f1.SetParLimits(1,0,10);
    f1.SetParLimits(2,0,10);

    h1d2ndptdy_stat = h1iy_stat.Clone("h1d2ndptdy_stat");
    h1d2ndptdy_stat.SetYTitle("#frac{d^{2}N}{dp_{T}dy} (GeV/c)^{-1}");
    for i in range(0,h1d2ndptdy_stat.GetNbinsX()):#convert invariant yield to d2ndptdy
        iy     = h1d2ndptdy_stat.GetBinContent(i+1);
        iy_err = h1d2ndptdy_stat.GetBinError(i+1);
        pT = h1d2ndptdy_stat.GetBinCenter(i+1);
        h1d2ndptdy_stat.SetBinContent(i+1, iy * pT * TMath.TwoPi());
        h1d2ndptdy_stat.SetBinError(i+1, iy * pT * TMath.TwoPi() * iy_err/iy);#keep relative unc. consistent

    g1d2ndptdy_stat = correct_bin_shift_x_graph(h1d2ndptdy_stat,f1)[-1];#last index
    g1iy_stat = convert_d2ndptdy_to_iy(g1d2ndptdy_stat);
    g1iy_stat.SetName("g1gamma_dir_stat_bin_shift_correction");
    g1iy_stat.SetTitle("direct photon after bin shift correction with stat. unc.");
    return g1iy_stat;

#_______________________________________________
if __name__ == "__main__":
    apply_bin_shift_correction_x(0,10);

