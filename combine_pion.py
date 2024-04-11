import numpy as np
import ROOT
from ROOT import TFile, TDirectory, TList, gStyle, gSystem, gROOT, TF1, TMath
from histo_manager import slice_histogram, rebin_histogram

def merge(filename):
    rootfile = TFile.Open(filename, "READ");
    outfile = TFile("charged_pion_PbPb5.02TeV_0090.root", "RECREATE");

    arr_cen = np.array([0,5,10,20,30,40,50,60,70,80,90], dtype=int);
    print(arr_cen);
    ncen = len(arr_cen);
    list_h1 = [];
    for icen in range(0, ncen-1):
        cen1 = arr_cen[icen];
        cen2 = arr_cen[icen+1];
        h1 = rootfile.Get("h1d2NdpTdy_Pion_Cen{0:02d}{1:02d}_stat".format(cen1, cen2));
        h1.Sumw2();
        list_h1.append(h1);

    h1all = list_h1[0].Clone("h1all");
    h1all.Sumw2();
    h1all.Reset();
    for icen in range(0, ncen-1):
        cen1 = arr_cen[icen];
        cen2 = arr_cen[icen+1];
        h1all.Add(list_h1[icen], cen2-cen1);

    for ipt in range(0, h1all.GetNbinsX()):
        d2ndptdy = h1all.GetBinContent(ipt+1);
        d2ndptdy_err = h1all.GetBinError(ipt+1);
        pt = h1all.GetBinCenter(ipt+1);
        h1all.SetBinContent(ipt+1, d2ndptdy/pt);
        h1all.SetBinError(ipt+1, d2ndptdy_err/pt);

    h1all.Scale(1/2.); #(pi+ + pi-)/2
    h1all.Scale(1/90.);
    h1all.Scale(1/TMath.TwoPi());
    h1all.SetYTitle("#frac{1}{2#pi} #frac{d^{2}N}{p_{T} dp_{T} dy} (GeV/c)^{-2}");
    outfile.WriteTObject(h1all);
    outfile.Close();

if __name__ == "__main__":
    filename = "FitPiKa_PbPb_5.02TeV.root";
    merge(filename);
