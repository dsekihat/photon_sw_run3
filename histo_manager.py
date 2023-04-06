import numpy as np
import math
import ROOT
from ROOT import TH1D, TH2D, TH3D, TMath, TH1F

#______________________________________________________________________
def rebin_histogram(h1, arrX, isdiff, is_syst=False):
    h1tmp = h1.Clone("h1tmp");
    h1rebin = h1tmp.Rebin(len(arrX)-1,"h1rebin",arrX);
    h1rebin.SetName("{0}_rebin".format(h1.GetName()));
    h1rebin.Sumw2();

    if is_syst:
        for i in range(0,h1rebin.GetNbinsX()):
            x_min = h1rebin.GetBinLowEdge(i+1);
            x_max = h1rebin.GetBinLowEdge(i+2);
            bin0 = h1.FindBin(x_min + 1e-6);
            bin1 = h1.FindBin(x_max - 1e-6);
            y = 0;
            err = 0;
            for j in range(bin0, bin1+1):
                y += h1.GetBinContent(j);
                err += h1.GetBinError(j) / h1.GetBinContent(j) * h1.GetBinContent(j);
            h1rebin.SetBinContent(i+1, y);
            h1rebin.SetBinError(i+1, err);
            #print("check 0 rel syst. = ", h1rebin.GetBinError(i+1) / h1rebin.GetBinContent(i+1));

    if isdiff and (h1rebin.Class() == TH1D.Class() or h1rebin.Class() == TH1F.Class() ) : #do you want differential histogram? e.g. dN/dpT, dN/dm.
        h1rebin.Scale(1.,"width");
    return h1rebin;
#______________________________________________________________________
def slice_histogram(h2,x0,x1,axis,isdiff):
    h1 = 0;
    delta = 1e-6;
    if "x" in axis.lower():
        bin0 = h2.GetYaxis().FindBin(x0 + delta);
        bin1 = h2.GetYaxis().FindBin(x1 - delta);
        h1 = h2.ProjectionX("h1prjx_{0}".format(h2.GetName()),bin0,bin1,"");
    elif "y" in axis.lower():
        bin0 = h2.GetXaxis().FindBin(x0 + delta);
        bin1 = h2.GetXaxis().FindBin(x1 - delta);
        h1 = h2.ProjectionY("h1prjy_{0}".format(h2.GetName()),bin0,bin1,"");

    if isdiff and h1.Class() == TH1D.Class(): #do you want differential histogram? e.g. dN/dpT, dN/dm.
        h1.Scale(1.,"width");
    return h1;
#______________________________________________________________________
def slice_profile(h2,x0,x1,axis,isdiff=False):
    h1 = 0;
    delta = 1e-6;
    if "x" in axis.lower():
        bin0 = h2.GetYaxis().FindBin(x0 + delta);
        bin1 = h2.GetYaxis().FindBin(x1 - delta);
        h1 = h2.ProfileX("h1prfx_{0}".format(h2.GetName()),bin0,bin1,"");
    elif "y" in axis.lower():
        bin0 = h2.GetXaxis().FindBin(x0 + delta);
        bin1 = h2.GetXaxis().FindBin(x1 - delta);
        h1 = h2.ProfileY("h1prfy_{0}".format(h2.GetName()),bin0,bin1,"");

    if isdiff and h1.Class() == TProfile.Class(): #do you want differential profile?
        h1.Scale(1.,"width");
    return h1;
#______________________________________________________________________
def get_bkg_subtracted(h1m_ULS_same, h1bkg):
    h1sig = h1m_ULS_same.Clone("h1sig");
    #h1sig.Sumw2();
    h1sig.Add(h1bkg,-1);
    return h1sig;
#______________________________________________________________________
def get_ratio(h1sig,h1bkg):
    h1r = h1sig.Clone("h1r");
    h1r.Reset();
    h1r.Divide(h1sig,h1bkg,1.,1.,"B");
    return h1r;
#______________________________________________________________________
def get_significance(h1sig,h1bkg):
    h1r = h1sig.Clone("h1r");
    h1r.Reset();

    n = h1r.GetNbinsX();
    for i in range(0,n):
        s = h1sig.GetBinContent(i+1);
        b = h1bkg.GetBinContent(i+1);
        s_err = h1sig.GetBinError(i+1);
        b_err = h1bkg.GetBinError(i+1);
        sig = s/sqrt(s + 2.*b);
        sig_err = sqrt( pow(pow(s+2*b,-1/2.) - s/2.*pow(s+2.*b, -3/2) ,2) *pow(s_err,2) + pow( 2*s * pow(s+2*b,-3/2.) ,2) * pow(b_err,2)  );
        h1r.SetBinContent(i+1,sig);
        h1r.SetBinError(i+1,sig_err);
    return h1r;

#______________________________________________________________________
#______________________________________________________________________
#______________________________________________________________________
#______________________________________________________________________
