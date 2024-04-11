import sys
import ROOT
from ROOT import TH1D, TH1F, TF1, TFitResultPtr

class NMFitter:
    def __init__(self, h1same, h1mix, f1sig_name, f1bkg_name):
        self.h1same  = h1same.Clone("h1same");
        if h1mix is not None:
            self.h1mix   = h1mix.Clone("h1mix");
        else:
            self.h1mix = None;
        self.h1ratio = None;
        self.h1bkg   = None;
        self.h1sig   = None;
        self.f1sig   = None;
        self.f1bkg   = None;
        self.f1total = None;

        if "cb" in f1sig_name.lower() or "crystalball" in f1sig_name.lower():
            self.f1sig = TF1("f1sig", "crystalball(0)", 0, 1);
            if f1bkg_name == "pol1":
                self.f1bkg   = TF1("f1bkg", "pol1(0)", 0, 1);
                self.f1total = TF1("f1total", "crystalball(0) + pol1(5)", 0, 1);
            elif f1bkg_name == "pol2":
                self.f1bkg   = TF1("f1bkg", "pol2(0)", 0,1);
                self.f1total = TF1("f1total", "crystalball(0) + pol2(5)", 0, 1);
            else:
                sys.exit("Not supoorted. exit.");
        self.f1sig.SetNpx(10000);
        self.f1bkg.SetNpx(10000);
        self.f1total.SetNpx(10000);
        print("initially, ", self.f1total.GetExpFormula(""));

    def set_parameters(self, mean, sigma, cb_alpha, cb_n, is_cb_alpha_fixed=False, is_cb_n_fixed=False):
        if self.f1sig.GetNpar() == 5: #crystal ball
            self.f1sig.SetParameter(1, mean);
            self.f1sig.SetParameter(2, sigma);
            self.f1sig.SetParameter(3, cb_alpha);
            self.f1sig.SetParameter(4, cb_n);
            self.f1total.SetParameter(1, mean);
            self.f1total.SetParameter(2, sigma);
            self.f1total.SetParameter(3, cb_alpha);
            self.f1total.SetParameter(4, cb_n);
            self.f1total.SetParameter(5, 1.0);
            self.f1total.SetParameter(6, -1.0);
            if is_cb_n_fixed:
                self.f1sig.FixParameter(4, cb_n);
                self.f1total.FixParameter(4, cb_n);
            if is_cb_alpha_fixed:
                self.f1sig.FixParameter(3, cb_alpha);
                self.f1total.FixParameter(3, cb_alpha);
        self.f1total.SetParLimits(0, 1e-3, 1e+6);
        self.f1sig  .SetParLimits(0,    1, 1e+6);
        self.f1total.SetParLimits(1, mean * 0.7, mean*1.5);
        self.f1sig.  SetParLimits(1, mean * 0.7, mean*1.5);
        self.f1total.SetParLimits(2, sigma * 0.2, sigma*5.0);
        self.f1sig.  SetParLimits(2, sigma * 0.2, sigma*5.0);

    def fit_same_only(self, opt1, opt2, fitmin, fitmax):
        self.h1sig = self.h1same.Clone("h1sig");
        fitresult = self.h1sig.Fit(self.f1sig, opt1, opt2, fitmin, fitmax);
        return [fitresult, self.h1sig, None, None, self.f1sig, None, None];

    def fit(self, opt1, opt2, fitmin, fitmax):
        if self.h1mix is None:
            return self.fit_same_only(opt1, opt2, fitmin, fitmax);

        self.h1ratio = self.h1same.Clone("h1ratio");
        self.h1ratio.Reset();
        self.h1ratio.Divide(self.h1same, self.h1mix, 1., 1., "G");

        bin_mean = self.h1ratio.FindBin(self.f1total.GetParameter(1));
        height = self.h1ratio.GetBinContent(bin_mean) - 1.0;

        self.f1total.SetParameter(0, height);
        self.h1ratio.Fit(self.f1total, opt1, opt2, fitmin, fitmax);

        self.h1bkg = self.h1mix.Clone("h1bkg");
        npar_sig = self.f1sig.GetNpar();
        npar_bkg = self.f1bkg.GetNpar();
        for ip in range(npar_bkg):
            self.f1bkg.FixParameter(ip, self.f1total.GetParameter(ip + npar_sig));
        self.h1bkg.Multiply(self.f1bkg);

        self.h1sig = self.h1same.Clone("h1sig");
        self.h1sig.Add(self.h1bkg, -1);
        height = self.h1sig.GetBinContent(bin_mean);
        self.f1sig.SetParameter(0, height);

        fitresult = self.h1sig.Fit(self.f1sig, opt1, opt2, fitmin, fitmax);
        return [fitresult, self.h1sig, self.h1bkg, self.h1ratio, self.f1sig, self.f1bkg, self.f1total]

