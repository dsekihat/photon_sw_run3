import os, sys, shutil
import numpy as np
import math
import array
import ctypes
from ROOT import TFile, TDirectory, THashList, TF1, TH1F, TF3, TH3F
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue
from histo_manager import slice_histogram, rebin_histogram, get_bkg_subtracted, get_ratio
from pair_analyzer import PairAnalyzer

class HBTAnalyzer(PairAnalyzer):
    def __init__(self):
        print("default constructor of HBTAnalyzer is called");
    def __init__(self, particle, filename, dirname, ismc):
        super().__init__(particle, filename, dirname, ismc);

    def set_arr_kt(self, arr_kt):
        super().set_arr_pt(arr_kt);

    def analyze_qinv_kt(self):
        print(sys._getframe().f_code.co_name)
        outlist = THashList();
        #outlist.SetOwner(True);
        outlist.SetName("outlist");

        hs_same = self.list_pair_ss_cut.FindObject("hs_q_same").Clone("hs_same");
        hs_mix  = self.list_pair_ss_cut.FindObject("hs_q_mix").Clone("hs_mix");
        hs_same.Sumw2();
        hs_mix .Sumw2();

        bin1_m = 0.0;
        bin2_m = 0.14;

        bin1_m = hs_same.GetAxis(1).FindBin(bin1_m + 1e-3);
        bin2_m = hs_same.GetAxis(1).FindBin(bin2_m - 1e-3);

        hs_same.GetAxis(1).SetRange(bin1_m, bin2_m);
        hs_mix .GetAxis(1).SetRange(bin1_m, bin2_m);
        h2same = hs_same.Projection(2,3);
        h2mix  = hs_mix .Projection(2,3);
        h2same.SetName("h2same");
        h2mix .SetName("h2mix");
        outlist.Add(h2same);
        outlist.Add(h2mix);

        nkt = len(self.arr_pt);
        h1lambda = TH1F("h1lambda","#lambda_{inv} vs. k_{T}"  ,nkt-1, self.arr_pt);
        h1lambda.SetXTitle("#it{k}_{T} (GeV/#it{c})");
        h1lambda.SetYTitle("#lambda_{inv}");
        h1lambda.Sumw2();
        h1R = TH1F("h1R","R_{inv} vs. k_{T}"  ,nkt-1, self.arr_pt);
        h1R.SetXTitle("#it{k}_{T} (GeV/#it{c})");
        h1R.SetYTitle("R_{inv} (fm)");
        h1R.Sumw2();

        for i in range(0, nkt-1):
            kt1 = self.arr_pt[i];
            kt2 = self.arr_pt[i+1];

            h1same = slice_histogram(h2same, kt1, kt2, "X", False);
            h1same.SetName("h1qinv_same_kt{0}".format(i));
            h1same.SetTitle("q_{{inv}}^{{same}}, {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            bw = h1same.GetBinWidth(1);
            h1same.SetXTitle("#it{q}_{inv} (GeV/#it{c})");
            h1same.SetYTitle("counts / {0:d} MeV/#it{{c}}".format(int(bw*1e+3)));
            h1mix  = slice_histogram(h2mix , kt1, kt2, "X", False);
            h1mix.SetName("h1qinv_mix_kt{0}".format(i));
            h1mix.SetTitle("q_{{inv}}^{{mix}}, {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            h1mix.SetXTitle("#it{q}_{inv} (GeV/#it{c})");
            h1mix.SetYTitle("counts / {0:d} MeV/#it{{c}}".format(int(bw*1e+3)));
            h1same.SetDirectory(0);
            h1mix .SetDirectory(0);
 
            npair_same = h1same.GetEntries();
            npair_mix  = h1mix.GetEntries();
            if npair_mix < 1e-6:
                continue;
            h1mix.Scale(npair_same/npair_mix);

            outlist.Add(h1same);
            outlist.Add(h1mix);

            h1cf = get_ratio(h1same, h1mix);
            h1cf.SetName("h1qinv_cf_kt{0}".format(i));
            h1cf.SetTitle("C(q), {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            h1cf.SetXTitle("#it{q}_{inv} (GeV/#it{c})");
            h1cf.SetYTitle("C(q)");
            h1cf .SetDirectory(0);
            outlist.Add(h1cf);

            f1 = TF1("f1_cf_kt{0}".format(i),"1 + [0] * exp(-[1]*[1]/0.197/0.197 * x*x)", 0, 0.1);
            f1.SetLineColor(kRed+1);
            f1.SetParNames("#lambda_{inv}","R_{inv} (fm)");
            f1.SetNpx(1000);
            f1.SetParameters(1,10);
            f1.SetParLimits(0,-10, +10);
            f1.SetParLimits(1,0.1,1e+3);
            h1cf.Fit(f1,"SME","",0,self.fit_max);
            outlist.Add(f1);

            la = f1.GetParameter(0);
            la_err = f1.GetParError(0);
            R = f1.GetParameter(1);
            R_err = f1.GetParError(1);

            h1lambda.SetBinContent(i+1, la);
            h1lambda.SetBinError(i+1, la_err);
            h1R.SetBinContent(i+1, R);
            h1R.SetBinError(i+1, R_err);
        outlist.Add(h1lambda);
        outlist.Add(h1R);
        return outlist;

    def analyze_3d_kt(self):
        outlist = THashList();
        #outlist.SetOwner(True);
        outlist.SetName("outlist");

        hs_same = self.list_pair_ss_cut.FindObject("hs_q_same") .Clone("hs_same");
        hs_mix  = self.list_pair_ss_cut.FindObject("hs_q_mix").Clone("hs_mix");
        hs_same.Sumw2();
        hs_mix .Sumw2();

        ndim = 4;
        arr_axis = array.array('i', [4,1,2,3]);
        hs_4_same = hs_same.Projection(ndim, arr_axis);
        hs_4_mix  = hs_mix .Projection(ndim, arr_axis);
        hs_4_same.SetName("hs_4_same");
        hs_4_mix .SetName("hs_4_mix");
        outlist.Add(hs_4_same);
        outlist.Add(hs_4_mix);
        hs_4_same.GetAxis(0).SetRange(0,0);
        hs_4_mix .GetAxis(0).SetRange(0,0);

        nkt = len(self.arr_pt);
        for i in range(0, nkt-1):
            kt1 = self.arr_pt[i];
            kt2 = self.arr_pt[i+1];

            bin0 = hs_4_same.GetAxis(0).FindBin(kt1 + 1e-3);
            bin1 = hs_4_same.GetAxis(0).FindBin(kt2 - 1e-3);
            #print(bin0, bin1);
            hs_4_same.GetAxis(0).SetRange(bin0, bin1);
            hs_4_mix .GetAxis(0).SetRange(bin0, bin1);

            #h3same = hs_4_same.Clone("hs_4_same_tmp").Projection(1,2,3);
            h3same = hs_4_same.Projection(1,2,3);
            h3same.SetName("h3_same_kt{0}".format(i));
            h3same.SetTitle("3d q^{{same}}, {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            #h3mix = hs_4_mix.Clone("hs_4_mix_tmp").Projection(1,2,3);
            h3mix = hs_4_mix.Projection(1,2,3);
            h3mix.SetName("h3_mix_kt{0}".format(i));
            h3mix.SetTitle("3d q^{{mix}}, {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            h3same.SetDirectory(0);
            h3mix .SetDirectory(0);

            h3same.RebinX(2);
            h3same.RebinY(2);
            h3same.RebinZ(2);
            h3mix.RebinX(2);
            h3mix.RebinY(2);
            h3mix.RebinZ(2);
 
            npair_same = h3same.GetEntries();
            npair_mix  = h3mix.GetEntries();
            if npair_mix < 1e-6:
                continue;
            h3mix.Scale(npair_same/npair_mix);
            outlist.Add(h3same);
            outlist.Add(h3mix);

            h3cf = get_ratio(h3same, h3mix);
            h3cf.SetName("h3_cf_kt{0}".format(i));
            h3cf.SetTitle("C(q), {0:2.1f} < #it{{k}}_{{T}} < {1:2.1f} GeV/#it{{c}}".format(kt1, kt2));
            h3cf .SetDirectory(0);

            #f3 = TF3("f3_cf_kt{0}".format(i),"1 + [0] * exp(-[1]*[1]*x*x/0.197/0.197 -[2]*[2]*y*y/0.197/0.197 -[3]*[3]*z*z/0.197/0.197)",-0.1,+0.1, -0.1,+0.1, -0.1,+0.1);
            f3 = TF3("f3_cf_kt{0}".format(i),"1 + [0] * exp(-[1]*[1]*x*x/0.197/0.197 -[2]*[2]*y*y/0.197/0.197 -[3]*[3]*z*z/0.197/0.197)", -self.fit_max, +self.fit_max, -self.fit_max, +self.fit_max, -self.fit_max, +self.fit_max);
            f3.SetParameters(1, 10, 10, 10)
            f3.SetNpx(int(1000*(self.fit_max - -self.fit_max)));
            f3.SetNpy(int(1000*(self.fit_max - -self.fit_max)));
            f3.SetNpz(int(1000*(self.fit_max - -self.fit_max)));
            f3.SetParLimits(0,   0, 1e+2);
            f3.SetParLimits(1, 0.1, 1e+2);
            f3.SetParLimits(2, 0.1, 1e+2);
            f3.SetParLimits(3, 0.1, 1e+2);
            h3cf.Fit(f3,"SME","");
            outlist.Add(h3cf);
            outlist.Add(f3);
        return outlist;
