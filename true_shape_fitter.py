import ROOT
from ROOT import TH1D, TH1F, TF1, TFitResultPtr

class TrueShapeFitter:
    def __init__(self, list_h1true):
        self.list_h1true = [];
        ntrue = len(list_h1true);
        for i in range(0, ntrue):
            h1tmp = list_h1true[i].Clone("{0}_clone".format(list_h1true[i].GetName()));
            self.list_h1true.append(h1tmp);

    def __call__(self, x, par):
        y = 0.0;
        for i in range(0, len(self.list_h1true)):
            bin_id = self.list_h1true[i].FindBin(x[0]);
            y += par[i] * self.list_h1true[i].GetBinContent(bin_id);
        return y;
