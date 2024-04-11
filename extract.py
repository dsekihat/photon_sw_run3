import numpy as np
import ROOT
from ROOT import TFile, TDirectory, TList, gStyle, gSystem, gROOT, TF1
from histo_manager import slice_histogram, rebin_histogram, get_bkg_subtracted, get_ratio, get_corrected_bkg_simple, get_R_factor, get_corrected_bkg
gStyle.SetPalette(55);

#_____________________________________________________________________
def extract_emc(filename):
    rootfile = TFile.Open(filename, "READ");
    rootfile.ls();
    rootdire = rootfile.Get("pi0eta-to-gammagamma");
    rootdire.ls();
    list_pair = rootdire.Get("Pair");
    list_pair_ss = list_pair.FindObject("EMCEMC");
    list_pair_ss_cut = list_pair_ss.FindObject("standard_standard");
    list_pair_ss_cut.ls();
    list_pair_ss_cut.SetName("EMC");

    outfile = TFile("emc_mgg.root", "RECREATE");
    outfile.WriteTObject(list_pair_ss_cut);
    outfile.Close();
    rootfile.Close();

#_____________________________________________________________________
def extract_lmee(filename):
    rootfile = TFile.Open(filename, "READ");
    rootfile.ls();
    rootdire = rootfile.Get("analysis-same-event-pairing");
    rootdire.ls();
    rootlist = rootdire.Get("output");
    rootlist.ls();
    rootlist_uls_same  = rootlist.FindObject("PairsBarrelSEPM_lmee_pp502TeV_PID_Corr_excludePairPhiV");
    rootlist_lspp_same = rootlist.FindObject("PairsBarrelSEPP_lmee_pp502TeV_PID_Corr_excludePairPhiV");
    rootlist_lsmm_same = rootlist.FindObject("PairsBarrelSEMM_lmee_pp502TeV_PID_Corr_excludePairPhiV");
    rootlist_uls_same.ls();

    h2_uls_same  = rootlist_uls_same.FindObject("Mass_QuadDCAsigXY");
    h2_lspp_same = rootlist_lspp_same.FindObject("Mass_QuadDCAsigXY");
    h2_lsmm_same = rootlist_lsmm_same.FindObject("Mass_QuadDCAsigXY");

    bin0 = h2_uls_same.GetYaxis().FindBin(0.0 + 1e-3);
    bin1 = h2_uls_same.GetYaxis().FindBin(1.2 - 1e-3);

    h1_uls_same  = h2_uls_same.ProjectionX("h1_uls_same"  , bin0, bin1, "");
    h1_lspp_same = h2_lspp_same.ProjectionX("h1_lspp_same", bin0, bin1, "");
    h1_lsmm_same = h2_lsmm_same.ProjectionX("h1_lsmm_same", bin0, bin1, "");

    h1_uls_same .SetDirectory(0);
    h1_lspp_same.SetDirectory(0);
    h1_lsmm_same.SetDirectory(0);
    ROOT.SetOwnership(h1_uls_same, False);
    ROOT.SetOwnership(h1_lspp_same, False);
    ROOT.SetOwnership(h1_lsmm_same, False);
    h1bkg = get_corrected_bkg_simple(1.0, 0.0, h1_lspp_same, h1_lsmm_same);
    h1sig = get_bkg_subtracted(h1_uls_same, h1bkg);

    h1sig.SetDirectory(0);
    h1bkg.SetDirectory(0);
    ROOT.SetOwnership(h1sig, False);
    ROOT.SetOwnership(h1bkg, False);

    h1sig.Draw("same");
    h1sig.Scale(1, "width");

    rootfile.Close();

#_____________________________________________________________________
def extract_phiv(filename):
    gStyle.SetPalette(55);
    rootfile = TFile.Open(filename, "READ");
    rootfile.ls();
    rootdire = rootfile.Get("dalitz-ee-qc-mc");
    rootdire.ls();
    rootlist = rootdire.Get("DalitzEE");
    rootlist.ls();
    list_cut = rootlist.FindObject("mee_0_120_tpconly_wo_phiv_lowB");
    #list_cut = rootlist.FindObject("mee_0_120_tpconly_lowB");
    list_cut.ls();

    arr_mee = np.array([0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10], dtype=float);
    nm = len(arr_mee);
    outfile = TFile("phiv_lmee.root", "RECREATE");

    hs_uls = list_cut.FindObject("hs_dilepton_uls");
    hs_lspp = list_cut.FindObject("hs_dilepton_lspp");
    hs_lsmm = list_cut.FindObject("hs_dilepton_lsmm");
    
    h2_m_phiv_uls = hs_uls.Projection(0,3);
    h2_m_phiv_lspp = hs_lspp.Projection(0,3);
    h2_m_phiv_lsmm = hs_lsmm.Projection(0,3);

    h2_m_phiv_uls.SetDirectory(0);
    h2_m_phiv_lspp.SetDirectory(0);
    h2_m_phiv_lsmm.SetDirectory(0);
    ROOT.SetOwnership(h2_m_phiv_uls, False);
    ROOT.SetOwnership(h2_m_phiv_lspp, False);
    ROOT.SetOwnership(h2_m_phiv_lsmm, False);
    #h2_m_phiv_uls.Draw("colz");
    #h2_m_phiv_lspp.Draw("colz");
    #h2_m_phiv_lsmm.Draw("colz");

    outfile.WriteTObject(h2_m_phiv_uls);
    outfile.WriteTObject(h2_m_phiv_lspp);
    outfile.WriteTObject(h2_m_phiv_lsmm);


    for im in range(0, nm-1):
        m1 = arr_mee[im];
        m2 = arr_mee[im+1];
        bin1 = h2_m_phiv_uls.GetYaxis().FindBin(m1 + 1e-3);
        bin2 = h2_m_phiv_uls.GetYaxis().FindBin(m2 - 1e-3);

        h1phiv_uls = h2_m_phiv_uls.ProjectionX("h1phiv_uls_m{0}".format(im), bin1, bin2, "");
        h1phiv_lspp = h2_m_phiv_lspp.ProjectionX("h1phiv_lspp_m{0}".format(im), bin1, bin2, "");
        h1phiv_lsmm = h2_m_phiv_lsmm.ProjectionX("h1phiv_lsmm_m{0}".format(im), bin1, bin2, "");
        outfile.WriteTObject(h1phiv_uls);
        outfile.WriteTObject(h1phiv_lspp);
        outfile.WriteTObject(h1phiv_lsmm);

        h1phiv_bkg = get_corrected_bkg_simple(1.0, 0.0, h1phiv_lspp, h1phiv_lsmm);
        h1phiv_bkg.SetName("h1phiv_bkg_m{0}".format(im));
        h1phiv_sig = get_bkg_subtracted(h1phiv_uls, h1phiv_bkg);
        h1phiv_sig.SetName("h1phiv_sig_m{0}".format(im));
        outfile.WriteTObject(h1phiv_bkg);
        outfile.WriteTObject(h1phiv_sig);

    rootfile.Close();
    outfile.Close();

#_____________________________________________________________________
def extract_phiv_2d(filename):
    rootfile = TFile.Open(filename, "READ");
    rootfile.ls();
    rootdire = rootfile.Get("dalitz-ee-qc-mc");
    rootdire.ls();
    rootlist1 = rootdire.Get("DalitzEE");
    rootlist1.ls();
    rootlist2 = rootlist1.FindObject("mee_0_120_tpconly_wo_phiv_lowB");
    rootlist2.ls();

    h2_pi0    = rootlist2.FindObject("hMvsPhiV_Pi0");
    h2_photon = rootlist2.FindObject("hMvsPhiV_Photon");
    h2_pi0   .Sumw2();
    h2_photon.Sumw2();

    h2_pi0.SetDirectory(0);
    h2_photon.SetDirectory(0);
    ROOT.SetOwnership(h2_pi0, False);
    ROOT.SetOwnership(h2_photon, False);

    h2_all = h2_pi0.Clone("h2_all");
    h2_all.Add(h2_photon, 1.);
    h2_ratio = h2_all.Clone("h2_ratio");
    h2_ratio.Reset();
    h2_ratio.Divide(h2_pi0, h2_all, 1. , 1., "B");
    h2_ratio.Draw("colz");
    h2_ratio.SetContour(100);
    h2_ratio.SetDirectory(0);
    ROOT.SetOwnership(h2_ratio, False);

    f1 = TF1("f1", "[0] + [1]*x", 0, 3.2);
    f1.SetParameters(-0.0280, 0.0185);
    f1.Draw("same");
    ROOT.SetOwnership(f1, False);

    rootfile.Close();

#_____________________________________________________________________
def extract_lmmumu(filename):
    rootfile = TFile.Open(filename, "READ");
    rootfile.ls();
    rootdire = rootfile.Get("dalitz-mumu-qc");
    rootdire.ls();
    rootlist1 = rootdire.Get("DalitzMuMu");
    rootlist1.ls();
    rootlist2 = rootlist1.FindObject("mmumu_0_1100_lowB");
    rootlist2.ls();

    list_ev = rootdire.Get("Event");
    list_ev.ls();
    h1z = list_ev.FindObject("hZvtx_after");
    nev = h1z.GetEntries();

    hs_uls_same  = rootlist2.FindObject("hs_dilepton_uls_same");
    hs_lspp_same = rootlist2.FindObject("hs_dilepton_lspp_same");
    hs_lsmm_same = rootlist2.FindObject("hs_dilepton_lsmm_same");
   
    hs_uls_mix  = rootlist2.FindObject("hs_dilepton_uls_mix");
    hs_lspp_mix = rootlist2.FindObject("hs_dilepton_lspp_mix");
    hs_lsmm_mix = rootlist2.FindObject("hs_dilepton_lsmm_mix");

    h1m_uls_same  = hs_uls_same.Projection(0);
    h1m_lspp_same = hs_lspp_same.Projection(0);
    h1m_lsmm_same = hs_lsmm_same.Projection(0);
    h1m_uls_same  .SetDirectory(0);
    h1m_lspp_same .SetDirectory(0);
    h1m_lsmm_same .SetDirectory(0);
    ROOT.SetOwnership(h1m_uls_same , False);
    ROOT.SetOwnership(h1m_lspp_same, False);
    ROOT.SetOwnership(h1m_lsmm_same, False);

    h1m_uls_mix  = hs_uls_mix.Projection(0);
    h1m_lspp_mix = hs_lspp_mix.Projection(0);
    h1m_lsmm_mix = hs_lsmm_mix.Projection(0);
    h1m_uls_mix  .SetDirectory(0);
    h1m_lspp_mix .SetDirectory(0);
    h1m_lsmm_mix .SetDirectory(0);
    ROOT.SetOwnership(h1m_uls_mix , False);
    ROOT.SetOwnership(h1m_lspp_mix, False);
    ROOT.SetOwnership(h1m_lsmm_mix, False);

    h1R = get_R_factor(h1m_uls_mix, None, h1m_lspp_mix, h1m_lsmm_mix);
    h1R.SetDirectory(0);
    ROOT.SetOwnership(h1R, False);
    #h1R.Draw("");

    h1bkg = get_corrected_bkg(h1R, h1m_lspp_same, h1m_lsmm_same);
    h1bkg.SetDirectory(0);
    ROOT.SetOwnership(h1bkg, False);
   
    h1sig = get_bkg_subtracted(h1m_uls_same, h1bkg); 
    h1sig.SetDirectory(0);
    ROOT.SetOwnership(h1sig, False);

    h1m_uls_same.Scale(1, "width");
    h1bkg       .Scale(1, "width");
    h1sig       .Scale(1, "width");
    h1m_uls_same.Scale(1/nev);
    h1bkg       .Scale(1/nev);
    h1sig       .Scale(1/nev);


    h1m_uls_same.Draw("");
    h1bkg.Draw("same");
    h1sig.Draw("same");



    rootfile.Close();

#_____________________________________________________________________
if __name__ == "__main__":
    #filename = "AnalysisResults_HL_85687.root";
    #extract_emc(filename);

    #filename = "AnalysisResults_HL_115176.root";
    #extract_lmee(filename);

    #filename = "AnalysisResults_HL_134157.root";
    #extract_phiv(filename);

    #filename = "AnalysisResults_HL_135197.root";
    #filename = "AnalysisResults_lowB_LIR.root";
    #filename = "AnalysisResults_HL_136258.root";
    #filename = "AnalysisResults_HL_137134.root";
    #extract_lmmumu(filename);

    filename = "AnalysisResults_HL_136507.root";
    extract_phiv_2d(filename);

