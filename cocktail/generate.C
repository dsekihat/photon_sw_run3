#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TF1.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>

using namespace std;

void generate_photon_cocktail(const int pdg, const int nev, const TF1 *f1pt, const TString suffix);

void generate()
{
  const int nev = 1e+6;

  TString filename = "~/alice/AliPhysics/PWG/Cocktail/parametrisations/pp_13TeV.root";
  TFile *rootfile = TFile::Open(filename, "READ");
  TDirectory *dir = (TDirectory*)rootfile->Get("13TeV_Comb");
  TF1 *f1pt_111 = (TF1*)dir->Get("111_pt"); //d2N/dpT/dy
  TF1 *f1pt_221 = (TF1*)dir->Get("221_pt"); //d2N/dpT/dy
  f1pt_111->SetNpx(1000);
  f1pt_221->SetNpx(1000);

  cout << f1pt_111->GetExpFormula("p") << endl;
  cout << f1pt_221->GetExpFormula("p") << endl;

  generate_photon_cocktail(111, nev, f1pt_111, "pp13TeV");
  generate_photon_cocktail(221, nev, f1pt_221, "pp13TeV");
//  generate_photon_cocktail(223, nev, f1pt_111, "pp13TeV");
//  generate_photon_cocktail(331, nev, f1pt_111, "pp13TeV");

  rootfile->Close();

}

void generate_photon_cocktail(const int pdg, const int nev, const TF1 *f1pt, const TString suffix)
{
  TFile *outfile = new TFile(Form("photon_cocktail_%s_%d.root", suffix.Data(), pdg),"RECREATE");

  TH1F *h1pt = new TH1F("h1pt","h1pt",200,0,20);
  h1pt->Sumw2();

  const double yacc = 0.9;
  const int ndaughters = 2;
  double mass_mother = 0.135;
  double mass_daughters[ndaughters] = {0,0};
  double br = 1.0;
  if(pdg == 111){
    mass_mother = 0.135;
    mass_daughters[0] = 0;
    mass_daughters[1] = 0;
    br = 0.98823; // 0.98823 +/- 0.00034
  }
  else if(pdg == 221){
    mass_mother = 0.548;
    mass_daughters[0] = 0;
    mass_daughters[1] = 0;
    br = 0.7196; // 0.7196 +/- 0.0030
  }
  else if(pdg == 223){
    mass_mother = 0.782;
    mass_daughters[0] = 0.135;
    mass_daughters[1] = 0;
    br = 0.0835; // 0.0835 +/- 0.007
  }
  else if(pdg == 331){
    mass_mother = 0.958;
    mass_daughters[0] = 0.775;
    mass_daughters[1] = 0;
    br = 0.295; // 0.295 +/- 0.004
  }
  printf("nev = %d, pdg = %d, mass_mother = %f, mass_daughters[0] = %f, mass_daughters[1] = %f, br = %f\n", nev, pdg, mass_mother, mass_daughters[0], mass_daughters[1], br);

  const double Ymax  = 1.5;
  const double pTmin = 0.0;
  const double pTmax = 20.0;

  TRandom3 *r3 = new TRandom3(pdg);
  TGenPhaseSpace event; //be careful! TGenPhaseSpace is valid only for scaler mesons
  TLorentzVector mother;
  for(int i=0;i<nev;i++){
    double pTmother = pTmin + (r3->Rndm() * (pTmax - pTmin));
    double ymother = Ymax * (2. * r3->Rndm() -1); // -Ymax < y < Ymax;
    double phimother = TMath::TwoPi() * r3->Rndm();

    double mT = TMath::Sqrt(pTmother*pTmother + mass_mother*mass_mother);
    double px = pTmother * TMath::Cos(phimother);
    double py = pTmother * TMath::Sin(phimother);
    double pz = mT * TMath::SinH(ymother);
    double E  = mT * TMath::CosH(ymother);
    mother.SetPxPyPzE(px,py,pz,E);
    event.SetDecay(mother, ndaughters, mass_daughters);

    double w = event.Generate(); // w = 1 always.
    TLorentzVector *g1 = event.GetDecay(0);
    TLorentzVector *g2 = event.GetDecay(1);
    double w_d2NdpTdy = f1pt->Eval(pTmother);
    double weight = br * w_d2NdpTdy;

    if(pdg == 111 || pdg == 221){
      if(abs(g1->Rapidity()) < yacc){
        h1pt->Fill(g1->Pt(), weight);
      }
      if(abs(g2->Rapidity()) < yacc){
        h1pt->Fill(g2->Pt(), weight);
      }
    }
    else if(pdg == 223 || pdg == 331){
      if(abs(g2->Rapidity()) < yacc){
        h1pt->Fill(g2->Pt(), weight);
      }
    }
  }//end of event loop

  outfile->WriteTObject(h1pt);
  outfile->WriteTObject(f1pt);
  outfile->Close();
  delete outfile;
  outfile = 0x0;
  delete r3;
  r3 = 0x0;
}
