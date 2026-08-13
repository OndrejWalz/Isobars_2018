// compare_svd_binning.C
//
// Loads hTrueTest and Unfolded_SVD_iter3 from three responses_embedding.root
// files (one per measured-binning choice), forms the ratio
//     ratio = hTrueTest / Unfolded_SVD_iter3
// for a given (R, centrality, ptlead) tag, and overlays the three ratios
// on a single canvas.
//
// Usage (inside ROOT):
//   root -l
//   .L compare_svd_binning.C+
//   compare_svd_binning();                 // uses defaults below
//
// or from the command line:
//   root -l -b -q 'compare_svd_binning.C("R0.2_CENT_0_10_ptlead0")'
//
// NOTE: the function name below (compare_svd_binning) must exactly
// match this file's base name (compare_svd_binning.C) -- when ROOT is
// invoked as `root file.C(args)`, it calls the function with the same
// name as the file, regardless of what other functions are defined
// inside it. If you rename this file, rename the function to match.
//
// Edit kFiles / kLabels below to point at your three files.

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"

#include <vector>
#include <string>
#include <iostream>

using std::string;
using std::vector;
using std::cout;
using std::endl;

// ------------------------- EDIT THESE -------------------------
//
// NOTE: using plain C-style arrays of const char* here (instead of
// vector<string> with brace-init) because some ROOT/CINT interpreter
// versions (e.g. on STAR RCF) choke on C++11 initializer-list syntax
// for vector<string> at file scope. These are converted to
// vector<string> inside the function body below.

// Paths to the three response files, one per measured-binning choice.
static const char* kFiles[3] = {
  "responses_embedding_binning1.root",
  "responses_embedding_binning2.root",
  "responses_embedding_binning3.root"
};

// Legend labels matching kFiles, in the same order.
static const char* kLabels[3] = {
  "Binning 1 ",
  "Binning 2 ",
  "Binning 3 "
};

// Output directory for the saved plot.
static const char* kOutDir = "svd_binning_comparison";

// -----------------------------------------------------------------

static const int kCols[3]  = { kBlack, kRed+1, kAzure+2 };
static const int kMarks[3] = { 20, 21, 22 };

// Same label style as the main unfolding macro: "0#minus10 %" etc.
static TString NiceCentLabel(const std::string& centToken)
{
  int a = -1, b = -1;
  if (sscanf(centToken.c_str(), "CENT_%d_%d", &a, &b) == 2 ||
      sscanf(centToken.c_str(), "MID_%d_%d",  &a, &b) == 2 ||
      sscanf(centToken.c_str(), "PERI_%d_%d", &a, &b) == 2) {
    return Form("%d#font[52]{#minus}%d %%", a, b);
  }
  return TString(centToken.c_str());
}

// Parses a tag of the form "R0.2_CENT_0_10_ptlead0" into its pieces,
// matching how the main unfolding macro builds tags
// (R + "_" + C + Form("_ptlead%.0f", cut)).
static bool ParseTag(const string& tag, double& Rval, string& centToken, double& ptleadCut)
{
  const size_t pos1 = tag.find('_');
  const size_t posPtlead = tag.find("_ptlead");
  if (pos1 == string::npos || posPtlead == string::npos || posPtlead <= pos1) {
    return false;
  }

  const string Rtok = tag.substr(0, pos1);
  if (sscanf(Rtok.c_str(), "R%lf", &Rval) != 1) return false;

  centToken = tag.substr(pos1 + 1, posPtlead - (pos1 + 1));

  const string cutTok = tag.substr(posPtlead + 7); // skip "_ptlead"
  ptleadCut = atof(cutTok.c_str());

  return true;
}

static void EnsureDir(const string& path) {
  if (gSystem->AccessPathName(path.c_str()))
    gSystem->mkdir(path.c_str(), /*recursive=*/true);
}

// The response files contain RooUnfoldResponse objects (saved alongside
// the histograms). Even though this macro never touches those objects,
// ROOT parses the file's TStreamerInfo on open, which references that
// class. Without the RooUnfold dictionary loaded, some ROOT/CINT builds
// segfault instead of just warning -- so load the library defensively
// before opening any file. Tries a couple of common library names/paths;
// safe to call even if RooUnfold is already loaded.
static void LoadRooUnfoldLib()
{
  static bool tried = false;
  if (tried) return;
  tried = true;

  // 1) explicit override via environment variable, e.g.:
  //      export ROOUNFOLD_LIB=/path/to/libRooUnfold.so
  //    set this to whatever path/name your normal unfolding setup uses
  //    (check your .rootlogon.C or the script you use to run
  //    unfold_embedding_inclusive.C for an existing gSystem->Load call
  //    or LD_LIBRARY_PATH entry pointing at RooUnfold).
  const char* envPath = gSystem->Getenv("ROOUNFOLD_LIB");
  if (envPath && envPath[0] != '\0') {
    if (gSystem->Load(envPath) >= 0) {
      cout << "[info] loaded RooUnfold from $ROOUNFOLD_LIB = " << envPath << endl;
      return;
    }
    cout << "[warn] $ROOUNFOLD_LIB is set to '" << envPath
         << "' but gSystem->Load() failed on it -- check the path." << endl;
  }

  // 2) fall back to common bare names already on the library search path
  const char* candidates[] = {
    "libRooUnfold.so",
    "libRooUnfold",
    "libRooUnfold.dylib"
  };

  for (int i = 0; i < 3; ++i) {
    if (gSystem->Load(candidates[i]) >= 0) {
      cout << "[info] loaded " << candidates[i] << endl;
      return;
    }
  }

  cout << "[warn] could not load libRooUnfold. Set the ROOUNFOLD_LIB "
       << "environment variable to its full path before running this "
       << "macro, e.g.:\n"
       << "         export ROOUNFOLD_LIB=/full/path/to/libRooUnfold.so\n"
       << "       Find the right path with something like:\n"
       << "         find / -iname 'libRooUnfold*' 2>/dev/null\n"
       << "       or check whatever setup script/.rootlogon.C you use "
       << "when running unfold_embedding_inclusive.C -- it must already "
       << "load this library successfully since that macro includes "
       << "RooUnfoldBayes.h / RooUnfoldSvd.h directly." << endl;
}

static vector<string> ToStringVector(const char* const* arr, int n)
{
  vector<string> v;
  v.reserve(n);
  for (int i = 0; i < n; ++i) v.push_back(string(arr[i]));
  return v;
}

// Loads hTrueTest and Unfolded_SVD_iter3 from a given file/tag,
// returns the ratio histogram hTrueTest / Unfolded_SVD_iter3.
// Caller owns the returned histogram (SetDirectory(0) already applied).
static TH1D* GetRatio(const string& filePath, const string& tag)
{
  LoadRooUnfoldLib();

  TFile* f = TFile::Open(filePath.c_str(), "READ");
  if (!f || f->IsZombie()) {
    cout << "[error] cannot open " << filePath << endl;
    return nullptr;
  }

  const string truePath = tag + "/hTrueTest";
  const string unfPath  = tag + "/Unfolded_SVD_iter3";

  TH1D* hTrue = dynamic_cast<TH1D*>(f->Get(truePath.c_str()));
  TH1D* hUnf  = dynamic_cast<TH1D*>(f->Get(unfPath.c_str()));

  if (!hTrue || !hUnf) {
    cout << "[error] missing " << (!hTrue ? truePath : unfPath)
         << " in " << filePath << endl;
    f->Close();
    delete f;
    return nullptr;
  }

  TH1D* ratio = (TH1D*)hTrue->Clone(("ratio_" + tag + "_" + filePath).c_str());
  ratio->SetDirectory(0);
  ratio->Divide(hUnf);

  f->Close();
  delete f;
  return ratio;
}

void compare_svd_binning(const string tag = "R0.2_CENT_0_10_ptlead0")
{
  gStyle->SetOptStat(0);
  const string outDir(kOutDir);
  EnsureDir(outDir);

  const vector<string> files  = ToStringVector(kFiles,  3);
  const vector<string> labels = ToStringVector(kLabels, 3);

  if (files.size() != labels.size() || files.empty()) {
    cout << "[error] kFiles and kLabels must be non-empty and the same size\n";
    return;
  }

  vector<TH1D*> ratios;
  ratios.reserve(files.size());

  for (size_t i = 0; i < files.size(); ++i) {
    TH1D* r = GetRatio(files[i], tag);
    ratios.push_back(r);
  }

  TCanvas* c = new TCanvas(("c_svdBinning_" + tag).c_str(), "", 800, 700);
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.04);
  c->SetTopMargin(0.06);
  c->SetBottomMargin(0.13);

 TLegend* leg = new TLegend(0.58, 0.15, 0.93, 0.35);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.030);

  TH1D* frame = nullptr;

  for (size_t i = 0; i < ratios.size(); ++i) {
    TH1D* r = ratios[i];
    if (!r) continue;

    const int col  = kCols[i % 3];
    const int mark = kMarks[i % 3];

    r->SetMarkerStyle(mark);
    r->SetMarkerColor(col);
    r->SetLineColor(col);
    r->SetMarkerSize(1.1);

    if (!frame) {
      frame = r;
      frame->SetTitle("");
      frame->GetXaxis()->SetTitle("#it{p}_{T, reco}^{corr} (GeV/#it{c})");
      frame->GetYaxis()->SetTitle("Truth / Unfolded (SVD, #it{k}_{reg}=3)");
      frame->GetYaxis()->SetRangeUser(0.7, 1.3);
      frame->GetXaxis()->SetTitleSize(0.045);
      frame->GetXaxis()->SetLabelSize(0.038);
      frame->GetYaxis()->SetTitleSize(0.045);
      frame->GetYaxis()->SetLabelSize(0.038);
      frame->GetYaxis()->SetTitleOffset(1.35);
      frame->Draw("E1");
    } else {
      r->Draw("E1 SAME");
    }

    leg->AddEntry(r, labels[i].c_str(), "lp");
  }

  if (!frame) {
    cout << "[error] no valid ratio histograms found for tag " << tag << endl;
    delete c;
    delete leg;
    return;
  }

  // reference lines at 1.0, +-5%, +-10%
  const double xmin = frame->GetXaxis()->GetXmin();
  const double xmax = frame->GetXaxis()->GetXmax();

  TLine* l0 = new TLine(xmin, 1.00, xmax, 1.00);
  l0->SetLineColor(kBlack);
  l0->SetLineStyle(1);
  l0->Draw();

  TLine* l1 = new TLine(xmin, 1.05, xmax, 1.05);
  l1->SetLineColor(kGray+1);
  l1->SetLineStyle(2);
  l1->Draw();

  TLine* l2 = new TLine(xmin, 0.95, xmax, 0.95);
  l2->SetLineColor(kGray+1);
  l2->SetLineStyle(2);
  l2->Draw();

  TLine* l3 = new TLine(xmin, 1.10, xmax, 1.10);
  l3->SetLineColor(kGray+2);
  l3->SetLineStyle(2);
  l3->Draw();

  TLine* l4 = new TLine(xmin, 0.90, xmax, 0.90);
  l4->SetLineColor(kGray+2);
  l4->SetLineStyle(2);
  l4->Draw();

  leg->Draw();

  {
    double Rval = 0.0, ptleadCut = 0.0;
    string centToken;
    const bool parsed = ParseTag(tag, Rval, centToken, ptleadCut);

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextSize(0.040);

    if (parsed) {
      lat.DrawLatex(0.16, 0.88, "p+p #otimes Au+Au  #sqrt{#it{s}_{NN}} = 200 GeV");
      lat.DrawLatex(0.16, 0.82,
                    Form("#it{R} = %.1f, %s", Rval, NiceCentLabel(centToken).Data()));
      lat.DrawLatex(0.16, 0.76, Form("#it{p}_{T}^{lead} #geq %.0f GeV/#it{c}", ptleadCut));
        lat.DrawLatex(0.16, 0.25, "THIS THESIS");
    } else {
      // fall back to raw tag if it doesn't match the expected pattern
      lat.DrawLatex(0.16, 0.88, "Au+Au  #sqrt{#it{s}_{NN}} = 200 GeV");
      lat.DrawLatex(0.16, 0.82, "Unfolding method: SVD, measured-binning comparison");
    }
  }

  const string outPdf = outDir + "/SVD_binningDependence_" + tag + ".pdf";
  const string outPng = outDir + "/SVD_binningDependence_" + tag + ".png";
  c->SaveAs(outPdf.c_str());
  c->SaveAs(outPng.c_str());

  cout << "Saved: " << outPdf << endl;

  // cleanup
  delete l0; delete l1; delete l2; delete l3; delete l4;
  delete leg;
  delete c;
  for (auto* r : ratios) delete r;
}