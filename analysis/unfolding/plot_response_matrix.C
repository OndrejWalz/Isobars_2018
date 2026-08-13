#include "TFile.h"
#include "TSystem.h"
#include <csignal>
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"

#include <iostream>

using std::cout;
using std::endl;

// Normalize each TRUTH column (y-axis) to unit probability, so that
// for a fixed truth bin j, sum_i P(reco=i | truth=j) = 1.
// This is the "response probability matrix" actually used in unfolding,
// as opposed to the raw yield-weighted 2D histogram.
static TH2D* MakeColumnNormalized(const TH2D* hRaw, const char* newName)
{
  TH2D* h = (TH2D*)hRaw->Clone(newName);
  h->SetDirectory(0);

  const int nx = h->GetNbinsX(); // reco
  const int ny = h->GetNbinsY(); // truth

  for (int jy = 1; jy <= ny; ++jy) {
    double colSum = 0.0;
    for (int ix = 1; ix <= nx; ++ix) colSum += h->GetBinContent(ix, jy);
    if (colSum <= 0.0) continue;
    for (int ix = 1; ix <= nx; ++ix) {
      h->SetBinContent(ix, jy, h->GetBinContent(ix, jy) / colSum);
      h->SetBinError  (ix, jy, h->GetBinError(ix, jy)   / colSum);
    }
  }
  return h;
}

// tag example: "R0.4_CENT_0_10_ptlead5"
void plot_response_matrix(const char* inFile  = "responses_embedding.root",
                           const char* tag     = "R0.4_CENT_0_10_ptlead5",
                           const char* outDir  = ".")
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1); // ROOT5-compatible "rainbow" palette (kBird enum is ROOT6-only)

  TFile* f = TFile::Open(inFile, "READ");
  if (!f || f->IsZombie()) {
    cout << "[error] cannot open " << inFile << endl;
    return;
  }

  TDirectory* d = f->GetDirectory(tag);
  if (!d) {
    cout << "[error] no directory named " << tag << " in file" << endl;
    f->Close();
    return;
  }

  // ---- raw (yield-weighted) response matrix, as filled during embedding ----
  TH2D* hRaw = dynamic_cast<TH2D*>(d->Get("hRespRecoVsTruth_full"));
  if (!hRaw) {
    cout << "[error] hRespRecoVsTruth_full not found in " << tag << endl;
    f->Close();
    return;
  }
  hRaw->SetDirectory(0);

  // ---- column-normalized probability matrix (P(reco|truth)) ----
  TH2D* hNorm = MakeColumnNormalized(hRaw, ("hRespNorm_" + TString(tag)).Data());

  // ===================== plot raw response matrix =====================
  {
    TCanvas* c1 = new TCanvas(("c_raw_" + TString(tag)).Data(), "", 900, 750);
    c1->SetRightMargin(0.14);
    c1->SetLeftMargin(0.12);
    c1->SetLogz();

    hRaw->SetTitle("");
    hRaw->GetXaxis()->SetTitle("p_{T}^{meas} (GeV/#it{c})");
    hRaw->GetXaxis()->SetRangeUser(0.0, 60.0);
    hRaw->GetYaxis()->SetRangeUser(0.0, 60.0);
    hRaw->GetYaxis()->SetTitle("p_{T}^{true} (GeV/#it{c})");
    hRaw->GetZaxis()->SetTitle("weighted yield");
    hRaw->Draw("COLZ");

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextSize(0.032);
    lat.DrawLatex(0.14, 0.93, Form("Response matrix (full stat.): %s", tag));

    c1->SaveAs(Form("%s/RespMatrix_raw_%s.pdf", outDir, tag));
    c1->SaveAs(Form("%s/RespMatrix_raw_%s.png", outDir, tag));
  }

  // ===================== plot column-normalized matrix =====================
  {
    TCanvas* c2 = new TCanvas(("c_norm_" + TString(tag)).Data(), "", 900, 750);
    c2->SetRightMargin(0.14);
    c2->SetLeftMargin(0.12);

    hNorm->SetTitle("");
    hNorm->GetXaxis()->SetTitle("p_{T}^{meas} (GeV/#it{c})");
    hNorm->GetYaxis()->SetTitle("p_{T}^{true} (GeV/#it{c})");
    hNorm->GetZaxis()->SetTitle("P(reco | truth)");
    hNorm->GetZaxis()->SetRangeUser(0.0, 1.0);
    hNorm->Draw("COLZ");

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextSize(0.032);
    lat.DrawLatex(0.14, 0.93, Form("Column-normalized response: %s", tag));

    c2->SaveAs(Form("%s/RespMatrix_norm_%s.pdf", outDir, tag));
    c2->SaveAs(Form("%s/RespMatrix_norm_%s.png", outDir, tag));
  }

  cout << "Raw response integral      : " << hRaw->Integral() << endl;

  // Deliberately do NOT call f->Close() here, and deliberately do NOT just
  // fall off the end of the function either. On this ROOT5/CINT build, the
  // "response" key in this directory holds a RooUnfoldResponse object with
  // no loaded dictionary. Anything that walks/cleans up this TFile's key
  // list -- an explicit f->Close(), or the atexit handler ROOT registers
  // that closes any TFiles still open when the interpreter session ends --
  // segfaults while touching that object. Falling through normally does
  // NOT avoid this: it just defers the crash to session exit, which is
  // easy to mistake for an unrelated failure.
  //
  // So instead we hard-kill the process right here, immediately after our
  // outputs are flushed to disk. We can't use exit()/gSystem->Exit()
  // (they run atexit handlers -- same crash) or _exit() (CINT on this
  // build fails to parse unistd.h). raise() with a numeric signal is
  // plain C89, which CINT parses fine; the macro SIGKILL from <csignal>
  // is NOT reliably resolved by CINT (it doesn't fully preprocess system
  // headers), so we use the literal value 9, which is SIGKILL on every
  // POSIX platform this runs on. SIGKILL can't be caught or handled by
  // anything, so nothing runs after this line -- including atexit cleanup.
  //
  // Practical consequence: the process exits via a signal, not a normal
  // return. In bash this shows up as exit status 137. If this macro is
  // driven from a loop over multiple tags, make sure that driver script
  // does not treat a signaled/non-zero exit as a fatal batch error.
  cout.flush();
  raise(9); // SIGKILL
}