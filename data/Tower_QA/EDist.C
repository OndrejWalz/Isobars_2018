// EDist.C -- mean tower energy <E> distribution for every run (all runs overlaid) and pooled over all runs
//
// Run (compiled):   root -l -b -q 'EDist.C+("Test_Tower_Isobars_full.root")'
//
// Input (made by StPicoTowerTest, merged):
//   MapTower_RunID_NHits_TowerID    x = run index, y = tower (bin number 1..4800), content = number of events with E_T > 0.2 GeV
//   MapTower_RunID_Energy_TowerID   same axes, content = sum of the energy of those hits
//   MapTower_RunID_NEvents          number of good events per run
// <E> of a tower in a run = sum(E) / NHits.  One entry per tower and run.
//
// Towers that are dead (NHits/NEvents < deadCut) or hot (> hotCut) are left out, as are towers with fewer than minHits hits
// (their <E> is too noisy).  Same cuts as in NHitsDist.C, so the same towers are removed.
//
// Outputs (prefix outTag):
//   <outTag>_pooled.pdf   <E> of all towers in all runs, linear and log y; Gaussian fit of the core, mu +- 3 sigma (red),
//                         0.135 % / 99.865 % quantiles (blue): equal distances from the median = symmetric
//   <outTag>_all.pdf, <outTag>_all_log.pdf   one histogram per run overlaid (blue = first run, red = last run)
//   <outTag>_2d.pdf       run index vs <E>, colour = number of towers
//   <outTag>_runs.txt     runId, runIndex, nEvents, mean <E> over the towers
//   <outTag>.root         all histograms

#include <vector>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
using namespace std;

static double Clamp01E(double x) { return x < 0. ? 0. : (x > 1. ? 1. : x); }

static void JetColorE(double x, double& r, double& g, double& b)
{
   r = Clamp01E(1.5 - fabs(4. * x - 3.));
   g = Clamp01E(1.5 - fabs(4. * x - 2.));
   b = Clamp01E(1.5 - fabs(4. * x - 1.));
}

static void DrawAllE(vector<TH1D*>& hs, bool logY, double ymaxData, const char* fileName, int firstRunIdx, int lastRunIdx)
{
   TCanvas* c = new TCanvas(logY ? "cAllLogE" : "cAllE", "<E> distributions of all runs", 1200, 800);
   c->SetLogy(logY ? 1 : 0);
   c->SetGrid(1, 1);
   const double ymax = logY ? 5. * ymaxData : 1.15 * ymaxData;
   const double ymin = logY ? 0.5 : 0.;
   for (size_t i = 0; i < hs.size(); ++i) {
      hs[i]->SetStats(0);
      hs[i]->SetLineWidth(1);
      hs[i]->SetMinimum(ymin);
      hs[i]->SetMaximum(ymax);
      hs[i]->Draw(i == 0 ? "hist" : "hist same");
   }
   TLatex* tx = new TLatex(0.50, 0.85, Form("%d runs:  blue = run index %d,  red = run index %d",
                                            (int)hs.size(), firstRunIdx, lastRunIdx));
   tx->SetNDC(kTRUE);
   tx->SetTextSize(0.035);
   tx->Draw();
   c->SaveAs(fileName);
}

void EDist(const char* inFile = "Test_Tower_Isobars_full.root",
           const char* runListFile = "/gpfs/mnt/gpfs01/star/pwg/walz/Isobars2018/data/Isobars_2018.txt",
           const char* outTag = "EDist",
           int    minEvents = 1,     // runs with fewer good events are skipped
           int    nBins = 200,         // bins of each distribution
           double xMax = 2.0,          // upper edge of the x axis in GeV (the overflow is reported)
           double deadCut = 0.02,      // NHits/NEvents below this: dead -> left out
           double hotCut  = 0.40,      // NHits/NEvents above this: hot  -> left out
           int    minHits = 0)        // towers with fewer hits in the run are left out (<E> too noisy)
{
   gStyle->SetOptStat(0);

   TFile* f = TFile::Open(inFile);
   if (!f || f->IsZombie()) { cout << "Cannot open " << inFile << endl; return; }
   TH2D* hN  = (TH2D*) f->Get("MapTower_RunID_NHits_TowerID");
   TH2D* hE  = (TH2D*) f->Get("MapTower_RunID_Energy_TowerID");
   TH1D* hEv = (TH1D*) f->Get("MapTower_RunID_NEvents");
   if (!hN || !hE || !hEv) { cout << "Missing MapTower_* histograms in " << inFile << endl; return; }

   vector<unsigned int> runId(1, 0);
   {
      ifstream in(runListFile);
      if (!in.is_open()) cout << "WARNING: cannot open run list " << runListFile << " - run numbers will be 0 in the text file." << endl;
      unsigned int runNum;
      while (in >> runNum) runId.push_back(runNum);
   }

   const int NT = hN->GetNbinsY();
   const int nRuns = hN->GetNbinsX() - 1;
   cout << "Input: " << nRuns << " runs, " << NT << " towers" << endl;

   TFile* fo = new TFile(Form("%s.root", outTag), "RECREATE");
   TH2D* h2 = new TH2D("hRunVsE", "<E> of the towers vs run;runIndex;<E> [GeV]", nRuns, 0.5, nRuns + 0.5, nBins, 0., xMax);
   TH1D* hAll = new TH1D("hE_all", "<E> of all towers in all runs (dead and hot removed);<E> [GeV];tower x run entries", nBins, 0., xMax);

   vector<TH1D*> hs;
   vector<int> idx;
   double vmax = 0.;
   long nAbove = 0, nDead = 0, nHot = 0, nFew = 0;
   int nSkipped = 0;
   FILE* ftxt = fopen(Form("%s_runs.txt", outTag), "w");
   fprintf(ftxt, "# runId runIndex nEvents meanE_overTowers\n");

   for (int r = 1; r <= nRuns; ++r) {
      const double ev = hEv->GetBinContent(r + 1);       // run index r is in bin r+1
      if (ev < minEvents) { ++nSkipped; continue; }

      TH1D* h = new TH1D(Form("hE_run%d", r), Form("<E> of the towers, run index %d;<E> [GeV];towers", r), nBins, 0., xMax);
      double sum = 0.;
      int nGood = 0;
      for (int t = 1; t <= NT; ++t) {
         const double nh = hN->GetBinContent(r + 1, t);
         const double frac = nh / ev;
         if (frac < deadCut) { ++nDead; continue; }
         if (frac > hotCut)  { ++nHot;  continue; }
         if (nh < minHits)   { ++nFew;  continue; }
         const double v = hE->GetBinContent(r + 1, t) / nh;
         h->Fill(v);
         hAll->Fill(v);
         h2->Fill((double)r, v);
         sum += v;
         ++nGood;
         if (v > vmax) vmax = v;
         if (v >= xMax) ++nAbove;
      }
      hs.push_back(h);
      idx.push_back(r);
      fprintf(ftxt, "%u %d %.0f %.5f\n", (r < (int)runId.size()) ? runId[r] : 0u, r, ev, nGood > 0 ? sum / nGood : 0.);
   }
   fclose(ftxt);

   if (hs.empty()) { cout << "No run has >= " << minEvents << " good events." << endl; return; }
   cout << "Runs drawn: " << hs.size() << ",  skipped (< " << minEvents << " events): " << nSkipped << endl;
   cout << "Left out tower-runs: dead " << nDead << ", hot " << nHot << ", fewer than " << minHits << " hits " << nFew << endl;
   cout << "Largest <E> found: " << vmax << " GeV;  entries at or above xMax=" << xMax << " (overflow): " << nAbove << endl;

   // ---- pooled distribution: Gaussian? symmetric? ----
   const double pMean = hAll->GetMean(), pRMS = hAll->GetRMS(), pSkew = hAll->GetSkewness();
   hAll->Fit("gaus", "Q", "", pMean - 2. * pRMS, pMean + 2. * pRMS);     // core only, so the tails stay visible
   TF1* gf = hAll->GetFunction("gaus");
   const double gMu = gf->GetParameter(1), gSig = fabs(gf->GetParameter(2));
   double pq[3] = {0.00135, 0.5, 0.99865}, qv[3] = {0., 0., 0.};
   hAll->GetQuantiles(3, qv, pq);
   cout << "Pooled: mean " << pMean << "  RMS " << pRMS << "  skewness " << pSkew << endl;
   cout << "Core Gaussian fit: mu " << gMu << "  sigma " << gSig << "   -> mu +- 3 sigma = [" << gMu - 3. * gSig << ", " << gMu + 3. * gSig << "]" << endl;
   cout << "Quantiles 0.135 % / 50 % / 99.865 %: " << qv[0] << " / " << qv[1] << " / " << qv[2]
        << "   (distance from median: -" << qv[1] - qv[0] << " / +" << qv[2] - qv[1] << ", equal if symmetric)" << endl;

   const int nr = (int)hs.size();
   double ymaxData = 0.;
   for (int i = 0; i < nr; ++i) {
      double rr, gg, bb;
      JetColorE(nr > 1 ? (double)i / (nr - 1) : 0., rr, gg, bb);
      new TColor(4000 + i, rr, gg, bb);
      hs[i]->SetLineColor(4000 + i);
      if (hs[i]->GetMaximum() > ymaxData) ymaxData = hs[i]->GetMaximum();
   }

   fo->cd();
   fo->Write();
   cout << "Saved histograms to " << outTag << ".root" << endl;

   TCanvas* c3 = new TCanvas("cPooledE", "pooled <E> distribution", 1400, 600);
   c3->Divide(2, 1);
   const double hmax = hAll->GetMaximum();
   for (int p = 1; p <= 2; ++p) {
      c3->cd(p);
      gPad->SetLogy(p == 2 ? 1 : 0);
      gPad->SetGrid(1, 1);
      const double y1 = (p == 2) ? 0.5 : 0., y2 = (p == 2) ? 5. * hmax : 1.15 * hmax;
      hAll->SetStats(0);
      hAll->SetLineColor(kBlack);
      hAll->SetMinimum(y1);
      hAll->SetMaximum(y2);
      hAll->Draw("hist");
      gf->SetRange(0., xMax);
      gf->SetLineColor(kRed);
      gf->Draw("same");
      TLine* l1 = new TLine(gMu - 3. * gSig, y1, gMu - 3. * gSig, 0.7 * y2);
      TLine* l2 = new TLine(gMu + 3. * gSig, y1, gMu + 3. * gSig, 0.7 * y2);
      TLine* l3 = new TLine(qv[0], y1, qv[0], 0.7 * y2);
      TLine* l4 = new TLine(qv[2], y1, qv[2], 0.7 * y2);
      l1->SetLineColor(kRed);  l2->SetLineColor(kRed);
      l3->SetLineColor(kBlue); l4->SetLineColor(kBlue);
      TLine* ls[4] = {l1, l2, l3, l4};
      for (int i = 0; i < 4; ++i) { ls[i]->SetLineStyle(2); ls[i]->SetLineWidth(2); ls[i]->Draw(); }
      TLegend* lg = new TLegend(0.55, 0.62, 0.89, 0.89);
      lg->SetHeader(Form("mean %.4f  RMS %.4f  skew %.2f", pMean, pRMS, pSkew));
      lg->AddEntry(hAll, "all towers, all runs", "l");
      lg->AddEntry(gf, Form("Gaussian fit of the core (#sigma=%.4f)", gSig), "l");
      lg->AddEntry(l1, "#mu #pm 3#sigma (Gaussian)", "l");
      lg->AddEntry(l3, "0.135 % / 99.865 % quantiles", "l");
      lg->Draw();
   }
   c3->SaveAs(Form("%s_pooled.pdf", outTag));

   DrawAllE(hs, false, ymaxData, Form("%s_all.pdf", outTag),     idx[0], idx[nr - 1]);
   DrawAllE(hs, true,  ymaxData, Form("%s_all_log.pdf", outTag), idx[0], idx[nr - 1]);

   TCanvas* c2 = new TCanvas("cMapE", "run vs <E>", 1200, 700);
   c2->SetLogz(1);
   c2->SetRightMargin(0.12);
   h2->SetMinimum(0.5);
   h2->Draw("COLZ");
   c2->SaveAs(Form("%s_2d.pdf", outTag));

   cout << "Wrote " << outTag << "_pooled.pdf, _all.pdf, _all_log.pdf, _2d.pdf, _runs.txt and " << outTag << ".root" << endl;
}