// NHitsDist.C -- <nHits> distribution of the BEMC towers for every run, all runs drawn in one graph
//
// Run (compiled, ROOT 5 or 6):
//   root -l -b -q 'NHitsDist.C+("merged.root")'
//
// Input (made by StPicoTowerTest, merged over the jobs):
//   MapTower_RunID_NHits_TowerID   x = run index (line number in the run list), y = tower index (tower ID = bin number 1..4800)
//                                  content = number of events in which the tower had E_T > 0.2 GeV
//   MapTower_RunID_NEvents         number of good events per run
//
// For every run with at least minEvents good events one histogram is filled with ONE ENTRY PER TOWER:
//   perEvent = true :  <nHits> = NHits / NEvents   (fraction of the run's events in which the tower fired; comparable between runs)
//   perEvent = false:  NHits                       (raw number of hits; depends on how many events the run has)
// All histograms are then drawn on top of each other, coloured by run index (blue = first run, red = last run).
//
// Outputs (prefix outTag):
//   <outTag>_all.pdf      all runs overlaid, linear y axis
//   <outTag>_all_log.pdf  the same with a logarithmic y axis (shows the tails: hot towers)
//   <outTag>_2d.pdf       the same information as a map: run index (x) vs <nHits> (y), colour = number of towers
//   <outTag>_runs.txt     run number, run index, number of events and the mean <nHits> over all towers
//   <outTag>.root         one histogram hNHits_run<index> per run, plus the 2D map hRunVsNHits
//
// Plain C++98 style so that it also compiles with ROOT 5 / gcc 4.8. Run it compiled (with the +).

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

static double Clamp01(double x) { return x < 0. ? 0. : (x > 1. ? 1. : x); }

// "jet" colour map: x = 0 -> blue, 0.5 -> green, 1 -> red
static void JetColor(double x, double& r, double& g, double& b)
{
   r = Clamp01(1.5 - fabs(4. * x - 3.));
   g = Clamp01(1.5 - fabs(4. * x - 2.));
   b = Clamp01(1.5 - fabs(4. * x - 1.));
}

static void DrawAll(vector<TH1D*>& hs, bool logY, double ymaxData, const char* fileName,
                    int firstRunIdx, int lastRunIdx, bool perEvent)
{
   TCanvas* c = new TCanvas(logY ? "cAllLog" : "cAll", "nHits distributions of all runs", 1200, 800);
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
   (void) perEvent;
   c->SaveAs(fileName);
}

void NHitsDist(const char* inFile = "Test_Tower_Isobars_full.root",
               const char* runListFile = "/gpfs/mnt/gpfs01/star/pwg/walz/Isobars2018/data/Isobars_2018.txt",
               const char* outTag = "NHitsDist",
               int    minEvents = 0,     // runs with fewer good events are skipped (their <nHits> is too noisy)
               int    nBins = 150,         // bins of each distribution
               double xMax = 0.15,         // upper edge of the x axis (use ~ the largest NHits if perEvent = false)
               bool   perEvent = true,     // true: <nHits> = NHits/NEvents, false: raw NHits
               double deadCut = 0.02,      // towers with NHits/NEvents below this are dead  -> left out of the distributions
               double hotCut  = 0.40)      // towers with NHits/NEvents above this are hot   -> left out of the distributions
{
   gStyle->SetOptStat(0);

   TFile* f = TFile::Open(inFile);
   if (!f || f->IsZombie()) { cout << "Cannot open " << inFile << endl; return; }
   TH2D* hN  = (TH2D*) f->Get("MapTower_RunID_NHits_TowerID");
   TH1D* hEv = (TH1D*) f->Get("MapTower_RunID_NEvents");
   if (!hN || !hEv) { cout << "Missing MapTower_* histograms in " << inFile << endl; return; }

   // run numbers (optional): run index = line number, starting at 1
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

   TFile* fo = new TFile(Form("%s.root", outTag), "RECREATE");     // histograms created below belong to this file
   TH2D* h2 = new TH2D("hRunVsNHits",
                       Form("%s of the towers vs run;runIndex;%s", perEvent ? "<nHits>" : "nHits", perEvent ? "<nHits> = NHits/NEvents" : "NHits"),
                       nRuns, 0.5, nRuns + 0.5, nBins, 0., xMax);

   // pooled distribution = projection on the y axis of the <nHits> vs tower ID plot, summed over all analysed runs
   TH1D* hAll = new TH1D("hNHits_all", Form("%s of all towers in all runs (dead and hot removed);%s;tower x run entries",
                         perEvent ? "<nHits>" : "nHits", perEvent ? "<nHits> = NHits / NEvents" : "NHits"),
                         nBins, 0., xMax);
   long nDead = 0, nHot = 0;

   vector<TH1D*> hs;
   vector<int> idx;
   double vmax = 0.;
   long nAbove = 0;
   int nSkipped = 0;
   FILE* ftxt = fopen(Form("%s_runs.txt", outTag), "w");
   fprintf(ftxt, "# runId runIndex nEvents meanNHitsOverTowers\n");

   for (int r = 1; r <= nRuns; ++r) {
      const double ev = hEv->GetBinContent(r + 1);       // input axes: run index r is in bin r+1
      if (ev < minEvents) { ++nSkipped; continue; }

      TH1D* h = new TH1D(Form("hNHits_run%d", r),
                         Form("%s of the towers, run index %d;%s;towers", perEvent ? "<nHits>" : "nHits", r,
                              perEvent ? "<nHits> = NHits / NEvents" : "NHits"),
                         nBins, 0., xMax);
      double sum = 0.;
      int nGood = 0;
      for (int t = 1; t <= NT; ++t) {
         const double nh = hN->GetBinContent(r + 1, t);
         const double frac = nh / ev;
         if (frac < deadCut) { ++nDead; continue; }
         if (frac > hotCut)  { ++nHot;  continue; }
         const double v = perEvent ? frac : nh;
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
   cout << "Largest value found: " << vmax << ";  tower entries at or above xMax=" << xMax << " (overflow, not drawn): " << nAbove << endl;
   if (nAbove > 0) cout << "  -> increase xMax if you want to see them (e.g. hot towers)." << endl;

   // ---- pooled distribution: Gaussian? symmetric? ----
   cout << "Left out: dead (NHits/NEvents < " << deadCut << ") " << nDead << " tower-runs,  hot (> " << hotCut << ") " << nHot << " tower-runs" << endl;
   const double pMean = hAll->GetMean(), pRMS = hAll->GetRMS(), pSkew = hAll->GetSkewness();
   hAll->Fit("gaus", "Q", "", pMean - 2. * pRMS, pMean + 2. * pRMS);     // fit the core only, so the tails stay visible
   TF1* gf = hAll->GetFunction("gaus");
   const double gMu = gf->GetParameter(1), gSig = fabs(gf->GetParameter(2));
   double pq[3] = {0.00135, 0.5, 0.99865}, qv[3] = {0., 0., 0.};           // 0.135 % / 99.865 % = where a Gaussian has -3 sigma / +3 sigma
   hAll->GetQuantiles(3, qv, pq);
   cout << "Pooled: mean " << pMean << "  RMS " << pRMS << "  skewness " << pSkew << endl;
   cout << "Core Gaussian fit: mu " << gMu << "  sigma " << gSig << "   -> mu +- 3 sigma = [" << gMu - 3. * gSig << ", " << gMu + 3. * gSig << "]" << endl;
   cout << "Quantiles 0.135 % / 50 % / 99.865 %: " << qv[0] << " / " << qv[1] << " / " << qv[2]
        << "   (distance from median: -" << qv[1] - qv[0] << " / +" << qv[2] - qv[1] << ", equal if symmetric)" << endl;

   // colours: one per run, jet map from the first to the last run
   const int n = (int)hs.size();
   double ymaxData = 0.;
   for (int i = 0; i < n; ++i) {
      double rr, gg, bb;
      JetColor(n > 1 ? (double)i / (n - 1) : 0., rr, gg, bb);
      new TColor(3000 + i, rr, gg, bb);
      hs[i]->SetLineColor(3000 + i);
      if (hs[i]->GetMaximum() > ymaxData) ymaxData = hs[i]->GetMaximum();
   }

   fo->cd();
   fo->Write();
   cout << "Saved histograms to " << outTag << ".root" << endl;

   DrawAll(hs, false, ymaxData, Form("%s_all.pdf", outTag),     idx[0], idx[n - 1], perEvent);
   DrawAll(hs, true,  ymaxData, Form("%s_all_log.pdf", outTag), idx[0], idx[n - 1], perEvent);

   TCanvas* c3 = new TCanvas("cPooled", "pooled nHits distribution", 1400, 600);
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

   TCanvas* c2 = new TCanvas("cMap", "run vs nHits", 1200, 700);
   c2->SetLogz(1);
   c2->SetRightMargin(0.12);
   h2->SetMinimum(0.5);
   h2->Draw("COLZ");
   c2->SaveAs(Form("%s_2d.pdf", outTag));

   cout << "Wrote " << outTag << "_pooled.pdf, " << outTag << "_all.pdf, " << outTag << "_all_log.pdf, " << outTag << "_2d.pdf, "
        << outTag << "_runs.txt and " << outTag << ".root" << endl;
}