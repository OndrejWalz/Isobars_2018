// TowerQA.C -- BEMC bad-tower QA built from the StPicoTowerTest output (merged.root)
//
// Run (compiled, ROOT 5 or 6):
//   root -l -b -q 'TowerQA.C+("merged.root","/gpfs/mnt/gpfs01/star/pwg/walz/Isobars2018/data/Isobars_2018.txt")'
//
// Input histograms (x = run index = line number in the run list, y = tower index, tower ID = bin number 1..4800):
//   MapTower_RunID_NHits_TowerID   number of events in which the tower had E_T > 0.2 GeV
//   MapTower_RunID_Energy_TowerID  sum of the tower energies for those hits
//   MapTower_RunID_NEvents         number of good events per run
//
// Procedure (per run), as in the reference study:
//   <nHits> = NHits / NEvents,   <E> = sum(E) / NHits
//   1. dead towers  : <nHits> < deadFactor * mu_all   (default deadFactor = 0.2; mu_all as below)
//                     (deadFactor <= 0: absolute cut <nHits> < deadCut instead)
//   2. hot towers   : <nHits> > hotFactor * mu_all, mu_all = mean <nHits> over ALL towers of the run (default hotFactor = 5)
//                     (hotFactor <= 0: absolute cut <nHits> > hotCut instead)
//   3. mu, sigma of <nHits> from the remaining towers
//   4. 3sigma towers: |<nHits> - mu| > 3 sigma                  -> excluded
//   5. mu_E, sigma_E of <E> from the remaining towers
//   6. 3sigmaE      : |<E> - mu_E| > 3 sigma_E                  -> excluded
//   7. every excluded tower is "bad" for that run
// EVERY threshold is a parameter of TowerQA() and is set by hand; nothing is tuned automatically:
//   deadFactor (or absolute deadCut), hotFactor (or absolute hotCut), nSigma (or fixed windows nhLo/nhHi, eLo/eHi), minEvents,
//   deadRunCut/hotRunCut/sigRunCut/sigERunCut (bad run if count > value; negative = criterion off),
//   towerMinBadRuns (+ optional towerBadFrac), twiceBadMask, towerRuleMode.
// Tower level, applied to the runs that survive the run-level cuts. A tower is "twice-bad" if it is flagged
// (dead, hot, 3sigma or 3sigmaE) in at least towerMinBadRuns (2) of those runs. towerRuleMode selects the rule:
//   0 (default, as in the reference map): a twice-bad tower is bad (red) in every kept run; a tower bad exactly ONCE is
//     considered good, and the RUN in which it is bad is excluded (grey).
//   1: a twice-bad tower is switched to good (green), but the runs in which it is bad are excluded (grey).
//      WARNING: towers that are bad in many runs (e.g. permanently dead ones) are present in every run, so this
//      rule excludes every run on real data.
//   twiceBadMask selects which categories count (1 dead + 2 hot + 4 3sigma + 8 3sigmaE; default 15 = all).
//   Bottom map: grey = excluded run (run-level cuts, the tower rule, or too few events), green = good tower, red = bad tower.
// A run is labelled bad if it has > 480 dead, > 3 hot, > 55 3sigma or > 15 3sigmaE towers (defaults = reference study), and towers that are
// bad in a large fraction of the remaining runs are labelled globally bad towers.
//
// Written in plain C++98 style so that it also compiles with ROOT 5 / gcc 4.8.

#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
using namespace std;

// ---- CORE-BEGIN (no ROOT dependence below this line, tested stand-alone) ----
enum { kGood = 0, kDead = 1, kHot = 2, kSig = 3, kSigE = 4 };

struct RunResult { int nDead, nHot, nSig, nSigE; double mu, sigma, muE, sigmaE, muAll, deadThr, hotThr; };

static double Median(vector<double> v)
{
   if (v.empty()) return 0.;
   sort(v.begin(), v.end());
   size_t n = v.size();
   return (n % 2) ? v[n/2] : 0.5 * (v[n/2 - 1] + v[n/2]);
}

static void MeanRms(const vector<double>& v, double& mean, double& rms)
{
   mean = 0.; rms = 0.;
   if (v.empty()) return;
   double s = 0.;
   for (size_t i = 0; i < v.size(); ++i) s += v[i];
   mean = s / v.size();
   double q = 0.;
   for (size_t i = 0; i < v.size(); ++i) q += (v[i] - mean) * (v[i] - mean);
   rms = (v.size() > 1) ? sqrt(q / (v.size() - 1)) : 0.;
}

// nh[t], en[t] for t = 1..NT (index 0 unused). flag[t] receives kGood/kDead/kHot/kSig/kSigE.
static void AnalyzeRun(const vector<double>& nh, const vector<double>& en,
                       double deadCut, double hotCut, double deadFactor, double hotFactor, double nSigma,
                       double nhLo, double nhHi, double eLo, double eHi,
                       unsigned char* flag, RunResult& r)
{
   // 3sigma step : if nhHi > nhLo, a tower is bad when <nHits> is outside [nhLo, nhHi]  (fixed window, nothing computed)
   //               otherwise bad when |<nHits> - mu| > nSigma * sigma                    (mu, sigma from the run)
   // 3sigmaE step: same with eLo/eHi for <E>.
   const bool fixedNh = (nhHi > nhLo), fixedE = (eHi > eLo);
   const int NT = (int)nh.size() - 1;
   r.nDead = r.nHot = r.nSig = r.nSigE = 0;
   r.mu = r.sigma = r.muE = r.sigmaE = 0.;

   // mu_all: mean <nHits> over ALL towers of this run (dead and hot ones included), as in the reference definition
   double sumAll = 0.;
   for (int t = 1; t <= NT; ++t) sumAll += nh[t];
   r.muAll  = (NT > 0) ? sumAll / NT : 0.;
   r.deadThr = (deadFactor > 0.) ? deadFactor * r.muAll : deadCut; // dead if <nHits> < deadThr
   r.hotThr  = (hotFactor  > 0.) ? hotFactor  * r.muAll : hotCut;  // hot  if <nHits> > hotThr

   vector<int> pool;
   for (int t = 1; t <= NT; ++t) {
      flag[t] = kGood;
      if (nh[t] < r.deadThr)     { flag[t] = kDead; ++r.nDead; }
      else if (nh[t] > r.hotThr) { flag[t] = kHot;  ++r.nHot;  }
      else pool.push_back(t);
   }

   vector<double> v;
   for (size_t i = 0; i < pool.size(); ++i) v.push_back(nh[pool[i]]);
   MeanRms(v, r.mu, r.sigma);

   vector<int> pool2;
   for (size_t i = 0; i < pool.size(); ++i) {
      int t = pool[i];
      bool outNh = fixedNh ? (nh[t] < nhLo || nh[t] > nhHi) : (fabs(nh[t] - r.mu) > nSigma * r.sigma);
      if (outNh) { flag[t] = kSig; ++r.nSig; }
      else pool2.push_back(t);
   }

   v.clear();
   for (size_t i = 0; i < pool2.size(); ++i) v.push_back(en[pool2[i]]);
   MeanRms(v, r.muE, r.sigmaE);

   for (size_t i = 0; i < pool2.size(); ++i) {
      int t = pool2[i];
      bool outE = fixedE ? (en[t] < eLo || en[t] > eHi) : (fabs(en[t] - r.muE) > nSigma * r.sigmaE);
      if (outE) { flag[t] = kSigE; ++r.nSigE; }
   }
}

// ---- CORE-END ----

// pointer to the per-tower flags of run r (helper: CINT cannot parse the cast inside a subscript)
static unsigned char* FlagPtr(vector<unsigned char>& f, int r, int NT)
{
   size_t off = (size_t)r * (size_t)(NT + 1);
   return &f[off];
}

// colour of a tower in a KEPT run of the bottom map: true = red (bad), false = green (good)
static bool IsRed(int mode, bool flagged, bool twiceBad)
{
   if (mode == 1) return flagged && !twiceBad;   // twice-bad towers are switched to good
   return flagged || twiceBad;                   // twice-bad towers are masked (bad) in every kept run
}

// Down-sampled copy of a +1/0/-1 map for DRAWING only (the full-resolution map stays in the ROOT file).
// A coarse cell is -1 if any fine cell in it is bad, +1 if any is good, else 0. Drawing the full
// 2000 x 4800 map costs ~10 million cells, which is far too much for a PDF (and can exhaust memory).
static TH2F* CoarseMap(TH2C* m, const char* name, int gx, int gy)
{
   const int nbx = m->GetNbinsX(), nby = m->GetNbinsY();
   const int nx = (nbx + gx - 1) / gx, ny = (nby + gy - 1) / gy;
   TH2F* c = new TH2F(name, m->GetTitle(), nx, 0.5, 0.5 + nx * gx, ny, 0.5, 0.5 + ny * gy);
   c->GetXaxis()->SetTitle("runIndex");
   c->GetYaxis()->SetTitle("towerID");
   for (int bx = 1; bx <= nbx; ++bx) {
      for (int by = 1; by <= nby; ++by) {
         double v = m->GetBinContent(bx, by);
         if (v == 0.) continue;
         int cx = (bx - 1) / gx + 1, cy = (by - 1) / gy + 1;
         double cur = c->GetBinContent(cx, cy);
         if (v < 0. || cur == 0.) c->SetBinContent(cx, cy, (cur < 0.) ? cur : v);   // -1 wins over +1
      }
   }
   return c;
}

static double LoadRun(TH2D* hN, TH2D* hE, TH1D* hEv, int r, vector<double>& nh, vector<double>& en)
{
   const int NT = hN->GetNbinsY();
   nh.clear(); en.clear();                       // (no vector::assign: CINT cannot call it)
   for (int t0 = 0; t0 <= NT; ++t0) { nh.push_back(0.); en.push_back(0.); }
   double ev = hEv->GetBinContent(r + 1);          // input axes: run index r is in bin r+1
   if (ev <= 0.) return ev;
   for (int t = 1; t <= NT; ++t) {
      double n = hN->GetBinContent(r + 1, t);
      nh[t] = n / ev;
      en[t] = (n > 0.) ? hE->GetBinContent(r + 1, t) / n : 0.;
   }
   return ev;
}

static TGraph* MakeGraph(const vector<double>& y, const unsigned char* flag, int cls, bool needPositive,
                         int color, int style, double size)
{
   vector<double> xs, ys;
   for (size_t t = 1; t < y.size(); ++t) {
      if (cls >= 0 && flag[t] != cls) continue;
      if (needPositive && y[t] <= 0.) continue;
      xs.push_back((double)t);
      ys.push_back(y[t]);
   }
   if (xs.empty()) return 0;
   TGraph* g = new TGraph((int)xs.size(), &xs[0], &ys[0]);
   g->SetMarkerColor(color);
   g->SetMarkerStyle(style);
   g->SetMarkerSize(size);
   return g;
}

static void HLine(double x1, double x2, double y, int color, int style)
{
   TLine* l = new TLine(x1, y, x2, y);
   l->SetLineColor(color);
   l->SetLineStyle(style);
   l->SetLineWidth(2);
   l->Draw();
}

// Four stacked panels (dead / hot / 3sigma / 3sigmaE towers per run) saved as one PDF page.
// withCuts = false: plain counts.  withCuts = true: the same plot plus one horizontal line per panel at the
// run-level cut (runs ABOVE the line are bad). Panels whose criterion is off (thr < 0) get no line.
static void DrawCountsPlot(TH1D** hp, const double* thr, bool withCuts, const char* fileName, int nRuns)
{
   TCanvas* c = new TCanvas(withCuts ? "cCuts" : "cCounts", withCuts ? "counts with cuts" : "counts", 1100, 1000);
   c->Divide(1, 4);
   for (int k = 0; k < 4; ++k) {
      c->cd(k + 1);
      hp[k]->Draw("hist");
      if (withCuts && thr[k] >= 0.) HLine(0.5, nRuns + 0.5, thr[k] , kRed, 1);
   }
   c->SaveAs(fileName);
}

void TowerQA(const char* inFile = "Test_Tower_Isobars_full.root",
             const char* runListFile = "/gpfs/mnt/gpfs01/star/pwg/walz/Isobars2018/data/Isobars_2018.txt",
             const char* outTag = "TowerQA",
             int    exampleRun = 958,     // run index for the example plots (<=0: a run near the middle)
             double deadCut = 0.02,      // absolute dead cut on <nHits>; only used if deadFactor <= 0
             double hotCut = 0.40,       // absolute hot cut on <nHits>; only used if hotFactor <= 0
             double nSigma = 3.0,        // sigma cut for <nHits> and <E>
             int    minEvents = 0,     // runs with fewer good events carry no information
             double deadRunCut = 480,    // a run is BAD if it has MORE than this many dead towers      (negative = off)
             double towerBadFrac = 0.0,  // optional extra condition for a bad tower: flagged in MORE than this fraction of the kept runs (0 = off; reference rule: just "at least twice")
             const char* ext = "pdf",    // image format for the example plot (counts, cuts and maps are always PDF)
             double hotRunCut = 2,       // bad run if it has MORE than this many hot towers      (negative = off)
             double sigRunCut = 125,      // bad run if it has MORE than this many 3sigma towers   (negative = off)
             double sigERunCut = 10,     // bad run if it has MORE than this many 3sigmaE towers  (negative = off)
             double nhLo = 0, double nhHi = 0,   // fixed <nHits> window for the 3sigma step (nhHi <= nhLo: use nSigma instead)
             double eLo = 0,  double eHi = 0,    // fixed <E> window in GeV for the 3sigmaE step (eHi <= eLo: use nSigma instead)
             double hotFactor = 0,             // tower is hot if <nHits> > hotFactor * mu_all  (<= 0: use absolute hotCut)
             double deadFactor = 0,           // tower is dead if <nHits> < deadFactor * mu_all (<= 0: use absolute deadCut)
             int    towerMinBadRuns = 2,         // a tower is "twice-bad" if flagged in at least this many runs that passed the run-level cuts
             int    twiceBadMask = 15,           // categories that count for the "bad in >= 2 runs" rule: 1 dead + 2 hot + 4 3sigma + 8 3sigmaE
             int    towerRuleMode = 0)           // 0: twice-bad towers red, runs with a once-bad tower excluded; 1: twice-bad towers green, their runs excluded
{
   gStyle->SetOptStat(0);

   TFile* f = TFile::Open(inFile);
   if (!f || f->IsZombie()) { cout << "Cannot open " << inFile << endl; return; }
   TH2D* hN  = (TH2D*) f->Get("MapTower_RunID_NHits_TowerID");
   TH2D* hE  = (TH2D*) f->Get("MapTower_RunID_Energy_TowerID");
   TH1D* hEv = (TH1D*) f->Get("MapTower_RunID_NEvents");
   if (!hN || !hE || !hEv) { cout << "Missing MapTower_* histograms in " << inFile << endl; return; }

   // run list: run index = line number (starting at 1), as in StPicoTowerTest::Init
   vector<unsigned int> runId(1, 0);
   {
      ifstream in(runListFile);
      if (!in.is_open()) { cout << "Cannot open run list " << runListFile << endl; return; }
      unsigned int runNum;
      while (in >> runNum) runId.push_back(runNum);
   }
   const int NT = hN->GetNbinsY();
   int nRuns = hN->GetNbinsX() - 1;
   if ((int)runId.size() - 1 != nRuns)
      cout << "WARNING: run list has " << runId.size() - 1 << " runs, histogram has " << nRuns << endl;
   if ((int)runId.size() - 1 < nRuns) nRuns = (int)runId.size() - 1;
   cout << "Input: " << nRuns << " runs, " << NT << " towers, " << hEv->GetEntries() << " good events" << endl;

   // ---------------- per-run tower classification ----------------
   vector<unsigned char> flags((size_t)(nRuns + 1) * (NT + 1), 0);
   vector<int> valid(nRuns + 1, 0), nD(nRuns + 1, 0), nH(nRuns + 1, 0), nS(nRuns + 1, 0), nSE(nRuns + 1, 0);
   vector<double> nEv(nRuns + 1, 0.), cD, cH, cS, cSE, nh, en;
   int nZeroEv = 0;
   for (int r = 1; r <= nRuns; ++r) {
      nEv[r] = LoadRun(hN, hE, hEv, r, nh, en);
      if (nEv[r] <= 0.) ++nZeroEv;
      if (nEv[r] < minEvents) continue;
      RunResult res;
      AnalyzeRun(nh, en, deadCut, hotCut, deadFactor, hotFactor, nSigma, nhLo, nhHi, eLo, eHi, FlagPtr(flags, r, NT), res);
      valid[r] = 1;
      nD[r] = res.nDead; nH[r] = res.nHot; nS[r] = res.nSig; nSE[r] = res.nSigE;
      cD.push_back(res.nDead); cH.push_back(res.nHot); cS.push_back(res.nSig); cSE.push_back(res.nSigE);
   }
   cout << "Runs with 0 good events: " << nZeroEv << ",  with 1-" << minEvents - 1 << " events: "
        << nRuns - nZeroEv - (int)cD.size() << ",  analysed (>= " << minEvents << "): " << cD.size() << endl;
   if (cD.empty()) { cout << "No run has >= " << minEvents << " events. Nothing to do." << endl; return; }

   // ---------------- run-level rule (all thresholds set by hand): bad run if count > threshold ----------------
   const char* nm[4] = {"dead", "hot", "3sigma", "3sigmaE"};
   double thr[4] = {deadRunCut, hotRunCut, sigRunCut, sigERunCut};
   vector<double>* cv[4] = {&cD, &cH, &cS, &cSE};
   cout << "Tower classification: dead <nHits> < ";
   if (deadFactor > 0.) cout << deadFactor << " x mu_all";
   else                 cout << deadCut;
   if (hotFactor > 0.) cout << ", hot <nHits> > " << hotFactor << " x mu_all";
   else                cout << ", hot <nHits> > " << hotCut;
   if (nhHi > nhLo) cout << ", 3sigma step: fixed <nHits> window [" << nhLo << ", " << nhHi << "]";
   else             cout << ", 3sigma step: |<nHits>-mu| > " << nSigma << " sigma";
   if (eHi > eLo)   cout << ", 3sigmaE step: fixed <E> window [" << eLo << ", " << eHi << "] GeV" << endl;
   else             cout << ", 3sigmaE step: |<E>-mu_E| > " << nSigma << " sigma_E" << endl;
   cout << "Bad towers per run over " << cD.size() << " analysed runs (for reference), and your run-level rule:" << endl;
   for (int k = 0; k < 4; ++k) {
      double mn = (*cv[k])[0], mx = (*cv[k])[0];
      for (size_t i = 0; i < cv[k]->size(); ++i) { if ((*cv[k])[i] < mn) mn = (*cv[k])[i]; if ((*cv[k])[i] > mx) mx = (*cv[k])[i]; }
      if (thr[k] >= 0.) cout << Form("  %-8s median %.0f  min %.0f  max %.0f   -> bad run if > %.0f", nm[k], Median(*cv[k]), mn, mx, thr[k]) << endl;
      else              cout << Form("  %-8s median %.0f  min %.0f  max %.0f   (criterion off)", nm[k], Median(*cv[k]), mn, mx) << endl;
   }

   vector<int> bad(nRuns + 1, 0), kept(nRuns + 1, 0);
   int nKept = 0, nBadRuns = 0, nLow = 0;
   for (int r = 1; r <= nRuns; ++r) {
      if (!valid[r]) { ++nLow; continue; }
      double cnt[4] = {(double)nD[r], (double)nH[r], (double)nS[r], (double)nSE[r]};
      int reason = 0;
      for (int k = 0; k < 4; ++k) if (thr[k] >= 0. && cnt[k] > thr[k]) reason |= (1 << k);
      bad[r] = reason;
      if (reason) ++nBadRuns; else { kept[r] = 1; ++nKept; }
   }
   if (nKept == 0)
      cout << "WARNING: every analysed run meets a bad-run criterion. Relax the run-level thresholds." << endl;

   // ---------------- tower level: towers that are bad in two or more runs ----------------
   // 1) in how many of the runs that passed the run-level cuts is each tower flagged?
   //    (only the categories selected by twiceBadMask count: bit 0 dead, 1 hot, 2 3sigma, 3 3sigmaE)
   vector<int> nBadT(NT + 1, 0);
   for (int r = 1; r <= nRuns; ++r) {
      if (!kept[r]) continue;
      const unsigned char* flagsRun = FlagPtr(flags, r, NT);
      for (int t = 1; t <= NT; ++t) if (flagsRun[t] && ((twiceBadMask >> (flagsRun[t] - 1)) & 1)) ++nBadT[t];
   }
   // 2) twice-bad towers: flagged in >= towerMinBadRuns of those runs. In the map they count as GOOD (green).
   vector<int> twiceBad;
   vector<unsigned char> isTwiceBad(NT + 1, 0);
   int nOnce = 0;
   for (int t = 1; t <= NT; ++t) {
      if (nBadT[t] == 1) ++nOnce;
      if (nKept > 0 && nBadT[t] >= 1 && nBadT[t] >= towerMinBadRuns && (double)nBadT[t] / nKept > towerBadFrac) {
         twiceBad.push_back(t);
         isTwiceBad[t] = 1;
      }
   }
   // 3) run exclusion
   //    mode 1: every run in which a twice-bad tower is bad is excluded
   //    mode 0: every run in which a tower is bad exactly ONCE is excluded (that tower is considered good)
   int nExclTower = 0;
   for (int r = 1; r <= nRuns; ++r) {
      if (!kept[r]) continue;
      const unsigned char* flagsRun = FlagPtr(flags, r, NT);
      bool hit = false;
      for (int t = 1; t <= NT && !hit; ++t) {
         if (!flagsRun[t] || !((twiceBadMask >> (flagsRun[t] - 1)) & 1)) continue;
         if (towerRuleMode == 1 && isTwiceBad[t]) hit = true;
         if (towerRuleMode == 0 && nBadT[t] == 1) hit = true;
      }
      if (hit) { kept[r] = 0; bad[r] |= (1 << 4); --nKept; ++nBadRuns; ++nExclTower; }
   }

   cout << "Tower rule " << towerRuleMode << ": towers bad in >= " << towerMinBadRuns << " runs: " << twiceBad.size()
        << (towerRuleMode == 1 ? " (shown as good)" : " (bad in every kept run)")
        << ";  towers bad in exactly one run: " << nOnce << endl;
   cout << "Runs excluded by the tower rule (" << (towerRuleMode == 1 ? "they contain a twice-bad tower" : "they contain a tower that is bad only once")
        << "): " << nExclTower << endl;
   cout << "Runs: " << nKept << " kept, " << nBadRuns << " excluded in total, " << nLow
        << " with < " << minEvents << " events (no information)" << endl;
   if (nKept == 0) {
      cout << "WARNING: every run is excluded by the tower rule." << endl;
      if (towerRuleMode == 1)
         cout << "         Towers that are bad in many runs (e.g. permanently dead ones) are present in every run, so rule 1 removes all runs."
                 " Use towerRuleMode = 0, or restrict the categories with twiceBadMask." << endl;
   }

   // ---------------- output histograms ----------------
   TFile* fo = new TFile(Form("%s.root", outTag), "RECREATE");
   const char* hn[4] = {"hDeadTowers", "hHotTowers", "h3sigmaTowers", "h3sigmaETowers"};
   const char* ht[4] = {"nHits/nEvents < 2%", "nHits/nEvents > 40%", "mean +/- 3sigma", "meanE +/- 3sigmaE"};
   const char* ya[4] = {"#DeadTowers", "#HotTowers", "#3sigmaTowers", "#3sigmaETowers"};
   TH1D* hc[4];
   for (int k = 0; k < 4; ++k) {
      hc[k] = new TH1D(hn[k], Form("%s;runIndex;%s", ht[k], ya[k]), nRuns, 0.5, nRuns + 0.5);
      hc[k]->SetFillColor(kBlue); hc[k]->SetLineColor(kBlue);
   }
   TH2C* map1 = new TH2C("hGoodTowers", "GoodTowers (all runs);runIndex;towerID",
                         nRuns, 0.5, nRuns + 0.5, NT, 0.5, NT + 0.5);
   TH2C* map2 = new TH2C("hGoodTowersAfterRunCut", "GoodTowers (excluded runs grey);runIndex;towerID",
                         nRuns, 0.5, nRuns + 0.5, NT, 0.5, NT + 0.5);
   TH1D* hTowerFrac = new TH1D("hTowerBadFraction", "fraction of kept runs in which the tower is bad;towerID;fraction",
                               NT, 0.5, NT + 0.5);
   for (int r = 1; r <= nRuns; ++r) {
      if (!valid[r]) continue;
      hc[0]->SetBinContent(r, nD[r]);  hc[1]->SetBinContent(r, nH[r]);
      hc[2]->SetBinContent(r, nS[r]);  hc[3]->SetBinContent(r, nSE[r]);
      const unsigned char* flagsRun = FlagPtr(flags, r, NT);
      for (int t = 1; t <= NT; ++t) {
         double val = flagsRun[t] ? -1. : 1.;
         map1->SetBinContent(r, t, val);
         // bottom map: excluded runs stay 0 (grey); in kept runs a flagged tower is red unless it is twice-bad (then green)
         if (kept[r]) map2->SetBinContent(r, t, IsRed(towerRuleMode, flagsRun[t] != 0, isTwiceBad[t] != 0) ? -1. : 1.);
      }
   }
   for (int t = 1; t <= NT; ++t) hTowerFrac->SetBinContent(t, nKept > 0 ? (double)nBadT[t] / nKept : 0.);

   // ---------------- text outputs ----------------
   {
      FILE* fs = fopen(Form("%s_runs_summary.txt", outTag), "w");
      fprintf(fs, "# runId runIndex nEvents status nDead nHot n3sigma n3sigmaE   (status: kept | BAD(reasons) | lowstat)\n");
      FILE* fb = fopen(Form("%s_bad_runs.list", outTag), "w");
      FILE* fLow = fopen(Form("%s_lowstat_runs.list", outTag), "w");
      FILE* fg = fopen(Form("%s_good_runs.list", outTag), "w");
      for (int r = 1; r <= nRuns; ++r) {
         string st = "kept";
         if (!valid[r]) st = "lowstat";
         else if (bad[r]) {
            st = "BAD(";
            for (int k = 0; k < 4; ++k) if (bad[r] & (1 << k)) { st += nm[k]; st += " "; }
            if (bad[r] & (1 << 4)) st += (towerRuleMode == 1 ? "twiceBadTower " : "onceBadTower ");
            st += ")";
         }
         fprintf(fs, "%u %d %.0f %s %d %d %d %d\n", runId[r], r, nEv[r], st.c_str(), nD[r], nH[r], nS[r], nSE[r]);
         if (!valid[r]) { fprintf(fLow, "%u\n", runId[r]); continue; }
         if (bad[r]) fprintf(fb, "%u\n", runId[r]); else fprintf(fg, "%u\n", runId[r]);
      }
      fclose(fs); fclose(fb); fclose(fLow); fclose(fg);

      FILE* ft = fopen(Form("%s_bad_towers.list", outTag), "w");
      fprintf(ft, "# towerID nBadRuns   (twice-bad towers: flagged in >= %d of the runs that passed the run-level cuts; rule %d: %s)\n", towerMinBadRuns, towerRuleMode,
              towerRuleMode == 1 ? "shown as good, their runs are excluded" : "bad in every kept run");
      for (int ig = 0; ig < (int)twiceBad.size(); ++ig)
         fprintf(ft, "%d %d\n", twiceBad[ig], nBadT[twiceBad[ig]]);
      fclose(ft);

      FILE* fp = fopen(Form("%s_bad_towers_perRun.txt", outTag), "w");
      fprintf(fp, "# runId nBad towerIDs...   (kept runs only; the red cells of the bottom map for that run)\n");
      for (int r = 1; r <= nRuns; ++r) {
         if (!kept[r]) continue;
         const unsigned char* flagsRun2 = FlagPtr(flags, r, NT);
         int n = 0;
         for (int t = 1; t <= NT; ++t) if (IsRed(towerRuleMode, flagsRun2[t] != 0, isTwiceBad[t] != 0)) ++n;
         fprintf(fp, "%u %d", runId[r], n);
         for (int t = 1; t <= NT; ++t) if (IsRed(towerRuleMode, flagsRun2[t] != 0, isTwiceBad[t] != 0)) fprintf(fp, " %d", t);
         fprintf(fp, "\n");
      }
      fclose(fp);
   }

   fo->cd();
   fo->Write();      // histograms are saved before any plotting
   cout << "Saved histograms to " << outTag << ".root" << endl;

   // ---------------- plots 1a/1b: counts per run, without and with the applied cuts (PDF) ----------------
   {
      TH1D* hp[4] = {hc[0], hc[1], hc[2], hc[3]};
      for (int k = 0; k < 4; ++k) {          // same y range in both files, and room for the cut line
         double ymax = 1.1 * hp[k]->GetMaximum();
         if (thr[k] >= 0. && 1.1 * thr[k] > ymax) ymax = 1.1 * thr[k];
         hp[k]->SetMinimum(0.);
         hp[k]->SetMaximum(ymax > 0. ? ymax : 1.);
      }
      DrawCountsPlot(hp, thr, false, Form("%s_counts.pdf", outTag), nRuns);   // counts only
      DrawCountsPlot(hp, thr, true,  Form("%s_cuts.pdf",   outTag), nRuns);   // counts + one cut line per panel
   }

   // ---------------- plot 2: good-tower maps ----------------
   {
      int pal[3] = {kRed, kGray, kGreen + 2};
      gStyle->SetPalette(3, pal);
      TCanvas* c = new TCanvas("cMaps", "good tower maps", 1400, 1000);
      c->Divide(1, 2);
      TH2C* mm[2] = {map1, map2};
      const char* cn[2] = {"hGoodTowers_draw", "hGoodTowersAfterRunCut_draw"};
      const int gx = (nRuns + 999) / 1000, gy = (NT + 799) / 800;     // <= 1000 x 800 cells
      for (int im = 0; im < 2; ++im) {
         c->cd(im + 1);
         gPad->SetRightMargin(0.12);
         TH2F* cm = CoarseMap(mm[im], cn[im], gx, gy);
         cm->SetMinimum(-1.5); cm->SetMaximum(1.5); cm->SetContour(3);
         cm->SetEntries(1.);   // a histogram with no entries is not drawn at all; this makes an all-grey map visible
         cm->Draw("COLZ");
      }
      c->SaveAs(Form("%s_maps.pdf", outTag));
   }

   // ---------------- plot 3: example run (nHits and mean energy vs tower) ----------------
   {
      int ex = exampleRun;
      if (ex <= 0 || ex > nRuns || !valid[ex]) {
         ex = -1;
         for (int d = 0; d < nRuns && ex < 0; ++d) {
            int a = nRuns / 2 + d, b = nRuns / 2 - d;
            if (a >= 1 && a <= nRuns && kept[a]) ex = a;
            else if (b >= 1 && b <= nRuns && kept[b]) ex = b;
         }
      }
      if (ex > 0) {
         vector<unsigned char> fe(NT + 1, 0);
         RunResult res;
         LoadRun(hN, hE, hEv, ex, nh, en);
         AnalyzeRun(nh, en, deadCut, hotCut, deadFactor, hotFactor, nSigma, nhLo, nhHi, eLo, eHi, &fe[0], res);
         cout << "Example run index " << ex << " (run " << runId[ex] << "): mu_all=" << res.muAll << " dead cut=" << res.deadThr << " hot cut=" << res.hotThr << " mu=" << res.mu << " sigma=" << res.sigma
              << " muE=" << res.muE << " sigmaE=" << res.sigmaE << endl;

         TCanvas* c = new TCanvas("cExample", "example run", 1000, 900);
         c->Divide(1, 2);

         c->cd(1);
         TH1F* fr1 = gPad->DrawFrame(0, 0, NT + 1, 1.2);
         fr1->SetTitle(Form("Example of nHits for run index %d (run %u);Tower ID;<nHits>", ex, runId[ex]));
         TLegend* l1 = new TLegend(0.12, 0.62, 0.36, 0.89);
         TGraph* g = MakeGraph(nh, &fe[0], -1, false, kBlack, 20, 0.4);       if (g) { g->Draw("P same"); l1->AddEntry(g, "nHits/NEvents", "p"); }
         g = MakeGraph(nh, &fe[0], kHot, false, kRed, 20, 0.6);               if (g) { g->Draw("P same"); l1->AddEntry(g, "Hot towers", "p"); }
         g = MakeGraph(nh, &fe[0], kDead, false, kBlue, 20, 0.6);             if (g) { g->Draw("P same"); l1->AddEntry(g, "Dead towers", "p"); }
         g = MakeGraph(nh, &fe[0], kSig, false, kGreen + 1, 20, 0.6);         if (g) { g->Draw("P same"); l1->AddEntry(g, "3sigma towers", "p"); }
         HLine(0, NT + 1, res.mu, kRed, 3);
         HLine(0, NT + 1, res.mu + nSigma * res.sigma, kGreen + 1, 3);
         HLine(0, NT + 1, res.mu - nSigma * res.sigma, kGreen + 1, 3);
         HLine(0, NT + 1, res.deadThr, kMagenta, 3);
         HLine(0, NT + 1, res.hotThr, kMagenta, 3);
         l1->Draw();

         c->cd(2);
         double ymax = res.muE + 4. * res.sigmaE;
         if (ymax < 1.) ymax = 1.;
         TH1F* fr2 = gPad->DrawFrame(0, 0, NT + 1, 1.5 * ymax);
         fr2->SetTitle("Example of mean energy;Tower ID;<E> [GeV]");
         TLegend* l2 = new TLegend(0.65, 0.62, 0.89, 0.89);
         g = MakeGraph(en, &fe[0], -1, true, kBlack, 20, 0.4);                if (g) { g->Draw("P same"); l2->AddEntry(g, "All towers", "p"); }
         g = MakeGraph(en, &fe[0], kHot, true, kRed, 20, 0.6);                if (g) { g->Draw("P same"); l2->AddEntry(g, "Hot towers", "p"); }
         g = MakeGraph(en, &fe[0], kSig, true, kGreen + 1, 20, 0.6);          if (g) { g->Draw("P same"); l2->AddEntry(g, "3sigma towers", "p"); }
         g = MakeGraph(en, &fe[0], kSigE, true, kMagenta, 20, 0.6);           if (g) { g->Draw("P same"); l2->AddEntry(g, "E 3sigma towers", "p"); }
         HLine(0, NT + 1, res.muE, kRed, 3);
         HLine(0, NT + 1, res.muE + nSigma * res.sigmaE, kMagenta, 3);
         HLine(0, NT + 1, res.muE - nSigma * res.sigmaE, kMagenta, 3);
         l2->Draw();
         c->SaveAs(Form("%s_example.%s", outTag, ext));
      }
   }

   cout << "Wrote " << outTag << ".root, " << outTag << "_counts.pdf, " << outTag << "_cuts.pdf, " << outTag << "_maps.pdf, " << outTag << "_example." << ext << " and "
        << outTag << "_{runs_summary.txt,bad_runs.list,good_runs.list,lowstat_runs.list,bad_towers.list,bad_towers_perRun.txt}" << endl;
}