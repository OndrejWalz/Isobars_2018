// TowerQA2.C -- bad tower / bad run QA, built step by step
//   root -l -b -q 'TowerQA2.C+("merged.root")'
//
// Input: MapTower_RunID_NHits_TowerID (x = run index in bin r+1, y = tower ID = bin number 1..4800)
//        MapTower_RunID_NEvents       (run index r in bin r+1)
//
// Before the QA: eta-phi distribution of the towers (TowerMap.h) -> <outTag>_etaphi.pdf.
// Edge towers (|eta| = 0.967, 240 towers): their <nHits> is scaled by 1/edgeFactor BEFORE the dead / hot cut,
// so the dead / hot cut, the <nHits> distribution, the 3sigma cut and the example plot all use the scaled values.
// <E_T> = sum(E_T) / NHits is computed from the scaled sum(E_T) and the scaled NHits: the factor cancels, so <E_T> is unchanged.
//
// Per-run procedure (one pass over the runs with >= minEvents events):
//   step 1: <nHits> = NHits / NEvents for every tower
//           tower is DEAD if <nHits> < deadCut (0.02), HOT if <nHits> > hotCut (0.40); both are excluded from the next steps
//   step 1b: towers listed in excludeTowers (e.g. "3407") are excluded by hand right after the dead / hot cut:
//           in every run in which they are not dead / hot they are marked "excluded by hand" and skip the next steps;
//           they do not count for the bad-run cuts of step 5 and are final bad towers in step 6
//   step 2: <nHits> distribution of the remaining good towers, all analysed runs together;
//           a tower with <nHits> outside the 0.135 % / 99.865 % quantiles of that distribution is bad:
//           "<nHits> 3sigma equivalent" tower (isSig); it is excluded from the next steps
//   step 3: <E> = sum(E) / NHits of the remaining good towers, distribution over all analysed runs
//   step 4: a tower with <E> outside the 0.135 % / 99.865 % quantiles of that distribution is bad:
//           "<E> 3sigma equivalent" tower (isSigE); it is excluded from the next steps
//   steps 2 and 4: cutMode chooses the cut window:  0 = the 0.135 % / 99.865 % quantiles,
//           1 = the fixed values nHitsCutLow / nHitsCutHigh and eCutLow / eCutHigh,
//           2 = histogram mean +- nSigmaHist (5) x histogram RMS (from the histogram itself)
//           (the quantiles are always computed and printed; they are drawn for comparison except with cutMode 2)
//   step 4 only: with eLandau = true the <E_T> window is taken from a Landau fit of the step 3 distribution instead:
//           [MPV - nSigmaLandau x sigma, MPV + nSigmaLandau x sigma] (MPV = Landau location, sigma = Landau width); overrides cutMode
//   step 4 only: with eLandauRings = true one Landau is fitted per |eta| ring (20 rings, finely binned <E_T> of each ring),
//           and every tower is cut with the window of its own ring: [MPV_ring - nSigmaLandau x sigma_ring, MPV_ring + ...];
//           overrides eLandau and cutMode for step 4 (a ring whose fit fails uses the global window) -> <outTag>_meanE_landau_rings.pdf
//   step 4 only: with eQuantileRings = true (and eLandauRings = false) every |eta| ring gets its own window from the
//           0.135 % / 99.865 % quantiles of that ring's <E_T> distribution -> <outTag>_meanE_quantile_rings_grid.pdf
//   step 5: a run is bad if it has more than deadRunCut (1000) dead, more than hotRunCut (2) hot, more than sigRunCut (200)
//           <nHits> 3sigma or more than sigERunCut (200) <E> 3sigma towers; bad runs are removed from the next steps
//   step 6: clean the status map (kept runs only): for every tower, count the kept runs in which it was bad
//           (dead / hot / <nHits> 3sigma / <E_T> 3sigma; histogram hTowerNBad), then
//           tower bad in MORE than towerCleanCut (10) kept runs -> bad tower: marked bad in every run
//           tower bad in 1..towerCleanCut kept runs              -> tower stays good; the runs in which it was bad are excluded
//           -> <outTag>_final_badtowers.txt, <outTag>_final_badruns.txt, hTowerStatusFinal, <outTag>_status_final.pdf
//   (next steps are added after step 6)

#include <vector>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TParameter.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLine.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TNamed.h"
#include "TowerMap.h"                 // TowerMap::Eta(softID), TowerMap::Phi(softID); keep TowerMap.h next to this macro
using namespace std;

void TowerQA2(const char* inFile = "Tower_Test_Isobars_full.root",
              const char* runListFile = "/gpfs/mnt/gpfs01/star/pwg/walz/Isobars2018/data/Isobars_2018.txt",
              const char* outTag = "TowerQA2",
              int    minEvents = 0,    // runs with fewer good events are not analysed
              double deadCut = 0.02,     // step 1: tower is dead if NHits/NEvents < deadCut
              double hotCut  = 0.40,     // step 1: tower is hot  if NHits/NEvents > hotCut
              int    nBins = 150,        // <nHits> distribution: bins
              double xMax  = 0.2,       // <nHits> distribution: upper edge
              int    nBinsE = 200,       // <E> distribution: bins
              double xMaxE  = 5.5,       // <E> distribution: upper edge [GeV]
              int    exampleRun = 1400,     // example plot: run number (e.g. 19098021) or run index (e.g. 500); 0 = first analysed run
              int    deadRunCut = 480,  // step 5: run is bad if it has MORE than this many dead towers
              int    hotRunCut  = 2,     //         ... or more than this many hot towers
              int    sigRunCut  = 40,   //         ... or more than this many <nHits> 3sigma towers
              int    sigERunCut = 25,    //         ... or more than this many <E> 3sigma towers
              double edgeFactor = 0.67,  // edge towers (|eta| = 0.967): <nHits> and <E> are divided by this (= scaled by 1/0.67); 1 = no scaling
              int    cutMode = 2,        // steps 2 and 4: 0 = quantiles, 1 = the fixed values below, 2 = hist mean +- nSigmaHist x hist RMS
              double nHitsCutLow  = 0.02,   // step 2 (cutMode 1): tower is <nHits> 3sigma if scaled <nHits> < nHitsCutLow ...
              double nHitsCutHigh = 0.12,   //                      ... or > nHitsCutHigh
              double eCutLow  = 0.2,        // step 4 (cutMode 1): tower is <E> 3sigma if <E> < eCutLow [GeV] ...
              double eCutHigh = 0.8,        //                      ... or > eCutHigh [GeV]
              double nSigmaHist = 5.0,       // steps 2 and 4 (cutMode 2): width of the window in units of the histogram RMS
              bool   eLandau = false,        // step 4: true = <E_T> window from a Landau fit (overrides cutMode for step 4 only)
              double nSigmaLandau = 7.0,     // step 4 (eLandau): window = Landau MPV +- nSigmaLandau x Landau sigma
              int    nBinsFine = 1000,       // Landau fit: finely binned <E_T> histogram, number of bins ...
              double xMaxFine  = 5.5,        //             ... over 0 .. xMaxFine [GeV] (default 1 MeV per bin)
              double landauFitLow  = 0.3,    // Landau fit range [GeV] ...
              double landauFitHigh = 0.8,
              int    towerCleanCut = 5,    // step 6: tower bad in MORE than this many kept runs -> bad in all runs; otherwise its bad runs are excluded
              bool   eLandauRings = false,  // step 4: true = one Landau fit per |eta| ring, each tower cut with its ring's MPV +- nSigmaLandau x sigma
              bool   eQuantileRings = false, // step 4: true = 0.135 % / 99.865 % quantiles per |eta| ring, each tower cut with its ring's window (eLandauRings wins if both)
              const char* excludeTowers = "")  // step 1b: tower IDs excluded by hand after the dead / hot cut (comma or space separated, "" = none)
{
   gStyle->SetOptStat(0);
   // ---------- input ----------
   TFile* f = TFile::Open(inFile);
   if (!f || f->IsZombie()) { cout << "Cannot open " << inFile << endl; return; }
   TH2D* hN  = (TH2D*) f->Get("MapTower_RunID_NHits_TowerID");
   TH2D* hE  = (TH2D*) f->Get("MapTower_RunID_Energy_TowerID");   // sum of the hit energies, same axes as hN
   TH1D* hEv = (TH1D*) f->Get("MapTower_RunID_NEvents");
   if (!hN || !hE || !hEv) { cout << "Missing MapTower_* histograms in " << inFile << endl; return; }

   vector<unsigned int> runId(1, 0);                          // run index = line number, starting at 1
   {
      ifstream in(runListFile);
      if (!in.is_open()) cout << "WARNING: cannot open run list " << runListFile << " - run numbers will be 0." << endl;
      unsigned int runNum;
      while (in >> runNum) runId.push_back(runNum);
   }

   const int NT = hN->GetNbinsY();                            // 4800
   const int nRuns = hN->GetNbinsX() - 1;
   cout << "Input: " << nRuns << " runs, " << NT << " towers" << endl;
   double nEvTot = 0.;                                        // events of all runs, before any QA step
   for (int r = 1; r <= nRuns; ++r) nEvTot += hEv->GetBinContent(r + 1);
   cout << "Total number of events before QA: " << Form("%.0f", nEvTot) << endl;

   // ---------- eta-phi distribution of the towers (TowerMap.h) ----------
   // edge towers: |eta| = 0.967 (the last ring of each half of the BEMC). Their <nHits> is multiplied by sc[t] = 1/edgeFactor
   // already before the dead / hot cut; sc[t] = 1 for all other towers. <E_T> uses scaled sum(E_T) / scaled NHits (factor cancels).
   vector<double> sc(NT + 1, 1.0);
   const double scEdge = (edgeFactor > 0.) ? 1. / edgeFactor : 1.;
   int nEdge = 0;
   TH2D* hEtaPhi = new TH2D("hTowerEtaPhi", "BEMC towers (value 2 = edge towers, |#eta| = 0.967);#eta;#phi [rad]", 80, -1.0, 1.0, 120, -3.14159265, 3.14159265);
   for (int t = 1; t <= NT; ++t) {
      double teta = 0., tphi = 0.;
      if (!TowerMap::Get(t, teta, tphi)) continue;
      const bool edge = fabs(fabs(teta) - 0.967) < 0.0005;
      if (edge) { sc[t] = scEdge; ++nEdge; }
      hEtaPhi->Fill(teta, tphi, edge ? 2. : 1.);
   }
   cout << "Edge towers (|eta| = 0.967): " << nEdge << " (expected 240), scale factor " << scEdge << " for <nHits>, sum(E_T) and NHits (cancels in <E_T>)" << endl;
   {
      TCanvas* cep = new TCanvas("cEtaPhi", "eta-phi distribution of the towers", 1000, 650);
      cep->SetRightMargin(0.12);
      hEtaPhi->SetStats(0);
      hEtaPhi->SetMinimum(0.);
      hEtaPhi->SetMaximum(2.);
      hEtaPhi->Draw("COLZ");
      cep->SaveAs(Form("%s_etaphi.pdf", outTag));
      cout << "Wrote " << outTag << "_etaphi.pdf" << endl;
   }

   // ---------- eta rings (20 rings in |eta|; towers 1..20 are one phi row) for the check of the edge scaling ----------
   const int NRING = 20;
   double ringEta[20];
   for (int i = 0; i < NRING; ++i) ringEta[i] = fabs(TowerMap::Eta(i + 1));
   vector<int> ringOf(NT + 1, -1);
   for (int t = 1; t <= NT; ++t) {
      const double ae = fabs(TowerMap::Eta(t));
      for (int i = 0; i < NRING; ++i) if (fabs(ae - ringEta[i]) < 0.0005) ringOf[t] = i;
   }
   vector<double> rnU(NRING, 0.), rnS(NRING, 0.), rnC(NRING, 0.);     // <nHits>: sum unscaled, sum scaled, count
   vector<double> reU(NRING, 0.), reS(NRING, 0.), reC(NRING, 0.);     // <E>: the same
   TCanvas* cring = 0;

   // ---------- per-run results ----------
   vector<char> analysed(nRuns + 1, 0);
   vector<vector<char> > isDead(nRuns + 1);                   // isDead[r][t], t = 1..NT
   vector<vector<char> > isHot(nRuns + 1);
   vector<vector<char> > isSig(nRuns + 1);                    // <nHits> 3sigma equivalent (outside the quantiles)
   vector<vector<char> > alive(nRuns + 1);                    // alive[r][t] = not dead, not hot, not isSig -> input of the next steps
   vector<int> nDead(nRuns + 1, 0);
   vector<int> nHot(nRuns + 1, 0);
   vector<int> nSig(nRuns + 1, 0);
   vector<vector<char> > isSigE(nRuns + 1);                   // <E> 3sigma equivalent (outside the <E> quantiles)
   vector<int> nSigE(nRuns + 1, 0);
   vector<vector<char> > isMan(nRuns + 1);                    // excluded by hand (step 1b)

   // towers excluded by hand (step 1b)
   vector<char> manEx(NT + 1, 0);
   {
      TString lst(excludeTowers);
      lst.ReplaceAll(",", " ");
      istringstream is(lst.Data());
      int tx;
      while (is >> tx) {
         if (tx >= 1 && tx <= NT) manEx[tx] = 1;
         else cout << "WARNING: tower " << tx << " in excludeTowers is outside 1.." << NT << " - ignored" << endl;
      }
   }

   TH1D* hAll = new TH1D("hNHits_good", "<nHits> of good towers, all analysed runs;<nHits> = NHits / NEvents;tower x run entries", nBins, 0., xMax);

   // ---------- step 1: dead and hot towers ----------
   int nAnalysed = 0;
   long nManTot = 0;                                          // tower-runs excluded by hand (step 1b)
   for (int r = 1; r <= nRuns; ++r) {
      const double ev = hEv->GetBinContent(r + 1);
      if (ev < minEvents) continue;
      analysed[r] = 1;
      ++nAnalysed;

      isDead[r] = vector<char>(NT + 1, 0);                    // (no vector::assign: CINT)
      isHot[r]  = vector<char>(NT + 1, 0);
      isSig[r]  = vector<char>(NT + 1, 0);
      isSigE[r] = vector<char>(NT + 1, 0);
      isMan[r]  = vector<char>(NT + 1, 0);
      alive[r]  = vector<char>(NT + 1, 0);
      for (int t = 1; t <= NT; ++t) {
         const double nh0 = (ev > 0) ? hN->GetBinContent(r + 1, t) / ev : 0.;  // unscaled <nHits>; empty run -> 0 (dead), avoids NaN in the sums
         const double nh  = nh0 * sc[t];                         // scaled <nHits> (edge towers x 1/edgeFactor), used from the dead / hot cut on
         if (nh < deadCut)      { isDead[r][t] = 1; ++nDead[r]; }
         else if (nh > hotCut)  { isHot[r][t]  = 1; ++nHot[r]; }
         else if (manEx[t])     { isMan[r][t]  = 1; ++nManTot; }      // step 1b: excluded by hand, skips all next steps
         else {
            alive[r][t] = 1;
            hAll->Fill(nh);
            if (ringOf[t] >= 0) {
               rnU[ringOf[t]] += nh0; rnS[ringOf[t]] += nh; rnC[ringOf[t]] += 1.;
               const double nhitsR = hN->GetBinContent(r + 1, t);           // <E_T> ring sums: the same tower-runs as <nHits> (before the <nHits> 3sigma cut)
               if (nhitsR > 0) {
                  const double mERu = hE->GetBinContent(r + 1, t) / nhitsR;                          // unscaled <E_T>
                  const double mERs = (hE->GetBinContent(r + 1, t) * sc[t]) / (nhitsR * sc[t]);      // scaled sum(E_T) / scaled NHits (factor cancels)
                  reU[ringOf[t]] += mERu; reS[ringOf[t]] += mERs; reC[ringOf[t]] += 1.;
               }
            }
         }
      }
   }
   cout << "Runs analysed: " << nAnalysed << " (skipped < " << minEvents << " events: " << nRuns - nAnalysed << ")" << endl;
   cout << "Step 1b: towers excluded by hand: \"" << excludeTowers << "\"  (" << nManTot << " tower-runs that were not dead / hot)" << endl;

   // ---------- <nHits> distribution of the good towers (all analysed runs) ----------
   const double pMean = hAll->GetMean(), pRMS = hAll->GetRMS(), pSkew = hAll->GetSkewness();
   double pq[3] = {0.00135, 0.5, 0.99865}, qv[3] = {0., 0., 0.};            // quantiles where a Gaussian has -3 sigma / median / +3 sigma
   hAll->GetQuantiles(3, qv, pq);
   cout << "Good towers: mean " << pMean << "  RMS " << pRMS << "  skewness " << pSkew << "  overflow " << hAll->GetBinContent(nBins + 1) << endl;
   cout << "Quantiles 0.135 % / 50 % / 99.865 %: " << qv[0] << " / " << qv[1] << " / " << qv[2]
        << "  (distance from median: -" << qv[1] - qv[0] << " / +" << qv[2] - qv[1] << ")" << endl;
   const TString modeName = (cutMode == 1) ? TString("manual") : (cutMode == 2) ? TString::Format("hist mean #pm %g RMS", nSigmaHist) : TString("quantiles");
   const char* modeTxt  = (cutMode == 1) ? "manual" : (cutMode == 2) ? "hist mean +- nSigmaHist x RMS" : "quantiles";
   const double nLo = (cutMode == 1) ? nHitsCutLow  : (cutMode == 2) ? pMean - nSigmaHist * pRMS : qv[0];   // cuts applied in step 2
   const double nHi = (cutMode == 1) ? nHitsCutHigh : (cutMode == 2) ? pMean + nSigmaHist * pRMS : qv[2];
   cout << "Step 2 cut on <nHits> (" << modeTxt << "): [" << nLo << ", " << nHi << "]" << endl;

   TCanvas* cd = new TCanvas("cNHitsGood", "<nHits> of good towers", 900, 650);
   cd->SetLogy(1);
   cd->SetGrid(1, 1);
   const double hmax = hAll->GetMaximum();
   const double y1 = 0.5, y2 = 5. * hmax;
   hAll->SetMinimum(y1);
   hAll->SetMaximum(y2);
   hAll->SetLineColor(kBlack);
   hAll->Draw("hist");
   TLine* l3 = new TLine(qv[0], y1, qv[0], 0.7 * y2);
   TLine* l4 = new TLine(qv[2], y1, qv[2], 0.7 * y2);
   l3->SetLineColor(kBlue); l4->SetLineColor(kBlue);
   TLine* ls[2] = {l3, l4};
   if (cutMode != 2) for (int i = 0; i < 2; ++i) { ls[i]->SetLineStyle(2); ls[i]->SetLineWidth(2); ls[i]->Draw(); }   // quantiles not drawn with cutMode 2
   TLegend* lg = new TLegend(0.55, 0.62, 0.89, 0.89);
   lg->SetHeader(Form("mean %.4f  RMS %.4f  skew %.2f", pMean, pRMS, pSkew));
   lg->AddEntry(hAll, "good towers, all runs", "l");
   if (cutMode != 2) lg->AddEntry(l3, "0.135 % / 99.865 % quantiles", "l");
   if (cutMode != 0) {                                                        // applied cuts (manual or hist mean +- n RMS) in green
      TLine* l5 = new TLine(nLo, y1, nLo, 0.7 * y2);
      TLine* l6 = new TLine(nHi, y1, nHi, 0.7 * y2);
      l5->SetLineColor(kGreen + 2);  l6->SetLineColor(kGreen + 2);
      l5->SetLineWidth(3);           l6->SetLineWidth(3);
      l5->Draw();  l6->Draw();
      lg->AddEntry(l5, Form("%s cuts [%.4f, %.4f]", modeName.Data(), nLo, nHi), "l");
   }
   lg->Draw();
   cd->SaveAs(Form("%s_nhits_good.pdf", outTag));
   cout << "Wrote " << outTag << "_nhits_good.pdf" << endl;

   // ---------- step 2: <nHits> 3sigma equivalent towers (outside the quantiles) ----------
   TH1D* hCut = new TH1D("hNHits_afterSig", "<nHits> of good towers after the 3sigma equivalent cut;<nHits> = NHits / NEvents;tower x run entries", nBins, 0., xMax);
   long nSigTot = 0, nAliveTot = 0;
   for (int r = 1; r <= nRuns; ++r) {
      if (!analysed[r]) continue;
      const double ev = hEv->GetBinContent(r + 1);
      for (int t = 1; t <= NT; ++t) {
         if (!alive[r][t]) continue;
         const double nh = hN->GetBinContent(r + 1, t) / ev * sc[t];      // scaled <nHits> (edge towers x 1/0.69)
         if (nh < nLo || nh > nHi) { isSig[r][t] = 1; alive[r][t] = 0; ++nSig[r]; ++nSigTot; }
         else { ++nAliveTot; hCut->Fill(nh); }
      }
   }
   cout << "Step 2: <nHits> 3sigma equivalent towers: " << nSigTot << " tower-runs outside [" << nLo << ", " << nHi
        << "];  good towers left: " << nAliveTot << " tower-runs" << endl;

   TCanvas* cc = new TCanvas("cNHitsAfter", "<nHits> after the 3sigma equivalent cut", 900, 650);
   cc->SetLogy(1);
   cc->SetGrid(1, 1);
   hCut->SetMinimum(0.5);
   hCut->SetMaximum(5. * hmax);
   hCut->SetLineColor(kBlack);
   hCut->Draw("hist");
   TLegend* lc = new TLegend(0.55, 0.75, 0.89, 0.89);
   lc->SetHeader(Form("mean %.4f  RMS %.4f  skew %.2f", hCut->GetMean(), hCut->GetRMS(), hCut->GetSkewness()));
   lc->AddEntry(hCut, "good towers after the cut", "l");
   lc->Draw();
   cc->SaveAs(Form("%s_nhits_after3sigma.pdf", outTag));
   cout << "After the cut: mean " << hCut->GetMean() << "  RMS " << hCut->GetRMS() << "  skewness " << hCut->GetSkewness() << endl;
   cout << "Wrote " << outTag << "_nhits_after3sigma.pdf" << endl;

   // ---------- step 3: <E> = sum(E) / NHits of the remaining good towers, all analysed runs ----------
   TH1D* hEfine = new TH1D("hMeanE_fine", "<E_{T}> of good towers, fine binning;<E_{T}> = #SigmaE_{T} / NHits [GeV];tower x run entries", nBinsFine, 0., xMaxFine);   // for the Landau fit
   vector<TH1D*> hRingE(NRING, (TH1D*)0);                       // <E_T> per |eta| ring, fine binning (per-ring Landau fits)
   for (int i = 0; i < NRING; ++i)
      hRingE[i] = new TH1D(Form("hMeanE_ring%d", i), Form("<E_{T}> of good towers, |#eta| = %.3f;<E_{T}> = #SigmaE_{T} / NHits [GeV];tower x run entries", ringEta[i]),
                           nBinsFine, 0., xMaxFine);
   TH1D* hEall = new TH1D("hMeanE_good", "<E_{T}> of good towers, all analysed runs;<E_{T}> = #SigmaE_{T} / NHits [GeV];tower x run entries", nBinsE, 0., xMaxE);
   for (int r = 1; r <= nRuns; ++r) {
      if (!analysed[r]) continue;
      for (int t = 1; t <= NT; ++t) {
         if (!alive[r][t]) continue;
         const double nhits = hN->GetBinContent(r + 1, t);
         if (nhits <= 0) continue;
         const double mE0 = (hE->GetBinContent(r + 1, t) * sc[t]) / (nhits * sc[t]);   // <E_T> = scaled sum(E_T) / scaled NHits (edge factor cancels)
         hEfine->Fill(mE0);
         if (ringOf[t] >= 0) hRingE[ringOf[t]]->Fill(mE0);
         hEall->Fill(mE0);
      }
   }
   const double eMean = hEall->GetMean(), eRMS = hEall->GetRMS(), eSkew = hEall->GetSkewness();
   double eq[3] = {0.00135, 0.5, 0.99865}, eqv[3] = {0., 0., 0.};
   hEall->GetQuantiles(3, eqv, eq);
   cout << "Step 3: <E> of good towers: mean " << eMean << "  RMS " << eRMS << "  skewness " << eSkew
        << "  entries " << hEall->GetEntries() << "  overflow " << hEall->GetBinContent(nBinsE + 1) << endl;
   cout << "Quantiles 0.135 % / 50 % / 99.865 %: " << eqv[0] << " / " << eqv[1] << " / " << eqv[2]
        << "  (distance from median: -" << eqv[1] - eqv[0] << " / +" << eqv[2] - eqv[1] << ")" << endl;
   // Landau fit of <E_T>: always done and drawn (<outTag>_meanE_landau.pdf); used as the step 4 cut only if eLandau = true
   // fitted on the finely binned histogram hEfine, only in [landauFitLow, landauFitHigh]
   double lMPV = 0., lSig = 0.;
   TF1* lfit = new TF1("lfitE", "landau", landauFitLow, landauFitHigh);
   {
      const int bPk = hEfine->GetMaximumBin();                                // start values: peak position and height, 10 MeV width
      lfit->SetParameters(1e5, 0.6, 0.5);
   }
   lfit->SetNpx(10000);                                                       // draw the narrow peak with enough points
   if (hEfine->Fit(lfit, "QRL0") != 0) lfit = 0;                               // R = fit range, 0 = do not draw yet; non-zero status = failed
   if (lfit) {
      lMPV = lfit->GetParameter(1);  lSig = fabs(lfit->GetParameter(2));
      cout << "Landau fit of <E_T>: MPV " << lMPV << "  sigma " << lSig << "  chi2/ndf " << lfit->GetChisquare() << " / " << lfit->GetNDF() << endl;
   } else cout << "WARNING: Landau fit of <E_T> failed" << (eLandau ? " - using cutMode for step 4." : ".") << endl;
   const bool useLandau = eLandau && (lfit != 0);
   const double eLo = useLandau ? lMPV - nSigmaLandau * lSig : (cutMode == 1) ? eCutLow  : (cutMode == 2) ? eMean - nSigmaHist * eRMS : eqv[0];   // cuts applied in step 4
   const double eHi = useLandau ? lMPV + nSigmaLandau * lSig : (cutMode == 1) ? eCutHigh : (cutMode == 2) ? eMean + nSigmaHist * eRMS : eqv[2];
   const TString eModeName = useLandau ? TString::Format("Landau MPV #pm %g#sigma", nSigmaLandau) : modeName;
   cout << "Step 4 cut on <E> (" << (useLandau ? "Landau MPV +- nSigmaLandau x sigma" : modeTxt) << "): [" << eLo << ", " << eHi << "] GeV" << endl;

   TCanvas* ce = new TCanvas("cMeanE", "<E_{T}> of good towers", 900, 650);
   ce->SetLogy(1);
   ce->SetGrid(1, 1);
   const double eY1 = 0.5, eY2 = 5. * hEall->GetMaximum();
   hEall->SetMinimum(eY1);
   hEall->SetMaximum(eY2);
   hEall->SetLineColor(kBlack);
   hEall->Draw("hist");
   TLine* el3 = new TLine(eqv[0], eY1, eqv[0], 0.7 * eY2);
   TLine* el4 = new TLine(eqv[2], eY1, eqv[2], 0.7 * eY2);
   el3->SetLineColor(kBlue); el4->SetLineColor(kBlue);
   TLine* els[2] = {el3, el4};
   if (cutMode != 2) for (int i = 0; i < 2; ++i) { els[i]->SetLineStyle(2); els[i]->SetLineWidth(2); els[i]->Draw(); }   // quantiles not drawn with cutMode 2
   TLegend* leg = new TLegend(0.55, 0.62, 0.89, 0.89);
   leg->SetHeader(Form("mean %.4f  RMS %.4f  skew %.2f", eMean, eRMS, eSkew));
   leg->AddEntry(hEall, "good towers, all runs", "l");
   if (cutMode != 2) leg->AddEntry(el3, "0.135 % / 99.865 % quantiles", "l");
   if (eLandauRings) leg->AddEntry((TObject*)0, "cut: per-#eta-ring Landau windows", "");
   else if (eQuantileRings) leg->AddEntry((TObject*)0, "cut: per-#eta-ring quantile windows", "");
   if (!eLandauRings && !eQuantileRings && (cutMode != 0 || useLandau)) {                  // applied cuts (manual or hist mean +- n RMS) in green
      TLine* el5 = new TLine(eLo, eY1, eLo, 0.7 * eY2);
      TLine* el6 = new TLine(eHi, eY1, eHi, 0.7 * eY2);
      el5->SetLineColor(kGreen + 2);  el6->SetLineColor(kGreen + 2);
      el5->SetLineWidth(3);           el6->SetLineWidth(3);
      el5->Draw();  el6->Draw();
      leg->AddEntry(el5, Form("%s cuts [%.3f, %.3f] GeV", eModeName.Data(), eLo, eHi), "l");
   }
   leg->Draw();
   ce->SaveAs(Form("%s_meanE_good.pdf", outTag));
   cout << "Wrote " << outTag << "_meanE_good.pdf" << endl;

   // ---------- Landau fit of <E_T>: own plot, linear y axis so the shape of the peak is visible ----------
   TCanvas* cLan = 0;
   if (lfit) {
      cLan = new TCanvas("cMeanELandau", "Landau fit of <E_{T}>", 900, 650);
      cLan->SetGrid(1, 1);
      cLan->SetLogy(1);
      TH1D* hL = hEfine;
      hL->SetTitle(Form("Landau fit of <E_{T}> of good towers (%.1f MeV bins, fit range %.2f-%.2f GeV);<E_{T}> = #SigmaE_{T} / NHits [GeV];tower x run entries",
                        1000. * xMaxFine / nBinsFine, landauFitLow, landauFitHigh));
      const double lyMax = 5. * hEfine->GetMaximum(), ly0 = 0.5;
      hL->SetMinimum(ly0);
      hL->SetMaximum(lyMax);
      hL->SetLineColor(kBlack);
      hL->Draw("hist");
      lfit->SetLineColor(kRed);
      lfit->SetLineWidth(2);
      lfit->Draw("same");                                                     // drawn only over the fit range
      TLine* lm = new TLine(lMPV, ly0, lMPV, lyMax);                            // MPV (Landau location)
      lm->SetLineColor(kRed);  lm->SetLineStyle(2);  lm->SetLineWidth(2);  lm->Draw();
      TLine* lw1 = new TLine(lMPV - nSigmaLandau * lSig, ly0, lMPV - nSigmaLandau * lSig, lyMax);   // MPV +- nSigmaLandau x sigma
      TLine* lw2 = new TLine(lMPV + nSigmaLandau * lSig, ly0, lMPV + nSigmaLandau * lSig, lyMax);
      lw1->SetLineColor(kGreen + 2);  lw2->SetLineColor(kGreen + 2);
      lw1->SetLineStyle(useLandau ? 1 : 2);  lw2->SetLineStyle(useLandau ? 1 : 2);
      lw1->SetLineWidth(2);  lw2->SetLineWidth(2);
      lw1->Draw();  lw2->Draw();
      TLegend* ll = new TLegend(0.50, 0.62, 0.89, 0.89);
      ll->AddEntry(hL, "good towers, all runs", "l");
      ll->AddEntry(lfit, Form("Landau fit, #chi^{2}/ndf = %.0f / %d", lfit->GetChisquare(), lfit->GetNDF()), "l");
      ll->AddEntry(lm, Form("MPV = %.3f GeV, #sigma = %.3f GeV", lMPV, lSig), "l");
      ll->AddEntry(lw1, Form("MPV #pm %g#sigma%s", nSigmaLandau, useLandau ? " (applied cut)" : " (not applied)"), "l");
      ll->Draw();
      cLan->SaveAs(Form("%s_meanE_landau.pdf", outTag));
      cout << "Wrote " << outTag << "_meanE_landau.pdf" << endl;
   }

   // ---------- per-ring Landau fits of <E_T> (eLandauRings): one Landau per |eta| ring, own cut window per ring ----------
   vector<double> rMPV(NRING, 0.), rSig(NRING, 0.);
   vector<double> wLo(NT + 1, eLo), wHi(NT + 1, eHi);            // step 4 window of every tower (global unless eLandauRings)
   TH1D* hRingMPV = new TH1D("hLandauRingMPV",   "per-ring Landau MPV of <E_{T}>;ring (1 = smallest |#eta|);MPV [GeV]",   NRING, 0.5, NRING + 0.5);
   TH1D* hRingSig = new TH1D("hLandauRingSigma", "per-ring Landau #sigma of <E_{T}>;ring (1 = smallest |#eta|);#sigma [GeV]", NRING, 0.5, NRING + 0.5);
   TCanvas* cLanR = 0;
   TCanvas* cLanRG = 0;                                          // one pad per ring
   TCanvas* cMPV = 0;                                            // MPV +- sigma vs |eta|
   if (eLandauRings) {
      const int rcol[20] = {kRed, kOrange + 7, kOrange, kYellow + 2, kSpring - 1, kGreen + 2, kTeal - 1, kCyan + 1, kAzure + 7, kBlue,
                            kViolet, kMagenta, kPink + 1, kRed + 2, kOrange + 3, kGreen + 4, kCyan + 3, kBlue + 3, kMagenta + 3, kGray + 2};
      vector<char> rOK(NRING, 0);
      vector<TF1*> rfit(NRING, (TF1*)0);
      double yMaxR = 0., xHiR = 0.;
      cout << "Per-ring Landau fits of <E_T> (window = MPV +- " << nSigmaLandau << " sigma, fit range " << landauFitLow << "-" << landauFitHigh << " GeV):" << endl;
      for (int i = 0; i < NRING; ++i) {
         TH1D* h = hRingE[i];
         if (h->GetEntries() < 10) { cout << Form("  ring %2d |eta| = %.3f: too few entries - global window used", i + 1, ringEta[i]) << endl; continue; }
         TF1* fr = new TF1(Form("lfitE_ring%d", i), "landau", landauFitLow, landauFitHigh);
         const int bPk = h->GetMaximumBin();                                  // start values from this ring's peak, 10 MeV width
         fr->SetParameters(h->GetBinContent(bPk) / 0.18, h->GetBinCenter(bPk), 0.01);
         fr->SetNpx(10000);
         if (h->Fit(fr, "QRL0") != 0) { cout << Form("  ring %2d |eta| = %.3f: fit failed - global window used", i + 1, ringEta[i]) << endl; continue; }
         rOK[i] = 1;  rfit[i] = fr;
         rMPV[i] = fr->GetParameter(1);  rSig[i] = fabs(fr->GetParameter(2));
         hRingMPV->SetBinContent(i + 1, rMPV[i]);  hRingSig->SetBinContent(i + 1, rSig[i]);
         cout << Form("  ring %2d |eta| = %.3f: MPV %.4f  sigma %.4f  chi2/ndf %.0f / %d  window [%.3f, %.3f] GeV", i + 1, ringEta[i], rMPV[i], rSig[i],
                      fr->GetChisquare(), fr->GetNDF(), rMPV[i] - nSigmaLandau * rSig[i], rMPV[i] + nSigmaLandau * rSig[i]) << endl;
         if (h->GetMaximum() > yMaxR) yMaxR = h->GetMaximum();
         if (rMPV[i] + nSigmaLandau * rSig[i] > xHiR) xHiR = rMPV[i] + nSigmaLandau * rSig[i];
      }
      for (int t = 1; t <= NT; ++t) {
         const int rg = ringOf[t];
         if (rg >= 0 && rOK[rg]) { wLo[t] = rMPV[rg] - nSigmaLandau * rSig[rg];  wHi[t] = rMPV[rg] + nSigmaLandau * rSig[rg]; }
      }
      // all ring distributions with their fits on one plot, log y axis
      cLanR = new TCanvas("cMeanELandauRings", "per-ring Landau fits of <E_{T}>", 1000, 700);
      cLanR->SetGrid(1, 1);
      cLanR->SetLogy();  
     
      const double xr = (xHiR > 0.) ? ((2. * xHiR < xMaxFine) ? 2. * xHiR : xMaxFine) : xMaxFine;
      TH1F* frR = gPad->DrawFrame(0., 0.5, xr, 20. * (yMaxR > 0. ? yMaxR : 1.));   // log y: lower edge > 0, extra room for the legend
      frR->SetTitle(Form("<E_{T}> per |#eta| ring with Landau fits (%.1f MeV bins);<E_{T}> = #SigmaE_{T} / NHits [GeV];tower x run entries", 1000. * xMaxFine / nBinsFine));
      TLegend* lR = new TLegend(0.55, 0.40, 0.89, 0.89);
      lR->SetNColumns(2);
      lR->SetTextSize(0.022);
      for (int i = 0; i < NRING; ++i) {
         if (!rOK[i]) continue;
         hRingE[i]->SetLineColor(rcol[i]);
         hRingE[i]->Draw("hist same");
         rfit[i]->SetLineColor(rcol[i]);  rfit[i]->SetLineStyle(2);  rfit[i]->SetLineWidth(2);
         rfit[i]->Draw("same");
         lR->AddEntry(hRingE[i], Form("|#eta| %.3f: MPV %.3f", ringEta[i], rMPV[i]), "l");
      }
      
      lR->Draw();
      cLanR->SaveAs(Form("%s_meanE_landau_rings.pdf", outTag));

      cout << "Wrote " << outTag << "_meanE_landau_rings.pdf" << endl;

      // the same, one ring per pad (4 x 5), each with its fit and its cut window
      cLanRG = new TCanvas("cMeanELandauRingsGrid", "per-ring Landau fits of <E_{T}>, one pad per ring", 1600, 1200);
      cLanRG->Divide(5, 4, 0.002, 0.002);
      for (int i = 0; i < NRING; ++i) {
         cLanRG->cd(i + 1);
         gPad->SetGrid(1, 1);
         gPad->SetLeftMargin(0.15);  gPad->SetRightMargin(0.03);
         TH1D* hp = (TH1D*) hRingE[i]->Clone(Form("hMeanE_ring%d_pad", i));   // clone: the zoom must not change the overlay plot
         hp->SetStats(0);
         hp->SetLineColor(kBlack);
         if (rOK[i]) {
            const double lo = rMPV[i] - nSigmaLandau * rSig[i], hi = rMPV[i] + nSigmaLandau * rSig[i];
            const double x0 = (lo - 0.1 > 0.) ? lo - 0.1 : 0., x1 = hi + 0.2;
            hp->GetXaxis()->SetRangeUser(x0, x1 < xMaxFine ? x1 : xMaxFine);
            hp->SetMaximum(1.15 * hp->GetMaximum());
            hp->SetTitle(Form("|#eta| = %.3f: MPV %.3f, #sigma %.4f, #chi^{2}/ndf %.0f/%d;<E_{T}> [GeV];entries",
                              ringEta[i], rMPV[i], rSig[i], rfit[i]->GetChisquare(), rfit[i]->GetNDF()));
            hp->Draw("hist");
            rfit[i]->SetLineColor(kRed);  rfit[i]->SetLineStyle(1);
            rfit[i]->Draw("same");
            const double ym = hp->GetMaximum();
            TLine* w1 = new TLine(lo, 0., lo, ym);
            TLine* w2 = new TLine(hi, 0., hi, ym);
            w1->SetLineColor(kGreen + 2);  w2->SetLineColor(kGreen + 2);
            w1->SetLineWidth(2);           w2->SetLineWidth(2);
            w1->Draw();  w2->Draw();
         } else {
            hp->SetTitle(Form("|#eta| = %.3f: no fit (global window);<E_{T}> [GeV];entries", ringEta[i]));
            hp->Draw("hist");
         }
      }
      cLanRG->SaveAs(Form("%s_meanE_landau_rings_grid.pdf", outTag));
      cout << "Wrote " << outTag << "_meanE_landau_rings_grid.pdf" << endl;

      // MPV +- sigma of every ring vs |eta|
      vector<double> gxE, gyE, gex, gey;
      for (int i = 0; i < NRING; ++i) if (rOK[i]) { gxE.push_back(ringEta[i]);  gyE.push_back(rMPV[i]);  gex.push_back(0.);  gey.push_back(rSig[i]); }
      if (!gxE.empty()) {
         cMPV = new TCanvas("cLandauMPVvsEta", "Landau MPV of <E_{T}> vs |eta|", 900, 650);
         cMPV->SetGrid(1, 1);
         TGraphErrors* gm = new TGraphErrors((int)gxE.size(), &gxE[0], &gyE[0], &gex[0], &gey[0]);
         gm->SetTitle("Landau MPV of <E_{T}> per |#eta| ring (error bar = Landau #sigma);|#eta|;MPV [GeV]");
         gm->SetMarkerStyle(20);  gm->SetMarkerColor(kBlue);  gm->SetLineColor(kBlue);
         gm->Draw("AP");
         cMPV->SaveAs(Form("%s_landau_mpv_eta.pdf", outTag));
         cout << "Wrote " << outTag << "_landau_mpv_eta.pdf" << endl;
      }
   }

   // ---------- per-ring quantiles of <E_T> (eQuantileRings): each ring gets its own 0.135 % / 99.865 % window ----------
   const bool useQRings = eQuantileRings && !eLandauRings;
   if (eQuantileRings && eLandauRings) cout << "WARNING: eQuantileRings and eLandauRings both set - the per-ring Landau windows are used." << endl;
   TCanvas* cQR = 0;                                             // one pad per ring
   TCanvas* cQRe = 0;                                            // window and median vs |eta|
   if (useQRings) {
      vector<double> qLo(NRING, eLo), qMed(NRING, 0.), qHi(NRING, eHi);
      vector<char> qOK(NRING, 0);
      cout << "Per-ring quantiles of <E_T> (0.135 % / 50 % / 99.865 %, " << 1000. * xMaxFine / nBinsFine << " MeV bins):" << endl;
      for (int i = 0; i < NRING; ++i) {
         if (hRingE[i]->GetEntries() < 10) { cout << Form("  ring %2d |eta| = %.3f: too few entries - global window used", i + 1, ringEta[i]) << endl; continue; }
         double pr[3] = {0.00135, 0.5, 0.99865}, qr[3] = {0., 0., 0.};
         hRingE[i]->GetQuantiles(3, qr, pr);
         qLo[i] = qr[0];  qMed[i] = qr[1];  qHi[i] = qr[2];  qOK[i] = 1;
         cout << Form("  ring %2d |eta| = %.3f: window [%.3f, %.3f] GeV, median %.3f GeV, entries %.0f", i + 1, ringEta[i], qLo[i], qHi[i], qMed[i], hRingE[i]->GetEntries()) << endl;
      }
      for (int t = 1; t <= NT; ++t) {
         const int rg = ringOf[t];
         if (rg >= 0 && qOK[rg]) { wLo[t] = qLo[rg];  wHi[t] = qHi[rg]; }
      }
      // one ring per pad (4 x 5), log y, with the ring's window
      cQR = new TCanvas("cMeanEQuantileRingsGrid", "per-ring quantile windows of <E_{T}>, one pad per ring", 1600, 1200);
      cQR->Divide(5, 4, 0.002, 0.002);
      for (int i = 0; i < NRING; ++i) {
         cQR->cd(i + 1);
         gPad->SetGrid(1, 1);  gPad->SetLogy(1);
         gPad->SetLeftMargin(0.15);  gPad->SetRightMargin(0.03);
         TH1D* hp = (TH1D*) hRingE[i]->Clone(Form("hMeanE_ring%d_qpad", i));   // clone: own axis range per pad
         hp->SetStats(0);
         hp->SetLineColor(kBlack);
         hp->SetMinimum(0.5);
         if (qOK[i]) {
            hp->GetXaxis()->SetRangeUser(0., 5.5);                              // fixed range 0 - 5.5 GeV for all rings
            hp->SetMaximum(5. * hp->GetMaximum());
            hp->SetTitle(Form("|#eta| = %.3f: window [%.3f, %.3f];<E_{T}> [GeV];entries", ringEta[i], qLo[i], qHi[i]));
            hp->Draw("hist");
            const double ym = hp->GetMaximum();
            TLine* w1 = new TLine(qLo[i], 0.5, qLo[i], ym);
            TLine* w2 = new TLine(qHi[i], 0.5, qHi[i], ym);
            w1->SetLineColor(kGreen + 2);  w2->SetLineColor(kGreen + 2);
            w1->SetLineWidth(2);           w2->SetLineWidth(2);
            w1->Draw();  w2->Draw();
         } else {
            hp->SetTitle(Form("|#eta| = %.3f: too few entries (global window);<E_{T}> [GeV];entries", ringEta[i]));
            hp->Draw("hist");
         }
      }
      cQR->SaveAs(Form("%s_meanE_quantile_rings_grid.pdf", outTag));
      cout << "Wrote " << outTag << "_meanE_quantile_rings_grid.pdf" << endl;

      // lower edge, median and upper edge of every ring vs |eta|
      vector<double> gx, gl, gm, gh;
      for (int i = 0; i < NRING; ++i) if (qOK[i]) { gx.push_back(ringEta[i]);  gl.push_back(qLo[i]);  gm.push_back(qMed[i]);  gh.push_back(qHi[i]); }
      if (!gx.empty()) {
         cQRe = new TCanvas("cQuantileRingsVsEta", "per-ring quantile window of <E_{T}> vs |eta|", 900, 650);
         cQRe->SetGrid(1, 1);
         double ymn = gl[0], ymx = gh[0];
         for (size_t k = 0; k < gx.size(); ++k) { if (gl[k] < ymn) ymn = gl[k];  if (gh[k] > ymx) ymx = gh[k]; }
         TH1F* fq = gPad->DrawFrame(0., ymn - 0.1 * (ymx - ymn), 1.0, ymx + 0.2 * (ymx - ymn));
         fq->SetTitle("per-ring window of <E_{T}> (0.135 % / 99.865 % quantiles);|#eta|;<E_{T}> [GeV]");
         TGraph* g1 = new TGraph((int)gx.size(), &gx[0], &gl[0]);
         TGraph* g2 = new TGraph((int)gx.size(), &gx[0], &gm[0]);
         TGraph* g3 = new TGraph((int)gx.size(), &gx[0], &gh[0]);
         g1->SetMarkerStyle(22);  g1->SetMarkerColor(kGreen + 2);  g1->SetLineColor(kGreen + 2);
         g2->SetMarkerStyle(20);  g2->SetMarkerColor(kBlue);       g2->SetLineColor(kBlue);
         g3->SetMarkerStyle(23);  g3->SetMarkerColor(kGreen + 2);  g3->SetLineColor(kGreen + 2);
         g1->Draw("PL same");  g2->Draw("PL same");  g3->Draw("PL same");
         TLegend* lq = new TLegend(0.14, 0.72, 0.50, 0.89);
         lq->AddEntry(g3, "upper edge (99.865 %)", "pl");
         lq->AddEntry(g2, "median", "pl");
         lq->AddEntry(g1, "lower edge (0.135 %)", "pl");
         lq->Draw();
         cQRe->SaveAs(Form("%s_quantile_rings_eta.pdf", outTag));
         cout << "Wrote " << outTag << "_quantile_rings_eta.pdf" << endl;
      }
   }

   // ---------- check of the edge scaling: mean <nHits> and <E> per eta ring, before and after scaling ----------
   // <nHits> and <E_T>: tower-runs that are not dead / hot / excluded by hand (step 1), before the 3sigma cuts; all analysed runs
   {
      vector<double> px, nU, nS, eU, eS;
      for (int i = 0; i < NRING; ++i) {
         if (rnC[i] <= 0. || reC[i] <= 0.) continue;
         px.push_back(ringEta[i]);
         nU.push_back(rnU[i] / rnC[i]);  nS.push_back(rnS[i] / rnC[i]);
         eU.push_back(reU[i] / reC[i]);  eS.push_back(reS[i] / reC[i]);
      }
      const int np = (int)px.size();
      if (np < 2) cout << "Ring check skipped: fewer than 2 eta rings have good towers." << endl;
      else {
      cring = new TCanvas("cRingMeans", "mean <nHits> and <E_{T}> per eta ring", 900, 900);
      cring->Divide(1, 2);
      double nMaxR = 0., eMaxR = 0.;
      for (int i = 0; i < np; ++i) {
         if (nU[i] > nMaxR) nMaxR = nU[i];  if (nS[i] > nMaxR) nMaxR = nS[i];
         if (eU[i] > eMaxR) eMaxR = eU[i];  if (eS[i] > eMaxR) eMaxR = eS[i];
      }
      for (int p = 1; p <= 2; ++p) {
         cring->cd(p);
         gPad->SetGrid(1, 1);
         TH1F* frr = gPad->DrawFrame(0., 0., 1.0, 1.3 * (p == 1 ? nMaxR : eMaxR));
         frr->SetTitle(p == 1 ? "mean <nHits> per #eta ring;|#eta|;mean <nHits>" : "mean <E_{T}> per #eta ring;|#eta|;mean <E_{T}> [GeV]");
         TGraph* gu = new TGraph(np, &px[0], p == 1 ? &nU[0] : &eU[0]);
         TGraph* gs2 = new TGraph(np, &px[0], p == 1 ? &nS[0] : &eS[0]);
         gu->SetMarkerStyle(20);   gu->SetMarkerSize(1.0);   gu->SetMarkerColor(kBlack);  gu->SetLineColor(kBlack);
         gs2->SetMarkerStyle(24);  gs2->SetMarkerSize(1.3);  gs2->SetMarkerColor(kRed);   gs2->SetLineColor(kRed);
         gu->Draw("PL same");
         gs2->Draw("PL same");
         TLegend* lr = new TLegend(0.14, 0.72, 0.60, 0.89);
         lr->AddEntry(gu, "before scaling", "pl");
         lr->AddEntry(gs2, Form("after scaling (edge towers x %.3f)", scEdge), "pl");
         lr->Draw();
      }
      cring->SaveAs(Form("%s_ring_means.pdf", outTag));
      cout << "Edge ring (|eta| = " << px[np - 1] << "): <nHits> " << nU[np - 1] << " -> " << nS[np - 1] << " (next ring " << nU[np - 2]
           << "),  <E> " << eU[np - 1] << " -> " << eS[np - 1] << " (next ring " << eU[np - 2] << ")" << endl;
      cout << "Wrote " << outTag << "_ring_means.pdf" << endl;
      }
   }

   // ---------- step 4: <E> 3sigma equivalent towers (outside the <E> quantiles) ----------
   TH1D* hEcut = new TH1D("hMeanE_afterSig", "<E_{T}> of good towers after the 3sigma equivalent cut;<E_{T}> = #SigmaE_{T} / NHits [GeV];tower x run entries", nBinsE, 0., xMaxE);
   long nSigETot = 0, nAliveETot = 0;
   for (int r = 1; r <= nRuns; ++r) {
      if (!analysed[r]) continue;
      for (int t = 1; t <= NT; ++t) {
         if (!alive[r][t]) continue;
         const double nhits = hN->GetBinContent(r + 1, t);
         if (nhits <= 0) continue;
         const double meanE = (hE->GetBinContent(r + 1, t) * sc[t]) / (nhits * sc[t]);   // as in step 3 (edge factor cancels)
         if (meanE < wLo[t] || meanE > wHi[t]) { isSigE[r][t] = 1; alive[r][t] = 0; ++nSigE[r]; ++nSigETot; }
         else { ++nAliveETot; hEcut->Fill(meanE); }
      }
   }
   if (eLandauRings || useQRings) cout << "Step 4: <E> 3sigma equivalent towers: " << nSigETot << " tower-runs outside their ring's " << (eLandauRings ? "Landau" : "quantile") << " window";
   else              cout << "Step 4: <E> 3sigma equivalent towers: " << nSigETot << " tower-runs outside [" << eLo << ", " << eHi << "] GeV";
   cout << ";  good towers left: " << nAliveETot << " tower-runs" << endl;

   TCanvas* cf = new TCanvas("cMeanEAfter", "<E_{T}> after the 3sigma equivalent cut", 900, 650);
   cf->SetLogy(1);
   cf->SetGrid(1, 1);
   hEcut->SetMinimum(0.5);
   hEcut->SetMaximum(eY2);
   hEcut->SetLineColor(kBlack);
   hEcut->Draw("hist");
   TLegend* lf = new TLegend(0.55, 0.75, 0.89, 0.89);
   lf->SetHeader(Form("mean %.4f  RMS %.4f  skew %.2f", hEcut->GetMean(), hEcut->GetRMS(), hEcut->GetSkewness()));
   lf->AddEntry(hEcut, "good towers after the cuts", "l");
   lf->Draw();
   cf->SaveAs(Form("%s_meanE_after3sigma.pdf", outTag));
   cout << "After the cut: <E> mean " << hEcut->GetMean() << "  RMS " << hEcut->GetRMS() << "  skewness " << hEcut->GetSkewness() << endl;
   cout << "Wrote " << outTag << "_meanE_after3sigma.pdf" << endl;

   // ---------- step 5: bad runs ----------
   // a run is bad if it has more than deadRunCut dead, more than hotRunCut hot, more than sigRunCut <nHits> 3sigma
   // or more than sigERunCut <E> 3sigma towers; bad runs are removed from the next steps
   vector<char> isBadRun(nRuns + 1, 0);
   int nBadRuns = 0, nBadDead = 0, nBadHot = 0, nBadSig = 0, nBadSigE = 0;
   FILE* fbad = fopen(Form("%s_badruns.txt", outTag), "w");
   fprintf(fbad, "# runId runIndex nEvents nDead nHot nSigNHits nSigE   (reasons: D dead, H hot, N <nHits> 3sigma, E <E> 3sigma)\n");
   for (int r = 1; r <= nRuns; ++r) {
      if (!analysed[r]) continue;
      const bool bD = nDead[r] > deadRunCut, bH = nHot[r] > hotRunCut, bN = nSig[r] > sigRunCut, bE = nSigE[r] > sigERunCut;
      if (bD) ++nBadDead;
      if (bH) ++nBadHot;
      if (bN) ++nBadSig;
      if (bE) ++nBadSigE;
      if (!(bD || bH || bN || bE)) continue;
      isBadRun[r] = 1;
      ++nBadRuns;
      fprintf(fbad, "%u %d %.0f %d %d %d %d   %s%s%s%s\n", (r < (int)runId.size()) ? runId[r] : 0u, r, hEv->GetBinContent(r + 1),
              nDead[r], nHot[r], nSig[r], nSigE[r], bD ? "D" : "", bH ? "H" : "", bN ? "N" : "", bE ? "E" : "");
   }
   fclose(fbad);
   cout << "Step 5: bad runs: " << nBadRuns << " of " << nAnalysed << " analysed runs, good runs left: " << nAnalysed - nBadRuns << endl;
   cout << "  (runs failing each criterion, a run can fail several: dead > " << deadRunCut << ": " << nBadDead << ",  hot > " << hotRunCut << ": " << nBadHot
        << ",  <nHits> 3sigma > " << sigRunCut << ": " << nBadSig << ",  <E> 3sigma > " << sigERunCut << ": " << nBadSigE << ")" << endl;
   cout << "Wrote " << outTag << "_badruns.txt" << endl;

   // ---------- step 6: clean the status map (kept runs only) ----------
   // count for every tower the kept runs in which it is bad
   vector<int> cntDead(NT + 1, 0), cntHot(NT + 1, 0), cntSig(NT + 1, 0), cntSigE(NT + 1, 0), cntMan(NT + 1, 0);
   int nKept = 0;
   for (int r = 1; r <= nRuns; ++r) {
      if (!analysed[r] || isBadRun[r]) continue;
      ++nKept;
      for (int t = 1; t <= NT; ++t) {
         if (isDead[r][t])      ++cntDead[t];
         else if (isHot[r][t])  ++cntHot[t];
         else if (isSig[r][t])  ++cntSig[t];
         else if (isSigE[r][t]) ++cntSigE[t];
         else if (isMan[r][t])  ++cntMan[t];
      }
   }
   TH1D* hTowerBad = new TH1D("hTowerNBad", "number of kept runs in which the tower is bad;towerID;N runs", NT, 0.5, NT + 0.5);
   for (int t = 1; t <= NT; ++t) hTowerBad->SetBinContent(t, cntDead[t] + cntHot[t] + cntSig[t] + cntSigE[t] + cntMan[t]);

   // tower bad in MORE than towerCleanCut kept runs -> bad tower, marked bad in every run
   // tower bad in 1..towerCleanCut kept runs        -> tower stays good, the runs in which it was bad are excluded
   // ponytail: single pass on these counts; runs excluded here do not lower the counts of other towers (iterate if that matters)
   vector<char> towerBadFinal(NT + 1, 0);
   vector<char> runCleaned(nRuns + 1, 0);                     // run excluded by step 6
   int nTowFinal = 0, nRunClean = 0;
   for (int t = 1; t <= NT; ++t)
      if (manEx[t] || cntDead[t] + cntHot[t] + cntSig[t] + cntSigE[t] > towerCleanCut) { towerBadFinal[t] = 1; ++nTowFinal; }   // excluded by hand: always a final bad tower

   // before the clean-up: towers bad in EVERY kept run vs the towers the clean-up adds or keeps
   {
      int nAll = 0, nAdded = 0, nFew = 0;
      for (int t = 1; t <= NT; ++t) {
         const int nb = cntDead[t] + cntHot[t] + cntSig[t] + cntSigE[t] + cntMan[t];
         if (nb >= nKept && nKept > 0) ++nAll;                // bad in every kept run (already a bad tower before step 6)
         else if (towerBadFinal[t])    ++nAdded;              // bad in > towerCleanCut but not all kept runs -> made bad in every run by step 6
         else if (nb > 0)              ++nFew;                // bad in 1..towerCleanCut kept runs -> stays good, its bad runs are excluded
      }
      cout << "Step 6: before the clean-up: towers bad in all " << nKept << " kept runs: " << nAll << endl;
      cout << "        added by the clean-up (bad in > " << towerCleanCut << " but not all kept runs): " << nAdded
           << "  -> final bad towers: " << nAll + nAdded << endl;
      cout << "        towers bad in 1.." << towerCleanCut << " kept runs (stay good, their bad runs are excluded): " << nFew << endl;
   }

   FILE* ffR = fopen(Form("%s_final_badruns.txt", outTag), "w");
   fprintf(ffR, "# final bad runs: step 5 bad runs and runs excluded in step 6 (towerCleanCut = %d)\n", towerCleanCut);
   fprintf(ffR, "# runId runIndex nEvents step reason   (step 5: D dead, H hot, N <nHits> 3sigma, E <E_T> 3sigma;  step 6: towers bad in <= %d runs that were bad in this run)\n", towerCleanCut);
   double nEvFinal = 0.;                                      // events in the final good runs
   double nEvClean = 0.;                                      // events in the runs excluded by step 6
   for (int r = 1; r <= nRuns; ++r) {
      if (!analysed[r]) continue;
      const unsigned int rid = (r < (int)runId.size()) ? runId[r] : 0u;
      if (isBadRun[r]) {
         fprintf(ffR, "%u %d %.0f 5 %s%s%s%s\n", rid, r, hEv->GetBinContent(r + 1),
                 nDead[r] > deadRunCut ? "D" : "", nHot[r] > hotRunCut ? "H" : "", nSig[r] > sigRunCut ? "N" : "", nSigE[r] > sigERunCut ? "E" : "");
         continue;
      }
      TString why;                                            // towers that cause the exclusion of this run
      for (int t = 1; t <= NT; ++t) {
         if (towerBadFinal[t]) continue;
         if (isDead[r][t] || isHot[r][t] || isSig[r][t] || isSigE[r][t]) why += TString::Format("%s%d", why.Length() ? "," : "towers ", t);
      }
      if (why.Length()) {
         runCleaned[r] = 1;  ++nRunClean;  nEvClean += hEv->GetBinContent(r + 1);
         fprintf(ffR, "%u %d %.0f 6 %s\n", rid, r, hEv->GetBinContent(r + 1), why.Data());
      } else nEvFinal += hEv->GetBinContent(r + 1);
   }
   fclose(ffR);

   FILE* ffT = fopen(Form("%s_final_badtowers.txt", outTag), "w");
   fprintf(ffT, "# final bad towers: bad in more than %d of the %d kept runs (step 5), marked bad in every run\n", towerCleanCut, nKept);
   fprintf(ffT, "# + towers excluded by hand (excludeTowers = \"%s\"), bad in every run\n", excludeTowers);
   fprintf(ffT, "# towerID nDead nHot nSigNHits nSigE nExclHand nBadTotal\n");
   for (int t = 1; t <= NT; ++t)
      if (towerBadFinal[t]) fprintf(ffT, "%d %d %d %d %d %d %d\n", t, cntDead[t], cntHot[t], cntSig[t], cntSigE[t], cntMan[t],
                                    cntDead[t] + cntHot[t] + cntSig[t] + cntSigE[t] + cntMan[t]);
   fclose(ffT);

   // short summary: final bad towers, each counted once
   //   bad in ALL kept runs           -> under its main reason (dead / hot / <nHits> 3sigma / <E_T> 3sigma, the kind it was bad for in most runs)
   //   bad in > towerCleanCut but not all kept runs -> fifth category "badInMoreThanNRuns"
   {
      int nMain[6] = {0, 0, 0, 0, 0, 0};                      // dead, hot, <nHits> 3sigma, <E_T> 3sigma, bad in > towerCleanCut (not all) runs, excluded by hand
      const char* mainName[6] = {"dead", "hot", "nHits3sigma", "ET3sigma", "badInMoreThanNRuns", "excludedByHand"};
      TString ids;
      for (int t = 1; t <= NT; ++t) {
         if (!towerBadFinal[t]) continue;
         const int c[4] = {cntDead[t], cntHot[t], cntSig[t], cntSigE[t]};
         int m = 0;
         for (int k = 1; k < 4; ++k) if (c[k] > c[m]) m = k;
         if (c[0] + c[1] + c[2] + c[3] < nKept) m = 4;               // not bad in every kept run
         if (manEx[t]) m = 5;                                         // excluded by hand
         ++nMain[m];
         ids += TString::Format("%d %s\n", t, mainName[m]);
      }
      FILE* fs = fopen(Form("%s_summary.txt", outTag), "w");
      fprintf(fs, "# bad towers: bad in more than %d of the %d kept runs -> bad in every run\n", towerCleanCut, nKept);
      fprintf(fs, "# each tower is counted once: towers bad in all kept runs under the reason they were bad for in most runs,\n");
      fprintf(fs, "# towers bad in more than %d but not all kept runs under \"bad in > %d runs\"\n", towerCleanCut, towerCleanCut);
      fprintf(fs, "bad towers total (> %d runs): %d\n", towerCleanCut, nTowFinal);
      fprintf(fs, "  dead:                %d\n", nMain[0]);
      fprintf(fs, "  hot:                 %d\n", nMain[1]);
      fprintf(fs, "  <nHits> 3sigma:      %d\n", nMain[2]);
      fprintf(fs, "  <E_T> 3sigma:        %d\n", nMain[3]);
      fprintf(fs, "  bad in > %d runs:    %d\n", towerCleanCut, nMain[4]);
      fprintf(fs, "  excluded by hand:    %d   (excludeTowers = \"%s\")\n", nMain[5], excludeTowers);
      fprintf(fs, "# towerID mainReason\n%s", ids.Data());
      fclose(fs);
      cout << "Wrote " << outTag << "_summary.txt" << endl;
   }

   cout << "Step 6: final bad towers (bad in > " << towerCleanCut << " of " << nKept << " kept runs): " << nTowFinal << " of " << NT << endl;
   cout << "        runs excluded because of towers bad in <= " << towerCleanCut << " runs: " << nRunClean
        << " (" << Form("%.0f", nEvClean) << " events);  final good runs: " << nKept - nRunClean << " of " << nAnalysed << endl;
   cout << "Total number of events after QA (final good runs): " << Form("%.0f", nEvFinal) << " of " << Form("%.0f", nEvTot) << endl;
   cout << "Wrote " << outTag << "_final_badtowers.txt, " << outTag << "_final_badruns.txt" << endl;

   // ---------- run characteristics: number of bad towers of each kind per run ----------
   TH1D* hRun[4];
   hRun[0] = new TH1D("hRunDead", Form("nHits/nEvents < %.0f%%;runIndex;#DeadTowers", 100. * deadCut), nRuns, 0.5, nRuns + 0.5);
   hRun[1] = new TH1D("hRunHot",  Form("nHits/nEvents > %.0f%%;runIndex;#HotTowers", 100. * hotCut), nRuns, 0.5, nRuns + 0.5);
   hRun[2] = new TH1D("hRunSig",  "mean nhits sigma;runIndex;#sigma nHits", nRuns, 0.5, nRuns + 0.5);
   hRun[3] = new TH1D("hRunSigE", "mean energy sigma;runIndex;#sigma E_{T}", nRuns, 0.5, nRuns + 0.5);
   for (int r = 1; r <= nRuns; ++r) {
      if (!analysed[r]) continue;
      hRun[0]->SetBinContent(r, nDead[r]);
      hRun[1]->SetBinContent(r, nHot[r]);
      hRun[2]->SetBinContent(r, nSig[r]);
      hRun[3]->SetBinContent(r, nSigE[r]);
   }
   TCanvas* cr = new TCanvas("cRunChar", "run characteristics", 900, 1100);
   cr->Divide(1, 4);
   for (int i = 0; i < 4; ++i) {
      cr->cd(i + 1);
      gPad->SetGrid(1, 1);
      hRun[i]->SetStats(0);
      hRun[i]->SetFillColor(kBlue);
      hRun[i]->SetLineColor(kBlue);
      hRun[i]->Draw("hist");
      const double rc = (i == 0) ? deadRunCut : (i == 1) ? hotRunCut : (i == 2) ? sigRunCut : sigERunCut;   // step 5 cut of this panel
      TLine* lrc = new TLine(0.5, rc, nRuns + 0.5, rc);
      lrc->SetLineColor(kRed);  lrc->SetLineStyle(2);  lrc->SetLineWidth(2);
      lrc->Draw();
   }
   cr->SaveAs(Form("%s_runs.pdf", outTag));
   cout << "Wrote " << outTag << "_runs.pdf" << endl;

   // ---------- run characteristics of the kept (good) runs only ----------
   TH1D* hRunG[4];
   hRunG[0] = new TH1D("hRunDeadGood", Form("good runs: nHits/nEvents < %.0f%%;runIndex;#DeadTowers", 100. * deadCut), nRuns, 0.5, nRuns + 0.5);
   hRunG[1] = new TH1D("hRunHotGood",  Form("good runs: nHits/nEvents > %.0f%%;runIndex;#HotTowers", 100. * hotCut), nRuns, 0.5, nRuns + 0.5);
   hRunG[2] = new TH1D("hRunSigGood",  "mean nhits sigma;runIndex;#sigma nHits", nRuns, 0.5, nRuns + 0.5);
   hRunG[3] = new TH1D("hRunSigEGood", "mean energy sigma;runIndex;#sigma E_{T}", nRuns, 0.5, nRuns + 0.5);
   for (int r = 1; r <= nRuns; ++r) {
      if (!analysed[r] || isBadRun[r]) continue;
      hRunG[0]->SetBinContent(r, nDead[r]);
      hRunG[1]->SetBinContent(r, nHot[r]);
      hRunG[2]->SetBinContent(r, nSig[r]);
      hRunG[3]->SetBinContent(r, nSigE[r]);
   }
   TCanvas* crg = new TCanvas("cRunCharGood", "run characteristics of the good runs", 900, 1100);
   crg->Divide(1, 4);
   for (int i = 0; i < 4; ++i) {
      crg->cd(i + 1);
      gPad->SetGrid(1, 1);
      hRunG[i]->SetStats(0);
      hRunG[i]->SetFillColor(kBlue);
      hRunG[i]->SetLineColor(kBlue);
      hRunG[i]->Draw("hist");
      const double rcg = (i == 0) ? deadRunCut : (i == 1) ? hotRunCut : (i == 2) ? sigRunCut : sigERunCut;
      TLine* lrg = new TLine(0.5, rcg, nRuns + 0.5, rcg);
      lrg->SetLineColor(kRed);  lrg->SetLineStyle(2);  lrg->SetLineWidth(2);
      lrg->Draw();
   }
   crg->SaveAs(Form("%s_runs_good.pdf", outTag));
   cout << "Wrote " << outTag << "_runs_good.pdf" << endl;

   // ---------- example run: the same procedure drawn for one run (tower by tower) ----------
   // The cuts are the global ones found above (dead, hot, and the <nHits> and <E> quantiles of all runs).
   TCanvas* cx = 0;
   {
      int er = 0;                                                 // run index of the example run
      for (int r = 1; r <= nRuns && r < (int)runId.size(); ++r) if (exampleRun > 0 && runId[r] == (unsigned int)exampleRun) er = r;
      if (er == 0 && exampleRun > 0 && exampleRun <= nRuns) er = exampleRun;      // a small number is taken as run index
      if (er == 0 && exampleRun <= 0) { for (int r = 1; r <= nRuns && er == 0; ++r) if (analysed[r]) er = r; }
      if (er == 0 || !analysed[er]) {
         cout << "Example run " << exampleRun << " not found among the analysed runs - no example plot." << endl;
      } else {
         const double evx = hEv->GetBinContent(er + 1);
         const unsigned int rn = (er < (int)runId.size()) ? runId[er] : 0u;
         // classes for the <nHits> panel: 0 good, 1 dead, 2 hot, 3 <nHits> 3sigma equivalent
         vector<vector<double> > gx(5), gy(5);
         vector<double> ex, ey, exs, eys;                          // <E> panel: good, <E> 3sigma equivalent
         for (int t = 1; t <= NT; ++t) {
            const double nhits = hN->GetBinContent(er + 1, t);
            int cl = 0;
            if (isDead[er][t]) cl = 1; else if (isHot[er][t]) cl = 2; else if (isSig[er][t]) cl = 3; else if (isMan[er][t]) cl = 4;
            const double nh0 = (evx > 0) ? nhits / evx : 0.;                             // empty run -> 0, as in step 1
            const double nhx = nh0 * sc[t];                                              // scaled for all towers, as in the dead / hot cut
            gx[cl].push_back(t);  gy[cl].push_back(nhx);
            if (cl == 0 && nhits > 0) {
               const double mE = (hE->GetBinContent(er + 1, t) * sc[t]) / (nhits * sc[t]);   // as in step 3 (edge factor cancels)
               if (isSigE[er][t]) { exs.push_back(t); eys.push_back(mE); }
               else               { ex.push_back(t);  ey.push_back(mE); }
            }
         }
         cout << "Example run " << rn << " (index " << er << ", " << evx << " events): dead " << nDead[er] << ", hot " << nHot[er]
              << ", <nHits> 3sigma " << nSig[er] << ", <E> 3sigma " << nSigE[er] << endl;

         const char* badTag = isBadRun[er] ? " (BAD RUN)" : "";
         cx = new TCanvas("cExample", Form("example run %u", rn), 1300, 900);
         cx->Divide(1, 2);
         const int mcol[5] = {kBlack, kBlue, kRed, kGreen + 2, kOrange + 7};
         const char* mlab[5] = {"nHits/NEvents", "Dead towers", "Hot towers", "<nHits> 3sigma towers", "excluded by hand"};

         cx->cd(1);
         TH1F* fr1 = gPad->DrawFrame(0, 0, NT, xMax);
         fr1->SetTitle(Form("run %u%s:  <nHits> of the towers;Tower ID;<nHits>", rn, badTag));
         TLegend* lx1 = new TLegend(0.12, 0.60, 0.40, 0.89);
         for (int c = 0; c < 5; ++c) {
            if (gx[c].empty()) continue;
            TGraph* g = new TGraph((int)gx[c].size());
            for (size_t i = 0; i < gx[c].size(); ++i) g->SetPoint((int)i, gx[c][i], gy[c][i]);
            g->SetMarkerStyle(20);  g->SetMarkerSize(0.5);  g->SetMarkerColor(mcol[c]);
            g->Draw("P same");
            lx1->AddEntry(g, mlab[c], "p");
         }
         const double hl[4] = {deadCut, hotCut, nLo, nHi};
         for (int i = 0; i < 4; ++i) {
            TLine* lh = new TLine(0, hl[i], NT, hl[i]);
            lh->SetLineStyle(3);  lh->SetLineWidth(2);  lh->SetLineColor(i < 2 ? kMagenta : kGreen + 2);
            lh->Draw();
         }
         lx1->Draw();

         cx->cd(2);
         TH1F* fr2 = gPad->DrawFrame(0, 0, NT, xMaxE);
         fr2->SetTitle(Form("run %u%s:  <E_{T}> of the good towers;Tower ID;<E_{T}> [GeV]", rn, badTag));
         TLegend* lx2 = new TLegend(0.12, 0.72, 0.40, 0.89);
         if (!ex.empty()) {
            TGraph* ge = new TGraph((int)ex.size());
            for (size_t i = 0; i < ex.size(); ++i) ge->SetPoint((int)i, ex[i], ey[i]);
            ge->SetMarkerStyle(20);  ge->SetMarkerSize(0.5);  ge->SetMarkerColor(kBlack);
            ge->Draw("P same");
            lx2->AddEntry(ge, "<E_{T}> of good towers", "p");
         }
         if (!exs.empty()) {
            TGraph* gs = new TGraph((int)exs.size());
            for (size_t i = 0; i < exs.size(); ++i) gs->SetPoint((int)i, exs[i], eys[i]);
            gs->SetMarkerStyle(20);  gs->SetMarkerSize(0.5);  gs->SetMarkerColor(kGreen + 2);
            gs->Draw("P same");
            lx2->AddEntry(gs, "<E_{T}> 3sigma towers", "p");
         }
         if (eLandauRings || useQRings) {                          // window differs per tower (ring): draw it tower by tower
            TGraph* gwl = new TGraph(NT);
            TGraph* gwh = new TGraph(NT);
            for (int t = 1; t <= NT; ++t) { gwl->SetPoint(t - 1, t, wLo[t]);  gwh->SetPoint(t - 1, t, wHi[t]); }
            gwl->SetLineColor(kGreen + 2);  gwh->SetLineColor(kGreen + 2);
            gwl->Draw("L same");  gwh->Draw("L same");
         } else {
            const double el[2] = {eLo, eHi};
            for (int i = 0; i < 2; ++i) {
               TLine* le = new TLine(0, el[i], NT, el[i]);
               le->SetLineStyle(3);  le->SetLineWidth(2);  le->SetLineColor(kGreen + 2);
               le->Draw();
            }
         }
         lx2->Draw();
         cx->SaveAs(Form("%s_example_%u.pdf", outTag, rn));
         cout << "Wrote " << outTag << "_example_" << rn << ".pdf" << endl;
      }
   }

   // ---------- store the results in <outTag>.root ----------
   // tower status map: x = run index (bin r = run index r), y = tower ID (bin t = tower t)
   //   0 good (or run not analysed), 1 bad tower in this run (dead, hot, <nHits> 3sigma or <E_T> 3sigma)
   TH2C* hStatus = new TH2C("hTowerStatus", "tower status;runIndex;towerID", nRuns, 0.5, nRuns + 0.5, NT, 0.5, NT + 0.5);
   // same map after excluding the bad runs: 0 good tower in a good run, 1 bad tower, 5 = run excluded (bad run or not analysed)
   TH2C* hStatusGood = new TH2C("hTowerStatusGood", "tower status, bad runs excluded;runIndex;towerID", nRuns, 0.5, nRuns + 0.5, NT, 0.5, NT + 0.5);
   for (int r = 1; r <= nRuns; ++r) {
      const bool excluded = (!analysed[r] || isBadRun[r]);
      for (int t = 1; t <= NT; ++t) {
         const int code = (analysed[r] && (isDead[r][t] || isHot[r][t] || isSig[r][t] || isSigE[r][t] || isMan[r][t])) ? 1 : 0;
         if (code > 0) hStatus->SetBinContent(r, t, code);
         const int codeG = excluded ? 5 : code;
         if (codeG > 0) hStatusGood->SetBinContent(r, t, codeG);
      }
   }
   // fixed colour range 0..5 (stored in the histograms): bad tower (1) = blue, excluded run (5) = yellow in the default palette
   hStatus->SetMinimum(0);      hStatus->SetMaximum(5);
   hStatusGood->SetMinimum(0);  hStatusGood->SetMaximum(5);
   // final map (step 6): 0 good, 1 bad tower (bad in every run), 5 run excluded (step 5 bad run, step 6 or not analysed)
   TH2C* hStatusFinal = new TH2C("hTowerStatusFinal", "final tower status;runIndex;towerID", nRuns, 0.5, nRuns + 0.5, NT, 0.5, NT + 0.5);
   for (int r = 1; r <= nRuns; ++r) {
      const bool excluded = (!analysed[r] || isBadRun[r] || runCleaned[r]);
      for (int t = 1; t <= NT; ++t) {
         const int codeF = excluded ? 5 : (towerBadFinal[t] ? 1 : 0);
         if (codeF > 0) hStatusFinal->SetBinContent(r, t, codeF);
      }
   }
   hStatusFinal->SetMinimum(0);  hStatusFinal->SetMaximum(5);
   TCanvas* csf = new TCanvas("cStatusFinal", "final tower status", 1000, 700);
   TH2C* hStatusFinalDraw = (TH2C*) hStatusFinal->Clone("hTowerStatusFinal_draw");   // draw a copy: hTowerStatusFinal is written untouched, like hTowerStatusGood
   hStatusFinalDraw->SetStats(0);
   hStatusFinalDraw->Draw("COL");
   csf->SaveAs(Form("%s_status_final.pdf", outTag));
   csf->SaveAs(Form("%s_status_final.svg", outTag));
   csf->SaveAs(Form("%s_status_final.png", outTag));                 // light-weight image, opens anywhere
   cout << "Wrote " << outTag << "_status_final.pdf, svg and png" << endl;
   // the same information as two small 1D histograms (the 2D map has ~10M bins, too big for web viewers such as the VS Code ROOT viewer)
   TH1C* hFinalTower = new TH1C("hFinalBadTowers", "final bad towers (1 = bad in every run);towerID;bad", NT, 0.5, NT + 0.5);
   TH1C* hFinalRun   = new TH1C("hFinalRunStatus", "final run status (0 = good, 5 = excluded);runIndex;status", nRuns, 0.5, nRuns + 0.5);
   for (int t = 1; t <= NT; ++t) if (towerBadFinal[t]) hFinalTower->SetBinContent(t, 1);
   for (int r = 1; r <= nRuns; ++r) if (!analysed[r] || isBadRun[r] || runCleaned[r]) hFinalRun->SetBinContent(r, 5);

   TFile* fo = new TFile(Form("%s.root", outTag), "RECREATE");
   TTree* tr = new TTree("runs", "run characteristics");
   unsigned int tRunId = 0;
   int tIdx = 0, tDead = 0, tHot = 0, tSig = 0, tSigE = 0, tBad = 0;
   double tEv = 0.;
   tr->Branch("runId", &tRunId, "runId/i");
   tr->Branch("runIndex", &tIdx, "runIndex/I");
   tr->Branch("nEvents", &tEv, "nEvents/D");
   tr->Branch("nDead", &tDead, "nDead/I");
   tr->Branch("nHot", &tHot, "nHot/I");
   tr->Branch("nSigNHits", &tSig, "nSigNHits/I");
   tr->Branch("nSigE", &tSigE, "nSigE/I");
   tr->Branch("badRun", &tBad, "badRun/I");             // 1 = bad run (step 5)
   for (int r = 1; r <= nRuns; ++r) {
      if (!analysed[r]) continue;
      tRunId = (r < (int)runId.size()) ? runId[r] : 0u;
      tIdx = r;  tEv = hEv->GetBinContent(r + 1);
      tDead = nDead[r];  tHot = nHot[r];  tSig = nSig[r];  tSigE = nSigE[r];  tBad = isBadRun[r];
      tr->Fill();
   }
   fo->cd();
   tr->Write();
   for (int i = 0; i < 4; ++i) hRun[i]->Write();       // number of dead / hot / <nHits> 3sigma / <E> 3sigma towers per run
   hAll->Write();                                      // <nHits> of the good towers after step 1
   hCut->Write();                                      // <nHits> after step 2
   hEall->Write();                                     // <E> of the good towers after step 2
   hEfine->Write();                                    // the same, fine binning, with the Landau fit attached
   hEcut->Write();                                     // <E> after step 4
   hStatus->Write();
   hStatusGood->Write();                               // status map: 1 bad tower, 5 run excluded
   hStatusFinal->Write();                              // FINAL map (step 6): 1 bad tower in every run, 5 run excluded
   hFinalTower->Write();                               // final bad towers, 1D (the rows of the final map)
   hFinalRun->Write();                                 // final run status, 1D (the columns of the final map)
   hTowerBad->Write();                                 // per tower: number of kept runs in which it is bad
   hEtaPhi->Write();                                   // eta-phi map of the towers (2 = edge towers)
   for (int i = 0; i < 4; ++i) hRunG[i]->Write();      // run characteristics of the good runs only
   cr->Write();
   crg->Write();
   if (cx) cx->Write();                                // example run
   if (cLan) cLan->Write();
   if (cLanR) {                                        // per-ring Landau fits
      cLanR->Write();
      if (cLanRG) cLanRG->Write();
      if (cMPV) cMPV->Write();
      for (int i = 0; i < NRING; ++i) hRingE[i]->Write();
      hRingMPV->Write();  hRingSig->Write();
   }                                // Landau fit of <E_T>
   if (cring) cring->Write();                          // mean <nHits> and <E> per eta ring, before / after the edge scaling
   // cuts that were used
   TParameter<double>("deadCut", deadCut).Write();
   TParameter<double>("hotCut", hotCut).Write();
   TParameter<double>("minEvents", minEvents).Write();
   TParameter<double>("nHitsQuantileLow", qv[0]).Write();
   TParameter<double>("nHitsQuantileHigh", qv[2]).Write();
   TParameter<double>("meanEQuantileLow", eqv[0]).Write();
   TParameter<double>("meanEQuantileHigh", eqv[2]).Write();
   TParameter<double>("cutMode", cutMode).Write();                   // 0 quantiles, 1 manual, 2 hist mean +- nSigmaHist x RMS
   TParameter<double>("nSigmaHist", nSigmaHist).Write();
   TParameter<double>("eLandau", useLandau ? 1. : 0.).Write();      // 1 = step 4 window from the Landau fit
   TParameter<double>("nSigmaLandau", nSigmaLandau).Write();
   TParameter<double>("landauMPV", lMPV).Write();
   TParameter<double>("landauSigma", lSig).Write();
   TParameter<double>("nHitsCutLow", nLo).Write();                   // cuts actually applied in step 2
   TParameter<double>("nHitsCutHigh", nHi).Write();
   TParameter<double>("meanECutLow", eLo).Write();                   // cuts actually applied in step 4
   TParameter<double>("meanECutHigh", eHi).Write();
   TParameter<double>("deadRunCut", deadRunCut).Write();
   TParameter<double>("hotRunCut", hotRunCut).Write();
   TParameter<double>("sigRunCut", sigRunCut).Write();
   TParameter<double>("sigERunCut", sigERunCut).Write();
   TParameter<double>("edgeFactor", edgeFactor).Write();
   TParameter<double>("towerCleanCut", towerCleanCut).Write();
   TParameter<double>("eLandauRings", eLandauRings ? 1. : 0.).Write();   // 1 = step 4 windows from per-ring Landau fits
   TNamed("excludeTowers", excludeTowers).Write();                      // towers excluded by hand (step 1b)
   TParameter<double>("eQuantileRings", useQRings ? 1. : 0.).Write();    // 1 = step 4 windows from per-ring quantiles
   if (cQR) cQR->Write();
   if (cQRe) cQRe->Write();
   fo->Close();
   cout << "Wrote " << outTag << ".root" << endl;

   // ---------- summary ----------
   FILE* ft = fopen(Form("%s_runs.txt", outTag), "w");
   fprintf(ft, "# runId runIndex nEvents nDead nHot nSigNHits nSigE badRun\n");
   for (int r = 1; r <= nRuns; ++r)
      if (analysed[r]) fprintf(ft, "%u %d %.0f %d %d %d %d %d\n", (r < (int)runId.size()) ? runId[r] : 0u, r, hEv->GetBinContent(r + 1), nDead[r], nHot[r], nSig[r], nSigE[r], (int)isBadRun[r]);
   fclose(ft);
   cout << "Wrote " << outTag << "_runs.txt" << endl;
}