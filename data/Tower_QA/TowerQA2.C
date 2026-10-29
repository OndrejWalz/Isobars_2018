// TowerQA2.C -- bad tower / bad run QA, built step by step
//   root -l -b -q 'TowerQA2.C+("merged.root")'
//
// Input: MapTower_RunID_NHits_TowerID (x = run index in bin r+1, y = tower ID = bin number 1..4800)
//        MapTower_RunID_NEvents       (run index r in bin r+1)
//
// Before the QA: eta-phi distribution of the towers (TowerMap.h) -> <outTag>_etaphi.pdf.
// Edge towers (|eta| = 0.967, 240 towers): after the dead / hot step their <nHits> and their <E> are scaled by 1/edgeFactor
// (1/0.69), so everything from step 2 on (distributions, quantile cuts, example plot) uses the scaled values.
//
// Per-run procedure (one pass over the runs with >= minEvents events):
//   step 1: <nHits> = NHits / NEvents for every tower
//           tower is DEAD if <nHits> < deadCut (0.02), HOT if <nHits> > hotCut (0.40); both are excluded from the next steps
//   step 2: <nHits> distribution of the remaining good towers, all analysed runs together;
//           a tower with <nHits> outside the 0.135 % / 99.865 % quantiles of that distribution is bad:
//           "<nHits> 3sigma equivalent" tower (isSig); it is excluded from the next steps
//   step 3: <E> = sum(E) / NHits of the remaining good towers, distribution over all analysed runs
//   step 4: a tower with <E> outside the 0.135 % / 99.865 % quantiles of that distribution is bad:
//           "<E> 3sigma equivalent" tower (isSigE); it is excluded from the next steps
//   step 5: a run is bad if it has more than deadRunCut (1000) dead, more than hotRunCut (2) hot, more than sigRunCut (200)
//           <nHits> 3sigma or more than sigERunCut (200) <E> 3sigma towers; bad runs are removed from the next steps
//   step 6: list of bad towers (<outTag>_badtowers.txt): for every tower the number of kept (good) runs in which it was
//           dead / hot / <nHits> 3sigma / <E> 3sigma; towers bad in at least towerMinBad kept runs are written to the file
//   step 7: clean the tower status (kept runs only): a tower that is bad (dead / hot / 3sigma / 3sigmaE) in MORE than towerFinalCut (10)
//           runs is a bad tower (bad in all runs); a tower bad in 1..towerFinalCut runs is kept, only those tower-in-run entries are bad.
//           Lists: <outTag>_final_badruns.txt, _final_badtowers.txt, _final_badtowerruns.txt; final status map hTowerStatusFinal in the root file
//   (next steps are added after step 7)

#include <vector>
#include <cstdio>
#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TParameter.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLine.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TowerMap.h"                 // TowerMap::Eta(softID), TowerMap::Phi(softID); keep TowerMap.h next to this macro
using namespace std;

void TowerQA2(const char* inFile = "Tower_Test_Isobars_full.root",
              const char* runListFile = "/gpfs/mnt/gpfs01/star/pwg/walz/Isobars2018/data/Isobars_2018.txt",
              const char* outTag = "TowerQA2",
              int    minEvents = 1,    // runs with fewer good events are not analysed
              double deadCut = 0.02,     // step 1: tower is dead if NHits/NEvents < deadCut
              double hotCut  = 0.40,     // step 1: tower is hot  if NHits/NEvents > hotCut
              int    nBins = 150,        // <nHits> distribution: bins
              double xMax  = 0.15,       // <nHits> distribution: upper edge
              int    nBinsE = 200,       // <E> distribution: bins
              double xMaxE  = 2.0,       // <E> distribution: upper edge [GeV]
              int    exampleRun = 0,     // example plot: run number (e.g. 19098021) or run index (e.g. 500); 0 = first analysed run
              int    deadRunCut = 500,  // step 5: run is bad if it has MORE than this many dead towers
              int    hotRunCut  = 2,     //         ... or more than this many hot towers
              int    sigRunCut  = 40,   //         ... or more than this many <nHits> 3sigma towers
              int    sigERunCut = 15,    //         ... or more than this many <E> 3sigma towers
              int    towerMinBad = 2,    // step 6: a tower is listed as bad if it is bad in at least this many kept runs
              double edgeFactor = 0.69,  // edge towers (|eta| = 0.967): <nHits> and <E> are divided by this (= scaled by 1/0.69); 1 = no scaling
              int    towerFinalCut = 10) // step 7: a tower bad in MORE than this many kept runs is a bad tower; bad in fewer: only those tower-in-run entries are bad
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

   // ---------- eta-phi distribution of the towers (TowerMap.h) ----------
   // edge towers: |eta| = 0.967 (the last ring of each half of the BEMC). After the dead / hot step their <nHits> and <E>
   // are multiplied by sc[t] = 1/edgeFactor (1/0.69); sc[t] = 1 for all other towers.
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
   cout << "Edge towers (|eta| = 0.967): " << nEdge << " (expected 240), scale factor " << scEdge << endl;
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

   TH1D* hAll = new TH1D("hNHits_good", "<nHits> of good towers, all analysed runs;<nHits> = NHits / NEvents;tower x run entries", nBins, 0., xMax);

   // ---------- step 1: dead and hot towers ----------
   int nAnalysed = 0;
   for (int r = 1; r <= nRuns; ++r) {
      const double ev = hEv->GetBinContent(r + 1);
      if (ev < minEvents) continue;
      analysed[r] = 1;
      ++nAnalysed;

      isDead[r] = vector<char>(NT + 1, 0);                    // (no vector::assign: CINT)
      isHot[r]  = vector<char>(NT + 1, 0);
      isSig[r]  = vector<char>(NT + 1, 0);
      isSigE[r] = vector<char>(NT + 1, 0);
      alive[r]  = vector<char>(NT + 1, 0);
      for (int t = 1; t <= NT; ++t) {
         const double nh = hN->GetBinContent(r + 1, t) / ev;  // <nHits>
         if (nh < deadCut)      { isDead[r][t] = 1; ++nDead[r]; }
         else if (nh > hotCut)  { isHot[r][t]  = 1; ++nHot[r]; }
         else {                                                 // dead / hot are judged on the unscaled <nHits>; the good towers are scaled (edge towers x 1/0.69)
            alive[r][t] = 1;
            hAll->Fill(nh * sc[t]);
            if (ringOf[t] >= 0) { rnU[ringOf[t]] += nh; rnS[ringOf[t]] += nh * sc[t]; rnC[ringOf[t]] += 1.; }
         }
      }
   }
   cout << "Runs analysed: " << nAnalysed << " (skipped < " << minEvents << " events: " << nRuns - nAnalysed << ")" << endl;

   // ---------- <nHits> distribution of the good towers (all analysed runs) ----------
   const double pMean = hAll->GetMean(), pRMS = hAll->GetRMS(), pSkew = hAll->GetSkewness();
   hAll->Fit("gaus", "Q", "", pMean - 50. * pRMS, pMean + 50. * pRMS);      // whole range
   TF1* gf = hAll->GetFunction("gaus");
   const double gMu = gf->GetParameter(1), gSig = fabs(gf->GetParameter(2));
   double pq[3] = {0.00135, 0.5, 0.99865}, qv[3] = {0., 0., 0.};            // quantiles where a Gaussian has -3 sigma / median / +3 sigma
   hAll->GetQuantiles(3, qv, pq);
   cout << "Good towers: mean " << pMean << "  RMS " << pRMS << "  skewness " << pSkew << "  overflow " << hAll->GetBinContent(nBins + 1) << endl;
   cout << "Gaussian fit: mu " << gMu << "  sigma " << gSig << "  -> mu +- 3 sigma = [" << gMu - 3. * gSig << ", " << gMu + 3. * gSig << "]" << endl;
   cout << "Quantiles 0.135 % / 50 % / 99.865 %: " << qv[0] << " / " << qv[1] << " / " << qv[2]
        << "  (distance from median: -" << qv[1] - qv[0] << " / +" << qv[2] - qv[1] << ")" << endl;

   TCanvas* cd = new TCanvas("cNHitsGood", "<nHits> of good towers", 900, 650);
   cd->SetLogy(1);
   cd->SetGrid(1, 1);
   const double hmax = hAll->GetMaximum();
   const double y1 = 0.5, y2 = 5. * hmax;
   hAll->SetMinimum(y1);
   hAll->SetMaximum(y2);
   hAll->SetLineColor(kBlack);
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
   lg->AddEntry(hAll, "good towers, all runs", "l");
   lg->AddEntry(gf, Form("Gaussian fit (#sigma=%.4f)", gSig), "l");
   lg->AddEntry(l1, "#mu #pm 3#sigma (Gaussian)", "l");
   lg->AddEntry(l3, "0.135 % / 99.865 % quantiles", "l");
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
         if (nh < qv[0] || nh > qv[2]) { isSig[r][t] = 1; alive[r][t] = 0; ++nSig[r]; ++nSigTot; }
         else { ++nAliveTot; hCut->Fill(nh); }
      }
   }
   cout << "Step 2: <nHits> 3sigma equivalent towers: " << nSigTot << " tower-runs outside [" << qv[0] << ", " << qv[2]
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
   TH1D* hEall = new TH1D("hMeanE_good", "<E> of good towers, all analysed runs;<E> = #SigmaE / NHits [GeV];tower x run entries", nBinsE, 0., xMaxE);
   for (int r = 1; r <= nRuns; ++r) {
      if (!analysed[r]) continue;
      for (int t = 1; t <= NT; ++t) {
         if (!alive[r][t]) continue;
         const double nhits = hN->GetBinContent(r + 1, t);
         if (nhits <= 0) continue;
         const double mE0 = hE->GetBinContent(r + 1, t) / nhits;         // unscaled <E>
         hEall->Fill(mE0 );                                       // scaled <E> (edge towers x 1/0.69) *sc[t] turned off for now
         if (ringOf[t] >= 0) { reU[ringOf[t]] += mE0; reS[ringOf[t]] += mE0 ; reC[ringOf[t]] += 1.; } // mE0 *sc[t] turned off for now
      }
   }
   const double eMean = hEall->GetMean(), eRMS = hEall->GetRMS(), eSkew = hEall->GetSkewness();
   hEall->Fit("gaus", "Q", "", eMean - 50. * eRMS, eMean + 50. * eRMS);     // whole range
   TF1* egf = hEall->GetFunction("gaus");
   const double eMu = egf->GetParameter(1), eSig = fabs(egf->GetParameter(2));
   double eq[3] = {0.00135, 0.5, 0.99865}, eqv[3] = {0., 0., 0.};
   hEall->GetQuantiles(3, eqv, eq);
   cout << "Step 3: <E> of good towers: mean " << eMean << "  RMS " << eRMS << "  skewness " << eSkew
        << "  entries " << hEall->GetEntries() << "  overflow " << hEall->GetBinContent(nBinsE + 1) << endl;
   cout << "Gaussian fit: mu " << eMu << "  sigma " << eSig << "  -> mu +- 3 sigma = [" << eMu - 3. * eSig << ", " << eMu + 3. * eSig << "]" << endl;
   cout << "Quantiles 0.135 % / 50 % / 99.865 %: " << eqv[0] << " / " << eqv[1] << " / " << eqv[2]
        << "  (distance from median: -" << eqv[1] - eqv[0] << " / +" << eqv[2] - eqv[1] << ")" << endl;

   TCanvas* ce = new TCanvas("cMeanE", "<E> of good towers", 900, 650);
   ce->SetLogy(1);
   ce->SetGrid(1, 1);
   const double eY1 = 0.5, eY2 = 5. * hEall->GetMaximum();
   hEall->SetMinimum(eY1);
   hEall->SetMaximum(eY2);
   hEall->SetLineColor(kBlack);
   hEall->Draw("hist");
   egf->SetRange(0., xMaxE);
   egf->SetLineColor(kRed);
   egf->Draw("same");
   TLine* el1 = new TLine(eMu - 3. * eSig, eY1, eMu - 3. * eSig, 0.7 * eY2);
   TLine* el2 = new TLine(eMu + 3. * eSig, eY1, eMu + 3. * eSig, 0.7 * eY2);
   TLine* el3 = new TLine(eqv[0], eY1, eqv[0], 0.7 * eY2);
   TLine* el4 = new TLine(eqv[2], eY1, eqv[2], 0.7 * eY2);
   el1->SetLineColor(kRed);  el2->SetLineColor(kRed);
   el3->SetLineColor(kBlue); el4->SetLineColor(kBlue);
   TLine* els[4] = {el1, el2, el3, el4};
   for (int i = 0; i < 4; ++i) { els[i]->SetLineStyle(2); els[i]->SetLineWidth(2); els[i]->Draw(); }
   TLegend* leg = new TLegend(0.55, 0.62, 0.89, 0.89);
   leg->SetHeader(Form("mean %.4f  RMS %.4f  skew %.2f", eMean, eRMS, eSkew));
   leg->AddEntry(hEall, "good towers, all runs", "l");
   leg->AddEntry(egf, Form("Gaussian fit (#sigma=%.4f)", eSig), "l");
   leg->AddEntry(el1, "#mu #pm 3#sigma (Gaussian)", "l");
   leg->AddEntry(el3, "0.135 % / 99.865 % quantiles", "l");
   leg->Draw();
   ce->SaveAs(Form("%s_meanE_good.pdf", outTag));
   cout << "Wrote " << outTag << "_meanE_good.pdf" << endl;

   // ---------- check of the edge scaling: mean <nHits> and <E> per eta ring, before and after scaling ----------
   // <nHits>: towers that are not dead / hot (step 1);  <E>: towers that are also inside the <nHits> quantiles (step 2); all analysed runs
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
      cring = new TCanvas("cRingMeans", "mean <nHits> and <E> per eta ring", 900, 900);
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
         frr->SetTitle(p == 1 ? "mean <nHits> per #eta ring;|#eta|;mean <nHits>" : "mean <E> per #eta ring;|#eta|;mean <E> [GeV]");
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
   TH1D* hEcut = new TH1D("hMeanE_afterSig", "<E> of good towers after the 3sigma equivalent cut;<E> = #SigmaE / NHits [GeV];tower x run entries", nBinsE, 0., xMaxE);
   long nSigETot = 0, nAliveETot = 0;
   for (int r = 1; r <= nRuns; ++r) {
      if (!analysed[r]) continue;
      for (int t = 1; t <= NT; ++t) {
         if (!alive[r][t]) continue;
         const double nhits = hN->GetBinContent(r + 1, t);
         if (nhits <= 0) continue;
         const double meanE = hE->GetBinContent(r + 1, t) / nhits * sc[t];   // scaled <E>
         if (meanE < eqv[0] || meanE > eqv[2]) { isSigE[r][t] = 1; alive[r][t] = 0; ++nSigE[r]; ++nSigETot; }
         else { ++nAliveETot; hEcut->Fill(meanE); }
      }
   }
   cout << "Step 4: <E> 3sigma equivalent towers: " << nSigETot << " tower-runs outside [" << eqv[0] << ", " << eqv[2]
        << "] GeV;  good towers left: " << nAliveETot << " tower-runs" << endl;

   TCanvas* cf = new TCanvas("cMeanEAfter", "<E> after the 3sigma equivalent cut", 900, 650);
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

   // ---------- step 6: bad towers, counted over the kept (good) runs ----------
   vector<int> cntDead(NT + 1, 0), cntHot(NT + 1, 0), cntSig(NT + 1, 0), cntSigE(NT + 1, 0);
   int nKept = 0;
   for (int r = 1; r <= nRuns; ++r) {
      if (!analysed[r] || isBadRun[r]) continue;
      ++nKept;
      for (int t = 1; t <= NT; ++t) {
         if (isDead[r][t])      ++cntDead[t];
         else if (isHot[r][t])  ++cntHot[t];
         else if (isSig[r][t])  ++cntSig[t];
         else if (isSigE[r][t]) ++cntSigE[t];
      }
   }
   TH1D* hTowerBad = new TH1D("hTowerNBad", "number of kept runs in which the tower is bad;towerID;N runs", NT, 0.5, NT + 0.5);
   int nBadTowers = 0;
   FILE* ftw = fopen(Form("%s_badtowers.txt", outTag), "w");
   fprintf(ftw, "# bad towers: bad (dead, hot, <nHits> 3sigma or <E> 3sigma) in at least %d of the %d kept runs\n", towerMinBad, nKept);
   fprintf(ftw, "# towerID nDead nHot nSigNHits nSigE nBadTotal\n");
   for (int t = 1; t <= NT; ++t) {
      const int nbad = cntDead[t] + cntHot[t] + cntSig[t] + cntSigE[t];
      hTowerBad->SetBinContent(t, nbad);
      if (nbad < towerMinBad) continue;
      ++nBadTowers;
      fprintf(ftw, "%d %d %d %d %d %d\n", t, cntDead[t], cntHot[t], cntSig[t], cntSigE[t], nbad);
   }
   fclose(ftw);
   cout << "Step 6: bad towers (bad in >= " << towerMinBad << " of " << nKept << " kept runs): " << nBadTowers << " of " << NT << endl;
   cout << "Wrote " << outTag << "_badtowers.txt" << endl;

   // ---------- step 7: clean the tower status (kept runs only) ----------
   // tower bad in MORE than towerFinalCut kept runs  -> bad tower (bad in all runs)
   // tower bad in 1..towerFinalCut kept runs          -> the tower stays good, only these tower-in-run entries are marked bad
   vector<char> towerBadFinal(NT + 1, 0);
   int nTowFinal = 0;
   FILE* ffT = fopen(Form("%s_final_badtowers.txt", outTag), "w");
   fprintf(ffT, "# final bad towers: bad (dead, hot, <nHits> 3sigma or <E> 3sigma) in more than %d of the %d kept runs\n", towerFinalCut, nKept);
   fprintf(ffT, "# towerID nDead nHot nSigNHits nSigE nBadTotal\n");
   for (int t = 1; t <= NT; ++t) {
      const int nbad = cntDead[t] + cntHot[t] + cntSig[t] + cntSigE[t];
      if (nbad <= towerFinalCut) continue;
      towerBadFinal[t] = 1;
      ++nTowFinal;
      fprintf(ffT, "%d %d %d %d %d %d\n", t, cntDead[t], cntHot[t], cntSig[t], cntSigE[t], nbad);
   }
   fclose(ffT);

   long nTowRun = 0;                                          // tower-in-run entries that stay bad (tower itself is not bad)
   FILE* ffR = fopen(Form("%s_final_badtowerruns.txt", outTag), "w");
   fprintf(ffR, "# bad towers in runs: towers bad in 1..%d kept runs, listed for those runs only (kept runs only)\n", towerFinalCut);
   fprintf(ffR, "# runId runIndex towerID reason   (D dead, H hot, N <nHits> 3sigma, E <E> 3sigma)\n");
   for (int r = 1; r <= nRuns; ++r) {
      if (!analysed[r] || isBadRun[r]) continue;
      for (int t = 1; t <= NT; ++t) {
         if (towerBadFinal[t]) continue;
         const char* why = isDead[r][t] ? "D" : isHot[r][t] ? "H" : isSig[r][t] ? "N" : isSigE[r][t] ? "E" : 0;
         if (!why) continue;
         fprintf(ffR, "%u %d %d %s\n", (r < (int)runId.size()) ? runId[r] : 0u, r, t, why);
         ++nTowRun;
      }
   }
   fclose(ffR);

   FILE* ffB = fopen(Form("%s_final_badruns.txt", outTag), "w");   // runs removed in step 5
   fprintf(ffB, "# final bad runs (step 5)\n# runId runIndex nEvents nDead nHot nSigNHits nSigE   (reasons: D dead, H hot, N <nHits> 3sigma, E <E> 3sigma)\n");
   for (int r = 1; r <= nRuns; ++r) {
      if (!analysed[r] || !isBadRun[r]) continue;
      fprintf(ffB, "%u %d %.0f %d %d %d %d   %s%s%s%s\n", (r < (int)runId.size()) ? runId[r] : 0u, r, hEv->GetBinContent(r + 1), nDead[r], nHot[r], nSig[r], nSigE[r],
              nDead[r] > deadRunCut ? "D" : "", nHot[r] > hotRunCut ? "H" : "", nSig[r] > sigRunCut ? "N" : "", nSigE[r] > sigERunCut ? "E" : "");
   }
   fclose(ffB);
   cout << "Step 7: final bad towers (bad in > " << towerFinalCut << " of " << nKept << " kept runs): " << nTowFinal << " of " << NT
        << ";  bad tower-in-run entries of the other towers: " << nTowRun << endl;
   cout << "Wrote " << outTag << "_final_badruns.txt, _final_badtowers.txt, _final_badtowerruns.txt" << endl;

   // ---------- run characteristics: number of bad towers of each kind per run ----------
   TH1D* hRun[4];
   hRun[0] = new TH1D("hRunDead", Form("nHits/nEvents < %.0f%%;runIndex;#DeadTowers", 100. * deadCut), nRuns, 0.5, nRuns + 0.5);
   hRun[1] = new TH1D("hRunHot",  Form("nHits/nEvents > %.0f%%;runIndex;#HotTowers", 100. * hotCut), nRuns, 0.5, nRuns + 0.5);
   hRun[2] = new TH1D("hRunSig",  "mean #pm 3sigma;runIndex;#SigmaTowers", nRuns, 0.5, nRuns + 0.5);
   hRun[3] = new TH1D("hRunSigE", "meanE #pm 3sigmaE;runIndex;#SigmaETowers", nRuns, 0.5, nRuns + 0.5);
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
   hRunG[2] = new TH1D("hRunSigGood",  "good runs: mean #pm 3sigma;runIndex;#SigmaTowers", nRuns, 0.5, nRuns + 0.5);
   hRunG[3] = new TH1D("hRunSigEGood", "good runs: meanE #pm 3sigmaE;runIndex;#SigmaETowers", nRuns, 0.5, nRuns + 0.5);
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
         vector<vector<double> > gx(4), gy(4);
         vector<double> ex, ey, exs, eys;                          // <E> panel: good, <E> 3sigma equivalent
         for (int t = 1; t <= NT; ++t) {
            const double nhits = hN->GetBinContent(er + 1, t);
            int cl = 0;
            if (isDead[er][t]) cl = 1; else if (isHot[er][t]) cl = 2; else if (isSig[er][t]) cl = 3;
            const double nhx = (cl == 1 || cl == 2) ? nhits / evx : nhits / evx * sc[t];   // dead / hot: unscaled; the others: scaled
            gx[cl].push_back(t);  gy[cl].push_back(nhx);
            if (cl == 0 && nhits > 0) {
               const double mE = hE->GetBinContent(er + 1, t) / nhits * sc[t];
               if (isSigE[er][t]) { exs.push_back(t); eys.push_back(mE); }
               else               { ex.push_back(t);  ey.push_back(mE); }
            }
         }
         cout << "Example run " << rn << " (index " << er << ", " << evx << " events): dead " << nDead[er] << ", hot " << nHot[er]
              << ", <nHits> 3sigma " << nSig[er] << ", <E> 3sigma " << nSigE[er] << endl;

         const char* badTag = isBadRun[er] ? " (BAD RUN)" : "";
         cx = new TCanvas("cExample", Form("example run %u", rn), 1300, 900);
         cx->Divide(1, 2);
         const int mcol[4] = {kBlack, kBlue, kRed, kGreen + 2};
         const char* mlab[4] = {"nHits/NEvents", "Dead towers", "Hot towers", "<nHits> 3sigma towers"};

         cx->cd(1);
         TH1F* fr1 = gPad->DrawFrame(0, 0, NT, xMax);
         fr1->SetTitle(Form("run %u%s:  <nHits> of the towers;Tower ID;<nHits>", rn, badTag));
         TLegend* lx1 = new TLegend(0.12, 0.60, 0.40, 0.89);
         for (int c = 0; c < 4; ++c) {
            if (gx[c].empty()) continue;
            TGraph* g = new TGraph((int)gx[c].size());
            for (size_t i = 0; i < gx[c].size(); ++i) g->SetPoint((int)i, gx[c][i], gy[c][i]);
            g->SetMarkerStyle(20);  g->SetMarkerSize(0.5);  g->SetMarkerColor(mcol[c]);
            g->Draw("P same");
            lx1->AddEntry(g, mlab[c], "p");
         }
         const double hl[4] = {deadCut, hotCut, qv[0], qv[2]};
         for (int i = 0; i < 4; ++i) {
            TLine* lh = new TLine(0, hl[i], NT, hl[i]);
            lh->SetLineStyle(3);  lh->SetLineWidth(2);  lh->SetLineColor(i < 2 ? kMagenta : kGreen + 2);
            lh->Draw();
         }
         lx1->Draw();

         cx->cd(2);
         TH1F* fr2 = gPad->DrawFrame(0, 0, NT, xMaxE);
         fr2->SetTitle(Form("run %u%s:  <E> of the good towers;Tower ID;<E> [GeV]", rn, badTag));
         TLegend* lx2 = new TLegend(0.12, 0.72, 0.40, 0.89);
         if (!ex.empty()) {
            TGraph* ge = new TGraph((int)ex.size());
            for (size_t i = 0; i < ex.size(); ++i) ge->SetPoint((int)i, ex[i], ey[i]);
            ge->SetMarkerStyle(20);  ge->SetMarkerSize(0.5);  ge->SetMarkerColor(kBlack);
            ge->Draw("P same");
            lx2->AddEntry(ge, "<E> of good towers", "p");
         }
         if (!exs.empty()) {
            TGraph* gs = new TGraph((int)exs.size());
            for (size_t i = 0; i < exs.size(); ++i) gs->SetPoint((int)i, exs[i], eys[i]);
            gs->SetMarkerStyle(20);  gs->SetMarkerSize(0.5);  gs->SetMarkerColor(kGreen + 2);
            gs->Draw("P same");
            lx2->AddEntry(gs, "<E> 3sigma towers", "p");
         }
         const double el[2] = {eqv[0], eqv[2]};
         for (int i = 0; i < 2; ++i) {
            TLine* le = new TLine(0, el[i], NT, el[i]);
            le->SetLineStyle(3);  le->SetLineWidth(2);  le->SetLineColor(kGreen + 2);
            le->Draw();
         }
         lx2->Draw();
         cx->SaveAs(Form("%s_example_%u.pdf", outTag, rn));
         cout << "Wrote " << outTag << "_example_" << rn << ".pdf" << endl;
      }
   }

   // ---------- store the results in <outTag>.root ----------
   // tower status map: x = run index (bin r = run index r), y = tower ID (bin t = tower t)
   //   0 good (or run not analysed), 1 dead, 2 hot, 3 <nHits> 3sigma equivalent, 4 <E> 3sigma equivalent
   TH2C* hStatus = new TH2C("hTowerStatus", "tower status;runIndex;towerID", nRuns, 0.5, nRuns + 0.5, NT, 0.5, NT + 0.5);
   // same map after excluding the bad runs: 0 good tower in a good run, 1-4 as above, 5 = run excluded (bad run or not analysed)
   TH2C* hStatusGood = new TH2C("hTowerStatusGood", "tower status, bad runs excluded;runIndex;towerID", nRuns, 0.5, nRuns + 0.5, NT, 0.5, NT + 0.5);
   // final map (step 7): 0 good, 1-4 bad tower in this run (tower itself ok), 5 = run excluded, 6 = bad tower (bad in more than towerFinalCut runs)
   TH2C* hStatusFinal = new TH2C("hTowerStatusFinal", "final tower status;runIndex;towerID", nRuns, 0.5, nRuns + 0.5, NT, 0.5, NT + 0.5);
   for (int r = 1; r <= nRuns; ++r) {
      const bool excluded = (!analysed[r] || isBadRun[r]);
      for (int t = 1; t <= NT; ++t) {
         int code = 0;
         if (analysed[r]) {
            if (isDead[r][t])      code = 1;
            else if (isHot[r][t])  code = 2;
            else if (isSig[r][t])  code = 3;
            else if (isSigE[r][t]) code = 4;
         }
         if (code > 0) hStatus->SetBinContent(r, t, code);
         const int codeG = excluded ? 5 : code;
         if (codeG > 0) hStatusGood->SetBinContent(r, t, codeG);
         const int codeF = excluded ? 5 : (towerBadFinal[t] ? 6 : code);
         if (codeF > 0) hStatusFinal->SetBinContent(r, t, codeF);
      }
   }

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
   hAll->Write();                                      // <nHits> of the good towers after step 1 (with the Gaussian fit)
   hCut->Write();                                      // <nHits> after step 2
   hEall->Write();                                     // <E> of the good towers after step 2 (with the Gaussian fit)
   hEcut->Write();                                     // <E> after step 4
   hStatus->Write();
   hStatusGood->Write();                               // status map with the bad runs excluded (code 5)
   hStatusFinal->Write();                              // FINAL tower / run status map (0 good, 1-4 bad in this run, 5 run excluded, 6 bad tower)
   hTowerBad->Write();                                 // per tower: number of kept runs in which it is bad
   hEtaPhi->Write();                                   // eta-phi map of the towers (2 = edge towers)
   for (int i = 0; i < 4; ++i) hRunG[i]->Write();      // run characteristics of the good runs only
   cr->Write();
   crg->Write();
   if (cx) cx->Write();                                // example run
   if (cring) cring->Write();                          // mean <nHits> and <E> per eta ring, before / after the edge scaling
   // cuts that were used
   TParameter<double>("deadCut", deadCut).Write();
   TParameter<double>("hotCut", hotCut).Write();
   TParameter<double>("minEvents", minEvents).Write();
   TParameter<double>("nHitsQuantileLow", qv[0]).Write();
   TParameter<double>("nHitsQuantileHigh", qv[2]).Write();
   TParameter<double>("meanEQuantileLow", eqv[0]).Write();
   TParameter<double>("meanEQuantileHigh", eqv[2]).Write();
   TParameter<double>("deadRunCut", deadRunCut).Write();
   TParameter<double>("hotRunCut", hotRunCut).Write();
   TParameter<double>("sigRunCut", sigRunCut).Write();
   TParameter<double>("sigERunCut", sigERunCut).Write();
   TParameter<double>("towerMinBad", towerMinBad).Write();
   TParameter<double>("edgeFactor", edgeFactor).Write();
   TParameter<double>("towerFinalCut", towerFinalCut).Write();
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