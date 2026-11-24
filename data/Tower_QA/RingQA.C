// RingQA.C -- look at the <E_T> = sum(E_T) / NHits of the towers of ONE |eta| ring, tower by tower and run by run,
//             to find which towers / runs make the high-<E_T> entries (e.g. the 1 - 5.5 GeV bump in the |eta| = 0.325 ring)
//   root -l -b -q 'RingQA.C+("merged.root")'
//
// Input: the same file as TowerQA2.C (MapTower_RunID_NHits_TowerID, MapTower_RunID_Energy_TowerID, MapTower_RunID_NEvents,
//        run index r in x bin r+1, tower ID = y bin), the same run list and TowerMap.h next to this macro.
//
// Only the dead / hot cut of TowerQA2 step 1 is applied (towers with <nHits> < deadCut or > hotCut are skipped);
// the <nHits> 3sigma cut is NOT applied, so the high-<E_T> entries are seen even if TowerQA2 removes some of them earlier.
//
// Output:
//   <outTag>.pdf (several pages):
//     1. <E_T> of the ring: all towers / towers with high entries / the other towers (log y)
//     2. number of high-<E_T> runs per tower of the ring (vs tower ID)
//     3. <E_T> vs run index, all towers of the ring (2D)  -> is the problem limited to a run period?
//     4. eta-phi position of the ring towers, colour = number of high-<E_T> runs  -> one area / one module?
//     5. the nTop worst towers: <E_T> and <nHits> vs run index
//   <outTag>_towers.txt : per tower of the ring: number of analysed runs, number of high runs, first / last high run, mean <E_T>
//   <outTag>_entries.txt: every high entry (runId, runIndex, towerID, <E_T>, <nHits>, NHits, nEvents)
//   <outTag>.root       : all histograms and graphs

#include <vector>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TowerMap.h"
using namespace std;

void RingQA(const char* inFile = "Isobars_full.root",
            const char* runListFile = "/gpfs/mnt/gpfs01/star/pwg/walz/Isobars2018/data/Isobars_2018.txt",
            const char* outTag = "RingQA_eta0325",
            double ringAbsEta = 0.325,   // |eta| of the ring to look at
            double eHigh   = 1.0,        // <E_T> above this [GeV] counts as a high entry
            int    minEvents = 0,        // runs with fewer events are skipped (as in TowerQA2)
            double deadCut = 0.02,       // as TowerQA2 step 1
            double hotCut  = 0.40,
            int    nTop    = 6,          // number of worst towers drawn on the last page
            int    nBinsE  = 1000,       // <E_T> histograms: bins ...
            double xMaxE   = 5.5)        //                    ... over 0 .. xMaxE [GeV]
{
   gStyle->SetOptStat(0);
   TFile* f = TFile::Open(inFile);
   if (!f || f->IsZombie()) { cout << "Cannot open " << inFile << endl; return; }
   TH2D* hN  = (TH2D*) f->Get("MapTower_RunID_NHits_TowerID");
   TH2D* hE  = (TH2D*) f->Get("MapTower_RunID_Energy_TowerID");
   TH1D* hEv = (TH1D*) f->Get("MapTower_RunID_NEvents");
   if (!hN || !hE || !hEv) { cout << "Missing MapTower_* histograms in " << inFile << endl; return; }

   vector<unsigned int> runId(1, 0);                           // run index = line number, starting at 1
   {
      ifstream in(runListFile);
      if (!in.is_open()) cout << "WARNING: cannot open run list " << runListFile << " - run numbers will be 0." << endl;
      unsigned int runNum;
      while (in >> runNum) runId.push_back(runNum);
   }
   const int NT = hN->GetNbinsY();
   const int nRuns = hN->GetNbinsX() - 1;

   // ---------- towers of the ring ----------
   vector<int> ringTowers;
   vector<double> tEta(NT + 1, 0.), tPhi(NT + 1, 0.);
   for (int t = 1; t <= NT; ++t) {
      double e = 0., p = 0.;
      if (!TowerMap::Get(t, e, p)) continue;
      tEta[t] = e;  tPhi[t] = p;
      if (fabs(fabs(e) - ringAbsEta) < 0.0005) ringTowers.push_back(t);
   }
   cout << "Ring |eta| = " << ringAbsEta << ": " << ringTowers.size() << " towers (expected 240), " << nRuns << " runs" << endl;
   if (ringTowers.empty()) { cout << "No tower found at this |eta| - check ringAbsEta." << endl; return; }
   const int NR = (int) ringTowers.size();

   // ---------- per tower-run <E_T> ----------
   TH1D* hAll   = new TH1D("hET_ring",      Form("<E_{T}>, |#eta| = %.3f;<E_{T}> = #SigmaE_{T} / NHits [GeV];tower x run entries", ringAbsEta), nBinsE, 0., xMaxE);
   TH1D* hBadT  = new TH1D("hET_badTowers", "towers with high entries", nBinsE, 0., xMaxE);
   TH1D* hGoodT = new TH1D("hET_otherTowers", "the other towers", nBinsE, 0., xMaxE);
   TH2D* hRunE  = new TH2D("hET_vs_run", Form("<E_{T}> vs run, |#eta| = %.3f;run index;<E_{T}> [GeV]", ringAbsEta), nRuns, 0.5, nRuns + 0.5, 110, 0., xMaxE);
   TH1D* hNHigh = new TH1D("hNHighRuns", Form("runs with <E_{T}> > %.2f GeV per tower, |#eta| = %.3f;tower ID;number of runs", eHigh, ringAbsEta), NT, 0.5, NT + 0.5);
   TH1D* hHighRun = new TH1D("hNHighTowers", Form("towers of the ring with <E_{T}> > %.2f GeV per run;run index;number of towers", eHigh), nRuns, 0.5, nRuns + 0.5);

   vector<int> nAna(NT + 1, 0), nHigh(NT + 1, 0), firstHigh(NT + 1, 0), lastHigh(NT + 1, 0);
   vector<double> sumE(NT + 1, 0.);
   FILE* fe = fopen(Form("%s_entries.txt", outTag), "w");
   fprintf(fe, "# entries of the |eta| = %.3f ring with <E_T> > %.3f GeV (dead / hot towers skipped)\n", ringAbsEta, eHigh);
   fprintf(fe, "# runId runIndex towerID eta phi <E_T>[GeV] <nHits> NHits nEvents\n");
   long nHighTot = 0, nEntries = 0;
   for (int r = 1; r <= nRuns; ++r) {
      const double ev = hEv->GetBinContent(r + 1);
      if (ev < minEvents || ev <= 0) continue;
      for (int k = 0; k < NR; ++k) {
         const int t = ringTowers[k];
         const double nhits = hN->GetBinContent(r + 1, t);
         const double nh = nhits / ev;
         if (nh < deadCut || nh > hotCut || nhits <= 0) continue;      // dead / hot as in TowerQA2 step 1
         const double et = hE->GetBinContent(r + 1, t) / nhits;
         ++nAna[t];  sumE[t] += et;  ++nEntries;
         hAll->Fill(et);
         hRunE->Fill(r, et);
         if (et > eHigh) {
            ++nHigh[t];  ++nHighTot;
            if (!firstHigh[t]) firstHigh[t] = r;
            lastHigh[t] = r;
            hHighRun->Fill(r);
            fprintf(fe, "%u %d %d %.3f %.3f %.4f %.5f %.0f %.0f\n", (r < (int)runId.size()) ? runId[r] : 0u, r, t, tEta[t], tPhi[t], et, nh, nhits, ev);
         }
      }
   }
   fclose(fe);
   // split the distribution into towers with / without high entries (second pass, now that nHigh is known)
   for (int r = 1; r <= nRuns; ++r) {
      const double ev = hEv->GetBinContent(r + 1);
      if (ev < minEvents || ev <= 0) continue;
      for (int k = 0; k < NR; ++k) {
         const int t = ringTowers[k];
         const double nhits = hN->GetBinContent(r + 1, t);
         const double nh = nhits / ev;
         if (nh < deadCut || nh > hotCut || nhits <= 0) continue;
         const double et = hE->GetBinContent(r + 1, t) / nhits;
         if (nHigh[t] > 0) hBadT->Fill(et); else hGoodT->Fill(et);
      }
   }

   // ---------- per-tower list, sorted by number of high runs ----------
   vector<int> order;
   for (int k = 0; k < NR; ++k) { const int t = ringTowers[k]; hNHigh->SetBinContent(t, nHigh[t]); if (nHigh[t] > 0) order.push_back(t); }
   for (size_t i = 0; i < order.size(); ++i)                       // simple sort by nHigh, descending
      for (size_t j = i + 1; j < order.size(); ++j)
         if (nHigh[order[j]] > nHigh[order[i]]) { const int tmp = order[i]; order[i] = order[j]; order[j] = tmp; }

   FILE* ftw = fopen(Form("%s_towers.txt", outTag), "w");
   fprintf(ftw, "# towers of the |eta| = %.3f ring with <E_T> > %.3f GeV in at least one run, sorted by the number of such runs\n", ringAbsEta, eHigh);
   fprintf(ftw, "# towerID eta phi nRunsAnalysed nRunsHigh firstHighRunId lastHighRunId meanET[GeV]\n");
   for (size_t i = 0; i < order.size(); ++i) {
      const int t = order[i];
      fprintf(ftw, "%d %.3f %.3f %d %d %u %u %.4f\n", t, tEta[t], tPhi[t], nAna[t], nHigh[t],
              (firstHigh[t] < (int)runId.size()) ? runId[firstHigh[t]] : 0u, (lastHigh[t] < (int)runId.size()) ? runId[lastHigh[t]] : 0u,
              nAna[t] ? sumE[t] / nAna[t] : 0.);
   }
   fclose(ftw);

   int nRunsWithHigh = 0;
   for (int r = 1; r <= nRuns; ++r) if (hHighRun->GetBinContent(r) > 0) ++nRunsWithHigh;
   cout << "Entries analysed: " << nEntries << ";  entries with <E_T> > " << eHigh << " GeV: " << nHighTot
        << " (" << (nEntries ? 100. * nHighTot / nEntries : 0.) << " %)" << endl;
   cout << "Towers with at least one high entry: " << order.size() << " of " << NR << ";  runs with at least one high entry: " << nRunsWithHigh << endl;
   cout << "Worst towers (towerID: high runs / analysed runs, eta, phi):" << endl;
   for (size_t i = 0; i < order.size() && i < 15; ++i)
      cout << Form("  %4d: %4d / %4d   eta %+.3f  phi %+.3f", order[i], nHigh[order[i]], nAna[order[i]], tEta[order[i]], tPhi[order[i]]) << endl;

   // ---------- plots ----------
   const TString pdf = Form("%s.pdf", outTag);
   TCanvas* c1 = new TCanvas("cRingET", "<E_T> of the ring", 1000, 700);
   c1->SetLogy(1);  c1->SetGrid(1, 1);
   hAll->SetMinimum(0.5);  hAll->SetMaximum(5. * hAll->GetMaximum());
   hAll->SetLineColor(kBlack);   hAll->Draw("hist");
   hBadT->SetLineColor(kRed);    hBadT->Draw("hist same");
   hGoodT->SetLineColor(kBlue);  hGoodT->Draw("hist same");
   TLine* lh = new TLine(eHigh, 0.5, eHigh, hAll->GetMaximum());
   lh->SetLineStyle(2);  lh->SetLineColor(kGreen + 2);  lh->SetLineWidth(2);  lh->Draw();
   TLegend* l1 = new TLegend(0.50, 0.70, 0.89, 0.89);
   l1->AddEntry(hAll, "all towers of the ring", "l");
   l1->AddEntry(hBadT, Form("the %d towers with entries > %.2f GeV", (int)order.size(), eHigh), "l");
   l1->AddEntry(hGoodT, "the other towers", "l");
   l1->AddEntry(lh, Form("threshold %.2f GeV", eHigh), "l");
   l1->Draw();
   c1->Print(pdf + "(");

   TCanvas* c2 = new TCanvas("cRingNHigh", "high runs per tower", 1000, 700);
   c2->SetGrid(1, 1);
   int tMin = NT, tMax = 1;
   for (int k = 0; k < NR; ++k) { if (ringTowers[k] < tMin) tMin = ringTowers[k];  if (ringTowers[k] > tMax) tMax = ringTowers[k]; }
   hNHigh->GetXaxis()->SetRangeUser(tMin - 20, tMax + 20);
   hNHigh->SetFillColor(kRed);  hNHigh->SetLineColor(kRed);
   hNHigh->Draw("hist");
   c2->Print(pdf);

   TCanvas* c3 = new TCanvas("cRingETvsRun", "<E_T> vs run", 1200, 700);
   c3->SetLogz(1);  c3->SetRightMargin(0.12);
   hRunE->Draw("COLZ");
   TLine* lh3 = new TLine(0.5, eHigh, nRuns + 0.5, eHigh);
   lh3->SetLineStyle(2);  lh3->SetLineColor(kRed);  lh3->Draw();
   c3->Print(pdf);

   TCanvas* c4 = new TCanvas("cRingEtaPhi", "position of the ring towers", 1000, 500);
   c4->SetGrid(1, 1);
   TH1F* fr4 = gPad->DrawFrame(-TMath::Pi(), -1.05, TMath::Pi(), 1.05);
   fr4->SetTitle(Form("ring towers (grey) and towers with <E_{T}> > %.2f GeV (red, size ~ number of runs);#phi [rad];#eta", eHigh));
   TGraph* gAllT = new TGraph(NR);
   for (int k = 0; k < NR; ++k) gAllT->SetPoint(k, tPhi[ringTowers[k]], tEta[ringTowers[k]]);
   gAllT->SetMarkerStyle(24);  gAllT->SetMarkerColor(kGray + 1);  gAllT->Draw("P same");
   for (size_t i = 0; i < order.size(); ++i) {
      TGraph* gb = new TGraph(1);
      gb->SetPoint(0, tPhi[order[i]], tEta[order[i]]);
      gb->SetMarkerStyle(20);  gb->SetMarkerColor(kRed);
      gb->SetMarkerSize(0.6 + 2.4 * nHigh[order[i]] / (double) nHigh[order[0]]);
      gb->Draw("P same");
   }
   c4->Print(pdf);

   // the nTop worst towers: <E_T> and <nHits> vs run index
   const int nShow = (int) order.size() < nTop ? (int) order.size() : nTop;
   TCanvas* c5 = new TCanvas("cRingWorst", "worst towers vs run", 1400, 300 * (nShow > 0 ? nShow : 1));
   if (nShow > 0) {
      c5->Divide(2, nShow);
      for (int iw = 0; iw < nShow; ++iw) {                       // (not "i": CINT already has i as size_t above)
         const int t = order[iw];
         TGraph* gE = new TGraph();  TGraph* gH = new TGraph();
         for (int r = 1; r <= nRuns; ++r) {
            const double ev = hEv->GetBinContent(r + 1);
            if (ev < minEvents || ev <= 0) continue;
            const double nhits = hN->GetBinContent(r + 1, t);
            gH->SetPoint(gH->GetN(), r, nhits / ev);
            if (nhits > 0) gE->SetPoint(gE->GetN(), r, hE->GetBinContent(r + 1, t) / nhits);
         }
         c5->cd(2 * iw + 1);  gPad->SetGrid(1, 1);
         gE->SetTitle(Form("tower %d (#eta %+.3f, #phi %+.3f): <E_{T}>;run index;<E_{T}> [GeV]", t, tEta[t], tPhi[t]));
         gE->SetMarkerStyle(20);  gE->SetMarkerSize(0.4);  gE->SetMarkerColor(kRed);
         gE->Draw("AP");
         TLine* le = new TLine(400, eHigh, 2200, eHigh);  le->SetLineStyle(2); le->SetLineColor(kMagenta); le->Draw();
         c5->cd(2 * iw + 2);  gPad->SetGrid(1, 1);
         gH->SetTitle(Form("tower %d: <nHits>;run index;<nHits> = NHits / NEvents", t));
         gH->SetMarkerStyle(20);  gH->SetMarkerSize(0.4);  gH->SetMarkerColor(kBlue);
         gH->Draw("AP");
         TLine* ld = new TLine(0.5, deadCut, 2200, deadCut);  ld->SetLineStyle(2);  ld->SetLineColor(kMagenta);  ld->Draw();
      }
   }
   c5->Print(pdf + ")");
   cout << "Wrote " << pdf << ", " << outTag << "_towers.txt, " << outTag << "_entries.txt" << endl;

   TFile* fo = new TFile(Form("%s.root", outTag), "RECREATE");
   hAll->Write();  hBadT->Write();  hGoodT->Write();  hRunE->Write();  hNHigh->Write();  hHighRun->Write();
   c1->Write();  c2->Write();  c3->Write();  c4->Write();  c5->Write();
   fo->Close();
   cout << "Wrote " << outTag << ".root" << endl;
}