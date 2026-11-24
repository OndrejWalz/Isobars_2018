// plotTowerMap.C -- draw the BEMC tower map in eta-phi
// Usage (TowerMap.h in the same directory):
//   root -l plotTowerMap.C                 // signed eta, colour = softID
//   root -l 'plotTowerMap.C(false)'        // without labels
//   root -l 'plotTowerMap.C(true,true)'    // |eta|: West+East merged into one plot,
//                                          // bin = sum of the two softIDs,
//                                          // labels = "WestID/EastID" of each row
// Note: no ?: with '::' anywhere -- CINT (ROOT 5) mis-parses them.

#include "TowerMap.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"
#include <cstdio>

void plotTowerMap(bool labels = true, bool absEta = false)
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);   // rainbow; kViridis needs ROOT >= 6.04

    int nEta = 40;  double etaMin = -1.;
    if (absEta) { nEta = 20; etaMin = 0.; }

    TH2D *h = new TH2D("hTowerMap", "", nEta, etaMin, 1., 120, -TMath::Pi(), TMath::Pi());
    TH2D *hCount = (TH2D*)h->Clone("hTowerCount");
    h->GetYaxis()->SetTitle("#phi [rad]");
    if (absEta) {
        h->SetTitle("BEMC tower map, |#eta| (West + East)");
        h->GetXaxis()->SetTitle("|#eta|");
        h->GetZaxis()->SetTitle("softID_{West} + softID_{East}");
    } else {
        h->SetTitle("BEMC tower map   (West #eta>0: softID 1-2400, East #eta<0: 2401-4800)");
        h->GetXaxis()->SetTitle("#eta");
        h->GetZaxis()->SetTitle("softID");
    }

    for (int id = 1; id <= TowerMap::kNTowers; ++id) {
        double eta, phi;
        if (!TowerMap::Get(id, eta, phi)) continue;
        double x = eta;
        if (absEta) x = TMath::Abs(eta);
        hCount->Fill(x, phi);
        h->Fill(x, phi, id);
    }

    // sanity check: 1 tower per bin (signed eta) or 2 per bin (|eta|)
    int expected = 1;
    if (absEta) expected = 2;
    int bad = 0;
    for (int i = 1; i <= nEta; ++i)
        for (int j = 1; j <= 120; ++j)
            if (hCount->GetBinContent(i, j) != expected) ++bad;
    printf("bins with != %d tower(s): %d\n", expected, bad);

    TCanvas *c = new TCanvas("cTowerMap", "BEMC tower map", 1400, 1000);
    c->SetRightMargin(0.13);

    h->Draw("COLZ");
    if (!absEta) {
        TLine *l = new TLine(0, -TMath::Pi(), 0, TMath::Pi());
        l->SetLineColor(kRed);
        l->SetLineWidth(2);
        l->Draw();
    }

    if (labels) {
        // first tower of each 20-tower eta row, placed at its lowest-|eta| bin
        TLatex t;
        t.SetTextSize(0.006);
        t.SetTextAlign(22);
        t.SetTextColor(kWhite);
        for (int id = 1; id <= TowerMap::kNTowers; id += 20) {
            double eta = TowerMap::Eta(id), phi = TowerMap::Phi(id);
            if (!absEta) {
                if (eta > 0) t.DrawLatex(0.025, phi, Form("%d", id));
                else         t.DrawLatex(-0.025, phi, Form("%d", id));
                continue;
            }
            if (eta < 0) continue;   // West rows only; find East row at same phi
            int idE = -1;
            for (int k = TowerMap::kNTowers / 2 + 1; k <= TowerMap::kNTowers; k += 20)
                if (TMath::Abs(TowerMap::Phi(k) - phi) < 1e-3) { idE = k; break; }
            t.DrawLatex(0.025, phi, Form("%d/%d", id, idE));
        }
    }

    if (absEta) c->SaveAs("bemc_tower_map_abseta_root.pdf");
    else        c->SaveAs("bemc_tower_map_root.pdf");
}