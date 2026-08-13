// divideJetSpectra.C
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>
#include <TList.h>
#include <TDirectory.h>
#include <iostream>
#include <vector>
#include <TPaveText.h>

// Helper: convert internal centrality string to a nice label
TString centLabel(TString cent) {
    if (cent == "CENT_0_10")  return "0-10%";
    if (cent == "MID_20_40")  return "20-40%";
    if (cent == "PERI_60_80") return "60-80%";
    return cent;
}

// Helper: convert internal radius string to a nice label
TString radLabel(TString R) {
    if (R == "R0.2") return "R = 0.2";
    if (R == "R0.3") return "R = 0.3";
    if (R == "R0.4") return "R = 0.4";
    return R;
}

// Helper: read total number of events from data_merged_<species>.root
// Structure: file -> "QA_histograms" (TList) -> "hevents" (TH1)
// (The "trees" directory is a filesystem directory holding this file on
// disk, not a TDirectory inside the ROOT file — resolved by the caller.)
// NOTE: this assumes hevents is a simple event-count histogram and uses the
// full integral as the event count. If hevents is a cutflow histogram with
// multiple stages, change this to read a specific bin instead
// (e.g. hEvents->GetBinContent(1)).
double getNEvents(TString fileName) {
    TFile *f = TFile::Open(fileName, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open events file: " << fileName << "\n";
        return -1.0;
    }

    TList *qaList = dynamic_cast<TList*>(f->Get("QA_histograms"));
    if (!qaList) {
        std::cerr << "Cannot find top-level TList 'QA_histograms' in " << fileName << "\n";
        f->Close();
        return -1.0;
    }

    TH1 *hEvents = dynamic_cast<TH1*>(qaList->FindObject("hevents"));
    if (!hEvents) {
        std::cerr << "Cannot find 'hevents' in QA_histograms in " << fileName << "\n";
        f->Close();
        return -1.0;
    }

    double nEvents = hEvents->Integral();
    std::cout << "  " << fileName << " -> nEvents (hevents integral) = " << nEvents << "\n";

    f->Close();
    return nEvents;
}

void divideJetSpectra(TString fileNum       = "hists_data_Ru.root",
                      TString fileDen        = "hists_data_Zr.root",
                      TString fileOut        = "ratio_Ru_over_Zr.root",
                      TString labelNum       = "Ru",
                      TString labelDen       = "Zr",
                      TString pdfDir         = "ratio_Ru_over_Zr",
                      TString fileEventsNum  = "../../trees/data_merged_Ru.root",
                      TString fileEventsDen  = "../../trees/data_merged_Zr.root") {

    // ──────────────────────────────────────────
    // 1.  Input files
    // ──────────────────────────────────────────
    TFile *fNum = TFile::Open(fileNum, "READ");
    TFile *fDen = TFile::Open(fileDen, "READ");

    if (!fNum || fNum->IsZombie()) { std::cerr << "Cannot open numerator file: "   << fileNum << "\n"; return; }
    if (!fDen || fDen->IsZombie()) { std::cerr << "Cannot open denominator file: " << fileDen << "\n"; return; }

    std::cout << "Numerator   : " << fileNum  << "  (" << labelNum << ")\n";
    std::cout << "Denominator : " << fileDen  << "  (" << labelDen << ")\n";
    std::cout << "Output      : " << fileOut  << "\n";
    std::cout << "PDF folder  : " << pdfDir   << "/\n\n";

    // ──────────────────────────────────────────
    // 1b.  Event counts for per-event normalization
    // ──────────────────────────────────────────
    std::cout << "Reading event counts:\n";
    double nEventsNum = getNEvents(fileEventsNum);
    double nEventsDen = getNEvents(fileEventsDen);

    if (nEventsNum <= 0 || nEventsDen <= 0) {
        std::cerr << "\nInvalid event count(s), aborting. "
                   << "nEventsNum=" << nEventsNum
                   << " nEventsDen=" << nEventsDen << "\n";
        return;
    }
    std::cout << "\n";

    // ──────────────────────────────────────────
    // 2.  Dimension arrays
    // ──────────────────────────────────────────
    std::vector<TString> radii;
    radii.push_back("R0.2");
    radii.push_back("R0.3");
    radii.push_back("R0.4");

    std::vector<TString> centralities;
    centralities.push_back("CENT_0_10");
    centralities.push_back("MID_20_40");
    centralities.push_back("PERI_60_80");

    std::vector<int> ptLeads;
    ptLeads.push_back(0);
    ptLeads.push_back(5);
    ptLeads.push_back(7);
    //ptLeads.push_back(9);

    std::vector<int> colors;
    colors.push_back(kBlack);
    colors.push_back(kRed+1);
    colors.push_back(kBlue+1);
    colors.push_back(kGreen+2);

    std::vector<int> markers;
    markers.push_back(20);
    markers.push_back(21);
    markers.push_back(22);
    markers.push_back(23);

    // ──────────────────────────────────────────
    // 3.  Output file
    // ──────────────────────────────────────────
    TFile *fOut = TFile::Open(fileOut, "RECREATE");

    // ──────────────────────────────────────────
    // 4.  Loop and divide
    //     Bins with fewer than minEntries in either
    //     histogram are zeroed out before dividing
    //     to suppress low-statistics noise.
    //     Each spectrum is normalized per-event before
    //     dividing, so the ratio reflects relative yield,
    //     not just spectrum shape.
    // ──────────────────────────────────────────
    const double minEntries = 1.0;   // tune this threshold as needed

    for (int ir = 0; ir < (int)radii.size(); ++ir) {
        TString R = radii[ir];
        for (int ic = 0; ic < (int)centralities.size(); ++ic) {
            TString cent = centralities[ic];
            for (int ip = 0; ip < (int)ptLeads.size(); ++ip) {
                int ptl = ptLeads[ip];

                TString hName = Form("hSpec_recoPtCorr_ptlead%d_%s_%s",
                                     ptl, R.Data(), cent.Data());
                TString path  = Form("%s/%s/%s", R.Data(), cent.Data(), hName.Data());

                TH1 *hNum = dynamic_cast<TH1*>(fNum->Get(path));
                TH1 *hDen = dynamic_cast<TH1*>(fDen->Get(path));

                if (!hNum) { std::cerr << "Missing in " << labelNum << ": " << path << "\n"; continue; }
                if (!hDen) { std::cerr << "Missing in " << labelDen << ": " << path << "\n"; continue; }

                TString ratioName  = Form("ratio_ptlead%d_%s_%s",
                                          ptl, R.Data(), cent.Data());

                TH1 *hRatioNum = dynamic_cast<TH1*>(hNum->Clone(ratioName));
                TH1 *hRatioDen = dynamic_cast<TH1*>(hDen->Clone("hDenTemp"));

                // Zero out bins below the minimum entries threshold
                // in both histograms to suppress low-statistics bins
                for (int ib = 1; ib <= hRatioNum->GetNbinsX(); ++ib) {
                    bool numLow = hRatioNum->GetBinContent(ib) < minEntries;
                    bool denLow = hRatioDen->GetBinContent(ib) < minEntries;
                    if (numLow || denLow) {
                        hRatioNum->SetBinContent(ib, 0.0);
                        hRatioNum->SetBinError(ib,   0.0);
                        hRatioDen->SetBinContent(ib, 0.0);
                        hRatioDen->SetBinError(ib,   0.0);
                    }
                }

                if (hRatioNum->GetSumw2N() == 0) hRatioNum->Sumw2();
                if (hRatioDen->GetSumw2N() == 0) hRatioDen->Sumw2();

                // Per-event normalization: turn raw counts into a per-event yield
                // BEFORE dividing, so the ratio reflects relative production rate
                // between the two species, not just spectrum shape.
                hRatioNum->Scale(1.0 / nEventsNum);
                hRatioDen->Scale(1.0 / nEventsDen);

                hRatioNum->Divide(hRatioDen);

                // Set axis range to positive pT only
                hRatioNum->GetXaxis()->SetRangeUser(0.0, hRatioNum->GetXaxis()->GetXmax());

                // No title on the histogram itself
                hRatioNum->SetTitle("");

                // Create directories one level at a time
                TString dirR  = Form("%s", R.Data());
                TString dirRC = Form("%s/%s", R.Data(), cent.Data());
                if (!fOut->GetDirectory(dirR))  fOut->mkdir(dirR);
                if (!fOut->GetDirectory(dirRC)) {
                    TDirectory *dR = fOut->GetDirectory(dirR);
                    dR->mkdir(cent);
                }

                fOut->cd(dirRC);
                hRatioNum->Write();

                delete hRatioDen;
                std::cout << "Saved: " << dirRC << "/" << ratioName << "\n";
            }
        }
    }

// ──────────────────────────────────────────
    // 5.  One canvas per (radius, centrality)
    // ──────────────────────────────────────────
    fOut->cd();
    gStyle->SetOptTitle(0);

    for (int ir = 0; ir < (int)radii.size(); ++ir) {
        TString R = radii[ir];
        for (int ic = 0; ic < (int)centralities.size(); ++ic) {
            TString cent = centralities[ic];

            TString cName = Form("c_%s_%s", R.Data(), cent.Data());
            TCanvas *c = new TCanvas(cName, cName, 800, 600);
            c->SetLeftMargin(0.12);
            c->SetBottomMargin(0.12);

            TLegend *leg = new TLegend(0.64, 0.65, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.035);
            // No SetHeader — info box used instead

            TH1  *hLast = NULL;
            bool  first = true;

            for (int ip = 0; ip < (int)ptLeads.size(); ++ip) {
                int ptl = ptLeads[ip];

                TString ratioName = Form("ratio_ptlead%d_%s_%s",
                                         ptl, R.Data(), cent.Data());
                TString path      = Form("%s/%s/%s",
                                         R.Data(), cent.Data(), ratioName.Data());

                TH1 *h = dynamic_cast<TH1*>(fOut->Get(path));
                if (!h) continue;

                h->SetLineColor(colors[ip]);
                h->SetMarkerColor(colors[ip]);
                h->SetMarkerStyle(markers[ip]);
                h->SetMarkerSize(0.8);
                h->SetLineWidth(2);
                h->GetYaxis()->SetRangeUser(0.5, 2.0);
                h->GetXaxis()->SetRangeUser(0.0, h->GetXaxis()->GetXmax());
                h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                h->GetXaxis()->SetTitleSize(0.045);
                h->GetYaxis()->SetTitle(Form("Ratio %s / %s",
                                              labelNum.Data(), labelDen.Data()));
                h->GetYaxis()->SetTitleSize(0.045);
                h->SetTitle("");

                h->Draw(first ? "E1" : "E1 SAME");
                hLast = h;
                first = false;

                leg->AddEntry(h, Form("p_{T}^{lead} > %d GeV/c", ptl), "lp");
            }

            // Reference line at 1
            if (hLast) {
                double xMin = 0.0;
                double xMax = hLast->GetXaxis()->GetXmax();
                TLine *refLine = new TLine(xMin, 1.0, xMax, 1.0);
                refLine->SetLineStyle(2);
                refLine->SetLineColor(kGray+1);
                refLine->SetLineWidth(2);
                refLine->Draw("SAME");
            }

            // ── Info box: system + kinematics ──
            TPaveText *info = new TPaveText(0.64, 0.15, 0.88, 0.30, "NDC");
            info->SetBorderSize(0);
            info->SetFillStyle(0);
            info->SetTextAlign(12);
            info->SetTextSize(0.03);
            info->AddText(Form("%s/%s  at  #sqrt{s_{NN}} = 200 GeV",
                                labelNum.Data(), labelDen.Data()));
            info->AddText(Form("%s,  %s",
                                radLabel(R).Data(), centLabel(cent).Data()));
            info->Draw();

            // ── "THIS THESIS" stamp ──
            TPaveText *thesis = new TPaveText(0.64, 0.30, 0.78, 0.36, "NDC");
            thesis->SetBorderSize(0);
            thesis->SetFillStyle(0);
            thesis->SetTextAlign(12);
            thesis->SetTextSize(0.03);
            thesis->SetTextFont(62);
            thesis->AddText("THIS THESIS");
            thesis->Draw();

            leg->Draw();
            c->Write();

            TString rSafe = R;
            rSafe.ReplaceAll(".", "p");
            TString pdfName = Form("%s/ratio_%s_%s.pdf",
                                    pdfDir.Data(), rSafe.Data(), cent.Data());
            c->SaveAs(pdfName);
            std::cout << "PDF saved: " << pdfName << "\n";

            delete c;
        }
    }

    // ──────────────────────────────────────────
    // 6.  Cleanup
    // ──────────────────────────────────────────
    fOut->Close();
    fNum->Close();
    fDen->Close();

    std::cout << "\nDone! Output written to " << fileOut << "\n";
    std::cout << "PDFs saved in " << pdfDir << "/\n";
}