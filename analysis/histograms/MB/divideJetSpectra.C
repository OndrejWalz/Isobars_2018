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
#include <TObjArray.h>
#include <TObjString.h>

// Helper: convert internal centrality string to a nice label
TString centLabel(TString cent) {
    if (cent == "CENT_0_10")  return "0-10%";
    if (cent == "MID_20_40")  return "20-40%";
    if (cent == "PERI_60_80") return "60-80%";
    return cent;
}

// Helper: convert internal radius string to a nice label
TString radLabel(TString R) {
    if (R == "R0.2") return "#it{R} = 0.2";
    if (R == "R0.3") return "#it{R} = 0.3";
    if (R == "R0.4") return "#it{R} = 0.4";
    return R;
}

// Helper: read total number of events from data_merged_<species>.root
// Structure: file -> "QA_histograms"
// (TList) -> "hevents" (TH1)
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

    TH1 *hEventStat1 = dynamic_cast<TH1*>(qaList->FindObject("hEventStat1"));
    if (!hEventStat1) {
        std::cerr << "Cannot find 'hEventStat1' in QA_histograms in " << fileName << "\n";
        f->Close();
        return -1.0;
    }

    double nEvents = hEventStat1->GetBinContent(7);
    std::cout << "  " << fileName << " -> nEvents (hEventStat1 bin 7) = " << nEvents << "\n";

    f->Close();
    return nEvents;
}

// Helper: parse a comma-separated string of bin edges (e.g. "0,1,2,5,10")
// into a vector<double>. Returns an empty vector if the string is empty
// or malformed (fewer than 2 edges).
std::vector<double> parseBinEdges(TString edgesStr) {
    std::vector<double> edges;
    if (edgesStr.Length() == 0) return edges;

    TObjArray *toks = edgesStr.Tokenize(",");
    for (int i = 0; i < toks->GetEntries(); ++i) {
        TString tok = ((TObjString*)toks->At(i))->GetString();
        tok.ReplaceAll(" ", "");
        if (tok.Length() == 0) continue;
        edges.push_back(tok.Atof());
    }
    delete toks;

    if (edges.size() < 2) {
        std::cerr << "WARNING: customBinEdges parsed to fewer than 2 edges, ignoring: \""
                   << edgesStr << "\"\n";
        edges.clear();
    }
    return edges;
}

void divideJetSpectra(TString fileNum       = "hists_data_Ru.root",
                      TString fileDen        = "hists_data_Zr.root",
                      TString fileOut        = "ratio_Ru_over_Zr.root",
                      TString labelNum       = "Ru",
                      TString labelDen       = "Zr",
                      TString pdfDir         = "ratio_Ru_over_Zr",
                      TString fileEventsNum  = "../../../trees/data_merged_Ru.root",
                      TString fileEventsDen  = "../../../trees/data_merged_Zr.root",
                      // Comma-separated list of new bin edges, e.g.
                      // "0,1,2,3,5,7,10,15,20". Leave empty ("") for no
                      // rebinning. Each edge MUST coincide exactly with an
                      // existing bin boundary in the original histograms --
                      // ROOT's variable-width Rebin() can only merge whole
                      // original bins together, it cannot split one.
                      TString customBinEdges = "") {

    // -- CRITICAL FIX: without this, every Clone()/new TH1 silently
    //    auto-registers itself into whatever TDirectory is "current"
    //    (gDirectory) at the moment it's created. Below, hRatioDen is
    //    Clone()'d while fOut (or one of its subdirectories, since the
    //    loop calls fOut->cd(dirRC) each pass and never cd's back) is the
    //    current directory, and is then manually delete'd -- leaving that
    //    directory holding a DANGLING POINTER to freed memory. When
    //    fOut->Close() walks its object list, this can double-free /
    //    crash / silently corrupt the output. AddDirectory(kFALSE) makes
    //    all subsequent histogram objects directory-agnostic in memory,
    //    so manual delete is safe and Write() is the only way they get
    //    attached to a file.
    TH1::AddDirectory(kFALSE);

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
    std::cout << "PDF folder  : " << pdfDir   << "/\n";

    std::vector<double> binEdges = parseBinEdges(customBinEdges);
    if (!binEdges.empty()) {
        std::cout << "Rebinning   : custom edges (" << binEdges.size() - 1 << " bins): ";
        for (size_t i = 0; i < binEdges.size(); ++i) {
            std::cout << binEdges[i] << (i + 1 < binEdges.size() ? ", " : "\n");
        }
    } else {
        std::cout << "Rebinning   : none\n";
    }
    std::cout << "\n";

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
    // Ncoll normalization values
    // ──────────────────────────────────────────


    // <Ncoll> per centrality, same order as `centralities` above.
    // Set both to empty {} to disable Ncoll normalization entirely.
    std::vector<double> ncollNum; // Ru+Ru
    ncollNum.push_back(356.0);   // 0-10%
    ncollNum.push_back(117.7);   // 20-40%
    ncollNum.push_back(11.85);    // 60-80%

    std::vector<double> ncollDen; // Zr+Zr
    ncollDen.push_back(351.5);
    ncollDen.push_back(111.2);
    ncollDen.push_back(10.71);

    bool doNcollNorm = !ncollNum.empty() && !ncollDen.empty();

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

                // -- NOTE: if hNum/hDen turn out to be integer-typed
                //    (TH1I/TH1C/TH1S) rather than TH1F/TH1D, Clone()
                //    preserves that integer type, and the Scale() call
                //    below (per-event normalization) will silently
                //    truncate every bin to 0. Check hNum->ClassName();
                //    if it's an integer type, clone into an explicit
                //    TH1D the way divide_hrefmult.C does, instead of
                //    using Clone() directly.
                TString ratioName  = Form("ratio_ptlead%d_%s_%s",
                                          ptl, R.Data(), cent.Data());

                TH1 *hRatioNum = dynamic_cast<TH1*>(hNum->Clone(ratioName));
                TH1 *hRatioDen = dynamic_cast<TH1*>(hDen->Clone("hDenTemp"));

                // Apply custom rebinning, if requested, BEFORE the low-stat
                // zeroing/threshold check and BEFORE per-event Scale() --
                // that way the stat threshold applies to the already-merged
                // (higher-count) bins, which is what you want.
                if (!binEdges.empty()) {
                    // TH1::Rebin(nbins, newname, xbins) rebins the histogram
                    // in place and returns the same pointer when newname=""
                    // is not used this way -- here we DO pass a name because
                    // Rebin's variable-width overload requires one; passing
                    // the same name as the original just replaces it in
                    // memory (this is a local TH1*, not registered in any
                    // TDirectory, since AddDirectory(kFALSE) is set above).
                    TH1 *hRatioNumRebinned = hRatioNum->Rebin(
                        (Int_t)binEdges.size() - 1, ratioName, &binEdges[0]);
                    TH1 *hRatioDenRebinned = hRatioDen->Rebin(
                        (Int_t)binEdges.size() - 1, "hDenTemp", &binEdges[0]);

                    if (!hRatioNumRebinned || !hRatioDenRebinned) {
                        std::cerr << "ERROR: Rebin failed for " << path
                                  << " -- customBinEdges likely don't align "
                                  << "with this histogram's existing bin edges. "
                                  << "Skipping.\n";
                        delete hRatioNum;
                        delete hRatioDen;
                        continue;
                    }
                    // NOTE: ROOT's Rebin() doesn't always return null on
                    // misaligned edges -- it can also just print a warning
                    // to stderr (e.g. "Bin edges array not sorted" or
                    // "Bins in the new histogram are not defined in terms
                    // of the old ones") and still return a histogram, which
                    // may not be what you intended. Watch stderr on the
                    // first run with a new customBinEdges list, and spot
                    // check the resulting bin edges against the originals.

                    // Rebin() with a new name returns a NEW histogram object;
                    // the originals still need to be freed separately.
                    delete hRatioNum;
                    delete hRatioDen;
                    hRatioNum = hRatioNumRebinned;
                    hRatioDen = hRatioDenRebinned;
                }

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
                          // Per-event (and optionally per-Ncoll) normalization, applied
                // BEFORE dividing, so the ratio reflects relative yield.
                double scaleNum = 1.0 / nEventsNum;
                double scaleDen = 1.0 / nEventsDen;
                if (doNcollNorm) {
                    scaleNum /= ncollNum[ic];
                    scaleDen /= ncollDen[ic];
                }
                hRatioNum->Scale(scaleNum);
                hRatioDen->Scale(scaleDen);

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

                // Safe now: hRatioNum/hRatioDen were never auto-attached
                // to fOut or any subdirectory (AddDirectory(kFALSE) above),
                // so deleting them here can't leave a dangling pointer
                // behind in any TDirectory's object list.
                delete hRatioNum;
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
            std::vector<TH1*> loadedHists; // track so we can free after drawing

            for (int ip = 0; ip < (int)ptLeads.size(); ++ip) {
                int ptl = ptLeads[ip];

                TString ratioName = Form("ratio_ptlead%d_%s_%s",
                                         ptl, R.Data(), cent.Data());
                TString path      = Form("%s/%s/%s",
                                         R.Data(), cent.Data(), ratioName.Data());

                TH1 *h = dynamic_cast<TH1*>(fOut->Get(path));
                if (!h) continue;
                loadedHists.push_back(h);

                h->SetLineColor(colors[ip]);
                h->SetMarkerColor(colors[ip]);
                h->SetMarkerStyle(markers[ip]);
                h->SetMarkerSize(0.8);
                h->SetLineWidth(2);
                h->GetYaxis()->SetRangeUser(0.7, 1.5);
                h->GetXaxis()->SetRangeUser(0.0, h->GetXaxis()->GetXmax());
                h->GetXaxis()->SetTitle("#it{p}_{T, reco}^{corr} (GeV/#it{c})");
                h->GetXaxis()->SetTitleSize(0.045);
                h->GetYaxis()->SetTitle(doNcollNorm
                    ? "[1/(N_{evt}#timesN_{coll}) dN/d#it{p}_{T}]_{Ru+Ru} / [1/(N_{evt}#timesN_{coll}) dN/d#it{p}_{T}]_{Zr+Zr}"
                    : "[1/N_{evt} dN/d#it{p}_{T}]_{Ru+Ru} / [1/N_{evt} dN/d#it{p}_{T}]_{Zr+Zr}");
                h->GetYaxis()->SetTitleSize(0.045);
                h->SetTitle("");

                h->Draw(first ? "E1" : "E1 SAME");
                hLast = h;
                first = false;

                leg->AddEntry(h, Form("p_{T}^{lead} > %d GeV/c", ptl), "lp");
            }

            // Reference line at 1
            TLine *refLine = NULL;
            if (hLast) {
                double xMin = 0.0;
                double xMax = hLast->GetXaxis()->GetXmax();
                refLine = new TLine(xMin, 1.0, xMax, 1.0);
                refLine->SetLineStyle(2);
                refLine->SetLineColor(kGray+1);
                refLine->SetLineWidth(2);
                refLine->Draw("SAME");
            }

            // ── Info box: system + kinematics ──

            if(centLabel(cent)=="0-10%") {
                TPaveText *info = new TPaveText(0.16, 0.72, 0.55, 0.82, "NDC");
                info->SetBorderSize(0);
                info->SetFillStyle(0);
                info->SetTextAlign(12); // Left-aligned, vertically centered
                info->SetTextSize(0.03);
                info->SetTextFont(42);
                info->AddText("STAR Ru+Ru/Zr+Zr (MB) #sqrt{s_{NN}} = 200 GeV");
                info->AddText(Form("%s,  %s", radLabel(R).Data(), centLabel(cent).Data()));
                info->Draw();
             
                TPaveText *thesis = new TPaveText(0.16, 0.82, 0.40, 0.88, "NDC");
                thesis->SetBorderSize(0);
                thesis->SetFillStyle(0);
                thesis->SetTextAlign(12);
                thesis->SetTextSize(0.035);
                thesis->SetTextFont(42); // Bold font for a stamp
                thesis->AddText("THIS THESIS");
                thesis->Draw();
       
            }

            if (centLabel(cent)=="60-80%") {
                TPaveText *thesis = new TPaveText(0.56, 0.3, 0.5, 0.3, "NDC");
                thesis->SetBorderSize(0);
                thesis->SetFillStyle(0);
                thesis->SetTextAlign(12);
                thesis->SetTextSize(0.035);
                thesis->SetTextFont(42); // Bold font for a stamp
                thesis->AddText("THIS THESIS");
                thesis->Draw();


                TPaveText *info = new TPaveText(0.5, 0.15, 0.55, 0.25, "NDC");
                info->SetBorderSize(0);
                info->SetFillStyle(0);
                info->SetTextAlign(12); // Left-aligned, vertically centered
                info->SetTextSize(0.03);
                info->SetTextFont(42);
                info->AddText("STAR Ru+Ru/Zr+Zr (MB) #sqrt{s_{NN}} = 200 GeV");
                info->AddText(Form("%s,  %s", radLabel(R).Data(), centLabel(cent).Data()));
                info->Draw();
             
 
            }

                leg->Draw();
                c->Write();

                TString rSafe = R;
                rSafe.ReplaceAll(".", "p");
                TString pdfName = Form("%s/ratio_%s_%s.pdf",
                                        pdfDir.Data(), rSafe.Data(), cent.Data());
                c->SaveAs(pdfName);
                std::cout << "PDF saved: " << pdfName << "\n";

                // Only delete the canvas itself — TPad::~TPad deletes every primitive
                // that was Draw()'n onto it (leg, info, thesis, refLine, and the
                // histograms in loadedHists), so deleting them individually here would
                // be a double free.
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