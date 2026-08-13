// divide_hrefmult.C
//
// Opens two ROOT files, retrieves the "hrefmult" histogram from the
// "QA_histograms" list in each, normalizes each to its own number of
// entries, divides the normalized histogram1 by the normalized histogram2,
// and saves the result to a new ROOT file + PNG/PDF plots.
//
// Note: QA_histograms was written with TObject::kSingleKey, so it is a
// TList (not a TDirectory) even though it displays like a folder in a
// TBrowser -- retrieve it with TFile::Get() + TList::FindObject().
//
// NOTE (centrality shading): the ratio plot now gets the same shaded
// 0-10% and 60-80% centrality regions used for hRefMultCorr in
// make_qa_pdfs.C, per the calibration table. The boundaries differ
// slightly between Zr+Zr and Ru+Ru, so the correct column is picked
// automatically from the input file names (or can be forced) via the
// new `species` argument. Since the ratio is Ru / Zr by convention,
// species detection prefers file1Name (the numerator) and falls back
// to file2Name, defaulting to "Zr" if neither name is informative.
//
// NOTE (low-stats marker): a dashed vertical line + label is drawn at the
// refMult above which per-bin statistics permanently drop below 2% of
// each file's total entries (computed separately for file1 and file2;
// the more conservative/lower of the two refMult values is used, since
// the ratio only makes sense while both inputs have decent statistics).
//
// Usage (from within ROOT or via the bash wrapper):
//   root4star -l -b -q 'divide_hrefmult.C("file1.root", "file2.root")'
//
// Optional third argument sets the output file base name (default: "hrefmult_ratio")
//   root4star -l -b -q 'divide_hrefmult.C("file1.root", "file2.root", "my_ratio")'
//
// Optional fourth argument forces which species' centrality boundaries
// are used for the shaded regions ("Ru" or "Zr"); default is auto-detect.
//   root4star -l -b -q 'divide_hrefmult.C("file1.root", "file2.root", "my_ratio", "Ru")'
//
// Optional fifth argument rebins hrefmult (merging this many original bins
// into one) before normalizing/dividing, e.g. to reduce bin-to-bin noise
// in the ratio; default is 16.
//   root4star -l -b -q 'divide_hrefmult.C("file1.root", "file2.root", "my_ratio", "Ru", 4)'

#include <TFile.h>
#include <TList.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TString.h>
#include <algorithm>
#include <iostream>

// -------------------- centrality class boundaries (refMultCorr >) --------------------
// (ported from make_qa_pdfs.C)
struct CentRegion {
  double min;
  double max;
  const char* label;
  Color_t fillColor;
  Color_t textColor;
};

// 60-80% is bounded by the 80% threshold (7.0/8.0) and the 60% threshold (25.0/26.0)
static const CentRegion kCentRegionsZr[] = {
  {213.0, 9999.0, "0-10%",  kRed-9,   kRed+2},
  {7.0,   25.0,   "60-80%", kAzure-9, kBlue+2}
};

static const CentRegion kCentRegionsRu[] = {
  {216.0, 9999.0, "0-10%",  kRed-9,   kRed+2},
  {8.0,   26.0,   "60-80%", kAzure-9, kBlue+2}
};

static const int kNumCentRegions = 2;

// Bins whose (normalized) content falls below this fraction of the total
// entries are considered "low statistics". Default: 2%.
static const double LOW_STATS_FRACTION = 0.01;

// Scan a normalized histogram (bin content = fraction of total entries) from
// high x to low x and find the upper edge of the last (highest-x) bin whose
// content is still >= fracThreshold. Everything above that edge is treated
// as the low-statistics tail. This is robust against the low-refMult bins
// near zero, which can also dip below threshold but aren't part of the tail.
static double FindLowStatsEdge(TH1* hNorm, double fracThreshold = LOW_STATS_FRACTION)
{
  if (!hNorm) return 0.0;

  const int nb = hNorm->GetNbinsX();
  int lastGoodBin = 0;
  for (int i = 1; i <= nb; ++i) {
    if (hNorm->GetBinContent(i) >= fracThreshold) lastGoodBin = i;
  }

  if (lastGoodBin == 0)  return hNorm->GetXaxis()->GetXmin(); // never reached threshold
  if (lastGoodBin == nb) return hNorm->GetXaxis()->GetXmax(); // stays above threshold throughout

  return hNorm->GetXaxis()->GetBinUpEdge(lastGoodBin);
}

// Draw a vertical dashed line at x = xcut spanning the current pad's
// y-range, with a short text label near the top of the pad.
static void DrawVLineX(double xcut, const char* label = 0,
                        Color_t lineColor = kGray+2, Style_t lineStyle = 3)
{
  if (!gPad) return;

  gPad->Update();

  double y1, y2;
  if (gPad->GetLogy()) {
    y1 = pow(10.0, gPad->GetUymin());
    y2 = pow(10.0, gPad->GetUymax());
  } else {
    y1 = gPad->GetUymin();
    y2 = gPad->GetUymax();
  }

  TLine* ln = new TLine(xcut, y1, xcut, y2);
  ln->SetLineColor(lineColor);
  ln->SetLineWidth(2);
  ln->SetLineStyle(lineStyle);
  ln->Draw("SAME");

  if (label && label[0]) {
    TLatex lx;
    lx.SetTextFont(42);
    lx.SetTextSize(0.032);
    lx.SetTextColor(lineColor);
    lx.SetTextAlign(11); // left-bottom, so text sits just right of the line
    double yLabel = gPad->GetLogy() ? y2 * 0.5 : y1 + 0.92 * (y2 - y1);
    lx.DrawLatex(xcut + 0.01 * (gPad->GetUxmax() - gPad->GetUxmin()), yLabel, label);
  }
}

// Shades the area under the histogram for specific centrality regions
// (ported from make_qa_pdfs.C, unchanged in behavior)
static void DrawCentralityRegions(TH1* hMain, const CentRegion* regions, int nRegions, const TString& drawOpt)
{
  if (!gPad || !hMain || !regions || nRegions <= 0) return;

  gPad->Update();
  double y2 = gPad->GetLogy() ? pow(10.0, gPad->GetUymax()) : gPad->GetUymax();

  TLatex lx;
  lx.SetTextFont(42);
  lx.SetTextSize(0.035);
  lx.SetTextAlign(22); // center horizontally & vertically

  for (int i = 0; i < nRegions; ++i) {
    // Clone to make a shaded copy for just this region
    TH1* hShade = (TH1*)hMain->Clone(Form("%s_shade_%d", hMain->GetName(), i));

    for (int b = 1; b <= hShade->GetNbinsX(); ++b) {
      double xC = hShade->GetXaxis()->GetBinCenter(b);
      // Zero out bins outside this region
      if (xC < regions[i].min || xC >= regions[i].max) {
        hShade->SetBinContent(b, 0);
        hShade->SetBinError(b, 0);
      }
    }

    hShade->SetFillColor(regions[i].fillColor);
    hShade->SetFillStyle(1001); // solid fill
    hShade->SetLineColor(regions[i].fillColor); // hide outline
    hShade->SetLineWidth(1);
    hShade->Draw("HIST SAME"); // Draw shading

    // Draw label above the shaded region
    double xLabel;
    if (regions[i].max > 1000) {
      xLabel = regions[i].min + 60.0;
    } else {
      // Shift 60-80% more to the right (closer to max, or add a fixed offset)
      if (TString(regions[i].label).Contains("60-80")) {
        xLabel = regions[i].max + 45.0; // Pushes it near the right edge of the 60-80% band
      } else {
        xLabel = (regions[i].min + regions[i].max) / 2.0;
      }
    }
    double yLabel = gPad->GetLogy() ? y2 * 0.15 : y2 * 0.82;

    lx.SetTextColor(regions[i].textColor);
    lx.DrawLatex(xLabel, yLabel, regions[i].label);
  }

  // Redraw the main histogram outline/points so they stay strictly ON TOP of the new shading
  hMain->Draw(Form("%s SAME", drawOpt.Data()));
}

// Auto-detect species ("Ru" or "Zr") from a file name; returns "" if not recognized.
static TString SpeciesFromFileName(const char* fname)
{
  TString fn(fname);
  if (fn.Contains("Ru", TString::kIgnoreCase)) return "Ru";
  if (fn.Contains("Zr", TString::kIgnoreCase)) return "Zr";
  return "";
}

void divide_hrefmult(const char *file1Name, const char *file2Name,
                      const char *outBaseName = "hRefMultCorr_ratio",
                      const char *species = 0,
                      Int_t nRebin = 42) {

  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);

  // -- Figure out which centrality-region boundaries to shade with.
  //    Prefer an explicit override, then auto-detect from file1 (numerator),
  //    then file2 (denominator), defaulting to Zr if nothing is recognized.
  TString spec = species ? TString(species) : TString();
  if (spec.IsNull()) spec = SpeciesFromFileName(file1Name);
  if (spec.IsNull()) spec = SpeciesFromFileName(file2Name);
  if (spec.IsNull()) spec = "Zr";
  TString specUp = spec; specUp.ToUpper();
  spec = (specUp == "RU") ? "Ru" : "Zr";

  const CentRegion* centRegions = (spec == "Ru") ? kCentRegionsRu : kCentRegionsZr;
  std::cout << "Using " << spec << "+" << spec << " centrality boundaries for shading" << std::endl;

  TFile *f1 = TFile::Open(file1Name, "READ");
  TFile *f2 = TFile::Open(file2Name, "READ");

  if (!f1 || f1->IsZombie()) {
    std::cerr << "ERROR: could not open file " << file1Name << std::endl;
    return;
  }
  if (!f2 || f2->IsZombie()) {
    std::cerr << "ERROR: could not open file " << file2Name << std::endl;
    return;
  }

  TList *list1 = (TList *)f1->Get("QA_histograms");
  TList *list2 = (TList *)f2->Get("QA_histograms");

  if (!list1) {
    std::cerr << "ERROR: 'QA_histograms' list not found in " << file1Name << std::endl;
    return;
  }
  if (!list2) {
    std::cerr << "ERROR: 'QA_histograms' list not found in " << file2Name << std::endl;
    return;
  }

  TH1 *h1 = (TH1 *)list1->FindObject("hRefMultCorr");
  TH1 *h2 = (TH1 *)list2->FindObject("hRefMultCorr");

  if (!h1) {
    std::cerr << "ERROR: 'hRefMultCorr' not found in " << file1Name << std::endl;
    return;
  }
  if (!h2) {
    std::cerr << "ERROR: 'hRefMultCorr' not found in " << file2Name << std::endl;
    return;
  }

  Double_t n1 = h1->GetEntries();
  Double_t n2 = h2->GetEntries();

  std::cout << "File 1: " << file1Name << "  -> hRefMultCorr entries: " << n1 << std::endl;
  std::cout << "File 2: " << file2Name << "  -> hRefMultCorr entries: " << n2 << std::endl;

  if (n1 <= 0) {
    std::cerr << "ERROR: hRefMultCorr in " << file1Name << " has zero entries, cannot normalize" << std::endl;
    return;
  }
  if (n2 <= 0) {
    std::cerr << "ERROR: hRefMultCorr in " << file2Name << " has zero entries, cannot normalize" << std::endl;
    return;
  }

  // -- Build TH1D copies of h1/h2 explicitly before normalizing. h1/h2 are
  //    TH1I (integer bin content); calling Scale() directly on a TH1I
  //    SILENTLY TRUNCATES the normalized (fractional) values to integers,
  //    e.g. a bin with content 3 scaled by 1/1000 becomes 0 instead of
  //    0.003. Copying into TH1D first preserves the true floating-point
  //    normalized values.
  TH1D *h1n = new TH1D("hRefMultCorr_file1_norm", "hRefMultCorr file1 (normalized);refMult;fraction of entries",
                        h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
  TH1D *h2n = new TH1D("hRefMultCorr_file2_norm", "hRefMultCorr file2 (normalized);refMult;fraction of entries",
                        h2->GetNbinsX(), h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());
  h1n->SetDirectory(0);
  h2n->SetDirectory(0);
  h1n->Sumw2();
  h2n->Sumw2();

  // -- Copy bin contents/errors (including under/overflow) before scaling.
  //    If Sumw2 was never called on h1/h2, GetBinError() falls back to
  //    sqrt(content), which is the correct counting-statistics error here.
  for (Int_t i = 0; i <= h1->GetNbinsX() + 1; i++) {
    h1n->SetBinContent(i, h1->GetBinContent(i));
    h1n->SetBinError(i, h1->GetBinError(i));
  }
  for (Int_t i = 0; i <= h2->GetNbinsX() + 1; i++) {
    h2n->SetBinContent(i, h2->GetBinContent(i));
    h2n->SetBinError(i, h2->GetBinError(i));
  }

  // -- Optionally rebin (merge neighboring bins) before normalizing. Doing
  //    this before Scale() means Rebin() sums the still-integer-valued
  //    contents/errors correctly; scaling afterwards then normalizes the
  //    coarser bins. Rebin() preserves the total integral, so normalizing
  //    by n1/n2 (the unrebinned entry counts) is still correct. If nRebin
  //    doesn't evenly divide the bin count, ROOT folds the remainder into
  //    the last bin and prints a warning.
  if (nRebin > 1) {
    h1n->Rebin(nRebin);
    h2n->Rebin(nRebin);
    std::cout << "Rebinning hrefmult by a factor of " << nRebin
              << " (" << h1n->GetNbinsX() << " bins after rebin)" << std::endl;
  }

  h1n->Scale(1.0 / n1);
  h2n->Scale(1.0 / n2);

  // -- Build the ratio histogram as a TH1D explicitly, same reasoning as above.
  //    Use h1n's (possibly rebinned) axis so it matches h1n/h2n exactly.
  TH1D *hRatio = new TH1D("hRefMultCorr_ratio", "hRefMultCorr ratio, normalized (Ru / Zr);refMult;ratio",
                           h1n->GetNbinsX(), h1n->GetXaxis()->GetXmin(), h1n->GetXaxis()->GetXmax());
  hRatio->SetDirectory(0);
  hRatio->Sumw2();

  hRatio->Divide(h1n, h2n);

  // -- Save ratio histogram to a new ROOT file
  TString outRootName = TString::Format("%s.root", outBaseName);
  TFile *fOut = new TFile(outRootName, "RECREATE");
  hRatio->Write();
  h1n->Write();
  h2n->Write();
  fOut->Close();
  std::cout << "Saved ratio histogram to: " << outRootName << std::endl;




  // -- Save PNG and PDF plots
  TCanvas *cRatio = new TCanvas("cRatio", "hrefmult ratio", 800, 600);
  hRatio->SetTitle("");
  hRatio->SetStats(0);
  hRatio->GetYaxis()->SetRangeUser(0.8, 1.31);
  hRatio->GetXaxis()->SetTitle("#it{N}_{charged} refMult");
  hRatio->GetXaxis()->SetTitleSize(0.045);
  hRatio->GetYaxis()->SetTitleSize(0.045);
  hRatio->GetYaxis()->SetTitle("[1/#it{N}_{evt} Counts]_{Ru+Ru} / [1/#it{N}_{evt} Counts]_{Zr+Zr}");
  TString ratioDrawOpt = "E";
  hRatio->Draw(ratioDrawOpt.Data());

  cRatio->Update();
  DrawCentralityRegions(hRatio, centRegions, kNumCentRegions, ratioDrawOpt);


  // -- Guide horizontal line at y = 1 ──
  TLine *line1 = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
  line1->SetLineColor(kBlack);
  line1->SetLineStyle(2); // dashed
  line1->SetLineWidth(1);
  line1->Draw("SAME");

  TPaveText *thesis = new TPaveText(0.11, 0.83, 0.40, 0.88, "NDC");
  thesis->SetBorderSize(0);
  thesis->SetFillStyle(0);
  thesis->SetTextAlign(12);
  thesis->SetTextSize(0.035);
  thesis->SetTextFont(42); // Bold font
  thesis->AddText("THIS THESIS");
  thesis->Draw();

  TPaveText *info = new TPaveText(0.1, 0.72, 0.55, 0.82, "NDC");
  info->SetBorderSize(0);
  info->SetFillStyle(0);
  info->SetTextAlign(12); // Left-aligned, vertically centered
  info->SetTextSize(0.029);
  info->SetTextFont(42);
  info->AddText(" Ru+Ru/Zr+Zr (HT2) #sqrt{s_{NN}} = 200 GeV");
  info->AddText(" #it{E}_{T}^{trig} > 4.2 GeV");
  info->Draw();
  // -- Low-stats marker: find, per file, the refMult above which that
  //    file's bins permanently fall below LOW_STATS_FRACTION of its total
  //    entries, then mark the more conservative (lower) of the two, since
  //    the ratio is only meaningful while *both* inputs have decent stats.
  double lowStatsX1 = FindLowStatsEdge(h1n);
  double lowStatsX2 = FindLowStatsEdge(h2n);
  double lowStatsX  = std::min(lowStatsX1, lowStatsX2);

  std::cout << "Low-stats edge (< " << (LOW_STATS_FRACTION * 100.0) << "% of total entries per bin): "
            << "File 1 -> refMult = " << lowStatsX1 << ",  "
            << "File 2 -> refMult = " << lowStatsX2 << ",  "
            << "using refMult = " << lowStatsX << std::endl;

  DrawVLineX(lowStatsX, Form("->low stat. (<%.0f%% of all events)", LOW_STATS_FRACTION * 100.0));

  TString outPngName = TString::Format("%s.png", outBaseName);
  TString outPdfName = TString::Format("%s.pdf", outBaseName);
  //cRatio->SaveAs(outPngName);
  cRatio->SaveAs(outPdfName);
  std::cout << "Saved plots to: " << outPngName << " and " << outPdfName << std::endl;

  f1->Close();
  f2->Close();
}