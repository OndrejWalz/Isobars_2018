// make_qa_pdfs.C
//
// Produces QA PDFs + Jet PDFs.
// Adds a vertical red dashed line at the jet-area cut in:
//   (a) 1D area histograms
//   (b) 2D reco_pt_corr vs reco_area histograms
//
// NOTE (NEF change): hNEF is still produced per (R, centrality) exactly as
// before. In addition, at the end of the run, the three R-histograms for a
// given centrality are combined into a single overlay plot (three labeled,
// differently-colored curves + a legend, NOT summed into one curve).
//
// NOTE (centrality shading): hRefMultCorr now gets shaded regions for the
// 0-10% and 60-80% centrality classes, per the calibration table. The boundaries
// differ slightly between Zr+Zr and Ru+Ru, so the correct column is picked
// automatically (or can be forced) via the new `species` argument.

#include "TFile.h"
#include "TDirectoryFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TObject.h"
#include "TList.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TLatex.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"


#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <cstring>
#include <iterator>
#include <cmath>

using namespace std;

// Convert directory name -> human label for header.
// Return "" if not recognized.
static TString CentralityLabelFromDir(const std::string& cname)
{
  if (cname == "CENT_0_10")  return "0-10 %";
  if (cname == "MID_20_40")  return "20-40 %";
  if (cname == "PERI_60_80") return "60-80 %";
  return "";
}

// -------------------- cuts (VISUAL ONLY) --------------------
const double RECO_DUMMY_CUT = -500.0;
const double pt_reco_corr_cut = 15.0;

// Jet-area cut lines (draw only; do not apply)
const double CUT_AREA_02 = 0.07;  // R = 0.2
const double CUT_AREA_03 = 0.20;  // R = 0.3
const double CUT_AREA_04 = 0.40;  // R = 0.4

// -------------------- centrality class boundaries (refMultCorr >) --------------------
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

// Set at runtime in make_qa_pdfs() based on the chosen/detected species.
// Used by DrawHeader()/DrawHeaderSplit() instead of a hardcoded string.
static TString gSpeciesLabel = "Zr+Zr  (MB)";

// -------------------- utils --------------------
static TString SafeName(const char* name) {
  TString s(name);
  s.ReplaceAll("/", "_");
  s.ReplaceAll(" ", "_");
  return s;
}

static void NormalizeProb(TH1* h) {
  if (!h) return;
  const double I = h->Integral();
  if (I > 0) h->Scale(1.0 / I);
}

static double AreaCutForR(const std::string& rname) {
  if (rname == "R0.2") return CUT_AREA_02;
  if (rname == "R0.3") return CUT_AREA_03;
  if (rname == "R0.4") return CUT_AREA_04;
  return -1.0;
}

// Draw a vertical red dashed line at x = xcut spanning current pad y-range.
static void DrawVLineX(double xcut)
{
  if (!gPad) return;
  if (xcut == 0) return;

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
  ln->SetLineColor(kRed);
  ln->SetLineWidth(3);
  ln->SetLineStyle(2);
  ln->Draw("SAME");
}

// Shades the area under the histogram for specific centrality regions
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
    double yLabel = gPad->GetLogy() ? y2 * 0.15 : y2 * 0.85;


    lx.SetTextColor(regions[i].textColor);
    lx.DrawLatex(xLabel, yLabel, regions[i].label);
  }

  // Redraw the main histogram outline/points so they stay strictly ON TOP of the new shading
  hMain->Draw(Form("%s SAME", drawOpt.Data()));
}

// Draws explicit small vertical tick marks at every bin edge (and the rightmost edge)
static void DrawBinEdgeTicks(TH1* h)
{
  if (!gPad || !h) return;
  gPad->Update();

  TAxis* axis = h->GetXaxis();
  double y0 = gPad->GetUymin();
  double yMax = gPad->GetUymax();
  
  // Calculate a clear tick height (2% of axis span)
  double tickHeight = gPad->GetLogy() ? (yMax / y0) * 0.02 : (yMax - y0) * 0.02;
  double y1 = gPad->GetLogy() ? y0 * 1.05 : y0 + tickHeight;

  // Loop through all bin lower edges, plus the final upper edge
  for (int i = 1; i <= h->GetNbinsX() + 1; ++i) {
    double xEdge = axis->GetBinLowEdge(i);
    TLine* tick = new TLine(xEdge, y0, xEdge, y1);
    tick->SetLineColor(kBlack);
    tick->SetLineWidth(1);
    tick->Draw("SAME");
  }
}


static TH2* TransposeTH2(TH2* h, const char* newName)
{
  if (!h) return 0;

  TAxis* xa = h->GetXaxis();
  TAxis* ya = h->GetYaxis();

  TH2D* ht = new TH2D(newName, "",
                      ya->GetNbins(), ya->GetXmin(), ya->GetXmax(),
                      xa->GetNbins(), xa->GetXmin(), xa->GetXmax());
  ht->Sumw2();

  for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
      // new x = old y, new y = old x
      ht->SetBinContent(iy, ix, h->GetBinContent(ix, iy));
      ht->SetBinError(iy, ix, h->GetBinError(ix, iy));
    }
  }

  return ht;
}

// -------------------- palette (push to yellow/green, avoid red) --------------------
static void SetYellowGreenPalette() {
  const Int_t NRGBs = 4;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = {0.00, 0.35, 0.70, 1.00};

  Double_t red[NRGBs]   = {0.05, 0.00, 0.20, 1.00};
  Double_t green[NRGBs] = {0.10, 0.80, 0.95, 1.00};
  Double_t blue[NRGBs]  = {0.50, 0.75, 0.25, 0.00};

  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

// -------------------- palette (sunset) --------------------
static void SetSunsetPalette() {
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = {0.00, 0.25, 0.50, 0.75, 1.00};

  Double_t red[NRGBs]   = {0.10, 0.40, 0.75, 0.95, 1.00};
  Double_t green[NRGBs] = {0.00, 0.05, 0.20, 0.55, 0.95};
  Double_t blue[NRGBs]  = {0.30, 0.45, 0.35, 0.10, 0.20};

  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}


// -------------------- palette (blue/gold) --------------------
static void SetBlueGoldPalette() {
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = {0.00, 0.25, 0.50, 0.75, 1.00};

  Double_t red[NRGBs]   = {0.02, 0.00, 0.10, 0.65, 1.00};
  Double_t green[NRGBs] = {0.02, 0.35, 0.65, 0.80, 0.90};
  Double_t blue[NRGBs]  = {0.20, 0.55, 0.55, 0.30, 0.10};

  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

// -------------------- style --------------------
static void BasicStyle() {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetTitle(0);

  // remove ticks on top/right
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  gStyle->SetErrorX(0);        // remove horizontal error bars on TH1
  gStyle->SetEndErrorSize(2);  // smaller caps (E1). Still useful elsewhere.
  SetYellowGreenPalette();
}

// -------------------- header text --------------------
enum EHeaderCorner { kTopRight = 0, kTopLeft = 1 };

static void DrawHeader(EHeaderCorner corner = kTopRight,
                       const TString& centLabel = "")
{
  if (!gPad) return;

  const double yTop = 1.0 - gPad->GetTopMargin() - 0.02;
  const double dy   = 0.050;

  double x = 0.0;
  int align = 13; // left-top

  if (corner == kTopRight) {
    x = 1.0 - gPad->GetRightMargin() - 0.02;
    align = 33; // right-top
  } else {
    x = gPad->GetLeftMargin() + 0.02;
    align = 13; // left-top
  }

  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(align);

  tx.SetTextSize(0.040);

  if (centLabel.Length() > 0) {
    tx.DrawLatex(x, yTop,
      Form(" %s %s  #sqrt{#it{s}_{NN}} = 200 GeV", gSpeciesLabel.Data(), centLabel.Data()));
  } else {
    tx.DrawLatex(x, yTop,
      Form(" %s  #sqrt{#it{s}_{NN}} = 200 GeV", gSpeciesLabel.Data()));
  }

  tx.SetTextSize(0.034);
 tx.DrawLatex(x, yTop - dy , "THIS THESIS");
 // tx.DrawLatex(x, yTop - dy, "#it{k}_{T}, #it{R} = 0.3");

}

// -------------------- axis cosmetics helpers --------------------
static void ApplyAxisStyle1D(TH1* h)
{
  if (!h) return;

  const double lab = 0.045;  // label size
  const double tit = 0.050;  // title size

  if (h->GetXaxis()) {
    h->GetXaxis()->SetLabelSize(lab);
    h->GetXaxis()->SetTitleSize(tit);
    h->GetXaxis()->SetTitleOffset(1.10);
  }
  if (h->GetYaxis()) {
    h->GetYaxis()->SetLabelSize(lab);
    h->GetYaxis()->SetTitleSize(tit);
    h->GetYaxis()->SetTitleOffset(1.35);
  }
}

static void ApplyAxisStyle2D(TH2* h)
{
  if (!h) return;

  const double lab = 0.045;
  const double tit = 0.050;

  if (h->GetXaxis()) {
    h->GetXaxis()->SetLabelSize(lab);
    h->GetXaxis()->SetTitleSize(tit);
    h->GetXaxis()->SetTitleOffset(1.10);
  }
  if (h->GetYaxis()) {
    h->GetYaxis()->SetLabelSize(lab);
    h->GetYaxis()->SetTitleSize(tit);
    h->GetYaxis()->SetTitleOffset(1.35);
  }
  if (h->GetZaxis()) {
    h->GetZaxis()->SetLabelSize(lab);
    h->GetZaxis()->SetTitleSize(tit);
    h->GetZaxis()->SetTitleOffset(1.10);
  }
}

static double SmallestPositiveBin(TH1* h)
{
  double minPos = 0.0;
  if (!h) return minPos;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double c = h->GetBinContent(i);
    if (c > 0 && (minPos <= 0 || c < minPos)) minPos = c;
  }
  return minPos;
}

// -------------------- saving helpers --------------------
static void Save1D(TH1* hist,
                   const TString& outPdf,
                   const char* xTitle,
                   const char* yTitle,
                   bool logY,
                   const char* drawOpt,
                   EHeaderCorner headerCorner,
                   double vline_x  = 0.0,
                   const TString& centLabel = "",
                   int cw = 800, int ch = 600,
                   double vline_x2 = 0.0,
                   const CentRegion* centRegions = 0,
                   int nCentRegions = 0,
                   bool drawEdgeTicks = false,
                   bool hideAxisTicks = false)
{
  if (!hist) return;

  TH1* h = (TH1*)hist->Clone(Form("%s_plot", hist->GetName()));
  if (!h) return;
  h->SetDirectory(0);
  h->SetTitle("");

  if (h->GetXaxis()) h->GetXaxis()->SetTitle(xTitle ? xTitle : "");
  if (h->GetYaxis()) h->GetYaxis()->SetTitle(yTitle ? yTitle : "");

  ApplyAxisStyle1D(h);

  TCanvas c(Form("c_%s", h->GetName()), "", cw, ch);

  c.SetLeftMargin(0.15);
  c.SetBottomMargin(0.15);
  c.SetTopMargin(0.06);
  c.SetRightMargin(0.02);

  c.SetTickx(0);
  c.SetTicky(0);

  if (logY) {
    c.SetLogy();
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.2);
    h->SetLineWidth(1);
  }

  if (logY) {
    const double minPos = SmallestPositiveBin(h);
    if (minPos > 0) h->SetMinimum(0.5 * minPos);
  }

  if (hideAxisTicks && h->GetXaxis()) {
    h->GetXaxis()->SetTickLength(0); 
  }

  TString opt = (drawOpt && drawOpt[0]) ? drawOpt : "E1";
  if (logY && opt.Contains("E")) {
    opt = "E0 P";
  }
  h->Draw(opt.Data());

  c.Update();

  if (vline_x  != 0) DrawVLineX(vline_x);
  if (vline_x2 != 0) DrawVLineX(vline_x2);
  
  if (centRegions && nCentRegions > 0) DrawCentralityRegions(h, centRegions, nCentRegions, opt);
  if (drawEdgeTicks) DrawBinEdgeTicks(h);

  DrawHeader(headerCorner, centLabel);

  c.SaveAs(outPdf.Data());
  delete h;
}

static void Save2D(TH2* hist,
                   const TString& outPdf,
                   const char* xTitle,
                   const char* yTitle,
                   bool logZ,
                   EHeaderCorner headerCorner,
                   double vline_x = -1.0,
                   const TString& centLabel = "",
                   int cw = 800, int ch = 600)
{
  if (!hist) return;

  //SetSunsetPalette();

  SetYellowGreenPalette();

  TH2* h = (TH2*)hist->Clone(Form("%s_plot", hist->GetName()));
  if (!h) return;
  h->SetDirectory(0);
  h->SetTitle("");

  if (h->GetXaxis()) h->GetXaxis()->SetTitle(xTitle ? xTitle : "");
  if (h->GetYaxis()) h->GetYaxis()->SetTitle(yTitle ? yTitle : "");

  ApplyAxisStyle2D(h);

  TCanvas c(Form("c_%s", h->GetName()), "", cw, ch);

  c.SetLeftMargin(0.15);
  c.SetBottomMargin(0.15);
  c.SetTopMargin(0.06);
  c.SetRightMargin(0.16);

  c.SetTickx(0);
  c.SetTicky(0);

  if (logZ) c.SetLogz();


  //h->GetXaxis()->SetRangeUser(0.0, 550.0);
  //h->GetYaxis()->SetRangeUser(0.0, 90.0);
  h->Draw("COLZ");
  c.Update();
  if (vline_x > 0) DrawVLineX(vline_x);
  DrawHeader(headerCorner, centLabel);

  c.SaveAs(outPdf.Data());
  delete h;
}

static void DrawHeaderSplit(const TString& centLabel = "")
{
  if (!gPad) return;

  const double yTop = 1.0 - gPad->GetTopMargin() - 0.02;
  const double dy   = 0.055;

  TLatex txL;
  txL.SetNDC();
  txL.SetTextFont(42);
  txL.SetTextAlign(13);
  txL.SetTextSize(0.045);
  txL.DrawLatex(gPad->GetLeftMargin() + 0.02, yTop, "THIS THESIS");
  txL.DrawLatex(gPad->GetLeftMargin() + 0.02, yTop - 1.1*dy, "#it{p}_{T, reco}^{corr} > 15 GeV/#it{c}");


  const double xR = 1.0 - gPad->GetRightMargin() - 0.02;
  TLatex txR;
  txR.SetNDC();
  txR.SetTextFont(42);
  txR.SetTextAlign(33);
  txR.SetTextSize(0.040);

  if (centLabel.Length() > 0) {
    txR.DrawLatex(xR, yTop, Form("%s %s", gSpeciesLabel.Data(), centLabel.Data()));
  } else {
    txR.DrawLatex(xR, yTop, gSpeciesLabel.Data());
  }
  txR.DrawLatex(xR, yTop - dy, "#sqrt{#it{s}_{NN}} = 200 GeV");

}

static void NEFStyleForR(const std::string& rname,
                         Color_t& color, Style_t& marker, TString& label)
{
  if (rname == "R0.2") {
    color = kRed;   marker = 20; label = "R = 0.2";
  } else if (rname == "R0.3") {
    color = kBlue;  marker = 21; label = "R = 0.3";
  } else if (rname == "R0.4") {
    color = kBlack; marker = 22; label = "R = 0.4";
  } else {
    color = kGray+2; marker = 24; label = rname.c_str();
  }
}

static void SaveNEFOverlay(TList* hists,
                           const TString& outPdf,
                           const TString& centLabel = "",
                           int cw = 800, int ch = 600)
{
  if (!hists || hists->GetSize() == 0) return;

  TCanvas c("c_hNEF_overlay", "", cw, ch);
  c.SetLeftMargin(0.15);
  c.SetBottomMargin(0.15);
  c.SetTopMargin(0.06);
  c.SetRightMargin(0.05);
  c.SetTickx(0);
  c.SetTicky(0);
  c.SetLogy();

  TLegend* leg = new TLegend(0.65, 0.2, 0.9, 0.40);
  leg->SetHeader("anti-#it{k}_{T}, Full jets");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.038);

  TObjArray clones;
  clones.SetOwner(kFALSE);

  TIter it(hists);
  TH1* horig = 0;

  while ((horig = (TH1*)it())) {
    std::string rname = horig->GetTitle();

    Color_t col;
    Style_t mk;
    TString label;
    NEFStyleForR(rname, col, mk, label);

    TH1* h = (TH1*)horig->Clone(Form("%s_ov", horig->GetName()));
    h->SetDirectory(0);
    h->SetTitle("");
    h->GetXaxis()->SetTitle("NEF (-)");
    h->GetYaxis()->SetTitle("Probability (a.u.)");
    ApplyAxisStyle1D(h);

    h->SetMarkerColor(col);
    h->SetLineColor(col);
    h->SetMarkerStyle(mk);
    h->SetMarkerSize(1.1);
    h->SetLineWidth(1);

    h->SetMinimum(1e-7);
    h->SetMaximum(1.0);

    leg->AddEntry(h, label.Data(), "P");
    clones.Add(h);
  }

  if (clones.GetEntriesFast() == 0) { delete leg; return; }

  for (int i = 0; i < clones.GetEntriesFast(); ++i) {
    TH1* h = (TH1*)clones.At(i);
    TString drawOpt = (i == 0) ? "P" : "P SAME";
    h->Draw(drawOpt.Data());
  }

  c.Update();
  DrawHeaderSplit(centLabel);
  leg->Draw();

  c.SaveAs(outPdf.Data());

  for (int i = 0; i < clones.GetEntriesFast(); ++i) delete clones.At(i);
  delete leg;
}

void make_qa_pdfs(const char *infile = "data_merged.root", const char* species = 0)
{
  BasicStyle();

  TString spec = species ? TString(species) : TString();
  if (spec.IsNull()) {
    TString fn(infile);
    if (fn.Contains("Ru", TString::kIgnoreCase))      spec = "Ru";
    else if (fn.Contains("Zr", TString::kIgnoreCase)) spec = "Zr";
    else                                               spec = "Zr";
  }
  TString specUp = spec; specUp.ToUpper();
  spec = (specUp == "RU") ? "Ru" : "Zr";

  gSpeciesLabel = Form("%s+%s  (MB)", spec.Data(), spec.Data());
  const CentRegion* centRegionsRefMult = (spec == "Ru") ? kCentRegionsRu : kCentRegionsZr;

  std::cout << "Species: " << spec << "  (header: " << gSpeciesLabel << ")\n";

  gSystem->mkdir("pdf/QA", kTRUE);
  gSystem->mkdir("pdf/Jets", kTRUE);

  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Error: cannot open input file " << infile << "\n";
    return;
  }

  std::cout << "Input : " << infile << "\n";

  {
    TObject* qaObj = fin->Get("QA_histograms");
    TList* qaList = qaObj ? dynamic_cast<TList*>(qaObj) : 0;

    if (!qaList) {
      std::cerr << "QA_histograms not found as TList. Got: "
                << (qaObj ? qaObj->ClassName() : "null") << "\n";
    } else {
      TH1* hcent9   = (TH1*)qaList->FindObject("hcent9");
      TH1* hzvertex = (TH1*)qaList->FindObject("hzvertex");
      TH1* hTrackPt = (TH1*)qaList->FindObject("hTrackPt");
      TH1* hTowE    = (TH1*)qaList->FindObject("hTowE");
      TH1* hDca     = (TH1*)qaList->FindObject("hTrackDca");
      TH1* hEvt0    = (TH1*)qaList->FindObject("hEventStat0");
      TH1* hEvt1    = (TH1*)qaList->FindObject("hEventStat1");
      TH1* hdeltaz_acc = (TH1*)qaList->FindObject("hdeltaz_acc");
      TH1* hdeltaz = (TH1*)qaList->FindObject("hdeltaz");
      TH1* hRefMultCorr = (TH1*)qaList->FindObject("hRefMultCorr");

      TH2* hTowEtaPhi   = (TH2*)qaList->FindObject("hTowEtaPhi");
      TH2* hTrackEtaPhi = (TH2*)qaList->FindObject("hTrackEtaPhi");
      TH2* hRhoVsRefMult = (TH2*)qaList->FindObject("hRhoVsRefMult");
      TH2* hRefMultVsbTofMatch = (TH2*)qaList->FindObject("hRefMultVsbTofMatch");
      TH2* hRefMultVsbTofMatch_ACC = (TH2*)qaList->FindObject("hRefMultVsbTofMatch_ACC");

      if (hcent9) Save1D(hcent9, "pdf/QA/hcent9.pdf",
                         "Centrality", "Counts", false, "HIST", kTopRight);

      if (hEvt0)  Save1D(hEvt0,  "pdf/QA/hEventStat0.pdf",
                         "", "Counts", false, "HIST", kTopRight, 
                         0.0, "", 800, 600, 0.0, 0, 0, true, true);

      if (hEvt1)  Save1D(hEvt1,  "pdf/QA/hEventStat1.pdf",
                         "", "#it{N}_{evt}", false, "HIST", kTopRight,
                         0.0, "", 800, 600, 0.0, 0, 0, true, true);
      if (hTrackPt) Save1D(hTrackPt, "pdf/QA/hTrackPt.pdf",
                           "#it{p}_{T}^{track} (GeV/#it{c})", "#it{N}_{evt}^{weighted}", true, "E1", kTopRight);
      if (hTowE)    Save1D(hTowE, "pdf/QA/hTowE.pdf",
                           "#it{E}_{T} (GeV)", "#it{N}_{evt}^{weighted}", true, "E1", kTopRight, 30.0);
      if (hDca) Save1D(hDca, "pdf/QA/hTrackDca.pdf",
                 "DCA (cm)", "#it{N}_{evt}", true, "E1", kTopRight,
                 1.0);
      if (hdeltaz_acc) Save1D(hdeltaz_acc, "pdf/QA/hdeltaz_acc.pdf",
                                 "#it{V}_{z,TPC} - #it{V}_{z,VPD} (cm)", "#it{N}_{evt}", false, "HIST", kTopRight);
      if (hdeltaz) Save1D(hdeltaz, "pdf/QA/hdeltaz.pdf",
                           "#it{V}_{z,TPC} - #it{V}_{z,VPD} (cm)", "#it{N}_{evt}", false, "HIST", kTopRight,
                           3.0, "", 800, 600, -3.0);
      if (hzvertex) Save1D(hzvertex, "pdf/QA/hzvertex.pdf",
                           "#it{v}_{z} (cm)", "#it{N}_{evt}", false, "E1", kTopRight,
                           -35.0, "", 800, 600, 25.0);

      if (hRefMultCorr) Save1D(hRefMultCorr, "pdf/QA/hRefMultCorr.pdf",
                               "#it{N}_{charged} refMult", "#it{N}_{evt}", true, "HIST", kTopRight,
                               0.0, "", 800, 600, 0.0, centRegionsRefMult, kNumCentRegions);

      if (hTowEtaPhi) {
        TH2* hTowEtaPhiT = (TH2*)TransposeTH2(hTowEtaPhi, "hTowEtaPhi_T");
        Save2D(hTowEtaPhiT, "pdf/QA/hTowEtaPhi.pdf",
               "#varphi_{tow} (rad)", "#eta_{tow} (-)", true, kTopRight);
        delete hTowEtaPhiT;
      }
      if (hTrackEtaPhi) {
        TH2* hTrackEtaPhiT = (TH2*)TransposeTH2(hTrackEtaPhi, "hTrackEtaPhi_T");
        Save2D(hTrackEtaPhiT, "pdf/QA/hTrackEtaPhi.pdf",
               "#varphi_{trk} (rad)", "#eta_{trk} (-)", true, kTopRight);
        delete hTrackEtaPhiT;
      }
      if (hRhoVsRefMult) {
        Save2D(hRhoVsRefMult, "pdf/QA/hRhoVsRefMult.pdf",
               "#it{N}_{charged} refMult", "#rho (GeV/#it{c} per unit area)", true, kTopRight);
        }
      if(hRefMultVsbTofMatch) {
        Save2D(hRefMultVsbTofMatch, "pdf/QA/hRefMultVsbTofMatch.pdf",
               "ToFMatched signal", "#it{N}_{charged} refMult", true, kTopRight);
      }
      if(hRefMultVsbTofMatch_ACC) {
        Save2D(hRefMultVsbTofMatch_ACC, "pdf/QA/hRefMultVsbTofMatch_ACC.pdf",
               "ToFMatched signal", "#it{N}_{charged} refMult", true, kTopRight);
        }
    }
  }

  TList* hNEFGroups = new TList();
  hNEFGroups->SetOwner(kTRUE);

  TIter itR(fin->GetListOfKeys());
  TKey* keyR = 0;

  while ((keyR = (TKey*)itR())) {
    if (std::string(keyR->GetClassName()) != "TDirectoryFile") continue;

    std::string rname = keyR->GetName();
    if (rname.empty() || rname[0] != 'R') continue;

    TDirectoryFile* dirR = (TDirectoryFile*)fin->Get(rname.c_str());
    if (!dirR) continue;

    TString pdfR = Form("pdf/Jets/%s", rname.c_str());
    gSystem->mkdir(pdfR.Data(), kTRUE);

    const double areaCutLine = AreaCutForR(rname);

    TIter itC(dirR->GetListOfKeys());
    TKey* keyC = 0;

    while ((keyC = (TKey*)itC())) {
      if (std::string(keyC->GetClassName()) != "TDirectoryFile") continue;

      std::string cname = keyC->GetName();
      
      if (cname != "CENT_0_10" && cname != "PERI_60_80") continue;

      TDirectoryFile* dirC = (TDirectoryFile*)dirR->Get(cname.c_str());
      if (!dirC) continue;

      TString centLabel = CentralityLabelFromDir(cname);

      TTree* tree = (TTree*)dirC->Get("JetTree");
      if (!tree) continue;

      TString pdfC = Form("%s/%s", pdfR.Data(), cname.c_str());
      gSystem->mkdir(pdfC.Data(), kTRUE);

      TString sel = Form("(reco_pt_corr > %g) && (reco_trigger_match==1)", pt_reco_corr_cut);
      //TString sel = Form("(reco_pt > %g) && (reco_trigger_match==1)", RECO_DUMMY_CUT);
      TString cut = Form("(centralityWeight) * (%s)", sel.Data());

      TString tag = Form("%s_%s", rname.c_str(), cname.c_str());

      TH1D* hAreaJet = new TH1D(Form("hArea_%s", tag.Data()), "", 120, 0.0, 1.0);
      TH1D* hNEFJet  = new TH1D(Form("hNEF_%s",  tag.Data()), "", 120, 0.0, 1.0);
      TH2D* hPtArea  = new TH2D(Form("h2_ptcorr_vs_area_%s", tag.Data()),
                                "", 120, 0.0, 1.0, 140, -20.0, 80.0);

      hAreaJet->Sumw2();
      hNEFJet->Sumw2();
      hPtArea->Sumw2();

      tree->Draw(Form("reco_area>>%s", hAreaJet->GetName()), cut.Data(), "goff");
      tree->Draw(Form("reco_neutral_fraction>>%s", hNEFJet->GetName()), cut.Data(), "goff");
      tree->Draw(Form("reco_pt_corr:reco_area>>%s", hPtArea->GetName()), cut.Data(), "goff");

      NormalizeProb(hNEFJet);

      Save1D(hAreaJet,
             Form("%s/%s.pdf", pdfC.Data(), SafeName(hAreaJet->GetName()).Data()),
             "A (-)", "Weighted entries", true, "E1", kTopRight,
             areaCutLine, centLabel);

      Save1D(hNEFJet,
             Form("%s/%s.pdf", pdfC.Data(), SafeName(hNEFJet->GetName()).Data()),
             "NEF (-)", "Probability (a.u.)", true, "E1", kTopLeft, -1.0, centLabel);

      TList* group = (TList*)hNEFGroups->FindObject(cname.c_str());
      if (!group) {
        group = new TList();
        group->SetName(cname.c_str());
        group->SetOwner(kTRUE);
        hNEFGroups->Add(group);
      }
      TH1D* hNEFforOverlay = (TH1D*)hNEFJet->Clone(Form("hNEF_ov_%s", tag.Data()));
      hNEFforOverlay->SetDirectory(0);
      hNEFforOverlay->SetTitle(rname.c_str());
      group->Add(hNEFforOverlay);

      Save2D(hPtArea,
             Form("%s/%s.pdf", pdfC.Data(), SafeName(hPtArea->GetName()).Data()),
             "A (-)", "#it{p}_{T,reco}^{corr} (GeV/#it{c})", true, kTopRight,
             areaCutLine, centLabel);

      delete hAreaJet;
      delete hNEFJet;
      delete hPtArea;
    }
  }

  {
    TIter itGroup(hNEFGroups);
    TList* group = 0;
    while ((group = (TList*)itGroup())) {
      std::string cname = group->GetName();
      TString centLabel = CentralityLabelFromDir(cname);

      TString outPdf = Form("pdf/Jets/hNEF_allR_%s.pdf", cname.c_str());
      SaveNEFOverlay(group, outPdf, centLabel);
    }
  }

  delete hNEFGroups;

  fin->Close();
  std::cout << "Done. PDFs under pdf/QA, pdf/Jets\n";
}