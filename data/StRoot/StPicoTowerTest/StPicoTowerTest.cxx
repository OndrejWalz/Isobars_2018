#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TLorentzVector.h"

#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoPrescales/StPicoPrescales.h"
#include "StPicoTowerTest.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoBTowHit.h"
#include "StCuts.h"

#include "TSystem.h"
#include "TEnv.h"

#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "BemcNewCalib.h"
//////Refit include lib
#include "PhysicalConstants.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorD.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TFile.h"
#include "TDirectory.h"
#include "StEvent/StDcaGeometry.h"
//
#include <bitset>

#include <set>
#include <vector>
#include <stdio.h>
#include <time.h>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>
using namespace std;

//-------------------------------------
// >>> EDIT THIS: absolute path of the run list (one run number per line) <<<
// Used when $RUNLIST is not set. Find it with:  find /pwg/walz -iname 'isobars_2018.txt'
static const char* const kRunListPath = "/gpfs/mnt/gpfs01/star/pwg/walz/Isobars2018/data/Isobars_2018.txt";
//-------------------------------------

ClassImp(StPicoTowerTest)

  StPicoTowerTest::StPicoTowerTest(char const * name, char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil,int pYear):
    StMaker(name),mPicoDstMaker(picoDstMaker), mGRefMultCorrUtil(grefmultCorrUtil),
    mOutFileName(outName),mOutputFile(NULL), mChain(NULL), mEventCounter(0),mYear(pYear){}

Int_t StPicoTowerTest::Init()
{
  mEmcPosition = new StEmcPosition();

  mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
  mOutputFile->cd();

  // ---- run list: $RUNLIST if set, otherwise try the usual locations ----
  {
    // Priority: $RUNLIST (override) > hardcoded kRunListPath
    std::string runListPath = gSystem->Getenv("RUNLIST") ? gSystem->Getenv("RUNLIST") : kRunListPath;
    std::ifstream runs(runListPath.c_str());
    if (runListPath.empty() || !runs.is_open()) {
      LOG_ERROR << "StPicoTowerTest - Cannot find run list file." << endm;
      LOG_ERROR << "  working directory : " << gSystem->WorkingDirectory() << endm;
      LOG_ERROR << "  tried             : " << runListPath << endm;
      LOG_ERROR << "  Fix: edit kRunListPath in StPicoTowerTest.cxx (and rebuild) or export RUNLIST=..." << endm;
      return kStFatal;   // stop now instead of processing millions of events for nothing
    }
    unsigned int run;
    int idx = 1;
    while (runs >> run) {
      mRunIndexMap[run] = idx++;
    }
    runs.close();
    LOG_INFO << "StPicoTowerTest - run list " << runListPath << " : " << mRunIndexMap.size() << " runs" << endm;
  }
  int nRuns = mRunIndexMap.size();
  if (nRuns == 0) {
    LOG_ERROR << "StPicoTowerTest - run list is empty!" << endm;
    return kStFatal;
  }

  MapTower_RunID_NHits_TowerID = new TH2D("MapTower_RunID_NHits_TowerID","MapTower_RunID_NHits_TowerID;runIndex;towerID", nRuns+1,0,nRuns+1,4800,0,4800);
  MapTower_RunID_Energy_TowerID = new TH2D("MapTower_RunID_Energy_TowerID","MapTower_RunID_Energy_TowerID_ET;runIndex;towerID", nRuns+1,0,nRuns+1,4800,0,4800);
  MapTower_RunID_Energy_TowerID_energy = new TH2D("MapTower_RunID_Energy_TowerID_energy","MapTower_RunID_Energy_TowerID_energy;runIndex;towerID", nRuns+1,0,nRuns+1,4800,0,4800);
  MapTower_RunID_NEvents = new TH1D("MapTower_RunID_NEvents", "MapTower_RunID_NEvents", nRuns+1,0,nRuns+1);
  vtxr = new TH2D("vtxr",";PVtx.x() [cm]; PVtx.y() [cm]",100,-3,3,100,-3,3);

  hRunPxl1 = new TH2D("hRunPxl1", "hRunPxl1", nRuns+1,0,nRuns+1,1500,0,1500);
  hRunPxl2 = new TH2D("hRunPxl2", "hRunPxl2", nRuns+1,0,nRuns+1,1500,0,1500);
  hRunIst = new TH2D("hRunIst", "hRunIst", nRuns+1,0,nRuns+1,1500,0,1500);
  hRunSsd = new TH2D("hRunSsd", "hRunSsd", nRuns+1,0,nRuns+1,1500,0,1500);

  mOutputFile->cd();

  StMaker* maker = GetMaker("Eread");
  cout << (maker ? "Maker is valid!" : "Maker is invalid!") << endl;
  mADCtoEMaker = dynamic_cast<StEmcADCtoEMaker*>(maker);
  cout << (mADCtoEMaker ? "mADCtoEMaker is valid!" : "mADCtoEMaker is invalid!") << endl;
  if (!mADCtoEMaker) {
    LOG_ERROR << "StPicoTowerTest - no StEmcADCtoEMaker named Eread" << endm;
    return kStFatal;
  }
  mTables = mADCtoEMaker->getBemcData()->getTables();

  return kStOK;
}
//-----------------------------------------------------------------------------
StPicoTowerTest::~StPicoTowerTest()
{
  // mGRefMultCorrUtil is not owned by this maker (it is created outside and may be nullptr),
  // so it is NOT deleted here.
  delete mEmcPosition;
}

//-----------------------------------------------------------------------------
Int_t StPicoTowerTest::Finish()
{
  LOG_INFO << " StPicoTowerTest - writing data and closing output file " <<endm;

  if (!mOutputFile) return kStOK;
  mOutputFile->cd();

  vtxr->Write();
  MapTower_RunID_NHits_TowerID->Write();
  MapTower_RunID_Energy_TowerID->Write();
  MapTower_RunID_Energy_TowerID_energy->Write();
  MapTower_RunID_NEvents->Write();

  hRunPxl1->Write();
  hRunPxl2->Write();
  hRunIst ->Write();
  hRunSsd ->Write();

  mOutputFile->Close();
  delete mOutputFile;
  mOutputFile = nullptr;

  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoTowerTest::Make()
{
  ++mEventCounter;

  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoTowerTest - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  picoDst = mPicoDstMaker->picoDst();

  if (!picoDst)
  {
    LOG_WARN << "StPicoTowerTest - No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  if (!event)
  {
    LOG_WARN << "StPicoTowerTest - No StPicoEvent (bad read?)! Skip! " << endm;
    return kStWarn;
  }

  // Periodic save so a crash/stall does not lose everything (overwrites previous cycle)
  if (mOutputFile && mEventCounter % 50000 == 0) {
    mOutputFile->Write(0, TObject::kOverwrite);
  }

  if(!(isGoodEvent()))   // minBias trigger + vertex cuts
  {
    return kStOK;
  }

  StThreeVectorF pVtx(event->primaryVertex().x(), event->primaryVertex().y(), event->primaryVertex().z());
  vtxr->Fill(pVtx.x(),pVtx.y());

  // ---- per-run caches (function-local statics: one maker instance per process) ----
  // Run-index lookup and BEMC pedestals only change when the run changes, so they are
  // recomputed once per run instead of once per event (pedestals: 4800 DB lookups/event before).
  static unsigned int sRun    = 0;
  static int          sRunIdx = -1;
  static Float_t      sPed[4800];
  const unsigned int runId = event->runId();
  if (runId != sRun) {
    sRun    = runId;
    sRunIdx = this->runIndex(runId);
    if (sRunIdx >= 0) {
      Float_t rms;
      for (int i = 0; i < 4800; ++i) {
        sPed[i] = 0;
        mTables->getPedestal(BTOW, i+1, 0, sPed[i], rms);
      }
    }
  }
  const int runIdx = sRunIdx;
  if (runIdx < 0) return kStOK;
  MapTower_RunID_NEvents->Fill(runIdx);

  // Tower coefficients for the different runs, parameters are saved in BemcNewCalib.h
  const Double_t *TowerCoeff = (runId <= 15094020) ? CPre : CLowMidHigh;

  // loop over towers
  for (int iTow = 0; iTow < 4800; iTow++){
    StPicoBTowHit *towHit = picoDst->btowHit(iTow);
    if (!towHit) continue;

    const int realtowID = towHit->numericIndex2SoftId(iTow);
    if (realtowID < 1 || realtowID > 4800) continue;   // guard the TowerCoeff index

    // Calibrated energy E = C*(ADC-Pedestal)
    //const Double_t energy = TowerCoeff[realtowID-1]*(towHit->adc() - sPed[iTow]);
    const Double_t energy = towHit->energy();  // already calibrated 
    // E_T = E/cosh(eta) <= E, so E <= ETmin can never pass E_T > 0.2.
    // Skip before the (expensive) position calculation. Towers with bad status are
    // deliberately NOT skipped: this is a QA study that needs to see them.
    if (energy <= 0.2) continue;

    // Eta corrected for the vertex position (tabulated eta is w.r.t. TPC centre)
    StThreeVectorF towerPosition = mEmcPosition->getPosFromVertex(pVtx, realtowID);
    const double E_T = energy/cosh(towerPosition.pseudoRapidity());
    if (E_T > 0.2){
      MapTower_RunID_NHits_TowerID->Fill(runIdx, iTow);
      MapTower_RunID_Energy_TowerID_energy->Fill(runIdx, iTow, energy);
      MapTower_RunID_Energy_TowerID->Fill(runIdx, iTow, E_T);
    }
  }

  // HFT hit counts
  int pxl1 = 0, pxl2 = 0, ist = 0, ssd = 0;
  UInt_t nTracks = picoDst->numberOfTracks();
  for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack){
    StPicoTrack* trk = picoDst->track(iTrack);
    if(!trk) continue;

    UInt_t map = trk->hftHitsMap();
    pxl1 += (map >> 0) & 0x1;       // bit 0
    pxl2 += (map >> 1) & 0x3;       // bits 1-2
    ist  += (map >> 3) & 0x3;       // bits 3-4
    ssd  += (map >> 5) & 0x3;       // bits 5-6
  }
  hRunPxl1->Fill(runIdx, pxl1);
  hRunPxl2->Fill(runIdx, pxl2);
  hRunIst ->Fill(runIdx, ist);
  hRunSsd ->Fill(runIdx, ssd);

  return kStOK;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//Correct tower eta for Vz position
//-----------------------------------------------------------------------------
Double_t StPicoTowerTest::vertexCorrectedEta(double eta, double vz) {
    double tower_theta = 2.0 * atan(exp(-eta));
    double z = 0.0;
    if (eta != 0.0) z = mBarrelRadius / tan(tower_theta);
    double z_diff = z - vz;
    double theta_corr = atan2(mBarrelRadius, z_diff);
    double eta_corr = -log(tan(theta_corr / 2.0));
    return eta_corr;
}

// bool StPicoTowerTest::isGoodEvent()
// {
//   StPicoEvent *event = (StPicoEvent *)picoDst->event();
//   return (isMBTrigger(mYear) &&
//       sqrt(event->primaryVertex().x()*event->primaryVertex().x()+event->primaryVertex().y()*event->primaryVertex().y()) < mycuts::vz &&
//       fabs(event->primaryVertex().z()) < mycuts::vz &&
//       fabs(event->primaryVertex().z() - event->vzVpd()) < mycuts::vzVpdVz);
// }
bool StPicoTowerTest::isGoodEvent()
{
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  const float vz = event->primaryVertex().z();
  return (isMBTrigger(mYear) &&
      vz > mycuts::vzMin && vz < mycuts::vzMax &&
      fabs(vz - event->vzVpd()) < mycuts::vzVpdVz);
}

bool StPicoTowerTest::isMBTrigger(int mYear)
{
  const std::set<int>* mbTriggers = nullptr;

  if(mYear ==2016){
      static const std::set<int> mbTriggers2016 = {
        520001, 520011, 520021, 520031, 520041, 520051,           // VPDMB-5-p-sst
        570002, 570001,                                           // VPDMB-5-nosst / VPDMB-5-sst (production 2)
        520201, 520211, 520221, 520231, 520241, 520251, 520261,   // BHT1*VPDMB-10
        520203,                                                   // BHT3
        520101, 520111, 520121, 520131, 520141,                   // central-5
        520007, 520017, 520027, 520037,                           // vpdmb-10
        520003, 520013, 520023, 520033, 520043,                   // VPDMB-5
        520802, 520812, 520822, 520832, 520842,                   // VPDMB-5-p-hlt
        520002, 520012, 520022, 520032, 520042                    // VPDMB-5-p-nosst
      };
      mbTriggers = &mbTriggers2016;
  }

  if(mYear ==2014){
      static const std::set<int> mbTriggers2014 = {
      450050, 450060,   //VPDMB-5-p-nobsmd-hlt
      450005, 450015, 450025,
      450201, //BHT1*VPDMB-30
      450202, //BHT2*VPDMB-30
      450203, //BHT3
      450014, //VPDMB-5-nobsmd
      450009, //VPDMB-5-p-nobsmd-ssd-hlt
      };
      mbTriggers = &mbTriggers2014;
  }

  if(mYear ==2017){
      static const std::set<int> mbTriggers2017 = {
      570214, //BHT1*VPD30
      570204, //BHT1*VPD100
      570205, 570215, //BHT2*BBCMB
      570201, //BHT3
      570001, //VPDMB-30
      570404, //JP0*VPDMB30
      570403, //JP1*VPDMB30
      570401, //JP2
      570402, //JP2*L2JetHigh
      };
      mbTriggers = &mbTriggers2017;
  }
  if(mYear ==2018){
      static const std::set<int> mbTriggers2018 = {
      600001, // VPDMB-30 run 2018
      600011,
      600021,
      600031,
      600203, //BHT2 run 2018
      600213
      };
      mbTriggers = &mbTriggers2018;
  }

  if (!mbTriggers) return false;   // unknown year: previously a null-pointer dereference

  StPicoEvent* event = static_cast<StPicoEvent*>(picoDst->event());
  return std::any_of(mbTriggers->begin(), mbTriggers->end(), [&](int trigger) { return event->isTrigger(trigger); });
}
//-----------------------------------------------------------------------------
bool StPicoTowerTest::isGoodTrack(StPicoTrack const * const trk) const
{
  bool HFTCondition = (trk->hasPxl1Hit() && trk->hasPxl2Hit()) && (trk->hasSstHit() || trk->hasIstHit());
  return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && HFTCondition;
}

int StPicoTowerTest::runIndex(unsigned int run) const
{
  auto it = mRunIndexMap.find(run);
  return (it == mRunIndexMap.end()) ? -1 : it->second;
}