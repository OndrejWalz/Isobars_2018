#ifndef __CINT__
#include "TString.h"
#include <iostream>
#include <fstream>
#include <memory>
#include <ctime>
#include <cstdio>
#include <string>
#include <vector>
#include <TROOT.h>
#include <TSystem.h>
#include <TEnv.h>
#include "TChain.h"
#include "StChain/StMaker.h"
#include "StChain/StChain.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StPicoTowerTest/StPicoTowerTest.h"
#include "St_db_Maker/St_db_Maker.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
using namespace std;

#else
class StChain;
#endif

StChain *chain;

void loadSharedLibraries();
void loadSharedAnalysisLibraries2();
void progres(double N, double D);

// Take files [chunk*perChunk, (chunk+1)*perChunk) from a file list and write them
// to a small list file. Returns the new list name ("" if the chunk is empty/unreadable).
string makeChunkList(const string& inList, int chunk, int perChunk)
{
   ifstream in(inList.c_str());
   if (!in.is_open()) {
      cout << "\033[0;31mCannot open file list: " << inList << "\033[0m" << endl;
      return "";
   }
   vector<string> files;
   string line;
   while (getline(in, line)) {
      if (!line.empty()) files.push_back(line);
   }
   size_t first = (size_t)chunk * perChunk;
   size_t last  = first + perChunk;
   if (first >= files.size()) {
      cout << "Chunk " << chunk << " is beyond the end of the list (" << files.size() << " files)." << endl;
      return "";
   }
   if (last > files.size()) last = files.size();

   string outList = Form("chunk_%d.list", chunk);
   ofstream out(outList.c_str());
   for (size_t i = first; i < last; ++i) out << files[i] << "\n";
   out.close();
   cout << "Chunk " << chunk << ": files " << first << "-" << last - 1
        << " of " << files.size() << " -> " << outList << endl;
   return outList;
}

// Copy every remote file of a list to node-local scratch (once, up front) and return a
// list pointing at the local copies. Local reads are far faster than streaming over
// XRootD and avoid the dropped-connection stalls. Falls back to the remote path on failure.
string gStageDir = "";
string stageFilesLocally(const string& inList, const string& tag)
{
   ifstream in(inList.c_str());
   if (!in.is_open()) return inList;
   const char* tmp = gSystem->Getenv("TMPDIR");
   gStageDir = Form("%s/pico_%s_%d", tmp ? tmp : "/tmp", tag.c_str(), (int)gSystem->GetPid());
   gSystem->mkdir(gStageDir.c_str(), kTRUE);

   string outList = Form("%s/staged.list", gStageDir.c_str());
   ofstream out(outList.c_str());
   string line;
   int n = 0, ok = 0;
   while (getline(in, line)) {
      if (line.empty()) continue;
      ++n;
      if (line.find("root://") == string::npos) { out << line << "\n"; ++ok; continue; } // already local
      string local = Form("%s/%s", gStageDir.c_str(), gSystem->BaseName(line.c_str()));
      bool copied = false;
      for (int attempt = 0; attempt < 3 && !copied; ++attempt)
         copied = (gSystem->Exec(Form("xrdcp -f -s %s %s", line.c_str(), local.c_str())) == 0);
      if (copied) { out << local << "\n"; ++ok; }
      else { cout << "Staging failed, reading remotely: " << line << endl; out << line << "\n"; }
   }
   out.close();
   cout << "Staged " << ok << "/" << n << " files to " << gStageDir << endl;
   return outList;
}

// chunk < 0  : process the whole list (old behaviour)
// stage = true : copy the input files to local scratch first ($TMPDIR or /tmp), removed at the end
// chunk >= 0 : process only that chunk of 'filesPerChunk' files, output name gets _chunkN suffix
void RunPicoTowerTest(string pico="filelists/test_list.list",
 string outFileName="Test_Tower_test.root", int pYear = 2018, bool Testing = false,
 int chunk = -1, int filesPerChunk = 20, bool stage = false)
{
   // Fail fast on a bad XRootD server instead of stalling for ~30 min.
   // (Names are for the ROOT legacy XNet client; check them against your ROOT version.)
   gEnv->SetValue("XNet.RequestTimeout",  60);
   gEnv->SetValue("XNet.ConnectTimeout",  20);
   gEnv->SetValue("XNet.MaxRedirectCount", 4);
   gEnv->SetValue("XNet.ReconnectWait",    5);

   string SL_version = "Unknown";
   string RefMult = "grefmult";
   string Runcode;
   string prodID;
   if (pYear==2016){
      SL_version = "SL20c";
      Runcode = "Run16_AuAu200_VpdMB5";
      prodID = "P16ij";
   } else if (pYear==2014){
      SL_version = "SL22c";
      Runcode = "Run14_AuAu200_VpdMB5";
      prodID = "P16id";
   } else if (pYear==2017){
      SL_version = "SL22b";
   } else if (pYear==2018){
      SL_version = "SL20c";
      prodID = "P20ic";
   } else {
      cout << "\033[0;31m Not valid year.\033[0m"<<endl;
      exit(0);
   }

   Int_t nEntries = 100;

   if(Testing){
      pico=Form("../../TestLists/testPico_%.d.list",pYear);
      outFileName=Form("../../Test_Tower_%.d.root",pYear);
   }

   // Optional chunking: one small job -> one output file, merge later with hadd
   if (!Testing && chunk >= 0) {
      pico = makeChunkList(pico, chunk, filesPerChunk);
      if (pico.empty()) exit(0);
      size_t pos = outFileName.rfind(".root");
      string base = (pos == string::npos) ? outFileName : outFileName.substr(0, pos);
      outFileName = Form("%s_chunk%d.root", base.c_str(), chunk);
   }

   if (stage && !Testing) pico = stageFilesLocally(pico, (chunk >= 0) ? Form("c%d", chunk) : "all");

   string env_SL = getenv("STAR");
   if (env_SL.find(SL_version) == string::npos)
   {
      cout << "\033[0;31mEnvironment Star Library does not match the requested library: \033[0m" << "\033[0;32m" << SL_version << "\033[0m";
      cout << "\033[0;31m for run: \033[0m" << "\033[0;32m" << pYear << "\033[0m";
      cout << "\033[0;31m in RunPicoTowerTest.C. Exiting...\033[0m" << endl;
      exit(1);
   }

   // 1. Load core STAR framework dependencies FIRST
   gSystem->Load("libTable");
   gSystem->Load("St_base");
   gSystem->Load("StChain");
   gSystem->Load("StUtilities");

   // 2. Load standard HF / analysis helper libraries
   gROOT->LoadMacro("loadSharedHFLibraries.C");
   loadSharedHFLibraries();

   // 3. Load STAR system utilities
   gSystem->Load("StBTofUtil");

   // 4. Load your local compiled libraries
   gSystem->Load("StRefMultCorr");
   gSystem->Load("StPicoCuts");
   gSystem->Load("StPicoPrescales");
   gSystem->Load("StPicoTowerTest");

   chain = new StChain();

   StRefMultCorr *mRefmultCorrUtil = 0;   // was uninitialized (undefined behaviour in destructor); 0 not nullptr: CINT/ROOT5
   StPicoDstMaker* picoDstMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, pico.c_str(), "picoDstMaker");
   // Read only the branches this analysis uses (default reads all of them: TOF, MTD, BBC, EPD, ...)
   picoDstMaker->SetStatus("*", 0);
   picoDstMaker->SetStatus("Event",   1);
   picoDstMaker->SetStatus("Track",   1);   // hftHitsMap
   picoDstMaker->SetStatus("BTowHit", 1);
   StMessMgr *msg = StMessMgr::Instance();
   msg->SwitchOff("Could not make BEMC detector");
   St_db_Maker *dbMaker = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");
   StEmcADCtoEMaker *adc = new StEmcADCtoEMaker();
   StPicoTowerTest*  picoTowerTest = new StPicoTowerTest("picoTowerTest", outFileName.c_str(), picoDstMaker,mRefmultCorrUtil,pYear);

   // -------------- USER variables -------------------------
   picoTowerTest->setCutETmin(0.2);
   picoTowerTest->setHadronCorr(1.);
   picoTowerTest->setMaxDcaZHadronCorr(3.0); //cm, max DCA_z for global tracks used for hadronic correction

   if (chain->Init() != kStOK) {
      cout << "\033[0;31mchain->Init() failed. Exiting.\033[0m" << endl;
      exit(1);
   }

   // TChain::GetEntries() opens EVERY file in the list. Over XRootD with thousands of files
   // that takes tens of minutes, so only do it for small (chunked) or local (staged) inputs.
   bool knownTotal = true;
   if(!Testing){
      if (chunk >= 0 || stage) {
         nEntries = picoDstMaker->chain()->GetEntries();
      } else {
         nEntries = 2000000000;   // just loop until end of input (kStEOF)
         knownTotal = false;
      }
   }

   if (knownTotal) cout << "nEntries: " << nEntries << endl;
   else            cout << "nEntries: not counted (running until end of input)" << endl;

   time_t tStart = time(0);
   //--------------------------------------------------------
   for (int iEvent = 0; iEvent < nEntries; ++iEvent)
   {
      chain->Clear();
      int iret = chain->Make();
      if (Testing && iEvent%200==0) progres(iEvent,nEntries);
      if (iEvent % 10000 == 0) {                       // heartbeat so stalls are visible
         time_t now = time(0);
         double rate = (now > tStart) ? iEvent / double(now - tStart) : 0.;
         cout << "Event " << iEvent;
         if (knownTotal) cout << " / " << nEntries;
         cout << "  (" << Form("%.0f", rate) << " ev/s)  " << ctime(&now) << flush;
      }
      if (iret)
      {
         // kStEOF (end of input) and real errors both end here; Finish() below still writes output
         cout << "Stopping at event " << iEvent << ", return code " << iret << endl;
         break;
      }
   }
   if(Testing && nEntries%200!=0) progres(nEntries,nEntries);
   //--------------------------------------------------------

   chain->Finish();
   delete chain;
   if (!gStageDir.empty()) gSystem->Exec(Form("rm -rf %s", gStageDir.c_str()));
}

void progres(double citatel, double jmenovatel){
      int Ndilky=50;
      int proc=floor(citatel/jmenovatel*Ndilky);

      cout << "\r"<< flush;
      cout << "  │" << flush;
      for (int i=1; i<=proc; i++){
      cout <<"█"<<flush;
      }
      for (int j=proc+1; j<=Ndilky;j++){
      cout <<  "░"<<flush;
      }
      cout << "│ " << flush;
      if(citatel!=jmenovatel){
         cout << Form("Completed: %.2f ",citatel/jmenovatel*100.)<<"%"<<flush;
      }
      else{
         cout << Form("\033[1;32m Completed: %.2f \033[0m",citatel/jmenovatel*100.)<<"\033[1;32m%\033[0m"<<endl;
      }
return;
}