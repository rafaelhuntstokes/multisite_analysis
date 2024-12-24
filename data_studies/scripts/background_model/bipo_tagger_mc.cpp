#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <glob.h>
#include <iostream>
#include <fstream>
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/GeoUtils.hh>
#include <math.h>
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "THistPainter.h"
#include "TPaveStats.h"

void bipo_tagger_full(std::vector<std::string> filelist, double zOff, std::string outfname){
  
  double PROMPT_LOWER;
  double PROMPT_UPPER; 
  double DELAYED_LOWER;
  double DELAYED_UPPER;   
  double dR_cut;     
  double dT_low;    
  double dT_high;    
  double FV_CUT; 

  /* ### ANALYSIS CUTS ### */
  
  // implementing a GENERAL coincidence cut to get rid of (alpha, n), BiPo214 and BiPo212 together ...
  PROMPT_LOWER  = 350;    // FILL_ME_IN, MeV;
  PROMPT_UPPER  = 850;    // FILL_ME_IN, MeV;
  DELAYED_LOWER = 425;    // FILL_ME_IN, MeV;
  DELAYED_UPPER = 500;    // FILL_ME_IN, MeV;
  dR_cut        = 1000;   // FILL_ME_IN, mm;
  dT_low        = 400;   // FILL_ME_IN, ns;
  dT_high       = 800;// FILL_ME_IN, ns;
  FV_CUT        = 6000;   // FILL_ME_IN, mm;        // do this in some FV you want the efficiency for 

  /* ### ANALYSIS CUTS ### */

  // where to save the output ntuple
  std::string output_name = outfname;
  TFile *g = new TFile(output_name.c_str(), "RECREATE");
  std::string ntuple_name = "general_coincidence_cuts";
  
  // tree to save the coincidence pair information within the ntuple
  TTree *tag_info = new TTree("tag_info", "tag_info");

  // variables which handle the muon and high energy vetos
  ULong64_t cut_trig_time = 0;        // the event time according to 50 MHz clock that starts muon veto
  bool cut_status         = false;    // flag which says whether or not to skip event from consideration

  // variables used in the coincidence tagging + analysis
  double energy1, energy2, dR, dT, nhits_scaled1, nhits_scaled2, x1, x2, y1, y2, z1, z2;
  int nhits_cleaned1, nhits_cleaned2, nhits_raw1, nhits_raw2, gtid1, gtid2;
  
  /* 
  The variables to be saved in the coincidence ntuple output.
  */

  tag_info->Branch("energy_po", &energy1);
  tag_info->Branch("energy_bi", &energy2);
  tag_info->Branch("nhits_clean_po", &nhits_cleaned1);
  tag_info->Branch("nhits_clean_bi", &nhits_cleaned2);
  tag_info->Branch("nhits_raw_po", &nhits_raw1);
  tag_info->Branch("nhits_raw_bi", &nhits_raw2);
  tag_info->Branch("nhits_scaled_po", &nhits_scaled1);
  tag_info->Branch("nhits_scaled_bi", &nhits_scaled2);
  tag_info->Branch("gtid_po", &gtid1);
  tag_info->Branch("gtid_bi", &gtid2);
  tag_info->Branch("dR", &dR);
  tag_info->Branch("dT", &dT);
  tag_info->Branch("x_po", &x1);
  tag_info->Branch("y_po", &y1);
  tag_info->Branch("z_po", &z1);
  tag_info->Branch("x_bi", &x2);
  tag_info->Branch("y_bi", &y2);
  tag_info->Branch("z_bi", &z2);
  

  // begin coincidence tagging the MC files in the file list
  std::cout << "Number of files found: " << filelist.size() << std::endl;
  for (int i=0; i < filelist.size(); i++){
    std::string input_file  = filelist[i];
    TFile *f = TFile::Open(input_file.c_str());
    TNtuple *T = (TNtuple *) f->Get("output");
    
    double posX, posY, posZ, energy, nhitsScaled;
    int Nhits, NhitsRaw, mcIndex, ev, parent1, parent2, gtid;
    bool fit, fitValid;
    ULong64_t clock50, DCapplied, DCflagged;

    T->SetBranchAddress("posx", &posX);
    T->SetBranchAddress("posy", &posY);
    T->SetBranchAddress("posz", &posZ);
    T->SetBranchAddress("clockCount50", &clock50);
    T->SetBranchAddress("nhitsCleaned", &Nhits);
    T->SetBranchAddress("nhits", &NhitsRaw);
    T->SetBranchAddress("dcApplied", &DCapplied);
    T->SetBranchAddress("dcFlagged", &DCflagged);
    T->SetBranchAddress("scintFit", &fit);
    T->SetBranchAddress("fitValid", &fitValid);
    T->SetBranchAddress("mcIndex", &mcIndex);
    T->SetBranchAddress("parentpdg1", &parent1);
    T->SetBranchAddress("parentpdg2", &parent2);
    T->SetBranchAddress("energy", &energy);
    T->SetBranchAddress("evIndex", &ev);
    T->SetBranchAddress("eventID", &gtid);
    T->SetBranchAddress("correctedNhits", &nhitsScaled);

    // vector stores the i and j index of the events already paired up to remove duplicate pairings
    std::vector<int> paired; 

    std::cout << "Num. Entries: " << T->GetEntries() << std::endl;

    // deadtime is equal to (Number vetos x 20 s) + pileupTime + lone follower time
    bool veto_flag          = false;
    bool lone_follower_flag = false;
    double loneFollowerTime = 0;
    double pileupTime       = 0;
    int num_vetos           = 0;
    ULong64_t veto_start_time, lone_start_time;

    for (int iEntry = 0; iEntry < T->GetEntries(); iEntry++){
      T->GetEntry(iEntry);

      // check the reconstruction is valid (so we can rely on the energy reconstruction for highE veto)
      if (fit == false or fitValid == false){
        continue;
      }

      // looping backwards in time so skip first entry
      if (iEntry == 0){
        continue;
      }      

      nhits_cleaned1 = Nhits;
      nhits_scaled1  = nhitsScaled;
      nhits_raw1     = NhitsRaw;
      gtid1          = gtid;

      double rho1, r1;
      ULong64_t flagged1; 
      flagged1 = DCflagged;
      energy1 = energy;

      x1   = posX;
      y1   = posY;
      z1   = posZ - zOff;  // account for AV offset
      rho1 = pow(x1*x1 + y1*y1, 0.5);
      r1   = pow(x1*x1 + y1*y1 + z1*z1, .5);

      // check event is inside the FV
      if (r1 > FV_CUT){
        continue;
      }

      // energy cuts
      if (nhits_scaled1 < DELAYED_LOWER || nhits_cleaned1 > DELAYED_UPPER){
        continue;
      }
      
      ULong64_t clock1 = clock50;
      
      // all cuts passed and we've got a Po candidate; loop backwards to find Bi
      for (int jEntry = iEntry - 1; jEntry >= 0; jEntry--){
          T->GetEntry(jEntry);
          ULong64_t clock2;
          clock2 = clock50;
          
          // apply dT cut
          dT = (clock1 - clock2) * 20;
          if (dT > dT_high){
              break;
          }
          if (dT < dT_low){
              continue;
          }

          // check if event is already paired up
          bool pair = false;
          for (int iPair = 0; iPair < paired.size(); iPair++){
            // loop over every entry in iPair and see if already paired
            
            if (jEntry == paired[iPair]){
              pair = true;
              break;
            }
          }
          if (pair == true){
            // this bismuth has already been paired up with something! skip.
            continue;
          }
          // else it hasn't been paired already and we're free to continue our checks

          // check reconstruction valid
          if (fit == false){
            continue;
          }
          if (fitValid == false){
            continue;
          }
          ULong64_t flagged2; 
          flagged2 = DCflagged;
          
          double rho2, r2;
          nhits_cleaned2 = Nhits;
          nhits_raw2     = NhitsRaw;
          nhits_scaled2  = nhitsScaled;
          gtid2          = gtid;
          
          
          x2   = posX;
          y2   = posY;
          z2   = posZ - zOff;
          rho2 = pow(x2*x2 + y2*y2, 0.5);
          r2   = pow(x2*x2 + y2*y2 + z2*z2, .5);

          // check event is inside the FV
          if (r2 > FV_CUT){
            continue;
          }

          energy2 = energy;

          // check if nhit cuts satisfied for Bi214 candidate
          if (nhits_scaled2 < PROMPT_LOWER || nhits_scaled2 > PROMPT_UPPER){
            continue;
          }


          // apply dR cut
          dR = pow(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2), 0.5);
          if (dR > dR_cut){
            continue;
          }

          // all cuts passed and should be a unique BiPo214 pair
          paired.push_back(iEntry);
          paired.push_back(jEntry);

          tag_info->Fill();
          // exit Bi214 candidate search as we have already paired it up
          break; 
          }
      }
  }

  // save the output NTUPLE for this run
  g->cd();
  tag_info->Write();
  g->Close();
}


std::vector<std::string> globVector(const std::string& pattern){
  glob_t glob_result;
  glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
  std::vector<std::string> files;
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
    files.push_back(std::string(glob_result.gl_pathv[i]));
    std::cout<<std::string(glob_result.gl_pathv[i])<<std::endl;
  }
  globfree(&glob_result);
  return files;
}

int main(int argc, char* argv[]){

  // inputs are a run number and the name of the isotope to tag
  std::string runNum      = argv[1];
  std::string isotope     = argv[2];
  
  std::string inFileName = "/data/snoplus2/miniPROD_bismsb_newOptics_oldFitter7015/ntuples/" + isotope + "_" + runNum + "_" + "bismsb.ntuple.root"; 
  
  // where to save the output ntuple
  std::string outfname = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/powei_212/" + runNum + ".root";
  
  // load all the ntuple files corresponding to this isotope into a vector
  std::vector<std::string> filelist = globVector(inFileName);

  // if we find no files, skip the run!
  if (filelist.size() == 0){
    std::cout << "Found no valid files for: " << inFileName << std::endl;
    return 1;
  }

  // Begin run and load DB
  int runID = std::stoi(runNum);                                                                                                                                                                 
  RAT::DB *db = RAT::DB::Get();
  db->LoadDefaults();
  RAT::DS::Run run;
  run.SetRunID(runID);
  db->BeginOfRun(run);
  std::cout << "Run ID = " << runID << std::endl;

  // Get AV Shift for position reconstruction FV cuts to be accurate                                                                                                                                                       
  std::vector<double> AVPos = RAT::GeoUtil::UpdateAVOffsetVectorFromDB();
  double zOff = AVPos[2];
  std::cout << "AV Shift is: " << zOff << std::endl;

  // get the run time ... for no reason
  double start_day, start_sec, start_nsec, end_day, end_sec, end_nsec;
  long long start_time, end_time, duration;
  RAT::DBLinkPtr dblink = db->GetLink("RUN");
  start_day = dblink->GetD("start_day");
  start_sec = dblink->GetD("start_sec");
  start_nsec = dblink->GetD("start_nsc");
  end_day = dblink->GetD("stop_day");
  end_sec = dblink->GetD("stop_sec");
  end_nsec = dblink->GetD("stop_nsc");
  start_time = start_day * 24 * 3600 * pow(10, 9) + start_sec * pow(10, 9) + start_nsec;
  end_time = end_day * 24 * 3600 * pow(10, 9) + end_sec * pow(10, 9) + end_nsec;
  duration = end_time - start_time;
  std::cout << "Duration of run: " << duration * pow(10, -9) << " s" << std::endl;

  // run the general coincidence tagging!
  bipo_tagger_full(filelist, zOff, outfname);
  return 0;
}