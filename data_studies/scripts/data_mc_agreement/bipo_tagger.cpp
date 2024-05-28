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

double bipo_tagger_full(std::vector<std::string> filelist, std::string outdir, double zOff, int run_num) {

  int num_vetos      = 0;        // number of times the high E / muon veto is triggered        
  int failedFit      = 0;
  int failedEnergy   = 0;
  int failedBitmask  = 0;
  int failedCleaning = 0;  
  
  int num_tagged     = 0;
  
  double BI_LOWER;
  double BI_UPPER; 
  double PO_LOWER;
  double PO_UPPER;   
  double dR_cut;     
  double dT_low;    
  double dT_high;    
  double FV_CUT; 

  // ANALYSIS CUTS 
  BI_LOWER   = 1.25; // energy cuts MeV
  BI_UPPER   = 3.0; 
  PO_LOWER   = 0.7;
  PO_UPPER   = 1.1; 
  dR_cut     = 1000; // mm  
  dT_low     = 3690; // ns
  dT_high    = 1000000; // ns
  FV_CUT     = 6000; // mm 
  
  // variables which handle the muon and high energy vetos
  ULong64_t cut_trig_time = 0;        // the event time according to 50 MHz clock that starts muon veto
  bool cut_status         = false;    // flag which says whether or not to skip event from consideration
  
  // set up the output file location and name
  std::string output_name = outdir + "/" + std::to_string(run_num) + ".root";
  TFile *g = new TFile(output_name.c_str(), "RECREATE");
  
  // variables of interest to save
  double energy1, energy2, dR, dT, nhits_scaled1, nhits_scaled2, x1, x2, y1, y2, z1, z2, itr1, itr2, posFOM1, posFOM2;
  float HC_ratio1, HC_ratio2;
  int nhits_clean1, nhits_clean2, nhits_raw1, nhits_raw2, gtid1, gtid2;
  UInt_t posFOM1_hits, posFOM2_hits;
  // define the output TTree variables
  TTree *tag_info = new TTree("tag_info", "tag_info");
  
  tag_info->Branch("energy_po", &energy1);
  tag_info->Branch("energy_bi", &energy2);
  tag_info->Branch("nhits_clean_po", &nhits_clean1);
  tag_info->Branch("nhits_clean_bi", &nhits_clean2);
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
  tag_info->Branch("hc_po", &HC_ratio1);
  tag_info->Branch("hc_bi", &HC_ratio2);
  tag_info->Branch("itr_po", &itr1);
  tag_info->Branch("itr_bi", &itr2);
  tag_info->Branch("posFOM_po", &posFOM1);
  tag_info->Branch("posFOM_po_hits", &posFOM1_hits); // nhits to scale posFOM with
  tag_info->Branch("posFOM_bi", &posFOM2);
  tag_info->Branch("posFOM_po_hits", &posFOM2_hits);

  std::cout << "Number of files found: " << filelist.size() << std::endl;

  // create a TChain that contains the entries of all the events in this run
  TChain *chain = new TChain("output");
  for (int i=0; i < filelist.size(); i++){
    std::string input_file  = filelist[i];
    chain->AddFile(input_file.c_str());
  }
  std::cout << "TChain contains " << chain->GetEntries() << " entries." << std::endl;
  
  // define the variables extracted from each event in the TChain
  double posX, posY, posZ, energy, nhits_scaled, posFOM, itr;
  int nhits_clean, nhits_raw, ev, gtid;
  UInt_t posFOM_hits;
  Int_t triggerWord;
  bool fit, fitValid;
  ULong64_t clock50, DCapplied, DCflagged;
  chain->SetBranchAddress("posx", &posX);
  chain->SetBranchAddress("posy", &posY);
  chain->SetBranchAddress("posz", &posZ);
  chain->SetBranchAddress("energy", &energy);
  chain->SetBranchAddress("nhitsCleaned", &nhits_clean);
  chain->SetBranchAddress("nhits", &nhits_raw);
  chain->SetBranchAddress("correctedNhits", &nhits_scaled);
  chain->SetBranchAddress("dcApplied", &DCapplied);
  chain->SetBranchAddress("dcFlagged", &DCflagged);
  chain->SetBranchAddress("scintFit", &fit);
  chain->SetBranchAddress("fitValid", &fitValid);
  chain->SetBranchAddress("clockCount50", &clock50);
  chain->SetBranchAddress("evIndex", &ev);
  chain->SetBranchAddress("eventID", &gtid);
  chain->SetBranchAddress("itr", &itr);
  chain->SetBranchAddress("posFOM", &posFOM);
  chain->SetBranchAddress("posFOM2", &posFOM_hits);
  chain->SetBranchAddress("triggerWord", &triggerWord);
  
  // loop over the entries for the BiPo214 analysis
  int numEntries = chain->GetEntries();
  int numMuons   = chain->GetEntries("(dcFlagged&0x80)!=0x80");
  int numHighE   = chain->GetEntries("(nhitsCleaned>5000)");
  std::cout << "Total muons + highE events in file: " << numMuons << " + " << numHighE << " = " << numMuons+numHighE << std::endl;

  // vector stores the i and j index of the events already paired up to remove duplicate pairings
  std::vector<int> paired; 

  std::cout << "Num. Entries: " << numEntries << std::endl;

  // deadtime is equal to (Number vetos x 20 s) + pileupTime + lone follower time
  bool veto_flag          = false;
  bool lone_follower_flag = false;
  double loneFollowerTime = 0;
  double pileupTime       = 0;
  bool track_flag = false;
  ULong64_t veto_start_time, lone_start_time, previous_entry_time;

  for (int iEntry = 0; iEntry < chain->GetEntries(); iEntry++){
    chain->GetEntry(iEntry);
    
    // check if the event is flagged as an ORPHAN!
    if (triggerWord == 0){
        std::cout << "Orphan! Skipping GTID: " << gtid << std::endl;
        continue;
    }
    // check if event is high energy or tagged as a muon by DC
    if (nhits_clean > 5000 or ((DCflagged&0x80) != 0x80)){
        std::cout << "\nDeadtime Triggered - GTID: " << gtid << std::endl;
        std::cout << "NhitsCleaned: " << nhits_clean << std::endl;
        std::cout << "Muon Flag: " << ((DCflagged&0x80)!=0x80) << std::endl;
        std::cout << "Veto Flag Status: " << veto_flag << std::endl;
        // check if we are inside a lone muon follower veto window
        if (lone_follower_flag == true){
            // add the dT and switch off the muon follower veto
            ULong64_t deltaT = ((clock50-lone_start_time) & 0x7FFFFFFFFFF)*20.0;
            loneFollowerTime += deltaT;
            lone_follower_flag = false;
        }

        if (veto_flag == false){
            num_vetos++;
            veto_flag = true;
            std::cout << "Setting veto flag to " << veto_flag << std::endl;
            std::cout << "Number of Vetos: " << num_vetos << std::endl;
            veto_start_time = clock50;
            continue;
        }
        else{
            // we have pileup! Need to calculate the additional pileup time.
            ULong64_t deltaT = ((clock50-veto_start_time) & 0x7FFFFFFFFFF)*20.0;
            pileupTime += deltaT;
            std::cout << "Pileup: Added: " << deltaT * pow(10, -9) << " s." << std::endl;
            // reset the veto window
            veto_start_time = clock50;
            continue;
        }
    }

    // now handle the veto window for follower events from highE / muon veto
    if (veto_flag == true){
        ULong64_t deltaT = ((clock50-veto_start_time) & 0x7FFFFFFFFFF)*20.0;

        // check if event falls inside veto window
        if (deltaT < 2.0e10){
          continue;
        }
        else{
          // we are no longer inside the veto time window --> switch it off!
          veto_flag = false;
          std::cout << "Time since veto window activated: " << deltaT * pow(10, -9) << " s." << std::endl;
          std::cout << "Setting veto_flag to " << veto_flag << std::endl;
        }
    }

    // check if we have a lone muon follower (i.e. muon at end of previous run)
    // this can only trigger if there was a muon at the end of the previous run
    if ((DCflagged&0x4000)!=0x4000){
        std::cout << "Event tagged as lone follower!" << std::endl;
        if (lone_follower_flag == false){
          lone_follower_flag = true;
          lone_start_time = clock50;
          continue;
        } else{
          // we are within a lone follower veto period!
          ULong64_t deltaT = ((clock50 - lone_start_time) & 0x7FFFFFFFFFF)*20.0;
          loneFollowerTime += deltaT;
          lone_start_time = clock50;
          continue;
        }
    }
    else {
        // check if event is NOT flagged as muon follower, but muon follower flag is on
        if (lone_follower_flag == true){
            // but we haven't triggered the other blocks! --> we've left the muon follower deadtime window
            lone_follower_flag = false;
        }
    }

    // looping backwards in time so skip first entry
    if (iEntry == 0){
      // record the first event time to begin the retrigger veto
      previous_entry_time = clock50;
      std::cout << "iEntry == 0 so skipping." << std::endl;
      continue;
    }

    // find the time of this event to the previous one
    ULong64_t retrigger_cut_time = ((clock50 - previous_entry_time) & 0x7FFFFFFFFFF)*20.0; // in ns

    // reset the previous event time regardless of the outcome of retrigger cut
    previous_entry_time = clock50;

    // apply retrigger cut
    if (retrigger_cut_time < 500){
        // event is a retrigger!
        continue;
    }

    // check the reconstruction is valid (so we can rely on the energy reconstruction for highE veto)
    if (fit == false or fitValid == false){
      failedFit += 1;
      continue;
    }

    // int nhits1, nhitsRaw1, gtid1;
    nhits_clean1  = nhits_clean;
    nhits_scaled1 = nhits_scaled;
    nhits_raw1    = nhits_raw;
    gtid1         = gtid;

    double rho1, r1;
    ULong64_t flagged1; 
    flagged1 = DCflagged;
    energy1 = energy;

    if (((DCapplied & 0x2100000042C2) & DCflagged ) != (DCapplied & 0x2100000042C2)){
        continue;
      }

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
    if (energy1 < PO_LOWER || energy1 > PO_UPPER){
        continue;
    }
    
    // check hit cleaning cleanliness cut
    HC_ratio1 = float(nhits_clean1) / float(nhits_raw1);

    // save the position of merit details from scintFitter
    posFOM1 = posFOM;
    posFOM1_hits = posFOM_hits;

    // save the clock 50 counts to variable to work out dT to bi candidates
    ULong64_t clock1 = clock50;

    // all cuts passed and we've got a Po candidate; loop backwards to find Bi
    for (int jEntry = iEntry - 1; jEntry >= 0; jEntry--){
        chain->GetEntry(jEntry);
        ULong64_t clock2;
        clock2 = clock50;
        
        // apply dT cut
        dT = ((clock1 - clock2) & 0x7FFFFFFFFFF) * 20; // in ns
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
        failedFit += 1;
        continue;
      }
      if (fitValid == false){
        failedFit += 1;
        continue;
      }

      ULong64_t flagged2; 
      flagged2 = DCflagged;
      if (((DCapplied & 0x2100000042C2) & DCflagged ) != (DCapplied & 0x2100000042C2)) continue;
      
      double rho2, r2;
      nhits_clean2  = nhits_clean;
      nhits_raw2    = nhits_raw;
      nhits_scaled2 = nhits_scaled;
      gtid2         = gtid;
      
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
      // check if energy cuts satisfied for Bi214 candidate
      if (energy2 < BI_LOWER || energy2 > BI_UPPER){
        continue;
      }

      // apply dR cut
      dR = pow(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2), 0.5);
      if (dR > dR_cut){
        continue;
      }

      // check hit cleaning cleanliness cut
      HC_ratio2 = float(nhits_clean2) / float(nhits_raw2);
      
      // save the position of merit details from scintFitter
      posFOM2 = posFOM;
      posFOM2_hits = posFOM_hits;
      
      // all cuts passed and should be a unique BiPo214 pair
      paired.push_back(iEntry);
      paired.push_back(jEntry);

      tag_info->Fill();

      num_tagged += 1;
      break; 
    }
  }

  // save the output NTUPLE for this run
  g->cd();
  tag_info->Write();
  g->Close();

  double deadtime =  (20 * num_vetos) + (pileupTime + loneFollowerTime) * pow(10, -9); // in s
  std::cout << "Num Tagged: " << num_tagged << std::endl;
  std::cout << "PileupTime: " << pileupTime * pow(10, -9) << " s" << std::endl;
  std::cout << "Num. highE vetos: " << num_vetos << std::endl;
  std::cout << "highE deadtime: " << deadtime << " s" << std::endl;

  std::cout << "% failed fit: " << float(failedFit) / float(numEntries) << std::endl;
  std::cout << "% failed energy: " << float(failedEnergy) / float(numEntries) << std::endl;
  std::cout << "% failed cleanliness: " << float(failedCleaning)  / float(numEntries) << std::endl;

  return deadtime;
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

  std::string runNum      = argv[1];
  std::string out_dir     = argv[2];
  
  std::vector<int> tagList;
  int iRun;

  std::string inFile1, ratds1;
  std::string data_path  = "/data/snoplus3/SNOplusData/processing/fullFill/rat-7.0.8/ntuples";
  std::string ratds_path = "/data/snoplus3/SNOplusData/processing/fullFill/rat-7.0.8/ratds";
  if (std::stoi(runNum) < 307613){
    inFile1 = data_path  + "/" + "Analysis20R_r0000";
    ratds1  = ratds_path + "/" + "Analysis20R_r0000";
  }
  if (std::stoi(runNum) >= 307613){
      inFile1 = data_path + "/"  + "Analysis20_r0000";
      ratds1  = ratds_path + "/" + "Analysis20_r0000";
  }
  
  std::string inFile2 = "*";
  const std::string inFileName   = inFile1+runNum+inFile2;
  const std::string ratdsFileName = ratds1 + runNum + inFile2; 
  std::cout << inFileName << std::endl;
  std::cout << "Ratds Fname: " << ratdsFileName << std::endl;
  std::vector<std::string> filelist = globVector(inFileName);

  // check if the corresponding RATDS files exist
  std::vector<std::string> ratds_filelist = globVector(ratdsFileName);

  // check if the number of files in the ratds vector matches the ntuple
  if (filelist.size() != ratds_filelist.size()){
    std::cout << "Num. ntuple subruns: " << filelist.size() << std::endl;
    std::cout << "Num. RATDS subruns: " << ratds_filelist.size() << std::endl;
    std::cout << "Not tagging run due to discrepency!" << std::endl;
    return 0;
  }
  // if we find no subrun files, skip the run!
  if (filelist.empty()){
    std::cout << "No valid ntuple files!" << std::endl;
    return 0;
  } else {
    // Begin run and load DB
    int runID = std::stoi(runNum);                                                                                                                                                                 
    RAT::DB *db = RAT::DB::Get();
    db->LoadDefaults();
    RAT::DS::Run run;
    run.SetRunID(runID);
    db->BeginOfRun(run);
    std::cout << "Run ID = " << runID << std::endl;
    
    // Get AV Shift                                                                                                                                                                     
    std::vector<double> AVPos = RAT::GeoUtil::UpdateAVOffsetVectorFromDB();
    double zOff = AVPos[2];
    std::cout << "AV Shift is: " << zOff << std::endl;

    // get the run time
    double start_day, start_sec, start_nsec, end_day, end_sec, end_nsec;
    long long start_time, end_time, duration;
    RAT::DBLinkPtr dblink = db->GetLink("RUN");
    start_day = dblink->GetD("start_day");
    start_sec = dblink->GetD("start_sec");
    start_nsec = dblink->GetD("start_nsc");
    end_day = dblink->GetD("stop_day");
    end_sec = dblink->GetD("stop_sec");
    end_nsec = dblink->GetD("stop_nsc");
    start_time = start_day * 24 * 3600 + start_sec;
    end_time = end_day * 24 * 3600 + end_sec;
    duration = end_time - start_time;
    std::cout << "Duration of run: " << duration << " s" << std::endl;

    double deadtime = bipo_tagger_full(filelist, out_dir, zOff, runID);
    double adjusted_livetime = duration - deadtime; 

    std::cout << "Adjusted Livetime: " << adjusted_livetime << " s" << std::endl;
    
    // create output txt file containing the livetime
    std::ofstream livetime;
    livetime.open(out_dir + "/livetime/" + runNum + ".txt");
    livetime << adjusted_livetime;
    livetime.close();
    
    return 0;
  }
}
