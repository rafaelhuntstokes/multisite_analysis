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

double bipo_tagger_full(std::vector<std::string> filelist, double zOff, int run_num, int isotope) {

  std::cout << "ISOTOPE: " << isotope << std::endl;
  int num_vetos      = 0;        // number of times the high E / muon veto is triggered        
  int failedFit      = 0;
  int failedEnergy   = 0;
  int failedBitmask  = 0;
  int failedCleaning = 0;  
  
  int numTagged_6m    = 0;
  int numTagged_4m    = 0;
  int numTagged_donut = 0;
  int totalEntries    = 0;
  
  double BI_LOWER;
  double BI_UPPER; 
  double PO_LOWER;
  double PO_UPPER;   
  double dR_cut;     
  double dT_low;    
  double dT_high;    
  double FV_CUT; 

  bool bisMSB = false;
  // ANALYSIS CUTS 
  if (isotope == 212){
    BI_LOWER   = 0.6; // energy cuts MeV
    BI_UPPER   = 3.0;
    PO_LOWER   = 0.85;
    PO_UPPER   = 1.3;
    dR_cut     = 1000; // mm  
    dT_low     = 400;//3690; // ns
    dT_high    = 800;//1000000; // ns
    FV_CUT     = 6000; // mm 
    if (bisMSB == true){
      // bisMSB 
      BI_LOWER   = 1.4; // energy cuts MeV
      BI_UPPER   = 4.0; 
      PO_LOWER   = 1.36; 
      PO_UPPER   = 2.1;
    }

  }
  if (isotope == 214){
    std::cout << "ISOTOPE 214" << std::endl;
    BI_LOWER   = 1.25; // energy cuts MeV
    BI_UPPER   = 3.0; 
    PO_LOWER   = 0.7;
    PO_UPPER   = 1.1; 
    dR_cut     = 1000; // mm  
    dT_low     = 3690; // ns
    dT_high    = 1000000; // ns
    FV_CUT     = 6000; // mm 
    //if (run_num >= 311234){
    if (bisMSB == true){ 
      // bisMSB 
      BI_LOWER   = 2.0; // energy cuts MeV
      BI_UPPER   = 4.8;
      PO_LOWER   = 1.1; 
      PO_UPPER   = 1.8;
    }
    std::cout << "Cut values are: " << BI_LOWER << " --->" << BI_UPPER << std::endl;
  }
  // FV of the torus
  double MAJOR_R = 3500; // mm
  double TUBE_R  = 1500; // mm
  std::cout << "Cut values are: " << BI_LOWER << " --->" << BI_UPPER << std::endl;
  // variables which handle the muon and high energy vetos
  ULong64_t cut_trig_time = 0;        // the event time according to 50 MHz clock that starts muon veto
  bool cut_status         = false;    // flag which says whether or not to skip event from consideration

  // std::string output_name = "/data/snoplus3/hunt-stokes/clean_multisite/true_rate_efficiency_tag_output/output_" + std::to_string(isotope) + "_" + std::to_string(run_num) + ".root";
  std::string output_name = "/data/snoplus3/hunt-stokes/clean_multisite/bipo_" + std::to_string(isotope) + "_cleanSelec3/output_" + std::to_string(isotope) + "_" + std::to_string(run_num) + ".root";
  // std::string output_name = "/data/snoplus3/hunt-stokes/clean_multisite/bipo212_mc/output_" + std::to_string(run_num) + ".root";
  // std::string output_name = "/data/snoplus3/hunt-stokes/clean_bipo/tagged_ntuples/bisMSB_BHT_runs" + std::to_string(isotope) + "/output_bisMSB_" + std::to_string(bisMSB) + "_" + std::to_string(run_num) + ".root";
  // std::string output_name = "/data/snoplus3/hunt-stokes/clean_bipo/tagged_ntuples/bht_check_runs214/output_" + std::to_string(run_num) + ".root";
  TFile *g = new TFile(output_name.c_str(), "RECREATE");
  std::string ntuple_name = "ntuple" + std::to_string(isotope);
  // TNtuple *ntuple = new TNtuple(ntuple_name.c_str(), ntuple_name.c_str(), "dt:dr:x1:y1:z1:energy1:x2:y2:z2:energy2:gtid1:gtid2:nhits1:nhits2");
  // TNtuple *ntuple = new TNtuple(ntuple_name.c_str(), ntuple_name.c_str(), "dt:dr:x1:y1:z1:energy1:x2:y2:z2:energy2:gtid1:gtid2:nhits1:nhits2");
  double energy1, energy2, dR, dT, nhitsScaled1, nhitsScaled2, x1, x2, y1, y2, z1, z2;
  float HC_ratio1, HC_ratio2;
  int nhits1, nhits2, nhitsRaw1, nhitsRaw2, gtid1, gtid2;
//   TTree *tree = new TTree(ntuple_name.c_str(), "tree to hold the output of tagging");
  TTree *tree_6m          = new TTree("6m", "Data Information");
  TTree *tree_5p5m        = new TTree("5p5m", "Data Information");
  TTree *tree_5m          = new TTree("5m", "Data Information");
  TTree *tree_4p5m        = new TTree("4p5m", "Data Information");
  TTree *tree_4m          = new TTree("4m", "Data Information");
  TTree *tree_3p5m        = new TTree("3p5m", "Data Information");
  TTree *tree_3m          = new TTree("3m", "Data Information");
  
  tree_6m->Branch("energy_po", &energy1);
  tree_6m->Branch("nhits_po", &nhits1);
  tree_6m->Branch("nhits_bi", &nhits2);
  tree_6m->Branch("nhitsRaw_po", &nhitsRaw1);
  tree_6m->Branch("nhitsRaw_bi", &nhitsRaw2);
  tree_6m->Branch("nhitsScaled_po", &nhitsScaled1);
  tree_6m->Branch("nhitsScaled_bi", &nhitsScaled2);
  tree_6m->Branch("energy_bi", &energy2);
  tree_6m->Branch("gtid_po", &gtid1);
  tree_6m->Branch("gtid_bi", &gtid2);
  tree_6m->Branch("dR", &dR);
  tree_6m->Branch("dT", &dT);
  tree_6m->Branch("x_po", &x1);
  tree_6m->Branch("y_po", &y1);
  tree_6m->Branch("z_po", &z1);
  tree_6m->Branch("x_bi", &x2);
  tree_6m->Branch("y_bi", &y2);
  tree_6m->Branch("z_bi", &z2);

  tree_5p5m->Branch("energy_po", &energy1);
  tree_5p5m->Branch("nhits_po", &nhits1);
  tree_5p5m->Branch("nhits_bi", &nhits2);
  tree_5p5m->Branch("nhitsRaw_po", &nhitsRaw1);
  tree_5p5m->Branch("nhitsRaw_bi", &nhitsRaw2);
  tree_5p5m->Branch("nhitsScaled_po", &nhitsScaled1);
  tree_5p5m->Branch("nhitsScaled_bi", &nhitsScaled2);
  tree_5p5m->Branch("energy_bi", &energy2);
  tree_5p5m->Branch("gtid_po", &gtid1);
  tree_5p5m->Branch("gtid_bi", &gtid2);
  tree_5p5m->Branch("dR", &dR);
  tree_5p5m->Branch("dT", &dT);
  tree_5p5m->Branch("x_po", &x1);
  tree_5p5m->Branch("y_po", &y1);
  tree_5p5m->Branch("z_po", &z1);
  tree_5p5m->Branch("x_bi", &x2);
  tree_5p5m->Branch("y_bi", &y2);
  tree_5p5m->Branch("z_bi", &z2);

  tree_5m->Branch("energy_po", &energy1);
  tree_5m->Branch("nhits_po", &nhits1);
  tree_5m->Branch("nhits_bi", &nhits2);
  tree_5m->Branch("nhitsRaw_po", &nhitsRaw1);
  tree_5m->Branch("nhitsRaw_bi", &nhitsRaw2);
  tree_5m->Branch("nhitsScaled_po", &nhitsScaled1);
  tree_5m->Branch("nhitsScaled_bi", &nhitsScaled2);
  tree_5m->Branch("energy_bi", &energy2);
  tree_5m->Branch("gtid_po", &gtid1);
  tree_5m->Branch("gtid_bi", &gtid2);
  tree_5m->Branch("dR", &dR);
  tree_5m->Branch("dT", &dT);
  tree_5m->Branch("x_po", &x1);
  tree_5m->Branch("y_po", &y1);
  tree_5m->Branch("z_po", &z1);
  tree_5m->Branch("x_bi", &x2);
  tree_5m->Branch("y_bi", &y2);
  tree_5m->Branch("z_bi", &z2);

  tree_4p5m->Branch("energy_po", &energy1);
  tree_4p5m->Branch("nhits_po", &nhits1);
  tree_4p5m->Branch("nhits_bi", &nhits2);
  tree_4p5m->Branch("nhitsRaw_po", &nhitsRaw1);
  tree_4p5m->Branch("nhitsRaw_bi", &nhitsRaw2);
  tree_4p5m->Branch("nhitsScaled_po", &nhitsScaled1);
  tree_4p5m->Branch("nhitsScaled_bi", &nhitsScaled2);
  tree_4p5m->Branch("energy_bi", &energy2);
  tree_4p5m->Branch("gtid_po", &gtid1);
  tree_4p5m->Branch("gtid_bi", &gtid2);
  tree_4p5m->Branch("dR", &dR);
  tree_4p5m->Branch("dT", &dT);
  tree_4p5m->Branch("x_po", &x1);
  tree_4p5m->Branch("y_po", &y1);
  tree_4p5m->Branch("z_po", &z1);
  tree_4p5m->Branch("x_bi", &x2);
  tree_4p5m->Branch("y_bi", &y2);
  tree_4p5m->Branch("z_bi", &z2);

  tree_4m->Branch("energy_po", &energy1);
  tree_4m->Branch("nhits_po", &nhits1);
  tree_4m->Branch("nhits_bi", &nhits2);
  tree_4m->Branch("nhitsRaw_po", &nhitsRaw1);
  tree_4m->Branch("nhitsRaw_bi", &nhitsRaw2);
  tree_4m->Branch("nhitsScaled_po", &nhitsScaled1);
  tree_4m->Branch("nhitsScaled_bi", &nhitsScaled2);
  tree_4m->Branch("energy_bi", &energy2);
  tree_4m->Branch("gtid_po", &gtid1);
  tree_4m->Branch("gtid_bi", &gtid2);
  tree_4m->Branch("dR", &dR);
  tree_4m->Branch("dT", &dT);
  tree_4m->Branch("x_po", &x1);
  tree_4m->Branch("y_po", &y1);
  tree_4m->Branch("z_po", &z1);
  tree_4m->Branch("x_bi", &x2);
  tree_4m->Branch("y_bi", &y2);
  tree_4m->Branch("z_bi", &z2);

  tree_3p5m->Branch("energy_po", &energy1);
  tree_3p5m->Branch("nhits_po", &nhits1);
  tree_3p5m->Branch("nhits_bi", &nhits2);
  tree_3p5m->Branch("nhitsRaw_po", &nhitsRaw1);
  tree_3p5m->Branch("nhitsRaw_bi", &nhitsRaw2);
  tree_3p5m->Branch("nhitsScaled_po", &nhitsScaled1);
  tree_3p5m->Branch("nhitsScaled_bi", &nhitsScaled2);
  tree_3p5m->Branch("energy_bi", &energy2);
  tree_3p5m->Branch("gtid_po", &gtid1);
  tree_3p5m->Branch("gtid_bi", &gtid2);
  tree_3p5m->Branch("dR", &dR);
  tree_3p5m->Branch("dT", &dT);
  tree_3p5m->Branch("x_po", &x1);
  tree_3p5m->Branch("y_po", &y1);
  tree_3p5m->Branch("z_po", &z1);
  tree_3p5m->Branch("x_bi", &x2);
  tree_3p5m->Branch("y_bi", &y2);
  tree_3p5m->Branch("z_bi", &z2);

  tree_3m->Branch("energy_po", &energy1);
  tree_3m->Branch("nhits_po", &nhits1);
  tree_3m->Branch("nhits_bi", &nhits2);
  tree_3m->Branch("nhitsRaw_po", &nhitsRaw1);
  tree_3m->Branch("nhitsRaw_bi", &nhitsRaw2);
  tree_3m->Branch("nhitsScaled_po", &nhitsScaled1);
  tree_3m->Branch("nhitsScaled_bi", &nhitsScaled2);
  tree_3m->Branch("energy_bi", &energy2);
  tree_3m->Branch("gtid_po", &gtid1);
  tree_3m->Branch("gtid_bi", &gtid2);
  tree_3m->Branch("dR", &dR);
  tree_3m->Branch("dT", &dT);
  tree_3m->Branch("x_po", &x1);
  tree_3m->Branch("y_po", &y1);
  tree_3m->Branch("z_po", &z1);
  tree_3m->Branch("x_bi", &x2);
  tree_3m->Branch("y_bi", &y2);
  tree_3m->Branch("z_bi", &z2);

  std::cout << "Number of files found: " << filelist.size() << std::endl;

  // create a TChain that contains the entries of all the events in this run
  TChain *chain = new TChain("output");
  for (int i=0; i < filelist.size(); i++){
    std::string input_file  = filelist[i];
    chain->AddFile(input_file.c_str());
  }
  std::cout << "TChain contains " << chain->GetEntries() << " entries." << std::endl;
  
  // TFile *f = TFile::Open(input_file.c_str());
  // TNtuple *T = (TNtuple *) f->Get("output");
  
  double posX, posY, posZ, energy, nhitsScaled;
  int Nhits, NhitsRaw, mcIndex, ev, parent1, parent2, gtid;
  bool fit, fitValid;
  ULong64_t clock50, DCapplied, DCflagged;

  chain->SetBranchAddress("posx", &posX);
  chain->SetBranchAddress("posy", &posY);
  chain->SetBranchAddress("posz", &posZ);
  chain->SetBranchAddress("clockCount50", &clock50);
  chain->SetBranchAddress("nhitsCleaned", &Nhits);
  chain->SetBranchAddress("nhits", &NhitsRaw);
  chain->SetBranchAddress("dcApplied", &DCapplied);
  chain->SetBranchAddress("dcFlagged", &DCflagged);
  chain->SetBranchAddress("scintFit", &fit);
  chain->SetBranchAddress("fitValid", &fitValid);
  chain->SetBranchAddress("mcIndex", &mcIndex);
  chain->SetBranchAddress("parentpdg1", &parent1);
  chain->SetBranchAddress("parentpdg2", &parent2);
  chain->SetBranchAddress("energy", &energy);
  chain->SetBranchAddress("evIndex", &ev);
  chain->SetBranchAddress("eventID", &gtid);
  chain->SetBranchAddress("correctedNhits", &nhitsScaled);

  
  // loop over the entries for the BiPo214 analysis
  // note that in MC we have multiple EVs in each entry (re-triggers)
  int numEntries = chain->GetEntries();
  int numMuons   = chain->GetEntries("(dcFlagged&0x80)!=0x80");
  int numHighE   = chain->GetEntries("(nhitsCleaned>5000)");
  std::cout << "Total muons + highE events in file: " << numMuons << " + " << numHighE << " = " << numMuons+numHighE << std::endl;
  // num_vetos += numMuons;

  // vector stores the i and j index of the events already paired up to remove duplicate pairings
  std::vector<int> paired; 

  totalEntries += numEntries;
  std::cout << "Num. Entries: " << numEntries << std::endl;

  // deadtime is equal to (Number vetos x 20 s) + pileupTime + lone follower time
  bool veto_flag          = false;
  bool lone_follower_flag = false;
  double loneFollowerTime = 0;
  double pileupTime       = 0;
  // int num_vetos           = 0;
  ULong64_t veto_start_time, lone_start_time;

  for (int iEntry = 0; iEntry < chain->GetEntries(); iEntry++){
    chain->GetEntry(iEntry);

    // if (gtid == 13282390){
    //   std::cout << "\nPo-Wei's Event: " << gtid << std::endl;
    //   std::cout << "NhitsCleaned: " << Nhits << std::endl;
    //   std::cout << "Muon Flag: " << ((DCflagged&0x80)!=0x80) << std::endl;
    //   std::cout << "Veto Flag Status: " << veto_flag << std::endl;
    // }
    // check if event is high energy or tagged as a muon by DC
    if (Nhits > 5000 or ((DCflagged&0x80) != 0x80)){
        std::cout << "\nDeadtime Triggered - GTID: " << gtid << std::endl;
        std::cout << "NhitsCleaned: " << Nhits << std::endl;
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
    // check the reconstruction is valid (so we can rely on the energy reconstruction for highE veto)
    if (fit == false or fitValid == false){
      failedFit += 1;
      continue;
    }

    // looping backwards in time so skip first entry
    if (iEntry == 0){
      continue;
    }      

    // int nhits1, nhitsRaw1, gtid1;
    nhits1    = Nhits;
    nhitsScaled1 = nhitsScaled;
    nhitsRaw1 = NhitsRaw;
    gtid1     = gtid;

    // double x1, y1, z1, rho1, r1;
    double rho1, r1;
    ULong64_t flagged1; 
    flagged1 = DCflagged;
    energy1 = energy;
    
    if (((DCapplied & 0x2100000042C2) & DCflagged ) != (DCapplied & 0x2100000042C2)) continue;

    x1 = posX;
    y1 = posY;
    z1 = posZ - zOff;  // account for AV offset
    rho1 = pow(x1*x1 + y1*y1, 0.5);
    r1 = pow(x1*x1 + y1*y1 + z1*z1, .5);

    // check event is inside the FV
    if (r1 > FV_CUT){
      continue;
    }

    // energy cuts
    if (energy1 < PO_LOWER || energy1 > PO_UPPER){
      failedEnergy++;
      continue;
    }
    
    // check hit cleaning cleanliness cut
    // float HC_ratio1 = float(nhits1) / float(nhitsRaw1);
    HC_ratio1 = float(nhits1) / float(nhitsRaw1);
    // if (HC_ratio1 < 0.8){
    //   failedCleaning++;
    //   continue;
    // }
    ULong64_t clock1 = clock50;
    
    // all cuts passed and we've got a Po candidate; loop backwards to find Bi
    for (int jEntry = iEntry - 1; jEntry >= 0; jEntry--){
        chain->GetEntry(jEntry);
        ULong64_t clock2;
        clock2 = clock50;
        
        // apply dT cut
        // double dT = (clock1 - clock2) * 20; // dT in ns 
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
      
      // double x2, y2, z2, rho2, r2, energy2;
      double rho2, r2;
      // int nhits2, nhitsRaw2, gtid2;
      nhits2       = Nhits;
      nhitsRaw2    = NhitsRaw;
      nhitsScaled2 = nhitsScaled;
      gtid2        = gtid;
      
      
      x2 = posX;
      y2 = posY;
      z2 = posZ - zOff;
      rho2 = pow(x2*x2 + y2*y2, 0.5);
      r2 = pow(x2*x2 + y2*y2 + z2*z2, .5);

      // check event is inside the FV
      if (r2 > FV_CUT){
        continue;
      }

      energy2 = energy;

      // check if nhit cuts satisfied for Bi214 candidate
      if (energy2 < BI_LOWER || energy2 > BI_UPPER){
        continue;
      }

      // apply dR cut
      // double dR = pow(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2), 0.5);
      dR = pow(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2), 0.5);
      if (dR > dR_cut){
        continue;
      }
      // check hit cleaning cleanliness cut
      // float HC_ratio2 = float(nhits2) / float(nhitsRaw2);
      HC_ratio2 = float(nhits2) / float(nhitsRaw2);
      // if (HC_ratio2 < 0.8){
      //   continue;
      // }
      // all cuts passed and should be a unique BiPo214 pair
      paired.push_back(iEntry);
      paired.push_back(jEntry);

      // ntuple->Fill(dT, dR, x1, y1, z1, energy1, x2, y2, z2, energy2, gtid1, gtid2, HC_ratio1, HC_ratio2);
      // ntuple->Fill(dT, dR, x1, y1, z1, energy1, x2, y2, z2, energy2, gtid1, gtid2, nhits1, nhits2);
      // tree->Fill();
      // check if in 6m volume
        if (r2 <= 6000){
            numTagged_6m += 1;
            tree_6m->Fill();

            if (r2 <= 5500){
              tree_5p5m->Fill();

              if (r2 <= 5000){
                  tree_5m->Fill();

                  if (r2 <= 4500){
                      tree_4p5m->Fill();

                      if (r2 <= 4000){
                          tree_4m->Fill();

                          if (r2 <= 3500){
                              tree_3p5m->Fill();

                              if (r2 <= 3000){
                                  tree_3m->Fill();
                              }
                          }
                      }
                  }
              }
            }
          }

      // exit Bi214 candidate search as we have already paired it up
      break; 
    }
  }

  // save the output NTUPLE for this run
  g->cd();
  // ntuple->Write();
  tree_6m->Write();
  tree_5p5m->Write();
  tree_5m->Write();
  tree_4p5m->Write();
  tree_4m->Write();
  tree_3p5m->Write();
  tree_3m->Write();

  g->Close();
  double deadtime =  (20 * pow(10, 9) * num_vetos) + (pileupTime + loneFollowerTime); // in ns
  std::cout << "PileupTime: " << pileupTime << std::endl;
  std::cout << "Num. Tagged R <= 6m: " << numTagged_6m << std::endl;
  std::cout << "Num. Tagged R <= 4m: " << numTagged_4m << std::endl;
  std::cout << "Num. Tagged inside Donut: " << numTagged_donut << std::endl;
  std::cout << "Num. highE vetos: " << num_vetos << std::endl;
  std::cout << "highE deadtime: " << deadtime * pow(10, -9) << " s" << std::endl;

  std::cout << "% failed fit: " << float(failedFit) / float(totalEntries) << std::endl;
  std::cout << "% failed energy: " << float(failedEnergy) / float(totalEntries) << std::endl;
  std::cout << "% failed cleanliness: " << float(failedCleaning)  / float(totalEntries) << std::endl;  

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
  int         isotope     = std::stoi(argv[2]);
  
  // check if run appears on tag list (ie physics > 30 minutes)
  // std::ifstream runList("/home/hunt-stokes/multisite_analysis/simulated_list.txt");
  // std::ifstream runList("/home/hunt-stokes/multisite_analysis/directionality_list.txt");
  std::vector<int> tagList;
  int iRun;
  // while (runList >> iRun){
  //   tagList.push_back(iRun);
  // }
//   if (std::count(tagList.begin(), tagList.end(), std::stoi(runNum))){
//     std::cout << "Physics run, t > 30 mins. Proceeding ..." << std::endl;
//   } else {
//     std::cout << "Run not found on tag list! Will not attempt tagging." << std::endl;
//     return 0;
//   }

  std::string inFile1;
  std::string data_path = "/data/snoplus3/SNOplusData/processing/fullFill/rat-7.0.8/ntuples";
  if (std::stoi(runNum) < 307613){
    inFile1 = data_path + "/" + "Analysis20R_r0000";
    // std::cout << "Loading: " << inFile1 << std::endl;
  }
  if (std::stoi(runNum) >= 307613){
      inFile1 = data_path + "/" + "Analysis20_r0000";
  }
  
  // std::string inFile1 = "/data/snoplus3/SNOplusData/production/miniProd_RAT-7-0-14_ASCI_RATHS_newRecoordination/ntuple/simulationBiPo" + std::to_string(isotope) + "_";
  // std::string inFile1 = "/data/snoplus3/hunt-stokes/bipo_data/bipo_ntuples8/Analysis20_r0000";
  // std::string inFile1 = "/data/snoplus3/SNOplusData/processing/fullFill/rat-7.0.8/ntuples/Analysis20_r0000";
  std::string inFile2 = "*";
  const std::string inFileName = inFile1+runNum+inFile2;
  std::cout << inFileName << std::endl;
  std::vector<std::string> filelist = globVector(inFileName);

  // if we find no subrun files, skip the run!
  if (filelist.empty()){
    std::cout << "No valid files!" << std::endl;
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
    start_time = start_day * 24 * 3600 * pow(10, 9) + start_sec * pow(10, 9) + start_nsec;
    end_time = end_day * 24 * 3600 * pow(10, 9) + end_sec * pow(10, 9) + end_nsec;
    duration = end_time - start_time;
    std::cout << "Duration of run: " << duration * pow(10, -9) << " s" << std::endl;

    double deadtime = bipo_tagger_full(filelist, zOff, runID, isotope);
    // double adjusted_livetime = duration - deadtime; 

    // std::cout << "Adjusted Livetime: " << adjusted_livetime * pow(10, -9) << " s" << std::endl;
    
    // create output txt file containing the livetime
    std::ofstream livetime;
    livetime.open("/data/snoplus3/hunt-stokes/clean_multisite/bipo_214_cleanSelec3/completed/" + runNum + ".txt");
    livetime << "1";
    livetime.close();
    
    return 0;
  }
}
