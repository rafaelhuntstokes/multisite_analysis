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
  float pileupT      = 0;        // if a highE event resets clock within a deadtime window, find the residual time
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
  // implementing a GENERAL coincidence cut to get rid of (alpha, n), BiPo214 and BiPo212 together ...
  BI_LOWER = 0;
  BI_UPPER = 10;
  PO_LOWER = 0.2;
  PO_UPPER = 10;
  dR_cut   = 2000;
  dT_low   = 0;
  dT_high  = 4000000;
  FV_CUT   = 6000;

  // variables which handle the muon and high energy vetos
  ULong64_t cut_trig_time = 0;        // the event time according to 50 MHz clock that starts muon veto
  bool cut_status         = false;    // flag which says whether or not to skip event from consideration
  std::string output_name;

  output_name = "/data/snoplus3/hunt-stokes/clean_multisite/general_cut_bipo214/output_" + std::to_string(isotope) + "_" + std::to_string(run_num) + ".root";
  
  TFile *g = new TFile(output_name.c_str(), "RECREATE");
  std::string ntuple_name = "ntuple" + std::to_string(isotope);
  double energy1, energy2, dR, dT, nhitsScaled1, nhitsScaled2, x1, x2, y1, y2, z1, z2, itr;
  float HC_ratio1, HC_ratio2;
  int nhits1, nhits2, nhitsRaw1, nhitsRaw2, gtid1, gtid2;
  
  int num_unpaired_in_roi  = 0; // number of unpaired prompt events reconstructing into ROI
  int num_unpaired_removed = 0; // number of unpaired prompt events reconstructing into ROI removed by ITR cuts

  TTree *prompt_6m = new TTree("6m_prompt", "general_cut");
  TTree *prompt_5p5m = new TTree("5p5m_prompt", "general_cut");
  TTree *prompt_5m = new TTree("5m_prompt", "general_cut");
  TTree *prompt_4p5m = new TTree("4p5m_prompt", "general_cut");
  TTree *prompt_4m = new TTree("4m_prompt", "general_cut");
  TTree *prompt_3p5m = new TTree("3p5m_prompt", "general_cut");
  TTree *prompt_3m = new TTree("3m_prompt", "general_cut");

  TTree *delayed_6m = new TTree("6m_delayed", "general_cut");
  TTree *delayed_5p5m = new TTree("5p5m_delayed", "general_cut");
  TTree *delayed_5m = new TTree("5m_delayed", "general_cut");
  TTree *delayed_4p5m = new TTree("4p5m_delayed", "general_cut");
  TTree *delayed_4m = new TTree("4m_delayed", "general_cut");
  TTree *delayed_3p5m = new TTree("3p5m_delayed", "general_cut");
  TTree *delayed_3m = new TTree("3m_delayed", "general_cut");

  prompt_6m->Branch("prompt_energy_inside_ROI", &energy2);
  delayed_6m->Branch("delayed_energy_inside_ROI", &energy1);
  prompt_5p5m->Branch("prompt_energy_inside_ROI", &energy2);
  delayed_5p5m->Branch("delayed_energy_inside_ROI", &energy1);
  prompt_5m->Branch("prompt_energy_inside_ROI", &energy2);
  delayed_5m->Branch("delayed_energy_inside_ROI", &energy1);
  prompt_4p5m->Branch("prompt_energy_inside_ROI", &energy2);
  delayed_4p5m->Branch("delayed_energy_inside_ROI", &energy1);
  prompt_4m->Branch("prompt_energy_inside_ROI", &energy2);
  delayed_4m->Branch("delayed_energy_inside_ROI", &energy1);
  prompt_3p5m->Branch("prompt_energy_inside_ROI", &energy2);
  delayed_3p5m->Branch("delayed_energy_inside_ROI", &energy1);
  prompt_3m->Branch("prompt_energy_inside_ROI", &energy2);
  delayed_3m->Branch("delayed_energy_inside_ROI", &energy1);

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
    T->SetBranchAddress("itr", &itr);

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
      double rho1, r1, itr1;
      itr1 = itr;
      ULong64_t flagged1; 
      flagged1 = DCflagged;
      energy1 = energy;

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
          failedFit += 1;
          continue;
        }
        if (fitValid == false){
          failedFit += 1;
          continue;
        }
        ULong64_t flagged2; 
        flagged2 = DCflagged;
        
        double rho2, r2, itr2;
        itr2 = itr;
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
        dR = pow(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2), 0.5);
        if (dR > dR_cut){
          continue;
        }
        
        // all cuts passed and should be a unique BiPo214 pair
        paired.push_back(iEntry);
        paired.push_back(jEntry);

        bool prompt_roi  = false;
        bool delayed_roi = false;
        if (energy2 <= 6.0 and energy2 >= 2.0){
          prompt_roi = true;
        }
        if (energy1 <= 6.0 and energy1 >= 2.0){
          delayed_roi = true;
        }
        std::cout << "Delayed Energy: " << energy1 << std::endl;
        // check if in 6m volume
          if (r2 <= 6000){
              numTagged_6m += 1;

              if (prompt_roi == true){
                prompt_6m->Fill();
              }
              if (delayed_roi == true){
                delayed_6m->Fill();
              }

              if (r2 <= 5500){
                
                if (prompt_roi == true){
                  prompt_5p5m->Fill();
                }
                if (delayed_roi == true){
                  delayed_5p5m->Fill();
                }

                if (r2 <= 5000){

                    if (prompt_roi == true){
                      prompt_5m->Fill();
                    }
                    if (delayed_roi == true){
                      delayed_5m->Fill();
                    }

                    if (r2 <= 4500){

                        if (prompt_roi == true){
                          prompt_4p5m->Fill();
                        }
                        if (delayed_roi == true){
                          delayed_4p5m->Fill();
                        }

                        if (r2 <= 4000){

                            if (prompt_roi == true){
                              prompt_4m->Fill();
                            }
                            if (delayed_roi == true){
                              delayed_4m->Fill();
                            }

                            if (r2 <= 3500){

                                if (prompt_roi == true){
                                  prompt_3p5m->Fill();
                                }
                                if (delayed_roi == true){
                                  delayed_3p5m->Fill();
                                }

                                if (r2 <= 3000){

                                    if (prompt_roi == true){
                                      prompt_3m->Fill();
                                    }
                                    if (delayed_roi == true){
                                      delayed_3m->Fill();
                                    }
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

  // check how many prompts are left UNPAIRED and how many get removed by the ITR cut
  for (int ientry = 0; ientry < T->GetEntries(); ientry++ ){
    
    T->GetEntry(ientry);

    // check it's a prompt
    if (ev != 0){
      // not a prompt event
      continue;
    }

    // check the prompt is not paired up
     bool pair = false;
        for (int iPair = 0; iPair < paired.size(); iPair++){
          // loop over every entry in iPair and see if already paired
          
          if (ientry == paired[iPair]){
            pair = true;
            break;
          }
        }
        if (pair == true){
          // this bismuth has already been paired up with something! skip.
          continue;
        }
        // else it hasn't been paired already and we're free to continue our checks

    // check that the prompt reconstructs into the ROI
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
    
    
    nhits2       = Nhits;
    nhitsRaw2    = NhitsRaw;
    nhitsScaled2 = nhitsScaled;
    gtid2        = gtid;
    
    
    x2 = posX;
    y2 = posY;
    z2 = posZ - zOff;
    double r2 = pow(x2*x2 + y2*y2 + z2*z2, .5);

    // check event is inside the FV
    if (r2 > 4500){
      continue;
    }

    // check if nhit cuts satisfied for Bi214 candidate
    if (energy < BI_LOWER || energy > BI_UPPER){
      continue;
    }

    // we now have an UNPAIRED prompt that reconstructs into our ROI --> increment counter
    num_unpaired_in_roi++;
    std::cout << itr << std::endl;
    // check if itr cuts remove prompt event
    if (itr > 0.3 or itr < 0.18){
      num_unpaired_removed++;
    }
  }
  
  std::cout << "Num. Unpaired Prompts in ROI: " << num_unpaired_in_roi << std::endl;
  std::cout << "Num. Unpaired Prompts in ROI removed by ITR: " << num_unpaired_removed << std::endl;
  std::cout << "Fraction of unpaired prompts in ROI removed by ITR cuts: " << num_unpaired_removed / num_unpaired_in_roi << std::endl;
  }
  
  // create output txt file containing the number of unpaired prompts remaining and the number of removed prompts with ITR cuts
  std::ofstream itr_info;
  itr_info.open("/data/snoplus3/hunt-stokes/clean_multisite/itr_cuts_bipo214/" + std::to_string(run_num) + ".txt");
  itr_info << num_unpaired_in_roi << std::endl;
  itr_info << num_unpaired_removed << std::endl;
  itr_info.close();
  

  // save the output NTUPLE for this run
  g->cd();
  prompt_6m->Write();
  prompt_5p5m->Write();
  prompt_5m->Write();
  prompt_4p5m->Write();
  prompt_4m->Write();
  prompt_3p5m->Write();
  prompt_3m->Write();

  delayed_6m->Write();
  delayed_5p5m->Write();
  delayed_5m->Write();
  delayed_4p5m->Write();
  delayed_4m->Write();
  delayed_3p5m->Write();
  delayed_3m->Write();

  g->Close();
  double deadtime =  (20 * pow(10, 9) * num_vetos) + (pileupT * pow(10, 9)); // in ns
  std::cout << "Num. Tagged R <= 6m: " << numTagged_6m << std::endl; 

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
  std::ifstream runList("/home/hunt-stokes/multisite_analysis/simulated_list.txt");
  std::vector<int> tagList;
  int iRun;
  while (runList >> iRun){
    tagList.push_back(iRun);
  }

  std::string inFile1;
  // std::string data_path = "/data/snoplus3/SNOplusData/production/miniProd_RAT-7-0-14_ASCI_RATHS_newRecoordination/ntuple/simulationBiPo" + std::to_string(isotope) + "_";
  std::string data_path = "/data/snoplus3/SNOplusData/production/miniProd_RAT-7-0-14_ASCI_RATHS_newRecoordination/ntuple/simulationBiPo214_";
  // std::string data_path = "/data/snoplus3/SNOplusData/production/miniProd_RAT-7-0-14_ASCI_RATHS_newRecoordination/ntuple/simulationAlphaN_LAB_13C_";
  // if (std::stoi(runNum) < 307613){
  //   inFile1 = data_path + "/" + "Analysis20R_r0000";
  //   // std::cout << "Loading: " << inFile1 << std::endl;
  // }
  // if (std::stoi(runNum) >= 307613){
  //     inFile1 = data_path + "/" + "Analysis20_r0000";
  // }
  inFile1 = data_path;
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
    return 0;
  }
}
