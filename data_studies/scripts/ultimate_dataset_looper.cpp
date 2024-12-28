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

/*
This script is designed to speed up the analysis by limiting the amount of I/O condor has to handle. Particualrly, looping through the
dataset to do coincidence tagging and extract the energy spectrum before / after cuts.

A function for coincidence tagging and energy spectrum + livetime calculations is included.

For each coincidence tagging and each run, energy spectra of the tags, position, dR and dT distributions are saved in a TTree. The run-by-run
event GTIDs are also saved to a .txt file.

For the total energy spectrum, each run returns the deadtime subtracted livetime, which is saved in an individual text file. A TTree containing
the energy spectrum, positions, ITR and biPolikelihood212 and bipolikelihood214 (one day I will cut out the in-windows with this!).

Two TTrees are returned from the energy spectrum code: 1. the full energy spectrum without cuts; 2. the energy spectrum with GTIDs flagged by
coincidence tagging removed.
*/

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

void coincidence_tagger(TChain *chain, double AV_OFFSET, int run_num, std::vector<int> &gtid_vector){
    /*
    Function runs over the events in the TChain and tags coincidence pairs. 

    Using a general tagging logic to remove all kinds of coincidence events.
    */
    
    // setup the variables to hold event information
    bool fit, fitValid;
    int nhitsCleaned, gtid;
    ULong64_t clock50, DCapplied, DCflagged;
    Double_t ITRVal, energy, posX, posY, posZ, correctedNhits, inWindow212LL, inWindow214LL;

    chain->SetBranchAddress("posx", &posX);
    chain->SetBranchAddress("posy", &posY);
    chain->SetBranchAddress("posz", &posZ);
    chain->SetBranchAddress("clockCount50", &clock50);
    chain->SetBranchAddress("correctedNhits", &correctedNhits);
    chain->SetBranchAddress("nhitsCleaned", &nhitsCleaned);
    chain->SetBranchAddress("scintFit", &fit);
    chain->SetBranchAddress("fitValid", &fitValid);
    chain->SetBranchAddress("energy", &energy);
    chain->SetBranchAddress("dcApplied", &DCapplied);
    chain->SetBranchAddress("dcFlagged", &DCflagged);
    chain->SetBranchAddress("alphaBeta212", &inWindow212LL);
    chain->SetBranchAddress("alphaBeta214", &inWindow214LL);
    chain->SetBranchAddress("itr", &ITRVal);
    chain->SetBranchAddress("eventID", &gtid);

    // variables for applied cuts
    double DELTA_T_LOWER, DELTA_T_HIGHER, DELTA_R, FV_CUT, PROMPT_LOWER, PROMPT_HIGHER, DELAYED_LOWER, DELAYED_HIGHER;

    // create an output .root for this run
    std::string output_name;
    
    // define the general coincidence cuts to run over the dataset
    PROMPT_LOWER   = 0;
    PROMPT_HIGHER  = 10;
    DELAYED_LOWER  = 0.2;
    DELAYED_HIGHER = 10;
    DELTA_R        = 2000;
    DELTA_T_LOWER  = 0;
    DELTA_T_HIGHER = 4000000;
    FV_CUT         = 6000;

    output_name = "path_to_coincidence_cut_info";//"/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs_general_coincidence/output_" + std::to_string(run_num) + ".root";


    // define the variables for tracking the coincidence pair attributes
    double x1, x2, y1, y2, z1, z2, energy1, energy2, itr1, itr2, nhitsCorrected1, nhitsCorrected2, dR, dT;
    float inWindowLL212_1, inWindowLL212_2, inWindowLL214_1, inWindowLL214_2;
    int nhitsCleaned1, nhitsCleaned2, gtid1, gtid2;
    ULong64_t clock1, clock2;

    // create the output file
    TFile *g                = new TFile(output_name.c_str(), "RECREATE");
    // save a TTree inside the file
    TTree *tree             = new TTree("coincidence_pairs", "coincidence_pairs");
    // define branches to save all the relevant information about coincidence pairs
    tree->Branch("prompt_energy", &energy1);
    tree->Branch("delayed_energy", &energy2);
    tree->Branch("deltaR", &dR);
    tree->Branch("deltaT", &dT);
    tree->Branch("prompt_x", &x1);
    tree->Branch("prompt_y", &y1);
    tree->Branch("prompt_z", &z1);
    tree->Branch("delayed_x", &x2);
    tree->Branch("delayed_y", &y2);
    tree->Branch("delayed_z", &z2);
    tree->Branch("prompt_itr", &itr1);
    tree->Branch("delayed_itr", &itr2);
    tree->Branch("prompt_LL212", &inWindowLL212_1);
    tree->Branch("delayed_LL212", &inWindowLL212_2);
    tree->Branch("prompt_LL214", &inWindowLL214_1);
    tree->Branch("delayed_LL214", &inWindowLL214_2);
    tree->Branch("prompt_gtid", &gtid1);
    tree->Branch("delayed_gtid", &gtid2);

    // setup flags for the deadtime calculation (accounting for deadtime veto pileup)
    bool veto_flag          = false;
    bool lone_follower_flag = false;
    double loneFollowerTime = 0;
    double pileupTime       = 0;
    ULong64_t veto_start_time, lone_start_time;
    
    // vector stores the i and j index of the events already paired up to remove duplicate pairings
    std::vector<int> paired; 
    
    // begin loop over the DELAYED event
    for (int ientry = 0; ientry < chain->GetEntries(); ientry++){
        chain->GetEntry(ientry);
        
        // check if event is high energy or tagged as a muon by DC
        if (nhitsCleaned > 5000 or (DCflagged&0x80)!=0x80){
          
            // check if we are inside a lone muon follower veto window
            if (lone_follower_flag == true){
                // add the dT and switch off the muon follower veto
                ULong64_t deltaT = ((clock50-lone_start_time) & 0x7FFFFFFFFFF)*20.0;
                loneFollowerTime += deltaT;
                lone_follower_flag = false;
            }

            if (veto_flag == false){
                veto_flag = true;
                veto_start_time = clock50;
                continue;
            }
            else{
                // we have pileup! Need to calculate the additional pileup time.
                ULong64_t deltaT = ((clock50-veto_start_time) & 0x7FFFFFFFFFF)*20.0;
                pileupTime += deltaT;

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
        
        // looping backwards in time so skip first entry
        if (ientry == 0){
            continue;
        }

        // check data cleaning
        if (((DCapplied & 0x2100000042C2) & DCflagged ) != (DCapplied & 0x2100000042C2)) {
            continue;
        }

        // check reconstruction validity
        if (fit == false or fitValid == false){
            continue;
        }

        // define all the variables for the delayed event
        energy2         = energy;
        x2              = posX;
        y2              = posY;
        z2              = posZ - AV_OFFSET; // subtract the AV offset
        itr2            = ITRVal;
        clock2          = clock50;
        nhitsCleaned2   = nhitsCleaned;
        nhitsCorrected2 = correctedNhits;
        gtid2           = gtid;
        inWindowLL212_2 = inWindow212LL;
        inWindowLL214_2 = inWindow214LL;

        // calculate the radius
        double rSq = x2*x2 + y2*y2 + z2*z2;

        // apply delayed energy and FV cuts
        if (rSq > FV_CUT*FV_CUT){
            continue;
        }
        // std::cout << "Passed FV and DataCleaning!" << std::endl;
        if (energy2 < DELAYED_LOWER || energy2 > DELAYED_HIGHER){
            continue;
        }
        // std::cout << "Found delayed candidate!" << std::endl;
        // loop backwards in time to find the PROMPT event
        for (int jentry = ientry - 1; jentry >= 0; jentry--){
            chain->GetEntry(jentry);
        
            // apply DC cleaning and reconstruction cuts
            if (((DCapplied & 0x2100000042C2) & DCflagged ) != (DCapplied & 0x2100000042C2)) {
                continue;
            }

            // check reconstruction validity
            if (fit == false or fitValid == false){
                continue;
            }

            // define all variables for PROMPT event
            energy1         = energy;
            x1              = posX;
            y1              = posY;
            z1              = posZ - AV_OFFSET; // subtract the AV offset
            itr1            = ITRVal;
            clock1          = clock50;
            nhitsCleaned1   = nhitsCleaned;
            nhitsCorrected1 = correctedNhits;
            gtid1           = gtid;
            inWindowLL212_1 = inWindow212LL;
            inWindowLL214_1 = inWindow214LL;

            // check dT first as we want to break out of this loop fast if we go beyond the dT upper limit!
            dT = ((clock2 - clock1) & 0x7FFFFFFFFFF) * 20;
            if (dT > DELTA_T_HIGHER){
                break;
            }
            if (dT < DELTA_T_LOWER){
                continue;
            }

            // only proceed if the event is within the dT window ...

            // calculate the radius
            double rSq = x1*x1 + y1*y1 + z1*z1;

            // apply prompt energy and FV cuts
            if (rSq > FV_CUT*FV_CUT){
                continue;
            }

            if (energy1 < PROMPT_LOWER || energy1 > PROMPT_HIGHER){
                continue;
            }
            
            // apply the dR cut
            double dRsq = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2);
            if (dRsq > DELTA_R*DELTA_R){
                continue;
            }
            dR = pow(dRsq, 0.5);
            
            // check if the event has already been paired up
            bool pair = false;
            for (int iPair = 0; iPair < paired.size(); iPair++){
                // loop over every entry in iPair and see if already paired
          
                if (jentry == paired[iPair]){
                    pair = true;
                    break;
                }
            }
            
            if (pair == true){
                // this bismuth has already been paired up with something! skip.
                continue;
            }
            
            // all cuts passed and should be a unique BiPo214 pair
            paired.push_back(ientry);
            paired.push_back(jentry);

            // fill the output tree
            tree->Fill();

            // fill the vector of GTIDs
            gtid_vector.push_back(gtid1);
            gtid_vector.push_back(gtid2);

            // break out of the prompt candidate loop as we have paired up already
            // std::cout << "Tagged " << gtid_vector.size() << " Events." << std::endl;
            break;
        }
    }

    // write the tree to the file once the loops have finished
    g->cd();
    tree->Write();
    g->Close();
}

void extract_energy_spectrum(TChain *chain, int run_number, double AV_OFFSET, std::vector<int> &gtid_214, double duration){
    /*
    Function loops through the dataset and records the energy of every event that passes data
    cleaning, reconstruction checks and falls inside the ROI.

    The livetime of the run (adjusted for deadtime) is also saved.

    For each FV, two copies of the energy spectrum are returned:
        1. Containing all selected events
        2. Containing all selected events which are not tagged as part of a coincidence pair.
    */

    // ROI ENERGY CUTS
    int ENERGY_LOW  = 2.0;
    int ENERGY_HIGH = 6.0;

    // set the variables and branch addresses of the data and output
    bool fit, fitValid;
    int nhitsCleaned, gtid;
    ULong64_t clock50, DCapplied, DCflagged;
    Double_t ITRVal, energy, posX, posY, posZ, correctedNhits, inWindow212LL, inWindow214LL;

    chain->SetBranchAddress("posx", &posX);
    chain->SetBranchAddress("posy", &posY);
    chain->SetBranchAddress("posz", &posZ);
    chain->SetBranchAddress("clockCount50", &clock50);
    chain->SetBranchAddress("correctedNhits", &correctedNhits);
    chain->SetBranchAddress("nhitsCleaned", &nhitsCleaned);
    chain->SetBranchAddress("scintFit", &fit);
    chain->SetBranchAddress("fitValid", &fitValid);
    chain->SetBranchAddress("energy", &energy);
    chain->SetBranchAddress("dcApplied", &DCapplied);
    chain->SetBranchAddress("dcFlagged", &DCflagged);
    chain->SetBranchAddress("alphaBeta212", &inWindow212LL);
    chain->SetBranchAddress("alphaBeta214", &inWindow214LL);
    chain->SetBranchAddress("itr", &ITRVal);
    chain->SetBranchAddress("eventID", &gtid);

    // create the output file
    std::string output_name = "/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs_general_coincidence/output_" + std::to_string(run_number) + ".root";
    TFile *g                = new TFile(output_name.c_str(), "RECREATE");
    // save a TTree for full spectrum and clean spectrum, for each FV
    TTree *tree_full_6m          = new TTree("full_6m", "full_6m");
    TTree *tree_full_5p5m        = new TTree("full_5p5m", "full_5p5m");
    TTree *tree_full_5m          = new TTree("full_5m", "full_5m");
    TTree *tree_full_4p5m        = new TTree("full_4p5m", "full_4p5m");
    TTree *tree_full_4m          = new TTree("full_4m", "full_4m");
    TTree *tree_full_3p5m        = new TTree("full_3p5m", "full_3p5m");
    TTree *tree_full_3m          = new TTree("full_3m", "full_3m");
    TTree *tree_clean_6m         = new TTree("clean_6m", "clean_6m");
    TTree *tree_clean_5p5m       = new TTree("clean_5p5m", "clean_5p5m");
    TTree *tree_clean_5m         = new TTree("clean_5m", "clean_5m");
    TTree *tree_clean_4p5m       = new TTree("clean_4p5m", "clean_4p5m");
    TTree *tree_clean_4m         = new TTree("clean_4m", "clean_4m");
    TTree *tree_clean_3p5m       = new TTree("clean_3p5m", "clean_3p5m");
    TTree *tree_clean_3m         = new TTree("clean_3m", "clean_3m");
    
    // each tree contains the position, energy, ITR, and alphabeta in window classifer result
    tree_full_6m->Branch("energy", &energy);
    tree_full_6m->Branch("x", &posX);
    tree_full_6m->Branch("y", &posY);
    tree_full_6m->Branch("z", &posZ);
    tree_full_6m->Branch("nhitsClean", &nhitsCleaned);
    tree_full_6m->Branch("nhitsCorrected", &correctedNhits);
    tree_full_6m->Branch("ITR", &ITRVal);
    tree_full_6m->Branch("gtid", &gtid);
    tree_full_6m->Branch("runid", &run_number);
    tree_full_6m->Branch("alphaBeta212", &inWindow212LL);
    tree_full_6m->Branch("alphaBeta214", &inWindow214LL);
    tree_clean_6m->Branch("energy", &energy);
    tree_clean_6m->Branch("x", &posX);
    tree_clean_6m->Branch("y", &posY);
    tree_clean_6m->Branch("z", &posZ);
    tree_clean_6m->Branch("nhitsClean", &nhitsCleaned);
    tree_clean_6m->Branch("nhitsCorrected", &correctedNhits);
    tree_clean_6m->Branch("ITR", &ITRVal);
    tree_clean_6m->Branch("gtid", &gtid);
    tree_clean_6m->Branch("runid", &run_number);
    tree_clean_6m->Branch("alphaBeta212", &inWindow212LL);
    tree_clean_6m->Branch("alphaBeta214", &inWindow214LL);

    tree_full_5p5m->Branch("energy", &energy);
    tree_full_5p5m->Branch("x", &posX);
    tree_full_5p5m->Branch("y", &posY);
    tree_full_5p5m->Branch("z", &posZ);
    tree_full_5p5m->Branch("nhitsClean", &nhitsCleaned);
    tree_full_5p5m->Branch("nhitsCorrected", &correctedNhits);
    tree_full_5p5m->Branch("ITR", &ITRVal);
    tree_full_5p5m->Branch("gtid", &gtid);
    tree_full_5p5m->Branch("runid", &run_number);
    tree_full_5p5m->Branch("alphaBeta212", &inWindow212LL);
    tree_full_5p5m->Branch("alphaBeta214", &inWindow214LL);
    tree_clean_5p5m->Branch("energy", &energy);
    tree_clean_5p5m->Branch("x", &posX);
    tree_clean_5p5m->Branch("y", &posY);
    tree_clean_5p5m->Branch("z", &posZ);
    tree_clean_5p5m->Branch("nhitsClean", &nhitsCleaned);
    tree_clean_5p5m->Branch("nhitsCorrected", &correctedNhits);
    tree_clean_5p5m->Branch("ITR", &ITRVal);
    tree_clean_5p5m->Branch("gtid", &gtid);
    tree_clean_5p5m->Branch("runid", &run_number);
    tree_clean_5p5m->Branch("alphaBeta212", &inWindow212LL);
    tree_clean_5p5m->Branch("alphaBeta214", &inWindow214LL);

    tree_full_5m->Branch("energy", &energy);
    tree_full_5m->Branch("x", &posX);
    tree_full_5m->Branch("y", &posY);
    tree_full_5m->Branch("z", &posZ);
    tree_full_5m->Branch("nhitsClean", &nhitsCleaned);
    tree_full_5m->Branch("nhitsCorrected", &correctedNhits);
    tree_full_5m->Branch("ITR", &ITRVal);
    tree_full_5m->Branch("gtid", &gtid);
    tree_full_5m->Branch("runid", &run_number);
    tree_full_5m->Branch("alphaBeta212", &inWindow212LL);
    tree_full_5m->Branch("alphaBeta214", &inWindow214LL);
    tree_clean_5m->Branch("energy", &energy);
    tree_clean_5m->Branch("x", &posX);
    tree_clean_5m->Branch("y", &posY);
    tree_clean_5m->Branch("z", &posZ);
    tree_clean_5m->Branch("nhitsClean", &nhitsCleaned);
    tree_clean_5m->Branch("nhitsCorrected", &correctedNhits);
    tree_clean_5m->Branch("ITR", &ITRVal);
    tree_clean_5m->Branch("gtid", &gtid);
    tree_clean_5m->Branch("runid", &run_number);
    tree_clean_5m->Branch("alphaBeta212", &inWindow212LL);
    tree_clean_5m->Branch("alphaBeta214", &inWindow214LL);

    tree_full_4p5m->Branch("energy", &energy);
    tree_full_4p5m->Branch("x", &posX);
    tree_full_4p5m->Branch("y", &posY);
    tree_full_4p5m->Branch("z", &posZ);
    tree_full_4p5m->Branch("nhitsClean", &nhitsCleaned);
    tree_full_4p5m->Branch("nhitsCorrected", &correctedNhits);
    tree_full_4p5m->Branch("ITR", &ITRVal);
    tree_full_4p5m->Branch("gtid", &gtid);
    tree_full_4p5m->Branch("runid", &run_number);
    tree_full_4p5m->Branch("alphaBeta212", &inWindow212LL);
    tree_full_4p5m->Branch("alphaBeta214", &inWindow214LL);
    tree_clean_4p5m->Branch("energy", &energy);
    tree_clean_4p5m->Branch("x", &posX);
    tree_clean_4p5m->Branch("y", &posY);
    tree_clean_4p5m->Branch("z", &posZ);
    tree_clean_4p5m->Branch("nhitsClean", &nhitsCleaned);
    tree_clean_4p5m->Branch("nhitsCorrected", &correctedNhits);
    tree_clean_4p5m->Branch("ITR", &ITRVal);
    tree_clean_4p5m->Branch("gtid", &gtid);
    tree_clean_4p5m->Branch("runid", &run_number);
    tree_clean_4p5m->Branch("alphaBeta212", &inWindow212LL);
    tree_clean_4p5m->Branch("alphaBeta214", &inWindow214LL);

    tree_full_4m->Branch("energy", &energy);
    tree_full_4m->Branch("x", &posX);
    tree_full_4m->Branch("y", &posY);
    tree_full_4m->Branch("z", &posZ);
    tree_full_4m->Branch("nhitsClean", &nhitsCleaned);
    tree_full_4m->Branch("nhitsCorrected", &correctedNhits);
    tree_full_4m->Branch("ITR", &ITRVal);
    tree_full_4m->Branch("gtid", &gtid);
    tree_full_4m->Branch("runid", &run_number);
    tree_full_4m->Branch("alphaBeta212", &inWindow212LL);
    tree_full_4m->Branch("alphaBeta214", &inWindow214LL);
    tree_clean_4m->Branch("energy", &energy);
    tree_clean_4m->Branch("x", &posX);
    tree_clean_4m->Branch("y", &posY);
    tree_clean_4m->Branch("z", &posZ);
    tree_clean_4m->Branch("nhitsClean", &nhitsCleaned);
    tree_clean_4m->Branch("nhitsCorrected", &correctedNhits);
    tree_clean_4m->Branch("ITR", &ITRVal);
    tree_clean_4m->Branch("gtid", &gtid);
    tree_clean_4m->Branch("runid", &run_number);
    tree_clean_4m->Branch("alphaBeta212", &inWindow212LL);
    tree_clean_4m->Branch("alphaBeta214", &inWindow214LL);

    tree_full_3p5m->Branch("energy", &energy);
    tree_full_3p5m->Branch("x", &posX);
    tree_full_3p5m->Branch("y", &posY);
    tree_full_3p5m->Branch("z", &posZ);
    tree_full_3p5m->Branch("nhitsClean", &nhitsCleaned);
    tree_full_3p5m->Branch("nhitsCorrected", &correctedNhits);
    tree_full_3p5m->Branch("ITR", &ITRVal);
    tree_full_3p5m->Branch("gtid", &gtid);
    tree_full_3p5m->Branch("runid", &run_number);
    tree_full_3p5m->Branch("alphaBeta212", &inWindow212LL);
    tree_full_3p5m->Branch("alphaBeta214", &inWindow214LL);
    tree_clean_3p5m->Branch("energy", &energy);
    tree_clean_3p5m->Branch("x", &posX);
    tree_clean_3p5m->Branch("y", &posY);
    tree_clean_3p5m->Branch("z", &posZ);
    tree_clean_3p5m->Branch("nhitsClean", &nhitsCleaned);
    tree_clean_3p5m->Branch("nhitsCorrected", &correctedNhits);
    tree_clean_3p5m->Branch("ITR", &ITRVal);
    tree_clean_3p5m->Branch("gtid", &gtid);
    tree_clean_3p5m->Branch("runid", &run_number);
    tree_clean_3p5m->Branch("alphaBeta212", &inWindow212LL);
    tree_clean_3p5m->Branch("alphaBeta214", &inWindow214LL);

    tree_full_3m->Branch("energy", &energy);
    tree_full_3m->Branch("x", &posX);
    tree_full_3m->Branch("y", &posY);
    tree_full_3m->Branch("z", &posZ);
    tree_full_3m->Branch("nhitsClean", &nhitsCleaned);
    tree_full_3m->Branch("nhitsCorrected", &correctedNhits);
    tree_full_3m->Branch("ITR", &ITRVal);
    tree_full_3m->Branch("gtid", &gtid);
    tree_full_3m->Branch("runid", &run_number);
    tree_full_3m->Branch("alphaBeta212", &inWindow212LL);
    tree_full_3m->Branch("alphaBeta214", &inWindow214LL);
    tree_clean_3m->Branch("energy", &energy);
    tree_clean_3m->Branch("x", &posX);
    tree_clean_3m->Branch("y", &posY);
    tree_clean_3m->Branch("z", &posZ);
    tree_clean_3m->Branch("nhitsClean", &nhitsCleaned);
    tree_clean_3m->Branch("nhitsCorrected", &correctedNhits);
    tree_clean_3m->Branch("ITR", &ITRVal);
    tree_clean_3m->Branch("gtid", &gtid);
    tree_clean_3m->Branch("runid", &run_number);
    tree_clean_3m->Branch("alphaBeta212", &inWindow212LL);
    tree_clean_3m->Branch("alphaBeta214", &inWindow214LL);

    // deadtime calculation variables
    bool veto_flag          = false;
    bool lone_follower_flag = false;
    double loneFollowerTime = 0;
    double pileupTime       = 0;
    int num_vetos           = 0;
    ULong64_t veto_start_time, lone_start_time;

    // fuck me, that was a lot of branches. Now we do the loop over entries
    for (int ientry = 0; ientry < chain->GetEntries(); ientry++){
        chain->GetEntry(ientry);
        
        // flag to track if event has been tagged as 214, 212 or (alpha, N)
        bool tagged = false;

        // deadtime calculations
        // check if event is high energy or tagged as a muon by DC
        if (nhitsCleaned > 5000 or (DCflagged&0x80) != 0x80){
            std::cout << "\nWe have a high nhit or muon: GTID: " << gtid << std::endl;
            std::cout << "Muon flag result: " << (DCflagged&0x80!=0x80) << std::endl;
            std::cout << "NHits: " << nhitsCleaned << std::endl;
            // check if we are inside a lone muon follower veto window
            if (lone_follower_flag == true){
                std::cout << "Inside lone follower flag!" << std::endl;
                // add the dT and switch off the muon follower veto
                ULong64_t deltaT = ((clock50-lone_start_time) & 0x7FFFFFFFFFF)*20.0;
                loneFollowerTime += deltaT;
                std::cout << "Added " << deltaT << " ns to lone follower time." << std::endl;
                lone_follower_flag = false;
            }

            if (veto_flag == false){
                std::cout << "Setting high Nhits / muon veto to True." << std::endl;
                num_vetos++;
                veto_flag = true;
                veto_start_time = clock50;
                continue;
            }
            else{
                std::cout << "High Nhits / muon veto is True and nhits / muon veto triggered - we have pileup!" << std::endl;
                // we have pileup! Need to calculate the additional pileup time.
                ULong64_t deltaT = ((clock50-veto_start_time) & 0x7FFFFFFFFFF)*20.0;
                std::cout << "clock50 Ticks: " << clock50 << std::endl;
                std::cout << "Veto Start Ticks: " << veto_start_time << std::endl;
                std::cout << "Added " << deltaT * pow(10, -9) << " s to pileupTime.\n" << std::endl;
                pileupTime += deltaT;

                // reset the veto window
                veto_start_time = clock50;
                continue;
            }
        }

        // now handle the veto window for follower events from highE / muon veto
        if (veto_flag == true){
            
            ULong64_t deltaT = ((clock50-veto_start_time) & 0x7FFFFFFFFFF)*20.0;
            // std::cout << "high nhits / muon veto flag is true. deltaT since start of veto is: " << deltaT * pow(10, -9) << " s.";
            // check if event falls inside veto window
            if (deltaT < 2.0e10){
                // std::cout << "DeltaT since veto start is < 20 s. Skipping event." << std::endl;
                continue;
            
            }
            else{
                std::cout << "DeltaT since start of veto is > 20 s. Switching veto window flag to false." << std::endl;
              // we are no longer inside the veto time window --> switch it off!
              veto_flag = false;
            }
        }

        // check if we have a lone muon follower (i.e. muon at end of previous run)
        // this can only trigger if there was a muon at the end of the previous run
        if ((DCflagged&0x4000)!=0x4000){
            // std::cout << "We have a lone muon follower!" << std::endl;
            if (lone_follower_flag == false){
                std::cout << "Activating lone follower flag. " << std::endl;
                lone_follower_flag = true;
                lone_start_time = clock50;
                continue;
            } else{
            
                // we are within a lone follower veto period!
                ULong64_t deltaT = ((clock50 - lone_start_time) & 0x7FFFFFFFFFF)*20.0;
                // std::cout << "We are already within lone follower veto period. Added " << deltaT << " ns to lone follower time." << std::endl;
                loneFollowerTime += deltaT;
                lone_start_time = clock50;
                continue;
            }
        }
        else {
            // std::cout << "Event is not tagged as a muon follower." << std::endl;
            // check if event is NOT flagged as muon follower, but muon follower flag is on
            if (lone_follower_flag == true){
                std::cout << "lone follower flag was on. Switching off (no longer inside muon follower veto period)." << std::endl;
                // but we haven't triggered the other blocks! --> we've left the muon follower deadtime window
                lone_follower_flag = false;
            }
        }
        // std::cout << "Deadtime = " << pileupTime * pow(10, -9) << " + " << loneFollowerTime * pow(10, -9) << " s." << std::endl;
        // now we check if the event passes regular DC and reconstruction cuts
        if (((DCapplied & 0x2100000042C2) & DCflagged ) != (DCapplied & 0x2100000042C2)) {
                continue;
            }
        // check reconstruction validity
        if (fit == false or fitValid == false){
            continue;
        }

        // apply energy cuts
        if (energy < ENERGY_LOW || energy > ENERGY_HIGH){
            continue;
        }
        
        // work out the FV
        double z = posZ - AV_OFFSET;
        double rsq = posX*posX + posY*posY + z*z;
        
        // apply 6 m FV cut
        if (rsq > 6000*6000){
            continue;
        }

        // event is observed in the ROI!
        // check if it's a tagged coincidence event!
        bool bipo214 = false;
        bool bipo212 = false;
        bool alphaN  = false;
        if (std::count(gtid_214.begin(), gtid_214.end(), gtid)){
            bipo214 = true;
            tagged  = true;
        }

        // std::cout << "Coincidence Tagging\nBiPo214: " << bipo214 << "\nBiPo212: " << bipo212 << "\n(alpha, n): " << alphaN << std::endl;

        // work out which tree to fill for this event
        if (tagged == false){
            // fill full 6 m spectrum and clean 6 m spectrum
            tree_full_6m->Fill();
            tree_clean_6m->Fill();
        } else{
            // only fill full spectrum
            tree_full_6m->Fill();
        }

        if (rsq < 5500*5500){
            if (tagged == false){
                tree_full_5p5m->Fill();
                tree_clean_5p5m->Fill();
            } else{
                tree_full_5p5m->Fill();
            }
        }
        if (rsq < 5000*5000){
            if (tagged == false){
                tree_full_5m->Fill();
                tree_clean_5m->Fill();
            } else{
                tree_full_5m->Fill();
            }
        }
        if (rsq < 4500*4500){
            if (tagged == false){
                tree_full_4p5m->Fill();
                tree_clean_4p5m->Fill();
            } else{
                tree_full_4p5m->Fill();
            }
        }
        if (rsq < 4000*4000){
            if (tagged == false){
                tree_full_4m->Fill();
                tree_clean_4m->Fill();
            } else{
                tree_full_4m->Fill();
            }
        }
        if (rsq < 3500*3500){
            if (tagged == false){
                tree_full_3p5m->Fill();
                tree_clean_3p5m->Fill();
            } else{
                tree_full_3p5m->Fill();
            }
        }
        if (rsq < 3000*3000){
            if (tagged == false){
                tree_full_3m->Fill();
                tree_clean_3m->Fill();
            } else{
                tree_full_3m->Fill();
            }
        }
    }

    // completed loop through entries --> calculate livetime and write the output files
    g->cd();
    tree_full_6m->Write();
    tree_full_5p5m->Write();
    tree_full_5m->Write();
    tree_full_4p5m->Write();
    tree_full_4m->Write();
    tree_full_3p5m->Write();
    tree_full_3m->Write();
    tree_clean_6m->Write();
    tree_clean_5p5m->Write();
    tree_clean_5m->Write();
    tree_clean_4p5m->Write();
    tree_clean_4m->Write();
    tree_clean_3p5m->Write();
    tree_clean_3m->Write();
    g->Close();
    
    // write out deadtime to .txt file
    double deadtime = (20 * num_vetos) + (pileupTime + loneFollowerTime) * pow(10, -9); // in s
    std::cout << "Deadtime: " << deadtime << std::endl;
    std::cout << "pileUpTime: " << pileupTime << std::endl;
    std::cout << "LoneFollowerTime: " << loneFollowerTime << std::endl;
    double livetime = (duration * pow(10, -9)) - deadtime; // in s
    std::cout << "Livetime: " << livetime << std::endl;
    std::ofstream flivetime;
    flivetime.open("/data/snoplus3/hunt-stokes/clean_multisite/data_spectrum_runs_general_coincidence/livetime/" + std::to_string(run_number) + ".txt");
    flivetime << std::to_string(livetime);
    flivetime.close();
}

int main(int argc, char* argv[]){

    // input to the script is a run number
    std::string runNum      = argv[1];

    // find all the data ntuple subruns corresponding to this run number
    std::string inFile1;
    std::string data_path = "/data/snoplus3/SNOplusData/processing/fullFill/rat-7.0.8/ntuples";

    // be mindful of the reprocessing of earlier data
    if (std::stoi(runNum) < 307613){
        inFile1 = data_path + "/" + "Analysis20R_r0000";
    }
    if (std::stoi(runNum) >= 307613){
        inFile1 = data_path + "/" + "Analysis20_r0000";
    }
    std::string inFile2 = "*";
    const std::string inFileName = inFile1+runNum+inFile2;
    std::cout << inFileName << std::endl;

    // create a vector of all the subrun files
    std::vector<std::string> filelist = globVector(inFileName);

    // find the duration of the run
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

        // create a TChain that contains the entries of all the events in this run
        TChain *chain = new TChain("output");
        // add each subrun file to the TChain
        for (int iFile = 0; iFile < filelist.size(); iFile++){
            std::string fname = filelist.at(iFile);
            chain->AddFile(fname.c_str());
        }
        std::cout << "TChain contains " << chain->GetEntries() << " entries." << std::endl;
        
        // create a vector to save the GTIDs of each coincidence event
        std::vector<int> gtid_coincidence;

        // now run the coincidence tagging
        coincidence_tagger(chain, zOff, runID, gtid_coincidence);

        // now extract the energy spectrum before and after coincidence tagging
        extract_energy_spectrum(chain, runID, zOff, gtid_coincidence, duration);
        std::cout << "Extracted data spectra." << std::endl; 



    }
}