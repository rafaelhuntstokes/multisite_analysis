#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
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

void count_simulated(std::string input_file, double zOff, std::vector<int> &entire_spectrum, std::vector<int> &roi_spectrum, std::vector<int> &systematic_fvs){
    
    TFile *f = TFile::Open(input_file.c_str());
    TNtuple *T = (TNtuple *) f->Get("output");
    
    double posX, posY, posZ, mcR, mcX, mcY, mcZ, energy, itr, posFOM;
    UInt_t posFOM_hits;
    int Nhits, mcIndex, ev;
    int pdg;
    bool fit, fitValid;

    ULong64_t clock50;

    T->SetBranchAddress("posx", &posX);
    T->SetBranchAddress("posy", &posY);
    T->SetBranchAddress("posz", &posZ);
    T->SetBranchAddress("nhitsCleaned", &Nhits);
    T->SetBranchAddress("energy", &energy);
    T->SetBranchAddress("scintFit", &fit);
    T->SetBranchAddress("fitValid", &fitValid);
    T->SetBranchAddress("mcIndex", &mcIndex);
    T->SetBranchAddress("evIndex", &ev);
    T->SetBranchAddress("itr", &itr);
    T->SetBranchAddress("posFOM", &posFOM);
    T->SetBranchAddress("posFOM2", &posFOM_hits);
    // T->SetBranchAddress("mcPosr", &mcR);
    T->SetBranchAddress("mcPosx", &mcX);
    T->SetBranchAddress("mcPosy", &mcY);
    T->SetBranchAddress("mcPosz", &mcZ);
    T->SetBranchAddress("parentpdg1", &pdg);
    T->SetBranchAddress("clockCount50", &clock50);
    
    int full_numSimulated_6m   = 0;

    int roi_numSimulated_6m   = 0;
    int roi_numSimulated_5p5m = 0;
    int roi_numSimulated_5m   = 0;
    int roi_numSimulated_4p5m = 0;
    int roi_numSimulated_4m   = 0;
    int roi_numSimulated_3p5m = 0;
    int roi_numSimulated_3m   = 0;

    int systematic_fv_central = 0;
    int systematic_fv_plus    = 0;
    int systematic_fv_neg     = 0;
    double MAJOR_R = 3500.0;
    double TUBE_R = 1500.0;
    
    /*
    This script counts events that: 
    
    1. Trigger the detector
    2. Are prompt candidates (e.g. bismuth214 / bismuth212)
    3. Reconstruct into the region of interest
    
    This number is used to calculate the TAGGING EFFICIENCY. This number is
    not valid for working out the efficiency scaled rate of BiPo214/212. 
    You cannot make estimates of background contaminations from this.
    */
    bool in_roi;
    bool re_trig;
    ULong64_t previous_clock;
    for (int iEntry = 0; iEntry < T->GetEntries(); iEntry++){
        T->GetEntry(iEntry);

        // include all events that trigger the detector 
        // - negative numbers are simulated events that do not trigger
        // 0 are triggered events
        // ev > 0 are retriggers
        in_roi = true;
        re_trig = false;
        if (ev > 0){
            re_trig = true;
            // skip retriggers
            // continue;
        }

        // check that the most energetic parent of this event is the Po214 (before quenching Q-val = 7.8 MeV)
        // may also have done Bi214 --> Tl210 in which case parent is Bi214 
        // if (pdg != 1000842140){
        // if (pdg == 1000832120){
            // std::cout << "Most energetic parent of this event not Po212!\nIt is: " << pdg << std::endl;
            // continue;
        // }

        if (iEntry != 0){
            // check if the event would be removed by the retrigger cut
            if (((clock50 - previous_clock) & 0x7FFFFFFFFFF) * 20 < 500){
                // event removed by retrigger cut
                in_roi  = false;
            }
        }
        // update previous clock for next event
        previous_clock = clock50;
        
        // check reconstruction is valid --> if not the event would be rejected anyway
        if (fit == false){
            in_roi = false;
        }
        if (fitValid == false){
            in_roi = false;
        }

        double z        = mcZ - zOff;
        double reconZ   = posZ - zOff;
        double Rsq_mc   = mcX*mcX + mcY*mcY + z*z;
        double rReconSq = posX*posX + posY*posY + reconZ*reconZ;
        double Rsq      = rReconSq;

        double posFOM_normed = posFOM / posFOM_hits;
        // check the prompt event passes the energy ROI cuts
        if (energy < 3.5 || energy > 5.0){
            in_roi = false;
        }
        if (itr < 0.22 or itr > 0.3){
            in_roi = false;
        }
        if (posFOM_normed < 13.7){
            in_roi = false;
        }

        // check simulated inside the AV and not a retrigger
        if (Rsq_mc < 6000*6000 and re_trig == false){
            full_numSimulated_6m += 1;
        }

        // check if event reconstructs inside the AV + ROI and triggers the detector
        if (in_roi == true and rReconSq < 6000*6000){
            roi_numSimulated_6m += 1;
        
            // check if in each volume
            if (Rsq <= 5500*5500){
                if (in_roi == true and rReconSq < 5500*5500){
                    roi_numSimulated_5p5m += 1;
                }
                if (Rsq <= 5000*5000){
                    if (in_roi == true and rReconSq < 5000*5000){
                        roi_numSimulated_5m += 1;
                    }       
                    if (Rsq <= 4500*4500){
                        if (in_roi == true and rReconSq < 4500*4500){
                            roi_numSimulated_4p5m += 1;
                        }
                        if (Rsq <= 4000*4000){
                            if (in_roi == true and rReconSq < 4000*4000){
                                roi_numSimulated_4m += 1;
                            }
                            if (Rsq <= 3500*3500){
                                if (in_roi == true and rReconSq < 3500*3500){
                                    roi_numSimulated_3p5m += 1;
                                }
                                if (Rsq <= 3000*3000){
                                    if (in_roi == true and rReconSq < 3000*3000){
                                        roi_numSimulated_3m += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    
        if (in_roi == true and rReconSq < 4000*4000){
            systematic_fv_central += 1;
        }
        if (in_roi == true and rReconSq < (4000 - 40)*(4000 - 40)){
            systematic_fv_neg += 1;
        }
        if (in_roi == true and rReconSq < (4000 + 40)*(4000 + 40)){
            systematic_fv_plus += 1;
        }
    }
    entire_spectrum.at(0) = entire_spectrum.at(0) + full_numSimulated_6m;

    roi_spectrum.at(0) = roi_spectrum.at(0) + roi_numSimulated_3m;
    roi_spectrum.at(1) = roi_spectrum.at(1) + roi_numSimulated_3p5m;
    roi_spectrum.at(2) = roi_spectrum.at(2) + roi_numSimulated_4m;
    roi_spectrum.at(3) = roi_spectrum.at(3) + roi_numSimulated_4p5m;
    roi_spectrum.at(4) = roi_spectrum.at(4) + roi_numSimulated_5m;
    roi_spectrum.at(5) = roi_spectrum.at(5) + roi_numSimulated_5p5m;
    roi_spectrum.at(6) = roi_spectrum.at(6) + roi_numSimulated_6m;

    systematic_fvs.at(0) = systematic_fvs.at(0) + systematic_fv_neg;
    systematic_fvs.at(1) = systematic_fvs.at(1) + systematic_fv_central;
    systematic_fvs.at(2) = systematic_fvs.at(2) + systematic_fv_plus;   
}

inline bool exists_test0 (const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

int main(int argc, char* argv[]){

    std::string run_list_fname = argv[1];
    std::vector<int> runs;

    // load the run list and push back each entry to the runs vector to loop over
    std::fstream runList;
    runList.open(run_list_fname, std::ios::in);
    if (runList.is_open()){
        std::string run;
        while (std::getline(runList, run)){
            // std::cout << "Reading in run: " << run << std::endl;
            runs.push_back(std::stoi(run));
        }
    }
    std::cout << "Created the file list..." << std::endl;
    std::vector<int> entire_spectrum;
    entire_spectrum.push_back(0);

    std::vector<int> roi_spectrum;
    roi_spectrum.push_back(0);
    roi_spectrum.push_back(0);
    roi_spectrum.push_back(0);
    roi_spectrum.push_back(0);
    roi_spectrum.push_back(0);
    roi_spectrum.push_back(0);
    roi_spectrum.push_back(0);

    std::vector<int> systematic_fvs;
    systematic_fvs.push_back(0);
    systematic_fvs.push_back(0);
    systematic_fvs.push_back(0);
    
    int number_files = runs.size();
    for (int i = 0; i < runs.size(); i++){
        std::cout << "Run ID: " << i << std::endl;
        // std::string input_file = "/data/snoplus3/SNOplusData/production/miniProd_RAT-7-0-14_ASCI_RATHS_newRecoordination/ntuple/simulationgammas_2p2MeV_" + std::to_string(runs.at(i)) + ".ntuple.root";
        // std::string input_file = "/data/snoplus3/SNOplusData/production/miniProd_RAT-7-0-14_ASCI_RATHS_newRecoordination/ntuple/simulationAlphaN_LAB_13C_" + std::to_string(runs.at(i)) + ".ntuple.root";
        // std::string input_file = "/data/snoplus3/SNOplusData/production/miniProd_RAT-7-0-14_ASCI_RATHS_newRecoordination/ntuple/simulationBiPo214_" + std::to_string(runs.at(i)) + ".ntuple.root";
        // std::string input_file = "/data/snoplus3/SNOplusData/production/miniProd_RAT-7-0-14_ASCI_RATHS_newRecoordination/ntuple/simulationB8_solar_numu_" + std::to_string(runs.at(i)) + ".ntuple.root";
        std::string input_file = "/data/snoplus3/SNOplusData/production/miniProd_RAT-7-0-14_ASCI_RATHS_newRecoordination/ntuple/simulationTl208_" + std::to_string(runs.at(i)) + ".ntuple.root";
        bool file_exists = exists_test0(input_file);
        if (file_exists == false){
            std::cout << "File doesn't exist." << std::endl;
            number_files--;
            continue;
        }                                                                                                                                                           
        RAT::DB *db = RAT::DB::Get();
        db->LoadDefaults();
        RAT::DS::Run run;
        run.SetRunID(runs.at(i));
        db->BeginOfRun(run);
        std::vector<double> AVPos = RAT::GeoUtil::UpdateAVOffsetVectorFromDB();
        double zOff = AVPos[2];

        std::cout << "The AV offset for run " << runs.at(i) << " is " << zOff << " mm." << std::endl;

        // find the number of events simulated in each volume for this run
        
        count_simulated(input_file, zOff, entire_spectrum, roi_spectrum, systematic_fvs);

        // add the results to totals
        std::cout << "\nFull Spectrum Number: " << std::endl;

        for (int i =0; i < entire_spectrum.size(); i++){

            std::cout << entire_spectrum.at(i) << std::endl;
        }

        std::cout << "\nROI Spectrum Number: " << std::endl;
        for (int i =0; i < roi_spectrum.size(); i++){

            std::cout << roi_spectrum.at(i) << std::endl;
        }

        std::cout << "Systematic Study FVs: " << std::endl;
        for (int i = 0; i < systematic_fvs.size(); i++){
            std::cout << systematic_fvs.at(i) << std::endl;
        }

        
    }
    
    std::cout << "Finished Counting!" << std::endl; 
    std::cout << "Number of runs used: " << number_files << std::endl;   
    // std::string input_file = "/data/snoplus3/griddata/Prod_7_0_8_BiPo212/prod_used_for_tuning/prod2.root";
    
    
    return 0;   
}
