#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DB.hh>
#include <RAT/GeoUtils.hh>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector3.h>
#include <string>
#include <RAT/SunUtil.hh>
#include <glob.h>
#include <fstream>
#include "TNtuple.h"
#include <RAT/DU/DSReader.hh>
#include <TChain.h>
/*

This script is used to identify solar candidate events between 5 --> 12 MeV.

The identification is done using ntuple files for a given run number.

Deadtime, cleanliness, FV, energy, ITR etc. cuts are applied. Events which
pass every cut is considered a 'good' solar candidate and the run ID and GTID
is saved to an output TTRee.

In addition, the TTRee contains event level information (excepting, of course, time 
information) for preliminary checks on the selection quality. Energy, X, Y, Z and ITR 
are saved.

The deadtime corrected run time is returned and saved to a .txt file.

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

void identify_solar_candidates(int run_number, float fv_cut, float z_cut, std::string input_fpath, std::string output_fpath, double offset, double &deadtime){
    /*
    Function loads in a run_by_run test .root ntuple file and loops over all the events. 
    
    Every event that passes the FV, ITR, Z, Recon and energy cuts is saved to TTree.

    A run_by_run ntuple file is returned containing the [energy, position, ITR, GTID, RUNID] of each event and for each energy.
    */

    // define ntuple input file location and output TTree file locations
    std::string input_fname;
    std::cout << run_number << std::endl;
    if (run_number < 307613){
    std::cout << "Setting fname..." << std::endl;
        input_fname  = input_fpath + "/Analysis20R_r0000" + std::to_string(run_number) + "*.ntuple.root";
    } else {
        input_fname  = input_fpath + "/Analysis20_r0000" + std::to_string(run_number) + "*.ntuple.root";
    }
    std::cout << "input fname is: " << input_fname << std::endl;
    
    std::string output_fname = output_fpath + "/" + std::to_string(run_number) + ".root";

    // create the output ntuple with an output TTree containing tagged solar events
    TFile *output_file    = new TFile(output_fname.c_str(), "RECREATE");
    TTree *tagged_solar   = new TTree("tagged_solar", "tagged_solar");

    // define the variables to fill branches with
    double energy, x, y, z, itr;
    int nhits_cleaned, gtid;
    bool fit, fitvalid;
    
    tagged_solar->Branch("energy", &energy);
    tagged_solar->Branch("x", &x);
    tagged_solar->Branch("y", &y);
    tagged_solar->Branch("z", &z);
    tagged_solar->Branch("itr", &itr);
    tagged_solar->Branch("gtid", &gtid);
    tagged_solar->Branch("run_number", &run_number);
    tagged_solar->Branch("nhits_cleaned", &nhits_cleaned);

    // setup deadtime veto code
    bool veto_flag          = false;
    bool lone_follower_flag = false;
    double loneFollowerTime = 0;
    double pileupTime       = 0;
    int num_vetos           = 0;
    ULong64_t veto_start_time, lone_start_time, dcApplied, dcFlagged, clock50;
    int nhitsCleaned;
    Int_t fPass;

    // use glob function to create a filelist of each subrun ratds corresponding to this run
    std::vector<std::string> filelist = globVector(input_fname);
    
    // now we can TChain all the subrun files into a single chain to iterate over
    TChain *chain = new TChain("output");
    for (int ifile = 0; ifile < filelist.size(); ifile++){
        std::string fname = filelist.at(ifile);
        std::cout << "Adding file: " << fname << std::endl;
        chain->AddFile(fname.c_str());
    }
    
    // set the branch addresses of the TChain events
    chain->SetBranchAddress("posx", &x);
    chain->SetBranchAddress("posy", &y);
    chain->SetBranchAddress("posz", &z);
    chain->SetBranchAddress("clockCount50", &clock50);
    chain->SetBranchAddress("nhitsCleaned", &nhits_cleaned);
    chain->SetBranchAddress("scintFit", &fit);
    chain->SetBranchAddress("fitValid", &fitvalid);
    chain->SetBranchAddress("energy", &energy);
    chain->SetBranchAddress("dcApplied", &dcApplied);
    chain->SetBranchAddress("dcFlagged", &dcFlagged);
    chain->SetBranchAddress("itr", &itr);
    chain->SetBranchAddress("eventID", &gtid);
    
    std::cout << "Created a TChain for run " << run_number << " containing " << chain->GetEntries() << "." << std::endl;
    
    // begin loop over every entry
    for (int ientry = 0; ientry < chain->GetEntries(); ientry++){
        
        // all the variables defined above are now filled with this entry's values
        chain->GetEntry(ientry);

        // apply neck hotspot / weird stuff removal and deadtime cuts
        // check if high nhits or tagged as a muon by DC
        if (nhits_cleaned > 5000 or (dcFlagged&0x80)!=0x80){
            bool muon_veto = (dcFlagged&0x80)!=0x80;
            std::cout << "Nhits veto: " << nhits_cleaned << std::endl;
            std::cout << "Muon veto: " << muon_veto << std::endl;
            if (lone_follower_flag == true){
                ULong64_t deltaT = ((clock50-lone_start_time) & 0x7FFFFFFFFFF)*20.0;
                
                loneFollowerTime += deltaT;
                std::cout << "Adding " << deltaT * pow(10, -9) << " s to loneFollowerTime." << std::endl; 
                lone_follower_flag = false;
            }

            if (veto_flag == false){
                num_vetos++;
                veto_flag = true;
                veto_start_time = clock50;
                continue;
            }
            else{
                // we have pileup! Need to calculate the additional pileup time.
                ULong64_t deltaT = ((clock50-veto_start_time) & 0x7FFFFFFFFFF)*20.0;
                pileupTime += deltaT;
                std::cout << "Adding " << deltaT * pow(10, -9) << " s to pileupTime." << std::endl;
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
        if ((dcFlagged&0x4000)!=0x4000){
            if (lone_follower_flag == false){
            lone_follower_flag = true;
            lone_start_time = clock50;
            continue;
            } else{
            // we are within a lone follower veto period!
            ULong64_t deltaT = ((clock50 - lone_start_time) & 0x7FFFFFFFFFF)*20.0;
            loneFollowerTime += deltaT;
            std::cout << "Adding " << deltaT * pow(10, -9) << " s to loneFollowerTime." << std::endl;
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

        // check if basic DC cuts are passed
        if (((dcApplied & 0x2100000042C2) & dcFlagged ) != (dcApplied & 0x2100000042C2)) {
            continue;
        }

        // check the reconstruction succeeded
        // check reconstruction validity
        if (fit == false or fitvalid == false){
            continue;
        }
        
        
        // check if we pass the ITR cut
        if (itr < 0.18 or itr > 0.3){
            continue;
        }


        // convert to AV coordinates and apply selection cuts
        z = z - offset;
        double radius_sq = x*x + y*y + z*z;
        // std::cout << "R: " << event_position_recon.Mag() << std::endl;
        if (radius_sq > fv_cut*fv_cut){
            // failed FV cut!
            continue;
        }
        // std::cout << "Passed FV cut." << std::endl;
        if (z < z_cut){
            // not in the scintillator!
            continue;
        }

        // check that event energy fits into the large domain
        if (energy < 5.0 or energy > 12.0){
            // outside energy ROI!
            continue;
        }

        // passed all event selection cuts! Save the output to the TTree.
        std::cout << "Found a solar candidate!" << std::endl;
        tagged_solar->Fill();
    }

    std::cout << "Writing results to file." << std::endl;
    // finished loop over all events and TTrees are filled --> write them to the output file
    tagged_solar->Write();
    output_file->Close();

    // calculate the deadtime
    std::cout << "num vetos: " << num_vetos << std::endl;
    std::cout << "Pileup Time: " << pileupTime * pow(10, -9) << " s." << std::endl;
    std::cout << "Lone Follower Time: " << loneFollowerTime * pow(10, -9) << " s." << std::endl;
    deadtime = (20 * num_vetos) + (pileupTime + loneFollowerTime) * pow(10, -9); // in seconds
}

int main(int argc, char* argv[]){

    int run_number           = std::stoi(argv[1]);
    double fv_cut            = std::stod(argv[2]);
    double z_cut             = std::stod(argv[3]);
    std::string input_fpath  = argv[4];
    std::string output_fpath = argv[5];

    // find raw runtime of this run in seconds
    RAT::DB *db = RAT::DB::Get();
    db->LoadDefaults();
    RAT::DS::Run run;
    run.SetRunID(run_number);
    db->BeginOfRun(run);
    std::cout << "Run ID = " << run_number << std::endl;
    double start_day, start_sec, end_day, end_sec, deadtime;
    long long start_time, end_time, livetime;
    RAT::DBLinkPtr dblink = db->GetLink("RUN");
    start_day = dblink->GetD("start_day");
    start_sec = dblink->GetD("start_sec");
    end_day = dblink->GetD("stop_day");
    end_sec = dblink->GetD("stop_sec");
    start_time = start_day * 24 * 3600 + start_sec;
    end_time = end_day * 24 * 3600 + end_sec;
    
    std::vector<double> AVPos = RAT::GeoUtil::UpdateAVOffsetVectorFromDB();
    double offset = AVPos[2];
    std::cout << "AV Shift is: " << offset << " mm." << std::endl;
    
    // run the analysis!
    identify_solar_candidates(run_number, fv_cut, z_cut, input_fpath, output_fpath, offset, deadtime);
    
    // calculate deadtime adjusted livetime
    livetime = end_time - start_time - deadtime;
    std::cout << "Deadtime adjusted livetime is: " << livetime << " s." << std::endl;

    // read out this livetime into a file
    std::ofstream flivetime;
    std::string fpath = output_fpath + "/livetimes/" + std::to_string(run_number) + ".txt";
    flivetime.open(fpath);
    flivetime << std::to_string(livetime);
    flivetime.close();
}
