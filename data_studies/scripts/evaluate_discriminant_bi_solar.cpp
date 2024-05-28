#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/PMTInfo.hh>
#include <RAT/DU/LightPathCalculator.hh>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
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

double multisite_discriminant(std::vector<double> time_residuals, TH1D* pdf_B8_multi, TH1D* pdf_Tl208_multi){

    /*
    Loop through each time residual and find the loglikelihood ratio for B8 and Tl208.
    */

    double loglikelihood  = 0.0;
    int contributing_hits = 0;  // normalise loglikelihood by number of contributing hits 
    for (int ires = 0; ires < time_residuals.size(); ires++){
        
        double residual = time_residuals.at(ires);

        // find probability of being signal and background
        double prob_signal     = pdf_B8_multi->GetBinContent(pdf_B8_multi->FindBin(residual));
        double prob_background = pdf_Tl208_multi->GetBinContent(pdf_Tl208_multi->FindBin(residual));

        // check if we will have divide by zero error or log(0) error
        if (prob_background <= 0 or prob_signal <= 0){
            continue;
        }
        contributing_hits++;

        // likelihood ratio
        double ratio = prob_signal / prob_background;
        ratio = log(ratio);

        loglikelihood += ratio;
    }
    
    // normalise by number of hits
    loglikelihood /= contributing_hits;

    return loglikelihood;
}

void evaluate_discriminants(int run_number, float fv_cut, float z_cut, int event_type){
    /*
    Multisite discriminant, ITR and GTIDs are saved for each extracted Bi214 (event_type = 1)
    or solar candidate (event_type = 2).

    */

    std::string event_string;
    std::string tree_name;
    float energy_high, energy_low;
    if (event_type == 1){
        event_string = "bi214";
        tree_name = "bi214";
        energy_low = 1.25;
        energy_high = 3.0;
    }
    if (event_type == 2){
        event_string = "above_5MeV";
        tree_name = "above_5MeV";
        energy_low = 5.0;
        energy_high = 12.0;
    }
    if (event_type != 2 and event_type != 1){
        std::cout << "Event type " << event_type << " not recognised!." << std::endl;
        return;
    }

    // define input and output paths
    // std::string input_fname = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/" + event_string + "/reprocessed_7.0.15_ratds/" + std::to_string(run_number) + "*.root";    
    // std::string output_fname = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/" + event_string + "/reprocessed_7.0.15_discriminants_2.5_3.0_PDF/" + std::to_string(run_number) + ".root";
    
    std::string input_fname = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/" + event_string + "/bismuth214_extracted_ratds/" + std::to_string(run_number) + "*.root"; 
    std::string output_fname = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/" + event_string + "/discriminants_7.0.8_2.5_3.125_PDF/" + std::to_string(run_number) + ".root";
    
    // load the PDF files for each isotope
    std::string working_directory = "/data/snoplus3/hunt-stokes/multisite_clean/mc_studies";
    std::string pdf_B8_fname      = working_directory + "/run_by_run_pdf/B8_solar_nue/total.root";
    std::string pdf_Tl208_fname   = working_directory + "/run_by_run_pdf/Tl208/total.root";
    TFile *PDF_B8 = new TFile(pdf_B8_fname.c_str());
    TFile *PDF_Tl208 = new TFile(pdf_Tl208_fname.c_str());

    // using the same PDFs for both isotopes, even if the energy range doesn't exactly match the Bi214 / solar energy, should be okay
    TH1D *pdf_B8_multi = dynamic_cast<TH1D*>(PDF_B8->Get("multi_2.5_3.125"));
    TH1D *pdf_Tl208_multi = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_2.5_3.125"));

    // create the output ntuple with an output TTree containing discriminants for each energy PDF
    TFile *output_file    = new TFile(output_fname.c_str(), "RECREATE");

    TTree *info = new TTree(tree_name.c_str(), tree_name.c_str());
    
    // define the variables to fill branches with
    double energy, x, y, z, dlogL, itr;
    Int_t nhitsCleaned, gtid;

    info->Branch("nhitsCleaned", &nhitsCleaned);
    info->Branch("energy", &energy);
    info->Branch("itr", &itr);
    info->Branch("x", &x);
    info->Branch("y", &y);
    info->Branch("z", &z);
    info->Branch("dlogL", &dlogL);
    info->Branch("gtid", &gtid);

    // setup deadtime veto code
    bool veto_flag          = false;
    bool lone_follower_flag = false;
    double loneFollowerTime = 0;
    double pileupTime       = 0;
    int num_vetos           = 0;
    ULong64_t veto_start_time, lone_start_time, dcApplied, dcFlagged, clock50;
    Int_t fPass;

    // use glob function to create a filelist of each subrun ratds corresponding to this run
    std::vector<std::string> filelist = globVector(input_fname);
    for (int ifile = 0; ifile < filelist.size(); ifile++){

        /*
        From here we load the RATDS file for this run and evaluate the event by event discriminants!
        */

        std::cout << "Loading File: " << filelist.at(ifile) << std::endl;
        RAT::DU::DSReader dsReader(filelist.at(ifile).c_str());
        std::cout << "Loaded File." << std::endl;

        int num_entries = dsReader.GetEntryCount();
        std::cout << "Found " << num_entries << " events." << std::endl;

        // setup time residual calculator, Point3D and get PMT info
        RAT::DU::TimeResidualCalculator timeResCalc = RAT::DU::Utility::Get()->GetTimeResidualCalculator();
        size_t fPSUPSystemId = RAT::DU::Point3D::GetSystemId("innerPMT");
        size_t av_system_id  = RAT::DU::Point3D::GetSystemId("av");
        const RAT::DU::PMTInfo& pmt_info = RAT::DU::Utility::Get()->GetPMTInfo();

        
        // begin loop over every entry
        for (int ientry = 0; ientry < num_entries; ientry++){

            const RAT::DS::Entry& rDS = dsReader.GetEntry(ientry);

            // check event triggered detector
            if (rDS.GetEVCount() == 0){
                // did not trigger!
                continue;
            }

            // get primary trigger
            const RAT::DS::EV& rEV = rDS.GetEV(0);

            // apply neck hotspot / weird stuff removal and deadtime cuts
            nhitsCleaned = rEV.GetNhitsCleaned();
            fPass        = rEV.GetDataCleaningFlags().GetLatestPass();
            clock50      = rEV.GetClockCount50();
            gtid         = rEV.GetGTID();
            if(fPass >= 0){
                dcApplied = rEV.GetDataCleaningFlags().GetApplied(fPass).GetULong64_t(0);
                dcFlagged = rEV.GetDataCleaningFlags().GetFlags(fPass).GetULong64_t(0);
            }

            // check if high nhits or tagged as a muon by DC
            if (nhitsCleaned > 5000 or (dcFlagged&0x80)!=0x80){

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
            RAT::DS::FitVertex fVertex = rEV.GetDefaultFitVertex();
        
            // check that a position/time fit exists                 
            if( !fVertex.ContainsPosition() ) continue;
            if( !fVertex.ValidPosition() ) continue;
            if( !fVertex.ContainsTime() ) continue;
            if( !fVertex.ValidTime() ) continue;
            // std::cout << "Passed deadtime, DC and reconstruction checks." << std::endl;
            // check if we pass the ITR cut
            itr = rEV.GetClassifierResult("ITR:scintFitter").GetClassification("ITR");
            // if (itr < 0.18 or itr > 0.3){
            //     continue;
            // }
            // std:: cout << "Passed ITR cut." << std::endl;
            // get the reconstructed variables
            RAT::DU::Point3D event_position_recon(fPSUPSystemId);
            event_position_recon.SetXYZ(fPSUPSystemId, fVertex.GetPosition());
            energy            = fVertex.GetEnergy();
            double recon_time = fVertex.GetTime();

            // convert to AV coordinates and apply selection cuts
            event_position_recon.SetCoordinateSystem(av_system_id);
            x = event_position_recon.X();
            y = event_position_recon.Y();
            z = event_position_recon.Z();
            // std::cout << "R: " << event_position_recon.Mag() << std::endl;
            // if (event_position_recon.Mag() > fv_cut){
                // failed FV cut!
                // continue;
            // }
            // std::cout << "Passed FV cut." << std::endl;
            // if (z < z_cut){
                // not in the scintillator!
                // continue;
            // }
            // std::cout << "Passed position cuts." << std::endl;
            // convert back to PSUP Coords for time residual calculation
            event_position_recon.SetCoordinateSystem(fPSUPSystemId);

            // check that event energy fits into the large domain
            // if (energy < energy_low or energy > energy_high){
                // outside energy ROI!
                // continue;
            // }
            // std::cout << "Passed energy cut." << std::endl;
            // passed all event selection cuts! Now we evaluate the discriminants;

            // each discriminant involves the time residuals, so evaluate and store (temporarily) the time residuals and PMT positions
            std::vector<double> time_residuals;
            std::vector<double> pmt_x;
            std::vector<double> pmt_y;
            std::vector<double> pmt_z;
            double tres_recon; 

            const RAT::DS::CalPMTs& calibratedPMTs = rEV.GetCalPMTs();
            for (size_t ipmt = 0; ipmt < calibratedPMTs.GetCount(); ipmt++){
                const RAT::DS::PMTCal& pmtCal             = calibratedPMTs.GetPMT(ipmt);
                RAT::DU::Point3D pmt_point(fPSUPSystemId, pmt_info.GetPosition(pmtCal.GetID()));
                const RAT::DU::PMTCalStatus& PMTCalStatus = RAT::DU::Utility::Get()->GetPMTCalStatus();
                // check the hit is 'good'
                if (PMTCalStatus.GetHitStatus(pmtCal) != 0){
                    continue;
                }

                tres_recon = timeResCalc.CalcTimeResidual(pmtCal, event_position_recon, recon_time, true);
                
                time_residuals.push_back(tres_recon);
                pmt_x.push_back(pmt_point[0]);
                pmt_y.push_back(pmt_point[1]);
                pmt_z.push_back(pmt_point[2]);
            }
        
            // now have all the time residuals and PMT positions saved in vector --> can fit the direction and find multisite dLog(l)
            
            std::cout << "Evaluating multisite discriminant." << std::endl;
            
            dlogL         = multisite_discriminant(time_residuals, pdf_B8_multi, pdf_Tl208_multi);

            // fill the ntuple
            info->Fill();
        }
    }

    std::cout << "Writing results to file." << std::endl;
    
    // finished loop over all events and TTrees are filled --> write them to the output file
    output_file->cd();
    info->Write();
    output_file->Close();
}

int main(int argc, char* argv[]){

    int run_number           = std::stoi(argv[1]);
    double fv_cut            = std::stod(argv[2]);
    double z_cut             = std::stod(argv[3]);
    int event_type           = std::stoi(argv[4]);

    // find raw runtime of this run in seconds
    RAT::DB *db = RAT::DB::Get();
    db->LoadDefaults();
    RAT::DS::Run run;
    run.SetRunID(run_number);
    db->BeginOfRun(run);
    std::cout << "Run ID = " << run_number << std::endl;
    
    // run the analysis!
    evaluate_discriminants(run_number, fv_cut, z_cut, event_type);
}