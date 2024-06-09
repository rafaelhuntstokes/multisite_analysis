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

int record_retriggers(int run_number, std::vector<int> &flagged_ignore_evs){
    /*
    Function checks if event extracted into ROI for this run is a retrigger.
    */

    std::string inFile1;
    std::string data_path = "/data/snoplus3/SNOplusData/processing/fullFill/rat-7.0.8/ntuples";

    // be mindful of the reprocessing of earlier data
    if (run_number < 307613){
        inFile1 = data_path + "/" + "Analysis20R_r0000";
    }
    if (run_number >= 307613){
        inFile1 = data_path + "/" + "Analysis20_r0000";
    }
    std::string inFile2 = "*";
    const std::string inFileName = inFile1+std::to_string(run_number)+inFile2;
    std::cout << "\nChecking " << inFileName << " for retriggers." << std::endl;

    // create a vector of all the subrun files
    std::vector<std::string> filelist = globVector(inFileName);

    if (filelist.size() == 0){
        std::cout << "No valid ntuple files for this run!" << std::endl;
        return 1;
    }

    // obtain the gtids corresponding to this run
    std::string gtid_file = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis/gtids/" + std::to_string(run_number) + ".txt";
    std::ifstream infile(gtid_file);

    // read gtids to a vector
    std::vector<int> gtids;
    std::string line;

    while(std::getline(infile, line)){
        std::stringstream ss(line);
        int number;

        while (ss >> number){
            gtids.push_back(number);
        }
    }
    infile.close();

    if (gtids.size() == 0){
        std::cout << "No GTIDs for this run!" << std::endl;
        return 1;
    }
    std::cout << "\nFound " << gtids.size() << " gtids in file." << std::endl;

    // create a TChain out of all the NTUPLES to loop through
    TChain *chain = new TChain("output");
    // add each subrun file to the TChain
    for (int iFile = 0; iFile < filelist.size(); iFile++){
        std::string fname = filelist.at(iFile);
        chain->AddFile(fname.c_str());
    }
    std::cout << "\nCreated TChain with " << chain->GetEntries() << " entries." << std::endl;

    // define TChain branch variables
    int gtid;
    Int_t triggerWord;
    ULong64_t clock50, previousClock;
    chain->SetBranchAddress("clockCount50", &clock50);
    chain->SetBranchAddress("eventID", &gtid);

    // identify and remove events that are retriggers or orphans
    for (int ientry = 0; ientry < chain->GetEntries(); ientry++){
        chain->GetEntry(ientry);

        // first event can't be a retrigger (unless it's a retrigger from event at end of last run ... but not worried about that)
        if (ientry == 0){
            // set the previous event time to this event
            previousClock = clock50;
            continue;
        }

        // check if this event is within 500 ns of previous event
        if (((clock50 - previousClock) & 0x7FFFFFFFFFF) * 20 < 500){
            // event is a retrigger!

            if (std::find(gtids.begin(), gtids.end(), gtid) != gtids.end()){
                // event is a retrigger and within our ROI ! Need to ignore it
                std::cout << "Found GTID " << gtid << " within ROI is retrigger." << std::endl;
                flagged_ignore_evs.push_back(gtid);
            }
        }

        // reset the previous event clock
        previousClock = clock50;
    }
    return 0;
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

void evaluate_discriminants(int run_number, float fv_cut, float z_cut, std::vector<int> &flagged_ignore_evs){
    /*
    Function loads in a run_by_run test .root ratds file and loops over all the events. 
    
    Every event that passes the FV, Z, Recon and energy cuts is filtered by energy into a pre-defined energy bin.
    
    The multisite (dLog(L), Fisher, IQR) and directionality discriminants are calculted using previously created PDFs a given 
    energy.

    A run_by_run ntuple file is returned containing the [energy, position, discriminants] of each event and for each energy.
    */

    // define input and output paths
    std::string input_fname = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis/extracted_ratds/" + std::to_string(run_number) + "*.root";    
    std::string output_fname = "/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis3/processed_dataset/" + std::to_string(run_number) + ".root";
    
    // load the PDF files for each isotope
    std::string working_directory = "/data/snoplus3/hunt-stokes/multisite_clean/mc_studies";
    std::string pdf_B8_fname      = working_directory + "/run_by_run_pdf/full_analysis2_B8_solar_nue/total.root";
    std::string pdf_Tl208_fname   = working_directory + "/run_by_run_pdf/full_analysis2_Tl208/total.root";
    TFile *PDF_B8 = new TFile(pdf_B8_fname.c_str());
    TFile *PDF_Tl208 = new TFile(pdf_Tl208_fname.c_str());

    TH1D *multi_pdf_B8_2p5_5p0    = dynamic_cast<TH1D*>(PDF_B8->Get("multi_2.5_5.0"));
    TH1D *multi_pdf_Tl208_2p5_5p0 = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_2.5_5.0"));
    TH1D *multi_pdf_B8_2p5_3p0    = dynamic_cast<TH1D*>(PDF_B8->Get("multi_2.5_3.0"));
    TH1D *multi_pdf_Tl208_2p5_3p0 = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_2.5_3.0"));
    TH1D *multi_pdf_B8_3p0_3p5    = dynamic_cast<TH1D*>(PDF_B8->Get("multi_3.0_3.5"));
    TH1D *multi_pdf_Tl208_3p0_3p5 = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_3.0_3.5"));
    TH1D *multi_pdf_B8_3p5_4p0    = dynamic_cast<TH1D*>(PDF_B8->Get("multi_3.5_4.0"));
    TH1D *multi_pdf_Tl208_3p5_4p0 = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_3.5_4.0"));
    TH1D *multi_pdf_B8_4p0_4p5    = dynamic_cast<TH1D*>(PDF_B8->Get("multi_4.0_4.5"));
    TH1D *multi_pdf_Tl208_4p0_4p5 = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_4.0_4.5"));
    TH1D *multi_pdf_B8_4p5_5p0    = dynamic_cast<TH1D*>(PDF_B8->Get("multi_4.5_5.0"));
    TH1D *multi_pdf_Tl208_4p5_5p0 = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_4.5_5.0"));

    // split in half analysis
    TH1D *multi_pdf_B8_2p5_3p75    = dynamic_cast<TH1D*>(PDF_B8->Get("multi_2.5_3.75"));
    TH1D *multi_pdf_Tl208_2p5_3p75 = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_2.5_3.75"));
    TH1D *multi_pdf_B8_3p75_5p0    = dynamic_cast<TH1D*>(PDF_B8->Get("multi_3.75_5.0"));
    TH1D *multi_pdf_Tl208_3p75_5p0 = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_3.75_5.0"));

    // high energy analysis
    TH1D *multi_pdf_B8_3p5_5p0    = dynamic_cast<TH1D*>(PDF_B8->Get("multi_3.5_5.0"));
    TH1D *multi_pdf_Tl208_3p5_5p0 = dynamic_cast<TH1D*>(PDF_Tl208->Get("multi_3.5_5.0"));

    // create the output ntuple with an output TTree containing discriminants for each energy PDF
    TFile *output_file    = new TFile(output_fname.c_str(), "RECREATE");

    TTree *info_2p5_5p0 = new TTree("2p5_5p0", "2p5_5p0");
    TTree *info_2p5_3p0 = new TTree("2p5_3p0", "2p5_3p0");
    TTree *info_3p0_3p5 = new TTree("3p0_3p5", "3p0_3p5");
    TTree *info_3p5_4p0 = new TTree("3p5_4p0", "3p5_4p0");
    TTree *info_4p0_4p5 = new TTree("4p0_4p5", "4p0_4p5");
    TTree *info_4p5_5p0 = new TTree("4p5_5p0", "4p5_5p0");

    TTree *info_2p5_3p75 = new TTree("2p5_3p75", "2p5_3p75");
    TTree *info_3p75_5p0 = new TTree("3p75_5p0", "3p75_5p0");

    TTree *info_3p5_5p0 = new TTree("3p5_5p0", "3p5_5p0");
    
    // define the variables to fill branches with
    double energy, x, y, z, dlogL, itr, HC_ratio, posFOM, cos_theta_sun;
    Int_t nhitsCleaned, nhitsRaw, gtid;
    UInt_t posFOM_hits;

    info_2p5_5p0->Branch("energy", &energy);
    info_2p5_5p0->Branch("itr", &itr);
    info_2p5_5p0->Branch("x", &x);
    info_2p5_5p0->Branch("y", &y);
    info_2p5_5p0->Branch("z", &z);
    info_2p5_5p0->Branch("dlogL", &dlogL);
    info_2p5_5p0->Branch("nhitsCleaned", &nhitsCleaned);
    info_2p5_5p0->Branch("nhitsRaw", &nhitsRaw);
    info_2p5_5p0->Branch("HC_ratio", &HC_ratio);
    info_2p5_5p0->Branch("posFOM", &posFOM);
    info_2p5_5p0->Branch("posFOM_hits", &posFOM_hits);
    info_2p5_5p0->Branch("gtid", &gtid);

    info_2p5_3p0->Branch("energy", &energy);
    info_2p5_3p0->Branch("itr", &itr);
    info_2p5_3p0->Branch("x", &x);
    info_2p5_3p0->Branch("y", &y);
    info_2p5_3p0->Branch("z", &z);
    info_2p5_3p0->Branch("dlogL", &dlogL);
    info_2p5_3p0->Branch("nhitsCleaned", &nhitsCleaned);
    info_2p5_3p0->Branch("nhitsRaw", &nhitsRaw);
    info_2p5_3p0->Branch("HC_ratio", &HC_ratio);
    info_2p5_3p0->Branch("posFOM", &posFOM);
    info_2p5_3p0->Branch("posFOM_hits", &posFOM_hits);
    info_2p5_3p0->Branch("gtid", &gtid);

    info_3p0_3p5->Branch("energy", &energy);
    info_3p0_3p5->Branch("itr", &itr);
    info_3p0_3p5->Branch("x", &x);
    info_3p0_3p5->Branch("y", &y);
    info_3p0_3p5->Branch("z", &z);
    info_3p0_3p5->Branch("dlogL", &dlogL);
    info_3p0_3p5->Branch("nhitsCleaned", &nhitsCleaned);
    info_3p0_3p5->Branch("nhitsRaw", &nhitsRaw);
    info_3p0_3p5->Branch("HC_ratio", &HC_ratio);
    info_3p0_3p5->Branch("posFOM", &posFOM);
    info_3p0_3p5->Branch("posFOM_hits", &posFOM_hits);
    info_3p0_3p5->Branch("gtid", &gtid);

    info_3p5_4p0->Branch("energy", &energy);
    info_3p5_4p0->Branch("itr", &itr);
    info_3p5_4p0->Branch("x", &x);
    info_3p5_4p0->Branch("y", &y);
    info_3p5_4p0->Branch("z", &z);
    info_3p5_4p0->Branch("dlogL", &dlogL);
    info_3p5_4p0->Branch("nhitsCleaned", &nhitsCleaned);
    info_3p5_4p0->Branch("nhitsRaw", &nhitsRaw);
    info_3p5_4p0->Branch("HC_ratio", &HC_ratio);
    info_3p5_4p0->Branch("posFOM", &posFOM);
    info_3p5_4p0->Branch("posFOM_hits", &posFOM_hits);
    info_3p5_4p0->Branch("gtid", &gtid);

    info_4p0_4p5->Branch("energy", &energy);
    info_4p0_4p5->Branch("itr", &itr);
    info_4p0_4p5->Branch("x", &x);
    info_4p0_4p5->Branch("y", &y);
    info_4p0_4p5->Branch("z", &z);
    info_4p0_4p5->Branch("dlogL", &dlogL);
    info_4p0_4p5->Branch("nhitsCleaned", &nhitsCleaned);
    info_4p0_4p5->Branch("nhitsRaw", &nhitsRaw);
    info_4p0_4p5->Branch("HC_ratio", &HC_ratio);
    info_4p0_4p5->Branch("posFOM", &posFOM);
    info_4p0_4p5->Branch("posFOM_hits", &posFOM_hits);
    info_4p0_4p5->Branch("gtid", &gtid);

    info_4p5_5p0->Branch("energy", &energy);
    info_4p5_5p0->Branch("itr", &itr);
    info_4p5_5p0->Branch("x", &x);
    info_4p5_5p0->Branch("y", &y);
    info_4p5_5p0->Branch("z", &z);
    info_4p5_5p0->Branch("dlogL", &dlogL);
    info_4p5_5p0->Branch("nhitsCleaned", &nhitsCleaned);
    info_4p5_5p0->Branch("nhitsRaw", &nhitsRaw);
    info_4p5_5p0->Branch("HC_ratio", &HC_ratio);
    info_4p5_5p0->Branch("posFOM", &posFOM);
    info_4p5_5p0->Branch("posFOM_hits", &posFOM_hits);
    info_4p5_5p0->Branch("gtid", &gtid);

    info_2p5_3p75->Branch("energy", &energy);
    info_2p5_3p75->Branch("itr", &itr);
    info_2p5_3p75->Branch("x", &x);
    info_2p5_3p75->Branch("y", &y);
    info_2p5_3p75->Branch("z", &z);
    info_2p5_3p75->Branch("dlogL", &dlogL);
    info_2p5_3p75->Branch("nhitsCleaned", &nhitsCleaned);
    info_2p5_3p75->Branch("nhitsRaw", &nhitsRaw);
    info_2p5_3p75->Branch("HC_ratio", &HC_ratio);
    info_2p5_3p75->Branch("posFOM", &posFOM);
    info_2p5_3p75->Branch("posFOM_hits", &posFOM_hits);
    info_2p5_3p75->Branch("gtid", &gtid);

    info_3p75_5p0->Branch("energy", &energy);
    info_3p75_5p0->Branch("itr", &itr);
    info_3p75_5p0->Branch("x", &x);
    info_3p75_5p0->Branch("y", &y);
    info_3p75_5p0->Branch("z", &z);
    info_3p75_5p0->Branch("dlogL", &dlogL);
    info_3p75_5p0->Branch("nhitsCleaned", &nhitsCleaned);
    info_3p75_5p0->Branch("nhitsRaw", &nhitsRaw);
    info_3p75_5p0->Branch("HC_ratio", &HC_ratio);
    info_3p75_5p0->Branch("posFOM", &posFOM);
    info_3p75_5p0->Branch("posFOM_hits", &posFOM_hits);
    info_3p75_5p0->Branch("gtid", &gtid);

    info_3p5_5p0->Branch("energy", &energy);
    info_3p5_5p0->Branch("itr", &itr);
    info_3p5_5p0->Branch("x", &x);
    info_3p5_5p0->Branch("y", &y);
    info_3p5_5p0->Branch("z", &z);
    info_3p5_5p0->Branch("dlogL", &dlogL);
    info_3p5_5p0->Branch("nhitsCleaned", &nhitsCleaned);
    info_3p5_5p0->Branch("nhitsRaw", &nhitsRaw);
    info_3p5_5p0->Branch("HC_ratio", &HC_ratio);
    info_3p5_5p0->Branch("posFOM", &posFOM);
    info_3p5_5p0->Branch("posFOM_hits", &posFOM_hits);
    info_3p5_5p0->Branch("gtid", &gtid);

    // setup deadtime veto code
    bool veto_flag          = false;
    bool lone_follower_flag = false;
    double loneFollowerTime = 0;
    double pileupTime       = 0;
    int num_vetos           = 0;
    ULong64_t veto_start_time, lone_start_time, dcApplied, dcFlagged, clock50;
    Int_t fPass, triggerWord;

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
            
            // get primary trigger
            const RAT::DS::EV& rEV = rDS.GetEV(0);

            gtid        = rEV.GetGTID();
            triggerWord = rEV.GetTrigType();
            
            // check if event has been flagged as a retrigger!
            if (std::find(flagged_ignore_evs.begin(), flagged_ignore_evs.end(), gtid) != (flagged_ignore_evs.end())){
                std::cout << "Event " << gtid << " flagged as retrigger --> will not evaluate multisite discriminant." << std::endl;
                continue;
            }

            // check if the triggerWord is 0 --> i.e. event is an orphan
            if (triggerWord == 0){
                // orphan event, so skip!
                std::cout << "Event " << gtid << " is an orphan --> will not evaluate multisite discriminant." << std::endl;
                continue;
            }

            // apply neck hotspot / weird stuff removal and deadtime cuts
            nhitsCleaned = rEV.GetNhitsCleaned();
            nhitsRaw     = rEV.GetNhits();
            HC_ratio     = float(nhitsCleaned) / float(nhitsRaw);
            posFOM       = rEV.GetFitResult("scintFitter").GetFOM("PositionLogL");
            posFOM_hits  = rEV.GetFitResult("scintFitter").GetFOM("PositionSelectedNHit");
            fPass        = rEV.GetDataCleaningFlags().GetLatestPass();
            clock50      = rEV.GetClockCount50();
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
            if (event_position_recon.Mag() > fv_cut){
                // failed FV cut!
                continue;
            }
            // std::cout << "Passed FV cut." << std::endl;
            if (z < z_cut){
                // not in the scintillator!
                continue;
            }
            // std::cout << "Passed position cuts." << std::endl;
            // convert back to PSUP Coords for time residual calculation
            event_position_recon.SetCoordinateSystem(fPSUPSystemId);

            // check that event energy fits into the large domain
            if (energy < 2.5 or energy > 5.0){
                // outside energy ROI!
                continue;
            }
            // std::cout << "Passed energy cut." << std::endl;
            // passed all event selection cuts! Now we evaluate the discriminants

            // get the solar direction for this run
            // RAT::DS::UniversalTime event_time = rEV.GetUniversalTime();
            // TVector3 solar_dir                = RAT::SunDirection(event_time.GetDays(), event_time.GetSeconds(), event_time.GetNanoSeconds()).Unit();

            // check result exists for the directionality
            bool dir_flag = false;

            // work out which PDF to use for this energy
            TH2D *pdf_dir;
            TH1D *pdf_B8_multi;
            TH1D *pdf_Tl208_multi;

            // changed here for the Bi214 study !! and changed back for full study ...
            if (energy >= 2.5 and energy < 3.0){
                // includes the Bi214 tagged ROI and the normal binning
                // pdf_dir         = dir_pdf_2p5_3p125;
                pdf_B8_multi    = multi_pdf_B8_2p5_3p0;
                pdf_Tl208_multi = multi_pdf_Tl208_2p5_3p0;  
            }
            if (energy >= 3.0 and energy < 3.5){
                // pdf_dir         = dir_pdf_3p125_3p75;
                pdf_B8_multi    = multi_pdf_B8_3p0_3p5;
                pdf_Tl208_multi = multi_pdf_Tl208_3p0_3p5;  
            }
            if (energy >= 3.5 and energy < 4.0){
                // pdf_dir         = dir_pdf_3p75_4p375;
                pdf_B8_multi    = multi_pdf_B8_3p5_4p0;
                pdf_Tl208_multi = multi_pdf_Tl208_3p5_4p0;  
            }
            if (energy >= 4.0 and energy < 4.5){
                // includes case for E > 5 MeV solar candidate study ...
                // pdf_dir         = dir_pdf_4p375_5p0;
                pdf_B8_multi    = multi_pdf_B8_4p0_4p5;
                pdf_Tl208_multi = multi_pdf_Tl208_4p0_4p5;  
            }
            if (energy >= 4.5 and energy <= 5.0){
                // includes case for E > 5 MeV solar candidate study ...
                // pdf_dir         = dir_pdf_4p375_5p0;
                pdf_B8_multi    = multi_pdf_B8_4p5_5p0;
                pdf_Tl208_multi = multi_pdf_Tl208_4p5_5p0;  
            }

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
            cos_theta_sun = 100; //directionality_fitter(solar_dir, event_position_recon, pmt_x, pmt_y, pmt_z, time_residuals, pdf_dir); // return false if no in time residuals
            // only calculate multisite discriminants if cos_theta_sun fit succeeded
            if (cos_theta_sun != -999.9){
                std::cout << "Evaluating multisite discriminant." << std::endl;
                
                dlogL         = multisite_discriminant(time_residuals, pdf_B8_multi, pdf_Tl208_multi);
                // IQR           = 0; // not implemented yet
                // fisher        = 0; // not implemented yet

                // fill the ntuples
                if (energy >= 2.5 and energy < 3.0){
                    info_2p5_3p0->Fill();
                }
                if (energy >= 3.0 and energy < 3.5){
                    info_3p0_3p5->Fill();
                }
                if (energy >= 3.5 and energy < 4.0){
                    info_3p5_4p0->Fill();
                }
                if (energy >= 4.0 and energy < 4.5){
                    info_4p0_4p5->Fill();
                }
                if (energy >= 4.5 and energy <= 5.0){
                    info_4p5_5p0->Fill();
                }
            }

            if (energy >= 2.5 and energy <= 5.0){
                // now repeat it using the full 2.5 --> 5.0 PDF
                cos_theta_sun = 100;//directionality_fitter(solar_dir, event_position_recon, pmt_x, pmt_y, pmt_z, time_residuals, dir_pdf_2p5_5p0);
                // only calculate multisite discriminants if cos_theta_sun fit succeeded
                if (cos_theta_sun != -999.9){
                    dlogL         = multisite_discriminant(time_residuals, multi_pdf_B8_2p5_5p0, multi_pdf_Tl208_2p5_5p0);
                    // IQR           = 0; // not implemented yet
                    // fisher        = 0; // not implemented yet

                    info_2p5_5p0->Fill();
                }
            }

            // split in half analysis result
            if (energy >= 2.5 and energy < 3.75){
                std::cout << "Repeating for 2.5 --> 3.75 MeV PDF." << std::endl;
                // now repeat it using the full 2.5 --> 5.0 PDF
                cos_theta_sun = 100;//directionality_fitter(solar_dir, event_position_recon, pmt_x, pmt_y, pmt_z, time_residuals, dir_pdf_2p5_5p0);
                // only calculate multisite discriminants if cos_theta_sun fit succeeded
                if (cos_theta_sun != -999.9){
                    std::cout << "Evaluating multisite discriminant." << std::endl;
                    dlogL         = multisite_discriminant(time_residuals, multi_pdf_B8_2p5_3p75, multi_pdf_Tl208_2p5_3p75);
                    // IQR           = 0; // not implemented yet
                    // fisher        = 0; // not implemented yet

                    info_2p5_3p75->Fill();
                }
            }
            // split in half analysis result
            if (energy >= 3.75 and energy <= 5.0){
                std::cout << "Repeating for 3.75 --> 5.0 MeV PDF." << std::endl;
                // now repeat it using the full 2.5 --> 5.0 PDF
                cos_theta_sun = 100;//directionality_fitter(solar_dir, event_position_recon, pmt_x, pmt_y, pmt_z, time_residuals, dir_pdf_2p5_5p0);
                // only calculate multisite discriminants if cos_theta_sun fit succeeded
                if (cos_theta_sun != -999.9){
                    std::cout << "Evaluating multisite discriminant." << std::endl;
                    dlogL         = multisite_discriminant(time_residuals, multi_pdf_B8_3p75_5p0, multi_pdf_Tl208_3p75_5p0);
                    // IQR           = 0; // not implemented yet
                    // fisher        = 0; // not implemented yet

                    info_3p75_5p0->Fill();
                }
            }

            // high E analysis result
            if (energy >= 3.5 and energy <= 5.0){
                std::cout << "Repeating for 3.75 --> 5.0 MeV PDF." << std::endl;
                // now repeat it using the full 2.5 --> 5.0 PDF
                cos_theta_sun = 100;//directionality_fitter(solar_dir, event_position_recon, pmt_x, pmt_y, pmt_z, time_residuals, dir_pdf_2p5_5p0);
                // only calculate multisite discriminants if cos_theta_sun fit succeeded
                if (cos_theta_sun != -999.9){
                    std::cout << "Evaluating multisite discriminant." << std::endl;
                    dlogL         = multisite_discriminant(time_residuals, multi_pdf_B8_3p5_5p0, multi_pdf_Tl208_3p5_5p0);
                    // IQR           = 0; // not implemented yet
                    // fisher        = 0; // not implemented yet

                    info_3p5_5p0->Fill();
                }
            }
        }
    }

    std::cout << "Writing results to file." << std::endl;
    
    // finished loop over all events and TTrees are filled --> write them to the output file
    output_file->cd();
    info_2p5_5p0->Write();
    info_2p5_3p0->Write();
    info_3p0_3p5->Write();
    info_3p5_4p0->Write();
    info_4p0_4p5->Write();
    info_4p5_5p0->Write();
    info_2p5_3p75->Write();
    info_3p75_5p0->Write();
    info_3p5_5p0->Write();
    output_file->Close();
}

int main(int argc, char* argv[]){

    int run_number           = std::stoi(argv[1]);
    double fv_cut            = std::stod(argv[2]);
    double z_cut             = std::stod(argv[3]);

    // find raw runtime of this run in seconds
    RAT::DB *db = RAT::DB::Get();
    db->LoadDefaults();
    RAT::DS::Run run;
    run.SetRunID(run_number);
    db->BeginOfRun(run);
    std::cout << "Run ID = " << run_number << std::endl;
    
    // check for and record retrigger events in this run
    int retrig_return;
    std::vector<int> flagged_ignore_evs;
    retrig_return = record_retriggers(run_number, flagged_ignore_evs);
    if (retrig_return == 1){
        // no ntuple files for this run, or no GTIDs found within ROI, so skip this job
        std::cout << "No files or GTIDs so terminating job." << std::endl;
        return 0;
    }

    // run the analysis!
    evaluate_discriminants(run_number, fv_cut, z_cut, flagged_ignore_evs);

    // create the job completed flag output
    std::ofstream completed("/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis3/job_completed_flag/" + std::to_string(run_number) + ".txt");
    completed << "1";
    completed.close();

    // keep a log of events tagged as retriggers
    std::ofstream retriggers("/data/snoplus3/hunt-stokes/multisite_clean/data_studies/extracted_data/full_analysis3/retrigger_logs/" + std::to_string(run_number) + ".txt");
    for (int retrig : flagged_ignore_evs){
        retriggers << retrig << std::endl;
    }
    retriggers.close();
}