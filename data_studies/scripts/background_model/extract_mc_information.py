import numpy as np
import rat
import ROOT
from ROOT import RAT
import glob
import matplotlib.pyplot as plt
import argparse

def extract_mc_information(run_number, isotope):
    """
    Function loops over the Tl-208 MC ASCII simulated and extracts the: 

        1. Nhits
        2. Reconstructed Energy
        3. Reconstructed Position
        4. ITR
        5. Time Residuals
    """

    # save the reconstructed position, nhits, energy and ITR of each event
    posX           = []
    posY           = []
    posZ           = []
    posR           = []
    rho2           = []
    nhits_scaled   = []
    nhits_cleaned  = []
    recon_energy   = []
    itr            = []
    time_residuals = []
    counter        = 0

    for ientry, _ in rat.dsreader(f"/data/snoplus3/SNOplusData/production/miniProd_RAT-7-0-14_ASCI_RATHS_newRecoordination/ratds/simulation{isotope}_{run_number}.root"):
        
        # LPC setup and Point3D
        PMTCalStatus = RAT.DU.Utility.Get().GetPMTCalStatus()
        light_path = rat.utility().GetLightPathCalculator()
        group_velocity = rat.utility().GetGroupVelocity()
        pmt_info = rat.utility().GetPMTInfo()
        psup_system_id = RAT.DU.Point3D.GetSystemId("innerPMT")
        av_system_id = RAT.DU.Point3D.GetSystemId("av")
        stateCorrUtility = RAT.DU.Utility.Get().GetDetectorStateCorrection();

        # for this entry get the first triggered event
        if ientry.GetEVCount() == 0:
            continue
        reconEvent = ientry.GetEV(0)

        # apply normal reconstruction quality cuts
        fit_name = reconEvent.GetDefaultFitName()
        if not reconEvent.FitResultExists(fit_name):
            continue

        vertex = reconEvent.GetFitResult(fit_name).GetVertex(0)
        if (not vertex.ContainsPosition() or
            not vertex.ContainsTime() or
            not vertex.ValidPosition() or
            not vertex.ValidTime() or
            not vertex.ContainsEnergy() or
            not vertex.ValidEnergy()):
            continue

        # get the reconstructed variables
        reconPosition  = vertex.GetPosition() # returns in PSUP coordinates
        reconEnergy    = vertex.GetEnergy()        
        reconEventTime = vertex.GetTime()

        # apply the AV offset to the position
        event_point = RAT.DU.Point3D(psup_system_id, reconPosition)
        event_point.SetCoordinateSystem(av_system_id)
        radius = event_point.Mag()
        
        # save event position reconstruction in AV coordinates
        posX.append(event_point[0])
        posY.append(event_point[1])
        posZ.append(event_point[2])
        posR.append(radius)
        rho2.append((event_point[0]**2 + event_point[1]**2) / 6000**2)
        
        if radius > 6000:
            continue
        
        # convert back to PSUP coordinates
        event_point.SetCoordinateSystem(psup_system_id)

        # extract the time residuals for the event
        calibratedPMTs = reconEvent.GetCalPMTs()
        pmtCalStatus = rat.utility().GetPMTCalStatus()
        for j in range(calibratedPMTs.GetCount()):
            pmt = calibratedPMTs.GetPMT(j)
            if pmtCalStatus.GetHitStatus(pmt) != 0:
                continue
            
            pmt_point = RAT.DU.Point3D(psup_system_id, pmt_info.GetPosition(pmt.GetID()))
            light_path.CalcByPosition(event_point, pmt_point)
            inner_av_distance = light_path.GetDistInInnerAV()
            av_distance = light_path.GetDistInAV()
            water_distance = light_path.GetDistInWater()
            transit_time = group_velocity.CalcByDistance(inner_av_distance, av_distance, water_distance)
            residual_recon = pmt.GetTime() - transit_time - reconEventTime

            time_residuals.append(residual_recon)

        # save the events nhits, energy and ITR
        nhitscleaned = reconEvent.GetNhitsCleaned()
        itr_val       = reconEvent.GetClassifierResult("ITR:scintFitter").GetClassification("ITR")

        # apply the scaling factor to cleaned nhits to account for channel efficiency and coverage
        nhits_scaled_val  = stateCorrUtility.ApplyCorrectionNhits(nhitscleaned)
        nhits_scaled.append(nhits_scaled_val)
        recon_energy.append(reconEnergy)
        nhits_cleaned.append(nhitscleaned)
        itr.append(itr_val)
        counter += 1
        print("Completed: ", counter)
    
    # save the information for this run
    np.save(f"/data/snoplus3/hunt-stokes/clean_multisite/{isotope}_mc_information/position/posX_{run_number}.npy", posX)
    np.save(f"/data/snoplus3/hunt-stokes/clean_multisite/{isotope}_mc_information/position/posY_{run_number}.npy", posY)
    np.save(f"/data/snoplus3/hunt-stokes/clean_multisite/{isotope}_mc_information/position/posZ_{run_number}.npy", posZ)
    np.save(f"/data/snoplus3/hunt-stokes/clean_multisite/{isotope}_mc_information/position/posR_{run_number}.npy", posR)
    np.save(f"/data/snoplus3/hunt-stokes/clean_multisite/{isotope}_mc_information/position/rho2_{run_number}.npy", rho2)
    np.save(f"/data/snoplus3/hunt-stokes/clean_multisite/{isotope}_mc_information/energy/energy_{run_number}.npy", recon_energy)
    np.save(f"/data/snoplus3/hunt-stokes/clean_multisite/{isotope}_mc_information/itr/itr_{run_number}.npy", itr)
    np.save(f"/data/snoplus3/hunt-stokes/clean_multisite/{isotope}_mc_information/scaled_nhits/nhitsScaled_{run_number}.npy", nhits_scaled)
    np.save(f"/data/snoplus3/hunt-stokes/clean_multisite/{isotope}_mc_information/cleaned_nhits/nhitsCleaned_{run_number}.npy", nhits_cleaned)
    np.save(f"/data/snoplus3/hunt-stokes/clean_multisite/{isotope}_mc_information/time_residuals/residuals_{run_number}.npy", time_residuals)

parser = argparse.ArgumentParser()
parser.add_argument("run_number", type = int)
parser.add_argument("isotope", type = str)
args = parser.parse_args()

extract_mc_information(args.run_number, args.isotope)