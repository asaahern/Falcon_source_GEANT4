//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: SteppingAction.cc 71007 2013-06-09 16:14:59Z maire $
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include <G4UnitsTable.hh>

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4Timer.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
: G4UserSteppingAction()
{ 
  fEventNumber = -1;
}


SteppingAction::~SteppingAction()
{ ; }


void SteppingAction::UserSteppingAction(const G4Step* step)
{	

G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
                                              
if (eventNumber != fEventNumber) {
     fEventNumber = eventNumber;
}

//Tracks
G4Track* track = step->GetTrack();
	
//Runs
const G4Run* run = G4RunManager::GetRunManager()->GetCurrentRun();   //G4int runID = run->GetRunID();
	
// Events
const G4Event * event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
	

G4String ParticleName = track->GetDynamicParticle()->
                                 GetParticleDefinition()->GetParticleName();
G4String processName = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    
// ---------------- Primary Processes ------------------------
 
if (ParticleName == "gamma" && (processName=="compt" || processName=="phot")){ // !!!! No transport processes considered -->  processName!="Transportation"
	  
  	std::ofstream file;
  	G4String name="PrimaryProcesses_P"+ std::to_string(run->GetRunID())+ "bar.dat";
	
  	std::ifstream ifile(name.c_str());
  
  	if (ifile) {
  		// The file exists, and is open for input
  		ifile.close();
  		file.open(name.c_str(), std::ofstream::out | std::ofstream::app );
  
		file << processName <<", "<< event-> GetEventID() <<", "<< track->GetKineticEnergy() <<", "<<
			track->GetMomentumDirection().getTheta()/degree <<", "<< track->GetMomentumDirection().getPhi()/degree <<", "<<
			track->GetPosition().getX() <<", "<< track->GetPosition().getY() <<", "<< track->GetPosition().getZ() <<"\n";

		file.close();
	}
	
	
	else{
		file.open(name.c_str(), std::ofstream::out ); 
		file << "Process, eventID, Energy, Theta, Phi, X, Y, Z \n";
	
		file << processName <<", "<< event-> GetEventID() <<", "<< track->GetKineticEnergy() <<", "<<
			track->GetMomentumDirection().getTheta()/degree <<", "<< track->GetMomentumDirection().getPhi()/degree <<", "<<
			track->GetPosition().getX() <<", "<< track->GetPosition().getY() <<", "<< track->GetPosition().getZ() <<"\n";

    	file.close();
	}
}	  
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......