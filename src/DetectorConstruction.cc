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
/// \file /src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"
#include "XenonSD.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"


#include "G4PhysicalConstants.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4SystemOfUnits.hh"
#include "Materials.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include <G4NistMaterialBuilder.hh>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
 {
	// World
	world_half_side = 100*mm;

	// FATGEM
	fatgem_radius   = 20.*mm;
	fatgem_thickness= 2.3*mm; // half thickness has to be given!!!
	hole_radius		= 1.5*mm;
	hole_pitch		= 5.*mm;
	hole_voff		= sqrt((hole_pitch*hole_pitch) - (hole_pitch/2*hole_pitch/2));

	// Photomultipliers
	PMT_radius     = 12.5*mm;
	PMT_thickness  = 2.5*mm;  // half thickness has to be given!!!
	phc_thickness  = 0.05*mm; // half thickness has to be given!!!

	// source
	source_side_x = 12.*mm;
	source_side_y = 6.*mm;
	source_side_z = 0.5*mm;
	kapton_thickness = 0.035*mm; // half thickness has to be given!!!
	al_thickness = 0.025*mm; // half thickness has to be given!!!

	// Drifts
	drift_1    		= 15.*mm;
	drift_2			= 20.5*mm;

	// Teflon base
	base_radius    = 50.*mm;

	// booleans
	disable_FATGEM_base = false;
	disable_PMT_base   = false;

 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{

	// Option to switch on/off checking of volumes overlaps
	G4bool checkOverlaps = true;
	G4NistManager* nist = G4NistManager::Instance();


	// ---------- Definition of materials -------------
	
	// Xenon ------------------------------------------
	// gas density
	G4double pressure = 10.*bar;
	G4double temperature = 298.15 * kelvin;
	G4Material* xe_mat = nist->ConstructNewGasMaterial("GXe", "G4_Xe", temperature, pressure);
	G4double density_xe = (pressure / atmosphere) * 131.29 / (temperature / kelvin * 82.058); // g/cm^3

	// refractive index calculations
	const G4int nXe_entries = 100;
	G4double nXe_energy[nXe_entries];
	G4double refractiveIndex_Xe[nXe_entries];
	G4double energy_start = 1.8 * eV;
	G4double energy_finish = 8.9 * eV;

	for (G4int i = 0; i < nXe_entries; i++) {
		nXe_energy[i] = (energy_start + i * (energy_finish - energy_start) / nXe_entries);


		// Formula for the refractive index taken from
		// A. Baldini et al., "Liquid Xe scintillation calorimetry 
		// and Xe optical properties", arXiv:physics/0401072v1 [physics.ins-det]

		// The Lorentz-Lorenz equation (also known as Clausius-Mossotti equation)
		// relates the refractive index of a fluid with its density:
		// (n^2 - 1) / (n^2 + 2) = - A · d_M,     (1)
		// where n is the refractive index, d_M is the molar density and
		// A is the first refractivity viral coefficient:
		// A(E) = \sum_i^3 P_i / (E^2 - E_i^2),   (2)
		// with:
		G4double P[3] = { 71.23, 77.75, 1384.89 }; // [eV^3 cm3 / mole]
		G4double E[3] = { 8.4, 8.81, 13.2 };       // [eV]

		// Note.- Equation (1) has, actually, a sign difference with respect 
		// to the one appearing in the reference. Otherwise, it yields values
		// for the refractive index below 1.

		// Let's calculate the virial coefficient.
		// We won't use the implicit system of units of Geant4 because
		// it results in loss of numerical precision.

		G4double energy_ite = nXe_energy[i] / eV;

		G4double virial = 0.;
		for (G4int j = 0; j < 3; j++)
			virial = virial + P[j] / (energy_ite * energy_ite - E[j] * E[j]);

		G4double mol_density = density_xe / 131.29;
		G4double alpha = virial * mol_density;

		// Isolating now the n2 from equation (1) and taking the square root
		refractiveIndex_Xe[i] = (1. - 2 * alpha) / (1. + alpha);

		if (refractiveIndex_Xe[i] < 1.) {
			// "Non-physical refractive index for energy "
			refractiveIndex_Xe[i] = 1.;
		}

	}
	assert(sizeof(refractiveIndex_Xe) == sizeof(nXe_energy));

	G4MaterialPropertiesTable* Xe_MPT1 = new G4MaterialPropertiesTable();

	Xe_MPT1->AddProperty("RINDEX", nXe_energy, refractiveIndex_Xe, nXe_entries)
		->SetSpline(true);

	G4cout << "Xe G4MaterialPropertiesTable" << G4endl;

	xe_mat->SetMaterialPropertiesTable(Xe_MPT1);


	// PMT window ----------------------------------------
	G4Material* PMT_window_mat = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

	// define window refractive index (source: janis.com/Libraries/Window_Transmissions/FusedSilicaUVGrade_SiO2_TransmissionCurveDataSheet.sflb.ashx)
	G4double photonEnergy[] = { 7.3074*eV, 6.7149*eV, 6.2113*eV, 5.7941*eV, 4.4319*eV, 4.1121*eV, 3.4035*eV,
								3.0704*eV, 2.8505*eV, 2.2748*eV, 2.1141*eV, 2.108*eV, 1.9296*eV, 1.8928*eV, 1.441*eV };
	
	G4double refractiveIndex1[] = { 1.6150, 1.5750, 1.5500, 1.5337, 1.4940, 1.4872, 1.4745,
									1.4696, 1.4666, 1.4601, 1.4585, 1.4584, 1.4567, 1.4564, 1.4525 };

	const G4int nEntries = sizeof(photonEnergy) / sizeof(G4double);

	G4MaterialPropertiesTable* pmtMPT1 = new G4MaterialPropertiesTable();

	pmtMPT1->AddProperty("RINDEX", photonEnergy, refractiveIndex1, nEntries);
	PMT_window_mat->SetMaterialPropertiesTable(pmtMPT1);


	// PMT photocathode ----------------------------------
	G4Material* PMT_phc_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");


	// Plexiglass (for FATGEM and source) ----------------
	G4Material* PMMA_mat = nist->FindOrBuildMaterial("G4_PLEXIGLASS");


	// Teflon bases   ------------------------------------
	G4Material* teflon_mat = nist->FindOrBuildMaterial("G4_TEFLON");

	const G4int num3 = 4;
	// source: C. Silva et al. “Reflectance of Polytetrafluoroethylene (PTFE) for Xenon Scintillation Light”,
	// J. Appl. Phys. 107 (2010) 064902
	G4double ephoton_teflon[num3] = { 2.21*eV , 3.95*eV, 4.87*eV, 7.3*eV };
	G4double reflectivity_teflon[num3] = { 0.98, 0.93, 0.85, 0.61};
	//G4double reflectivity_teflon[num3] = { 1., 1., 1., 1.};

	// Kapton layer at source ---------------------------
	G4Material* kapton_mat = nist->FindOrBuildMaterial("G4_KAPTON");

	// Aluminum layer at source ---------------------------
	G4Material* al_mat = nist->FindOrBuildMaterial("G4_Al");


  // ----------- colors ------------
	G4VisAttributes * blue = new G4VisAttributes(G4Colour(0. ,0.8 ,0.9, 0.6));
	blue -> SetVisibility(true);
	blue -> SetForceSolid(true);
	G4VisAttributes * black = new G4VisAttributes(G4Colour(0.2 ,0.2 ,0.2));
	black -> SetVisibility(true);
	black -> SetForceSolid(true);
	G4VisAttributes * red = new G4VisAttributes(G4Colour(0.3 ,0.1 ,0., 0.5));
	red -> SetVisibility(true);
	red -> SetForceSolid(true);
	G4VisAttributes * white = new G4VisAttributes(G4Colour(0.95, 0.95, 0.95, 0.5));
	white -> SetVisibility(true);
	white -> SetForceSolid(true);
	G4VisAttributes * yellow = new G4VisAttributes(G4Colour(1. ,1. ,0.));
	yellow -> SetVisibility(true);
	yellow -> SetForceSolid(true);
	G4VisAttributes * gray = new G4VisAttributes(G4Colour(0.7 ,0.7 ,0.7));
	gray -> SetVisibility(true);
	gray -> SetForceSolid(true);



	// ------------- Volumes -------------

    // The experimental Hall
	G4Box* expHall_box = new G4Box("World", world_half_side, world_half_side, world_half_side);
	G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box, xe_mat, "World", 0, 0, 0);
	G4VPhysicalVolume* expHall_phys = new G4PVPlacement(0, G4ThreeVector(), expHall_log, "World", 0, false, 0, checkOverlaps);
/*
	// PMT
	// window
	G4Tubs* PMT_win = new G4Tubs("PMT_win", 0.*cm, PMT_radius, PMT_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* pmtWindow_log = new G4LogicalVolume(PMT_win, PMT_window_mat, "PMT_win", 0, 0, 0);
	G4VPhysicalVolume* pmtWindow_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -drift_2 -fatgem_thickness - PMT_thickness), pmtWindow_log, "PMT_win", expHall_log, false, 0, checkOverlaps);
	pmtWindow_log -> SetVisAttributes(blue);
	
	// photocathode
	G4Tubs* PMT_phc = new G4Tubs("PMT_phCa", 0.*cm, PMT_radius, phc_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* pmtPhc_log = new G4LogicalVolume(PMT_phc, PMT_phc_mat, "PMT_phCa", 0, 0, 0);
	G4VPhysicalVolume* pmtPhc_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -drift_2 -fatgem_thickness - 2*PMT_thickness - phc_thickness), pmtPhc_log, "PMT_phCa", expHall_log, false, 0, checkOverlaps);
	pmtPhc_log -> SetVisAttributes(black);
*/
	// Source
	G4Box* source = new G4Box("source", source_side_x, source_side_y, source_side_z);
	G4LogicalVolume* source_log = new G4LogicalVolume(source, PMMA_mat, "source", 0, 0, 0);
	G4VPhysicalVolume* source_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, drift_1 +fatgem_thickness +2*al_thickness +2*kapton_thickness +source_side_z), source_log, "source", expHall_log, false, 0, checkOverlaps);
	source_log -> SetVisAttributes(gray);

	// kapton layer
	G4Tubs* kap_layer = new G4Tubs("kap_layer", 0.*cm, base_radius, kapton_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* kap_layer_log = new G4LogicalVolume(kap_layer, kapton_mat, "kap_layer", 0, 0, 0);
	G4VPhysicalVolume* kap_layer_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, drift_1 +fatgem_thickness +2*al_thickness +kapton_thickness), kap_layer_log, "kap_layer", expHall_log, false, 0, checkOverlaps);
	kap_layer_log -> SetVisAttributes(red);	

	// aluminum layer
	G4Tubs* al_layer = new G4Tubs("al_layer", 0.*cm, base_radius, al_thickness, 0.*deg, 360.*deg);
	G4LogicalVolume* al_layer_log = new G4LogicalVolume(al_layer, al_mat, "al_layer", 0, 0, 0);
	G4VPhysicalVolume* al_layer_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, drift_1 +fatgem_thickness +al_thickness), al_layer_log, "al_layer", expHall_log, false, 0, checkOverlaps);
	al_layer_log -> SetVisAttributes(yellow);	

	// xenon volume
	G4Tubs* xe_volume = new G4Tubs("xe_volume", 0.*cm, base_radius, drift_1, 0.*deg, 360.*deg);
	G4LogicalVolume* xe_volume_log = new G4LogicalVolume(xe_volume, xe_mat, "xe_volume", 0, 0, 0);
	G4VPhysicalVolume* xe_volume_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, fatgem_thickness), xe_volume_log, "xe_volume", expHall_log, false, 0, checkOverlaps);
	xe_volume_log -> SetVisAttributes(blue);	

/*
	// FATGEM
	// volume
	G4Tubs* fatgem = new G4Tubs("fatgem", 0.*cm, fatgem_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* fatgem_log = new G4LogicalVolume(fatgem, PMMA_mat, "fatgem", 0, 0, 0);
	G4VPhysicalVolume* fatgem_phys = new G4PVPlacement(0, G4ThreeVector(), fatgem_log, "fatgem", expHall_log, false, 0, checkOverlaps);
	fatgem_log -> SetVisAttributes(red);

	// optical surface
	G4OpticalSurface* OS_fatgem = new G4OpticalSurface("fatgem_surface");
	G4LogicalSkinSurface* SS_fatgem = new G4LogicalSkinSurface("fatgem_surface", fatgem_log, OS_fatgem);

	OS_fatgem->SetType(dielectric_metal);
	OS_fatgem->SetFinish(ground);
	OS_fatgem->SetModel(glisur);
	OS_fatgem->SetPolish(0.4);

	const G4int num = 11;
	G4double ephoton_fatgem[num] = { 0.12*eV, 0.31*eV, 0.62*eV, 1.24*eV, 1.77*eV, 2.067*eV, 2.48*eV, 3.1*eV, 4.13*eV, 6.2*eV,  8.9*eV };
	G4double reflectivity_fatgem[num] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };

	G4MaterialPropertiesTable* R_fatgem = new G4MaterialPropertiesTable();
	R_fatgem->AddProperty("REFLECTIVITY", ephoton_fatgem, reflectivity_fatgem, num);
	OS_fatgem->SetMaterialPropertiesTable(R_fatgem);


	// holes
	//line 0
	G4Tubs* Hole1 = new G4Tubs("hole1", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole1_log = new G4LogicalVolume(Hole1, matXe, "hole1", 0, 0, 0);
	G4VPhysicalVolume* Hole1_phys = new G4PVPlacement(0, G4ThreeVector(-hole_pitch/2, 0, 0), Hole1_log, "hole1", fatgem_log, false, 0, checkOverlaps);
	
	G4Tubs* Hole2 = new G4Tubs("hole2", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole2_log = new G4LogicalVolume(Hole2, matXe, "hole2", 0, 0, 0);
	G4VPhysicalVolume* Hole2_phys = new G4PVPlacement(0, G4ThreeVector(hole_pitch/2, 0, 0), Hole2_log, "hole2", fatgem_log, false, 0, checkOverlaps);

	G4Tubs* Hole3 = new G4Tubs("hole3", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole3_log = new G4LogicalVolume(Hole3, matXe, "hole3", 0, 0, 0);
	G4VPhysicalVolume* Hole3_phys = new G4PVPlacement(0, G4ThreeVector(-hole_pitch*3/2, 0, 0), Hole3_log, "hole3", fatgem_log, false, 0, checkOverlaps);
	
	G4Tubs* Hole4 = new G4Tubs("hole4", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole4_log = new G4LogicalVolume(Hole4, matXe, "hole4", 0, 0, 0);
	G4VPhysicalVolume* Hole4_phys = new G4PVPlacement(0, G4ThreeVector(hole_pitch*3/2, 0, 0), Hole4_log, "hole4", fatgem_log, false, 0, checkOverlaps);

	G4Tubs* Hole5 = new G4Tubs("hole5", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole5_log = new G4LogicalVolume(Hole5, matXe, "hole1", 0, 0, 0);
	G4VPhysicalVolume* Hole5_phys = new G4PVPlacement(0, G4ThreeVector(-hole_pitch*5/2, 0, 0), Hole5_log, "hole5", fatgem_log, false, 0, checkOverlaps);
	
	G4Tubs* Hole6 = new G4Tubs("hole6", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole6_log = new G4LogicalVolume(Hole6, matXe, "hole6", 0, 0, 0);
	G4VPhysicalVolume* Hole6_phys = new G4PVPlacement(0, G4ThreeVector(hole_pitch*5/2, 0, 0), Hole6_log, "hole6", fatgem_log, false, 0, checkOverlaps);

	//line -1
	G4Tubs* Hole7 = new G4Tubs("hole7", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole7_log = new G4LogicalVolume(Hole7, matXe, "hole7", 0, 0, 0);
	G4VPhysicalVolume* Hole7_phys = new G4PVPlacement(0, G4ThreeVector(0, -hole_voff, 0), Hole7_log, "hole7", fatgem_log, false, 0, checkOverlaps);
	
	G4Tubs* Hole8 = new G4Tubs("hole8", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole8_log = new G4LogicalVolume(Hole8, matXe, "hole8", 0, 0, 0);
	G4VPhysicalVolume* Hole8_phys = new G4PVPlacement(0, G4ThreeVector(-hole_pitch, -hole_voff, 0), Hole8_log, "hole8", fatgem_log, false, 0, checkOverlaps);

	G4Tubs* Hole9 = new G4Tubs("hole9", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole9_log = new G4LogicalVolume(Hole9, matXe, "hole9", 0, 0, 0);
	G4VPhysicalVolume* Hole9_phys = new G4PVPlacement(0, G4ThreeVector(hole_pitch, -hole_voff, 0), Hole9_log, "hole9", fatgem_log, false, 0, checkOverlaps);
	
	G4Tubs* Hole10 = new G4Tubs("hole10", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole10_log = new G4LogicalVolume(Hole10, matXe, "hole10", 0, 0, 0);
	G4VPhysicalVolume* Hole10_phys = new G4PVPlacement(0, G4ThreeVector(-hole_pitch*2, -hole_voff, 0), Hole10_log, "hole10", fatgem_log, false, 0, checkOverlaps);

	G4Tubs* Hole11 = new G4Tubs("hole11", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole11_log = new G4LogicalVolume(Hole11, matXe, "hole11", 0, 0, 0);
	G4VPhysicalVolume* Hole11_phys = new G4PVPlacement(0, G4ThreeVector(hole_pitch*2, -hole_voff, 0), Hole11_log, "hole11", fatgem_log, false, 0, checkOverlaps);
	
	//line +1
	G4Tubs* Hole12 = new G4Tubs("hole12", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole12_log = new G4LogicalVolume(Hole12, matXe, "hole12", 0, 0, 0);
	G4VPhysicalVolume* Hole12_phys = new G4PVPlacement(0, G4ThreeVector(0, hole_voff, 0), Hole12_log, "hole12", fatgem_log, false, 0, checkOverlaps);
	
	G4Tubs* Hole13 = new G4Tubs("hole13", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole13_log = new G4LogicalVolume(Hole13, matXe, "hole13", 0, 0, 0);
	G4VPhysicalVolume* Hole13_phys = new G4PVPlacement(0, G4ThreeVector(-hole_pitch, hole_voff, 0), Hole13_log, "hole13", fatgem_log, false, 0, checkOverlaps);

	G4Tubs* Hole14 = new G4Tubs("hole14", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole14_log = new G4LogicalVolume(Hole14, matXe, "hole14", 0, 0, 0);
	G4VPhysicalVolume* Hole14_phys = new G4PVPlacement(0, G4ThreeVector(hole_pitch, hole_voff, 0), Hole14_log, "hole14", fatgem_log, false, 0, checkOverlaps);
	
	G4Tubs* Hole15 = new G4Tubs("hole15", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole15_log = new G4LogicalVolume(Hole15, matXe, "hole15", 0, 0, 0);
	G4VPhysicalVolume* Hole15_phys = new G4PVPlacement(0, G4ThreeVector(-hole_pitch*2, hole_voff, 0), Hole15_log, "hole10", fatgem_log, false, 0, checkOverlaps);

	G4Tubs* Hole16 = new G4Tubs("hole16", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole16_log = new G4LogicalVolume(Hole16, matXe, "hole16", 0, 0, 0);
	G4VPhysicalVolume* Hole16_phys = new G4PVPlacement(0, G4ThreeVector(hole_pitch*2, hole_voff, 0), Hole16_log, "hole16", fatgem_log, false, 0, checkOverlaps);

	//line-2
	G4Tubs* Hole17 = new G4Tubs("hole17", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole17_log = new G4LogicalVolume(Hole17, matXe, "hole17", 0, 0, 0);
	G4VPhysicalVolume* Hole17_phys = new G4PVPlacement(0, G4ThreeVector(-hole_pitch/2, -hole_voff*2, 0), Hole17_log, "hole17", fatgem_log, false, 0, checkOverlaps);
	
	G4Tubs* Hole18 = new G4Tubs("hole18", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole18_log = new G4LogicalVolume(Hole18, matXe, "hole18", 0, 0, 0);
	G4VPhysicalVolume* Hole18_phys = new G4PVPlacement(0, G4ThreeVector(hole_pitch/2, -hole_voff*2, 0), Hole18_log, "hole18", fatgem_log, false, 0, checkOverlaps);

	G4Tubs* Hole19 = new G4Tubs("hole19", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole19_log = new G4LogicalVolume(Hole19, matXe, "hole19", 0, 0, 0);
	G4VPhysicalVolume* Hole19_phys = new G4PVPlacement(0, G4ThreeVector(-hole_pitch*3/2, -hole_voff*2, 0), Hole19_log, "hole19", fatgem_log, false, 0, checkOverlaps);
	
	G4Tubs* Hole20 = new G4Tubs("hole20", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole20_log = new G4LogicalVolume(Hole20, matXe, "hole20", 0, 0, 0);
	G4VPhysicalVolume* Hole20_phys = new G4PVPlacement(0, G4ThreeVector(hole_pitch*3/2, -hole_voff*2, 0), Hole20_log, "hole20", fatgem_log, false, 0, checkOverlaps);

	//line+2
	G4Tubs* Hole21 = new G4Tubs("hole21", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole21_log = new G4LogicalVolume(Hole21, matXe, "hole21", 0, 0, 0);
	G4VPhysicalVolume* Hole21_phys = new G4PVPlacement(0, G4ThreeVector(-hole_pitch/2, hole_voff*2, 0), Hole21_log, "hole21", fatgem_log, false, 0, checkOverlaps);
	
	G4Tubs* Hole22 = new G4Tubs("hole22", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole22_log = new G4LogicalVolume(Hole22, matXe, "hole22", 0, 0, 0);
	G4VPhysicalVolume* Hole22_phys = new G4PVPlacement(0, G4ThreeVector(hole_pitch/2, hole_voff*2, 0), Hole22_log, "hole22", fatgem_log, false, 0, checkOverlaps);

	G4Tubs* Hole23 = new G4Tubs("hole23", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole23_log = new G4LogicalVolume(Hole23, matXe, "hole23", 0, 0, 0);
	G4VPhysicalVolume* Hole23_phys = new G4PVPlacement(0, G4ThreeVector(-hole_pitch*3/2, hole_voff*2, 0), Hole23_log, "hole23", fatgem_log, false, 0, checkOverlaps);
	
	G4Tubs* Hole24 = new G4Tubs("hole24", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole24_log = new G4LogicalVolume(Hole24, matXe, "hole24", 0, 0, 0);
	G4VPhysicalVolume* Hole24_phys = new G4PVPlacement(0, G4ThreeVector(hole_pitch*3/2, hole_voff*2, 0), Hole24_log, "hole24", fatgem_log, false, 0, checkOverlaps);

	//line -3
	G4Tubs* Hole25 = new G4Tubs("hole25", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole25_log = new G4LogicalVolume(Hole25, matXe, "hole25", 0, 0, 0);
	G4VPhysicalVolume* Hole25_phys = new G4PVPlacement(0, G4ThreeVector(0, -hole_voff*3, 0), Hole25_log, "hole25", fatgem_log, false, 0, checkOverlaps);
	
	G4Tubs* Hole26 = new G4Tubs("hole26", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole26_log = new G4LogicalVolume(Hole26, matXe, "hole26", 0, 0, 0);
	G4VPhysicalVolume* Hole26_phys = new G4PVPlacement(0, G4ThreeVector(-hole_pitch, -hole_voff*3, 0), Hole26_log, "hole26", fatgem_log, false, 0, checkOverlaps);

	G4Tubs* Hole27 = new G4Tubs("hole27", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole27_log = new G4LogicalVolume(Hole27, matXe, "hole27", 0, 0, 0);
	G4VPhysicalVolume* Hole27_phys = new G4PVPlacement(0, G4ThreeVector(hole_pitch, -hole_voff*3, 0), Hole27_log, "hole27", fatgem_log, false, 0, checkOverlaps);
	
	//line +3
	G4Tubs* Hole28 = new G4Tubs("hole28", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole28_log = new G4LogicalVolume(Hole28, matXe, "hole28", 0, 0, 0);
	G4VPhysicalVolume* Hole28_phys = new G4PVPlacement(0, G4ThreeVector(0, hole_voff*3, 0), Hole28_log, "hole28", fatgem_log, false, 0, checkOverlaps);
	
	G4Tubs* Hole29 = new G4Tubs("hole29", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole29_log = new G4LogicalVolume(Hole29, matXe, "hole29", 0, 0, 0);
	G4VPhysicalVolume* Hole29_phys = new G4PVPlacement(0, G4ThreeVector(-hole_pitch, hole_voff*3, 0), Hole29_log, "hole29", fatgem_log, false, 0, checkOverlaps);

	G4Tubs* Hole30 = new G4Tubs("hole30", 0, hole_radius, fatgem_thickness, 0. * deg, 360. * deg);
	G4LogicalVolume* Hole30_log = new G4LogicalVolume(Hole30, matXe, "hole30", 0, 0, 0);
	G4VPhysicalVolume* Hole30_phys = new G4PVPlacement(0, G4ThreeVector(hole_pitch, hole_voff*3, 0), Hole30_log, "hole30", fatgem_log, false, 0, checkOverlaps);
	
	Hole1_log -> SetVisAttributes(white);
	Hole2_log -> SetVisAttributes(white);
	Hole3_log -> SetVisAttributes(white);
	Hole4_log -> SetVisAttributes(white);
	Hole5_log -> SetVisAttributes(white);
	Hole6_log -> SetVisAttributes(white);
	Hole7_log -> SetVisAttributes(white);
	Hole8_log -> SetVisAttributes(white);
	Hole9_log -> SetVisAttributes(white);
	Hole10_log -> SetVisAttributes(white);
	Hole11_log -> SetVisAttributes(white);
	Hole12_log -> SetVisAttributes(white);
	Hole13_log -> SetVisAttributes(white);
	Hole14_log -> SetVisAttributes(white);
	Hole15_log -> SetVisAttributes(white);
	Hole16_log -> SetVisAttributes(white);
	Hole17_log -> SetVisAttributes(white);
	Hole18_log -> SetVisAttributes(white);
	Hole19_log -> SetVisAttributes(white);
	Hole20_log -> SetVisAttributes(white);
	Hole21_log -> SetVisAttributes(white);
	Hole22_log -> SetVisAttributes(white);
	Hole23_log -> SetVisAttributes(white);
	Hole24_log -> SetVisAttributes(white);
	Hole25_log -> SetVisAttributes(white);
	Hole26_log -> SetVisAttributes(white);
	Hole27_log -> SetVisAttributes(white);
	Hole28_log -> SetVisAttributes(white);
	Hole29_log -> SetVisAttributes(white);
	Hole30_log -> SetVisAttributes(white);
	
	// FATGEM base
	if (!disable_FATGEM_base) {

		// volume
		G4Tubs* FATGEM_base = new G4Tubs("FATGEM_base", fatgem_radius, base_radius, fatgem_thickness, 0. * deg, 360. * deg);
		G4LogicalVolume* FATGEM_base_log	= new G4LogicalVolume(FATGEM_base, teflon_mat, "FATGEM_base", 0, 0, 0);
		G4VPhysicalVolume* FATGEM_base_phys = new G4PVPlacement(0, G4ThreeVector(), FATGEM_base_log, "cathode", expHall_log, false, 0, checkOverlaps);
		FATGEM_base_log  -> SetVisAttributes(yellow);

		// optical surface
		G4OpticalSurface* op_FATGEM_base = new G4OpticalSurface("FATGEM_base_Surface");
		G4LogicalSkinSurface* FATGEM_base_Surface = new G4LogicalSkinSurface("FATGEM_base_Surface", FATGEM_base_log, op_FATGEM_base);

		G4MaterialPropertiesTable* FATGEM_base_ST2 = new G4MaterialPropertiesTable();
		FATGEM_base_ST2->AddProperty("REFLECTIVITY", ephoton_teflon, reflectivity_teflon, num3);

		op_FATGEM_base->SetMaterialPropertiesTable(FATGEM_base_ST2);
		op_FATGEM_base->SetType(dielectric_metal);
		op_FATGEM_base->SetFinish(ground);
		op_FATGEM_base->SetPolish(0.2);
		op_FATGEM_base->SetModel(glisur);
	}


	// PMT base
	if (!disable_PMT_base) {

		// volume
		G4Tubs* PMT_base = new G4Tubs("PMT_base", PMT_radius, base_radius, PMT_thickness, 0. * deg, 360. * deg);
		G4LogicalVolume* PMT_base_log	= new G4LogicalVolume(PMT_base, teflon_mat, "PMT_base", 0, 0, 0);
		G4VPhysicalVolume* PMT_base_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -drift_2 -fatgem_thickness - PMT_thickness), PMT_base_log, "PMT_base", expHall_log, false, 0, checkOverlaps);
		PMT_base_log  -> SetVisAttributes(yellow);

		// optical surface
		G4OpticalSurface* op_PMT_base = new G4OpticalSurface("PMT_base_Surface");
		G4LogicalSkinSurface* PMT_base_Surface = new G4LogicalSkinSurface("PMT_base_Surface", PMT_base_log, op_PMT_base);

		G4MaterialPropertiesTable* PMT_base_ST2 = new G4MaterialPropertiesTable();
		PMT_base_ST2->AddProperty("REFLECTIVITY", ephoton_teflon, reflectivity_teflon, num3);

		op_PMT_base->SetMaterialPropertiesTable(PMT_base_ST2);
		op_PMT_base->SetType(dielectric_metal);
		op_PMT_base->SetFinish(ground);
		op_PMT_base->SetPolish(0.2);
		op_PMT_base->SetModel(glisur);
	}


	// ------------- Surfaces --------------

	// PMT Quartz Window
	G4OpticalSurface* opPMT_window = new G4OpticalSurface("PMT_Surface");
	opPMT_window->SetType(dielectric_LUTDAVIS);
	opPMT_window->SetFinish(Polished_LUT);
	opPMT_window->SetModel(DAVIS);

	G4LogicalBorderSurface* PMT_Surface = new G4LogicalBorderSurface("PMT_BorderSurface1", pmtWindow_phys, expHall_phys, opPMT_window);
	G4OpticalSurface* PMTopticalSurface = dynamic_cast <G4OpticalSurface*>
		(PMT_Surface->GetSurface(pmtWindow_phys, expHall_phys)->
			GetSurfaceProperty());
	if (PMTopticalSurface) PMTopticalSurface->DumpInfo();

	// photocathode
	const G4int num2 = 2;
	G4double ephoton_abs[num2] = { 1.0 * eV , 10.0 * eV };
	G4double efficiency[num2] = { 1.0 , 1.0 };
	G4double reflectivity_ph[num2] = { 0.0, 0.0 };

	G4OpticalSurface* opPhCa = new G4OpticalSurface("PhCa_Surface");

	G4MaterialPropertiesTable* PhCaST2 = new G4MaterialPropertiesTable();
	PhCaST2->AddProperty("EFFICIENCY", ephoton_abs, efficiency, num2);
	PhCaST2->AddProperty("REFLECTIVITY", ephoton_abs, reflectivity_ph, num2);
	
	opPhCa->SetMaterialPropertiesTable(PhCaST2);
	opPhCa->SetType(dielectric_metal);
	opPhCa->SetFinish(polished);
	opPhCa->SetModel(glisur);

	G4LogicalSkinSurface* PhCa_Surface = new G4LogicalSkinSurface("PhCa_SkinSurface", pmtPhc_log, opPhCa);
*/

//always return the physical World
  return expHall_phys;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // ------------- Sensitive detectors ------------- 
  XenonSD* XeSD = new XenonSD("/XenonSD","XenonHitsCollection");  
  G4SDManager::GetSDMpointer()->AddNewDetector(XeSD);
  SetSensitiveDetector("xe_volume", XeSD, true); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......