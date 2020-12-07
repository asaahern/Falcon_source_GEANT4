#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "PhysicsList.hh"
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif


#include "G4UnitsTable.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//~ namespace {
  //~ void PrintUsage() {
    //~ G4cerr << " Usage: " << G4endl;
    //~ G4cerr << " OpNovice [-m macro ] [-u UIsession] [-t nThreads] [-r seed] "
           //~ << G4endl;
    //~ G4cerr << "   note: -t option is available only for multi-threaded mode."
           //~ << G4endl;
  //~ }
//~ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  
  G4String macro;
  G4String session;
  #ifdef G4MULTITHREADED
    G4int nThreads = 0;
  #endif

  if (argc<3){
        
    G4cout<<" "<<G4endl;
    G4cout<<"YOU MUST PROVIDE THE RUN NUMBER."<<G4endl;
    G4cout<<"Ex: ./SCXM (int)runNumber (int)UseVIS"<<G4endl;
    G4cout<<"0 - NoVisualization"<<G4endl;
    G4cout<<"1 - WithVisualization"<<G4endl;
    G4cout<<" "<<G4endl;        
    exit(1);
    }

  G4int RunNumber = atoi(argv[1]);
  G4int useVis=atoi(argv[2]);

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4long seeds[2];
  G4long systime = time(NULL);
  std::cout<<systime<<std::endl;
  seeds[0] = (long) systime+RunNumber;
  seeds[1] = (long) (systime*G4UniformRand()+RunNumber);
  // seeds[1] = (long) (systime*randSeed);
  G4Random::setTheSeeds(seeds);

  // Construct the default run manager

  #ifdef G4MULTITHREADED
    G4MTRunManager * runManager = new G4MTRunManager;
    if ( nThreads > 0 ) runManager->SetNumberOfThreads(nThreads);
    #else
    G4RunManager * runManager = new G4RunManager;
  #endif

  // Seed the random number generator manually
  //~ G4Random::setTheSeed(myseed);

  
  // Set mandatory initialization classes

  // Detector construction
  DetectorConstruction* detector =new DetectorConstruction();
  runManager-> SetUserInitialization(detector);
 
  // Physics list
  runManager-> SetUserInitialization(new PhysicsList()); //(detector)
  
  // User action initialization
  auto actionInitialization = new ActionInitialization();//(detector)
  runManager->SetUserInitialization(actionInitialization);

  // Initialize G4 kernel
  runManager->Initialize();
  runManager->SetRunIDCounter(RunNumber);
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/control/execute physics.mac");
  
  #ifdef G4VIS_USE
  // Initialize visualization
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
 
  switch (useVis){
    case 0:{
      G4String command = "/control/execute noVIS.mac";
      UImanager->ApplyCommand(command);
      break;
    }
    case 1:{
      G4VisManager* visManager = new G4VisExecutive;  
      visManager->Initialize();
      G4UIExecutive * ui = new G4UIExecutive(argc,argv,session);
      UImanager->ApplyCommand("/control/execute vis.mac");
      UImanager->ApplyCommand("/control/execute gui.mac");
      ui->SessionStart();
      break;
    }
    default:{
      G4cout<<"Ex: ./SCXM (int)runNumber (int)UseVIS"<<G4endl;
      G4cout<<"0 - NoVisualization"<<G4endl;
      G4cout<<"1 - WithVisualization"<<G4endl;
      G4cout<<" "<<G4endl;
      exit(1);
      break;
    }
  }
  #endif


  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  #ifdef G4VIS_USE
  #endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......