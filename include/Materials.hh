#ifndef Materials_h
#define Materials_h 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"

class Materials
{
  public:

    virtual ~Materials();
 
    static Materials* GetInstance();

    G4Material* GetMaterial(const G4String);
	G4MaterialPropertiesTable* mptGlass;
	
  private:
 
    Materials();

    void CreateMaterials();
    

  private:
    static Materials* fInstance;  	
	
};

#endif /*Materials_h*/
