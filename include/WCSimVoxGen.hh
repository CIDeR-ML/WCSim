#ifndef WCSimVoxGen_h
#define WCSimVoxGen_h

#include "WCSimDetectorMessenger.hh"

#include "G4SPSEneDistribution.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4RandomDirection.hh"

#include "globals.hh"
#include "jhfNtuple.h"
#include <fstream>

include <vector>

using namespace std;

class G4ParticleGun;
class G4Event;

class WCSimVoxGen
{
  public:
    WCSimVoxGen(WCSimDetectorConstruction* detector, G4int nGamma, G4double energy, G4double rRange[], G4double phiRange[], G4double zRange[]);
    ~WCSimVoxGen();

    // Initialise the AmBe generator
    void Initialise();
    
    // Set Particle Gun Properties of the neutron and gamma
    void GenerateVox(G4Event* anEvent);
    //void SetPositionBGOGeometry(G4double X, G4double Y, G4double Z) { BGOX=X, BGOY=Y, BGOZ=Z; }
    //G4ThreeVector GetPositionBGOGeometry() const { return G4ThreeVector(BGOX, BGOY, BGOZ); }

    std::vector<G4ThreeVector> GenRandomPosition();
    std::vector<G4ThreeVector> GenRandomDirection();
    G4double GenGammaEnergy();

  private:
    G4ParticleGun*        myVoxGun;
    G4SPSRandomGenerator* rGen;
    G4SPSEneDistribution* nEnergyDist;
    //G4SPSEneDistribution* nEnergyDistFE;
    //G4SPSEneDistribution* nEnergyDistSE;
    WCSimDetectorConstruction *myDetector;
    
    // Variables for the initialisation of Vox generator parameters
    G4double gEnergy;
    G4int nGamma;
    // ranges of the voxels in the cylindrical coordinate system
    G4double rRange[2];
    G4double phiRange[2];
    G4double zRange[2];

    std::vector<G4ThreeVector> pos_vox;
    std::vector<G4ThreeVector> dir_vox;

    static G4int nGammaOutcomes;
    static G4double gammaProbabilities[3];
    static G4int pdgids;
    G4double time;
    //G4double epsilon;
    //G4double BGOX, BGOY, BGOZ;

    // Variables for reading in the file
    /// Points to $WCSIM_BUILD_DIR/data/
    string wcsimdir;

};

#endif
