#ifndef WCSimVoxGen_h
#define WCSimVoxGen_h

#include "WCSimDetectorMessenger.hh"

#include "G4SPSEneDistribution.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4RandomDirection.hh"

#include "globals.hh"
#include "jhfNtuple.h"
#include <fstream>

#include "TH1D.h"
#include "TSpline.h"
#include <vector>

using namespace std;

class G4ParticleGun;
class G4Event;

class WCSimVoxGen
{
  public:
    WCSimVoxGen(WCSimDetectorConstruction* detector, G4double energy, G4int nphotons, G4double rRange[], G4double phiRange[], G4double zRange[], G4double phidir = -999., G4double thetadir = -999.,
                                                                                                                                                 G4double phidirsmear = 0., G4double thetadirsmear = 0.);
    ~WCSimVoxGen();

    // Initialise the AmBe generator
    void Initialise();
    
    // Set Particle Gun Properties of the neutron and gamma
    void GenerateVox(G4Event* anEvent);
    //void SetPositionBGOGeometry(G4double X, G4double Y, G4double Z) { BGOX=X, BGOY=Y, BGOZ=Z; }
    //G4ThreeVector GetPositionBGOGeometry() const { return G4ThreeVector(BGOX, BGOY, BGOZ); }

    void GenRandomPosition();
    void GenRandomDirection();
    void GenRandomDirection_wRange();
    G4double GenGammaEnergy();
    void SetNPhotonsPerEvent(G4int nPhotons);


    void Set_Direction_flag(){apply_dir_range = true;}
    void SetGammaEnergy(G4double energy){gEnergy = energy;}
    G4ThreeVector GetPrimaryPosition(){ return position; }
    G4ThreeVector GetPrimaryDirection(){ return direction; }
    G4int GetPDG(){ return pdgids; }
    TH1D* GethSpectrum(){return hSpectrum;}

  private:
    G4ParticleGun*        myVoxGun;
    //G4SPSRandomGenerator* rGen;
    G4SPSEneDistribution* nEnergyDist;
    //G4SPSEneDistribution* nEnergyDistFE;
    //G4SPSEneDistribution* nEnergyDistSE;
    WCSimDetectorConstruction *myDetector;

    bool apply_dir_range;
    // Variables for the initialisation of Vox generator parameters
    G4double gEnergy;
    G4int nPhotons;
    // ranges of the voxels in the cylindrical coordinate system
    G4double rRange[2];
    G4double phiRange[2];
    G4double zRange[2];
    G4double phiDir;
    G4double thetaDir;
    G4double phidirRange;
    G4double thetadirRange;

    G4ThreeVector position;
    G4ThreeVector direction;

    static G4int nGammaOutcomes;
    static G4double correctionFactor;
    //static G4double gammaProbabilities[3];
    static G4double gammaWavelengths[23];
    static G4double gammaSpectrum[23];
    static G4int pdgids;
    static G4double wavelength_binwidth;

    G4double hist_binedges[24];
    //G4double epsilon;
    //G4double BGOX, BGOY, BGOZ;
    G4double time;
    // Variables for reading in the file
    /// Points to $WCSIM_BUILD_DIR/data/
    string wcsimdir;

    TH1D* hSpectrum;
    TSpline3* spline;
};

#endif
