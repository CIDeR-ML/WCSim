#include "WCSimVoxGen.hh"
#include "WCSimPrimaryGeneratorMessenger.hh"
#include "WCSimPrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "G4SPSEneDistribution.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4ParticleTable.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4RandomTools.hh" // For random number generation

#include "TH1D.h"
#include "TSpline.h"
#include <vector>
#include <cmath>
#include <random>

using namespace std;

G4int    WCSimVoxGen::nGammaOutcomes = 21;
//copied from WCTE mPMT QE-wavelength table in WCSimPMTObject.cc
G4double WCSimVoxGen::gammaWavelengths[21] = { 300., 320., 340., 360., 380., 400., 420., 440., 460., 480., 500., 520., 540., 560., 580., 600., 620., 640., 660., 680., 700.};
G4double WCSimVoxGen::correctionFactor = 1.0;
G4double WCSimVoxGen::gammaSpectrum[21] = { .0787*correctionFactor, .1838*correctionFactor, .2401*correctionFactor, .2521*correctionFactor, .2695*correctionFactor, .2676*correctionFactor,
                                            .2593*correctionFactor, .2472*correctionFactor, .2276*correctionFactor, .1970*correctionFactor,  .1777*correctionFactor, .1547*correctionFactor,
                                            .1033*correctionFactor, .0727*correctionFactor, .0587*correctionFactor, .0470*correctionFactor, .0372*correctionFactor, .0285*correctionFactor,
                                            .0220*correctionFactor, .0130*correctionFactor, .0084*correctionFactor};

G4double WCSimVoxGen::wavelength_binwidth = WCSimVoxGen::gammaWavelengths[1] - WCSimVoxGen::gammaWavelengths[0];

//only considered gamma for now - J.Xia, Aug 24, 2024
G4int    WCSimVoxGen::pdgids = 22;

WCSimVoxGen::WCSimVoxGen(WCSimDetectorConstruction* detector, G4double energy, G4double rinputs[], G4double phiinputs[], G4double zinputs[]) :
                                           myDetector(detector),
                                           gEnergy(energy),
                                           rRange({rinputs[0], rinputs[1]}),
                                           phiRange({phiinputs[0], phiinputs[1]}),
                                           zRange({zinputs[0], zinputs[1]})
{
  wcsimdir = string(getenv("WCSIM_BUILD_DIR"))+"data/";
  
  // Create a histogram of the gamma spectrum
  // find the left and right most bin edges:
  for (int i = 0; i < nGammaOutcomes+1; i++){
    hist_binedges[i] = gammaWavelengths[i] - wavelength_binwidth/2.;
  }
  hist_binedges[nGammaOutcomes] = gammaWavelengths[nGammaOutcomes-1] + wavelength_binwidth/2.;

  hSpectrum = new TH1D("hSpectrum", "hSpectrum", nGammaOutcomes, hist_binedges[0], hist_binedges[nGammaOutcomes]);
  for (int i = 0; i < nGammaOutcomes; i++){
    hSpectrum->SetBinContent(i+1, gammaSpectrum[i] / hSpectrum->GetBinCenter(i+1));
  }
  hSpectrum->Scale(1./hSpectrum->Integral());
  hSpectrum->SetDirectory(0);
  spline = new TSpline3(hSpectrum);
  // Initialise
  this->Initialise();
}

WCSimVoxGen::~WCSimVoxGen(){
  
  // This needed to be deleted
  delete myVoxGun;
  delete hSpectrum;
  delete spline;
}

void WCSimVoxGen::Initialise(){
    myVoxGun     = new G4ParticleGun();
    time = 0.;
}

G4double WCSimVoxGen::GenGammaEnergy(){
    G4double energy = 0.*MeV;
    // Generate a random energy from the gamma spectrum
    G4double rand = G4UniformRand();
    G4double prob = 0.;
    G4double spec_range = hist_binedges[nGammaOutcomes] - hist_binedges[0];
    G4int nstep = 500;
    for (G4int i = 0; i < nstep; i++){
        G4double curr_lambda = hist_binedges[0]+i*spec_range/nstep;
        prob += spline->Eval(curr_lambda);
        if (rand < prob){
	  //G4double rand_lambda = hSpectrum->GetBinWidth(i+1) * G4UniformRand() + hSpectrum->GetBinLowEdge(i+1);
	  //energy = 1.24E-3/rand_lambda*MeV;
  	    energy = 1.24E-3/curr_lambda*MeV;
            break;
        }
    }
    return energy;
}

void WCSimVoxGen::GenRandomPosition(){
     // Generate a random position for the particle
     // (r, phi, z)
     double rand = G4UniformRand();
     G4double r = (rRange[1] - rRange[0]) * rand + rRange[0];
     G4double cos_phi_rad = fabs(cos(phiRange[1]) - cos(phiRange[0])) * rand + min(cos(phiRange[0]), cos(phiRange[1])); //already in radians
     G4double sin_phi_rad = std::sqrt(1.0 - cos_phi_rad * cos_phi_rad);
     G4double z = (zRange[1] - zRange[0]) * rand + zRange[0];
     position = G4ThreeVector(r*cos_phi_rad, r*sin_phi_rad, z); // because the rotation later will flip the sign
     if (myDetector->GetIsNuPrism()){
       position.rotateX(-CLHEP::pi / 2);	
     }
}

void WCSimVoxGen::GenRandomDirection(){
     double rand = G4UniformRand();
     double phi = 2.0 * CLHEP::pi * rand;

     // Generate random cosine of the polar angle θ in [-1, 1]
     double cosTheta = 2.0 * rand - 1.0;
     double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);

     // Convert to Cartesian coordinates
     double x = sinTheta * std::cos(phi);
     double y = sinTheta * std::sin(phi); // because the rotation later will flip the sign
     double z = cosTheta;

     direction = G4ThreeVector(x, y, z);
     
     if (myDetector->GetIsNuPrism()){
       direction.rotateX(-CLHEP::pi / 2);
     }
}

void WCSimVoxGen::GenerateVox(G4Event* anEvent){
    myVoxGun->SetParticleEnergy(gEnergy);
    myVoxGun->SetParticleTime(time);
    myVoxGun->SetParticleDefinition(G4Gamma::Definition());

    GenRandomPosition();
    GenRandomDirection();
    // Configure the final properties of the particle
    myVoxGun->SetParticlePosition(position);
    myVoxGun->SetParticleMomentumDirection(direction);
    // Generate the primary vertex for the particle
    myVoxGun->GeneratePrimaryVertex(anEvent);
}
