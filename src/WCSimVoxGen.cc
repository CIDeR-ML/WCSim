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
G4int    WCSimVoxGen::pdgids = 22;

WCSimVoxGen::WCSimVoxGen(WCSimDetectorConstruction* detector, G4int nGamma, G4double energy, G4double rRange[], G4double phiRange[], G4double zRange[]) :
                                           myDetector(detector),
                                           nGamma(nGamma),
                                           gEnergy(gammaEnergy),
                                           r(rRange),
                                           phi(phiRange),
                                           z(zRange)
  wcsimdir = string(getenv("WCSIMDIR_BUILD_DIR"))+"data/";

  // Initialise
  this->Initialise();
}

WCSimVoxGen::~WCSimVoxGen(){
  
  // This needed to be deleted
  delete myVoxGun;
  delete rGen;
  delete nEnergyDistGS;
  //delete nEnergyDistFE;
  //delete nEnergyDistSE;
}

void WCSimVoxGen::Initialise(){
    rGen          = new G4SPSRandomGenerator();
    myVoxGun     = new G4ParticleGun();

    //nEnergyDistGS = new G4SPSEneDistribution();
    //nEnergyDistGS->SetEnergyDisType("Arb");
    //nEnergyDistGS->ArbEnergyHistoFile(gs_path);
    //nEnergyDistGS->ArbInterpolate("Lin");
    //nEnergyDistGS->SetBiasRndm(rGen);

    time = 0.;
    //epsilon = 1e-6;
}

G4double WCSimVoxGen::GenGammaEnergy(){
    // Create a histogram of the gamma spectrum
    // find the left and right most bin edges:
    G4double wavelength_binwidth = gammaWavelengths[1] - gammaWavelengths[0];
    G4double hist_binedges[nGammaOutcomes+1];
    for (int i = 0; i < nGammaOutcomes+1; i++){
        hist_binedges[i] = gammaEnergies[i] - wavelength_binwidth/2.;
    }
    hist_binedges[nGammaOutcomes] = gammaEnergies[nGammaOutcomes-1] + wavelength_binwidth/2.;
    TH1D hSpectrum = new TH1D("hSpectrum", "hSpectrum", nGammaOutcomes, hist_binedges[0], hist_binedges[nGammaOutcomes]);
    for (int i = 0; i < nGammaOutcomes; i++){
        hSpectrum->SetBinContent(i+1, gammaSpectrum[i]);
    }
    hSpectrum->Scale(1./hSpectrum->Integral());
    // Generate a random energy from the gamma spectrum
    G4double energy = 0.*MeV;
    G4double rand = G4UniformRand();
    G4double prob = 0.;
    for (G4int i = 0; i < nGammaOutcomes; i++){
        prob += hSpectrum->GetBinContent(i+1);
        if (rand < prob){
            G4double rand_lambda = hSpectrum->GetBinWidth(i+1) * G4UniformRand() + hSpectrum->GetBinLowEdge(i+1);
            energy = 1.24E-3./rand_lambda*MeV;
            break;
        }
    }
    return energy;
}

std::vector<G4ThreeVector> WCSimVoxGen::GenRandomPosition(){
    for (G4int i=0; i<nGamma; i++){
        // Generate a random position for the particle
        // (r, phi, z)
        G4double r = (rRange[1] - rRange[0]) * G4UniformRand() + rRange[0];
        G4double phi = (phiRange[1] - phiRange[0]) * G4UniformRand() + phiRange[0];
        G4double z = (zRange[1] - zRange[0]) * G4UniformRand() + zRange[0];

        position = G4ThreeVector(r*cos(phi), r*sin(phi), z);
        if (myDetector->GetIsNuPrism()){
            position.rotateX(-90.*deg);
        }
        pos_vox.push_back(position);
    }

    return pos_vox;
}

std::vector<G4ThreeVector> WCSimVoxGen::GenRandomDirection(){
    for (G4int i=0; i<nGamma; i++){
        direction = rGen.fire();
        if (myDetector->GetIsNuPrism()){
            direction.rotateX(-90.*deg);
        }
        dir_vox.push_back(direction);
    }

    return dir_vox;
}

void WCSimVoxGen::GenerateVox(G4Event* anEvent){
    myVoxGun->SetParticleEnergy(gEnergy);
    myVoxGun->SetParticleTime(time);
    myVoxGun->SetParticleDefinition(G4Gamma::Definition());
    for (G4int iGamma = 0; iGamma < nGamma; iGamma++){
        // Configure the final properties of the particle
        myVoxGun->SetParticlePosition(pox_vox[iGamma]);
        myVoxGun->SetParticleMomentumDirection(dir_vox[iGamma]);
        // Generate the primary vertex for the particle
        myVoxGun->GeneratePrimaryVertex(anEvent);
    }
}
