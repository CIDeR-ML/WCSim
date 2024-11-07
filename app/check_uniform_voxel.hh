#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <unordered_map>

#include "TTree.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"

#include "WCSimRootOptions.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootEvent.hh"
#include "G4ThreeVector.hh"

//Variables
char* rootfilename;
char* voxdefname;
int verbose = 0;
int all_Chi2 = 0;
int print_usage = 0;
int skip_plotting = 0;
TH1D* hSpec;
TH1D* hVox[3];
TH2D* hDir;

std::vector<string> keys = {"r0", "r1", "phi0", "phi1", "z0", "z1", "phidir0", "phidir1", "thdir0", "thdir1", "criterion", "wrapup_file"};
std::unordered_map<std::string, double> variables;
std::unordered_map<std::string, double> check_values;
std::unordered_map<std::string, bool> range_values;

//Functions
double KS_test(std::vector<double> testvec, double start, double stop);
double Chi2_test(TH1D *hbase, TH1D *h);
void Load_Voxel_Definition(const char* filename);

//copied from WCTE mPMT QE-wavelength table in WCSimPMTObject.cc
int nGammaOutcomes = 23;
double gammaWavelengths[23] = { 280., 300., 320., 340., 360., 380., 400., 420., 440., 460., 480., 500., 520., 540., 560., 580., 600., 620., 640., 660., 680., 700., 720.};
double correctionFactor = 1.0;
double gammaSpectrum[23] = { 0., .0787*correctionFactor, .1838*correctionFactor, .2401*correctionFactor, .2521*correctionFactor, .2695*correctionFactor, .2676*correctionFactor, .2593*correctionFactor, .2472*correctionFactor, .2276*correctionFactor, .1970*correctionFactor,  .1777*correctionFactor, .1547*correctionFactor, .1033*correctionFactor, .0727*correctionFactor, .0587*correctionFactor, .0470*correctionFactor, .0372*correctionFactor, .0285*correctionFactor, .0220*correctionFactor, .0130*correctionFactor, .0084*correctionFactor, 0.};

double wavelength_binwidth = gammaWavelengths[1] - gammaWavelengths[0];
double hist_binedges[24];

std::string wrapup_file;
