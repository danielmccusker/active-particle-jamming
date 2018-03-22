/*
 
 TERNARY FLUID SIMULATION
 
 Originally written by Jack Treado for use on his senior thesis.
 
 COMPILE INSTRUCTIONS, plz read:
 
 The time evolution uses FFT (fast fourier transform) algorithms which are pretty hard to do in c++.
 In order to do FFTs, we use a package called fftw++, a library developed to use packing called FFTW3 in c++.
 The documentation for fftw++ is at http://fftwpp.sourceforge.net/
 The documentation for fftw3 is at http://www.fftw.org/
 In case you were wondering, FFTW stands for "fastest fourier transform in the west" ... yeah, painful, I know.
 
 FIRST THINGS FIRST! You need to install a few things from the fftw++ page in order to have the correct libraries on your home machine.
 Make sure the libraries are installed correctly, or else compilation is impossible.
 Luckily, the libraries are already on the Medusa cluster, but its good to have them for yourself too.
 
    IF YOU WANT TO COMPILE ON HOME MACHINE, you have to use the following command on the command line:
 
        g++ -I ../"headerlocation" ternaryFluidSimulation.cpp ../fftw++.cc -lfftw3
 
    - link to proper proper headers (accessed using #include right below this comment section) with -I command
    - also have to include the source file fftw++.cc after you include this source file
    - finally, link fftw3 library with -l command
 
    IF YOU WANT TO COMPILE ON MEDUSA CLUSTER, use:
    
        $ module load intel
        icpc -I ../header location ternaryFluidSimulation.cpp ../fftw++.cc -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
 
    - read docs on cluster first if you are unfamiliar with medusa, but long story short you need to compile with icpc, so you need to first load the intel module.
 
 
 
 Pro tip: edit on Xcode if on a mac, definitely use a wide wiewing window
 
 */


#include <stdio.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include "fftw++.h"
#include "Array.h"
#include <string>

#define N 50 // simulation grid is NxN in area

// namespaces used in FFTW package
using namespace std;
using namespace fftwpp;
using namespace utils;
using namespace Array;


// FUNCTION DEFINITIONS
void printArray(array2<double>, ofstream &);
void printXProfile(array2<double>, ofstream &);
double avgArray(array2<double>);
void packingFFT(array2<double>, array2<Complex>, array2<Complex>, rcfft2d&);
void unpackingFFT(array2<double>, array2<Complex>, array2<Complex>, crfft2d&);
void lapAndGrad(array2<Complex>, array2<Complex>, array2<Complex>, array2<Complex>,
                array2<Complex>, array2<Complex>, array2<Complex>, array2<Complex>);
void velocityField(array2<Complex>, array2<Complex>, array2<Complex>, array2<Complex>,
                   array2<Complex>, array2<Complex>, array2<Complex>, array2<Complex>);
void StructureFactor(array2<Complex>, ofstream&, array2<Complex>, ofstream&, double);
int delta(int, int);
double Oseen(int, int, double, double, double);


// CONSTANT DEFINITIONS
const double PI = 3.14159265359;   // why is this not in cmath

// chi parameters: see thesis for definition of critical chiIJ between chem species I and J
const double chiRN = 2.0;     // RNA - Nucleoplasm interaction
const double chiPN = 3.0;     // Protein - Nucleoplasm interaction
const double chiRP = 1.5;     // RNA - Protein interaction


// dimensionless quantities:
// NOTE time is scaled as t = tauD * tHat, so tHat is time we plot (as in units of tauD)
const double tauDVR = 1e2;  // protein diffusion/hydrodynamic ratio
const double tauDVP = 1e2;  // rna diffusion/hydrodynamic ratio
const double tauDT = 1.0;  // diffusion/transcription ratio
const double tauDP = 1.0;   // diffusion/protein ratio

//molecular weights
const double Nrna = 3;            // RNA molecular weight
const double Nprotein = 1;        // Protein molecular weight
const double Nnucleoplasm = 1;    // Nucleoplasm molecular weight
const double lambda = 3;          // surface energy term...scaled by kB*T*N^2
const double invNrn = ((1/(Nrna)) - (1/(Nnucleoplasm)));         // used in muR (chemical potential)
const double invNpn = ((1/(Nprotein)) - (1/(Nnucleoplasm)));     // used in muP (chemical potential)


int main()
{
    srand48(1000001);
    
    
    
    // Output file names: can create different strings using string variable definition
    
    string ms = "ms";
    string num = "-sum1-local";  // arbitrary, use to set up your file naming convention (I use ms* , number according to time scale choice)
    string underscore = "_";
    
    
    string dirLoc = "/Users/JackTreado/OlmstedGroup/summer/data/"; //comment out when not compiling on home machine
    //string dirLoc = "/home/jdt66/output/"; // (comment out when compiling on home machine)
    string groupName = ms + num + underscore;
    string rDataStr = "rData.dat";
    string rSFDataStr = "rSF.dat";
    string pDataStr = "pData.dat";
    string pSFDataStr = "pSF.dat";
    
    string groupLoc = dirLoc + groupName;
    string objRFile = groupLoc + rDataStr;
    string sfObjRFile = groupLoc + rSFDataStr;
    string objPFile = groupLoc + pDataStr;
    string sfObjPFile = groupLoc + pSFDataStr;
    
    
    ofstream objR(objRFile.c_str());
    ofstream sfObjR(sfObjRFile.c_str());
    
    ofstream objP(objPFile.c_str());
    ofstream sfObjP(sfObjPFile.c_str());
    
    
    double tEnd = 0;
    
    cout << "2D Multiple Species w/ Hydrodynamics and Structure Factor" << endl;
    cout << ms << num << endl;
    
    // fftw stuff
    fftw::maxthreads=get_max_threads();
    
    unsigned int Np = N/2 + 1;
    size_t align=sizeof(Complex);
    
    
    
    /*
     ------------------------------------------
     ARRAYS FOR FORWARD AND BACKWARD TRANSFORMS
     ------------------------------------------
     */
    
    
    // REAL SPACE ARRAYS
    // R is RNA, P is protein
    
    array2<double> R(N, N, align);
    array2<double> P(N, N, align);
    array2<double> muR(N, N, align);
    array2<double> muP(N, N, align);
    
    
    
    // RNA hydrodynamic real space arrays
    array2<double> L2R(N, N, align);
    array2<double> g1R(N, N, align);
    array2<double> g2R(N, N, align);
    array2<double> muGrad1R(N, N, align);
    array2<double> muGrad2R(N, N, align);
    
    array2<double> v1R(N, N, align);
    array2<double> v2R(N, N, align);
    array2<double> vgR(N, N, align);
    
    // Protein hydrodynamic real space arrays
    array2<double> L2P(N, N, align);
    array2<double> g1P(N, N, align);
    array2<double> g2P(N, N, align);
    array2<double> muGrad1P(N, N, align);
    array2<double> muGrad2P(N, N, align);
    
    array2<double> v1P(N, N, align);
    array2<double> v2P(N, N, align);
    array2<double> vgP(N, N, align);
    
    
    
    
    // Transcription localizer array set (real and FFT)
    array2<double> tSite(N, N, align);
    array2<Complex> F_tSite(N, Np, align);
    array2<Complex> F_tSite_full(N, N, align);
    rcfft2d ForwardT(N, N, tSite, F_tSite);
    crfft2d BackwardT(N, N, F_tSite, tSite);
    
    int x0 = N/2;  // x location, center of simulation
    int y0 = N/2;  // y location, center of simulation
    double sd = 1; // gaussian standard deviation
    
    // transcription-protein interaction
    array2<double> tp(N, N, align);
    array2<Complex> F_tp(N, Np, align);
    array2<Complex> F_tp_full(N, N, align);
    rcfft2d ForwardTP(N, N, tp, F_tp);
    crfft2d BackwardTP(N, N, F_tp, tp);
    
    // Nucleoplasm source
    array2<double> np(N, N, align);
    array2<Complex> F_np(N, Np, align);
    array2<Complex> F_np_full(N, N, align);
    rcfft2d Forwardnp(N, N, np, F_np);
    crfft2d Backwardnp(N, N, F_np, np);
    
    
    // FFTW3 arrays
    array2<Complex> F_R(N, Np, align);
    array2<Complex> F_P(N, Np, align);
    array2<Complex> F_muR(N, Np, align);
    array2<Complex> F_muP(N, Np, align);
    
    array2<Complex> F_L2R(N, Np, align);
    array2<Complex> F_g1R(N, Np, align);
    array2<Complex> F_g2R(N, Np, align);
    array2<Complex> F_muGrad1R(N, Np, align);
    array2<Complex> F_muGrad2R(N, Np, align);
    array2<Complex> F_v1R(N, Np, align);
    array2<Complex> F_v2R(N, Np, align);
    array2<Complex> F_vgR(N, Np, align);
    
    array2<Complex> F_L2P(N, Np, align);
    array2<Complex> F_g1P(N, Np, align);
    array2<Complex> F_g2P(N, Np, align);
    array2<Complex> F_muGrad1P(N, Np, align);
    array2<Complex> F_muGrad2P(N, Np, align);
    array2<Complex> F_v1P(N, Np, align);
    array2<Complex> F_v2P(N, Np, align);
    array2<Complex> F_vgP(N, Np, align);
    
    
    
    
    // Full Fourier Transform Arrays
    array2<Complex> F_R_full(N, N, align);
    array2<Complex> F_P_full(N, N, align);
    array2<Complex> F_muR_full(N, N, align);
    array2<Complex> F_muP_full(N, N, align);
    
    
    array2<Complex> F_L2R_full(N, N, align);
    array2<Complex> F_g1R_full(N, N, align);
    array2<Complex> F_g2R_full(N, N, align);
    array2<Complex> F_muGrad1R_full(N, N, align);
    array2<Complex> F_muGrad2R_full(N, N, align);
    array2<Complex> F_v1R_full(N, N, align);
    array2<Complex> F_v2R_full(N, N, align);
    array2<Complex> F_vgR_full(N, N, align);
    
    array2<Complex> F_L2P_full(N, N, align);
    array2<Complex> F_g1P_full(N, N, align);
    array2<Complex> F_g2P_full(N, N, align);
    array2<Complex> F_muGrad1P_full(N, N, align);
    array2<Complex> F_muGrad2P_full(N, N, align);
    array2<Complex> F_v1P_full(N, N, align);
    array2<Complex> F_v2P_full(N, N, align);
    array2<Complex> F_vgP_full(N, N, align);
    
    
    // fftw object instantiation
    // Forward is needed to transform to fourier space
    // Backward is needed for inverse transform
    
    
    rcfft2d ForwardR(N, N, R, F_R);
    rcfft2d ForwardP(N, N, P, F_P);
    rcfft2d ForwardmuR(N, N, muR, F_muR);
    rcfft2d ForwardmuP(N, N, muP, F_muP);
    
    
    rcfft2d Forward_L2R(N, N, L2R, F_L2R);
    rcfft2d Forward_g1R(N, N, g1R, F_g1R);
    rcfft2d Forward_g2R(N, N, g2R, F_g2R);
    rcfft2d Forward_muGrad1R(N, N, muGrad1R, F_muGrad1R);
    rcfft2d Forward_muGrad2R(N, N, muGrad2R, F_muGrad2R);
    rcfft2d Forward_v1R(N, N, v1R, F_v1R);
    rcfft2d Forward_v2R(N, N, v2R, F_v2R);
    rcfft2d Forward_vgR(N, N, vgR, F_vgR);
    
    rcfft2d Forward_L2P(N, N, L2P, F_L2P);
    rcfft2d Forward_g1P(N, N, g1P, F_g1P);
    rcfft2d Forward_g2P(N, N, g2P, F_g2P);
    rcfft2d Forward_muGrad1P(N, N, muGrad1P, F_muGrad1P);
    rcfft2d Forward_muGrad2P(N, N, muGrad2P, F_muGrad2P);
    rcfft2d Forward_v1P(N, N, v1P, F_v1P);
    rcfft2d Forward_v2P(N, N, v2P, F_v2P);
    rcfft2d Forward_vgP(N, N, vgP, F_vgP);
    
    
    
    crfft2d BackwardR(N, N, F_R, R);
    crfft2d BackwardP(N, N, F_P, P);
    crfft2d BackwardmuR(N, N, F_muR, muR);
    crfft2d BackwardmuP(N, N, F_muP, muP);
    
    crfft2d Backward_L2R(N, N, F_L2R, L2R);
    crfft2d Backward_g1R(N, N, F_g1R, g1R);
    crfft2d Backward_g2R(N, N, F_g2R, g2R);
    crfft2d Backward_muGrad1R(N, N, F_muGrad1R, muGrad1R);
    crfft2d Backward_muGrad2R(N, N, F_muGrad2R, muGrad2R);
    crfft2d Backward_v1R(N, N, F_v1R, v1R);
    crfft2d Backward_v2R(N, N, F_v2R, v2R);
    crfft2d Backward_vgR(N, N, F_vgR, vgR);
    
    crfft2d Backward_L2P(N, N, F_L2P, L2P);
    crfft2d Backward_g1P(N, N, F_g1P, g1P);
    crfft2d Backward_g2P(N, N, F_g2P, g2P);
    crfft2d Backward_muGrad1P(N, N, F_muGrad1P, muGrad1P);
    crfft2d Backward_muGrad2P(N, N, F_muGrad2P, muGrad2P);
    crfft2d Backward_v1P(N, N, F_v1P, v1P);
    crfft2d Backward_v2P(N, N, F_v2P, v2P);
    crfft2d Backward_vgP(N, N, F_vgP, vgP);
    
    
    
    /*
     ------------------------------------------------------------------------------------
     ------------------------------------------------------------------------------------
     
                                            Initialization
     
     ------------------------------------------------------------------------------------
     ------------------------------------------------------------------------------------
     */
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            P(i,j) = 0.3 + (drand48()/50);
            R(i,j) = 0.1 + (drand48()/50);
            
            tSite(i,j) = (1/(2*PI*sd*sd))*exp(-(((i-x0)*(i-x0))+((j-y0)*(j-y0))/(2*sd*sd))); // Normalized Gaussian
            
            // if gaussian is small, make 0
            // the gaussian localizes the transcription site, so it should only be non-zero near center of simulation
            if(tSite(i,j) < 1e-16)
                tSite(i,j) = 0;
            if(tSite(i,j) != tSite(i,j)){
                cout << "nan found in tSite at "<< i << " " << j << endl;
                cout << "i - x0 = " <<  i-x0 << endl;
                cout << "exponential = " << exp(-((i-x0))) << endl;
                return 0;
            }
        }
    }
    
    
    // total simulation time
    tEnd = 300;
    cout << "tEnd = " << tEnd << endl;
    
    // initialize variables in time evolution
    double kx = 0;
    double ky = 0;
    double kNorm = 0;
    double L2 = 0;
    double t = 0;
    int count = 0;
    double Nij = 0;  // 1 - R(i,j) - P(i,j), to be used in time evolution
    
    // Runge-Kutta variables (source: http://lpsa.swarthmore.edu/NumInt/NumIntFourth.html)
    double k1, k2, k3, k4;
    double y1, y2, y3, y4;
    
    // log variables
    double logR = 0;
    double logP = 0;
    double logNij = 0;
    
    // ------------------------------------ END INITIALIZATION ------------------------------------------
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*
     ----------------------------------------------------------------------------------------------------
     ----------------------------------------------------------------------------------------------------
     
                                            BEGIN TIME EVOLUTION
                                                    
     Sections:
     
     1. Adaptive Time Steps -- use different time steps to get good resolution at early times, coarse
     grain time step at later time steps to make more efficient
     
     2. Structure Factor Variables -- calculate structure factor
     
     3. Actual time evolution
        i.   laplacians and gradients in fourier space
        ii.  chemical potentials in real space
        iii. hydrodynamics in fourier space
        iv.  time evolution (euler scheme) in fourier space
     
     NOTE: printing of array, calculate the structure factor only happens with largest time step
     
     ----------------------------------------------------------------------------------------------------
     ----------------------------------------------------------------------------------------------------

     */
    
    
     
     // BEGIN ADAPTIVE TIME STEP INITIALIZATION
    
    
    
    int Nsteps = 3;             // number of different time steps
    int framesPerStep[Nsteps];  // number of total frames (outputs) per time step size
    double tPerFrame[Nsteps];   // number of itterations per frame
    int tPerStep[Nsteps];       // total amount of time given a single time step size
    double tDelta[Nsteps];      // itteration size
    int elapsed[Nsteps+1];      // amount of time elapsed, used in printing array
    int totalFrameNumber = 0;   // total number of frames
    
    
    
    // changing these values determines how long each time step is used
    tPerStep[0] = 2;
    tPerStep[1] = 8;
    tPerStep[2] = tEnd - (tPerStep[0] + tPerStep[1]);
    
    int transTime = tPerStep[2]/5;     // time at which transcription is turned on
    double trans = 0;           // switch to hold off on transcription until 3rd time step stage
    
    tDelta[0] = 0.0005; // early time step: small to see spinodal texture form
    tDelta[1] = 0.001; // median time step
    tDelta[2] = 0.0025;  // late time step: resolution isnt as important, so you can coarse grain the time step
    elapsed[0] = 0;
    
    
    for (int step = 0; step < Nsteps; step++) {
        // Earlier time steps
        if (step < Nsteps - 1) {
            tPerFrame[step] = 0.1;
            framesPerStep[step] = tPerStep[step]/tPerFrame[step];
        }
        // final time step
        else {
            tPerFrame[step] = 5;
            framesPerStep[step] = tPerStep[step]/tPerFrame[step];
        }
        
        elapsed[step+1] = elapsed[step] + tPerStep[step];
        
        totalFrameNumber += framesPerStep[step];
    }
    
    
    // END ADAPTIVE TIME STEP INITIALIZATION
    
    
    
    
    
    // BEGIN STRUCTURE FACTOR VARIABLES
    
    double da = ((PI)/N) + 0.001; // ensures that da gets good statistics, sqrt(2) takes radial distance into account
    int indMax = floor((PI + da)/da);
    int sfPrint = tPerFrame[Nsteps - 1];
    int sfPrintNumber = framesPerStep[Nsteps-1];
    
    // END STRUCTURE FACTOR VARIABLES
    
    
    // BEGIN Output system size and frame #
    
    objR << N << endl;
    objR << totalFrameNumber << endl;
    
    objP << N << endl;
    objP << totalFrameNumber << endl;
    
    sfObjR << sfPrintNumber << endl;
    sfObjP << sfPrintNumber << endl;
    sfObjR << indMax << endl;
    sfObjP << indMax << endl;
    
    // END initial outputs
    
    
    
    
    // BEGIN ACTUAL TIME EVOLUTION w/ adaptive time step
    
    for(int step = 0; step < Nsteps; step++){
        
        // BEGIN NEW TIME STEP SIZE
        
        cout << "\n\n------------------------------------------\n";
        cout << "------------------------------------------\n";
        cout << "Starting a new step!" << endl;
        cout << "step = " << step << endl;
        cout << "time per step = " << tPerStep[step] << endl;
        count = 0;
        
        while (count < framesPerStep[step]){
            
            // BEGIN NEW FRAME (output step when on last time step)
            t = 0;
            
            // print every frame itteration
            printArray(R, objR);
            printArray(P, objP);
            
            /*
            
            if(step > Nsteps-1){
                
                // print info for structure factor
                // if printing for movie purposes, should be commented out
                
                
                packingFFT(R, F_R, F_R_full, ForwardR);
                packingFFT(P, F_P, F_P_full, ForwardR);
                StructureFactor(F_R_full, sfObjR, F_P_full, sfObjP, da);
                unpackingFFT(R, F_R, F_R_full, BackwardR);
                unpackingFFT(P, F_P, F_P_full, BackwardR);
             
                
            }
            
             */
            
            cout << "t = " << elapsed[step] + count*tPerFrame[step] << endl;
            cout << "tDelta = " << tDelta[step] << endl;
            cout << "Average RNA volume fraction = " << avgArray(R) << endl;
            cout << "Average RNA velocity field = " << avgArray(vgR) << endl;
            cout << "Average Protein volume fraction = " << avgArray(P) << endl;
            cout << "Average Protein velocity field = " << avgArray(vgP) << endl;
            
            while (t < tPerFrame[step]) {
                //packingFFT(L2R, F_L2R, F_L2R_full, Forward_L2R);
                //packingFFT(L2P, F_L2P, F_L2P_full, Forward_L2P);
                
                // compute laplacians, gradients in k spce and transform back
                lapAndGrad(F_L2R_full, F_g1R_full, F_g2R_full, F_R_full,
                           F_L2P_full, F_g1P_full, F_g2P_full, F_P_full);
                
                unpackingFFT(L2R, F_L2R, F_L2R_full, Backward_L2R);
                unpackingFFT(g1R, F_g1R, F_g1R_full, Backward_g1R);
                unpackingFFT(g2R, F_g2R, F_g2R_full, Backward_g2R);
                
                unpackingFFT(L2P, F_L2P, F_L2P_full, Backward_L2P);
                unpackingFFT(g1P, F_g1P, F_g1P_full, Backward_g1P);
                unpackingFFT(g2P, F_g2P, F_g2P_full, Backward_g2P);
                
                for(int i = 0; i < N; i++){
                    for(int j = 0; j < N; j++){
                        Nij = 1 - R(i,j) - P(i,j);
                        logR = log(R(i,j));
                        logP = log(P(i,j));
                        logNij = log(Nij);
                        
                        if (Nij < 0 || R(i,j) < 0 || P(i,j) < 0) {
                            cout << "Negatives found" << endl;
                            cout << "Pij = " << P(i,j) << endl;
                            cout << "Rij = " << R(i,j) << endl;
                            cout << "Nij = " << Nij << endl;
                            return 0;
                        }
                        
                        /*
                        
                        if (Nij < 2*val){
                            logNij = log(val) + (1/val)*(Nij-val) - (1/(2*pow(val,2)))*pow(Nij-val,2);
                            //+ (1/(3*pow(val,3)))*pow(Nij-val,3) - (1/(4*pow(val,4)))*pow(Nij-val,4)
                            //+ (1/(5*pow(val,5)))*pow(Nij-val,5) - (1/(6*pow(val,6)))*pow(Nij-val,6) + (1/(4*pow(val,7)))*pow(Nij-val,7);
                        }
                        if (R(i,j) < 2*val){
                            logR = log(val) + (1/val)*(R(i,j)-val) - (1/(2*pow(val,2)))*pow(R(i,j)-val,2);
                            //+ (1/(3*pow(val,3)))*pow(R(i,j)-val,3) - (1/(4*pow(val,4)))*pow(R(i,j)-val,4)
                            //+ (1/(5*pow(val,5)))*pow(R(i,j)-val,5) - (1/(6*pow(val,6)))*pow(R(i,j)-val,6) + (1/(4*pow(val,7)))*pow(R(i,j)-val,7);
                        }
                        if (P(i,j) < 2*val){
                            logP = log(val) + (1/val)*(P(i,j)-val) - (1/(2*pow(val,2)))*pow(P(i,j)-val,2);
                            //+ (1/(3*pow(val,3)))*pow(P(i,j)-val,3) - (1/(4*pow(val,4)))*pow(P(i,j)-val,4)
                            //+ (1/(5*pow(val,5)))*pow(P(i,j)-val,5) - (1/(6*pow(val,6)))*pow(P(i,j)-val,6) + (1/(4*pow(val,7)))*pow(P(i,j)-val,7);
                        }
                        */
                         
                        /*
                        if (Nij != Nij || R(i,j) != R(i,j) || P(i,j) != P(i,j)){
                            cout << "nan found" << endl;
                            cout << "Pij = " << P(i,j) << endl;
                            cout << "Rij = " << R(i,j) << endl;
                            cout << "Nij = " << Nij << endl;
                            cout << "1 - Rij - Pij = " << 1 - R(i,j) - P(i,j) << endl;
                            return 0;
                        }
                        */
                        
                        muR(i,j) = logR/Nrna + invNrn - logNij/Nnucleoplasm + P(i,j)*(chiRP - chiPN) + chiRN*(1 - 2*R(i,j) - P(i,j)) - lambda*L2R(i,j);
                        muP(i,j) = logP/Nprotein + invNpn - logNij/Nnucleoplasm + R(i,j)*(chiRP - chiRN) + chiPN*(1 - R(i,j) - 2*P(i,j)) - lambda*L2P(i,j);
                        
                        muGrad1R(i,j) = muR(i,j)*g1R(i,j);
                        muGrad2R(i,j) = muR(i,j)*g2R(i,j);
                        
                        muGrad1P(i,j) = muP(i,j)*g1P(i,j);
                        muGrad2P(i,j) = muP(i,j)*g2P(i,j);
                        
                        
                        
                        tp(i,j) = P(i,j)*Nij*tSite(i,j); // real space Gamma*P for protein-transcription interaction term in time evolve
                        np(i,j) = Nij*tSite(i,j);
                        
                    }
                }
                
                // FFT of mu and phi
                packingFFT(R, F_R, F_R_full, ForwardR);
                packingFFT(P, F_P, F_P_full, ForwardP);
                packingFFT(muR, F_muR, F_muR_full, ForwardmuR);
                packingFFT(muP, F_muP, F_muP_full, ForwardmuP);
                packingFFT(tSite, F_tSite, F_tSite_full, ForwardT);
                packingFFT(np, F_np, F_np_full, Forwardnp);
                packingFFT(tp, F_tp, F_tp_full, ForwardTP);
                
                // FFT muGrads for velocity field
                packingFFT(muGrad1R, F_muGrad1R, F_muGrad1R_full, Forward_muGrad1R);
                packingFFT(muGrad2R, F_muGrad2R, F_muGrad2R_full, Forward_muGrad2R);
                
                packingFFT(muGrad1P, F_muGrad1P, F_muGrad1P_full, Forward_muGrad1P);
                packingFFT(muGrad2P, F_muGrad2P, F_muGrad2P_full, Forward_muGrad2P);
                
                // compute velocity field in k space
                velocityField(F_v1R_full, F_v2R_full, F_muGrad1R_full, F_muGrad2R_full,
                              F_v1P_full, F_v2P_full, F_muGrad1P_full, F_muGrad2P_full);
                
                // unpack velocity components into real space
                unpackingFFT(v1R, F_v1R, F_v1R_full, Backward_v1R);
                unpackingFFT(v2R, F_v2R, F_v2R_full, Backward_v2R);
                
                unpackingFFT(v1P, F_v1P, F_v1P_full, Backward_v1P);
                unpackingFFT(v2P, F_v2P, F_v2P_full, Backward_v2P);
                
                // compute hydrodynamic interaction field in real space
                for(int i = 0; i < N; i++){
                    for(int j = 0; j < N; j++){
                        vgR(i,j) = (v1R(i,j)*g1R(i,j) + v2R(i,j)*g2R(i,j));
                        vgP(i,j) = (v1P(i,j)*g1P(i,j) + v2P(i,j)*g2P(i,j));
                    }
                }
                
                // FFT hydrodynamic field for time evolution
                packingFFT(vgR, F_vgR, F_vgR_full, Forward_vgR);
                packingFFT(vgP, F_vgP, F_vgP_full, Forward_vgP);
                
                // BEGIN UPDATE EACH FOURIER COMPONENT OF F_R_full AND F_P_full
                for(int i = 0; i < N; i++){
                    for(int j = 0; j < N; j++){
                        // Wavevectors
                        kx = (PI*(i))/(N);
                        ky = (PI*(j))/(N);
                        
                        // If i and/or j is greater, wavevector has one or two negative components
                        if (i > N/2) {
                            kx = (PI*(i-N))/(N);
                        }
                        if (j > N/2) {
                            ky = (PI*(j-N))/(N);
                        }
                        
                        
                        kNorm = sqrt(kx*kx + ky*ky);     // wavevector norm (k^2)
                        L2 = -2*(1-cos(kNorm));          // Laplacian operator in k space
                        
                        
                        // Time update using explicit Euler method for A and B species
                        F_R_full(i,j) = F_R_full(i,j) + tDelta[step]*(L2*F_muR_full(i,j) + tauDT*F_np_full(i,j) + tauDP*F_tp_full(i,j) - tauDVR*F_vgR_full(i,j));
                        F_P_full(i,j) = F_P_full(i,j) + tDelta[step]*(L2*F_muP_full(i,j) - tauDVP*F_vgP_full(i,j));
                    }
                }
                // END UPDATE EACH FOURIER COMPONENT OF F_R_Full AND F_P_full
                
                
                // transform back to phi and mu
                unpackingFFT(R, F_R, F_R_full, BackwardR);
                unpackingFFT(muR, F_muR, F_muR_full, BackwardmuR);
                unpackingFFT(P, F_P, F_P_full, BackwardP);
                unpackingFFT(muP, F_muP, F_muP_full, BackwardmuP);
                unpackingFFT(tSite, F_tSite, F_tSite_full, BackwardT);
                
                // update time
                t = t + tDelta[step];
            }
            
            //update frame count
            count++;
        }
    }
    
    // print R and P one final time after all simulations
    printArray(R, objR);
    printArray(P, objP);
    
    // close files
    objR.close();
    objP.close();
    sfObjR.close();
    sfObjP.close();
    
    
    cout << "\nFOR THIS RUN: full time was " << tEnd << endl;
    cout << "chiRP = " << chiRP << endl;
    cout << "chiRN = " << chiRN << endl;
    cout << "chiPN = " << chiPN << endl;
    cout << "lambda = " << lambda << endl;
    cout << "tauDVP = " << tauDVP << endl;
    cout << "tauDVR = " << tauDVR << endl;
    cout << "tauDT = " << tauDT << endl;
    cout << "tauDP = " << tauDP << endl << endl;
    cout << "Protein Data: " << objPFile << endl;
    cout << "RNA Data: " << objRFile << endl;
    cout << "RNA Structure Factor Data: " << sfObjRFile << endl;
    cout << "Protein Structure Factor Data: " << sfObjPFile << endl;
    return 0;
}


//
//
//  *********    **    *    *******           **       **      *      *******   **    *
//  *            * *   *    *      *          * *     * *     * *        *      * *   *
//  ******       *  *  *    *       *         *  *   *  *    *****       *      *  *  *
//  *            *   * *    *      *          *   * *   *   *     *      *      *   * *
//  *********    *    **    *******           *    *    *  *       *  *******   *    **
//
//









/*
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
                    PRINT FUNCTIONS
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
 */

void printArray(array2<double> phi, ofstream &obj)
{
    cout << endl << "in printArray" << endl;
    obj << phi << endl;
    cout << "LEAVING printArray" << endl;
}

void printXProfile(array2<double> phi, ofstream &obj)
{
    cout << "in printXProfile" << endl;
    for(int i = 0; i < N; i++)
        obj << i << " " << phi(i, N/2) << endl;
    cout << "leaving printXprofile" << endl;
    
    
}

double avgArray(array2<double> phi)
{
    double sum = 0;
    double avg = 0;
    
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            sum += phi(i,j);
        }
    }
    
    avg = (sum)/(N*N);
    
    return avg;
}


/*
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
                    END PRINT FUNCTIONS
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
 */












/*
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
              PACKING AND UNPACKING FUNCTIONS
 
 Used to modify FFTW fourier transforms so they are NxN and
 can be used in time evolution. If this looks super con-
 fusing, its because the mapping between FFTW and a normal
 DFT is not obvious. If one wants to check out the docu-
 mentation on this crap, then head on over to 
 http://www.fftw.org/doc/Multi_002dDimensional-DFTs-of-Real-Data.html#Multi_002dDimensional-DFTs-of-Real-Data
 so you can learn about how FFTW actually does its business.
 
 IF you can figure out a simpler way to do this "packing",
 the simulation will be much faster and I will sing your praises
 forever. -JT
 
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
 */

void packingFFT(array2<double> phi, array2<Complex> F_phi, array2<Complex> F_phi_full, rcfft2d& Forward)
{
    //cout << " IN packingFFT" << endl;
    
    Forward.fft0(phi, F_phi);
    // cout << "Incomplete Fourier transform:\n" << F_phi << endl;
    int px = 0;
    int py = 0;
    
    
    
    // F_phi -> q, F_phi_full -> p
    // 1st: qx < N/2, qy < N/2 + 1
    for(int qx = 0; qx < N/2; qx++){
        for(int qy = 0; qy < (N/2) + 1; qy++){
            px = qx + (N/2);
            py = qy;
            F_phi_full(px, py) = F_phi(qx, qy);
        }
    }
    
    // 2nd: qx > N/2, qy < N/2 + 1
    for(int qx = N/2; qx < N; qx++){
        for(int qy = 0; qy < (N/2) + 1; qy++){
            px = qx - (N/2);
            py = qy;
            F_phi_full(px, py) = F_phi(qx, qy);
        }
    }
    
    int qx = 0;
    int qy = 0;
    // int pxTmp = 0;
    int pyTmp = 0;
    
    
    // 3rd qx < N/2 -> px > N/2: py > N/2 + 1, so py
    for(int px = 0; px < N/2 + 1; px++){
        for(int py = (N/2)+1; py < N; py++){
            if(px == 0){
                qy = N - py;
                qx = px + N/2;
                F_phi_full(px, py) = conj(F_phi(qx, qy));
            }
            else{
                qy = px;
                qx = py - N/2;
                F_phi_full(px, py) = F_phi(qx, qy);
            }
        }
    }
    
    for(int px = N/2 + 1; px < N; px++){
        for(int py = (N/2)+1; py < N; py++){
            qy = N - py;
            qx = N - px + (N/2);
            F_phi_full(px, py) = conj(F_phi(qx, qy));
        }
    }
    
    for(int px = N/2 + 2; px < N; px++){
        pyTmp = px;
        qx = (N - px)+ N/2;
        qy = N - pyTmp;
        if(qy < 0)
            cout << "found error" << endl;
        F_phi_full(px, pyTmp) = conj(F_phi(qx, qy));
    }
    // cout << endl << F_phi_full << endl;
    
    
    
    //cout << " Leaving packingFFT" << endl;
}

void unpackingFFT(array2<double> phi, array2<Complex> F_phi, array2<Complex> F_phi_full, crfft2d& Backward){
    //cout << " IN unpackingFFT" << endl;
    int qx = 0;
    int qy = 0;
    
    int Np = (N/2) + 1;
    for(int px = 0; px < N/2; px++){
        for(int py = 0; py < Np; py++){
            qx = px + N/2;
            qy = py;
            F_phi(qx, qy) = F_phi_full(px, py);
        }
    }
    
    for(int px = N/2; px < N; px++){
        for(int py = 0; py < Np; py++){
            qx = px - N/2;
            qy = py;
            F_phi(qx, qy) = F_phi_full(px, py);
        }
    }
    
    Backward.fft0Normalized(F_phi, phi);
    
    //cout << " Leaving unpackingFFT" << endl;
    
}

/*
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
              END PACKING AND UNPACKING FUNCTIONS
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
 */














/*
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
                    BEGIN MATH FUNCTIONS
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
 */


// homemand kronecker delta
int delta(int i, int j){
    if (i == j){
        return 1;
    }
    else
        return 0;
}


double Oseen(int i, int j, double kNorm, double ki, double kj){
    
    double value = 0;
    
    // 2D stokeslet
    
    if(kNorm == 0){
        value = 0;
        //cout << "kNorm = 0 at " << ki << " " << kj << endl;
    }
    else
        value = 1/((kNorm*kNorm))*(delta(i,j) - (ki*kj)/(kNorm*kNorm));
    
    
    //
    
    return value;
    
}



// function to compute fourier gradients and laplacians

void lapAndGrad(array2<Complex>F_L2R_full, array2<Complex> F_g1R_full, array2<Complex> F_g2R_full, array2<Complex> F_R_full,
                array2<Complex>F_L2P_full, array2<Complex> F_g1P_full, array2<Complex> F_g2P_full, array2<Complex> F_P_full){
    
    double kx = 0;
    double ky = 0;
    double kNorm = 0;
    double L2 = 0;
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            
            // Define wavevectors and norm
            kx = (PI*(i))/(N);
            ky = (PI*(j))/(N);
            
            // If i and/or j is greater, wavevector has one or two negative components
            
            if (i > N/2) {
                kx = (PI*(i-N))/(N);
            }
            if (j > N/2) {
                ky = (PI*(j-N))/(N);
            }
            
            kNorm = sqrt(kx*kx + ky*ky);
            L2 = -2*(1 - cos(kNorm));
            
            
            // trust me, this is an fourier-transformed gradient
            complex<double> com_one_R(-kx*F_R_full(i,j).imag(), kx*F_R_full(i,j).real());
            complex<double> com_two_R(-ky*F_R_full(i,j).imag(), ky*F_R_full(i,j).real());
            
            complex<double> com_one_P(-kx*F_P_full(i,j).imag(), kx*F_P_full(i,j).real());
            complex<double> com_two_P(-ky*F_P_full(i,j).imag(), ky*F_P_full(i,j).real());
            
            
            F_L2R_full(i,j) = L2*F_R_full(i,j);
            F_g1R_full(i,j) = com_one_R;
            F_g2R_full(i,j) = com_two_R;
            
            F_L2P_full(i,j) = L2*F_P_full(i,j);
            F_g1P_full(i,j) = com_one_P;
            F_g2P_full(i,j) = com_two_P;
            
        }
    }
    
}



// hydrodynamic coupling vector field in fourier space
void velocityField(array2<Complex> F_v1R_full, array2<Complex> F_v2R_full, array2<Complex> F_muGrad1R_full, array2<Complex> F_muGrad2R_full,
                   array2<Complex> F_v1P_full, array2<Complex> F_v2P_full, array2<Complex> F_muGrad1P_full, array2<Complex> F_muGrad2P_full){
    
    double kx = 0;
    double ky = 0;
    double kNorm = 0;
    
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            
            // Define wavevectors and norm
            kx = (PI*(i))/(N);
            ky = (PI*(j))/(N);
            
            // If i and/or j is greater, wavevector has one or two negative components
            if (i > N/2) {
                kx = (PI*(i-N))/(N);
            }
            if (j > N/2) {
                ky = (PI*(j-N))/(N);
            }
            
            kNorm = sqrt(kx*kx + ky*ky);
            
            
            F_v1R_full(i,j) = Oseen(0,0,kNorm,kx,kx)*F_muGrad1R_full(i,j) + Oseen(0,1,kNorm,kx,ky)*F_muGrad2R_full(i,j);
            F_v2R_full(i,j) = Oseen(1,0,kNorm,ky,kx)*F_muGrad1R_full(i,j) + Oseen(1,1,kNorm,ky,ky)*F_muGrad2R_full(i,j);
            
            F_v1P_full(i,j) = Oseen(0,0,kNorm,kx,kx)*F_muGrad1P_full(i,j) + Oseen(0,1,kNorm,kx,ky)*F_muGrad2P_full(i,j);
            F_v2P_full(i,j) = Oseen(1,0,kNorm,ky,kx)*F_muGrad1P_full(i,j) + Oseen(1,1,kNorm,ky,ky)*F_muGrad2P_full(i,j);
            
            
        }
    }
    
}

/*
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
                    END MATH FUNCTIONS
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
 */














/*
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
               BEGIN STRUCTURE FACTOR FUNCTION
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
 */


void StructureFactor(array2<Complex> F_R_full, ofstream& objR, array2<Complex> F_P_full, ofstream& objP, double da){
    
    double kx = 0;
    double ky = 0;
    double kNorm = 0;
    double sfR = 0;
    double sfP = 0;
    double count = 0;
    
    double a = 0;
    double aMax = PI/2 + da;
    int indMax = floor(aMax/da);
    
    for (int ind = 0; ind < indMax; ind++) {
        
        a = ind*da;
        sfR = 0;
        sfP = 0;
        
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                // Define wavevectors and norm
                kx = (PI*(i))/(N);
                ky = (PI*(j))/(N);
                
                // If i and/or j is greater, wavevector has one or two negative components
                if (i > N/2) {
                    kx = (PI*(i-N))/(N);
                }
                if (j > N/2) {
                    ky = (PI*(j-N))/(N);
                }
                
                kNorm = sqrt(kx*kx + ky*ky);
                
                
                if(a < kNorm && kNorm < a + da){
                    sfR += (F_R_full(i,j) * conj(F_R_full(i,j))).real();
                    sfP += (F_P_full(i,j) * conj(F_P_full(i,j))).real();
                }
            }
        }
        sfR = sfR/(N*N);
        sfP = sfP/(N*N);
        objR << a << " " << sfR << endl;
        objP << a << " " << sfP << endl;
    }
}

/*
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
              END STRUCTURE FACTOR FUNCTION
 
 ------------------------------------------------------------
 ------------------------------------------------------------
 
 */




