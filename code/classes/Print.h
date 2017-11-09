// Class containing all printing functions

// Uses Unix commands mkdir and cp

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sys/types.h>
#include <dirent.h>

using namespace std;

struct Print
{
    Print();
    ~Print();
    
    void init(string, string, string, int);
    void print_data(long int, double, double, double, double, double, double, double, double, double, double);
    void print_Ovito(int, int, int, double, double, double, int);
    void print_size(double);
    void print_v(double, double);
    void print_g(double, double);
    void print_MSD(int, double);
    void print_VAF(int, double);
    void print_GNF(double, double, double);
    void print_orientationCorr(double,double);
    void resize(double &, double &);
    void print_summary(string, int, double, long int, int, double, double, double, double, long int, double);
    
    ofstream data, size, OvitoVid, pairCorr, velDist, MSD, GNF, VAF, orient, summary;

    string run;
    string path;
};

Print::Print(){
    // This space is intentionally left blank
}

Print::~Print(){
    data.close();
    OvitoVid.close();
    pairCorr.close();
    velDist.close();
    summary.close();
    GNF.close();
    MSD.close();
    VAF.close();
    orient.close();
}

void Print::init(string location, string fullRun, string ID, int noCells){
// Create file directories, open output files
    
    run = ID;
    path = location+"output/"+fullRun+"/";
    
    int check = 0;
    const char *p = path.c_str();
    DIR* dir = opendir(p);
    if(!dir){
        check += system(("mkdir "+path).c_str());
        check += system(("cp "+location+"code/postprocessing/{GNFvis.gnu,MSDvis.gnu,p*.py} "+path).c_str());
        if(check < 0){
            cout << "Failed to create folder for the full run. Status 715\n";
            exit(715);
        }
    } else if(dir) closedir(dir);
    

    check += system(("mkdir "+path+run).c_str());
    check += system(("mkdir "+path+run+"/dat").c_str());
    check += system(("mkdir "+path+run+"/eps").c_str());
    check += system(("mkdir "+path+run+"/vid").c_str());
    
    if(check < 0){
        cout << "Something hasn't been copied. Status 716\n";
        exit(716);
    }
    
    data.open((path+run+"/dat/data.dat").c_str());
    OvitoVid.open((path+run+"/vid/ovito.txt").c_str());
    GNF.open((path+run+"/dat/GNF.dat").c_str());
    MSD.open((path+run+"/dat/MSD.dat").c_str());
    pairCorr.open((path+run+"/dat/pairCorr.dat").c_str());
    velDist.open((path+run+"/dat/velDist.dat").c_str());
    VAF.open((path+run+"/dat/vaf.dat").c_str());
    orient.open((path+run+"/dat/orientCorr.dat").c_str());
    summary.open((path+run+"/dat/summary.dat").c_str());

}

void Print::print_data(long int t, double xa, double ya, double phi, double Psi,
                         double x1, double y1, double x2, double y2, double x3, double y3){
    data << setprecision(8) << t << "\t" << xa << "\t" << ya << "\t" << "\t" << phi << "\t" << Psi << "\t"
        << x1 << "\t" << y1 << "\t" << x2 << "\t" << y2 << "\t" << x3 << "\t" << y3 << "\n";
}

void Print::print_Ovito(int k, int noCells, int cellIndex, double x, double y, double radius, int overlap){
// Print in "XYZ" file format, to be read by molecular dynamics visualizing software
    
    if(k==0){
    // Demarcate the beginning of each time step
        OvitoVid << noCells << endl;
        OvitoVid << "time step comment" << endl;
    }
    
    OvitoVid << cellIndex << "\t" << x << "\t" << y << "\t" << radius << "\t" << overlap << endl;
}

void Print::print_v(double v, double prob){
    velDist << v << "\t" << prob << "\n";
}

void Print::print_g(double r, double gr){
    pairCorr << r << "\t" << gr << "\n";
}

void Print::print_MSD(int t, double msd){
    MSD << t << "\t" << msd << "\t" << log(msd) << endl;//"\t" << log(1.0-msd) << endl;
}

void Print::print_VAF(int t, double vaf){
    VAF << t << "\t" << vaf << endl;
}

void Print::print_GNF(double Rad, double avg, double fluct){
// Use logs for gnuplot log-log plot
    GNF << Rad << " " << avg << "\t" << fluct << "\t" << log(avg) << "\t" << log(fluct) << endl;
}

void Print::print_orientationCorr(double r, double g){
    orient << r << "\t" << g << endl;
}

void Print::print_summary(string ID, int noCells, double L, long int numberOfSteps, int stepsPerTime, double C1, double C2, double rho, double seconds, long int resetCounter, double peak, double decay) {
    summary << "Run ID:                     " << "\t" << ID << endl;
    summary << "Number of cells:            " << "\t" << noCells << endl;
    summary << "Grid length:                " << "\t" << L << endl;
    summary << "Number of steps:            " << "\t" << numberOfSteps << endl;
    summary << "Steps per unit time:        " << "\t" << stepsPerTime << endl;
    summary << "lambda_s (self-propulsion): " << "\t" << C1 << endl;
    summary << "lambda_n (noise):           " << "\t" << C2 << endl;
    summary << "Rho (density):              " << "\t" << rho << endl;
    summary << "Simulation time (seconds):  " << "\t" << seconds << endl;
    summary << "Number of list refreshes:   " << "\t" << resetCounter << endl;
    summary << "Nearest neighbor correlation: " << "\t" << peak << endl;
}
