
// Class containing all printing functions

// Uses Unix commands mkdir and cp

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sys/types.h>
#include <dirent.h>

struct Print
{
    Print(string, string, string, int);
    ~Print();
    
    void print_COM(long int, vector<double>&);
    void print_orientation(long int, vector<double>&);
    void print_order(long int, double);
    void print_velDist(double, double);
    void print_corr(double,double);
    void print_dens(double,double);
    void print_corrNorm(double,double);
    void print_pairCorr(double,double);
    void print_autoCorr(int, double);
    void print_MSD(int, double);
    void print_fluct(double, double, double);
    void print_Ovito(int&, int&, int&, double&, int&, vector<double>&, vector<double>&);
    void print_summary(string, int, double, long int, int, double, double, double, double, long int, double, double, double, double);
    
    ofstream COM, orientation, order,
             corr, corrNorm, pairCorr, autoCorr,
             velDist, MSD, fluct, dens,
             OvitoVid, summary, summary2;

    string run;
    string path;
};

Print::Print(string location, string fullRun, string ID, int noCells){
    run = ID;
    path = location+"output/"+fullRun+"/";
    
    int check = 0;
    const char *p = path.c_str();
    DIR* dir = opendir(p);
    if(!dir){
        check += system(("mkdir "+path).c_str());
        check += system(("cp "+location+"code/postprocessing/* "+path).c_str());
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
        cout << "Directories not successfully created. Status 716\n";
        exit(716);
    }
    
    COM.open((path+run+"/dat/COM.dat").c_str());
    orientation.open((path+run+"/dat/orientation.dat").c_str());
    order.open((path+run+"/dat/order.dat").c_str());
    corr.open((path+run+"/dat/corr.dat").c_str());
    corrNorm.open((path+run+"/dat/corrNorm.dat").c_str());
    pairCorr.open((path+run+"/dat/pairCorr.dat").c_str());
    velDist.open((path+run+"/dat/velDist.dat").c_str());
    autoCorr.open((path+run+"/dat/autoCorr.dat").c_str());
    fluct.open((path+run+"/dat/fluct.dat").c_str());
    dens.open((path+run+"/dat/densDist.dat").c_str());
    MSD.open((path+run+"/dat/MSD.dat").c_str());
    OvitoVid.open((path+run+"/vid/ovito.txt").c_str());
    summary.open((path+run+"/dat/summary.dat").c_str());
    summary2.open((path+run+"/dat/summary2.dat").c_str());
}

Print::~Print(){
    COM.close();
    orientation.close();
    order.close();
    corr.close();
    corrNorm.close();
    pairCorr.close();
    autoCorr.close();
    velDist.close();
    fluct.close();
    dens.close();
    MSD.close();
    OvitoVid.close();
    summary.close();
    summary2.close();
}

void Print::print_COM(long int t, vector<double> &center)
{
    if(NDIM == 2) COM << t << "\t" << center[0] << "\t" << center[1] << endl;
    if(NDIM == 3) COM << t << "\t" << center[0] << "\t" << center[1] << "\t" << center[2] << endl;
}

void Print::print_orientation(long int t, vector<double> &orient)
{
    if(NDIM==2) orientation << t << "\t" << orient[0] << "\t" << orient[1] << endl;
    if(NDIM==3) orientation << t << "\t" << orient[0] << "\t" << orient[1] << "\t" << orient[2] << endl;
}

void Print::print_order(long int t, double o)
{
    order << t << "\t" << o << endl;
}

void Print::print_Ovito(int &k, int &noCells, int &cellIndex, double &radius, int &overlap,
                        vector<double> &x, vector<double> &v){
// Print in "XYZ" file format, to be read by molecular dynamics visualization software
    
    if(k==0){
    // Demarcate the beginning of each time step
        OvitoVid << noCells << endl;
        OvitoVid << "time step comment" << endl;
    }
    
    if(NDIM==2) {
        OvitoVid << cellIndex << "\t" << radius << "\t" << overlap << "\t" << x[0] << "\t"
                 << x[1] << "\t" << v[0] << "\t" << v[1] << endl;
    }
    
    if(NDIM==3) {
        OvitoVid << cellIndex << "\t" << radius << "\t" << overlap << "\t" << x[0] << "\t" << x[1]
                 << "\t" << x[2] << "\t" << v[0] << "\t" << v[1] << "\t" << v[2] << endl;
    }
}

void Print::print_corr(double r, double v){
    corr << r << "\t" << v << endl;
}

void Print::print_corrNorm(double r, double v){
    corrNorm << r << "\t" << v << endl;
}

void Print::print_pairCorr(double r, double gr){
    pairCorr << r << "\t" << gr << "\n";
}

void Print::print_velDist(double v, double prob){
    velDist << v << "\t" << prob << "\n";
}

void Print::print_autoCorr(int t, double vaf){
    autoCorr << t << "\t" << vaf << endl;
}

void Print::print_MSD(int t, double msd){
    MSD << t << "\t" << msd << "\t" << log(msd) << "\t" << log(1.0-msd) << endl;
}

void Print::print_fluct(double r, double avg, double f){
    fluct <<  avg << "\t" << f << endl;
}

void Print::print_dens(double density, double count){
    dens << density << "\t" << count << endl;
}

void Print::print_summary(string ID, int noCells, double L, long int numberOfSteps, int stepsPerTime, double C1, double C2, double rho, double seconds, long int resetCounter, double corr,
                          double binder, double order, double variance) {
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
    summary << "Nearest neighbor correlation: " << "\t" << corr << endl;
    summary << "Binder cumulant:            " << "\t" << binder << endl;
    summary << "Average order parameter:    " << "\t" << order << endl;
    summary << "Order parameter variance:   " << "\t" << variance << endl;
    summary2 << ID << endl;
    summary2 << noCells << endl;
    summary2 << L << endl;
    summary2 << numberOfSteps << endl;
    summary2 << stepsPerTime << endl;
    summary2 << C1 << endl;
    summary2 << C2 << endl;
    summary2 << rho << endl;
    summary2 << seconds << endl;
    summary2 << resetCounter << endl;
    summary2 << corr << endl;
    summary2 << binder << endl;
    summary2 << order << endl;
    summary2 << variance << endl;
}

