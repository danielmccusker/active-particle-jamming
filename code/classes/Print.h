// Class containing all printing functions

// Uses Unix commands mkdir and cp

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sys/types.h>
#include <dirent.h>

using namespace std;

const string location = "/Users/Daniel1/Desktop/ActiveMatterResearch/jamming-dynamics/";
//const string location = "/home/dmccusker/remote/jamming-dynamics/";

struct Print
{
    Print();
    ~Print();
    
    string init(string, string, int);
    void print_data(long int, double, double, double, double, double, double, double, double, double, double);
    void print_vid(double, double, double, int, double, int);
    void print_Ovito(int, double, double, double, double, double, double, int);
    void print_noCells(int);
    void print_size(double);
    void print_v(double, double);
    void print_g(double, double);
    void print_MSD(int, double);
    void print_VAF(int, double);
    void print_GNF(double, double);
    void resize(double &, double &);
    void print_summary(string, int, long int, int, double, double, double, double);
    
    ofstream data, size, JavaScriptVid, OvitoVid, pairCorr, velDist, MSD, GNF, VAF, summary;

    string run;
    string path;
    int N, sqrtN;
};

Print::Print(){
    // This space is intentionally left blank
}

Print::~Print(){
	JavaScriptVid << "]\n];";

    data.close();
    JavaScriptVid.close();
    OvitoVid.close();
    pairCorr.close();
    velDist.close();
    summary.close();
    GNF.close();
    MSD.close();
    VAF.close();
}

string Print::init(string fullRun, string ID, int noCells){
// Create file directories, open output files
    
    run = ID;
    path = location+"output/"+fullRun+"/";
    N = noCells;
    sqrtN = sqrt(N);
    
    int check = 0;
    
    const char *p = path.c_str();
    DIR* dir = opendir(p);
    if(!dir){
        check += system(("mkdir "+path).c_str());
        check += system(("cp "+location+"code/postprocessing/{GNFvis.gnu,MSDvis.gnu,p*.py} "+path).c_str());
        cout << "This message should only appear once" << endl;
        if(check < 0){
            cout << "Failed to create folder for the full run. Status 715\n";
            exit(715);
        }
    } else if(dir){
        cout << "this message should appear every other time" << endl;
        closedir(dir);
    }

    check += system(("mkdir "+path+run).c_str());
    check += system(("mkdir "+path+run+"/dat").c_str());
    check += system(("mkdir "+path+run+"/eps").c_str());
    check += system(("mkdir "+path+run+"/vid").c_str());
    check += system(("cp "+location+"code/postprocessing/jam.html "+path+run+"/vid/jam.html").c_str());
    check += system(("cp "+location+"code/postprocessing/jquery.js "+path+run+"/vid/jquery.js").c_str());
    
    if(check < 0){
        cout << "Something hasn't been copied. Status 716\n";
        exit(716);
    }
    
    data.open((path+run+"/dat/data.dat").c_str());
    size.open((path+run+"/vid/size.js").c_str());
    JavaScriptVid.open((path+run+"/vid/data.js").c_str());
    OvitoVid.open((path+run+"/vid/ovito.txt").c_str());
    GNF.open((path+run+"/dat/GNF.dat").c_str());
    MSD.open((path+run+"/dat/MSD.dat").c_str());
    pairCorr.open((path+run+"/dat/pairCorr.dat").c_str());
    velDist.open((path+run+"/dat/velDist.dat").c_str());
    VAF.open((path+run+"/dat/vaf.dat").c_str());
    summary.open((path+run+"/dat/summary.dat").c_str());
    
    JavaScriptVid << "var cells = [\n";

    return path;
}

void Print::print_size(double R){
// The particle sizes don't change over time. They are immediately written after which the file can be closed.
// The last line is the size of the box

    static int p2=0;
    
    if(p2==0)   size << "var size = [\n{'r': " << 20*R << "}";
    else        size << ",\n{'r': " << 20*R << "}";
    
    if(p2==N){
        size << "\n];";
        size.close();
    }
    ++p2;
}

void Print::print_data(long int t, double xa, double ya, double phi, double Psi,
                         double x1, double y1, double x2, double y2, double x3, double y3){
    data << setprecision(8) << t << "\t" << xa << "\t" << ya << "\t" << "\t" << phi << "\t" << Psi << "\t"
        << x1 << "\t" << y1 << "\t" << x2 << "\t" << y2 << "\t" << x3 << "\t" << y3 << "\n";
}

void Print::print_vid(double x, double y, double a, int v, double d, int o){
// Prints the information about the cells to the JavaScript video file
// When called, it adds the intrinsics of one cell being the
// x position, y position, direction, velocity (magnitude and direction), overlap

	static int p=0;			// p: cell counter	    
	static bool first = true;

    ++p;

    if(p==1 && first){
        JavaScriptVid << "[";
        first = false;
	}
    else if(p==1)
        JavaScriptVid << "],\n[";

    resize(x, y);
    int ix = x;
    int iy = y;

    JavaScriptVid << "{'x': " << ix << ", 'y': " << iy << ", 'a': " << -a << ", 'v': " << v << ", 'd': " << -d << ", 'o': " << o << "}";

    if(p==N){
        JavaScriptVid << "\n";
        p=0;
    }
    else
        JavaScriptVid << ",\n";
}

void Print::print_Ovito(int cellIndex, double radius, double x, double y, double phi, double speed, double direction, int overlap){
    OvitoVid << cellIndex << "\t" << x << "\t" << y << "\t" << radius << "\t" << cos(phi) << "\t" << sin(phi) << "\t" << speed << "\t" << speed*cos(direction) << "\t" << speed*sin(direction) << "\t" << overlap << endl;
}

void Print::print_noCells(int noCells){
    OvitoVid << noCells << endl;
    OvitoVid << "time step comment" << endl;
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

void Print::print_GNF(double avg, double fluct){
    GNF << avg << "\t" << fluct << "\t" << log(avg) << "\t" << log(fluct) << endl;//use logs for gnuplot
}

void Print::resize(double &x, double &y){
// Resizes and translates the cells for display in JavaScript video

    x += sqrtN+30;
    y += sqrtN-70;
    x *= 20;
    y *= -20;
}

void Print::print_summary(string ID, int noCells, long int numberOfSteps, int stepsPerTime, double C1, double C2, double rho, double seconds) {
    summary << "Run ID:                     " << "\t" << ID << endl;
    summary << "Number of cells:            " << "\t" << noCells << endl;
    summary << "Number of steps:            " << "\t" << numberOfSteps << endl;
    summary << "Steps per unit time:        " << "\t" << stepsPerTime << endl;
    summary << "lambda_s (self-propulsion): " << "\t" << C1 << endl;
    summary << "lambda_n (noise):           " << "\t" << C2 << endl;
    summary << "Rho (density):              " << "\t" << rho << endl;
    summary << "Simulation time (seconds):  " << "\t" << seconds << endl;
}
