// ============================
// Fixed:
// - 
// Added:
// - 02: velocity autocorrelation function in MSD print
// - 03: Quantify crystalization by comparing radii of neighbours
// - 04: print the input summaryeters (Dan)
// Changed:
// - 03: The order of output in print_data
// Removed:
// - 

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sys/types.h>
#include <dirent.h>
using namespace std;

const string location = "/Users/Daniel1/Desktop/ActiveMatterResearch/jamming-dynamics/";

struct Print
{
    Print();
    ~Print();
    
    string init(string, string, int);
    void print_data(long int, double, double, double, double, double, double, double, double, double, double, double, double);
    void print_positions(double x, double y);
    void print_vid(double, double, double, int, double, int);
    void print_Ovito(int, double, double, double, double, double, double, int);
    void print_noCells(int);
    void print_size(double);
    void print_v(int, double, double);
    void print_g(int, double, double);
    //void print_MSD(int, double, double, double);
    void print_MSD(int, double);
    void resize(double &, double &);
    void print_summary(string, int, long int, int, double, double, double);
    void print_GNF(double, double);
    
    ofstream data, positions, JavaScriptVid, OvitoVid, size, pairCorr, velDist, MSD, summary, GNF;

    string run;
    string path;
    int N, sqrtN;
};

Print::Print()
{
    // This space is intentionally left blank
}

Print::~Print()
{
	JavaScriptVid << "]\n];";

    data.close();
    positions.close();
    JavaScriptVid.close();
    OvitoVid.close();
    pairCorr.close();
    velDist.close();
    summary.close();
    GNF.close();
    MSD.close();
}

string Print::init(string fullRun, string ID, int noCells)
{
    run = ID;
    path = location+"output/"+fullRun+"/";
    const char *p = path.c_str();
    N = noCells;
    sqrtN = sqrt(N);
    int check = 0;
    DIR* dir = opendir(p); 
	cout << "the value of dir is " << dir << endl;
    
    if(!dir){
        check += system(("mkdir "+path).c_str());
        check += system(("mkdir "+path+run).c_str());
        cout << "This message should only appear once" << endl;
        if(check < 0){
            cout << "Failed to create folders. Status 715\n";
            exit(715);
        }
    }
    if(dir){
        cout << "this message should appear every other time" << endl;
        closedir(dir);
    }
    
    check += system(("cp "+location+"code/postprocessing/* "+path).c_str());
    if( check < 0 ){
        cout << "The postprocessing code was not copied to the output folder.";
    }
    check += system(("mkdir "+path+run+"/dat").c_str());
    check += system(("mkdir "+path+run+"/eps").c_str());
    check += system(("mkdir "+path+run+"/vid").c_str());
    if(check < 0){
        cout << "Failed to create folders. Status 715\n"; exit(715);
    }
    check += system(("cp "+location+"output/jam.html "+path+run+"/vid/jam.html").c_str());
        
    if(check < 0){
        cout << "Failed to copy files to vid/. Status 716\n"; exit(716);
    }
    data.open((path+run+"/dat/data.dat").c_str());
    positions.open((path+run+"/dat/positions.dat").c_str());
    JavaScriptVid.open((path+run+"/vid/data.js").c_str());
    OvitoVid.open((path+run+"/vid/ovito.txt").c_str());
    size.open((path+run+"/vid/size.js").c_str());
    summary.open((path+run+"/dat/summary.dat").c_str());
    GNF.open((path+run+"/dat/GNF.dat").c_str());
    MSD.open((path+run+"/dat/MSD.dat").c_str());

    if(data.fail() || JavaScriptVid.fail() || size.fail())
    {
        cout << "failed opening file(s)\n"; exit(717);
    }
    JavaScriptVid << "var cells = [\n";
    return path;
}

void Print::print_data(long int timeCounter, double xa, double ya, double vgem, double phi, double Psi,
                         double da2, double x1, double y1, double x2, double y2, double x3, double y3)
{
    data << setprecision(8) << timeCounter << "\t" << xa << "\t" << ya << "\t" << vgem << "\t" << phi << "\t" << Psi  << "\t"
         << da2 << "\t" << x1 << "\t" << y1 << "\t" << x2 << "\t" << y2 << "\t" << x3 << "\t" << y3 << "\n";
}

void Print::print_positions(double x, double y)
{
	positions << x << "\t" << y << endl;
}

void Print::print_vid(double x, double y, double a, int v, double d, int o)
// Prints the information about the cells to the data file
// When called, it adds the intrinsics of one cell being the
// x position, y position, direction, velocity (magnitude and direction), overlap
{
	static int p=0;			// p: cell counter	    
	static bool first = true;

    ++p;

    if(p==1 && first)
    {
        JavaScriptVid << "[";
        first = false;
	}
    else if(p==1)
        JavaScriptVid << "],\n[";

    resize(x, y);
    int ix = x;
    int iy = y;

    JavaScriptVid << "{'x': " << ix << ", 'y': " << iy << ", 'a': " << -a << ", 'v': " << v << ", 'd': " << -d << ", 'o': " << o << "}";

    if(p==N)
    {
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

void Print::print_size(double R)
// The particle sizes don't change over time. They are immediately written after which the file can be closed.
{
    static int p2=0;

	if(p2==0)
		size << "var size = [\n{'r': " << 20*R << "}";
	else
		size << ",\n{'r': " << 20*R << "}";

	if(p2==N)   // The last one is the size of the box
	{
		size << "\n];";
		size.close();
	}
	++p2;
}

void Print::print_v(int k, double v, double prob)
{
    if( k == 0 )
    {
        velDist.open((path+run+"/dat/velDist.dat").c_str());
    }
    velDist << v << "\t" << prob << "\n";
}

void Print::print_g(int k, double r, double gr)
{
    if( k == 1 )
    {
        pairCorr.open((path+run+"/dat/pairCorr.dat").c_str());
    }
    pairCorr << r << "\t" << gr << "\n";
}

//void Print::print_MSD(int tau, double msd, double err, double vaf)
//{
//    if( tau == 0 )
//    {
//        MSD.open((path+run+"/dat/MSD.dat").c_str());
//    }
//    MSD << tau << "\t" << msd << "\t" << log(msd) << "t" << log(1.-msd) << "\t" << err << "\t" << vaf << "\n";
//}

void Print::print_MSD(int t, double msd){
    MSD << t << "\t" << msd << "\t" << log(msd) << "\t" << log(1.0-msd) << endl;
}

void Print::resize(double &x, double &y)
// Resizes (and translates) the cells such that they are beautifully displayed
{
    x += sqrtN+30;
    y += sqrtN-70;
    x *= 20;
    y *= -20;
}

void Print::print_summary(string ID, int noCells, long int numberOfSteps, int stepsPerTime, double C1, double C2, double rho)
//outputs the summary
{
    summary << "Run ID:                     " << "\t" << ID << endl;
    summary << "Number of Cells:            " << "\t" << noCells << endl;
    summary << "Number of Steps:            " << "\t" << numberOfSteps << endl;
    summary << "Steps per unit time:        " << "\t" << stepsPerTime << endl;
    summary << "lambda_s (self-propulsion): " << "\t" << C1 << endl;
    summary << "lambda_n (noise):           " << "\t" << C2 << endl;
    summary << "rho (density):              " << "\t" << rho << endl;
    //summary << "simulation ran for " << seconds << " seconds." << endl;
}

void Print::print_GNF(double avg, double fluct)
{
    GNF << avg << "\t" << fluct << "\t" << log(avg) << "\t" << log(fluct) << endl;//use logs for gnuplot
}


