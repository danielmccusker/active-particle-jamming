// 1024 circular particles are put together in a box with double periodic boundary conditions. We will test for what values of the parameter
// setting the self-propulsion force, alignment torque and density the system will be jammed. We say a system is jammed when there are no
// internal rearrangements, that is, the set of neighbours does not change. Hence, a system can be jammed with a net movement. The output
// of the simulation is a file with: time, order, net orientation, average velocity, average position, position of three particles, and, in
// three different files: the velocity distribution, the pair correlation function, and the mean square displacement.
// In this file, the alignment is done by looking at the neighbours, rather than aligning with the velocity.
// This test file should work for all simulation length. The last 256 frames are recorded

// ============================
// Fixed:
// - 01: Relax the system as a passive particle system: CTalign = 0, -PI < cell[i].psi <= PI
// - 05: Lagtimes are divided by screenshotInterval
// - 05: Set the original position (x0,y0) using the real position instead of position in box
// - 05: Choose a 5 times smaller dr for the paircorrelation function
// - 06: MSD was not calculated at timeCounter = tau1, tau2, tau3, tau4
// - 06: atoi and atol are used correctly
// Added:
// - 01: GNU-plotting after simulation
// - 03: Mean square displacement
// - 06: Size simulation box output to config_2.dat
// - 06: Calculation of velocity auto correlation function
// - 07: Quantify crystalization by comparing radii of neighbours
// - 09: Average distance to neighbours to detect holes
// - 11: Output position, velocity, size and orientation 10 times during the run

// * Dan * //
// - 12: added print file for input parameters, improved function for pair correlation function
// - 13: added timer, circle overlap for calculating GNF

// Changed:
// - 01: Choose a ten times smaller value for dr
// - 01: screenshotInterval 512 -> 256
// - 02: Alignment and noise are combined into a single update mechanism [ psi_i = avg(psi_j) + lambda_n * PI * randuni() ]
// - 02: The binwidth dr depends on the displacement per step
// - 02: Uniform random distribution interval to (-PI, PI] instead of (-1, 1]
// - 02: "../cell/jam_01.h"
// - 02: Back to polydisperse particles
// - 02: screenshotInterval 256 -> 128
// - 03: "../cell/jam_02.h"
// - 03: "../print/jam_01.h"
// - 04: newSkinList incorporates movement of average position and saves the number of refreshes in config_2
// - 06: "../cell/jam_03.h"
// - 06: "../print/jam_02.h"
// - 07: Rearranged the output to data.dat
// - 07: "../cell/jam_04.h"
// - 07: "../print/jam_03.h"
// - 07: More output by gnuplot
// - 08: Cnoise depends on the number of steps per unit time (Undone in 09)
// - 09: Monodisperse particles that move at the same speed
// - 10: Polydisperse particles that move at the same speed
// - 11: More particles 1024 -> 4096

// * Dan * //
// - 12: "phi" to "rho" (density), ../cell/jam_04.h -> ../cell/cell_04.h, ../print/jam_03.h -> "../print/print_04.h"

// Removed:
// - 02: CTalign
// - 02: theta_norm
// - 02: uniform in distribution
// - 04: Removed the physics function and moved the contents all to find_nbr()
// - 07: Removed the nbr vector
// - 09: Removed adif2

// * Dan * //
// - 12: removed old correlation function


#include <vector>
#include <cmath>
#include <ctime>
#include <chrono>

#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include "../cell/cell_04.h"
#include "../print/print_04.h"

//Mersenne Twister pseudorandomnumber generator
boost::mt19937 gen(time(0));
boost::uniform_real<> unidist(-PI, PI);
boost::normal_distribution<> normdist(0, 1);
boost::variate_generator< boost::mt19937, boost::uniform_real<> > randuni(gen, unidist);
boost::variate_generator< boost::mt19937, boost::normal_distribution<> > randnorm(gen, normdist);

const int screenshotInterval = 128;

using namespace std;
using namespace std::chrono;

struct Engine
{
    Engine(string, string, int, long int, int, double, double, double);
    ~Engine();

    // Functions
    void start();
    void init();
    void relax();
    bool newSkinList();
    void find_nbrRegion();
    void find_nbr();
    double delta_norm(double);
    double cc_overlap(double R1, double R2, double r);
    double totalAreaIntersection(double R);
    void pairCorr();
    void graphics();
    void GNFinit();
    void densityFluctuations();
                                                // [++|--|--|--|--|--|--|--|--|--|-]
    void setMSD_and_update0();                  // [--+--|--|--|--|--|--|--|--|--|-]
    void setMSD_and_update();                   // [--|--+--|--|--|--|--|--|--|--|-]
    void setMSD_and_update(int);                // [--|--|--+--|--|--|--|--|--|--|-]
    void setMSD_and_update(int, int);           // [--|--|--|--+--|--|--|--|--|--|-]
    void setMSD_and_update(int, int, int);      // [--|--|--|--|--+--|--|--|--|--|-]
    void calcMSD_and_update1();                 // [--|++|--|--|--|--|--|--|--|--|-]
    void calcMSD_and_update2();                 // [--|--|++|--|--|--|--|--|--|--|-]
    void calcMSD_and_update3(int);              // [--|--|--|++|--|--|--|--|--|--|-]
    void calcMSD_and_update4(int, int);         // [--|--|--|--|++|--|--|--|--|--|-]
    void calcMSD_and_update5(int, int, int);    // [--|--|--|--|--|+++++++++++++++-]
    void calcMSD_and_update6(int, int, int);    // [--|--|--|--|--|--|--|--|--|--|+]
    void save_config();
    void save_config2();
    // Variables
    Cell cell[noCells]; //"cell" is array of "Cell" objects
    Print print;

    string run;
    string fullRun;
    string path;
    int nbrRegionRefresh;
    int film;

    long int totalSteps;
    long int countdown;
    int noCells;
    long int timeCounter;
    long int ts0, ts1, ts2, ts3, ts4;       // Times of the MSDs starting points
    long int tau0, tau1, tau2, tau3, tau4;  // Lagtimes in units of screenshotInterval starting at ts0, ts1, ...
    long int resetCounter;                  // Records the number of times the Skinlist is determined
    double deltat;
    double CFself;
    double CTnoise;
    double dens;
    double nbrRegionRadius2;
    double nbrRadius2;
    double listRefreshDistance2;
    double l;                               // Size of the box
    double lover2;
    double xavg;                            // Average position of all the particles at any time, based on the
    double yavg;                            // real positions as if there were no periodic boundary conditions
    double xavgold;                         // Old value of the average position of all particles used for 
    double yavgold;                         // newSkinlist
    double vcx, vcy;                        // Velocity of the center of mass
    double vg;                              // Average speed of the particles (not the speed of center of mass!)
    double ux, uy;                          // Sum of the x and y-components of the orientations of the particles
    double nbrdist2;                        // Average squared distance to nbr (small for crystalized with holes?)

    double x0[5*noCells];                   // Records the x and y positions where we start
    double y0[5*noCells];                   // measuring the MSD at five different points
    double vx0[5*noCells];                  // Records the x and y velocity for the velocity autocorrelation
    double vy0[5*noCells];                  // function at five different points in time
    double xnavg[5];                        // Records the average x and y position of all particle at the
    double ynavg[5];                        // time we start measuring the MSDs at five differents points
    vector<double> vaf;
    vector<double> MSD;
    vector<double> MSDerr;
    vector<int> MSDcounter;
    
    int rmsCounter;                         //counts time steps for calculating rms density fluctuations
    double rmsValue;
    double R;                               //radius of GNF intersection circle
    double numberOfGNFpoints;               //how many divisions of total run time for each density fluctuation measurement
    double GNFinterval;
};

Engine::Engine(string dir, string ID, int N, long int numberOfSteps, int stepsPerTime, double C1, double C2, double rho)
// Constructor: This is where variables get defined and initialized.
{
    fullRun = dir;
    run = ID;
    noCells = N;
    totalSteps = numberOfSteps;
    countdown = numberOfSteps+1;
    deltat    = 1./stepsPerTime;
    CFself    = C1;
    CTnoise   = C2;
    dens      = rho;

    nbrRegionRefresh     = 16;                            // Tune this value, preferably a power of two
    film                 = 256*screenshotInterval;
    timeCounter          = 0;
    ts0                  = 8192;
    while( 2 * ts0 < countdown / 5 ) { ts0 *= 2;}
    ts1                  = ts0 + ts0;
    ts2                  = ts1 + ts0;
    ts3                  = ts2 + ts0;
    ts4                  = ts3 + ts0;
    tau0 = -ts0/screenshotInterval;
    tau1 = -ts1/screenshotInterval;
    tau2 = -ts2/screenshotInterval;
    tau3 = -ts3/screenshotInterval;
    tau4 = -ts4/screenshotInterval;
    resetCounter         = 0;
    nbrRegionRadius2     = 20;                            // Tune this value
    nbrRadius2           = 2.8*2.8;
    listRefreshDistance2 = 0.25*(nbrRegionRadius2 + nbrRadius2 - 2*sqrt(nbrRegionRadius2*nbrRadius2));

    long int vecLength = (countdown - ts0)/screenshotInterval+1;
    vaf.reserve(vecLength);
    MSD.reserve(vecLength);
    MSDerr.reserve(vecLength);
    MSDcounter.reserve(vecLength);
    for(long int t=0; t<vecLength; ++t)
    {
        vaf.push_back(0);
        MSD.push_back(0);
        MSDerr.push_back(0);
        MSDcounter.push_back(0);
    }

}

Engine::~Engine()
{
    const int gmax = 60001;
    const int vmax = 150;
    int g[gmax] = {0};
    int veldist[vmax] = {0};
    double dr = 0.0001;
    double dv = CFself/100;
    
    for( int i=0; i<noCells; ++i)
    {
        double v = sqrt(cell[i].velx*cell[i].velx + cell[i].vely*cell[i].vely);
        int binv = v/dv;
        if(binv >= vmax)
        {
            cout << binv << "=149\n";
            binv = 149;
        }
        veldist[binv] += 1;
        
    }
    
    for( int k=0; k<vmax; ++k) print.print_v(k, k*dv, veldist[k] / (double)noCells);
    
    for(long int t=0; t<(timeCounter - ts0)/screenshotInterval+1; ++t)
    {
        long int norm = noCells*MSDcounter[t];
        double msdave = MSD[t]/norm;
        print.print_MSD(t*screenshotInterval, msdave, sqrt(MSDerr[t]/norm - msdave*msdave), (vaf[t]/norm)/2);
    }
    
    ofstream overWrite;
    overWrite.open((path+"/"+run+"/dat/config_2.dat").c_str());
    overWrite << resetCounter;
    overWrite.close();
}

void Engine::start()
// The system of particles in initialized. Also contains the main loop.
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now(); //start time
    
    path = print.init(fullRun, run, noCells);
    cout << path << endl;
    
    init();
    GNFinit();
    relax();

    while(countdown != 0)    // time loop
    {
        if(!(timeCounter%nbrRegionRefresh))
        {
            if(newSkinList())
                find_nbrRegion();
        }

        find_nbr();     // find nbrs

        graphics();     // output data
        
        densityFluctuations();

        ++timeCounter;
        --countdown;
        //cout << countdown << endl;
    }
    
    pairCorr();
    
    high_resolution_clock::time_point t2 = high_resolution_clock::now(); //end time
    
    auto duration = duration_cast<seconds>( t2 - t1 ).count();//elapsed time
    
    print.print_summary(run, noCells, timeCounter, 1./deltat, CFself, CTnoise, dens, duration); //print run summary
}

void Engine::init()
// Particle positions and orientations are initialized on a square lattice, initialize GNF measurements
{
    // Assign particle sizes and calculate boxsize
    double area = 0;
    for(int i=0; i<noCells; ++i)
    {
        double cellRad = 1. + randnorm()/10;	    // Use polydisperse particles
        double R2 = cellRad*cellRad;
        area += R2;
        cell[i].R  = cellRad;
        print.print_size(cellRad);
    }
    l = sqrt(PI * area / dens);
    lover2 = l/2.;
    print.print_size(l);

    // Place them in a square with the appropriate size L
    int sqrtnoCells = sqrt(noCells);
    double spacing = l / sqrtnoCells;
	xavg = 0;
	yavg = 0;
	xavgold = 0;
	yavgold = 0;

    for(int i=0; i<noCells; ++i)
    {
        cell[i].posx = -lover2 + spacing*(i%sqrtnoCells) + randnorm()/3.;
        cell[i].posy = -lover2 + spacing*(i/sqrtnoCells) + randnorm()/3.;
        cell[i].psi  = randuni();
        cell[i].cosp = cos(cell[i].psi);
        cell[i].sinp = sin(cell[i].psi);
        cell[i].cosp_new = cell[i].cosp;        // Average direction of particles in nbrhood also includes itself
        cell[i].sinp_new = cell[i].sinp;
    }
    
   
}

void Engine::relax()
// This function relaxes the system for 1 000 000 steps, slowly decreasing the activity to the final value
{
    int trelax = 1000000;
    double CFrelax = 0.05;
    double CFself_old = CFself;
    double CTnoise_old = CTnoise;
    timeCounter = 1;                // Prevent the programme from producing output
    CTnoise = 0;

    for( int t=0; t < trelax; ++t)
    {
        CFself = CFself_old + ((CFrelax - CFself_old)*(trelax - t))/trelax;

        if(!(t%nbrRegionRefresh))
        {
            if(newSkinList())
                find_nbrRegion();
        }

        find_nbr();     // find nbrs

        for( int i = 0; i<noCells; ++i)
        {
			cell[i].update(deltat, lover2, l);
            cell[i].psi = randuni();
        }
    }

    for(int i=0; i<noCells; ++i)
    {
        cell[i].xreal = cell[i].posx;
        cell[i].yreal = cell[i].posy;
    }
    timeCounter = 0;
    resetCounter = 0;
    CFself = CFself_old;
    CTnoise = CTnoise_old;
}

bool Engine::newSkinList()
// Checks whether particles have changed their positions far enough to force a skin list refresh
// Returns true if list needs to be refreshed
{
    bool refresh= false;
    for(int i=0; i<noCells; ++i)
    {
        double dx = delta_norm(cell[i].xreal - xavg - cell[i].xold + xavgold);
        double dy = delta_norm(cell[i].yreal - yavg - cell[i].yold + yavgold);
        if(dx*dx+dy*dy > listRefreshDistance2)
        {
            xavgold = 0;
            yavgold = 0;
            refresh = true;
            ++resetCounter;
            for(i=0; i<noCells; ++i)    // use same index so break is not needed.
            {
                cell[i].xold = cell[i].xreal;
                cell[i].yold = cell[i].yreal;
                xavgold += cell[i].xreal;
                yavgold += cell[i].yreal;
                cell[i].nbrRegion.clear();  // Empty the nbrRegion vectors before adding new indices
            }
            xavgold /= noCells;
            yavgold /= noCells;
        }
    }

    return refresh;
}

void Engine::find_nbrRegion()
// Fills nbrRegion with the indeces of the cells that are within radius 
// nbrRegionRadius of cell i. 
{
	for(int i=0; i < noCells; ++i)
	{
		for(int j=i+1; j < noCells; ++j)                  // Use symmetry
		{
			double deltax = delta_norm(cell[j].posx-cell[i].posx);
			double deltay = delta_norm(cell[j].posy-cell[i].posy);
			double r2 = deltax*deltax+deltay*deltay;

			if(r2 < nbrRegionRadius2)   // if j is in i's circle
			{
				cell[i].nbrRegion.push_back(j);
				cell[j].nbrRegion.push_back(i);
			}
        }
    }
}

void Engine::find_nbr()
// Determines the neighbours from the Verlet list
// Since we calculate the distance between neighbours anyway, we already calculate the forces due to overlap
// We use a simple spring-like potential, that is, the force is proportional to the linear overlap.
{
    if(timeCounter%screenshotInterval != 0)
    {
        for(int i=0; i < noCells; ++i)
        {
            int max = cell[i].nbrRegion.size();
            for(int k=0; k < max; ++k)
            {
                int j = cell[i].nbrRegion[k];
                if(j > i)       // Use symmetry; only use indexes j > i
                {
                    double deltax = delta_norm(cell[j].posx-cell[i].posx);
                    double deltay = delta_norm(cell[j].posy-cell[i].posy);
                    double r2 = deltax*deltax+deltay*deltay;

                    if(r2 < nbrRadius2)    // They are neighbours
                    {
                        double sumR = cell[i].R + cell[j].R;

                        if(r2 < (sumR*sumR))     // They overlap
                        {
                            double overlap = sumR / sqrt(r2) - 1;
                            cell[i].Fx -= overlap*deltax;
                            cell[i].Fy -= overlap*deltay;
                            cell[j].Fx += overlap*deltax;
                            cell[j].Fy += overlap*deltay;
                        }
                        cell[i].cosp_new += cell[j].cosp;
                        cell[i].sinp_new += cell[j].sinp;
                        cell[j].cosp_new += cell[i].cosp;
                        cell[j].sinp_new += cell[i].sinp;
                    }
                }
            }
            double FselfR = CFself*cell[i].R;
            cell[i].Fx += cell[i].cosp*FselfR;
            cell[i].Fy += cell[i].sinp*FselfR;
            cell[i].psi_new = atan2(cell[i].sinp_new, cell[i].cosp_new) + CTnoise*randuni();
        }
    }
    else
    {
        xavg = 0;
        yavg = 0;
        vcx = 0;
        vcy = 0;
        vg = 0;
        ux = 0;
        uy = 0;
        int nbrNormCounter = 0;

        for(int i=0; i < noCells; ++i)
        {
            int max = cell[i].nbrRegion.size();
            for(int k=0; k < max; ++k)
            {
                int j = cell[i].nbrRegion[k];
                if(j > i)       // Use symmetry; only use indexes j > i
                {
                    double deltax = delta_norm(cell[j].posx-cell[i].posx);
                    double deltay = delta_norm(cell[j].posy-cell[i].posy);
                    double r2 = deltax*deltax+deltay*deltay;

                    if(r2 < nbrRadius2)    // They are neighbours
                    {
                        double sumR = cell[i].R + cell[j].R;

                        nbrdist2 += r2;
                        nbrNormCounter += 1;

                        if(r2 < (sumR*sumR))     // They overlap
                        {
                            double overlap = sumR / sqrt(r2) - 1;
                            cell[i].Fx -= overlap*deltax;
                            cell[i].Fy -= overlap*deltay;
                            cell[j].Fx += overlap*deltax;
                            cell[j].Fy += overlap*deltay;
                            if(countdown <= film && timeCounter%screenshotInterval == 0)
                            {
                                cell[i].over -= 240*overlap;
                            }
                        }
                        cell[i].cosp_new += cell[j].cosp;
                        cell[i].sinp_new += cell[j].sinp;
                        cell[j].cosp_new += cell[i].cosp;
                        cell[j].sinp_new += cell[i].sinp;
                    }
                }
            }
            cell[i].Fx += cell[i].cosp*CFself;
            cell[i].Fy += cell[i].sinp*CFself;
            cell[i].psi_new = atan2(cell[i].sinp_new, cell[i].cosp_new) + CTnoise*randuni();
            xavg += cell[i].xreal;
            yavg += cell[i].yreal;
            vcx += cell[i].velx;
            vcy += cell[i].vely;
            vg += sqrt(cell[i].velx*cell[i].velx + cell[i].vely*cell[i].vely);
            ux += cell[i].cosp;
            uy += cell[i].sinp;
            
        }
        xavg /= noCells;
        yavg /= noCells;
        vcx /= noCells;
        vcy /= noCells;
        vg /= noCells;
        nbrdist2 /= nbrNormCounter;
    }
}

double Engine::delta_norm(double delta)
// Subtracts multiples of the boxsize until -L/2 < delta <= L/2
{
    int k=-1;
    if(delta < -lover2) k=1;
    while(delta < -lover2 || delta >= lover2) delta += k*l;
    
    return delta;
}

void Engine::pairCorr()
//calculates the pair correlation function at a certain time step.
{
    double g = 0.0; //value of the correlation function
    double r = 0; //current radius
    double dr = l/500;  //step size is 1/500 of box size
    double C = l*l/(2*PI); //geometric normalization factor that depends on the system dimension
    
    double drij = 0.0; //distance between a pair
    double dx = 0.0;
    double dy = 0.0;
    double sum = 0.0;
    
    int k = 1;
    
    while(r < lover2){
        sum = 0.0;
        for (int i=0; i < noCells; i++){ //iterate through pairs
            for(int j=i+1; j < noCells; j++){
                dx = delta_norm(cell[i].posx-cell[j].posx);
                dy = delta_norm(cell[i].posy-cell[j].posy);
                drij = sqrt(dx*dx + dy*dy);
                
                if(abs(drij-r) < dr)
                    sum = sum + 1.0; //add if a pair separation is in the interval [r, r+dr]
            }
        }
        sum = sum / noCells;
        g = (C/r)*sum; //surface area of shell, rad squared if in 3D
        print.print_g(k, r, g);
        r = r + dr;
        k++;
    }
}

void Engine::GNFinit()//just initialize GNF measurements
{
    numberOfGNFpoints = 10;
    R = lover2;
    GNFinterval = (double)totalSteps/(double)numberOfGNFpoints;
    rmsCounter = 0;
    rmsValue = 0;
}

double Engine::cc_overlap(double R1, double R2, double r)
// Two circles have radii R1 and R2 with center C1 and C2.
// Their centers are distance r < R1 + R2 apart such that we have intersections A and B.
// We calculate here the overlapping area A1 (closest to C1) and A2 (closest to C2).
// d(C1, AB) = (R1^2 - R2^2 + r^2)/(2r)    d(C2, AB) = (R2^2 - R1^2 + r^2)/(2r)
// d^2(A, C1C2) = R1^2 - d^2(C1, AB) = R2^2 - d^2(C2, AB) = d(B, C1C2)
// theta1 = arccos((R1^2 - R2^2 + r^2)/(2 R1 r))
// A1 = 2*theta1/(2 PI) * R1^2 * PI - d(C1, AB) * R1*sin(theta1)
//    = R1^2 * (theta1 - (R1^2 - R2^2 + r^2)/(2*r*R1) * sin(theta1))
{
    if ( R1 + R2 <= r )     // There is no overlap between C1 and C2
    {
        return 0;
    }
    if ( R2 >= R1 + r )     // Circle C1 lies completely in C2
    {
        return PI*R1*R1;
    }
    if ( R1 >= R2 + r )     // Circle C2 lies completely in C1
    {
        return PI*R2*R2;
    }
    
    double R12 = R1*R1;
    double R22 = R2*R2;
    double x = (R12 - R22 + r*r)/(2*r); // base of triangle at C1
    double theta = acos(x/R1);          // angle of triangle at C1
    double A = R12*theta - x*R1*sin(theta);
    
    x = r - x;                            // Other area (at C2-side)
    theta = acos(x/R2);
    A += R22*theta - x*R2*sin(theta);
    
    return A;
}

double Engine::totalAreaIntersection(double R)
{   //draw a circle of radius R at the center of the box, calculate total overlapping area with all cells
    double A = 0;
    double Atot = 0;
    double d = 0;
    
    for (int i=0; i<noCells; i++) {
        d = sqrt((cell[i].posx*cell[i].posx) + (cell[i].posy*cell[i].posy));
        A = cc_overlap(cell[i].R, R, d);
        Atot = Atot+A;
    }
    //cout << "total intersecting area: " << Atot << endl;
    return Atot;
    
}

void Engine::densityFluctuations(){
    
    rmsCounter++;
    double expectedA = dens*PI*R*R;

    if(rmsCounter < GNFinterval) { //add square of intersecting areas to calculate rms
        double A = totalAreaIntersection(R);
        //cout << "R: " << R << " total intersecting area: " << A << endl;
        rmsValue = rmsValue + (A-expectedA)*(A-expectedA);
        //cout << rmsValue << endl;
    } else { //store value of rms fluctuation and reset for next measurement, increase circle size
        print.print_GNF(dens*PI*R*R, sqrt(rmsValue/(double)rmsCounter) );
        //cout << "printing rms value: " << sqrt(rmsValue/(double)rmsCounter) << endl;
        rmsValue = 0;
        rmsCounter = 0;
        R = R/2;
    }
}

void Engine::graphics()
// Creates the output files and resets output parameters.
{
	if((timeCounter%screenshotInterval) != 0)
	{
        for(int i=0; i<noCells; i++)
        {
			cell[i].update(deltat, lover2, l);
		}
	}
	else
	{
        int n2 = 2*noCells, n3 = 3*noCells, n4 = 4*noCells;

        if( tau4 > 0 && countdown > film )      // [--|--|--|--|--|+++++++++++++++-]
            calcMSD_and_update5(n2, n3, n4);
        else if( tau0 < 0 )                     // [++|--|--|--|--|--|--|--|--|--|-]
        {
            for(int i=0; i<noCells; i++)
            {
                cell[i].update(deltat, lover2, l);
            }
        }
        else if( tau0 > 0 && tau1 < 0 )         // [--|++|--|--|--|--|--|--|--|--|-]
            calcMSD_and_update1();
        else if( tau1 > 0 && tau2 < 0 )         // [--|--|++|--|--|--|--|--|--|--|-]
            calcMSD_and_update2();
        else if( tau2 > 0 && tau3 < 0 )         // [--|--|--|++|--|--|--|--|--|--|-]
            calcMSD_and_update3(n2);
        else if( tau3 > 0 && tau4 < 0 )         // [--|--|--|--|++|--|--|--|--|--|-]
            calcMSD_and_update4(n2, n3);
        else if( tau4 > 0 && countdown <= film )// [--|--|--|--|--|--|--|--|--|--|+]
        {
            print.print_noCells(noCells);
            calcMSD_and_update6(n2, n3, n4);
        }
        else if( tau0 == 0 )                    // [--+--|--|--|--|--|--|--|--|--|-]
            setMSD_and_update0();
        else if( tau1 == 0 )                    // [--|--+--|--|--|--|--|--|--|--|-]
            setMSD_and_update();
        else if( tau2 == 0 )                    // [--|--|--+--|--|--|--|--|--|--|-]
            setMSD_and_update(n2);
        else if( tau3 == 0 )                    // [--|--|--|--+--|--|--|--|--|--|-]
            setMSD_and_update(n2, n3);
        else if( tau4 == 0 )                    // [--|--|--|--|--+--|--|--|--|--|-]
            setMSD_and_update(n2, n3, n4);
        else
            cout << "Something is wrong: " << timeCounter << endl;

		print.print_data(timeCounter, xavg, yavg, vg, sqrt(ux*ux+uy*uy)/noCells, atan2(uy, ux), nbrdist2,
                        cell[330].posx, cell[330].posy, cell[341].posx, cell[341].posy, cell[687].posx, cell[687].posy);
		if(timeCounter%65536==0)
		{
			save_config();
            if (timeCounter%ts0 == 0)
            {
                save_config2();
            }
		}

        ++tau0;
        ++tau1;
        ++tau2;
        ++tau3;
        ++tau4;
	}
}

void Engine::setMSD_and_update0()
// Sets the starting points for calculating the MSD and updates the positions
{
    for(int i=0; i<noCells; i++)
    {
        x0[i] = cell[i].xreal;
        y0[i] = cell[i].yreal;
        vx0[i] = cell[i].velx-vcx;
        vy0[i] = cell[i].vely-vcy;

        cell[i].update(deltat, lover2, l);
    }

    xnavg[0] = xavg;
    ynavg[0] = yavg;
}

void Engine::setMSD_and_update()
// Sets the starting points for calculating the MSD and updates the positions
{
    for(int i=0; i<noCells; i++)
    {
        x0[noCells+i] = cell[i].xreal;
        y0[noCells+i] = cell[i].yreal;
        vx0[noCells+i] = cell[i].velx-vcx;
        vy0[noCells+i] = cell[i].vely-vcy;

        double dxrel = (cell[i].xreal - x0[i] - (xavg - xnavg[0]));
        double dyrel = (cell[i].yreal - y0[i] - (yavg - ynavg[0]));
        double dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau0] += vx0[i]*(cell[i].velx - vcx) + vy0[i]*(cell[i].vely - vcy);
        MSD[tau0] += dr2;
        MSDerr[tau0] += dr2*dr2;

        cell[i].update(deltat, lover2, l);
    }
    MSDcounter[tau0] += 1;

    xnavg[1] = xavg;
    ynavg[1] = yavg;
}

void Engine::setMSD_and_update(int n2)
// Sets the starting points for calculating the MSD and updates the positions
{
    for(int i=0; i<noCells; i++)
    {
        x0[n2+i] = cell[i].xreal;
        y0[n2+i] = cell[i].yreal;
        vx0[n2+i] = cell[i].velx-vcx;
        vy0[n2+i] = cell[i].vely-vcy;

        double dxrel = (cell[i].xreal - x0[i] - (xavg - xnavg[0]));
        double dyrel = (cell[i].yreal - y0[i] - (yavg - ynavg[0]));
        double dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau0] += vx0[i]*(cell[i].velx - vcx) + vy0[i]*(cell[i].vely - vcy);
        MSD[tau0] += dr2;
        MSDerr[tau0] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[noCells+i] - (xavg - xnavg[1]));
        dyrel = (cell[i].yreal - y0[noCells+i] - (yavg - ynavg[1]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau1] += vx0[noCells+i]*(cell[i].velx - vcx) + vy0[noCells+i]*(cell[i].vely - vcy);
        MSD[tau1] += dr2;
        MSDerr[tau1] += dr2*dr2;
        cell[i].update(deltat, lover2, l);
    }
    MSDcounter[tau0] += 1;
    MSDcounter[tau1] += 1;


    xnavg[2] = xavg;
    ynavg[2] = yavg;
}

void Engine::setMSD_and_update(int n2, int n3)
// Sets the starting points for calculating the MSD and updates the positions
{
    for(int i=0; i<noCells; i++)
    {
        x0[n3+i] = cell[i].xreal;
        y0[n3+i] = cell[i].yreal;
        vx0[n3+i] = cell[i].velx-vcx;
        vy0[n3+i] = cell[i].vely-vcy;

        double dxrel = (cell[i].xreal - x0[i] - (xavg - xnavg[0]));
        double dyrel = (cell[i].yreal - y0[i] - (yavg - ynavg[0]));
        double dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau0] += vx0[i]*(cell[i].velx - vcx) + vy0[i]*(cell[i].vely - vcy);
        MSD[tau0] += dr2;
        MSDerr[tau0] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[noCells+i] - (xavg - xnavg[1]));
        dyrel = (cell[i].yreal - y0[noCells+i] - (yavg - ynavg[1]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau1] += vx0[noCells+i]*(cell[i].velx - vcx) + vy0[noCells+i]*(cell[i].vely - vcy);
        MSD[tau1] += dr2;
        MSDerr[tau1] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[n2+i] - (xavg - xnavg[2]));
        dyrel = (cell[i].yreal - y0[n2+i] - (yavg - ynavg[2]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau2] += vx0[n2+i]*(cell[i].velx - vcx) + vy0[n2+i]*(cell[i].vely - vcy);
        MSD[tau2] += dr2;
        MSDerr[tau2] += dr2*dr2;

        cell[i].update(deltat, lover2, l);
    }
    MSDcounter[tau0] += 1;
    MSDcounter[tau1] += 1;
    MSDcounter[tau2] += 1;

    xnavg[3] = xavg;
    ynavg[3] = yavg;
}

void Engine::setMSD_and_update(int n2, int n3, int n4)
// Sets the starting points for calculating the MSD and updates the positions
{
    for(int i=0; i<noCells; i++)
    {
        x0[n4+i] = cell[i].xreal;
        y0[n4+i] = cell[i].yreal;
        vx0[n4+i] = cell[i].velx-vcx;
        vy0[n4+i] = cell[i].vely-vcy;

        double dxrel = (cell[i].xreal - x0[i] - (xavg - xnavg[0]));
        double dyrel = (cell[i].yreal - y0[i] - (yavg - ynavg[0]));
        double dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau0] += vx0[i]*(cell[i].velx - vcx) + vy0[i]*(cell[i].vely - vcy);
        MSD[tau0] += dr2;
        MSDerr[tau0] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[noCells+i] - (xavg - xnavg[1]));
        dyrel = (cell[i].yreal - y0[noCells+i] - (yavg - ynavg[1]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau1] += vx0[noCells+i]*(cell[i].velx - vcx) + vy0[noCells+i]*(cell[i].vely - vcy);
        MSD[tau1] += dr2;
        MSDerr[tau1] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[n2+i] - (xavg - xnavg[2]));
        dyrel = (cell[i].yreal - y0[n2+i] - (yavg - ynavg[2]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau2] += vx0[n2+i]*(cell[i].velx - vcx) + vy0[n2+i]*(cell[i].vely - vcy);
        MSD[tau2] += dr2;
        MSDerr[tau2] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[n3+i] - (xavg - xnavg[3]));
        dyrel = (cell[i].yreal - y0[n3+i] - (yavg - ynavg[3]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau3] += vx0[n3+i]*(cell[i].velx - vcx) + vy0[n3+i]*(cell[i].vely - vcy);
        MSD[tau3] += dr2;
        MSDerr[tau3] += dr2*dr2;

        cell[i].update(deltat, lover2, l);
    }
    MSDcounter[tau0] += 1;
    MSDcounter[tau1] += 1;
    MSDcounter[tau2] += 1;
    MSDcounter[tau3] += 1;

    xnavg[4] = xavg;
    ynavg[4] = yavg;
}

void Engine::calcMSD_and_update1()
// [--|++|--|--|--|--|--|--|--|--|-]
// Calculates the MSD at starting n and lagtime tau and updates
// This function is overloaded because at early and late times there is only one MSD to be calculated
// but halfway through the simulation we average over five different starting points.
{
    for(int i=0; i<noCells; i++)
    {
        double dxrel = (cell[i].xreal - x0[i] - (xavg - xnavg[0]));
        double dyrel = (cell[i].yreal - y0[i] - (yavg - ynavg[0]));
        double dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau0] += vx0[i]*(cell[i].velx - vcx) + vy0[i]*(cell[i].vely - vcy);
        MSD[tau0] += dr2;
        MSDerr[tau0] += dr2*dr2;

        cell[i].update(deltat, lover2, l);
    }
    MSDcounter[tau0] += 1;
}

void Engine::calcMSD_and_update2()
// [--|--|++|--|--|--|--|--|--|--|-]
// Calculates the MSD at starting n and lagtime tau and updates
// This function is overloaded because at early and late times there is only one MSD to be calculated
// but halfway through the simulation we average over five different starting points.
{
    for(int i=0; i<noCells; i++)
    {
        double dxrel = (cell[i].xreal - x0[i] - (xavg - xnavg[0]));
        double dyrel = (cell[i].yreal - y0[i] - (yavg - ynavg[0]));
        double dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau0] += vx0[i]*(cell[i].velx - vcx) + vy0[i]*(cell[i].vely - vcy);
        MSD[tau0] += dr2;
        MSDerr[tau0] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[noCells+i] - (xavg - xnavg[1]));
        dyrel = (cell[i].yreal - y0[noCells+i] - (yavg - ynavg[1]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau1] += vx0[noCells+i]*(cell[i].velx - vcx) + vy0[noCells+i]*(cell[i].vely - vcy);
        MSD[tau1] += dr2;
        MSDerr[tau1] += dr2*dr2;

        cell[i].update(deltat, lover2, l);
    }
    MSDcounter[tau0] += 1;
    MSDcounter[tau1] += 1;
}

void Engine::calcMSD_and_update3(int n2)
// [--|--|--|++|--|--|--|--|--|--|-]
// Calculates the MSD at starting n and lagtime tau and updates
// This function is overloaded because at early and late times there is only one MSD to be calculated
// but halfway through the simulation we average over five different starting points.
{
    for(int i=0; i<noCells; i++)
    {
        double dxrel = (cell[i].xreal - x0[i] - (xavg - xnavg[0]));
        double dyrel = (cell[i].yreal - y0[i] - (yavg - ynavg[0]));
        double dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau0] += vx0[i]*(cell[i].velx - vcx) + vy0[i]*(cell[i].vely - vcy);
        MSD[tau0] += dr2;
        MSDerr[tau0] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[noCells+i] - (xavg - xnavg[1]));
        dyrel = (cell[i].yreal - y0[noCells+i] - (yavg - ynavg[1]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau1] += vx0[noCells+i]*(cell[i].velx - vcx) + vy0[noCells+i]*(cell[i].vely - vcy);
        MSD[tau1] += dr2;
        MSDerr[tau1] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[n2+i] - (xavg - xnavg[2]));
        dyrel = (cell[i].yreal - y0[n2+i] - (yavg - ynavg[2]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau2] += vx0[n2+i]*(cell[i].velx - vcx) + vy0[n2+i]*(cell[i].vely - vcy);
        MSD[tau2] += dr2;
        MSDerr[tau2] += dr2*dr2;

        cell[i].update(deltat, lover2, l);
    }
    MSDcounter[tau0] += 1;
    MSDcounter[tau1] += 1;
    MSDcounter[tau2] += 1;
}

void Engine::calcMSD_and_update4(int n2, int n3)
// [--|--|--|--|++|--|--|--|--|--|-]
// Calculates the MSD at starting n and lagtime tau and updates
// This function is overloaded because at early and late times there is only one MSD to be calculated
// but halfway through the simulation we average over five different starting points.
{
    for(int i=0; i<noCells; i++)
    {
        double dxrel = (cell[i].xreal - x0[i] - (xavg - xnavg[0]));
        double dyrel = (cell[i].yreal - y0[i] - (yavg - ynavg[0]));
        double dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau0] += vx0[i]*(cell[i].velx - vcx) + vy0[i]*(cell[i].vely - vcy);
        MSD[tau0] += dr2;
        MSDerr[tau0] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[noCells+i] - (xavg - xnavg[1]));
        dyrel = (cell[i].yreal - y0[noCells+i] - (yavg - ynavg[1]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau1] += vx0[noCells+i]*(cell[i].velx - vcx) + vy0[noCells+i]*(cell[i].vely - vcy);
        MSD[tau1] += dr2;
        MSDerr[tau1] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[n2+i] - (xavg - xnavg[2]));
        dyrel = (cell[i].yreal - y0[n2+i] - (yavg - ynavg[2]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau2] += vx0[n2+i]*(cell[i].velx - vcx) + vy0[n2+i]*(cell[i].vely - vcy);
        MSD[tau2] += dr2;
        MSDerr[tau2] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[n3+i] - (xavg - xnavg[3]));
        dyrel = (cell[i].yreal - y0[n3+i] - (yavg - ynavg[3]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau3] += vx0[n3+i]*(cell[i].velx - vcx) + vy0[n3+i]*(cell[i].vely - vcy);
        MSD[tau3] += dr2;
        MSDerr[tau3] += dr2*dr2;

        cell[i].update(deltat, lover2, l);
    }
    MSDcounter[tau0] += 1;
    MSDcounter[tau1] += 1;
    MSDcounter[tau2] += 1;
    MSDcounter[tau3] += 1;
}

void Engine::calcMSD_and_update5(int n2, int n3, int n4)
// [--|--|--|--|--|+++++++++++++++-]
// Calculates the MSD at starting n and lagtime tau and updates
// This function is overloaded because at early and late times there is only one MSD to be calculated
// but halfway through the simulation we average over five different starting points.
{
    for(int i=0; i<noCells; i++)
    {
        double dxrel = (cell[i].xreal - x0[i] - (xavg - xnavg[0]));
        double dyrel = (cell[i].yreal - y0[i] - (yavg - ynavg[0]));
        double dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau0] += vx0[i]*(cell[i].velx - vcx) + vy0[i]*(cell[i].vely - vcy);
        MSD[tau0] += dr2;
        MSDerr[tau0] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[noCells+i] - (xavg - xnavg[1]));
        dyrel = (cell[i].yreal - y0[noCells+i] - (yavg - ynavg[1]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau1] += vx0[noCells+i]*(cell[i].velx - vcx) + vy0[noCells+i]*(cell[i].vely - vcy);
        MSD[tau1] += dr2;
        MSDerr[tau1] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[n2+i] - (xavg - xnavg[2]));
        dyrel = (cell[i].yreal - y0[n2+i] - (yavg - ynavg[2]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau2] += vx0[n2+i]*(cell[i].velx - vcx) + vy0[n2+i]*(cell[i].vely - vcy);
        MSD[tau2] += dr2;
        MSDerr[tau2] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[n3+i] - (xavg - xnavg[3]));
        dyrel = (cell[i].yreal - y0[n3+i] - (yavg - ynavg[3]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau3] += vx0[n3+i]*(cell[i].velx - vcx) + vy0[n3+i]*(cell[i].vely - vcy);
        MSD[tau3] += dr2;
        MSDerr[tau3] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[n4+i] - (xavg - xnavg[4]));
        dyrel = (cell[i].yreal - y0[n4+i] - (yavg - ynavg[4]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau4] += vx0[n4+i]*(cell[i].velx - vcx) + vy0[n4+i]*(cell[i].vely - vcy);
        MSD[tau4] += dr2;
        MSDerr[tau4] += dr2*dr2;

        cell[i].update(deltat, lover2, l);
    }
    MSDcounter[tau0] += 1;
    MSDcounter[tau1] += 1;
    MSDcounter[tau2] += 1;
    MSDcounter[tau3] += 1;
    MSDcounter[tau4] += 1;
}

void Engine::calcMSD_and_update6(int n2, int n3, int n4)
// [--|--|--|--|--|--|--|--|--|--|+]
// Calculates the MSD at starting n and lagtime tau and updates
// This function is overloaded because at early and late times there is only one MSD to be calculated
// but halfway through the simulation we average over five different starting points.
{
    for(int i=0; i<noCells; i++)
    {
        double dxrel = (cell[i].xreal - x0[i] - (xavg - xnavg[0]));
        double dyrel = (cell[i].yreal - y0[i] - (yavg - ynavg[0]));
        double dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau0] += vx0[i]*(cell[i].velx - vcx) + vy0[i]*(cell[i].vely - vcy);
        MSD[tau0] += dr2;
        MSDerr[tau0] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[noCells+i] - (xavg - xnavg[1]));
        dyrel = (cell[i].yreal - y0[noCells+i] - (yavg - ynavg[1]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau1] += vx0[noCells+i]*(cell[i].velx - vcx) + vy0[noCells+i]*(cell[i].vely - vcy);
        MSD[tau1] += dr2;
        MSDerr[tau1] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[n2+i] - (xavg - xnavg[2]));
        dyrel = (cell[i].yreal - y0[n2+i] - (yavg - ynavg[2]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau2] += vx0[n2+i]*(cell[i].velx - vcx) + vy0[n2+i]*(cell[i].vely - vcy);
        MSD[tau2] += dr2;
        MSDerr[tau2] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[n3+i] - (xavg - xnavg[3]));
        dyrel = (cell[i].yreal - y0[n3+i] - (yavg - ynavg[3]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau3] += vx0[n3+i]*(cell[i].velx - vcx) + vy0[n3+i]*(cell[i].vely - vcy);
        MSD[tau3] += dr2;
        MSDerr[tau3] += dr2*dr2;
        dxrel = (cell[i].xreal - x0[n4+i] - (xavg - xnavg[4]));
        dyrel = (cell[i].yreal - y0[n4+i] - (yavg - ynavg[4]));
        dr2 = dxrel*dxrel + dyrel*dyrel;
        vaf[tau4] += vx0[n4+i]*(cell[i].velx - vcx) + vy0[n4+i]*(cell[i].vely - vcy);
        MSD[tau4] += dr2;
        MSDerr[tau4] += dr2*dr2;

        double velmag = sqrt(cell[i].velx*cell[i].velx + cell[i].vely*cell[i].vely);
        double veldir = atan2(cell[i].vely, cell[i].velx);
        print.print_vid(cell[i].posx, cell[i].posy, cell[i].psi, 50*velmag, veldir, cell[i].over);
        print.print_Ovito(i, cell[i].R, cell[i].posx, cell[i].posy, cell[i].psi, 50*velmag, veldir, cell[i].over);
        print.print_positions(cell[i].posx, cell[i].posy);
        cell[i].over = 240;

        cell[i].update(deltat, lover2, l);
    }
    MSDcounter[tau0] += 1;
    MSDcounter[tau1] += 1;
    MSDcounter[tau2] += 1;
    MSDcounter[tau3] += 1;
    MSDcounter[tau4] += 1;
}

void Engine::save_config()
// Two files that keep track of the size, position and directions of all particles while the simulation is running.
// If the simulation crashes, at least one of the files can be used to reinitiate the simulation without losing 
// already gathered data.
{
    static int phase = 0;
    ofstream config;
    
    if(phase)
        config.open((path+"/"+run+"/dat/config_1.dat").c_str());
    else
        config.open((path+"/"+run+"/dat/config_2.dat").c_str());
    
    if(config.fail())
    {
        cout << "Failed opening " << path+"/" << run << "/dat/config.dat";
        exit(717);
    }
    
    config << timeCounter << "\n";
    
    for(int i=0; i<noCells; ++i)
    {
        config << setprecision(8) << cell[i].posx << "\t" << cell[i].posy << "\t" << cell[i].psi << "\n";
    }
    
    config.close();
    phase = 1 - phase;
}

void Engine::save_config2()
// File that keeps track of the dynamics
{
    static int fileNumber = 0;
    string stringFileNumber = boost::lexical_cast<string>(fileNumber);
    ofstream config;
    
    config.open((path+"/"+run+"/dat/config_"+stringFileNumber+"1.dat").c_str());
    
    if(config.fail())
    {
        cout << "Failed opening " << path+"/" << run << "/dat/config.dat";
        exit(717);
    }
    
    config << timeCounter << "\n";
    
    for(int i=0; i<noCells; ++i)
    {
        config << setprecision(8) << cell[i].posx << "\t" << cell[i].posy << cell[i].velx << "\t" << cell[i].vely << "\t" << cell[i].R << "\t" << cell[i].psi << "\n";
    }
    
    config.close();
    ++fileNumber;
}

int main(int argc, char *argv[])
{
    string ID = "";
    string dir = "";
    int N = 0;
    if(argc != 8)
    {
        cout<< "Incorrect number of arguments. Need:" << endl
            << "- full run ID" << endl
	    << "- ID" << endl << "- number of cells" << endl << "- number of iterations" << endl
            << "- iterations per unit time" << endl << "- \\lambda_s" << endl
            << "- \\lambda_n" << endl << "- \\rho" << endl
            << "Program exit status (1)" << endl;
        return 1;
    }
    else
    {
        dir = argv[1];
        ID             = argv[2];
        N              = argv[3];
        long int steps = atol(argv[3]);
        int stepsPTime = atoi(argv[4]);
        double l_s     = atof(argv[5]);
        double l_n     = atof(argv[6]);
        double rho     = atof(argv[7]);

        Engine engine(dir, ID, N, steps, stepsPTime, l_s, l_n, rho);
        engine.start();
    }
    
    int check = system(("/home/dmccusker/remote/jamming-dynamics/code/plot/plot_jam_act_03.gnu "+ID+ " " +dir+" "+steps).c_str());
    if( check != 0 )
        cout << "An error occurred while plotting (one of) the graphs\n";
    return 0;

}
