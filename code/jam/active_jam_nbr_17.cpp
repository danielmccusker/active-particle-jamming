#include <vector>
#include <cmath>
#include <ctime>
#include <chrono>
#include <time.h>

#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include "../classes/Cell.h"
#include "../classes/Box.h"
#include "../classes/Print.h"
//#include "../classes/Correlations.h"
//#include "../classes/Geometry.h"
//#include "../classes/Physics.h"
//#include "../classes/Fluctuations.h"

#define PI 3.14159265
#define sqrt2 1.41421356

using namespace std;
using namespace std::chrono;

// Mersenne Twister pseudo-random number generator
std::chrono::time_point<std::chrono::high_resolution_clock> t1 = std::chrono::high_resolution_clock::now();
boost::mt19937 gen(std::chrono::duration_cast<std::chrono::nanoseconds>(t1.time_since_epoch()).count());
boost::uniform_real<> unidist(-PI, PI);
boost::normal_distribution<> normdist(0, 1);
boost::variate_generator< boost::mt19937, boost::uniform_real<> > randuni(gen, unidist);
boost::variate_generator< boost::mt19937, boost::normal_distribution<> > randnorm(gen, normdist);

const string location = "/Users/Daniel1/Desktop/ActiveMatterResearch/jamming-dynamics/";
//const string location = "/home/dmccusker/remote/jamming-dynamics/";

struct Engine
{
    // Engine
    Engine(string, string, long int, long int, double, double, double);
    ~Engine();
    
    int N;                                  // Number of cells
    double CFself;                          // Self-propulsion force
    double CTnoise;                         // Noise parameter
    double dens;                            // "Packing fraction" / average density
    string run;                             // Run with current set of parameters
    string fullRun;                         // Entire set of runs for a given particle number
    vector<Cell> cell;                      // "cell" is vector of "Cell" type objects
    
    long int totalSteps;
    long int countdown;
    long int t;                             // Current time
    int nSkip;                              // Print data every nSkip steps
    long int resetCounter;                  // Records the number of times the Verlet skin list is refreshed
    double dt;                              // Time step, in units of cell-cell repulsion time
    int film;                               // Make a video of the last time steps
    
    double xavg0, yavg0;                    // Initial position of center of mass
    double xavg, yavg;                      // Current position of center of mass, with no PBC
    double xavgold, yavgold;                // Stores values of average positions for Verlet list skin refresh
    
    double vcx, vcy;
    double u0x, u0y;                              // Initial order parameter
    double u0avg;
    
    void start();
    void relax();
    void assignCellsToGrid();
    void buildVerletLists();
    bool newSkinList();
    void neighborInteractions();
    void calculateCenterOfMass();
    void calculateInitialOrientation();
    void saveOldPositions();
    
    // Geometry
    
    vector<Box> grid;                       // Stores topology of simulation area
    vector<vector<int>> boxPairs;           // List of box pairs that are separated by less than a correlation cut-off
    void defineGrid();
    void defineCells();
    double delta_norm(double);

    double L;                               // Length of the simulation area
    double Lover2;                          // Read: "L-over-two", so we don't have to calculate L/2 every time we need it
    double lp;                              // Length of a box
    int b;                                  // Number of boxes in one dimension
    int b2;                                 // Total number of boxes
    double rn;                              // Radius that defines particle's interaction neighborhood
    double rs;                              // Verlet skin radius
    double rn2;                             // Square distances, to avoid calculating too many square roots
    double rs2;
    
    // Physics functions
    double calculateOrderParameter();
    double calculateOrientation();
    void calculateInitialVelocities();
    void calculateCOMvelocity();
    void MSD(int t);
    void velDist(int init);
    double orderAvg;
    double order2Avg;
    double order4Avg;
    double binder;
    double variance;
    double pressure;
    double pressureAvg;
    double pressure2Avg;
    double pressureVariance;
    
    // Fluctuations
    void densityFluctuations();
    void GNFinit();
    int GNFcounter;                         //counts time steps for calculating rms density fluctuations
    double fluctuationValue;
    double Rad;                              //radius of GNF intersection circle
    double RadMax, RadMin, RadInt;
    double numberOfGNFpoints;               //how many divisions of total run time for each density fluctuation measurement
    double timeOneDensityMeasurement;
    double cc_overlap(double R1, double R2, double r);
    
    // Correlation functions
    
    const int timeAvg = 100;                // Number of averages of static physics functions
    const int tCorrelation = 100;           // Number of time steps for calculating time-dependent correlations
    const int cutoff = 150;
    int tau;                                // Starting time for time-dependent correlations
    int tCorrCounter;
    
    void spatialCorrelations(int init);
    void VAF(int init, int t);
    void OAF(int init, int t);
    vector<double> velocityCorrelationValues;
    vector<double> velocityCorrelationCorrectedValues;
    vector<double> angleCorrelationValues;
    vector<double> angleCorrelationCorrectedValues;
    vector<double> counts;
    vector<double> pairCorrelationValues;
    vector<double> orderAutocorrelationValues;
    vector<double> velocityAutocorrelationValues;
    vector<double> velocityDistributionValues;
    
    Print print;                            // Print class has methods for printing data
    void printPhysicsFunctions();
};

Engine::Engine(string dir, string ID, long int n, long int steps, double l_s, double l_n, double rho){
    // Constructor
    
    fullRun     = dir;
    run         = ID;
    N           = n;
    totalSteps  = steps;
    countdown   = steps+1;
    dt          = 0.1;              // Fix 10 steps per unit time
    CFself      = l_s;
    CTnoise     = l_n;
    dens        = rho;
    
    nSkip = 100;
    film  = 100*nSkip;
    t = 0;
    resetCounter = 0;
    tau = totalSteps/timeAvg;
    tCorrCounter = tCorrelation;
    
    rn = 2.8;                       // Neighborhood radius is fixed at 2.8.
    rs = 1.5*rn;					// Optimize skin radius choice for the number of particles:
    rn2 = rn*rn;                    // More particles -> assignCellsToGrid takes more time.
    rs2 = rs*rs;                    // A larger value of rs then reduces the number of calls
                                    // to assignCellsToGrid.
    orderAvg = 0.0;
    order2Avg = 0.0;
    order4Avg = 0.0;
    binder = 0.0;
    variance = 0.0;
    pressure = 0.0;
    pressureAvg = 0.0;
    pressure2Avg = 0.0;
    pressureVariance = 0.0;
    
    vcx = 0.0;
    vcy = 0.0;
    u0x = 0.0;
    u0y = 0.0;
    u0avg = 0.0;
    
    pairCorrelationValues.reserve(200);
    velocityCorrelationValues.reserve(200);
    orderAutocorrelationValues.reserve(200);
    velocityAutocorrelationValues.reserve(200);
    angleCorrelationValues.reserve(200);
    angleCorrelationCorrectedValues.reserve(200);
}

Engine::~Engine(){
    // Destructor
    
    for (int i=0; i<N;  i++) cell[i].VerletList.clear();
    for (int j=0; j<b2; j++) grid[j].CellList.clear();
}

void Engine::start(){
    
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    
    print.init(location, fullRun, run, N);
    defineCells();
    defineGrid();
    GNFinit();
    
    assignCellsToGrid();
    buildVerletLists();
    
    relax();
    
    int init = 0;
    int init2 = 0;
    while(countdown != 0){                                          // Time loop
        
        if( newSkinList() ){
            assignCellsToGrid();
            buildVerletLists();
        }
        
        neighborInteractions();
        
        for(int i=0; i<N; i++) cell[i].update(dt, Lover2, L);      // Integrate equations of motion
        
        calculateCenterOfMass();
        
        if(t%nSkip == 0){                                          // Print data every nSkip steps
            
            double orientation = calculateOrientation();
            double order = calculateOrderParameter();
            double order2 = order*order;
            orderAvg+=order;
            order2Avg+=order2;
            order4Avg+=order2*order2;
            
            pressure/=2*PI*(double)N;
            pressureAvg+=pressure;
            pressure2Avg+=(pressure*pressure);
            
            densityFluctuations();
            MSD(t);
            print.print_data(t, xavg, yavg, order, orientation,
                             cell[0].xreal, cell[0].yreal, cell[N/3].xreal, cell[N/3].yreal, cell[2*N/3].xreal, cell[2*N/3].yreal, pressure);
            // Time, center of mass position, order parameter, system orientation, three cell trajectories
            
            if(countdown<film)                  // Print video for the last part of the run
            {
                int k=0;
                for(int i=0; i<N; i++){
                    print.print_Ovito(k, N, i, cell[i].posx, cell[i].posy, cell[i].R, cell[i].over, cell[i].velx, cell[i].vely);
                    cell[i].over = 240;
                    k=1;
                }
            }
        }
        
        if( t%tau == 0 && t!=0 )                // Calculate static physics functions and initialize time correlation functions
        {
            saveOldPositions();                 // Use the grid to build list of pairs: refresh box numbers
            assignCellsToGrid();
            buildVerletLists();
            
            spatialCorrelations(init);
            velDist(init);
            init++;
            
            calculateInitialOrientation();
            calculateCOMvelocity();
            calculateInitialVelocities();
            tCorrCounter=0;
        }
        
        if( tCorrCounter<tCorrelation )         // Calculate time correlation functions
        {
            calculateCOMvelocity();
            VAF(init2, tCorrCounter);
            OAF(init2, tCorrCounter);
            init2++;
            tCorrCounter++;
        }
        
        t++;
        countdown--;
    }
    
    printPhysicsFunctions();
    orderAvg  /= ((double)totalSteps/(double)nSkip);
    order2Avg /= ((double)totalSteps/(double)nSkip);
    order4Avg /= ((double)totalSteps/(double)nSkip);
    
    pressureAvg  /= ((double)totalSteps/(double)nSkip);
    pressure2Avg /= ((double)totalSteps/(double)nSkip);
    pressureVariance = pressure2Avg - pressureAvg*pressureAvg;
    binder = 1.0 - order4Avg/(3.0*order2Avg*order2Avg);
    variance = order2Avg - orderAvg*orderAvg;
    
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<seconds>( t2 - t1 ).count();
    
    print.print_summary(run, N, L, t, 1./dt, CFself, CTnoise, dens, duration, resetCounter, velocityCorrelationValues[0],
                        binder, orderAvg, variance, pressureAvg, pressureVariance);
}

void Engine::defineCells(){
    
    double area = 0;
    
    for(int i=0; i<N; i++){                             // Create vector of cell objects and assign their radii
        Cell *newCell = new Cell;
        cell.push_back(*newCell);
        cell[i].index = i;
        
        double cellRad = 1. + randnorm()/10;            // Cell radii follow a Gaussian distribution
        double R2 = cellRad*cellRad;
        area += R2;                                     // Simulation area depends on cell sizes
        cell[i].R  = cellRad;
    }
    
    L = sqrt(PI * area / dens);                         // L depends on the given density and calculated area
    Lover2 = L/2.0;
    
    int sqrtN = sqrt(N);
    double spacing = L / sqrtN;
    
    xavg0 = 0;
    yavg0 = 0;
    xavgold = 0;
    yavgold = 0;
    
    for(int i=0; i<N; i++){                             // Set up initial cell configuration
        cell[i].posx    = -Lover2 + spacing*(i%sqrtN) + randnorm()/3.;
        cell[i].posy    = -Lover2 + spacing*(i/sqrtN) + randnorm()/3.;
        cell[i].xreal   = cell[i].posx;
        cell[i].yreal   = cell[i].posy;
        cell[i].x0      = cell[i].posx;
        cell[i].y0      = cell[i].posy;
        cell[i].xold    = cell[i].posx;
        cell[i].yold    = cell[i].posy;
        cell[i].psi     = randuni();
        
        if(cell[i].posx >= Lover2)      cell[i].posx -= L;
        else if(cell[i].posx < -Lover2) cell[i].posx += L;
        if(cell[i].posy >= Lover2)      cell[i].posy -= L;
        else if(cell[i].posy < -Lover2) cell[i].posy += L;
        if(cell[i].psi >= PI)           cell[i].psi  -= 2*PI;
        else if(cell[i].psi < -PI)      cell[i].psi  += 2*PI;
        
        cell[i].cosp     = cos(cell[i].psi);
        cell[i].sinp     = sin(cell[i].psi);
        cell[i].cosp_new = cell[i].cosp;                // Average direction of particles in nbrhood also includes itself
        cell[i].sinp_new = cell[i].sinp;
        
        xavg0 += cell[i].xreal;
        yavg0 += cell[i].yreal;
    }
    
    xavg0 = xavg0/N;
    yavg0 = yavg0/N;
    xavg    = xavg0;
    yavg    = yavg0;
    xavgold = xavg0;
    yavgold = yavg0;
}

void Engine::defineGrid()
{
    
//index=z*width*height+y*width+x

    
    lp = 2*rn;                                          // lp is ~at least~ the assigned neighbor region diameter.
    b = static_cast<int>(floor(L/lp));                  // It can be a little bit bigger such that we have
    b2 = b*b;                                           // an integer number b2 of equally-sized boxes.
    lp = L/floor(L/lp);
    
    for (int k=0; k<b2; k++) {
        Box *newBox = new Box;
        grid.push_back(*newBox);
        grid[k].serial_index = k;
    }
    
    for(int i=0; i<b; i++){                             // Start in bottom left of simulation box and move to the right
        for(int j=0; j<b; j++){
            grid[i+(j*b)].vector_index[0] = i;
            grid[i+(j*b)].vector_index[1] = j;
            
            grid[i+(j*b)].xmin = -Lover2 + i*L/b;
            grid[i+(j*b)].xmax = -Lover2 + (i+1)*L/b;
            grid[i+(j*b)].center[0] = (grid[i+(j*b)].xmin + grid[i+(j*b)].xmax) / 2.;
            
            grid[i+(j*b)].ymin = -Lover2 + j*L/b;
            grid[i+(j*b)].ymax = -Lover2 + (j+1)*L/b;
            grid[i+(j*b)].center[1] = (grid[i+(j*b)].ymin + grid[i+(j*b)].ymax) / 2.;
        }
    }
    
    for(int k=0; k<b2; k++){                            // Find each box's adjacent boxes
        for (int j=0; j<3; j++) {
            for (int i=0; i<3; i++){
                int p = i;
                int q = j;
                
                if (grid[k].vector_index[0] == 0   && i==0) p=p+b;
                if (grid[k].vector_index[0] == b-1 && i==2) p=p-b;
                if (grid[k].vector_index[1] == 0   && j==0) q=q+b;
                if (grid[k].vector_index[1] == b-1 && j==2) q=q-b;
                
                grid[k].neighbors[i+(j*3)] = k + (p-1) + (q-1)*b;
            }
        }
    }
    
    boxPairs.reserve(b*(b+1)/2);                       // Max number of possible box pairs
    
    for (int p=0; p<b2; p++)
    {
        for(int q=p; q<b2; q++)
        {
            double deltax = delta_norm(grid[p].center[0]-grid[q].center[0]);
            double deltay = delta_norm(grid[p].center[1]-grid[q].center[1]);
            double boxDist2 = deltax*deltax+deltay*deltay;
            double boxCutoff = cutoff+(sqrt2*lp);
            if (boxDist2 < boxCutoff*boxCutoff)
            {
                vector<int> temp;
                temp.assign(2,0);
                temp[0] = grid[p].serial_index;
                temp[1] = grid[q].serial_index;
                boxPairs.push_back(temp);
            }
        }
    }
                
}

void Engine::relax()
// Relax the system as passive particles for 10000*unit relaxation time to allow many rearrangements
// Then allow to thermalize for 1,000,000 steps, slowly increasing activity to final value
{
    int trelax = (int)(1000.0/dt);
    int tthermalize = 1e4;
    
    double CFself_old = CFself;
    
    CFself = 0;
    
    for(int t_=0; t_<trelax+tthermalize; t_++)
    {
        if(t_<trelax)
        {
            CFself = 0;
        }
        else
        {
            CFself = CFself_old - (tthermalize - (t_-trelax))*CFself_old/tthermalize;
        }
       
        if( newSkinList() )
        {
            assignCellsToGrid();
            buildVerletLists();
        }
        
        neighborInteractions();
        
        for(int i=0; i<N; i++)
        {
            cell[i].update(dt, Lover2, L);
        }
        
        calculateCenterOfMass();
    }
    
    CFself = CFself_old;
    
    // Start cells at their current location in the grid after relaxing and thermalizing
    resetCounter = 0;
    xavg0 = 0;
    yavg0 = 0;
    xavgold = 0;
    yavgold = 0;
    for(int i=0; i<N; i++){
        cell[i].xreal = cell[i].posx;
        cell[i].yreal = cell[i].posy;
        cell[i].x0    = cell[i].posx;
        cell[i].y0    = cell[i].posy;
        cell[i].xold  = cell[i].posx;
        cell[i].yold  = cell[i].posy;
        xavg0 += cell[i].xreal;
        yavg0 += cell[i].yreal;
    }
    xavg0 = xavg0/N;
    yavg0 = yavg0/N;
    xavg    = xavg0;
    yavg    = yavg0;
    xavgold = xavg0;
    yavgold = yavg0;
}

void Engine::calculateCenterOfMass(){
    xavg = 0;
    yavg = 0;
    for(int i=0; i<N; i++)
    {
        xavg += cell[i].xreal;
        yavg += cell[i].yreal;
    }
    xavg = xavg/N;
    yavg = yavg/N;
}

void Engine::saveOldPositions(){
    xavgold = 0;
    yavgold = 0;
    for(int i=0; i<N; i++)
    {
        cell[i].xold = cell[i].xreal;
        cell[i].yold = cell[i].yreal;
        xavgold += cell[i].xreal;
        yavgold += cell[i].yreal;
    }
    xavgold /= N;
    yavgold /= N;
}

void Engine::assignCellsToGrid()
{
    for (int j=0; j<b2; j++) grid[j].CellList.clear();
    
    for (int i=0; i<N; i++)
    {
        double r2 = lp*lp/2;            // Circle of radius lp*√2/2 around each box center touches the box's corners
        for (int j=0; j<b2; j++)
        {
            double dx = cell[i].posx - grid[j].center[0];
            double dy = cell[i].posy - grid[j].center[1];
            double d2 = dx*dx+dy*dy;
            if(d2 < r2)                 // Find box to which the cell belongs
            {
                r2 = d2;
                cell[i].box = j;
            }
        }
        grid[cell[i].box].CellList.push_back(i);
    }
}

void Engine::buildVerletLists()
{
    for(int i=0; i<N; i++)
    {
        cell[i].VerletList.clear();
        for(int m=0; m<9; m++)
        {
            int p = grid[cell[i].box].neighbors[m];
            int max = grid[p].CellList.size();
            for(int k=0; k < max; k++)
            {
                int j = grid[p].CellList[k];
                if(j > i)
                {
                    double dx = delta_norm(cell[j].posx-cell[i].posx);
                    double dy = delta_norm(cell[j].posy-cell[i].posy);
                    if( dx*dx+dy*dy < rs2 )
                    {
                        cell[i].VerletList.push_back(j);
                        cell[j].VerletList.push_back(i);
                    }
                }
            }
        }
    }
}

bool Engine::newSkinList()
{
// Compare the two largest displacements to see if a skin refresh is required
// Refresh if any particle may have entered any other particle's neighborhood
    
    bool refresh = false;
    double largest2 = 0.;
    double second2 = 0.;
    for(int i=0; i<N; i++){
        double dx = delta_norm(cell[i].xreal - cell[i].xold - xavg + xavgold);
        double dy = delta_norm(cell[i].yreal - cell[i].yold - yavg + yavgold);
        double d2 = dx*dx+dy*dy;
        if(d2 > largest2)     { second2 = largest2; largest2 = d2; }
        else if(d2 > second2) { second2 = d2; }
    }
    
    if( ( sqrt(largest2)+sqrt(second2) ) > (rs-rn) ){
        resetCounter++;
        saveOldPositions();
        refresh = true;
    }
    
    return refresh;
}

void Engine::neighborInteractions()
// * Most physics happens here *
{
    pressure = 0.0;
    for(int i=0; i<N; i++){
        int max = cell[i].VerletList.size();
        for(int k=0; k < max; k++){                                 // Check each cell's Verlet list for neighbors
            int j = cell[i].VerletList[k];
            if(j > i){                                              // Symmetry reduces calculations by half
                double dx = delta_norm(cell[j].posx-cell[i].posx);
                double dy = delta_norm(cell[j].posy-cell[i].posy);
                double d2 = dx*dx+dy*dy;
                if(d2 < rn2){                                       // They're neighbors
                    double sumR = cell[i].R + cell[j].R;
                    if( d2 < (sumR*sumR) ){                         // They also overlap
                        double overlap = sumR / sqrt(d2) - 1;
                        double forceX = overlap*dx;
                        double forceY = overlap*dy;
                        cell[i].Fx -= forceX;                       // Spring repulsion force
                        cell[i].Fy -= forceY;
                        cell[j].Fx += forceX;
                        cell[j].Fy += forceY;
                        pressure+=2.0*sqrt(forceX*forceX+forceY*forceY);
                        if(countdown <= film && t%nSkip == 0){
                            cell[i].over -= 240*abs(overlap);
                            cell[j].over -= 240*abs(overlap);
                        }
                    }
                    cell[i].cosp_new += cell[j].cosp;               // Add up orientations of neighbors
                    cell[i].sinp_new += cell[j].sinp;
                    cell[j].cosp_new += cell[i].cosp;
                    cell[j].sinp_new += cell[i].sinp;
                }
            }
        }
        cell[i].Fx += cell[i].cosp*cell[i].R*CFself;                // Self-propulsion force
        cell[i].Fy += cell[i].sinp*cell[i].R*CFself;
        cell[i].psi_new = atan2(cell[i].sinp_new, cell[i].cosp_new) + CTnoise*randuni();
    }
}

void Engine::MSD(int t){
    double MSD = 0.0;
    for (int i=0; i<N; i++) {
        double dx = (cell[i].xreal - cell[i].x0 - xavg + xavg0);
        double dy = (cell[i].yreal - cell[i].y0 - yavg + yavg0);
        MSD += dx*dx+dy*dy;
    }
    MSD /= N;
    print.print_MSD(t, MSD);
}

void Engine::calculateCOMvelocity(){
    vcx = 0.0;
    vcy = 0.0;
    for(int i=0; i<N; i++){
        vcx+=cell[i].velx;
        vcy+=cell[i].vely;
    }
    vcx /= (double)N;
    vcy /= (double)N;
}

void Engine::calculateInitialVelocities()
{
    for(int i=0; i<N; i++){
        cell[i].vx0 = cell[i].velx - vcx;
        cell[i].vy0 = cell[i].vely - vcy;
    }
}

void Engine::VAF(int init, int t){
    double v = 0.0;
    
    if(init==0)
    {
        velocityAutocorrelationValues.assign(tCorrelation,0);
    }
    
    for (int i=0; i<N; i++) {
        v += (cell[i].velx-vcx)*cell[i].vx0 + (cell[i].vely-vcy)*cell[i].vy0;
    }
    velocityAutocorrelationValues[t] += v/(CFself*CFself*(double)N);
}

double Engine::calculateOrderParameter(){
    double ox = 0.0;
    double oy = 0.0;
    for (int i=0; i<N; i++) {
        ox+=cell[i].velx;
        oy+=cell[i].vely;
    }
    return sqrt(ox*ox+oy*oy)/((double)N*CFself);
}

double Engine::calculateOrientation(){
    double ox = 0.0;
    double oy = 0.0;
    for (int i=0; i<N; i++) {
        ox+=cell[i].velx;
        oy+=cell[i].vely;
    }
    return atan2(oy, ox);
}

void Engine::calculateInitialOrientation(){
    double ox = 0.0;
    double oy = 0.0;
    for (int i=0; i<N; i++) {
        ox+=cell[i].velx;
        oy+=cell[i].vely;
    }
    u0x=ox/N;
    u0y=oy/N;
    u0avg+=sqrt(u0x*u0x+u0y*u0y);
}

void Engine::OAF(int init, int t){
    
    if(init==0)
    {
        orderAutocorrelationValues.assign(tCorrelation,0);
    }
    
    double ox = 0.0;
    double oy = 0.0;
    for (int i=0; i<N; i++) {
        ox+=cell[i].velx;
        oy+=cell[i].vely;
    }
    
    orderAutocorrelationValues[t] += (ox*u0x+oy*u0y)/N;
    
}

void Engine::velDist(int init)
{
    const double vmax = 2.0*CFself;
    const double dv = CFself/50.0;
    int noBins = (int)ceil(vmax/dv);
    vector<double> temp;
    temp.assign(noBins, 0);
    
    if(init == 0) velocityDistributionValues.assign(noBins, 0);
    
    for(int i=0; i<N; i++){
        double v = sqrt(cell[i].velx*cell[i].velx + cell[i].vely*cell[i].vely);
        int binv = (int)floor(v/dv);
        if(binv < noBins) temp[binv] += 1.0;
    }
    
    for(int k=0; k<noBins; k++){
        temp[k] /= (double)N;
        velocityDistributionValues[k]+=temp[k];
    }
}

void Engine::spatialCorrelations(int init)
{
    const double dr_v = 2.0;
    const double dr_p = 0.01;
    
    const int np = (int)ceil(cutoff/dr_p);
    const int nv = (int)ceil(cutoff/dr_v);
    
    vector<double> velTemp;
    vector<double> velCorrTemp;
    vector<double> angleTemp;
    vector<double> angleCorrTemp;
    vector<double> pairTemp;
    vector<double> counts;
    
    velTemp.assign(nv,0);
    velCorrTemp.assign(nv,0);
    angleTemp.assign(nv,0);
    angleCorrTemp.assign(nv,0);
    pairTemp.assign(np,0);
    counts.assign(nv,0);
    
    calculateCOMvelocity();
    
    if(init == 0)
    {
        velocityCorrelationValues.assign(nv, 0);
        velocityCorrelationCorrectedValues.assign(nv, 0);
        angleCorrelationValues.assign(nv,0);
        angleCorrelationCorrectedValues.assign(nv,0);
        pairCorrelationValues.assign(np, 0);
    }
   
    for (int p=0; p<b2; p++)
    {
        for(int q=p; q<b2; q++)
        {
            double deltax = delta_norm(grid[p].center[0]-grid[q].center[0]);
            double deltay = delta_norm(grid[p].center[1]-grid[q].center[1]);
            double boxDist2 = deltax*deltax+deltay*deltay;
            double boxCutoff = cutoff+(sqrt2*lp);
            if (boxDist2 < boxCutoff*boxCutoff)
            {
                int maxp = grid[p].CellList.size();
                for (int m=0; m<maxp; m++)
                {
                    int maxq = grid[q].CellList.size();
                    for (int n=0; n<maxq; n++)
                    {
                        int i = grid[p].CellList[m];
                        int j = grid[q].CellList[n];

                        if( p!=q || (p==q && j>i) ) // Avoid double-counting pairs in the same box
                        {
                            double dx = delta_norm(cell[j].posx-cell[i].posx);
                            double dy = delta_norm(cell[j].posy-cell[i].posy);
                            double r = sqrt(dx*dx+dy*dy);

                            int binp = (int)floor(r/dr_p);
                            int binv = (int)floor(r/dr_v);

                            if(binp < np)   // Exclude any pairs beyond cutoff, normalize
                            {
                                 pairTemp[binp] += 1.0/r;
                            }

                            if(binv < nv)
                            {
                                double vxi = cell[i].velx;
                                double vyi = cell[i].vely;
                                double vxj = cell[j].velx;
                                double vyj = cell[j].vely;
                            
                                double magi = sqrt(vxi*vxi+vyi*vyi);
                                double magj = sqrt(vxj*vxj+vyj*vyj);
                            
                                double angle = cos(cell[i].psi - cell[j].psi);
                                angleTemp[binv] += angle;
                                
                                double diff = sqrt( (vxi-vxj)*(vxi-vxj) + (vyi-vyj)*(vyi-vyj) )/(magi+magj);
                                angleCorrTemp[binv] += angle*(1.0-diff);
                                
                                double vel = vxi*vxj+vyi*vyj;
                                velTemp[binv] += vel;
                                
                                double velCOM = (vxi-vcx)*(vxj-vcx)+(vyi-vcy)*(vyj-vcy);
                                velCorrTemp[binv] += velCOM;
                                
                                counts[binv] += 1.0;
                            }
                        }
                    }
                }
            }
        }
    }
    
    for(int k=0; k<nv; k++)
    {
        angleTemp[k]/=counts[k];
        angleCorrTemp[k]/=counts[k];
        velTemp[k]/=(CFself*CFself*counts[k]);
        velCorrTemp[k]/=(CFself*CFself*counts[k]);
        
        angleCorrelationValues[k]+=angleTemp[k];
        angleCorrelationCorrectedValues[k]+=angleCorrTemp[k];
        velocityCorrelationValues[k]+=velTemp[k];
        velocityCorrelationCorrectedValues[k]+=velCorrTemp[k];
    }
    double norm = 1/(2*PI*dr_p*(double)N*rho);
    for(int k=0; k<np; k++)
    {
        pairCorrelationValues[k]+=norm*pairTemp[k];
    }
}

void Engine::GNFinit(){
// Start with a circle the size of one neighborhood
// Doesn't reach L/2 because then cells leaving the circle could immediately reenter on the other side
    
    numberOfGNFpoints = 10;
    RadMin = 3.0;
    RadMax = Lover2*4.0/5.0;
    
    Rad = RadMin;
    RadInt = pow(RadMax/RadMin, 1./(double)(numberOfGNFpoints-1));
    timeOneDensityMeasurement = (double)totalSteps/((double)nSkip*(double)numberOfGNFpoints);
    GNFcounter = 0;
    fluctuationValue = 0;
}

void Engine::densityFluctuations(){
// Draw a circle at the center of mass and calculate its intersection with the cells
    
    double expectedA = dens*PI*Rad*Rad;
    if(GNFcounter < timeOneDensityMeasurement) {
        double A = 0.0;
        for (int i=0; i<N; i++)
        {
            double dx = delta_norm(cell[i].xreal-xavg);
            double dy = delta_norm(cell[i].yreal-yavg);
            double d2 = dx*dx+dy*dy;
            if((cell[i].R+Rad)*(cell[i].R+Rad) >= d2)
            // If the circles have a nonzero overlap
            {
                A += cc_overlap(cell[i].R, Rad, sqrt(d2));
            }
        }
        fluctuationValue+=(A-expectedA)*(A-expectedA);
    } else {
        // Store value of rms fluctuation and reset for next measurement, increase circle size
        fluctuationValue=sqrt(fluctuationValue/(double)GNFcounter);
        print.print_GNF(Rad, expectedA, fluctuationValue);
        fluctuationValue = 0;
        GNFcounter = 0;
        Rad = Rad*RadInt;
    }
    GNFcounter++;
}

double Engine::cc_overlap(double cellR, double measurementR, double r)
// Formula for calculating area intersection between two circles with centers separated by r
// In our case, measurement circle radius is always larger than the cell radius
{
    if (measurementR>=cellR+r)  return PI*cellR*cellR;  // Cell lies completely in measurement circle
    else
    {
        double R12 = cellR*cellR;
        double R22 = measurementR*measurementR;
        double x = (R12 - R22 + r*r)/(2.0*r);           // base of triangle at cell center
        double theta = acos(x/cellR);                      // angle of triangle at cell center
        double A = R12*theta - x*cellR*sin(theta);
        x = r - x;                                      // Other area (at C2-side)
        theta = acos(x/measurementR);
        A += R22*theta - x*measurementR*sin(theta);

        return A;
    }
}

void Engine::printPhysicsFunctions(){
// X-values depending on k use conversion factors defined in the appropriate physics function
    
    for(int k=0; k<velocityCorrelationValues.size(); k++)
    {
        print.print_velCorr(2.0*(k+1), velocityCorrelationValues[k]/timeAvg);
    }
    for(int k=0; k<velocityCorrelationCorrectedValues.size(); k++)
    {
        print.print_velCorrC(2.0*(k+1), velocityCorrelationCorrectedValues[k]/timeAvg);
    }
    for(int k=0; k<angleCorrelationValues.size(); k++)
    {
        print.print_angleCorr(2.0*(k+1), velocityCorrelationValues[k]/timeAvg);
    }
    for(int k=0; k<angleCorrelationCorrectedValues.size(); k++)
    {
        print.print_angleCorrC(2.0*(k+1), velocityCorrelationCorrectedValues[k]/timeAvg);
    }
    for(int k=0; k<pairCorrelationValues.size(); k++)
    {
        print.print_pairCorr(0.01*k, pairCorrelationValues[k]/timeAvg);
    }
    for(int t=0; t<velocityAutocorrelationValues.size(); t++)
    {
        print.print_VAF(t, velocityAutocorrelationValues[t]/timeAvg);
    }
    for(int t=0; t<orderAutocorrelationValues.size(); t++)
    {
        print.print_OAF(t, (orderAutocorrelationValues[t]/timeAvg)-(u0avg/timeAvg)*(u0avg/timeAvg));
    }
    for(int k=0; k<velocityDistributionValues.size(); k++)
    {
        print.print_velDist(k*CFself/50.0, velocityDistributionValues[k]/timeAvg);
    }
}

double Engine::delta_norm(double delta)
 // Subtracts multiples of the box size to account for periodic boundary conditions
{
    int k=-1;
    if(delta < -Lover2) k=1;
    while(delta < -Lover2 || delta >= Lover2) delta += k*L;
    
    return delta;
}

int main(int argc, char *argv[]){
    
    string dir = "";
    string ID = "";
    long int n = 0;
    long int steps = 0;
    double l_s = 0;
    double l_n = 0;
    double rho = 0;
    
    if(argc != 8){
        cout    << "Incorrect number of arguments. Need: " << endl
        << "- full run ID" << endl
        << "- single run ID" << endl
        << "- number of cells" << endl
        << "- number of steps" << endl
        << "- \\lambda_s" << endl
        << "- \\lambda_n" << endl
        << "- \\rho" << endl
        << "Program exit status (1)" << endl;
        return 1;
    } else {
        
        dir     = argv[1];
        ID      = argv[2];
        n       = atol(argv[3]);
        steps   = atol(argv[4]);
        l_s     = atof(argv[5]);
        l_n     = atof(argv[6]);
        rho     = atof(argv[7]);
        
        Engine engine(dir, ID, n, steps, l_s, l_n, rho);
        engine.start();
    }
    
    int check = system((location+"code/plot/plot_jam_act_03.gnu "+ID+ " " +dir+" "+std::to_string(steps)).c_str());
    if( check != 0 ) cout << "An error occurred while plotting (one of) the graphs\n";
    
    return 0;
}

