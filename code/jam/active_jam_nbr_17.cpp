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

#define PI 3.14159265

using namespace std;
using namespace std::chrono;

// Mersenne Twister pseudo-random number generator
std::chrono::time_point<std::chrono::high_resolution_clock> t1 = std::chrono::high_resolution_clock::now();
boost::mt19937 gen(std::chrono::duration_cast<std::chrono::nanoseconds>(t1.time_since_epoch()).count());
boost::uniform_real<> unidist(-PI, PI);
boost::normal_distribution<> normdist(0, 1);
boost::variate_generator< boost::mt19937, boost::uniform_real<> > randuni(gen, unidist);
boost::variate_generator< boost::mt19937, boost::normal_distribution<> > randnorm(gen, normdist);

//const string location = "/Users/Daniel1/Desktop/ActiveMatterResearch/jamming-dynamics/";
const string location = "/home/dmccusker/remote/jamming-dynamics/";

struct Engine
{
    Engine(string, string, long int, long int, double, double, double);
    ~Engine();
    
    // Simulation functions
    void defineGrid();
    void defineCells();
    void start();
    void relax();
    void GNFinit();
    
    void assignCellsToGrid();
    void buildVerletLists();
    bool newSkinList();
    void neighborInteractions(int k);
    void calculateCenterOfMass();
    void saveOldPositions();
    
    double delta_norm(double);
    double cc_overlap(double R1, double R2, double r);
    
    // Statistical physics functions
    double calculateOrderParameter();
    double calculateOrientation();
    void calculateInitialVelocities();
    void calculateCOMvelocity();
    void densityFluctuations();
    void pairCorr();
    void MSD(int t);
    void VAF(int t, double, double);
    void OAF(int t, double);
    void velDist(int init);
    void spatialCorrelations(int init);
    void printPhysicsFunctions();
    vector<double> velocityCorrelationValues;
    vector<double> pairCorrelationValues;
    vector<double> orderAutocorrelationValues;
    vector<double> velocityAutocorrelationValues;
    vector<double> velocityDistributionValues;
    
    // Order parameter variance and Binder cumulant
    double orderAvg;
    double order2Avg;
    double order4Avg;
    double binder;
    double variance;
    
    // Variables
    Print print;                            // Print class has methods for printing data
    string run;                             // Run with current set of parameters
    string fullRun;                         // Entire set of runs for a given particle number
    
    int film;                               // Make a video of the last time steps
    long int totalSteps;
    long int countdown;
    long int t;                             // Current time
    int nSkip;                              // Print data every nSkip steps
    long int resetCounter;                  // Records the number of times the Verlet skin list is refreshed
    double dt;                              // Time step, in units of cell-cell repulsion time
    const int timeAvg = 1;                  // Number of time steps for averaging static physics functions
    const int tCorrelation = 20;            // Number of time steps for calculating time-dependent correlations
    int tau;                                // Starting time for time-dependent correlations
    int tCorrCounter;
    
    int N;                                  // Number of cells
    vector<Cell> cell;                      // "cell" is vector of "Cell" type objects
    vector<Box> grid;                       // Stores topology of simulation area
    
    double CFself;                          // Self-propulsion force
    double CTnoise;                         // Noise parameter
    double dens;                            // "Packing fraction" / average density
    
    double L;                               // Length of the simulation area
    double Lover2;                          // Read: "L-over-two", so we don't have to calculate L/2 every time we need it
    double lp;                              // Length of a box
    int b;                                  // Number of boxes in one dimension
    int b2;                                 // Total number of boxes
    double rn;                              // Radius that defines particle's interaction neighborhood
    double rs;                              // Verlet skin radius
    double rn2;                             // Square distances, to avoid calculating too many square roots
    double rs2;
    
    double xavg0, yavg0;                    // Initial position of center of mass
    double xavg, yavg;                      // Current position of center of mass, with no PBC
    double xavgold, yavgold;                // Stores values of average positions for Verlet list skin refresh
    
    double vcx, vcy;
    double vcx0, vcy0;                      // Initial velocity of center of mass
    double u0;                              // Initial order parameter
    double u0avg;
    
    int GNFcounter;                         //counts time steps for calculating rms density fluctuations
    double fluctuationValue;
    double Rad;                              //radius of GNF intersection circle
    double RadMax, RadMin, RadInt;
    double numberOfGNFpoints;               //how many divisions of total run time for each density fluctuation measurement
    double timeOneDensityMeasurement;
    
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
    
    vcx = 0.0;
    vcy = 0.0;
    vcx0 = 0.0;
    vcy0 = 0.0;
    u0 = 0.0;
    u0avg = 0.0;
    
    orderAutocorrelationValues.assign(tCorrelation,0);
    velocityAutocorrelationValues.assign(tCorrelation,0);
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
    
    // Start cells at their current location in the grid after relaxing and thermalizing
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
    
    int init = 0;
    
    while(countdown != 0){                                          // Time loop
        
        if( newSkinList() ){
            assignCellsToGrid();
            buildVerletLists();
        }
        
        neighborInteractions(1);
        
        for(int i=0; i<N; i++) cell[i].update(dt, Lover2, L);      // Integrate equations of motion
        
        calculateCenterOfMass();
        
        if(t%nSkip == 0){                                          // Print data every nSkip steps
            
            double orientation = calculateOrientation();
            double order = calculateOrderParameter();
            double order2 = order*order;
            orderAvg+=order;
            order2Avg+=order2;
            order4Avg+=order2*order2;
            
            densityFluctuations();
            MSD(t);
            print.print_data(t, xavg, yavg, order, orientation,
                             cell[0].xreal, cell[0].yreal, cell[N/3].xreal, cell[N/3].yreal, cell[2*N/3].xreal, cell[2*N/3].yreal);
            // Time, center of mass position, order parameter, system orientation, three cell trajectories
            
            if(countdown<film){
            // Print video for the last part of the run
                int k=0;
                for(int i=0; i<N; i++){
                    print.print_Ovito(k, N, i, cell[i].posx, cell[i].posy, cell[i].R, cell[i].over, cell[i].velx, cell[i].vely);
                    cell[i].over = 240;
                    k=1;
                }
            }
        }
        
        if( t%tau == 0 && t!=0 ){
        // Calculate static physics functions and initialize time correlation functions
            
            spatialCorrelations(init);
            velDist(init);
            init++;
            
            u0 = calculateOrderParameter();
            u0avg += u0;
            calculateInitialVelocities();
            
            tCorrCounter=0;
        }
        
        if( tCorrCounter<tCorrelation ){
        // Calculate time correlation functions
            VAF(tCorrCounter, vcx0, vcy0);
            OAF(tCorrCounter, u0);
            tCorrCounter++;
        }
        
        t++;
        countdown--;
    }
    
    u0avg = u0avg*u0avg/timeAvg;
    
    printPhysicsFunctions();
    orderAvg  /= ((double)totalSteps/(double)nSkip);
    order2Avg /= ((double)totalSteps/(double)nSkip);
    order4Avg /= ((double)totalSteps/(double)nSkip);
    
    
    binder = 1.0 - order4Avg/(3.0*order2Avg*order2Avg);
    variance = order2Avg - orderAvg*orderAvg;
    
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<seconds>( t2 - t1 ).count();
    
    print.print_summary(run, N, L, t, 1./dt, CFself, CTnoise, dens, duration, resetCounter, velocityCorrelationValues[0],
                        binder, orderAvg, variance);
}

void Engine::relax(){
// Relax the system for 1,000,000 steps, slowly decreasing the activity to the final value

    int throwaway = totalSteps/10;
    int trelax = 1e6;
    double CFrelax = 0.05;
    double CFself_old = CFself;                         // Store parameter values
    double CTnoise_old = CTnoise;
    t = 1;                                              // Prevent the program from producing output
    CTnoise = 0;
    
    for(int t_=0; t_<trelax; t_++){
        CFself = CFself_old + ((CFrelax - CFself_old)*(trelax - t_))/trelax;
        
        if( newSkinList() ){
            assignCellsToGrid();
            buildVerletLists();
        }
        
        neighborInteractions(0);
        
        for(int i=0; i<N; i++){
            cell[i].update(dt, Lover2, L);
            cell[i].psi = randuni();
        }
        
        calculateCenterOfMass();
        
    }
    
    CFself = CFself_old;
    CTnoise = CTnoise_old;
    
    for(int t_=0; t_<throwaway; t_++){
        
        if( newSkinList() ){
            assignCellsToGrid();
            buildVerletLists();
        }
        
        neighborInteractions(1);
        
        for(int i=0; i<N; i++) cell[i].update(dt, Lover2, L);
        
        calculateCenterOfMass();
        
    }
    
    t = 0;
    resetCounter = 0;
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
    Lover2 = L/2.;
    
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

void Engine::defineGrid(){
    
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
}

void Engine::calculateCenterOfMass(){
    xavg = 0;
    yavg = 0;
    for(int i=0; i<N; i++){
        xavg += cell[i].xreal;
        yavg += cell[i].yreal;
    }
    xavg = xavg/N;
    yavg = yavg/N;
}

void Engine::saveOldPositions(){
    xavgold = 0;
    yavgold = 0;
    for(int i=0; i<N; i++){
        cell[i].xold = cell[i].xreal;
        cell[i].yold = cell[i].yreal;
        xavgold += cell[i].xreal;
        yavgold += cell[i].yreal;
    }
    xavgold /= N;
    yavgold /= N;
}

void Engine::assignCellsToGrid(){
// At the beginning of the simulation, for every particle, we search every box to find which box
// the cell is in. Afterwards, we search only the neighboring boxes of the cell's current box to see
// if it moved boxes.  The particles are moving slowly enough that they are not travelling more than
// one box length per skin refresh, relative to the center of mass.
    
    for (int j=0; j<b2; j++) grid[j].CellList.clear();
    
    //if(k==0){                               // Search all boxes
        for (int i=0; i<N; i++){
            double r2 = lp*lp/2;            // Circle of radius lp*√2/2 around each box center touches the box's corners
            for (int j=0; j<b2; j++){
                double dx = cell[i].posx - grid[j].center[0];
                double dy = cell[i].posy - grid[j].center[1];
                double d2 = dx*dx+dy*dy;
                if(d2 < r2)  { r2 = d2; cell[i].box = j; }      // Find box to which the cell belongs
            }
            grid[cell[i].box].CellList.push_back(i);            // Add this cell to the box's list
        }
//    } else if(k==1 && refresh<lp){          // Search only neighboring boxes of the old box
//        for (int i=0; i<N; i++){
//            int oldBox = cell[i].box;
//            double r2 = lp*lp/2;
//            for (int m=0; m<9; m++){
//                int j = grid[oldBox].neighbors[m];
//                double dx = delta_norm(cell[i].xreal - (xavg + grid[j].center[0]));
//                double dy = delta_norm(cell[i].yreal - (yavg + grid[j].center[1]));
//                double d2 = dx*dx+dy*dy;
//                if(d2 < r2)  { r2 = d2; cell[i].box = j; }
//            }
//            grid[cell[i].box].CellList.push_back(i);
//        }
//    } else {
//        cout << "The cells are moving too fast or something is wrong with the dynamics. exit 100"
//        << endl;
//        exit(100);
//    }
}

void Engine::buildVerletLists(){
    
    for(int i=0; i<N; i++){
        cell[i].VerletList.clear();
        for(int m=0; m<9; m++){
            int p = grid[cell[i].box].neighbors[m];                 // Get the indices of the 9 boxes to search
            int max = grid[p].CellList.size();
            for(int k=0; k < max; k++){                             // Iterate through the cell list of each of these 9 boxes
                int j = grid[p].CellList[k];
                if(j > i){                                          // Symmetry reduces calcuations by half
                    double dx = delta_norm(cell[j].posx-cell[i].posx);
                    double dy = delta_norm(cell[j].posy-cell[i].posy);
                    if( dx*dx+dy*dy < rs2 ){
                        cell[i].VerletList.push_back(j);            // Add to cell's Verlet list
                        cell[j].VerletList.push_back(i);
                    }
                }
            }
        }
    }
}

bool Engine::newSkinList(){
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

void Engine::neighborInteractions(int k){
// Orientational interaction between neighbors and soft repulsion force between overlapping neighbors
// k = 0: while relaxing, neglect orientational interaction
    
// ~Most physics happens here~
    
    if (k==0){
        for(int i=0; i<N; i++){
            int max = cell[i].VerletList.size();
            for(int k=0; k < max; k++){
                int j = cell[i].VerletList[k];
                if(j > i){
                    double dx = delta_norm(cell[j].posx-cell[i].posx);
                    double dy = delta_norm(cell[j].posy-cell[i].posy);
                    double d2 = dx*dx+dy*dy;
                    if(d2 < rn2){
                        double sumR = cell[i].R + cell[j].R;
                        if( d2 < (sumR*sumR) ){
                            double overlap = sumR / sqrt(d2) - 1;
                            double forceX = overlap*dx;
                            double forceY = overlap*dy;
                            cell[i].Fx -= forceX;
                            cell[i].Fy -= forceY;
                            cell[j].Fx += forceX;
                            cell[j].Fy += forceY;
                        }
                    }
                }
            }
            cell[i].Fx += cell[i].cosp*cell[i].R*CFself;
            cell[i].Fy += cell[i].sinp*cell[i].R*CFself;
            cell[i].psi_new = randuni();
        }
    } else if (k==1){
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

void Engine::VAF(int t, double vcx0, double vcy0){
    double VAF = 0.0;
    calculateCOMvelocity();
    for (int i=0; i<N; i++) {
        double vx0 = cell[i].vx0 - vcx0;
        double vy0 = cell[i].vy0 - vcy0;
        double num  = (cell[i].velx-vcx)*vx0 + (cell[i].vely-vcy)*vy0;
        double norm = vx0*vx0 + vy0*vy0;
        VAF += num/norm;
    }
    velocityAutocorrelationValues[t]+=(VAF/(double)N);
}

void Engine::OAF(int t, double u0){
    
    double order = calculateOrderParameter();
    
    orderAutocorrelationValues[t]+=order*u0;
    
}

void Engine::velDist(int init){
    // Histogram; velocity distribution graphed from 0 to 2 times average velocity
    
    double vmax = 2.0*CFself;
    const double dv = CFself/100;
    int noBins = (int)ceil(vmax/dv);
    
    if(init == 0) velocityDistributionValues.assign(noBins, 0);
    
    for(int i=0; i<N; i++){
        double v = sqrt(cell[i].velx*cell[i].velx + cell[i].vely*cell[i].vely);
        int binv = (int)floor(v/dv);
        if(binv >= noBins) velocityDistributionValues[noBins-1] += 1.0;
        else velocityDistributionValues[binv] += 1.0;
    }
    
    for(int k=0; k<noBins; k++) velocityDistributionValues[k] /= (double)N;
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

void Engine::calculateInitialVelocities(){
    vcx0 = 0.0;
    vcy0 = 0.0;
    for(int i=0; i<N; i++){
        double velx = cell[i].velx;
        double vely = cell[i].vely;
        cell[i].vx0 = velx;
        cell[i].vy0 = vely;
        vcx0+=velx;
        vcy0+=vely;
    }
    vcx0 /= (double)N;
    vcy0 /= (double)N;
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

void Engine::spatialCorrelations(int init){
// For pair correlations, increment radius by 1/10 of average cell radius.
// For velocity correlations, bin cells into intervals of 2 (average cell diameter)
// We calculate these correlations in the same function to reduce the number of pairwise calculations
    
    double r0 = 1.0;                    // Initial radius
    double r = r0;
    double inc = 0.1;
    double bin = 2.0;
    double dr = inc;                    // Start incrementing for pair correlations
    double C = 1.0/(2.0*PI*double(N));  // Geometric normalization that depends on the system dimension
    double pairCutoff = 10;
    double pairCutoff2 = pairCutoff*pairCutoff;
    double velCutoff = 3.0*Lover2/8.0;
    double velCutoff2 = velCutoff*velCutoff;
    
    int noPairPoints = (pairCutoff-r0)/inc;
    int noVelPoints = (velCutoff-r0)/bin;
    
    if (init==0){
        pairCorrelationValues.assign(noPairPoints,0);
        velocityCorrelationValues.assign(noVelPoints,0);
    }
    
    calculateCOMvelocity();             // Measure cell velocities relative to COM motion
    
    int k1 = 0;
    int k2 = 0;
    
    while(r <= velCutoff+bin){
        double pairVal = 0.0;
        double velVal = 0.0;
        double avgvelxi = 0.0;
        double avgvelyi = 0.0;
        double avgvelxj = 0.0;
        double avgvelyj = 0.0;
        bool cellDiameter = (abs(fmod((r-r0),bin)-bin) < 1e-3); // True if current radius is multiple of bin
    
        for (int i=0; i < N; i++){
            for(int j=i+1; j < N; j++){
                double dx = delta_norm(cell[i].posx-cell[j].posx);
                double dy = delta_norm(cell[i].posy-cell[j].posy);
                double d2 = dx*dx+dy*dy;
                if(d2 < pairCutoff2 && abs((sqrt(d2)-r)) < dr){
                    pairVal+=1.0;
                }
                if( cellDiameter &&  d2 < velCutoff2 && abs(sqrt(d2)-r) < bin){
                    double vxi = cell[i].velx - vcx;
                    double vyi = cell[i].vely - vcy;
                    double vxj = cell[j].velx - vcx;
                    double vyj = cell[j].vely - vcy;
                    velVal += vxi*vxj+vyi*vyj;
                    avgvelxi += vxi;
                    avgvelyi += vyi;
                    avgvelxj += vxj;
                    avgvelyj += vyj;
                }
            }
        }
        
        double norm = C/(double)r;
        pairVal = norm*pairVal;
        velVal = norm*( velVal - norm*(avgvelxi*avgvelxj + avgvelyi*avgvelyj) )/(CFself*CFself);
        
        if(k1 < noPairPoints) { pairCorrelationValues[k1] += pairVal; k1++; }
        if(k2 < noVelPoints && cellDiameter) { velocityCorrelationValues[k2] += velVal; k2++; }
        
        r+=dr;
        // Increase dr once we've reached the end of the pair correlation calculation
        if (abs(pairCutoff - (r-r0)) < 1e-6)   dr = bin;
        
    }
}

void Engine::GNFinit(){
// Start with a circle the size of one neighborhood
// Doesn't reach L/2 because then cells leaving the circle could immediately reenter on the other side
    
    numberOfGNFpoints = 10;
    RadMin = 3.0;
    RadMax = Lover2*9.0/10.0;
    
    Rad = RadMin;
    RadInt = pow(RadMax/RadMin, 1./(numberOfGNFpoints-1));
    timeOneDensityMeasurement = (double)totalSteps/((double)nSkip*(double)numberOfGNFpoints);
    GNFcounter = 0;
    fluctuationValue = 0;
}

void Engine::densityFluctuations(){
// Draw a circle at the center of mass and calculate its intersection with the cells
    
    double expectedA = dens*PI*Rad*Rad;
    if(GNFcounter < timeOneDensityMeasurement) {
        double A = 0.0;
        for (int i=0; i<N; i++) {
            double dx = delta_norm(cell[i].xreal-xavg);
            double dy = delta_norm(cell[i].yreal-yavg);
            double d2 = dx*dx+dy*dy;
            if((cell[i].R+Rad)*(cell[i].R+Rad) >= d2){
            // If the circles have a nonzero overlap
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

double Engine::cc_overlap(double cellR, double measurementR, double r){
// Formula for calculating area intersection between two circles with centers separated by r
// In our case, measurement circle radius is always larger than the cell radius
    
    if (measurementR>=cellR+r)  return PI*cellR*cellR;  // Cell lies completely in measurement circle
    else {
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

double Engine::delta_norm(double delta){
    // Subtracts multiples of the box size to account for PBC
    
    int k=-1;
    if(delta < -Lover2) k=1;
    while(delta < -Lover2 || delta >= Lover2) delta += k*L;
    
    return delta;
}

void Engine::printPhysicsFunctions(){
// X-values depending on k use conversion factors defined in the appropriate physics function
    
    for(int k=0; k<velocityCorrelationValues.size(); k++)
        print.print_velCorr(2*k+1, velocityCorrelationValues[k]/timeAvg);
    
    for(int k=0; k<pairCorrelationValues.size(); k++)
        print.print_pairCorr(0.1*k+1, pairCorrelationValues[k]/timeAvg);
    
    for(int k=0; k<orderAutocorrelationValues.size(); k++)
        print.print_OAF(k, (orderAutocorrelationValues[k]-u0avg)/timeAvg);
    
    for(int k=0; k<velocityAutocorrelationValues.size(); k++)
        print.print_VAF(k, velocityAutocorrelationValues[k]/timeAvg);
    
    for(int k=0; k<velocityDistributionValues.size(); k++)
        print.print_velDist((double)k*2.0*CFself/200.0, velocityDistributionValues[k]/timeAvg);
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
        
        dir             = argv[1];
        ID              = argv[2];
        n               = atol(argv[3]);
        steps           = atol(argv[4]);
        l_s             = atof(argv[5]);
        l_n             = atof(argv[6]);
        rho             = atof(argv[7]);
        
        Engine engine(dir, ID, n, steps, l_s, l_n, rho);
        engine.start();
    }
    
    int check = system((location+"code/plot/plot_jam_act_03.gnu "+ID+ " " +dir+" "+std::to_string(steps)).c_str());
    if( check != 0 ) cout << "An error occurred while plotting (one of) the graphs\n";
    
    return 0;
}

