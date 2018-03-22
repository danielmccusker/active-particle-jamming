#define NDIM 2
#define PI 3.14159265
#define PI2 6.28318531
#define sqrt2 1.41421356
#define sqrt3 1.73205081

#include <vector>
#include <cmath>
#include <ctime>
#include <chrono>
#include <time.h>
#include <iostream>

using namespace std;
using namespace std::chrono;

#include "../classes/Cell.h"
#include "../classes/Box.h"
#include "../classes/Print.h"
#include "../classes/Fluctuations.h"
#include "../classes/Correlations.h"
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>

const string location = "/Users/Daniel1/Desktop/ActiveMatterResearch/jamming-dynamics/";
//const string location = "/home/dmccusker/remote/jamming-dynamics/";

// Mersenne Twister pseudo-random number generator. Cell radii are drawn from a Gaussian
// distribution, while noise values are drawn from a normal distribution.
std::chrono::time_point<std::chrono::high_resolution_clock> t1 = std::chrono::high_resolution_clock::now();
boost::mt19937 gen(std::chrono::duration_cast<std::chrono::nanoseconds>(t1.time_since_epoch()).count());
boost::uniform_real<> unidist(-PI, PI);
boost::normal_distribution<> normdist(0, 1);
boost::variate_generator< boost::mt19937, boost::uniform_real<> > randuni(gen, unidist);
boost::variate_generator< boost::mt19937, boost::normal_distribution<> > randnorm(gen, normdist);

struct Engine
{
    Engine(string, string, long int, long int, double, double, double);
    ~Engine();
    
    int N;                                  // Number of cells
    double CFself;                          // Self-propulsion force
    double CTnoise;                         // Noise parameter
    double dens;                            // Packing fraction
    string run;                             // Run with current set of parameters
    string fullRun;                         // Entire set of runs
    
    //const int trelax = (int)(1000.0/dt);
    //const int tthermalize = 1e6;
    const int trelax = 1000;
    const int tthermalize = 1000;
    
    long int totalSteps;
    long int countdown;
    long int t;
    double dt;                              // Time step, in units of cell-cell repulsion time
    long int resetCounter;                  // Records the number of times the Verlet skin list is refreshed
    const int nSkip = 100;                  // Print data every nSkip steps
    const int film  = nSkip*100;                    // Make a video of the last time steps
   
    const int timeAvg = 10;                 // Number of instances to average correlation functions
    const int tCorrelation = 100;           // Number of time steps of auto-correlation function
    const int cutoff = 10;                  // Cutoff distance for spatial correlation functions
    
    void start();
    void topology();
    void initCells();
    void assignCellsToGrid();
    void buildVerletLists();
    void relax();
    bool newSkinList();
    void calculate_next_positions();
    void neighborInteractions();
    void calculateCOM();
    void saveOldPositions();
    double calculateOrderParameter();
    vector<double> calculateSystemOrientation();
    void print_video(Print&);
    double delta_norm(double);
    
    double COM[NDIM];                       // Current position of center of mass, with no PBC
    double COM_old[NDIM];                   // Stores old center of mass value for Verlet list skin refresh
    double orderAvg;
    double order2Avg;
    double order4Avg;
    double binder;
    double variance;
    
    vector<Cell> cell;                      // "cell" is vector of "Cell" type objects
    vector<Box> grid;                       // Stores topology of simulation area
    vector<vector<int>> boxPairs;           // List of box pairs that are separated by less than a correlation cut-off
    
    double L;                               // Length of the simulation area
    double Lover2;                          // Read: "L-over-two", so we don't have to calculate L/2 every time we need it
    double lp;                              // Length of one box in the grid
    int b;                                  // Number of boxes in one dimension
    int nbox;                               // Total number of boxes
    const int boxnb = pow(3,NDIM)           // Number of boxes neighboring each other: 9 in 2D, 27 in 3D
    double rn;                              // Radius that defines particle's interaction neighborhood
    double rs;                              // Verlet skin radius
    double rn2;                             // Square distances, to avoid calculating too many square roots
    double rs2;
    
    double MSD(int t);
    
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
    
    t = 0;
    resetCounter = 0;
    
    rn = 2.5;               // Neighborhood radius is fixed at 2.5.
    rs = 1.5*rn;		    // Optimize skin radius choice for the number of particles and speed.
    rn2 = rn*rn;
    rs2 = rs*rs;
    
    orderAvg = 0.0;
    order2Avg = 0.0;
    order4Avg = 0.0;
    binder = 0.0;
    variance = 0.0;
    
    for(int k=0; k<NDIM; k++){
        COM[k] = 0.0;
        COM_old[k] = 0.0;
    }
}

Engine::~Engine()
// Destructor
{
    for (int i=0; i<N; i++) cell[i].VerletList.clear();
    for (int j=0; j<nbox; j++) grid[j].CellList.clear();
}

void Engine::start(){
    
    int tCorrCounter = 0;
    
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    
    initCells();
    topology();
    
    Print print(location, fullRun, run, N);
    Fluctuations fluct(L, totalSteps, nSkip, dens);
    Correlations corr(L, dens, cutoff, tCorrelation, N, CFself);
    
    assignCellsToGrid();
    buildVerletLists();
    
    relax();
    
    while(countdown != 0){
        
        calculate_next_positions();
        
        // Print data every nSkip steps.
        
        if(t%nSkip == 0){
            
            fluct.measureFluctuations(cell, COM, print);
            
            double order = calculateOrderParameter();
            vector<double> orientation = calculateSystemOrientation();
            
            double order2 = order*order;
            orderAvg+=order;
            order2Avg+=order2;
            order4Avg+=order2*order2;
            
            print.print_COM(t, COM);
            print.print_order(t, order);
            print.print_orientation(t, orientation);
            print.print_MSD(t, MSD(t));
            
            //if(countdown<film) print_video(print);
        }
        
        // Calculate static correlation functions and initialize autocorrelation function.
        
        if( t%(totalSteps/timeAvg) == 0 && t!=0 )
        {
            assignCellsToGrid();
            buildVerletLists();
            
            corr.orientation0 = calculateSystemOrientation;
            corr.spatialCorrelations(boxPairs, grid, cell);
            corr.velDist(cell);
            
            tCorrCounter = 0;
        }
        
        // Calculate autocorrelation function.
        
        if( tCorrCounter < tCorrelation )
        {
            corr.autocorrelation( tCorrCounter, calculateSystemOrientation() );
            tCorrCounter++;
        }
        
        t++;
        countdown--;
    }
    
    double nsteps = (double)totalSteps/(double)nSkip;
    orderAvg  /= nsteps;
    order2Avg /= nsteps;
    order4Avg /= nsteps;
    
    binder = 1.0 - order4Avg/(3.0*order2Avg*order2Avg);
    variance = order2Avg - orderAvg*orderAvg;
    
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<seconds>( t2 - t1 ).count();

    corr.printCorrelations(timeAvg, print);
    print.print_summary(run, N, L, t, 1./dt, CFself, CTnoise, dens, duration, resetCounter,
                        corr.correlationValues[0], binder, orderAvg, variance);
}

void Engine::initCells(){
    
    double volume = 0;
    
    // Create vector of cell objects and assign their radii.
    
    for(int i=0; i<N; i++){
        Cell *newCell = new Cell(L, dt);
        cell.push_back(*newCell);
        cell[i].index = i;
        
        double cellRad = 1. + randnorm()/10;
        cell[i].R  = cellRad;
        volume += pow(cellRad,NDIM);
    }
    
    // Simulation area/volume depends on cell sizes.
    
    if(NDIM==2) L = cbrt( 4.0 * PI * volume / (3.0*dens) );
    if(NDIM==3) L = sqrt( PI * volume / dens );
    Lover2 = L/2.0;
    
    // Set up initial cell configuration in an hexagonal lattice.
    
    int rootN = 0;
    if(NDIM==2) rootN = sqrt(N);
    if(NDIM==3) rootN = cbrt(N);
    double spacing = L/rootN;
    
    for (int i=0; i<N; i++) {
        int j = i/rootN;
        cell[i].pos[0] = -Lover2 + spacing*(i%cbrtN) + randnorm()/3.;
        cell[i].pos[1] = -Lover2 + spacing*j + randnorm()/3.;
        if( i%2 == 0 ) { cell[i].pos[1] += 1.0; }
        
        cell[i].phi = randuni();
        cell[i].theta = PI/2.0;
        
        if(NDIM==3) {
            int k = i/(rootN*rootN);
            cell[i].pos[1] = -Lover2 + spacing*j + randnorm()/3. - L*k ;
            cell[i].pos[2] = -Lover2 + spacing*k + randnorm()/3.;
            if( k%2 == 0 ) { cell[i].pos[0] += 1.0; }
            else           { cell[i].pos[1] += 1.0; }
            cell[i].theta = (randuni() + PI)/2.0;
        }
    }
}

void Engine::topology()
// lp is ~at least~ the assigned neighbor region diameter. It can be a little bit bigger such that
// we have an integer number of equally-sized boxes.
// Dynamics are incorrect if b<2 (small number of particles).
{
    lp = 2*rn;
    b = static_cast<int>(floor(L/lp));
    nbox = pow(b,NDIM);
    lp = L/floor(L/lp);
    
    for (int k=0; k<nbox; k++) {
        Box *newBox = new Box;
        grid.push_back(*newBox);
        grid[k].serial_index = k;
    }
    
    // Label box indices, and find each box's neighbors. Including itself, each box has 9 neighbors
    // in 2D and 27 in 3D. Then, build a list of boxes that are separated by less than the cutoff
    // radius. Include the entire diagonal length of the boxes, not just their center-center
    // distance, so add sqrt2(3)*lp to the cutoff distance.
    
    boxPairs.reserve(b*(b+1)/2);
    
    if(NDIM==2){
        for(int i=0; i<b; i++){
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
        for(int k=0; k<b2; k++){
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
        for (int p=0; p<b2; p++){
            for(int q=p; q<b2; q++){
                double deltax = delta_norm(grid[p].center[0]-grid[q].center[0]);
                double deltay = delta_norm(grid[p].center[1]-grid[q].center[1]);
                double boxDist2 = deltax*deltax+deltay*deltay;
                double boxCutoff = cutoff+(sqrt2*lp);
                if (boxDist2 < boxCutoff*boxCutoff){
                    vector<int> temp;
                    temp.assign(2,0);
                    temp[0] = grid[p].serial_index;
                    temp[1] = grid[q].serial_index;
                    boxPairs.push_back(temp);
                }
            }
        }
    }
    
    if(NDIM==3){
        for(int i=0; i<b; i++){
            for(int j=0; j<b; j++){
                for(int k=0; k<b; k++){
                    int p = i+(j*b)+(k*b*b);
                    grid[p].vector_index[0] = i;
                    grid[p].vector_index[1] = j;
                    grid[p].vector_index[2] = k;

                    grid[p].min[0] = -Lover2 + i*lp;
                    grid[p].min[1] = -Lover2 + j*lp;
                    grid[p].min[2] = -Lover2 + k*lp;
                    grid[p].max[0] = grid[p].min[0] + lp;
                    grid[p].max[1] = grid[p].min[1] + lp;
                    grid[p].max[2] = grid[p].min[2] + lp;
                    
                    for (int m=0; m<NDIM; m++) {
                        grid[p].center[m] = (grid[p].min[m] + grid[p].max[m]) / 2.;
                    }
                }
            }
        }
        for(int p=0; p<nbox; p++){
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++){
                    for (int k=0; k<3; k++){
                        int a1 = i;
                        int a2 = j;
                        int a3 = k;
                        
                        if (grid[p].vector_index[0] == 0   && i==0) a1=a1+b;
                        if (grid[p].vector_index[0] == b-1 && i==2) a1=a1-b;
                        if (grid[p].vector_index[1] == 0   && j==0) a2=a2+b;
                        if (grid[p].vector_index[1] == b-1 && j==2) a2=a2-b;
                        if (grid[p].vector_index[2] == 0   && k==0) a3=a3+b;
                        if (grid[p].vector_index[2] == b-1 && k==2) a3=a3-b;
                        
                        grid[p].neighbors[i+(j*3)+(k*9)] = p + (a1-1) + (a2-1)*b + (a3-1)*b*b;
                    }
                }
            }
        }
        for (int p=0; p<nbox; p++){
            for(int q=p; q<nbox; q++){
                double boxDist2 = 0.0;
                for (int k=0; k<NDIM; k++){
                    double dr = delta_norm(grid[p].center[k]-grid[q].center[k]);
                    boxDist2+=dr*dr;
                }
                double boxCutoff = cutoff+(sqrt3*lp);
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
}

void Engine::relax()
// Relax the system as passive particles to allow many rearrangements.
// Then, allow to thermalize, slowly increasing activity to final value.
{
    double CFself_old = CFself;
    
    CFself = 0;
    
    for(int t_=0; t_<trelax; t_++) { calculate_next_positions(); }
    
    for(int t_=0; t_<tthermalize; t_++) {
        CFself = CFself_old - (tthermalize - t_)*CFself_old/tthermalize;
        calculate_next_positions();
    }
    
    CFself = CFself_old;
    
    // Translate to center of mass frame.
    // Start cells at their current location in the grid after relaxing and thermalizing.
    
    calculateCOM();
    
    for(int i=0; i<N; i++){
        for(int k=0; k<NDIM; k++) {
            cell[i].pos[k] = delta_norm(cell[i].pos[k]-COM[k]);
            cell[i].pos_real[k] = cell[i].pos[k];
            cell[i].pos0[k] = cell[i].pos[k];
        }
    }
    
    saveOldPositions();
    resetCounter = 0;
}

void Engine::saveOldPositions()
{
    for(int k=0; k<NDIM; k++) {
        COM_old[k] = COM[k];
        for(int i=0; i<N; i++){
            cell[i].pos_old[k]  = cell[i].pos[k];
        }
    }
}

void Engine::calculateCOM(){
    
    for(int k=0; k<NDIM; k++) COM[k] = 0.0;

    for(int i=0; i<N; i++){
        for(int k=0; k<NDIM; k++){
            COM[k] += cell[i].pos_real[k];
        }
    }
    
    for(int k=0; k<NDIM; k++) COM[k] /= N;
    
}

void Engine::assignCellsToGrid()
// A cell can be, at most, a distance of lp*√3/2 (3D) or lp*√2/2 (2D) from its closest box center.
// Starting at this distance, search through all boxes to find which box the cell is in.
{
    for (int j=0; j<nbox; j++) grid[j].CellList.clear();
    
    for (int i=0; i<N; i++)
    {
        double r2 = lp*lp*0.25*NDIM;
        for (int j=0; j<nbox; j++)
        {
            double d2 = 0;
            for (int k=0; k<NDIM; k++)
            {
                double dr = cell[i].pos[k] - grid[j].center[k];
                d2 += dr*dr;
            }
            if(d2 < r2) { r2 = d2; cell[i].box = j; }
        }
        grid[cell[i].box].CellList.push_back(i);
    }
}

void Engine::buildVerletLists()
{
    for(int i=0; i<N; i++)
    {
        cell[i].VerletList.clear();
        for(int m=0; m<nboxnb; m++)
        {
            int p = grid[cell[i].box].neighbors[m];
            int max = grid[p].CellList.size();
            for(int k=0; k<max; k++)
            {
                int j = grid[p].CellList[k];
                if(j > i)
                {
                    double d2 = 0;
                    for (int l=0; l<NDIM; l++) {
                        double dr = delta_norm(cell[j].pos[l] - cell[i].pos[l]);
                        d2 += dr*dr;
                    }
                    if( d2 < rs2 )
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
// Compare the two largest displacements to see if a skin refresh is required
// Refresh if any particle may have entered any other particle's neighborhood
{
    calculateCOM();
    
    bool refresh = false;
    double largest2 = 0.;
    double second2 = 0.;
    for(int i=0; i<N; i++){
        double d2 = 0;
        for (int k=0; k<NDIM; k++) {
            double dr = delta_norm(cell[i].pos[k] - cell[i].pos_old[k] - COM[k] + COM_old[k]);
            d2 += dr*dr;
        }
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
// *** Most physics happens here *** //
// Calculate spring repulsion force and neighbor orientational interactions.
{
    for(int i=0; i<N; i++){
        int max = cell[i].VerletList.size();
        for(int k=0; k<max; k++){                               // Check each cell's Verlet list for neighbors
            int j = cell[i].VerletList[k];
            if(j > i){                                          // Symmetry reduces calculations by half
                double dr[NDIM] = {0.0};
                double d2 = 0.0;
                
                for (int l=0; l<NDIM; l++){
                    dr[l] = delta_norm(cell[j].pos[l]-cell[i].pos[l]);
                    d2 += dr[l]*dr[l];
                }
                
                if(d2 < rn2){                                   // They're neighbors
                    double sumR = cell[i].R + cell[j].R;
                    
                    if( d2 < sumR*sumR ){                       // They also overlap
                        double overlap = sumR / sqrt(d2) - 1;
                        
                        for (int l=0; l<NDIM; l++)
                        {                                       // Spring repulsion force, watch the sign
                            cell[i].F[l] -= overlap*dr[l];
                            cell[j].F[l] += overlap*dr[l];
                        }
                        
                        if(countdown <= film && t%nSkip == 0){
                            cell[i].over -= 240*abs(overlap);
                            cell[j].over -= 240*abs(overlap);
                        }
                    }
                
                    cell[i].cosp_new += cell[j].cosp;            // Add up orientations of neighbors
                    cell[i].sinp_new += cell[j].sinp;
                    cell[j].cosp_new += cell[i].cosp;
                    cell[j].sinp_new += cell[i].sinp;
                    
                    if(NDIM==3){
                        cell[i].cost_new += cell[j].cost;
                        cell[i].sint_new += cell[j].sint;
                        cell[j].cost_new += cell[i].cost;
                        cell[j].sint_new += cell[i].sint;
                    }
                }
            }
        }
        if(NDIM==2){
            cell[i].phi_new   = atan2(cell[i].sinp_new, cell[i].cosp_new) + CTnoise*randuni();
        }
        if(NDIM==3){
            double delta = CTnoise*randuni();
            double alpha = randuni() + PI;
            cell[i].theta_new = atan2(cell[i].sint_new, cell[i].cost_new) + delta*sin(alpha);
            cell[i].phi_new   = atan2(cell[i].sinp_new, cell[i].cosp_new) + delta*cos(alpha);
        }
    }
}

double Engine::MSD(int t){
    double MSD = 0.0;
    for (int i=0; i<N; i++) {
        for (int k=0; k<NDIM; k++) {
            double dr = cell[i].pos_real[k] - cell[i].pos0[k] - COM[k];
            MSD += dr*dr;
        }
    }
    return MSD/N;
}

double Engine::calculateOrderParameter()
{
    double x = 0.0, y = 0.0, z = 0.0;
    
    for (int i=0; i<N; i++) {
        x += cell[i].cosp*cell[i].sint;
        y += cell[i].sinp*cell[i].sint;
        z += cell[i].cost;
    }
   
    return sqrt(x*x+y*y+z*z)/(double)N;
}

vector<double> Engine::calculateSystemOrientation()
{
    vector<double> orientation(NDIM,0);
    
    for (int i=0; i<N; i++) {
        orientation[0] += cell[i].cosp*cell[i].sint;
        orientation[1] += cell[i].sinp*cell[i].sint;
        if(NDIM==3)  orientation[2] += cell[i].cost;
    }
    
    for (int k=0; k<NDIM; k++ ) orientation[k] /= (double)N;
    
    return orientation;
}

void Engine::print_video(Print &print)
{
    int k=0;
    for(int i=0; i<N; i++){
        print.print_Ovito(k, N, i, cell[i].R, cell[i].over,
                          cell[i].pos[0], cell[i].pos[1], cell[i].pos[2],
                          cell[i].vel[0], cell[i].vel[1], cell[i].vel[2]);
        cell[i].over = 240;
        k=1;
    }
}

void Engine::calculate_next_positions()
{
    if( newSkinList() ) {
        assignCellsToGrid();
        buildVerletLists();
    }
    
    neighborInteractions();
    
    for(int i=0; i<N; i++) cell[i].update();
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
    
    int npoints = 0;
    
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
        npoints = steps/engine.nSkip;
    }
    
    int check = system(("python "+location+"code/plot/singleRunPlots.py "+ID+ " " +dir+" "+std::to_string(npoints)+" "+std::to_string(NDIM)).c_str());
    if( check != 0 ) cout << "An error occurred while plotting (one of) the graphs\n";
    
    return 0;
}
