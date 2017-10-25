
// N circular particles are put together in a 2D box with periodic boundary conditions. We will study
// for which values of the self-propulsion force, orientational noise, and density the system is
// jammed. We say a system is jammed when there are no internal rearrangements, that is, the set of
// neighbors does not change. In such a case, each particle is "caged" by its neighbors. The system
// can therefore be jammed even with a global translation.

// The simulation uses a modified Verlet list method. A Verlet list keeps track of every particle
// in a "skin radius" around each particle, which are considered "possible" neighbors. Neighbors are
// particles separated by less than 2.8 times the average particle radius. The Verlet lists are
// refreshed when the particles have moved far enough that they may have entered the neighbor radius
// from outside the skin radius. When refreshing the Verlet lists, rather than searching globally
// for all pairs to see if they are within a skin radius, we assign each cell to the simulation area
// partition ("box"), and search all pairs inside this box and the eight surrounding boxes. To determine
// when to refresh the list, we measure displacements relative to the center of mass. Therefore, we
// calculate the system's center of mass at every time step.

// Be sure to change the directory in the gnuplot file and the print class based on whether running
// locally or remotely on a computing cluster

#include <vector>
#include <cmath>
#include <ctime>
#include <chrono>

#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include "../classes/Cell.h"
#include "../classes/Box.h"
#include "../classes/Print.h"

#define PI 3.14159265

//Mersenne Twister pseudo-random number generator
boost::mt19937 gen(time(0));
boost::uniform_real<> unidist(-PI, PI);
boost::normal_distribution<> normdist(0, 1);
boost::variate_generator< boost::mt19937, boost::uniform_real<> > randuni(gen, unidist);
boost::variate_generator< boost::mt19937, boost::normal_distribution<> > randnorm(gen, normdist);

using namespace std;
using namespace std::chrono;

struct Engine
{
    Engine(string, string, long int, long int, int, double, double, double);
    ~Engine();
    
// Functions
    void defineGrid();
    void defineCells();
    void GNFinit();
    void start();
    void relax();
    
    void assignCellsToGrid(int,double);
    void buildVerletLists();
    double newSkinList();
    void neighborInteractions();
    void calculateCenterOfMass();
    
// Statistical physics functions
    void densityFluctuations();
    void pairCorr();
    double MSD();
    double VAF();
    void velDist();
    
    double delta_norm(double);
    double cc_overlap(double R1, double R2, double r);
    double totalAreaIntersection(double R);
    
// Variables
    Print print;                            // Print class has methods for printing data
    string path;                            // Change this in the Print class for local or remote runs
    string run;                             // Run with current set of parameters
    string fullRun;                         // Entire set of runs for a given particle number
    
    int film;                               // Make a video of the last time steps
    long int totalSteps;
    long int countdown;
    long int t;                             // Current time
    int nSkip;                              // Print data every nSkip steps
    long int resetCounter;                  // Records the number of times the Verlet skin list is refreshed
    double dt;                              // Time step, in units of cell-cell repulsion time

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
    double vcx0, vcy0;                      // Initial velocity of center of mass
    double xavg, yavg;                      // Current position of center of mass, with no PBC
    double xavgold, yavgold;                // Stores values of average positions for Verlet list skin refresh
    
    double vcx, vcy;                        // Current velocity of center of mass
    double vavg;                            // Average particle speed
    double ux, uy;                          // Average x- and y-components of particle orientation
    double nbrdist2;                        // Average squared distance to a neighbor (small for crystalized with holes?)
    
    int rmsCounter;                         //counts time steps for calculating rms density fluctuations
    double rmsValue;
    double R;                               //radius of GNF intersection circle
    double numberOfGNFpoints;               //how many divisions of total run time for each density fluctuation measurement
    double GNFinterval;
};

Engine::Engine(string dir, string ID, long int n, long int steps, int stepsPerTime, double l_s, double l_n, double rho){
// Constructor
    
    fullRun     = dir;
    run         = ID;
    N           = n;
    totalSteps  = steps;
    countdown   = steps+1;
    dt          = 1./stepsPerTime;
    CFself      = l_s;
    CTnoise     = l_n;
    dens        = rho;

    nSkip = 100;
    film  = 256*nSkip;
    t = 0;
    resetCounter = 0;
}

Engine::~Engine(){
// Destructor
    
    for (int i=0; i<N;  i++){
        cell[i].VerletList.clear();
    }
    for (int j=0; j<b2; j++){
        grid[j].CellList.clear();
    }

}

void Engine::start(){
    
// Start time
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    
    path = print.init(fullRun, run, N);
    defineCells();
    defineGrid();
    assignCellsToGrid(0,0);
    buildVerletLists();
    GNFinit();
    
    relax();
    
    while(countdown != 0){                      // Time loop
        
        double refresh = newSkinList();
        if( refresh > 0 ){
            assignCellsToGrid(1, refresh);
            buildVerletLists();
        }
        
        neighborInteractions();
        
        for(int i=0; i<N; i++) cell[i].update(dt, Lover2, L);      // Integrate equations of motion
            
        calculateCenterOfMass();
        
        if(t%nSkip == 0){                                          // Print data every nSkip steps
            
            ux = 0;
            uy = 0;
            vcx = 0;
            vcy = 0;
            for(int i=0; i<N; i++){
                ux += cell[i].cosp;
                uy += cell[i].sinp;
                vcx += cell[i].velx;
                vcy += cell[i].vely;
            }
            ux = ux/N;
            uy = uy/N;
            vcx = vcx/N;
            vcy = vcy/N;
            
            // Time, center of mass position, order parameter, system orientation, three cell trajectories
            print.print_data(t, xavg, yavg, sqrt(ux*ux+uy*uy), atan2(uy, ux),
                             cell[0].posx, cell[0].posy, cell[N/3].posx, cell[N/3].posy, cell[2*N/3].posx, cell[2*N/3].posy);
            print.print_VAF(t, VAF());
            print.print_MSD(t, MSD());
            
            densityFluctuations();
            
            if(countdown<film){                 // Print video for the last part of the run
                print.print_noCells(N);
                for(int i=0; i<N; i++){
                    double velmag = sqrt(cell[i].velx*cell[i].velx + cell[i].vely*cell[i].vely);
                    double veldir = atan2(cell[i].vely, cell[i].velx);
                    print.print_vid(cell[i].posx, cell[i].posy, cell[i].psi, 50*velmag, veldir, cell[i].over);
                    print.print_Ovito(i, cell[i].R, cell[i].posx, cell[i].posy, cell[i].psi, velmag, veldir, cell[i].over);
                    cell[i].over = 240;
                }
            }
        }
        
        t++;
        countdown--;
        //cout << countdown << endl;
    }
    
    //velDist();
    //pairCorr();

// End time
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
// Elapsed time
    auto duration = duration_cast<seconds>( t2 - t1 ).count();
    
    print.print_summary(run, N, t, 1./dt, CFself, CTnoise, dens, duration);
    
}

void Engine::relax(){
// Relax the system for 1,000,000 steps, slowly decreasing the activity to the final value
// Also get initial velocities
    
    int trelax = 1e5;
    double CFrelax = 0.05;
    double CFself_old = CFself;                         // Store parameter values
    double CTnoise_old = CTnoise;
    t = 1;                                              // Prevent the program from producing output
    CTnoise = 0;
    
    for(int t_=0; t_<trelax; t_++){
        CFself = CFself_old + ((CFrelax - CFself_old)*(trelax - t_))/trelax;
        
        double refresh = newSkinList();
        if( refresh > 0 ){
            assignCellsToGrid(1, refresh);
            buildVerletLists();
        }
        
        neighborInteractions();
        
        for(int i=0; i<N; i++){
            cell[i].update(dt, Lover2, L);
            cell[i].psi = randuni();
        }
        
        calculateCenterOfMass();
        
    }
    
    xavg0 = 0;
    yavg0 = 0;
    xavgold = 0;
    yavgold = 0;
    vcx0 = 0;
    vcy0 = 0;
    
    for(int i=0; i<N; i++){
        cell[i].xreal = cell[i].posx;                   // Start cells at their current location in the grid
        cell[i].yreal = cell[i].posy;
        cell[i].x0    = cell[i].posx;
        cell[i].y0    = cell[i].posy;
        xavg0 += cell[i].xreal;
        yavg0 += cell[i].yreal;
        
        cell[i].vx0 = cell[i].velx;                     // Initial velocities
        cell[i].vy0 = cell[i].vely;
        vcx0 += cell[i].vx0;
        vcy0 += cell[i].vy0;
    }
    
    xavg0 = xavg0/N;                                    // Reset center of mass
    yavg0 = yavg0/N;
    xavg    = xavg0;
    yavg    = yavg0;
    xavgold = xavg0;
    yavgold = yavg0;
    vcx0 = vcx0/N;                                      // Initial velocity of center of mass
    vcy0 = vcy0/N;
    
    t = 0;
    resetCounter = 0;
    CFself = CFself_old;
    CTnoise = CTnoise_old;
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
        print.print_size(cellRad);
    }
    
    L = sqrt(PI * area / dens);                         // L depends on the given density and calculated area
    Lover2 = L/2.;
    print.print_size(L);
    
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
        
        cell[i].FselfR = CFself*cell[i].R;

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
    
    rn = 2.8;
    rs = 1.5*rn;
    rn2 = rn*rn;
    rs2 = rs*rs;
    
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

void Engine::assignCellsToGrid(int k, double refresh){
// At the beginning of the simulation, for every particle, we search every box to find which box
// the cell is in. Afterwards, we search only the neighboring boxes of the cell's current box to see
// if it moved boxes.  The particles are moving slowly enough that they are not travelling more than
// one box length per skin refresh.
    
    if(k>0 && lp>refresh){
        for (int i=0; i<N; i++){
            double r2 = lp*lp/2;        // Circle of radius lp*√2/2 around each box center touches the box's corners
            for (int m=0; m<9; m++){
                int j = grid[cell[i].box].neighbors[m];
                double dx = cell[i].posx - grid[j].center[0];
                double dy = cell[i].posy - grid[j].center[1];
                double d2 = dx*dx+dy*dy;
                if(d2 < r2)  { r2 = d2; cell[i].box = j; }      // Find box to which the cell belongs
            }
            grid[cell[i].box].CellList.push_back(i);            // Add this cell to the box's list
        }
    } else if(k==0){
        for (int i=0; i<N; i++){
            double r2 = lp*lp/2;
            for (int j=0; j<b2; j++){
                double dx = cell[i].posx - grid[j].center[0];
                double dy = cell[i].posy - grid[j].center[1];
                double d2 = dx*dx+dy*dy;
                if(d2 < r2)  { r2 = d2; cell[i].box = j; }
            }
            grid[cell[i].box].CellList.push_back(i);
        }
    } else {
        cout << "The cells are moving too fast or something is wrong with the dynamics. exit 100"
             << endl;
        exit(100);
    }
}

void Engine::buildVerletLists(){
    for(int i=0; i<N; i++){
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

double Engine::newSkinList(){
// Compare the two largest displacements to see if a skin refresh is required
// Refresh if any particle may have entered any other particle's neighborhood
    
    double refresh = -1.0;
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
        xavgold = 0;
        yavgold = 0;
        refresh = sqrt(largest2)+sqrt(second2);
        resetCounter++;
        for(int i=0; i<N; i++){
            cell[i].xold = cell[i].xreal;
            cell[i].yold = cell[i].yreal;
            xavgold += cell[i].xreal;
            yavgold += cell[i].yreal;
            cell[i].VerletList.clear();
        }
        xavgold /= N;
        yavgold /= N;
        
        for (int j=0; j<b2; j++)
            grid[j].CellList.clear();
    }
    return refresh;
}

void Engine::neighborInteractions(){
// Orientational interaction between neighbors and soft repulsion force between overlapping neighbors
    
// ~Most physics happens here~
    
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
        cell[i].Fx += cell[i].cosp*cell[i].FselfR;                  // Self-propulsion force
        cell[i].Fy += cell[i].sinp*cell[i].FselfR;
        cell[i].psi_new = atan2(cell[i].sinp_new, cell[i].cosp_new) + CTnoise*randuni();
    }
}

double Engine::MSD(){
    double MSD = 0.0;
    for (int i=0; i<N; i++) {
        double dx = (cell[i].xreal - cell[i].x0 - xavg + xavg0);
        double dy = (cell[i].yreal - cell[i].y0 - yavg + yavg0);
        MSD += dx*dx+dy*dy;
    }
    MSD /= N;
    return MSD;
}

double Engine::VAF(){
    double VAF = 0.0;
    for (int i=0; i<N; i++) {
        double num = (cell[i].velx-vcx)*(cell[i].vx0-vcx0) + (cell[i].vely-vcy)*(cell[i].vy0-vcy0);
        double norm = cell[i].vx0*cell[i].vx0 + cell[i].vy0*cell[i].vy0;
        VAF += num/norm;
    }
    VAF /= N;
    return VAF;
}

void Engine::velDist(){
// Histogram; velocity distribution graphed from 0 to 2 times average velocity
    
    double vmax = 2*CFself;
    const int noBins = 200;
    int dist[noBins] = {0};
    double dv = vmax/(double)noBins;

    for(int i=0; i<N; i++){
        double v = sqrt(cell[i].velx*cell[i].velx + cell[i].vely*cell[i].vely);
        int binv = v/dv;
        dist[binv] += 1;
    }

    for(int k=0; k<noBins; k++)
        print.print_v(k*dv, dist[k] / (double)N);
    
}

double Engine::delta_norm(double delta){
// Subtracts multiples of the boxsize until -L/2 < delta <= L/2 to account for PBC

    int k=-1;
    if(delta < -Lover2) k=1;
    while(delta < -Lover2 || delta >= Lover2) delta += k*L;
    
    return delta;
}

void Engine::pairCorr(){
// Calculates the pair correlation function at a certain time step.

    double g = 0.0;                             // Value of the correlation function
    double r = 0.1;                             // Current radius
    double dr = 0.1;                            // Step size is 1/10 of average cell radius
    double C = L*L/(2*PI);                      // Geometric normalization that depends on the system dimension
    
    while(r < Lover2){
        double sum = 0.0;
        for (int i=0; i < N; i++){
            for(int j=i+1; j < N; j++){
                double dx = delta_norm(cell[i].posx-cell[j].posx);
                double dy = delta_norm(cell[i].posy-cell[j].posy);
                double d2 = dx*dx+dy*dy;
                if(abs((sqrt(d2)-r)) < dr)
                    sum = sum + 1.0;            // Add if a pair separation is in the interval [r, r+dr)
            }
        }
        sum = sum / N;
        g = (C/r)*sum;                          // Surface area of shell, rad squared if in 3D
        print.print_g(r, g);
        r = r + dr;
    }
}

void Engine::GNFinit(){
// Initialize GNF measurements
    
    numberOfGNFpoints = 10;
    R = Lover2;
    GNFinterval = (double)totalSteps/(nSkip*(double)numberOfGNFpoints);
    rmsCounter = 0;
    rmsValue = 0;
}

double Engine::cc_overlap(double R1, double R2, double r){
// Formula for calculating area intersection between two circles C1 and C2 with radii R1 and R2
    
    if      ( R1 + R2 <= r )    return 0;               // There is no overlap between C1 and C2
    else if ( R2 >= R1 + r )    return PI*R1*R1;        // Circle C1 lies completely in C2
    else if ( R1 >= R2 + r )    return PI*R2*R2;        // Circle C2 lies completely in C1
    else {
        double R12 = R1*R1;
        double R22 = R2*R2;
        double x = (R12 - R22 + r*r)/(2*r);             // base of triangle at C1
        double theta = acos(x/R1);                      // angle of triangle at C1
        double A = R12*theta - x*R1*sin(theta);
        
        x = r - x;                                      // Other area (at C2-side)
        theta = acos(x/R2);
        A += R22*theta - x*R2*sin(theta);
        
        return A;
    }
}

double Engine::totalAreaIntersection(double R){
// Draw a circle of radius R at the center of the box, and calculate its total overlapping area with all cells
    
    double Atot = 0;
    for (int i=0; i<N; i++) {
        double d = sqrt((cell[i].posx*cell[i].posx) + (cell[i].posy*cell[i].posy));
        double A = cc_overlap(cell[i].R, R, d);
        Atot = Atot+A;
    }
    return Atot;
}

void Engine::densityFluctuations(){
// Print density fluctuations as a function of average density
    
    rmsCounter++;
    if(rmsCounter < GNFinterval) {
        double A = totalAreaIntersection(R);
        double expectedA = dens*PI*R*R;
        rmsValue = rmsValue + (A-expectedA)*(A-expectedA);
    } else {
    // Store value of rms fluctuation and reset for next measurement, decrease circle size
        print.print_GNF(dens*PI*R*R, sqrt(rmsValue/(double)rmsCounter) );
        rmsValue = 0;
        rmsCounter = 0;
        R = R/2;
    }
}

int main(int argc, char *argv[])
{
    string dir = "";
    string ID = "";
    long int n = 0;
    long int steps = 0;
    int stepsPerTime = 0;
    double l_s = 0;
    double l_n = 0;
    double rho = 0;
    
    if(argc != 9){
        cout    << "Incorrect number of arguments. Need: " << endl
                << "- full run ID" << endl
                << "- single run ID" << endl
                << "- number of cells" << endl
                << "- number of steps" << endl
                << "- steps per unit time" << endl
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
        stepsPerTime    = atoi(argv[5]);
        l_s             = atof(argv[6]);
        l_n             = atof(argv[7]);
        rho             = atof(argv[8]);

        Engine engine(dir, ID, n, steps, stepsPerTime, l_s, l_n, rho);
        engine.start();
    }
    
    int check = system(("/Users/Daniel1/Desktop/ActiveMatterResearch/jamming-dynamics/code/plot/plot_jam_act_03.gnu "+ID+ " " +dir+" "+std::to_string(steps)).c_str());
    if( check != 0 ) cout << "An error occurred while plotting (one of) the graphs\n";
    
    return 0;
}
