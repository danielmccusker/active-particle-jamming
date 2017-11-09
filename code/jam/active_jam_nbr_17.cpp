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
    
    // Functions
    void defineGrid();
    void defineCells();
    void start();
    void relax();
    void GNFinit();
    
    void assignCellsToGrid(int,double);
    void buildVerletLists();
    double newSkinList();
    void neighborInteractions();
    void calculateCenterOfMass();
    void saveOldPositions();
    
    // Statistical physics functions
    void densityFluctuations();
    void pairCorr();
    double MSD();
    double VAF();
    void velDist();
    vector<double> orientationCorrelation();
    
    double delta_norm(double);
    double cc_overlap(double R1, double R2, double r);
    
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
    dt          = 0.1;           // Fix 10 steps per unit time
    CFself      = l_s;
    CTnoise     = l_n;
    dens        = rho;
    
    nSkip = 100;
    film  = 300*nSkip;
    t = 0;
    resetCounter = 0;
    
    rn = 2.8;                       // Neighborhood radius is fixed.
    rs = 1.5*rn;					// Optimize skin radius choice for the number of particles:
    rn2 = rn*rn;                    // More particles -> assignCellsToGrid takes more time.
    rs2 = rs*rs;                    // A larger value of rs then reduces the number of calls
                                    // to assignCellsToGrid.

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
    
    assignCellsToGrid(0,0);
    buildVerletLists();
    
    relax();
    
    while(countdown != 0){                                          // Time loop
        
        double refresh = newSkinList();
        if( refresh > 0 ){
            assignCellsToGrid(0, refresh);
            buildVerletLists();
        }
        
        neighborInteractions();
        
        for(int i=0; i<N; i++) cell[i].update(dt, Lover2, L);      // Integrate equations of motion
        
        calculateCenterOfMass();
        
        if(t%nSkip == 0){                                          // Print data every nSkip steps
            
            densityFluctuations();
            
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
            
            //if(countdown<film){                 // Print video for the last part of the run
            if(t<100000){
                int k=0;
                for(int i=0; i<N; i++){
                    print.print_Ovito(k, N, i, cell[i].posx, cell[i].posy, cell[i].R, cell[i].over);
                    cell[i].over = 240;
                    k=1;
                }
            }
        }
        
        t++;
        countdown--;
    }
    
    //velDist();
    //pairCorr();
    vector<double> o = orientationCorrelation();
    
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<seconds>( t2 - t1 ).count();
    
    print.print_summary(run, N, L, t, 1./dt, CFself, CTnoise, dens, duration, resetCounter, o[0]);
}

void Engine::relax(){
    // Relax the system for 1,000,000 steps, slowly decreasing the activity to the final value
    // Also get initial velocities
    
    int trelax = 1e6;
    double CFrelax = 0.05;
    double CFself_old = CFself;                         // Store parameter values
    double CTnoise_old = CTnoise;
    t = 1;                                              // Prevent the program from producing output
    CTnoise = 0;
    
    for(int t_=0; t_<trelax; t_++){
        CFself = CFself_old + ((CFrelax - CFself_old)*(trelax - t_))/trelax;
        
        double refresh = newSkinList();
        if( refresh > 0 ){
            assignCellsToGrid(0, refresh);
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
        cell[i].xold  = cell[i].posx;
        cell[i].yold  = cell[i].posy;
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

void Engine::assignCellsToGrid(int k, double refresh){
// At the beginning of the simulation, for every particle, we search every box to find which box
// the cell is in. Afterwards, we search only the neighboring boxes of the cell's current box to see
// if it moved boxes.  The particles are moving slowly enough that they are not travelling more than
// one box length per skin refresh.
    for (int j=0; j<b2; j++) grid[j].CellList.clear();
    
    if(k==0){                               // Search all boxes
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
    } else if(k==1 && refresh<lp){          // Search only neighboring boxes of the old box
        for (int i=0; i<N; i++){
            int oldBox = cell[i].box;
            double r2 = lp*lp/2;
            for (int m=0; m<9; m++){
                int j = grid[oldBox].neighbors[m];
                double dx = delta_norm(cell[i].posx - grid[j].center[0]);
                double dy = delta_norm(cell[i].posy - grid[j].center[1]);
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
        resetCounter++;
        saveOldPositions();
        refresh = sqrt(largest2)+sqrt(second2);
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
        cell[i].Fx += cell[i].cosp*cell[i].R*CFself;                  // Self-propulsion force
        cell[i].Fy += cell[i].sinp*cell[i].R*CFself;
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
    const double dv = 0.01;
    int noBins = (int)ceil(1./dv);
    vector<int> dist(noBins);
    
    for(int i=0; i<N; i++){
        double v = sqrt(cell[i].velx*cell[i].velx + cell[i].vely*cell[i].vely);
        int binv = (int)floor(v/dv);
        dist[binv] += 1;
    }
    
    for(int k=0; k<noBins; k++)
        print.print_v(k*dv, dist[k] / (double)N);
    
}

double Engine::delta_norm(double delta){
// Subtracts multiples of the box size to account for PBC
    
    int k=-1;
    if(delta < -Lover2) k=1;
    while(delta < -Lover2 || delta >= Lover2) delta += k*L;
    
    return delta;
}

void Engine::pairCorr(){
// Calculates the pair correlation function at a certain time step.
// Peaks where a pair separation is in the interval [r, r+dr)
    
    double g = 0.0;                             // Value of the correlation function
    double r = 0.1;                             // Current radius
    double dr = 0.1;                            // Step size is 1/10 of average cell radius
    double C = L*L/(2*PI);                      // Geometric normalization that depends on the system dimension
    double cutoff = 20;
    double cutoff2 = cutoff*cutoff;
    
    while(r < cutoff){
        int sum = 0;
        for (int i=0; i < N; i++){
            for(int j=i+1; j < N; j++){
                double dx = delta_norm(cell[i].posx-cell[j].posx);
                double dy = delta_norm(cell[i].posy-cell[j].posy);
                double d2 = dx*dx+dy*dy;
                if(d2 < cutoff2 && abs((sqrt(d2)-r)) < dr)  sum++;
            }
        }
        g = (C/r)*(double)sum / N;                          // Surface area of shell, rad squared if in 3D
        print.print_g(r, g);
        r = r + dr;
    }
}

vector<double> Engine::orientationCorrelation(){
    vector<double> values(2);
    int r = 3;                          // Start with nearest neighbors
    int dr = 2;                         // Bin size is the average diameter
    int max = (int)(Lover2/4.0);
    double C = 1.0/(2.0*PI*double(N));
    double peak = 0.0;
    double currentVal = 1.0;
    int k=0;                            // Store the value of the first iteration
    
    // Measure characteristic decay length, at which the correlation = peak value/e
    while(r < max){
        currentVal = 0.0;
        double avgsini = 0.0;
        double avgsinj = 0.0;
        double avgcosi = 0.0;
        double avgcosj = 0.0;
        for (int i=0; i < N; i++){
            for(int j=i+1; j < N; j++){
                double dx = delta_norm(cell[i].posx-cell[j].posx);
                double dy = delta_norm(cell[i].posy-cell[j].posy);
                double d2 = dx*dx+dy*dy;
                if(abs(sqrt(d2)-r) < dr){
                    currentVal += (cell[i].cosp*cell[j].cosp)+(cell[i].sinp*cell[j].sinp);
                    avgsini += cell[i].sinp;
                    avgsinj += cell[j].sinp;
                    avgcosi += cell[i].cosp;
                    avgcosj += cell[j].cosp;
                }
            }
        }
        double norm = C/(double)r;
        currentVal = norm*( currentVal - norm*(avgsini*avgsinj - avgcosi*avgcosj) );
        print.print_orientationCorr(r, currentVal);
        r = r + dr;
        if(k==0) peak = currentVal;
        k++;
    }
    values[0] = peak;
    values[1] = (r-1)/2;
    return values;
}

void Engine::GNFinit(){
    
    numberOfGNFpoints = 10;
    RadMin = 3;                 // Start with a circle the size of one neighborhood
    RadMax = Lover2/2;          // Don't reach L/2 because then cells leaving the circle can immediately reenter on the other side
    
    Rad = RadMin;
    RadInt = pow(RadMax/RadMin, 1./(numberOfGNFpoints-1));
    timeOneDensityMeasurement = (double)totalSteps/(nSkip*(double)numberOfGNFpoints);
    GNFcounter = 0;
    fluctuationValue = 0;
}

void Engine::densityFluctuations(){
// Draw nine circles in the grid and measure density fluctuations, average
// Circles move with the center of mass
    
    double expectedA = dens*PI*Rad*Rad;
    if(GNFcounter < timeOneDensityMeasurement) {
        double A = 0.0;
        for (int i=0; i<N; i++) {                               // For every cell, calculate its intersection with nine circles
            double dx = cell[i].posx-L*(double)q/3.0;
            double dy = cell[i].posy-L*(double)p/3.0;
            double d2 = dx*dx+dy*dy;
            if((cell[i].R+Rad)*(cell[i].R+Rad) >= d2)   // If the circles have a nonzero overlap
                A += cc_overlap(cell[i].R, Rad, sqrt(d2));
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
    
    //if      ( cellR + measurementR <= r ) return 0;
    // There is no overlap between circles
    //if ( cellR >= measurementR + r ) return PI*measurementR*measurementR;
    // Measurement circle lies completely in cell
    
    if (measurementR>=cellR+r)  return PI*cellR*cellR;    // Cell lies completely in measurement circle
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
    
    //int check = system(("/home/dmccusker/remote/jamming-dynamics/code/plot/plot_jam_act_03.gnu "+ID+ " " +dir+" "+std::to_string(steps)).c_str());
    int check = system(("/Users/Daniel1/Desktop/ActiveMatterResearch/jamming-dynamics/code/plot/plot_jam_act_03.gnu "+ID+ " " +dir+" "+std::to_string(steps)).c_str());
    if( check != 0 ) cout << "An error occurred while plotting (one of) the graphs\n";
    
    return 0;
}
