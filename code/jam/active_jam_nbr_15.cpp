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
#include "../classes/Cell.h"
#include "../print/print_04.h"
#include "../classes/Box.h"

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
    void buildVerletLists();
    void start();
    void defineCells();
    void relax();
    bool newSkinList();
    void assignCellsToGrid();
    void neighborInteractions();
    void find_nbr();
    double delta_norm(double);
    double cc_overlap(double R1, double R2, double r);
    double totalAreaIntersection(double R);
    void pairCorr();
    void graphics();
    void GNFinit();
    void densityFluctuations();
    
    double MSD();
    double MSDerror();
    double VAF();
    
    void save_config();
    void save_config2();
    
    // Variables
    Print print;                            // Print class has methods for printing all kinds of data
    string path;                            // Change this in the Print class for local or remote runs
    string run;                             // Run with current set of parameters
    string fullRun;                         // Entire set of runs for a given particle number
    
    int film;                               // Make a video of the last time steps
    long int totalSteps;
    long int countdown;
    long int t;                             // Current time
    long int resetCounter;                  // Records the number of times the Verlet skin list is refreshed
    double dt;                              // Time step, in units of cell-cell repulsion time

    int N;                                  // Number of cells
    vector<Cell> cell;                      // "cell" is vector of "Cell" type objects
    vector<Box> grid;                       // Stores topology of simulation area
    
    double L;                               // Length of the simulation area
    double Lover2;                          // Read: "L-over-two", so we don't have to calculate L/2 every time we need it
    double lp;                              // Length of a box
    int b;                                  // Number of boxes in one direction
    int b2;                                 // Total number of boxes
    double rn;                              // Radius that defines particle's interaction neighborhood
    double rs;                              // Verlet skin radius
    double rn2;                             // Square distances, to avoid calculating too many square roots
    double rs2;
    
    double CFself;                          // Self-propulsion force
    double CTnoise;                         // Noise parameter
    double dens;                            // "Packing fraction" / ~average~ density (recall: giant density fluctuations)
    
    double xavg0, yavg0;                    // Initial position of center of mass
    double vcx0, vcy0;                      // Initial velocity of center of mass
    double xavg, yavg;                      // Current position of center of mass, with no PBC
    double xavgold, yavgold;                // Stores values of average positions for Verlet list skin refresh
                                            //The four above are measured in terms of "real" positions
    double vcx, vcy;                        // Current velocity of center of mass
    double vg;                              // Average particle speed
    double ux, uy;                          // Average x- and y-components of particle orientation
        //normalize ux and uy
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

    film                 = 256;
    t          = 0;
    resetCounter         = 0;
}

Engine::~Engine(){
// Destructor
    
    const int gmax = 60001;
    const int vmax = 150;
    int g[gmax] = {0};
    int veldist[vmax] = {0};
    double dr = 0.0001;
    double dv = CFself/100;
    
    for( int i=0; i<N; ++i)
    {
        double v = sqrt(cell[i].velx*cell[i].velx + cell[i].vely*cell[i].vely);
        int binv = v/dv;
        if(binv >= vmax)
        {
           // cout << binv << "=149\n";
            binv = 149;
        }
        veldist[binv] += 1;
        
    }
    
    //for( int k=0; k<vmax; ++k) print.print_v(k, k*dv, veldist[k] / (double)N);
    
    for (int i=0; i<N; i++) {
        cell[i].clearVerletList();
    }
    for (int j=0; j<b2; j++) {
        grid[j].clear_cell_list();
    }
    
//    for(long int t=0; t<(timeCounter - ts0)/screenshotInterval+1; ++t)
//    {
//        long int norm = N*MSDcounter[t];
//        double msdave = MSD[t]/norm;
//        print.print_MSD(t*screenshotInterval, msdave, sqrt(MSDerr[t]/norm - msdave*msdave), (vaf[t]/norm)/2);
//    }
    
    ofstream overWrite;
    overWrite.open((path+"/"+run+"/dat/config_2.dat").c_str());
    overWrite << resetCounter;
    overWrite.close();
}


void Engine::start(){
    
    //high_resolution_clock::time_point t1 = high_resolution_clock::now();               // Start time
    
    path = print.init(fullRun, run, N);
    
    defineCells();
    defineGrid();
    assignCellsToGrid();
    buildVerletLists();
    GNFinit();
    //relax();
    
    while(countdown != 0){

        if( newSkinList() ){
            assignCellsToGrid();
            buildVerletLists();
        }
        
        neighborInteractions();
        print.print_noCells(N);
        xavg = 0;
        yavg = 0;
        for(int i=0; i<N; i++){
            cell[i].update(dt, Lover2, L);
            xavg+=cell[i].xreal;
            yavg+=cell[i].yreal;
            double velmag = sqrt(cell[i].velx*cell[i].velx + cell[i].vely*cell[i].vely);
            double veldir = atan2(cell[i].vely, cell[i].velx);
            print.print_Ovito(i, cell[i].R, cell[i].posx, cell[i].posy, cell[i].psi, 50*velmag, veldir, cell[i].over);
        }
        xavg = xavg/N;
        yavg = yavg/N;
        print.print_MSD(t, MSD());
        //graphics();
        //densityFluctuations();
        t++;
        countdown--;
        
    }
    
    pairCorr();
    
   // high_resolution_clock::time_point t2 = high_resolution_clock::now();                //end time
    
    //auto duration = duration_cast<seconds>( t2 - t1 ).count();                          //elapsed time
    
    print.print_summary(run, N, t, 1./dt, CFself, CTnoise, dens);   //print run summary
    
//    
//    double velmag = sqrt(cell[i].velx*cell[i].velx + cell[i].vely*cell[i].vely);
//    double veldir = atan2(cell[i].vely, cell[i].velx);
//    print.print_vid(cell[i].posx, cell[i].posy, cell[i].psi, 50*velmag, veldir, cell[i].over);
//    print.print_Ovito(i, cell[i].R, cell[i].posx, cell[i].posy, cell[i].psi, 50*velmag, veldir, cell[i].over);
//    print.print_positions(cell[i].posx, cell[i].posy);
//    cell[i].over = 240;
}

void Engine::defineCells(){
    
    double area = 0;
    
    for(int i=0; i<N; i++){                             // Create vector of cells and assign their radii
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
    //vcx0 = 0;
    //vcy0 = 0;
    
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
        //vcx0 += cell[i].velx;
        //vcy0 += cell[i].vely;
    }
    
    xavg0 = xavg0/N;
    yavg0 = yavg0/N;
    //vcx0 = vcx0/N;
    //vcy0 = vcy0/N;
    
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

void Engine::assignCellsToGrid(){                           // Could maybe optimize to not search every box every time
    for (int j=0; j<b2; j++) grid[j].clear_cell_list();
    for (int i=0; i<N; i++){
        for (int j=0; j<b2; j++){
            double R2 = lp*lp/2;          // Circle of radius lp*√2/2 around each box center touches the box's corners
            double dx = cell[i].posx - grid[j].center[0];
            double dy = cell[i].posy - grid[j].center[1];
            double d2 = dx*dx+dy*dy;
            if(d2 < R2)  { R2 = d2; cell[i].box = j; }      // Find box to which the cell belongs
        }
        grid[cell[i].box].add_cell(cell[i]);                // Add this cell to the box's list
    }
}

void Engine::buildVerletLists(){
    for(int i=0; i<N; i++){
        // check again to use symmetry
        for(int m=0; m<9; m++){
            int p = grid[cell[i].box].neighbors[m];         // Get the indices of the 9 boxes to search
            grid[p].check_box(&cell[i], rs2, L, Lover2);      // Add to cell's Verlet list
        }
    }
}

bool Engine::newSkinList(){
// Compare the two largest displacements to see if a skin refresh is required
    
    bool refresh= false;
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
        refresh = true;
        resetCounter++;
        for(int i=0; i<N; i++){
            cell[i].xold = cell[i].xreal;
            cell[i].yold = cell[i].yreal;
            xavgold += cell[i].xreal;
            yavgold += cell[i].yreal;
            cell[i].clearVerletList();
        }
        xavgold /= N;
        yavgold /= N;
    }
    return refresh;
}

void Engine::neighborInteractions(){
// Orientational interaction between neighbors and soft repulsion force between overlapping neighbors
    
// ~Most physics happens here~
    
    for(int i=0; i<N; i++){
        double noise = CTnoise*randuni();
        cell[i].checkVerletList(rn2, L, Lover2, noise);
    }
}

void Engine::relax(){
// Relax the system for 1,000,000 steps, slowly decreasing the activity to the final value
    
    int trelax = 1e3;
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
        
        neighborInteractions();
        
        for(int i=0; i<N; i++){
            cell[i].update(dt, Lover2, L);
            cell[i].psi = randuni();
        }
    }
    
    for(int i=0; i<N; i++){
        cell[i].xreal = cell[i].posx;
        cell[i].yreal = cell[i].posy;
    }
    
    t = 0;
    resetCounter = 0;
    CFself = CFself_old;
    CTnoise = CTnoise_old;
    cout << "stop relaxation " << endl;
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

double Engine::MSDerror(){
    double MSDerr = 0.0;
    for (int i=0; i<N; i++) {
        double dx = (cell[i].xreal - cell[i].x0 - xavg + xavg0);
        double dy = (cell[i].yreal - cell[i].y0 - yavg + yavg0);
        MSDerr += (dx*dx+dy*dy)*(dx*dx+dy*dy);
    }
    return MSDerr;
}

double Engine::VAF(){
    double VAF = 0.0;
    for (int i=0; i<N; i++) {
        double num = (cell[i].velx-vcx)*cell[i].vx0 + (cell[i].vely-vcy)*cell[i].vy0;
        double norm = cell[i].vx0*cell[i].vx0 + cell[i].vy0*cell[i].vy0;
        VAF += num/norm;
    }
    VAF /= N;
    return VAF;
}

double Engine::delta_norm(double delta){
// Subtracts multiples of the boxsize until -L/2 < delta <= L/2

    int k=-1;
    if(delta < -Lover2) k=1;
    while(delta < -Lover2 || delta >= Lover2) delta += k*L;
    
    return delta;
}

void Engine::pairCorr()
//calculates the pair correlation function at a certain time step.
{
    double g = 0.0; //value of the correlation function
    double r = 0; //current radius
    double dr = L/500;  //step size is 1/500 of box size
    double C = L*L/(2*PI); //geometric normalization factor that depends on the system dimension
    
    double drij = 0.0; //distance between a pair
    double dx = 0.0;
    double dy = 0.0;
    double sum = 0.0;
    
    int k = 1;
    
    while(r < Lover2){
        sum = 0.0;
        for (int i=0; i < N; i++){ //iterate through pairs
            for(int j=i+1; j < N; j++){
                dx = delta_norm(cell[i].posx-cell[j].posx);
                dy = delta_norm(cell[i].posy-cell[j].posy);
                drij = sqrt(dx*dx + dy*dy);
                
                if(abs(drij-r) < dr) sum = sum + 1.0; //add if a pair separation is in the interval [r, r+dr]
            }
        }
        sum = sum / N;
        g = (C/r)*sum; //surface area of shell, rad squared if in 3D
        print.print_g(k, r, g);
        r = r + dr;
        k++;
    }
}

void Engine::GNFinit(){ // Initialize GNF measurements
    numberOfGNFpoints = 10;
    R = Lover2;
    GNFinterval = (double)totalSteps/(double)numberOfGNFpoints;
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
    //cout << "total intersecting area: " << Atot << endl;
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

void Engine::graphics()
// Creates the output files and resets output parameters.
{
//
//
//		print.print_data(timeCounter, xavg, yavg, vg, sqrt(ux*ux+uy*uy)/N, atan2(uy, ux), nbrdist2,
//                        cell[0].posx, cell[0].posy, cell[N/3].posx, cell[N/3].posy, cell[2*N/3].posx, cell[2*N/3].posy);
//		if(timeCounter%65536==0)
//		{
//			save_config();
//            if (timeCounter%ts0 == 0)
//            {
//                save_config2();
//            }
//		}
//
//
//	}
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
    
    config << t << "\n";
    
    for(int i=0; i<N; ++i)
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
    
    config << t << "\n";
    
    for(int i=0; i<N; ++i)
    {
        config << setprecision(8) << cell[i].posx << "\t" << cell[i].posy << cell[i].velx << "\t" << cell[i].vely << "\t" << cell[i].R << "\t" << cell[i].psi << "\n";
    }
    
    config.close();
    ++fileNumber;
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
    
    //int check = system(("/Users/Daniel1/Desktop/ActiveMatterResearch/jamming-dynamics/code/plot/plot_jam_act_03.gnu "+ID+ " " +dir+" "+std::to_string(steps)).c_str());
    //if( check != 0 )
      //cout << "An error occurred while plotting (one of) the graphs\n";
    cout << "at the end " << endl;
    return 0;
}
