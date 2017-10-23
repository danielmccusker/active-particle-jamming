#define PI 3.14159265

#include <vector>

using namespace std;

struct Cell
{
    Cell();
    
    double delta_norm(double delta, double L, double Lover2);
    void update(double, double, double);
    void addToVerletList(Cell c);
    void clearVerletList();
    void printVerletList();
    void checkVerletList(double, double, double, double);
    
    double R;                       // Cell radius
    double FselfR;                  // Self-propulsion parameter
    
    int over;    					// Overlap that cell has with neighbours: 240(blue, no overlap) to 0(red, big overlap) [Hue]
    int index;                      // Cell number
    int box;                        // Box number to which the cell belongs
    int VerletListSize;
    
    double posx, posy;              // Position in the grid
    double xreal, yreal;            // Real position if it weren't for periodic boundary conditions
    double x0, y0;                  // Initial positions
    double xold, yold;              // Old position to determine displacement during a skin interval
    double psi;                     // Orientation
    double cosp, sinp;              // x-y projection of the orientation
    double velx, vely;              // Velocity
    double vx0, vy0;                // Initial velocity
    double Fx, Fy;                  // Force
    
    double cosp_new, sinp_new;      // For updating
    double psi_new;
    
    Cell *nextBoxMember, *nextVerletMember;
    Cell *VerletFirst, *VerletLast, *temp;
    
};

Cell::Cell(){
    box = -1;
    index = -1;
    
    posx = -100;
    posy = -100;
    xold = -100;
    yold = -100;
    xreal = -100;
    yreal = -100;
    x0 = -100;
    y0 = -100;
    psi = 0;
    psi_new = 0;
    cosp = 0;
    sinp = 0;
    cosp_new = 0;
    sinp_new = 0;
    
    velx = 0;
    vely = 0;
    Fx  = 0;
    Fy  = 0;
    over = 240;
    
    nextVerletMember = NULL;
    nextBoxMember = NULL;
    VerletFirst = NULL;
    VerletLast = NULL;
    temp = NULL;
    VerletListSize = 0;
};

void Cell::update(double dt, double Lover2, double L){
// Update particle positions and orientations from the equations of motion
// F_res = \zeta R v ,          with \zeta = 1
// New orientation is set by direction of the average of the neighbours and itself + a noise term

    velx = Fx/R;
    vely = Fy/R;
    
    double dx = velx * dt;
    double dy = vely * dt;

    posx += dx;
    posy += dy;
    psi = psi_new;
    xreal += dx;
    yreal += dy;

    cosp = cos(psi);
    sinp = sin(psi);

    if(posx >= Lover2)      posx -= L;
    else if(posx < -Lover2) posx += L;
    if(posy >= Lover2)      posy -= L;
    else if(posy < -Lover2) posy += L;
    if(psi >= PI)           psi -= 2*PI;
    else if(psi < -PI)      psi += 2*PI;

    Fx  = 0;
    Fy  = 0;
    
    cosp_new = cosp; // The average direction of particles in r_n also includes itself
    sinp_new = sinp;
}

void Cell::addToVerletList(Cell c){
    if(VerletFirst != NULL){
        temp = new Cell(c);
        temp->nextVerletMember = VerletFirst;
        VerletFirst = temp;
        VerletListSize++;
    } else {
        VerletLast = new Cell(c);
        VerletLast->nextVerletMember = NULL;
        VerletFirst = VerletLast;
        VerletListSize++;
    }
    //cout << "adding cell " << c.index << " " << " to the Verlet list of cell " << index << endl;
}

void Cell::clearVerletList(){
    temp = VerletFirst;
    while(temp != NULL){
        temp = temp->nextVerletMember;
        delete VerletFirst;
        VerletFirst = temp;
        VerletListSize--;
    }
}

void Cell::printVerletList(){
    temp = VerletFirst;
    cout << VerletListSize << " cells in the Verlet list: ";
    while(temp != NULL){
        cout << temp->index << " ";
        temp = temp->nextVerletMember;
    }
    cout << endl;
}

void Cell::checkVerletList(double rn2, double L, double Lover2, double noise){
    temp = VerletFirst;
    while(temp != NULL){
        //if((temp->index)>index){                                // Symmetry reduces calculations by half
            double dx = delta_norm((temp->posx)-posx,L,Lover2);
            double dy = delta_norm((temp->posy)-posy,L,Lover2);
            double d2 = dx*dx+dy*dy;
            if(d2 < rn2){                                       // They're neighbors
                double sumR = R + temp->R;
                if(d2 < (sumR*sumR)){                           // They also overlap
                    double overlap = sumR / sqrt(d2) - 1;
                    Fx          -= overlap*dx;  // Sign of overlap and dx is really important
                    Fy          -= overlap*dy;
                   // temp->Fx    += overlap*dx;
                   // temp->Fy    += overlap*dy;
                    over -= 240*overlap;
                }
                cosp_new += temp->cosp;
                sinp_new += temp->sinp;
               // temp->cosp_new += cosp;
                //temp->sinp_new += sinp;
            }
        //}
        temp = temp->nextVerletMember;
    }
    Fx += cosp*FselfR;                                //check self propulsion and radius
    Fy += sinp*FselfR;
    psi_new = atan2(sinp_new, cosp_new) + noise;
}

double Cell::delta_norm(double delta, double L, double Lover2)
// Subtracts multiples of the boxsize until -L/2 < delta <= L/2
{
    int k=-1;
    if(delta < -Lover2) k=1;
    while(delta < -Lover2 || delta >= Lover2) delta += k*L;
    
    return delta;
}

