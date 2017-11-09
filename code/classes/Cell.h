#define PI 3.14159265

#include <vector>

using namespace std;

struct Cell
{
    Cell();
    
    void update(double, double, double);
    
    double R;                       // Cell radius
    
    int over;    					// Overlap that cell has with neighbours: 240(blue, no overlap) to 0(red, big overlap) [Hue]
    int index;                      // Cell number
    int box;                        // Box number to which the cell belongs
    
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
    
    vector<int> VerletList;
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
    over = 0;
    
    VerletList.reserve(27);
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
