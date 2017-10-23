// ============================
// Fixed:
// - 
// Added:
// - 02: xreal and yreal for mean squared displacement
// Changed:
// - 01: Update mechanism for orientation
// Removed:
// - 01: Torque, omga and R2
// - 03: Removed x0, x1, ..., 0, y1, ...

#define PI 3.14159265

#include <vector>

using namespace std;

struct Cell
{
    Cell();
    
    void update(double, double, double);

    int over;    					// Overlap that cell has with neighbours: 240(blue, no overlap) to 0(red, big overlap) [Hue]

    double posx, posy, psi;
    double velx, vely;
    double Fx, Fy;
    double cosp_new, sinp_new, psi_new;
    double R;
    double cosp, sinp;
    double xold, yold;              // Old position to determine displacement during interval
    double xreal, yreal;            // Real position if it weren't for periodic boundary conditions

    vector<int> nbrRegion;          // indices of potentially neighbouring cells
};

Cell::Cell()
{
    xold = -10;
    yold = -10;
    velx = 0;
    vely = 0;
    Fx  = 0;
    Fy  = 0;
    over = 240;

    nbrRegion.reserve(15);
}

void Cell::update(double dt, double lover2, double l)
// Update particle positions and orientations from the Equations of Motion
// F_res = \zeta R v ,          with \zeta = 1
// New orientation is set by direction of the average of the neighbours and itself + a noise term
{
    velx = Fx / R;
    vely = Fy / R;
    
    double dx = velx * dt;
    double dy = vely * dt;

    posx += dx;
    posy += dy;
    psi = psi_new;
    xreal += dx;
    yreal += dy;

    cosp = cos(psi);
    sinp = sin(psi);

	if(posx >= lover2) posx -= l;
	else if(posx < -lover2) posx += l;
	if(posy >= lover2) posy -= l;
	else if(posy < -lover2) posy += l;
    if(psi >= PI) psi -= 2*PI;
    else if(psi < -PI) psi += 2*PI;

    Fx  = 0;
    Fy  = 0;
    cosp_new = cosp;            // The average direction of all particles in the neighbourhood also includes itself
    sinp_new = sinp;
}
