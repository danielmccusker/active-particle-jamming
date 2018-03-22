//#include <iostream>
//using std::cout;
//using std::endl;
//using std::cin;

struct Cell
{
    Cell();
    
    void update(double);
    void PBC();
    void periodicAngles();
    double get_speed();
    
    double L, Lover2, dt;
    
    double R;                       // Cell radius
    
    int over;    					// Overlap that cell has with neighbours: 240(blue, no overlap) to 0(red, big overlap) [Hue]
    int index;                      // Cell number
    int box;                        // Box number to which the cell belongs
    
    vector<double> pos;             // Position in the grid
    vector<double> pos_real;        // Real position if it weren't for periodic boundary conditions
    vector<double> pos0;            // Initial positions
    vector<double> pos_old;         // Old position to determine displacement during a skin interval
    vector<double> vel;             // Velocity
    vector<double> F;               // Force
    
    double phi, theta;              // Orientation in x-y plane, angle from z-axis
    double cosp, sinp,
           cost, sint;
    double x_new, y_new, z_new;     // Projections of the cell's orientation for adding
    double phi_new, theta_new;
    
    vector<int> VerletList;
};

Cell::Cell()
{
    L = -1.0;
    Lover2 = -1.0;
    dt = -1.0;
    
    box = -1;
    index = -1;
    over = 0;
    
    pos.assign(NDIM,-100);
    pos_old.assign(NDIM,-100);
    pos0.assign(NDIM,-100);
    pos_real.assign(NDIM,-100);
    vel.assign(NDIM,0);
    F.assign(NDIM,0);
    
    phi = theta = 0.0;
    cosp = sinp = cost = sint = 0.0;
    x_new = y_new = z_new = 0.0;
    phi_new=theta_new=0.0;
    
    if(NDIM==2){
        theta = PI/2.0;
        cost = 0.0;
        sint = 1.0;
        VerletList.reserve(20);
    }
    
    if(NDIM==3)
    {
        VerletList.reserve(40);
    }
};

void Cell::update(double CFself)
// Update particle positions and orientations from the equations of motion.
{
    if(NDIM==2)
    {
        periodicAngles();
        
        cosp = cos(phi);
        sinp = sin(phi);
        
        F[0] += cosp*CFself*R;   // Self-propulsion force
        F[1] += sinp*CFself*R;
        
        x_new = cosp;                // The average direction of particles in the neighborhood
        y_new = sinp;                // also includes itself
    }
    
    if(NDIM==3)
    {
        periodicAngles();
        
        cosp = cos(phi);
        sinp = sin(phi);
        cost = cos(theta);
        sint = sin(theta);
        
        F[0] += sint*cosp*CFself*R;
        F[1] += sint*sinp*CFself*R;
        F[2] += cost*CFself*R;
        
        phi_new = phi;
        theta_new = theta;
        
        x_new = cosp*sint;
        y_new = sinp*sint;
        z_new = cost;
    }
    
    for(int k=0; k<NDIM; k++)
    {
        vel[k]       = F[k] / R;
        double dr    = vel[k]*dt;
        pos[k]      += dr;
        pos_real[k] += dr;
        F[k] = 0;
    }
    
    PBC();
}

void Cell::periodicAngles()
{
    if(phi >= PI)           phi  -= PI2;
    else if(phi < -PI)      phi  += PI2;
    if(theta >= PI)         theta = PI2-theta;
    else if(theta < 0)      theta = -theta;
}

void Cell::PBC()
{
    for(int k=0; k<NDIM; k++)
    {
        if(pos[k] >= Lover2)      pos[k] -= L;
        else if(pos[k] < -Lover2) pos[k] += L;
    }
}

double Cell::get_speed()
{
    double speed = 0.0;
    for(int k=0; k<NDIM; k++) speed+=vel[k]*vel[k];
    return sqrt(speed);
}
