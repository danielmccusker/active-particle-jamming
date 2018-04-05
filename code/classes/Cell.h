//#include <iostream>
//using std::cout;
//using std::endl;
//using std::cin;

struct Cell
{
    Cell();
    
    void update(double&);
    void PBC();
    void periodicAngles();
    double get_speed();
    
    double L, Lover2, dt;
    
    double R;                       // Cell radius
    double Rinv; 					// 1/R
    const double Zinv = (9*PI)/16;  // Proportionality constant for Stoke's law in 3D
    
    int over;    					// Overlap that cell has with neighbours: 240(blue, no overlap) to 0(red, big overlap) [Hue]
    int index;                      // Cell number
    int box;                        // Box number to which the cell belongs
    
    vector<double> x;            	// Position in the grid
    vector<double> x_real;        	// Real position if it weren't for periodic boundary conditions
    vector<double> x0;            	// Initial positions
    vector<double> x_old;         	// Old position to determine displacement during a skin interval
	
    double vx;						// Velocity
    double vy;
    double vz;
	
    double Fx;               		// Force
    double Fy;
    double Fz;
    
    double phi, theta;              // Self-propulsion vector in x-y plane, angle from z-axis
    double cosp, sinp,
           cost, sint;
    double x_new, y_new, z_new;     // Projections of the cell's self-propulsion vector
    
    vector<int> VerletList;
};

Cell::Cell()
{
    L = -1.0;
    Lover2 = -1.0;
    dt = -1.0;
    R = -1.0;
	
    Fx = 0.0;
    Fy = 0.0;
    Fz = 0.0;
    vx = 0.0;
    vy = 0.0;
    vz = 0.0;
    
    box = -1;
    index = -1;
    over = 0;
    
    x.assign(NDIM,-100);
    x_old.assign(NDIM,-100);
    x0.assign(NDIM,-100);
    x_real.assign(NDIM,-100);
    
    phi = theta = 0.0;
    cosp = sinp = cost = sint = 0.0;
    x_new = y_new = z_new = 0.0;
    
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

void Cell::update(double &CFself)
// Update particle positions and orientations from the equations of motion:
// F_i = 6*pi*eta*R_i in 2D
// F_i = (32/3)*eta*R_i in 3D
// Let eta = 1/(6pi) in 2D, then in 3D the proportionality constant is 16/(9*pi)
{
    if(NDIM==2)
    {
        periodicAngles();
        
        cosp = cos(phi);
        sinp = sin(phi);
        
        Fx += cosp*CFself*R;         // Self-propulsion force
        Fy += sinp*CFself*R;
        
        x_new = cosp;                // The average direction of particles in the neighborhood
        y_new = sinp;                // also includes itself
		
        vx = Fx*Rinv;
        vy = Fy*Rinv;
		
        double dx = vx*dt;
        double dy = vy*dt;
		
		x[0] += dx;
		x[1] += dy;
		x_real[0] += dx;
		x_real[1] += dy;
		
		Fx = 0.0;
		Fy = 0.0;
    }
    
    if(NDIM==3)
    {
        periodicAngles();
        
        cosp = cos(phi);
        sinp = sin(phi);
        cost = cos(theta);
        sint = sin(theta);
        
        Fx += sint*cosp*CFself*R;
        Fy += sint*sinp*CFself*R;
        Fz += cost*CFself*R;
        
        x_new = cosp*sint;
        y_new = sinp*sint;
        z_new = cost;
		
        vx = Fx*Zinv*Rinv;
        vy = Fy*Zinv*Rinv;
        vz = Fz*Zinv*Rinv;
		
        double dx = vx*dt;
        double dy = vy*dt;
        double dz = vz*dt;
		
		x[0] += dx;
		x[1] += dy;
		x[2] += dz;
		x_real[0] += dx;
		x_real[1] += dy;
		x_real[2] += dz;
		
		Fx = 0.0;
		Fy = 0.0;
		Fz = 0.0;
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
        if(x[k] >= Lover2)      x[k] -= L;
        else if(x[k] < -Lover2) x[k] += L;
    }
}

double Cell::get_speed()
{
    if(NDIM==2) 		return sqrt(vx*vx+vy*vy);
    else if(NDIM==3) 	return sqrt(vx*vx+vy*vy+vz*vz);
}
