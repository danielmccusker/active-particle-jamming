
// ********************************************************************************
// **** Methods for calculating density fluctuations                            ***
// ********************************************************************************

// *** Gets NDIM from definition in main ***

struct Fluctuations
{
    Fluctuations(double, double, double, double);
    
    double overlap(double, double, double);
    void measureFluctuations(vector<Cell>&, vector<double>&, Print&);
    double delta_norm(double);
    void density_distribution(vector<Cell>&, vector<Box>&);
    void print_density_distribution(int, Print&);
    
    vector<double> distribution;
    
    int nPoints;
    double min_radius;
    double max_radius;
    double current_radius;
    double rad_interval;
    double time_interval;
    int counter;
    double current_value;
    double dens;
    
    double Lover2, L;
};

Fluctuations::Fluctuations(double L_, double totalSteps, double nSkip, double dens_)
{
    L = L_;
    Lover2 = L/2.0;
    dens = dens_;
    
    nPoints = 10;
    min_radius = 3.0;
    max_radius = Lover2;
    current_radius = min_radius;
    rad_interval = pow( max_radius/min_radius, 1./((double)nPoints+1.0) );
    time_interval = (double)totalSteps/((double)nSkip*(double)nPoints);
    counter = 0;
    current_value = 0;
    
    distribution.assign(50,0.0);
};

void Fluctuations::measureFluctuations(vector<Cell> &cell, vector<double> &COM, Print &print)
// Draw a measurement circle around the center of mass and calculate total intersecting volume
// Calculate fluctuation from the expected intersecting volume
{
    double expectedV = 0.0;
    if (NDIM == 2) expectedV = dens*PI*current_radius*current_radius;
    if (NDIM == 3) expectedV = 4.0*dens*PI*current_radius*current_radius*current_radius/3.0;
    int N = cell.size();
    
    if(counter < time_interval) {
        double V = 0.0;
        for (int i=0; i<N; i++)
        {
            double sumR = cell[i].R+current_radius;
            double d2 = 0.0;
            for (int k=0; k<NDIM; k++) {
                double dr = delta_norm(cell[i].pos_real[k]-COM[k]);
                d2 += dr*dr;
            }
            
            if( d2 <= sumR*sumR )
            {   // If the circles have a nonzero overlap
                V += overlap(cell[i].R, current_radius, sqrt(d2));
            }
        }
        current_value+=(V-expectedV)*(V-expectedV);
    }
    else
    {   // Store value, increase measurement region size
        current_value = sqrt(current_value/(double)counter);
        print.print_fluct(current_radius, expectedV, current_value);
        current_value = 0;
        counter = -1;
        current_radius = current_radius*rad_interval;
    }
    counter++;
}

double Fluctuations::overlap(double r, double R, double d)
{
    if (R>=r+d)
    {   // If cell lies completely in measurement circle
        if(NDIM == 2) return PI*r*r;
        if(NDIM == 3) return 4.0*PI*r*r*r/3.0;
    }
    
    else
    {   // If there is a partial intersection
        if(NDIM == 2)
        {
            double R12 = r*r;
            double R22 = R*R;
            double x = (R12 - R22 + d*d)/(2.0*d);
            double theta = acos(x/r);
            double A = R12*theta - x*r*sin(theta);
            x = d - x;
            theta = acos(x/R);
            A += R22*theta - x*R*sin(theta);
            
            return A;
        }
        
        if(NDIM == 3)
        {
            double x2 = (r+R-d)*(r+R-d)/(12.0*d);
            double y2 = d*d + 2.0*d*r - 3.0*r*r + 2.0*d*R + 6*r*R - 3*R*R;
            return PI*x2*y2;
        }
    }
}

void Fluctuations::density_distribution( vector<Cell> &cell, vector<Box> &grid )
// Count cell centers in each spatial partition and plot the frequencies of densities in a histogram
{
    int size = grid.size();
    int delta = 1;
    
    vector<double> counts(size,0.0);
    
    for (int j=0; j<size; j++) {
        double count = grid[j].CellList.size();
        counts[j] = count;
    }
    
    for (int j=0; j<size; j++) {
        int bin = (int)floor(counts[j]/delta);
        if(bin < distribution.size()) distribution[bin]+=1.0;
    }
}

void Fluctuations::print_density_distribution(int timeAvg, Print &print)
{
    for(int k=0; k<distribution.size(); k++)
        print.print_dens(k, distribution[k]/timeAvg);
}

double Fluctuations::delta_norm(double delta)
// Subtracts multiples of the box size to account for periodic boundary conditions
{
    int k=-1;
    if(delta < -Lover2) k=1;
    while(delta < -Lover2 || delta >= Lover2) delta += k*L;
    
    return delta;
}
