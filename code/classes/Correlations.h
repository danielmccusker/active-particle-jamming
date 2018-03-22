
// ********************************************************************************
// **** Methods for calculating correlation functions and velocity distribution ***
// ********************************************************************************

// *** Gets NDIM from definition in main ***

struct Correlations
{
    Correlations(double, double, double, double, long int, double);
    
    void spatialCorrelations(vector<vector<int>>&, vector<Box>&, vector<Cell>&);
    void autocorrelation(int, vector<double>& );
    void velDist(vector<Cell>&);
    
    void printCorrelations(int, Print&);
    double delta_norm(double);
    
    double L, Lover2, dens;
    long int N;
    
    int correlation_time;  // Time step cut-off for autocorrelation
    double cutoff;
    
    double dr_c;
    double dr_p;
    double dv;
    int np;
    int nc;
    int noBins;
    double norm;
    
    vector<double> orientation0;
    vector<double> correlationValues;
    vector<double> correlationValuesNorm;
    vector<double> counts;
    vector<double> pairCorrelationValues;
    vector<double> autocorrelationValues;
    vector<double> velocityDistributionValues;
    
};

Correlations::Correlations(double L_, double dens_, double cut, double time, long int N_, double CFself_)
{
    N = N_;
    L = L_;
    Lover2 = L/2.0;
    dens = dens_;
    
    cutoff = cut;
    correlation_time = time;
    
    dr_c = 2.0;
    dr_p = 0.01;
    dv = CFself_/50.0;
    np = (int)ceil(cutoff/dr_p);
    nc = (int)ceil(cutoff/dr_c);
    noBins = 100;

    if(NDIM == 2) norm = 2*L*L/(2*PI*dr_p*(double)(N*N));
    if(NDIM == 3) norm = 2*L*L*L/(4*PI*dr_p*(double)(N*N));
    
    orientation0.assign(NDIM,0);
    
    correlationValues.assign(nc, 0);
    correlationValuesNorm.assign(nc,0);
    pairCorrelationValues.assign(np, 0);
    autocorrelationValues.assign(correlation_time,0);
    velocityDistributionValues.assign(noBins,0);
}

void Correlations::spatialCorrelations
    (vector<vector<int>> &boxPairs, vector<Box> &grid, vector<Cell> &cell)
// We normalize the velocity correlations by the number of counts in the bin size. The pair
// correlation normalization is geometric and depends on the system dimension.
{
    vector<double> pairTemp(np,0.0);
    vector<double> corrTemp(nc,0.0);
    vector<double> counts(nc,0.0);
    
    vector<double> tiavg(nc,0.0);
    vector<double> tjavg(nc,0.0);
    vector<double> piavg(nc,0.0);
    vector<double> pjavg(nc,0.0);

    for (int a=0; a<boxPairs.size(); a++)
    {
        int p = boxPairs[a][0];
        int q = boxPairs[a][1];
        int maxp = grid[p].CellList.size();
        int maxq = grid[q].CellList.size();
        
        for (int m=0; m<maxp; m++)
        {
            for (int n=0; n<maxq; n++)
            {
                int i = grid[p].CellList[m];
                int j = grid[q].CellList[n];
                
                if( p!=q || (p==q && j>i) ) // Avoid double-counting pairs in the same box
                {
                    double r = 0.0;
                    for(int k=0; k<NDIM; k++){
                        double dk = delta_norm(cell[j].pos[k]-cell[i].pos[k]);
                        r += dk*dk;
                    }
                    
                    r = sqrt(r);
                    
                    int binp = (int)floor(r/dr_p);
                    int binc = (int)floor(r/dr_c);
                    
                    if(binp < np)   // Exclude any pairs beyond cutoff, normalize
                    {
                        if(NDIM == 2) pairTemp[binp] += 1.0/r;
                        if(NDIM == 3) pairTemp[binp] += 1.0/(r*r);
                    }
                    
                    if(binc < nc)
                    {
                        double pi = cell[i].phi;
                        double pj = cell[j].phi;
                        
                        if(NDIM == 2)
                        {
                            piavg[binc] += pi;
                            pjavg[binc] += pj;
                            corrTemp[binc] += cos(pi-pj);
                        }
                        if(NDIM == 3)
                        {
                            double ti = cell[i].theta;
                            double tj = cell[j].theta;
                            
                            piavg[binc] += pi;
                            pjavg[binc] += pj;
                            tiavg[binc] += ti;
                            tjavg[binc] += tj;
                            corrTemp[binc] += ( cos(ti)*cos(tj) + cos(pi-pj)*sin(ti)*sin(tj) ) ;
                            //cout << ( cos(ti)*cos(tj) + cos(pi-pj)*sin(ti)*sin(tj) ) << endl;
                        }
                        counts[binc] += 1.0;
                    }
                }
            }
        }
    }
    
    for(int k=0; k<nc; k++)
    {
        corrTemp[k] /= counts[k];
        //cout << k << " " << corrTemp[k] << endl;
        correlationValues[k] += corrTemp[k];
        
        piavg[k] /= counts[k];
        pjavg[k] /= counts[k];
        tiavg[k] /= counts[k];
        tjavg[k] /= counts[k];
        
        //cout << k << " " << piavg[k] << " " << tiavg[k] << " " << pjavg[k] << " " << tjavg[k] << endl;
        
        double offset = 0.0;
        if (NDIM==2)
            offset = cos(piavg[k]-pjavg[k]);
        if (NDIM==3)
            offset = cos(tiavg[k])*cos(tjavg[k]) + cos(piavg[k]-pjavg[k])*sin(tiavg[k])*sin(tjavg[k]);
        //cout << offset << endl;
        correlationValuesNorm[k] += ( corrTemp[k] - offset );
    }
    
    for(int k=0; k<np; k++) { pairTemp[k] *= norm;  pairCorrelationValues[k]+=pairTemp[k]; }
}

void Correlations::autocorrelation(int t, vector<double>&orientation)
{
    double value = 0.0;
    for(int k=0; k<NDIM; k++) { value += orientation0[k]*orientation[k]; }
    autocorrelationValues[t] += value;
}

void Correlations::velDist(vector<Cell> &cell)
{
    for(int i=0; i<N; i++){
        double v = cell[i].get_speed();
        int binv = (int)floor(v/dv);
        if(binv < noBins) velocityDistributionValues[binv] += 1.0/(double)N;
    }
}

void Correlations::printCorrelations(int timeAvg, Print &print)
{
    for(int k=0; k<correlationValues.size(); k++)
        print.print_corr(dr_c*(k+1), correlationValues[k]/timeAvg);
    
    for(int k=0; k<correlationValuesNorm.size(); k++)
        print.print_corrNorm(dr_c*(k+1), correlationValuesNorm[k]/timeAvg);
    
    for(int k=0; k<pairCorrelationValues.size(); k++)
        print.print_pairCorr(0.01*k, pairCorrelationValues[k]/timeAvg);
    
    for(int t=0; t<autocorrelationValues.size(); t++)
        print.print_autoCorr(t, autocorrelationValues[t]/timeAvg);
    
    for(int k=0; k<velocityDistributionValues.size(); k++)
        print.print_velDist(k*dv, velocityDistributionValues[k]/timeAvg);
}

double Correlations::delta_norm(double delta)
// Subtracts multiples of the box size to account for periodic boundary conditions
{
    int k=-1;
    if(delta < -Lover2) k=1;
    while(delta < -Lover2 || delta >= Lover2) delta += k*L;
    
    return delta;
}
