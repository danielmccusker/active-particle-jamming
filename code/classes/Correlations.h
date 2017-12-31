#include <vector>

struct Correlations
{
    Correlations();
    
    vector<vector<double>> buildFullPairList(double cutoff);
    void velocityCorrelation(vector<vector<double>> list);
    void pairCorrelation(vector<vector<double>> list)
}

vector<vector<double>> Correlations::buildFullPairList(double cutoff, int N){
    vector<vector<double>> list;
    list.reserve( N*(N-1)/2 );
    int listCount = 0;
    
    //    for (int p=0; p<b2; p++)
    //    {
    //        for(int q=0; q<b2; q++)                         // Check cells within the same box
    //        {
    //            double dx = delta_norm(grid[p].center[0]-grid[q].center[0]);
    //            double dy = delta_norm(grid[p].center[1]-grid[q].center[1]);
    //            double boxDist2 = dx*dx+dy*dy;
    //            double boxCutoff = cutoff+(sqrt2*lp);
    //            if (boxDist2 < boxCutoff*boxCutoff)
    //            {
    //                int maxp = grid[p].CellList.size();
    //
    //                for (int m=0; m<maxp; m++)
    //                {
    //                    int maxq = grid[q].CellList.size();
    //
    //                    for (int n=0; n<maxq; n++)
    //                    {
    //                        int i = grid[p].CellList[m];
    //                        int j = grid[q].CellList[n];
    //
    //                        if(j > i)
    //                        {
    //                            double deltax = delta_norm(cell[j].posx-cell[i].posx);
    //                            double deltay = delta_norm(cell[j].posy-cell[i].posy);
    //                            double d2 = deltax*deltax+deltay*deltay;
    //                            vector<double> entry(3);
    //                            entry[0] = (double)i;
    //                            entry[1] = (double)j;
    //                            entry[2] = d2;
    //                            list.push_back(entry);
    //                            listCount++;
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }
    
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            if(j > i)
            {
                double deltax = delta_norm(cell[j].posx-cell[i].posx);
                double deltay = delta_norm(cell[j].posy-cell[i].posy);
                double d2 = deltax*deltax+deltay*deltay;
                vector<double> entry(3);
                entry[0] = (double)i;
                entry[1] = (double)j;
                entry[2] = d2;
                list.push_back(entry);
                listCount++;
            }
        }
    }
    
    cout << "length of list " << listCount << " " << list.size() << endl;
    return list;
}

void Correlations::velocityCorrelation(vector<vector<double>> list)
{
    double r0 = 2.0;                    // Initial radius
    double r = r0;
    double dr = 2.0;
    double velCutoff = Lover2; //120;
    double velCutoff2 = velCutoff*velCutoff;
    int listsize = list.size();
    
    while(r < velCutoff-dr)
    {
        double velVal = 0.0;
        double velValCorrected = 0.0;
        double norm = 1.0/(2.0*PI*r);  // Geometric normalization that depends on the system dimension
        double avgvelxi = 0.0;
        double avgvelxj = 0.0;
        double avgvelyi = 0.0;
        double avgvelyj = 0.0;
        
        for (int k=0; k<listsize; k++)
        {
            double d2 = list[k][2];
            
            if(d2 < velCutoff2 && abs(sqrt(d2)-r) < dr)
            {
                int i = (int)list[k][0];
                int j = (int)list[k][1];
                double vxi = cell[i].velx;
                double vyi = cell[i].vely;
                double vxj = cell[j].velx;
                double vyj = cell[j].vely;
                avgvelxi+=vxi;
                avgvelxj+=vxj;
                avgvelyi+=vyi;
                avgvelyj+=vyj;
                velVal += vxi*vxj+vyi*vyj;
            }
        }
        
        avgvelxi = norm*avgvelxi/(CFself*(double)N);
        avgvelyi = norm*avgvelyi/(CFself*(double)N);
        avgvelxj = norm*avgvelxj/(CFself*(double)N);
        avgvelyj = norm*avgvelyj/(CFself*(double)N);
        velVal = norm*velVal/(CFself*CFself*(double)N);
        velValCorrected = velVal - (avgvelxi*avgvelxj+avgvelyi*avgvelyj);
        
        velocityCorrelationValues.push_back(velVal);
        velocityCorrelationCorrectedValues.push_back(velValCorrected);
        r+=dr;
    }
}

void Correlations::pairCorrelation(vector<vector<double>> list)
{
    double r0 = 1.0;
    double r = r0;
    double dr = 0.1;
    double pairCutoff = 10;
    double pairCutoff2 = pairCutoff*pairCutoff;
    
    while(r < pairCutoff-dr)
    {
        double pairVal = 0.0;
        double norm = 1.0/(2.0*PI*r);  // Geometric normalization that depends on the system dimension
        
        for (int k=0; k<list.size(); k++)
        {
            double d2 = list[k][2];
            
            if(d2 < pairCutoff2 && abs((sqrt(d2)-r)) < dr)
            {
                pairVal+=1.0;
            }
        }
        
        pairVal = norm*pairVal/(dens*(double)N);
        pairCorrelationValues.push_back(pairVal);
        r+=dr;
    }
}
