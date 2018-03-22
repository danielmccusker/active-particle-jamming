
// *** One spatial partition of the simulation area ***

struct Box
{
    Box();
    
    int serial_index;
    vector<int> vector_index;
    vector<double> min, max;
    vector<double> center;
    vector<int> neighbors;
    vector<int> CellList;
};

Box::Box(){
    serial_index = -1;
    vector_index.assign(NDIM,-1);
    min.assign(NDIM,0);
    max.assign(NDIM,0);
    center.assign(NDIM,0);
    
    if(NDIM==2) neighbors.assign(9,0);
    if(NDIM==3) neighbors.assign(27,0);

    CellList.reserve(50);
}

