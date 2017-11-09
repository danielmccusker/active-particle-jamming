using namespace std;

struct Box
{
    Box();
    
    int serial_index;
    int vector_index[2];
    double xmin, xmax, ymin, ymax;
    double center[2];
    int neighbors[9];
    
    vector<int> CellList;
};

Box::Box(){
    serial_index = -1;
    vector_index[0] = -1;
    vector_index[1] = -1;
   
    CellList.reserve(50);
}

