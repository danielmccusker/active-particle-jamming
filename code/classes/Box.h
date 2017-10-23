using namespace std;

struct Box
{
    Box();
    
   // void add_cell(Cell c);
   // void clear_cell_list();
   // void print_neighbors();
    //void print_cell_list();
   // void cell_positions();
   // double delta_norm(double delta, double L, double Lover2);
   // void check_box(Cell* c, double rs2, double L, double Lover2);
    
    int serial_index;
    int vector_index[2];
    double xmin, xmax, ymin, ymax;
    double center[2];
    int neighbors[9];
    
    //Cell *boxFirst, *boxLast, *temp;
    
    vector<int> CellList;
};

Box::Box(){
    serial_index = -1;
    vector_index[0] = -1;
    vector_index[1] = -1;
    //boxFirst = NULL;
    //boxLast = NULL;
    //temp = NULL;
    CellList.reserve(50);
}

//void Box::add_cell(Cell c){
//// Adds cell to front of linked list
//    
//    if(boxFirst != NULL){
//        temp = new Cell(c);
//        temp->nextBoxMember = boxFirst;
//        boxFirst = temp;
//        cell_count++;
//    } else{
//        boxLast = new Cell(c);
//        boxLast->nextBoxMember = NULL;
//        boxFirst = boxLast;
//        cell_count++;
//    }
//}
//
//void Box::clear_cell_list(){
//    temp = boxFirst;
//    while(temp != NULL){
//        temp = boxFirst->nextBoxMember;
//        delete boxFirst;
//        boxFirst = temp;
//        cell_count--;
//    }
//}
//
//void Box::print_neighbors(){
//    for(int i=0; i<9; i++) cout << neighbors[i] << " ";
//    cout << endl;
//}
//
//void Box::print_cell_list(){
//    temp = boxFirst;
//    while(temp != NULL){
//        cout << temp->index << " ";
//        temp = temp->nextBoxMember;
//    }
//    cout << endl;
//}

//void Box::check_box(Cell* c, double rs2, double L, double Lover2){
//// Iterate through this box's members to check for cells to add to c's Verlet list
//    
//    temp = boxFirst;
//    while(temp != NULL){
//        //if ( (temp->index) > (c->index) ){                      // Symmetry reduces calculations by half
//            double dx = delta_norm(c->posx-(temp->posx),L,Lover2);
//            double dy = delta_norm(c->posy-(temp->posy),L,Lover2);
//            if(dx*dx+dy*dy < rs2 && (c->index)!=(temp->index)){ // Don't include the cell itself in the list
//                c->addToVerletList(*temp);
//                //temp->addToVerletList(*c);
//            }
//        //}
//        temp = temp->nextBoxMember;
//    }
//}
//
//double Box::delta_norm(double delta, double L, double Lover2)
//// Subtracts multiples of the boxsize until -L/2 < delta <= L/2
//{
//    int k=-1;
//    if(delta < -Lover2) k=1;
//    while(delta < -Lover2 || delta >= Lover2) delta += k*L;
//    
//    return delta;
//}
