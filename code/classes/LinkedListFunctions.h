// Not currently implemented in the code
// In case there is a situation where linked lists are preferable to vectors for Verlet list

//*****************
//*** Cell class **
//*****************

void addToVerletList(Cell c);
void clearVerletList();
void printVerletList();
void checkVerletList(double, double, double, double);

Cell *nextBoxMember, *nextVerletMember;
Cell *VerletFirst, *VerletLast, *temp;

nextVerletMember = NULL;
nextBoxMember = NULL;
VerletFirst = NULL;
VerletLast = NULL;
temp = NULL;

void Cell::addToVerletList(Cell c){
    if(VerletFirst != NULL){
        temp = new Cell(c);
        temp->nextVerletMember = VerletFirst;
        VerletFirst = temp;
        VerletListSize++;
    } else {
        VerletLast = new Cell(c);
        VerletLast->nextVerletMember = NULL;
        VerletFirst = VerletLast;
        VerletListSize++;
    }
    //cout << "adding cell " << c.index << " " << " to the Verlet list of cell " << index << endl;
}

void Cell::clearVerletList(){
    temp = VerletFirst;
    while(temp != NULL){
        temp = temp->nextVerletMember;
        delete VerletFirst;
        VerletFirst = temp;
        VerletListSize--;
    }
}

void Cell::printVerletList(){
    temp = VerletFirst;
    cout << VerletListSize << " cells in the Verlet list: ";
    while(temp != NULL){
        cout << temp->index << " ";
        temp = temp->nextVerletMember;
    }
    cout << endl;
}

void Cell::checkVerletList(double rn2, double L, double Lover2, double noise){
    temp = VerletFirst;
    while(temp != NULL){
        //if((temp->index)>index){                                // Symmetry reduces calculations by half
            double dx = delta_norm((temp->posx)-posx,L,Lover2);
            double dy = delta_norm((temp->posy)-posy,L,Lover2);
            double d2 = dx*dx+dy*dy;
            if(d2 < rn2){                                       // They're neighbors
                double sumR = R + temp->R;
                if(d2 < (sumR*sumR)){                           // They also overlap
                    double overlap = sumR / sqrt(d2) - 1;
                    Fx          -= overlap*dx;  // Sign of overlap and dx is really important
                    Fy          -= overlap*dy;
                    //temp->Fx    += overlap*dx;
                    //temp->Fy    += overlap*dy;
                    //over -= 240*overlap;
                }
                cosp_new += temp->cosp;
                sinp_new += temp->sinp;
               // temp->cosp_new += cosp;
               // temp->sinp_new += sinp;
            }
       // }
        temp = temp->nextVerletMember;
    }
    Fx += cosp*FselfR;                                //check self propulsion and radius
    Fy += sinp*FselfR;
    psi_new = atan2(sinp_new, cosp_new) + noise;
}

// ************
// Box class **
// ************

 void add_cell(Cell c);
 void clear_cell_list();
 void print_neighbors();
 void print_cell_list();
 void cell_positions();
 double delta_norm(double delta, double L, double Lover2);
 void check_box(Cell* c, double rs2, double L, double Lover2);

Cell *boxFirst, *boxLast, *temp;

boxFirst = NULL;
boxLast = NULL;
temp = NULL;

void Box::add_cell(Cell c){
// Adds cell to front of linked list

    if(boxFirst != NULL){
        temp = new Cell(c);
        temp->nextBoxMember = boxFirst;
        boxFirst = temp;
        cell_count++;
    } else{
        boxLast = new Cell(c);
        boxLast->nextBoxMember = NULL;
        boxFirst = boxLast;
        cell_count++;
    }
}

void Box::clear_cell_list(){
    temp = boxFirst;
    while(temp != NULL){
        temp = boxFirst->nextBoxMember;
        delete boxFirst;
        boxFirst = temp;
        cell_count--;
    }
}

void Box::print_neighbors(){
    for(int i=0; i<9; i++) cout << neighbors[i] << " ";
    cout << endl;
}

void Box::print_cell_list(){
    temp = boxFirst;
    while(temp != NULL){
        cout << temp->index << " ";
        temp = temp->nextBoxMember;
    }
    cout << endl;
}

void Box::check_box(Cell* c, double rs2, double L, double Lover2){
// Iterate through this box's members to check for cells to add to c's Verlet list

    temp = boxFirst;
    while(temp != NULL){
        //if ( (temp->index) > (c->index) ){                      // Symmetry reduces calculations by half
            double dx = delta_norm(c->posx-(temp->posx),L,Lover2);
            double dy = delta_norm(c->posy-(temp->posy),L,Lover2);
            if(dx*dx+dy*dy < rs2 && (c->index)!=(temp->index)){ // Don't include the cell itself in the list
                c->addToVerletList(*temp);
                //temp->addToVerletList(*c);
            }
        //}
        temp = temp->nextBoxMember;
    }
}

double Box::delta_norm(double delta, double L, double Lover2)
// Subtracts multiples of the boxsize until -L/2 < delta <= L/2
{
    int k=-1;
    if(delta < -Lover2) k=1;
    while(delta < -Lover2 || delta >= Lover2) delta += k*L;

    return delta;
}
