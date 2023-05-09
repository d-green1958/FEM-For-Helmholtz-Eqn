#include "FiniteElementMethod.h"

using namespace std;
using namespace FEM;

int main()
{
    // define the wavenumebr
    double k = 10;


    // create parameters for the mesh:

    // mesh will have element edges at 0,0.25,1
    vector<double> Edges{0, 0.25, 1};

    // the first element will have 4 nodes, the second will have 2
    vector<int> Nodes_Per_Element{4,2};

    // now initalise the mesh made of Nodes
    Mesh<Node> mesh(Edges, Nodes_Per_Element);

    // lets make some Dirichlet Conditions:

    // pin the 0th node only
    vector<int> Dir_Nodes{0};

    // pin it with the value 2 + 3i
    complex<double> condition(2,3);
    vector<complex<double>> Dir_cond {condition};

    // make a Neumann condition:

    // Neumann condition at final node (node 5 ( 4 in cpp))
    vector<int> Neu_Nodes{4};

    // give the sommerfeld conditions (u'= 0 + 0i + (0 + ik)u)
    complex<double> const_cond(0,0); // no condition on the constants
    complex<double> u_cond(0,k); // u conditions have ik coeffcients

    // ready to build the mesh!
    mesh.build();

    // self explanatory
    mesh.apply_Dirichelet_BCs(Dir_cond,Dir_Nodes);
    mesh.apply_Neumann_BCs({const_cond},{u_cond},Neu_Nodes);
    mesh.assign_eq_nums();

    // give the mesh some initial values
    mesh.apply_homogenous_ICs();
    
    // initialise the matrices and vectors
    mesh.initialise_res_and_jac();

    // display mesh to check things
    mesh.display(true); // true => gives more information about mesh

    // complete an iteration 
    mesh.next_iteration(k);

    // display again to see what happened
    mesh.display(true);
    // looks good to me
    
    return 0;
}