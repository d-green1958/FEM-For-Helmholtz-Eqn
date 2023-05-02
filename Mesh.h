#ifndef MESH
#define MESH
#include <iomanip>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <typeinfo>
#include "LinearSolvers.h"
#include "elements.h"
#include "nodes.h"

// collection of nodes and elements linked together with member functions
// to build the mesh, display information about the mesh, calcualte the 
// results and write the solution to file. NODE_T is type of node that the 
// mesh will use. Currently there is only one node type.
template <class NODE_T>
class Mesh
{
private:
    // #########################################################
    // PRIVATE MEMBER VARIABLES

    // number of nodes and elements in total
    int Num_Nodes;
    int Num_Elements;

    // outer boundary of PML
    double X_PML;

    // edges of the elements in the mesh
    vector<double> Element_Edge_Pos;

    // numbers of nodes per element
    vector<int> Nodes_Per_Element;

    // residulas vector
    vector<complex<double>> Residuals_Vec;

    // jacobian matrix
    vector<vector<complex<double>>> Jacobian_Mat;

    // true if assign_eq_nums() has been used
    bool global_eqn_nums_assigned = false;

    // variables used to override the gauss quadrature process
    bool override_gauss_quadrature = false;
    int n_gauss;

    // is the system tridiagonal? (assume it is)
    bool tridiagonal = true;

public:
    // #########################################################
    // PUBLIC MEMBER VARIABLES

    // pointers to elements
    vector<Element *> List_Of_Elements;

    // pointers to nodes (global coordinates)
    vector<Node *> List_Of_Nodes;

    // #########################################################
    // CONSTRUCTOR FUNCTIONS

    // default constructor
    Mesh() {}

    // Contructor for the mesh class
    // arguements:
    // Element_Edge_Pos_ : the edge positions of each element
    // Nodes_Per_Element_ : the nodes per element
    Mesh(vector<double> Element_Edge_Pos_, vector<int> Nodes_Per_Element_)
    {
        if (Nodes_Per_Element_.size() != Element_Edge_Pos_.size() - 1)
        {
            throw invalid_argument("Incorrect sizes!");
        }

        Nodes_Per_Element = Nodes_Per_Element_;
        Element_Edge_Pos = Element_Edge_Pos_;
        Num_Elements = Element_Edge_Pos.size() - 1;
        Num_Nodes = 1;
        for (int e = 0; e < Num_Elements; e++)
        {
            Num_Nodes += Nodes_Per_Element[e] - 1;
        }
    }

    // #########################################################

    ~Mesh() = default;

    // #########################################################
    // SETUP FUNCTIONS:

    // build function used to create the mesh and link the elements and nodes
    // Arguments:
    // - n : number of nodes per element
    // - Element_Edge_Pos : positions of element edges
    void build()
    {
        // create pointers to nodes
        for (int i = 0; i < Num_Nodes; i++)
        {
            List_Of_Nodes.push_back(new NODE_T);
        }

        // initialise global node number
        int j_global = 0;
        // node spacing inside element
        double h;
        // number of nodes in element e
        int N;
        // LH boundary of element
        double L_edge;
        // RH boundary of element
        double R_edge;

        // create elements of type given by Nodes_Per_Element
        for (int e = 0; e < Num_Elements; e++)
        {
            if (Nodes_Per_Element[e] == 2)
            {
                List_Of_Elements.push_back(new TwoNodeElement);
            }
            else if (Nodes_Per_Element[e] == 3)
            {
                List_Of_Elements.push_back(new ThreeNodeElement);
                tridiagonal = false;
            }
            else if (Nodes_Per_Element[e] == 4)
            {
                List_Of_Elements.push_back(new FourNodeElement);
                tridiagonal = false;
            }
            else
            {
                throw invalid_argument("n Must Be In 1,2,3,4!");
            }

            // use Element_Edge_Pos to define nodal positions
            L_edge = Element_Edge_Pos[e];
            R_edge = Element_Edge_Pos[e + 1];
            h = (R_edge - L_edge) / (Nodes_Per_Element[e] - 1);

            // defining X for each node
            for (int j_local = 0; j_local < Nodes_Per_Element[e] - 1; j_local++)
            {
                List_Of_Nodes[j_global]->X = L_edge + (j_local)*h;
                j_global += 1;
            }
        }

        // assign position of final node (missed out by loop)
        List_Of_Nodes[j_global]->X = R_edge;

        for (int i = 0; i < Num_Nodes; i++)
        {
            List_Of_Nodes[i]->global_eqn_number = 0;
        }

        // assign elements with pointers to nodes
        for (int e = 0; e < Num_Elements; e++)
        {
            List_Of_Elements[e]->Node_Pointer.resize(Nodes_Per_Element[e]);
            for (int j = 0; j < Nodes_Per_Element[e]; j++)
            {
                List_Of_Elements[e]->Node_Pointer[j] = List_Of_Nodes[lookup_node_num(j, e)];
            }
        }
        cout << "MESH LINKED & BUILT! \t D_total = [" << List_Of_Nodes[0]->X << " , "
             << List_Of_Nodes[Num_Nodes - 1]->X << "]" << endl;
    }

    // lookup function used to find global element number from
    // elememt number and local node number
    int lookup_node_num(int j_local, int element_number)
    {
        // return (element_number * (n - 1) + j_local);
        int temp = 0;
        for (int e = 0; e < element_number; e++)
        {
            temp += List_Of_Elements[e]->Node_Pointer.size() - 1;
        }
        return temp + j_local;
    }

    // function to resize the vector of global residuals and jacobian matrix
    void initialise_res_and_jac()
    {
        int Num_Unknowns = get_nou();

        Residuals_Vec.resize(Num_Unknowns);
        Jacobian_Mat.resize(Num_Unknowns);

        for (int i = 0; i < Num_Unknowns; i++)
        {
            Jacobian_Mat[i].resize(Num_Unknowns);
        }
        cout << "JACOBIAN, RESIDUALS INITIALISED" << endl;

        for (int e = 0; e < Num_Elements; e++)
        {
            List_Of_Elements[e]->initialise_element_res_and_jac();
        }
    }

    // function used to pin nodes with dirichelet BCs
    // arguments:
    // - Dir_BCs : the values the nodes are pinned with
    // - Boundary_Nodes : the index of the pinned nodes
    void apply_Dirichelet_BCs(vector<complex<double>> Dir_BCs,
                              vector<int> Boundary_Nodes)
    {
        if (global_eqn_nums_assigned)
        {
            throw logic_error("Eqn Numbers Have Already Been Assigned!");
        }

        for (int i = 0; i < Boundary_Nodes.size(); i++)
        {
            List_Of_Nodes[Boundary_Nodes[i]]->U = Dir_BCs[i];
            List_Of_Nodes[Boundary_Nodes[i]]->global_eqn_number = -1;
        }
    }

    // function to impose neumann (flux conditions):
    // u'(x_N) = (a+bi) + (c+di)u
    // arguments:
    // - Const_Coeffcient : a+bi
    // - U_Coeffcient : c+di
    // - Flux_Nodes : the nodes that the conditions are imposed on
    void apply_Neumann_BCs(complex<double> &Const_Coeffcient,
                           complex<double> &U_Coeffcient, vector<int> &Flux_Nodes)
    {
        for (int i = 0; i < Flux_Nodes.size(); i++) //
        {
            List_Of_Nodes[Flux_Nodes[i]]->flux_node = 1;
            List_Of_Nodes[Flux_Nodes[i]]->Const_Coef = Const_Coeffcient;
            List_Of_Nodes[Flux_Nodes[i]]->U_Coef = U_Coeffcient;
        }
    }

    // Imposes intial conditions (guess) on the unpinned nodes.
    void apply_ICs(const vector<complex<double>> &ICs)
    {
        int unpinned_node = 0;
        for (int i = 0; i < Num_Nodes; i++)
        {
            if (List_Of_Nodes[i]->global_eqn_number != -1)
            {
                List_Of_Nodes[i]->U = ICs[unpinned_node];
                unpinned_node += 1;
            }
        }
        cout << "INITIAL CONDITIONS APPLIED" << endl;
    }

    // apply homogenous ICs to every element that isnt already pinned
    void apply_homogenous_ICs()
    {
        for (int i = 0; i < Num_Nodes; i++)
        {
            // check if pinned
            if (List_Of_Nodes[i]->global_eqn_number != -1)
            {
                // not pinned so set U to 0
                List_Of_Nodes[i]->U = 0;
            }
        }
    }

    // assign the global equation numbers assuming all pinned
    // nodes have already been pinned (Dirichlet conditions already imposed)
    void assign_eq_nums()
    {
        int unpinned_nodes = 0;
        for (int i = 0; i < Num_Nodes; i++)
        {
            if (List_Of_Nodes[i]->global_eqn_number != -1)
            {
                List_Of_Nodes[i]->global_eqn_number = unpinned_nodes;
                unpinned_nodes += 1;
            }
        }
        global_eqn_nums_assigned = true;
    }

    // #########################################################
    // CALCULATION FUNCTIONS:

    // get the number of unkowns (the number of unpinned nodes)
    int get_nou()
    {
        if (global_eqn_nums_assigned)
        {
            int No_Of_Unknowns = 0;
            for (int i = 0; i < Num_Nodes; i++)
            {
                if (List_Of_Nodes[i]->global_eqn_number != -1)
                {
                    No_Of_Unknowns += 1;
                }
            }
            if (No_Of_Unknowns > 8000)
            {
                throw logic_error("Possible Memory Error, Careful!");
            }

            return No_Of_Unknowns;
        }
        else
        {
            throw logic_error("Eqn Numbers Not Assigned!");
        }
    }

    // function used to calculate the next iteration and update the values accordingly
    void next_iteration(double k, int repeats = 1)
    {
        if (!global_eqn_nums_assigned)
        {
            throw logic_error("Eqn Numbers Not Assigned!");
        }

        int eq_num_j;
        int eq_num_k;
        int i_dof;
        int j_dof;

        for (int rep = 0; rep < repeats; rep++)
        {

            // now clear the residuals and jacobian vectors (from previous run)
            for (int i = 0; i < Residuals_Vec.size(); i++)
            {
                Residuals_Vec[i] = 0;
            }

            for (int i = 0; i < Jacobian_Mat.size(); i++)
            {
                for (int j = 0; j < Jacobian_Mat[0].size(); j++)
                {
                    Jacobian_Mat[i][j] = 0;
                }
            }

            for (int e = 0; e < Num_Elements; e++)
            {
                List_Of_Elements[e]->clear_jac_and_res();
            }

            // update the elemental jacobian and residuals for all elements
            for (int e = 0; e < Num_Elements; e++)
            {
                List_Of_Elements[e]->update_residual_and_jacobian(k, X_PML, override_gauss_quadrature, n_gauss);
            }

            cout << "ITERATION COMPLETE" << endl;

            // add contributions to global res and jac
            // loop over elements
            for (int e = 0; e < Num_Elements; e++)
            {
                // now loop over nodes
                i_dof = 0;
                for (int j_local = 0; j_local < List_Of_Elements[e]->Node_Pointer.size(); j_local++)
                {
                    eq_num_j = List_Of_Elements[e]->Node_Pointer[j_local]->global_eqn_number;
                    if (eq_num_j != -1)
                    {
                        // add elemental contribution to residuals vector
                        Residuals_Vec[eq_num_j] += List_Of_Elements[e]->Element_Residuals_Vec[i_dof];

                        if (List_Of_Elements[e]->Node_Pointer[j_local]->flux_node == 1)
                        {
                            Residuals_Vec[eq_num_j] -= List_Of_Elements[e]->Node_Pointer[j_local]->U_Coef *
                                                       List_Of_Elements[e]->Node_Pointer[j_local]->U;

                            Residuals_Vec[eq_num_j] -= List_Of_Elements[e]->Node_Pointer[j_local]->Const_Coef;
                        }

                        j_dof = 0;
                        for (int k_local = 0; k_local < List_Of_Elements[e]->Node_Pointer.size(); k_local++)
                        {
                            eq_num_k = List_Of_Elements[e]->Node_Pointer[k_local]->global_eqn_number;
                            if (eq_num_k != -1)
                            {
                                // add elemental contribution to jacobian matrix
                                Jacobian_Mat[eq_num_j][eq_num_k] +=
                                    List_Of_Elements[e]->Element_Jacobian_Mat[i_dof][j_dof];

                                // check if element is a flux boundary element
                                if (List_Of_Elements[e]->Node_Pointer[k_local]->flux_node == 1 &&
                                    List_Of_Elements[e]->Node_Pointer[j_local]->flux_node == 1)
                                {
                                    Jacobian_Mat[eq_num_j][eq_num_k] -=
                                        List_Of_Elements[e]->Node_Pointer[k_local]->U_Coef;
                                }

                                j_dof += 1;
                            }
                        }
                        i_dof += 1;
                    }
                }
            }

            cout << "CONTRIBUTIONS ADDED" << endl;

            // NOT VERY EFFCIENT CHANGE THIS (BELOW)
            vector<complex<double>> Negative_Residuals_Vec = Residuals_Vec;
            for (int i = 0; i < Residuals_Vec.size(); i++)
            {
                Negative_Residuals_Vec[i] *= -1.0;
            }

            // now use a solver to solve the linear system
            vector<complex<double>> solution;

            if (tridiagonal)
            {
                cout << "TRIDIAGONAL => THOMAAS SOLVER" << endl;
                // exploiting the sparsness of the matrix
                solution = LinearSolvers::Thomas_Solver(Jacobian_Mat, Negative_Residuals_Vec);

                // solution = LU_Solver(Jacobian_Mat, Negative_Residuals_Vec);
            }
            else
            {
                cout << "LU SOLVER" << endl;
                solution = LinearSolvers::LU_Solver(Jacobian_Mat, Negative_Residuals_Vec);
            }

            for (int i = 0; i < Num_Nodes; i++)
            {
                if (List_Of_Nodes[i]->global_eqn_number != -1)
                {
                    List_Of_Nodes[i]->U += solution[List_Of_Nodes[i]->global_eqn_number];
                }
            }

            cout << "VALUES UPDATED" << endl;
        }
    }

    // #########################################################
    // DISPLAY FUNCTIONS:

    // function to display mesh information
    void display(bool FULL_LIST = false)
    {
        cout << "Mesh contains " << Num_Nodes << " nodes." << endl;
        cout << "Mesh contains " << Num_Elements << " elements" << endl;

        if (FULL_LIST == true)
        {
            int node_counter = 0;
            for (int e = 0; e < List_Of_Elements.size(); e++)
            {
                cout << "Element: " << e << endl;
                cout << "#######################" << endl;
                for (int i = 0; i < List_Of_Elements[e]->Node_Pointer.size(); i++)
                {
                    cout << "Node: " << i << " Global Node: " << node_counter
                         << " Global Eq Num: " << List_Of_Elements[e]->Node_Pointer[i]->global_eqn_number
                         << " X: " << List_Of_Elements[e]->Node_Pointer[i]->X << " U: "
                         << List_Of_Elements[e]->Node_Pointer[i]->U << "Flux Node: " << List_Of_Elements[e]->Node_Pointer[i]->flux_node << endl;

                    node_counter += 1;
                }
                cout << "#######################" << endl;
                node_counter -= 1;
            }
        }
    }

    // display list of elements and some info about them
    void display_elements()
    {
        for (int e = 0; e < List_Of_Elements.size(); e++)
        {
            cout << "Element: " << e << " D_e: [" << List_Of_Elements[e]->Node_Pointer[0]->X
                 << " , " << List_Of_Elements[e]->Node_Pointer.back()->X << "] n: "
                 << List_Of_Elements[e]->Node_Pointer.size() << endl;
        }
    }

    // display list of the nodes and some info about them
    void display_nodes()
    {
        for (int i = 0; i < Num_Nodes; i++)
        {
            cout << "Node: " << i << " X: " << List_Of_Nodes[i]->X << " U: "
                 << List_Of_Nodes[i]->U << " Global Eqn Num: " << List_Of_Nodes[i]->global_eqn_number
                 << " Flux Node: " << List_Of_Nodes[i]->flux_node << endl;
        }
    }

    // print residual vector
    void display_residuals()
    {
        cout << "Residuals:" << endl;
        for (int i = 0; i < Residuals_Vec.size(); i++)
        {
            cout << Residuals_Vec[i] << endl;
        }
    }

    // print jacobian
    void display_jacobian()
    {
        cout << "Jacobian:" << endl;
        for (int i = 0; i < Jacobian_Mat.size(); i++)
        {
            for (int j = 0; j < Jacobian_Mat[0].size(); j++)
            {
                cout << Jacobian_Mat[i][j] << " , ";
            }
            cout << endl;
        }
    }

    // #########################################################
    // GET SOLUTION FUNCTIONS:

    // returns interpolated value of u(x) at a global coordiante x
    complex<double> interpolated_u_glb_coor(double x)
    {
        double L_bound;
        double R_bound;
        double local_coord;
        for (int e = 0; e < Num_Elements; e++)
        {
            L_bound = List_Of_Elements[e]->Node_Pointer[0]->X;
            R_bound = List_Of_Elements[e]->Node_Pointer.back()->X;

            if (L_bound <= x && x <= R_bound)
            {
                local_coord = (2 * x - L_bound - R_bound) / (R_bound - L_bound);

                return List_Of_Elements[e]->interpolated_u(local_coord);
            }
        }
        cout << "ERROR: OUT OF MESH!" << endl;
        throw invalid_argument("Attempting To Interpolate Outside Of Mesh!");

        complex<double> empty(0, 0);
        return empty;
    }

    // return value of u'(x) at a global coordinate x
    complex<double> interpolated_der_u_glb_coor(double x)
    {
        double L_bound;
        double R_bound;
        double local_coord;
        for (int e = 0; e < Num_Elements; e++)
        {
            L_bound = List_Of_Elements[e]->Node_Pointer[0]->X;
            R_bound = List_Of_Elements[e]->Node_Pointer.back()->X;

            if (L_bound <= x && x <= R_bound)
            {
                local_coord = (2 * x - L_bound - R_bound) / (R_bound - L_bound);

                return List_Of_Elements[e]->der_u_wrt_x(local_coord);
            }
        }
        cout << "ERROR: OUT OF MESH!" << endl;

        complex<double> empty(0, 0);
        return empty;
    }

    // #########################################################
    // WRITING FUNCTIONS:

    // write the nodal values to a file
    void write_nodal_values(ofstream &write_file)
    {
        // setting accuracy
        write_file << setw(20) << setprecision(18) << fixed;

        // writing nodal values
        for (int i = 0; i < Num_Nodes; i++)
        {
            write_file << List_Of_Nodes[i]->X << ", " << List_Of_Nodes[i]->U.real()
                       << ", " << List_Of_Nodes[i]->U.imag() << endl;
        }
    }

    // function to take a file and write interpolated solution at
    // equally spaced positions.
    void write_interpolated_values(ofstream &write_file, double L_Bound,
                                   double R_Bound, int Num_Points)
    {
        // setting accuracy
        write_file << setw(20) << setprecision(18) << fixed;

        // finding values to interpolate at
        double h = (R_Bound + L_Bound) / Num_Points;

        // writing nodal values
        complex<double> to_write;
        for (int i = 0; i < Num_Points; i++)
        {
            to_write = interpolated_u_glb_coor(L_Bound + i * h);
            write_file << L_Bound + i * h << " , " << to_write.real()
                       << " , " << to_write.imag() << endl;
        }
    }

    // #########################################################
    // TESTING QUADRATURE FUNCTIONS:

    // function used to overreide the gaussian quadrature points used.
    // it is assumed that every element will be override with the same number of
    void override_gauss_points(int n_gauss_)
    {
        override_gauss_quadrature = true;
        n_gauss = n_gauss_;

        if (n_gauss_ < 1 || n_gauss_ > gauss_values::GAUSS_POINTS.size())
        {
            throw invalid_argument("n out of range!");
        }
    }

    // #########################################################
    // PML FUNCTIONS:

    // set the elements that are in the PML
    void set_PML_elements(vector<int> PML_elements)
    {
        for (vector<int>::iterator ptr = PML_elements.begin(); ptr < PML_elements.end(); ptr++)
        {
            List_Of_Elements[*ptr]->PML_Element = true; // tell element it is a PML element
        }

        // find the PML domain
        double start_of_PML = List_Of_Elements[PML_elements[0]]->Node_Pointer[0]->X;
        X_PML = List_Of_Elements[PML_elements.back()]->Node_Pointer.back()->X;

        // print the PML domain
        cout << "PLM NODES SET \t D_PLM = [" << start_of_PML
             << " , " << X_PML << "]" << endl;
    }
};

#endif
