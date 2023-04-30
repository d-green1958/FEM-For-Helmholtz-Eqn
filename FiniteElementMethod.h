#ifndef FINITE_ELEMENT_METHOD
#define FINITE_ELEMENT_METHOD
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
#include "Mesh.h"
#include "PML.h"

using namespace std;

// some useful functions and prebuilt operations
namespace FEM
{
    // create a vector of element edge posisions used for mesh production
    // arguments:
    // - h : node spacing
    // - N : number of elements
    // - n : number of nodes per element (assumed constant)
    vector<double> form_equal_spacing(double h, int N, int n)
    {
        vector<double> Element_Edge_Pos;

        // now form vector
        for (int i = 0; i <= N; i++)
        {
            Element_Edge_Pos.push_back(i * (n - 1) * h);
        }

        return Element_Edge_Pos;
    }

    // runs convergence test over the domain approx [0,1]
    // arguments:
    // - <> : element type
    // - h : node spacing
    // - n : number of nodes per element (assumed constant)
    vector<double> CONVERGENCE_TEST(double h, double k, int n, bool sup = true)
    {

        // suppress cout
        if (sup == true)
        {
            cout.setstate(ios_base::failbit);
        }

        // need to use domain defined by number of elements
        int Num_Elements = ceil(1.0 / (h * (n - 1)));
        int Num_Nodes = Num_Elements * (n - 1) + 1;

        // Element edge positions
        vector<double> Element_Edge_Pos = form_equal_spacing(h, Num_Elements, n);
        vector<int> Nodes_Per_Element(Element_Edge_Pos.size() - 1, n);

        // Dirichelet Conditions
        vector<int> Dir_Boundary_Nodes{0};
        complex<double> start(0.0, 0.5);
        vector<complex<double>> Dir_Conditions{start};

        // Initial conditions
        vector<complex<double>> ICs(Num_Nodes, 0);

        // Neumann Conditions
        vector<int> Neu_Flux_Nodes{Num_Nodes - 1};
        complex<double> Const_Coef(0, 0);
        complex<double> U_Coef(0, k);

        // now perform the setup
        Mesh<Node> mesh(Element_Edge_Pos, Nodes_Per_Element);

        mesh.build(); // build the mesh
        // mesh.apply_ICs(ICs); // impose the ICs
        mesh.apply_homogenous_ICs();
        mesh.apply_Dirichelet_BCs(Dir_Conditions, Dir_Boundary_Nodes); // pin nodes
        mesh.assign_eq_nums();                                         // give equations numbers
        mesh.apply_Neumann_BCs(Const_Coef, U_Coef, Neu_Flux_Nodes);    // impose flux nodes
        mesh.initialise_res_and_jac();                                 // initialise the residuals and jacobians

        mesh.next_iteration(k);

        // now calcualte error:

        
        // taking approximation to integral over 101 values
        double test_spacing = 0.00001;
        int test_points = 100001;

        cout << "CALCULATING EXACT SOL" << endl;
        vector<complex<double>> Exact_Sol(test_points);
        for (int i = 0; i < test_points; i++)
        {
            Exact_Sol[i] = {-0.5 * sin(k * i * test_spacing),
                            0.5 * cos(k * i * test_spacing)};
        }

        vector<complex<double>> Exact_Der_Sol(test_points);
        for (int i = 0; i < test_points; i++)
        {
            Exact_Der_Sol[i] = {-0.5 * k * cos(k * i * test_spacing),
                                -0.5 * k * sin(k * i * test_spacing)};
        }

        // ||e||_0
        cout << "CALCULATING L0" << endl;
        double sum_2 = 0;
        for (int i = 0; i < test_points; i++)
        {
            sum_2 += abs(Exact_Sol[i] - mesh.interpolated_u_glb_coor(test_spacing * i)) *
                     abs(Exact_Sol[i] - mesh.interpolated_u_glb_coor(test_spacing * i)) *
                     test_spacing;
        }

        cout << "CALCULATING LE" << endl;
        // ||e||_E
        double sum_E = 0;
        for (int i = 0; i < test_points; i++)
        {
            sum_E += k * k * abs(Exact_Sol[i] - mesh.interpolated_u_glb_coor(test_spacing * i)) *
                         abs(Exact_Sol[i] - mesh.interpolated_u_glb_coor(test_spacing * i)) * test_spacing +
                     abs(Exact_Der_Sol[i] - mesh.interpolated_der_u_glb_coor(test_spacing * i)) *
                         abs(Exact_Der_Sol[i] - mesh.interpolated_der_u_glb_coor(test_spacing * i)) *
                         test_spacing;
        }

        cout.clear();
        cout << "################################## \n" << endl;
        cout << "k = " << k << endl;
        cout << "n = " << n << endl;
        cout << "L_2 Error = " << sqrt(sum_2) << endl;
        cout << "L_E Error = " << sqrt(sum_E) << endl;
        cout << "\n##################################" << endl;

        return {sqrt(sum_2), sqrt(sum_E)};
    }

    // returns L2 error from a single PML test run
    // arguments:
    // - N_PML : number of elements in PML
    // - h_PML : spacing between nodes in PML
    // - k : wavenumber (NOTE INT TYPE)
    // - n_pml : the number of nodes used per element in the PML
    double PML_CONVERGENCE_TEST(double N_PML, double h_PML, int k, int n_PML ,bool sup = true)
    {
        // suppress cout
        if (sup == true)
        {
            cout.setstate(std::ios_base::failbit);
        }

        // use 400 3 node elements in the bulk
        double h_comp = 0.001;
        double N_comp = 500;
        int n = 3;

        // form bulk elements vector
        vector<double> Element_Edge_Pos = FEM::form_equal_spacing(h_comp, N_comp, n);
         vector<int> Nodes_Per_Element(Element_Edge_Pos.size() - 1, n);

        // now add PML elements
        for (int i = 1; i <= N_PML; i++)
        {
            Element_Edge_Pos.push_back(Element_Edge_Pos[N_comp] +
                                       i * h_PML * (n_PML - 1));

            Nodes_Per_Element.push_back(n_PML);
        }


        Mesh<Node> mesh(Element_Edge_Pos, Nodes_Per_Element);
        mesh.build();

        // vector of the indicies of the PML elements
        int N = N_PML + N_comp;
        vector<int> PML_Elements;
        for (int i = N_comp; i < N; i++)
        {
            PML_Elements.push_back(i);
        }
        mesh.set_PML_elements(PML_Elements);

        // Dirichlet conditions
        int Num_Nodes = (n-1)*N_comp  + N_PML*(n_PML-1) + 1;
        vector<int> Dir_Boundary_Nodes{0, Num_Nodes-1};
        complex<double> start(0.0, 0.5);
        complex<double> end(0.0, 0.0);
        vector<complex<double>> Dir_Conditions{start, end};

        // calculate
        mesh.apply_homogenous_ICs();
        mesh.apply_Dirichelet_BCs(Dir_Conditions, Dir_Boundary_Nodes);
        mesh.assign_eq_nums();
        mesh.initialise_res_and_jac();
        mesh.next_iteration(k);

        mesh.display(true);


        // now calulate errors
        // using L^2 error

        double test_spacing = 0.0001;
        int test_points = 10001;

        // taking approximation to integral over 10001 values
        cout << "CALCULATING EXACT SOL" << endl;
        vector<complex<double>> Exact_Sol(test_points);
        for (int i = 0; i < test_points; i++)
        {
            Exact_Sol[i] = {-0.5 * sin(k * i * test_spacing),
                            0.5 * cos(k * i * test_spacing)};
        }

        cout << "CALCULATING L2 ERROR" << endl;
        double sum_2 = 0;
        for (int i = 0; i < test_points; i++)
        {
            sum_2 += abs(Exact_Sol[i] - mesh.interpolated_u_glb_coor(test_spacing * i)) *
                     abs(Exact_Sol[i] - mesh.interpolated_u_glb_coor(test_spacing * i)) *
                     test_spacing;
        }

        cout.clear();
        cout << "################################## \n"
             << endl;
        cout << "k = " << k << endl;
        cout << "n_PML = " << n_PML << endl;
        cout << "N_PML = " << N_PML << endl;
        cout << "delta_PML = " << N_PML * (n_PML - 1) * h_PML << endl;
        cout << "k * delta_PML = " << k * N_PML * (n_PML- 1) * h_PML << endl;
        cout << "L2 Error = " << sqrt(sum_2) << endl;
        cout << "\n##################################" << endl;

        ofstream outfile("values.csv");
        mesh.write_nodal_values(outfile);

        return sqrt(sum_2);
    }


    // return L2 error after completing test over domain [0,1].
    // importantly, allows for change in number of quadrature points
    // which is assumed the same over the entire mesh.
    // arguments:
    // - h : spacing
    // - k : wavenumber
    // - n : nodes per element
    // - n_gauss : number of gauss points used
    double GAUSSIAN_QUADRATURE_TEST(double h, double k, int n, int N, int n_gauss, bool sup = true)
    {
        // suppress cout
        if (sup == true)
        {
            cout.setstate(ios_base::failbit);
        }

        // need to use domain defined by number of elements
        int Num_Elements = N;
        int Num_Nodes = Num_Elements * (n - 1) + 1;

        cout << Num_Elements << "\t" << Num_Nodes << endl;

        // Element edge positions
        vector<double> Element_Edge_Pos = form_equal_spacing(h, Num_Elements, n);
        vector<int> Nodes_Per_Element(Element_Edge_Pos.size() - 1, n);

        // Dirichelet Conditions
        vector<int> Dir_Boundary_Nodes{0};
        complex<double> start(0.0, 0.5);
        vector<complex<double>> Dir_Conditions{start};

        // Initial conditions
        vector<complex<double>> ICs(Num_Nodes, 0);

        // Neumann Conditions
        vector<int> Neu_Flux_Nodes{Num_Nodes - 1};
        complex<double> Const_Coef(0, 0);
        complex<double> U_Coef(0, k);

        // now perform the setup
        Mesh<Node> mesh(Element_Edge_Pos, Nodes_Per_Element);

        mesh.build(); // build the mesh
        mesh.override_gauss_points(n_gauss);
        mesh.apply_homogenous_ICs();
        mesh.apply_Dirichelet_BCs(Dir_Conditions, Dir_Boundary_Nodes); // pin nodes
        mesh.assign_eq_nums();                                         // give equations numbers
        mesh.apply_Neumann_BCs(Const_Coef, U_Coef, Neu_Flux_Nodes);    // impose flux nodes
        mesh.initialise_res_and_jac();                                 // initialise the residuals and jacobians

        mesh.next_iteration(k);

        ofstream out_file("values_gauss_quad.csv");
        mesh.write_nodal_values(out_file);

        // now calcualte error
        // ||e||_0
        // taking approximation to integral over 101 values
        double L = mesh.List_Of_Nodes.back()->X;
        int test_points = 100001;
        double test_spacing = L/100000;

        vector<complex<double>> Exact_Sol(test_points);
        for (int i = 0; i < test_points; i++)
        {
            Exact_Sol[i] = {-0.5 * sin(k * i * test_spacing),
                            0.5 * cos(k * i * test_spacing)};
        }

        double sum_2 = 0;
        for (int i = 0; i < test_points; i++)
        {
            sum_2 += abs(Exact_Sol[i] - mesh.interpolated_u_glb_coor(test_spacing * i)) *
                     abs(Exact_Sol[i] - mesh.interpolated_u_glb_coor(test_spacing * i)) *
                     test_spacing;
        }

        cout.clear();
        cout << "################################## \n"
             << endl;
        cout << "k = " << k << endl;
        cout << "n = " << n << endl;
        cout << "N = " << N << endl;
        cout << "n_gauss " << n_gauss << endl;
        cout << "L2 Error = " << sqrt(sum_2) << endl;
        cout << "\n##################################" << endl;
        return sqrt(sum_2);
    }


    // function to quickly remind me of the order of operations
    void help()
    {
        cout << "\n \n";
        cout << "The correct order is:" << endl;
        cout << "1. Initialise a mesh giving it Element_Edge_Pos and Nodes_Per_Element" << endl;
        cout << "2. Build the mesh" << endl;
        cout << "3. Apply the intial conditions" << endl;
        cout << "4. Apply the Dirichlet conditions" << endl;
        cout << "5. Assign Eqn Numbers" << endl;
        cout << "6. Apply Neumann Conditions / set PML elements" << endl;
        cout << "7. Initalise the mesh Res and Jac" << endl;
        cout << "8. Calulate the iteration" << endl;
        cout << "(order is only a guide, can be changed a little)" << endl;

    }

}

#endif
