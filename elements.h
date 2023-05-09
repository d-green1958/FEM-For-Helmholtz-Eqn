#ifndef ELEMENTS
#define ELEMENTS
#include "nodes.h"
#include "gauss_values.h"
#include "PML.h"
#include "shape_functions.h"


using namespace std;

// Element class (pure virtual so requires derived class to use)
// contains pointers to nodes and some funcitons to calculate elemental
// Jacobian and residuals as well as some useful functions to interpolate u
// and u' at a point in the element. The shape functions, derivatives of
// shape functions, gauss points and gauss weights must be overriden.
class Element
{
public:
    // The list of pointers to nodes containted within the element
    vector<Node *> Node_Pointer;

    // elemental residuals and jacobian
    vector<complex<double>> Element_Residuals_Vec;
    vector<vector<complex<double>>> Element_Jacobian_Mat;

    // define if element is a PML element
    bool PML_Element = false;

    // #########################################################
    // CALCULATION FUNCTIONS:

    // loops over gauss points to calculate the residual at these points then update the
    // elemental residual vector
    // arguements:
    // - k : wavenumber
    // - X_PML : position of outerboundary PML layer
    // - override_gauss_points : bool to override the gaussian quadrature points (false by default)
    // - n_gauss : number of quadrature points we override with
    void update_residual_and_jacobian(double k, double X_PML, bool override_gauss_points = false, int n_gauss = 0)
    {

        vector<double> gauss_points;
        vector<double> gauss_weights;

        if (override_gauss_points)
        {
            gauss_points = gauss_values::GAUSS_POINTS[n_gauss - 1];
            gauss_weights = gauss_values::GAUSS_WEIGHTS[n_gauss - 1];
        }
        else
        {
            gauss_points = get_gauss_points();
            gauss_weights = get_gauss_weights();
        }

        // loop over gauss points
        for (int point = 0; point < gauss_points.size(); point++)
        {
            // counters for degrees of freedom
            int i_dof = 0;
            int j_dof = 0;

            for (int j_local = 0; j_local < Node_Pointer.size(); j_local++)
            {
                // check if node is pinned
                if (Node_Pointer[j_local]->global_eqn_number != -1)
                {

                    // update the vector of residuals
                    Element_Residuals_Vec[i_dof] +=
                        calc_residual_component(j_local, gauss_points[point], gauss_weights[point], k, X_PML);

                    j_dof = 0;
                    // update the jacobian matrix
                    for (int k_local = 0; k_local < Node_Pointer.size(); k_local++)
                    {
                        // checkif node is pinned
                        if (Node_Pointer[k_local]->global_eqn_number != -1)
                        {
                            Element_Jacobian_Mat[i_dof][j_dof] +=
                                calc_jacobian_component(j_local, k_local, gauss_points[point],
                                                        gauss_weights[point], k, X_PML);

                            // add degree of freedom
                            j_dof += 1;
                        }
                    }

                    // add degree of freedom
                    i_dof += 1;
                }
            }
        }
    }

    // calculate an residual vector component at a given gauss point
    virtual complex<double> calc_residual_component(int j_local, double gauss_point, double gauss_weight, double k, double X_PML)
    {

        if (PML_Element == false)
        {
            return (der_u_wrt_x(gauss_point) * der_shape_function_wrt_x(j_local, gauss_point) -
                    k * k * interpolated_u(gauss_point) * shape_function(j_local, gauss_point)) *
                   mapping_jacobian(gauss_point) * gauss_weight;
        }
        else if (PML_Element == true)
        {
            return (der_u_wrt_x(gauss_point) * der_shape_function_wrt_x(j_local, gauss_point) / PML::gamma(interpolated_x(gauss_point), X_PML, k) -
                    PML::gamma(interpolated_x(gauss_point), X_PML, k) * k * k * interpolated_u(gauss_point) * shape_function(j_local, gauss_point)) *
                   mapping_jacobian(gauss_point) * gauss_weight;
        }
        else
        {
            cout << "ERROR IN CALCULATING RESIDUAL COMPONENT" << endl;
            return 0;
        }
    }

    // calculate an elemental jacobian component at a given gauss point.
    virtual complex<double> calc_jacobian_component(int j_local, int k_local, double gauss_point, double gauss_weight,
                                                    double k, double X_PML)
    {

        if (PML_Element == false)
        {
            return (der_shape_function_wrt_x(k_local, gauss_point) * der_shape_function_wrt_x(j_local, gauss_point) -
                    k * k * shape_function(k_local, gauss_point) * shape_function(j_local, gauss_point)) *
                   mapping_jacobian(gauss_point) * gauss_weight;
        }
        else if (PML_Element == true)
        {
            return (der_shape_function_wrt_x(k_local, gauss_point) * der_shape_function_wrt_x(j_local, gauss_point) / PML::gamma(interpolated_x(gauss_point), X_PML, k) -
                    PML::gamma(interpolated_x(gauss_point), X_PML, k) * k * k * shape_function(k_local, gauss_point) * shape_function(j_local, gauss_point)) *
                   mapping_jacobian(gauss_point) * gauss_weight;
        }
        else
        {
            cout << "ERROR IN CALCULATING JACOBIAN COMPONENT" << endl;
            return 0;
        }
    }

    // function used to find the number of pinned nodes in the element. Used to define the size of the vectors
    int get_dof()
    {
        int unpinned_nodes = 0;
        for (int i = 0; i < Node_Pointer.size(); i++)
        {
            if (Node_Pointer[i]->global_eqn_number != -1)
            {
                unpinned_nodes += 1;
            }
        }
        return unpinned_nodes;
    }

    // empty the jacobian and residual matrix (for a second iteration)
    void clear_jac_and_res()
    {
        for (int i = 0; i < Element_Residuals_Vec.size(); i++)
        {
            Element_Residuals_Vec[i] = 0;
        }

        for (int i = 0; i < Element_Jacobian_Mat.size(); i++)
        {
            for (int j = 0; j < Element_Jacobian_Mat[i].size(); j++)
            {
                Element_Jacobian_Mat[i][j] = 0;
            }
        }
    }

    // #########################################################
    // SHAPE FUNCTIONS:

    // given the local node number j_local and the value s
    // evaluate the shape_function. (THIS IS A 2 NODE ELEMENT (CHANGE THIS))
    virtual double shape_function(int j_local, double s) = 0;

    // derivative of local shape function wrt to local coordinate s
    virtual double der_shape_function_wrt_s(int j_local, double s) = 0;

    // derivative of shape function wrt to global coordinate x
    double der_shape_function_wrt_x(int j_local, double s)
    {
        double denominator = 0;
        for (int i = 0; i < Node_Pointer.size(); i++)
        {
            denominator += Node_Pointer[i]->X * der_shape_function_wrt_s(i, s);
        }

        return der_shape_function_wrt_s(j_local, s) / denominator;
    }

    // #########################################################
    // GAUSSIAN QUADRATURE:

    // provide the required gauss points
    virtual vector<double> get_gauss_points() = 0;

    // provide the required gauss wights
    virtual vector<double> get_gauss_weights() = 0;

    // #########################################################
    // INTERPOLATION FUNCTIONS:

    // converts the local coordinate into the gloabl coordinate
    double interpolated_x(double s)
    {
        double temp = 0;
        for (int i = 0; i < Node_Pointer.size(); i++)
        {
            temp += Node_Pointer[i]->X * shape_function(i, s);
        }
        return temp;
    }

    // function to take local coordinate and find global coordinate
    double x2s(double x)
    {
        double L_Bound = Node_Pointer[0]->X;
        double R_Bound = Node_Pointer[Node_Pointer.size() - 1]->X;

        if (L_Bound <= x && x <= R_Bound)
        {
            return (2 * x - L_Bound - R_Bound) / (R_Bound - L_Bound);
        }
        cout << "OUTSIDE OF ELEMENT" << endl;
        return -100;
    }

    // interpolates u(x(s)) taking s as an argument
    complex<double> interpolated_u(double s)
    {
        complex<double> temp(0, 0);
        for (int i = 0; i < Node_Pointer.size(); i++)
        {
            temp += Node_Pointer[i]->U * shape_function(i, s);
        }
        return temp;
    }

    // finds derivative of u wrt to global coordinate, takes local coordinate as argument
    complex<double> der_u_wrt_x(double s)
    {
        complex<double> numerator(0, 0);
        double denomenator = 0;

        for (int i = 0; i < Node_Pointer.size(); i++)
        {
            numerator += der_shape_function_wrt_s(i, s) * Node_Pointer[i]->U;
            denomenator += der_shape_function_wrt_s(i, s) * Node_Pointer[i]->X;
        }

        return numerator / denomenator;
    }

    // Jacobain between local and global coordianates
    double mapping_jacobian(double s)
    {
        double temp = 0;
        for (int i = 0; i < Node_Pointer.size(); i++)
        {
            temp += Node_Pointer[i]->X * der_shape_function_wrt_s(i, s);
        }
        return temp;
    }

    // #########################################################
    // SETUP FUNCTIONS:

    // initialise the elemental jacobian matric and residuals vector of the correct size
    void initialise_element_res_and_jac()
    {
        int DOF = get_dof();
        Element_Residuals_Vec.resize(DOF);
        Element_Jacobian_Mat.resize(DOF);
        for (int i = 0; i < DOF; i++)
        {
            Element_Jacobian_Mat[i].resize(DOF);
        }
    }
};

// Two node element
class TwoNodeElement : public Element
{
    vector<double> get_gauss_points()
    {
        return {-1.0 / sqrt(3), +1.0 / sqrt(3)};
    }

    vector<double> get_gauss_weights()
    {
        return {1, 1};
    }

    // given the local node number j_local and the value s
    // evaluate the shape_function.
    double shape_function(int j_local, double s)
    {
        if (j_local == 0)
        {
            return 0.5 * (1 - s);
        }
        if (j_local == 1)
        {
            return 0.5 * (1 + s);
        }
        cout << "ERROR IN SHAPE FUNCTION" << endl;
        return -1;
    }

    // derivative of local shape function wrt to local coordinate s
    double der_shape_function_wrt_s(int j_local, double s)
    {
        if (j_local == 0)
        {
            return -0.5;
        }
        if (j_local == 1)
        {
            return 0.5;
        }
        cout << "ERROR IN SHAPE FUNCTION" << endl;
        return -1;
    }
};

// Three node element
class ThreeNodeElement : public Element
{
    vector<double> get_gauss_points()
    {
        return {-1.0 * sqrt(3.0 / 5.0), 0.0, 1.0 * sqrt(3.0 / 5.0)};
    }

    vector<double> get_gauss_weights()
    {
        return {5.0 / 9, 8.0 / 9.0, 5.0 / 9.0};
    }

    // given the local node number j_local and the value s
    // evaluate the shape_function.
    double shape_function(int j_local, double s)
    {
        if (j_local == 0)
        {
            return 0.5 * s * (s - 1);
        }
        if (j_local == 1)
        {
            return (1 + s) * (1 - s);
        }
        if (j_local == 2)
        {
            return 0.5 * s * (s + 1);
        }
        cout << "ERROR IN SHAPE FUNCTION" << endl;
        cout << "j_local = " << j_local << endl;
        cout << "s = " << s << endl;
        return -1;
    }

    // derivative of local shape function wrt to local coordinate s
    double der_shape_function_wrt_s(int j_local, double s)
    {
        if (j_local == 0)
        {
            return s - 0.5;
        }
        if (j_local == 1)
        {
            return -2 * s;
        }
        if (j_local == 2)
        {
            return s + 0.5;
        }
        cout << "ERROR IN DER SHAPE FUNCTION" << endl;
        return -1;
    }
};

// Four node element
class FourNodeElement : public Element
{
    vector<double> get_gauss_points()
    {
        return {-0.861136, -0.339981, +0.339981, +0.861136};
    }

    vector<double> get_gauss_weights()
    {
        return {0.347855, 0.652145, 0.652145, 0.347855};
    }

    // given the local node number j_local and the value s
    // evaluate the shape_function.
    double shape_function(int j_local, double s)
    {
        if (j_local == 0)
        {
            return (s - 1) * (s - 1.0 / 3) * (s + 1.0 / 3) * -9.0 / 16;
        }
        if (j_local == 1)
        {
            return (s - 1) * (s + 1) * (s - 1.0 / 3) * 27.0 / 16;
        }
        if (j_local == 2)
        {
            return (s - 1) * (s + 1) * (s + 1.0 / 3) * -27.0 / 16;
        }
        if (j_local == 3)
        {
            return (s + 1) * (s - 1.0 / 3) * (s + 1.0 / 3) * 9.0 / 16;
        }
        cout << "ERROR IN SHAPE FUNCTION" << endl;
        cout << "j_local = " << j_local << endl;
        cout << "s = " << s << endl;
        return -1;
    }

    // derivative of local shape function wrt to local coordinate s
    double der_shape_function_wrt_s(int j_local, double s)
    {
        if (j_local == 0)
        {
            return (-27.0 * s * s + 18 * s + 1) / 16;
        }
        if (j_local == 1)
        {
            return (81.0 * s * s - 18 * s - 27) / 16;
        }
        if (j_local == 2)
        {
            return (-81.0 * s * s - 18 * s + 27) / 16;
        }
        if (j_local == 3)
        {
            return (27.0 * s * s + 18 * s - 1) / 16;
        }
        cout << "ERROR IN DER SHAPE FUNCTION" << endl;
        return -1;
    }
};

// a generalised element with constructor to give it the number of basis functions
// it requires. (possibly a better method than using predefined n although may be slower
// either way, not worth the work required to implement this...)
class GeneralElement : public Element
{
public:
    int n;
    GeneralElement(int n_)
    {
        n = n_;
    }

    vector<double> get_gauss_points()
    {
        return gauss_values::GAUSS_POINTS[n - 1];
    }

    vector<double> get_gauss_weights()
    {
        return gauss_values::GAUSS_WEIGHTS[n - 1];
    }

    double shape_function(int j_local, double s)
    {
        return shape_functions::shape_function(j_local, s, n);
    }

    double der_shape_function_wrt_s(int j_local, double s)
    {
        return shape_functions::der_shape_function_wrt_s(j_local, s, n);
    }
};

#endif
