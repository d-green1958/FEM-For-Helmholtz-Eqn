#ifndef PML_funcs
#define PML_funcs
#include<complex>

using namespace std;

namespace PML
{
    // used to calculate gamma function in PML layer,
    // assumed the PML layer is attached to the right hand
    // of the 1D domain.
    complex<double> gamma(double x, double X_PML, double k)
    {
        return complex<double>(1, 1.0 / (k * (X_PML - x)));
    }
}


#endif

