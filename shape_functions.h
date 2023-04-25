#ifndef SHAPE_FUNCTIONS
#define SHAPE_FUCNTIONS

using namespace std;

// namespace for shape functions and their 1st derivatives all in local coordinates
namespace shape_functions
{
    // shape functions for n=2,3,4
    // arguments:
    // - j_local : the local equation number
    // - s : the local coordinate
    // - n : the order of the shape function
    double shape_function(int j_local, double s, int n)
    {
        if (n == 2)
        {
            if (j_local == 0)
            {
                return 0.5 * (1 - s);
            }
            if (j_local == 1)
            {
                return 0.5 * (1 + s);
            }
        }

        if (n == 3)
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
        }

        if (n == 4)
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
        }

        throw invalid_argument("argument in shape functions is wrong!");
        return -10;
    }

    // derivative of shape functions wrt s for n=2,3,4
    // arguments:
    // - j_local : the local equation number
    // - s : the local coordinate
    // - n : the order of the shape function
    double der_shape_function_wrt_s(int j_local, double s, int n)
    {
        if (n == 2)
        {
            if (j_local == 0)
            {
                return -0.5;
            }
            if (j_local == 1)
            {
                return 0.5;
            }
        }

        if (n == 3)
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
        }

        if (n == 4)
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
        }

        throw invalid_argument("argument in der shape functions is wrong!");
        return -10;
    }
}

#endif