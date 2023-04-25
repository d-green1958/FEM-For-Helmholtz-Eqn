#ifndef NODES
#define NODES
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <typeinfo>

using namespace std;

// General Node Class
class Node
{
public:
    // position of node (1D)
    double X;

    // complex valued U (nodal value)
    complex<double> U;

    // local and global equation numbers
    int global_eqn_number;

    // defines if the node is a flux boundary node
    int flux_node = 0;
    // 0 if not a flux node (default)
    // 1 if flux node

    // conditions for flux element in the form:
    // u'(x_N) = (a+bi) + (c+di) u(x_N)
    complex<double> Const_Coef; // a+bi
    complex<double> U_Coef;     // c+di

    // overload << operator
    friend ostream &operator<<(ostream &os, Node node)
    {
        return os << "X: " << node.X << " U: " << node.U << " Global Eqn Num: " << node.global_eqn_number << " Flux Element: " << node.flux_node;
    }
};



#endif
