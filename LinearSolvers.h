#ifndef LINEAR_SOLVERS
#define LINEAR_SOLVERS
#include <iomanip>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <typeinfo>

using namespace std;

namespace LinearSolvers
{
    // linear equation solver (tridiagonal) using thomas method
    vector<complex<double>> Thomas_Solver(const vector<vector<complex<double>>> &A,
                                          const vector<complex<double>> &B)
    {
        int n = B.size();
        vector<complex<double>> sol(n);
        vector<complex<double>> Y(n - 1);
        vector<complex<double>> P(n);

        if (n != A[0].size())
        {
            cout << "ERROR IN SOLVER" << endl;
        }

        // forward phase
        //  1st row by hand
        Y[0] = A[0][1] / A[0][0];
        P[0] = B[0] / A[0][0];

        // then the rest
        for (int i = 1; i < n; i++)
        {
            P[i] = (B[i] - A[i][i - 1] * P[i - 1]) / (A[i][i] - A[i][i - 1] * Y[i - 1]);
            if (i < n - 1)
            {
                Y[i] = A[i][i + 1] / (A[i][i] - A[i][i - 1] * Y[i - 1]);
            }
        }

        // backwards phase
        sol[n - 1] = P[n - 1];
        for (int i = n - 2; i > -1; i--)
        {
            sol[i] = P[i] - Y[i] * sol[i + 1];
        }

        return sol;
    }

    // LU solver for general matrix (doesnt use sparsness!)
    vector<complex<double>> LU_Solver(const vector<vector<complex<double>>> &A,
                                      const vector<complex<double>> &B)
    {
        int n = B.size();

        // for solution storage
        vector<complex<double>> sol(n);

        // The lower and upper parts of the matrix
        vector<vector<complex<double>>> L(n);
        vector<vector<complex<double>>> U(n);

        // resize the matrix
        for (int i = 0; i < n; i++)
        {
            L[i].resize(n);
            U[i].resize(n);
        }

        // form L and U
        for (int i = 0; i < n; i++)
        {
            // upper - U
            for (int k = i; k < n; k++)
            {
                // sum of L(i, j) * U(j, k)
                complex<double> sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (L[i][j] * U[j][k]);
                U[i][k] = A[i][k] - sum;
            }

            // lower - L
            for (int k = i; k < n; k++)
            {
                if (i == k)
                    L[i][i] = 1; // set diagonal to 1
                else
                {
                    // sum of L(k, j) * U(j, i)
                    complex<double> sum = 0;
                    for (int j = 0; j < i; j++)
                        sum += (L[k][j] * U[j][i]);

                    L[k][i] = (A[k][i] - sum) / U[i][i];
                }
            }
        }

        vector<complex<double>> Y(n);

        // now need to solve the system Ly = B
        for (int i = 0; i < n; i++)
        {
            complex<double> sum(0, 0);
            for (int j = 0; j <= i; j++)
            {
                sum += L[i][j] * Y[j];
            }

            Y[i] = B[i] - sum;
        }

        // now solve for system Ux = y
        for (int i = n - 1; i >= 0; i--)
        {
            complex<double> sum(0, 0);
            for (int j = n - 1; j > i; j--)
            {
                sum += U[i][j] * sol[j];
            }

            sol[i] = (Y[i] - sum) / U[i][i];
        }

        return sol;
    }
}

#endif