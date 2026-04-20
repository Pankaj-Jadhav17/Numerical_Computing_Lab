#include "../include/Leastsquares.hpp"
#include <stdexcept>
#include <cmath>
using namespace std;

// Constructor
// Chain: LeastSquares() --> Matrix(n, 2)
//   Matrix(n, 2) sets:  rows = n,  cols = 2,  data[n][2] all zeros
// Then fills:
//   data[i][0] = xs[i]   (x value)
//   data[i][1] = fs[i]   (f/y value)
// NO own members — everything lives in Matrix's rows, cols, data

LeastSquares::LeastSquares(const vector<double>& xs, const vector<double>& fs)
    : Matrix((int)xs.size(), 2)     // rows = n, cols = 2
{
    if (xs.size() != fs.size() || xs.empty())
        throw invalid_argument("LeastSquares: x and f must be non-empty and same length.");

    // Fill Matrix data[][] with data points
    for (int i = 0; i < rows; ++i) {     // rows = N (from Matrix)
        data[i][0] = xs[i];               // col 0 = xi
        data[i][1] = fs[i];               // col 1 = fi
    }
}

// fit(): Builds and solves the 2x2 normal equation system
//
// From notes (image 1):
//   ∂S/∂a = 0  =>  a·Σxi²  + b·Σxi  = Σfi·xi
//   ∂S/∂b = 0  =>  a·Σxi   + b·N    = Σfi
//
// Written as  [A]{coeff} = {rhs}:
//   A    = [ Σxi²  Σxi ]
//          [ Σxi   N   ]
//
//   rhs  = [ Σfi·xi ]
//          [ Σfi    ]
//
// Solve using Matrix::inverse():
//   {a, b} = A⁻¹ · rhs
//
// Uses Matrix members ONLY:
//   rows        → N  (number of points)
//   data[i][0]  → xi
//   data[i][1]  → fi

vector<double> LeastSquares::fit() const
{
    // Accumulate sums directly from Matrix data[][]
    double sumX  = 0.0;   // Σxi
    double sumX2 = 0.0;   // Σxi²
    double sumF  = 0.0;   // Σfi
    double sumFX = 0.0;   // Σfi·xi
    int    N     = rows;  // rows = number of data points (from Matrix)

    for (int i = 0; i < N; ++i) {
        double xi = data[i][0];   // xi from Matrix col 0
        double fi = data[i][1];   // fi from Matrix col 1
        sumX  += xi;
        sumX2 += xi * xi;
        sumF  += fi;
        sumFX += fi * xi;
    }

    // Build 2x2 normal equation matrix A using Matrix constructor
    //   A = [ Σxi²   Σxi ]
    //       [ Σxi    N   ]
    Matrix A(2, 2);
    A(0, 0) = sumX2;          // Σxi²
    A(0, 1) = sumX;           // Σxi
    A(1, 0) = sumX;           // Σxi
    A(1, 1) = (double)N;      // N  (total points)

    // Check determinant before inverting (using Matrix::determinant())
    if (abs(A.determinant()) < 1e-12)
        throw runtime_error("LeastSquares: Normal equation matrix is singular. "
                            "Check for duplicate or collinear x-values.");

    // Invert A using Matrix::inverse()
    Matrix Ainv = A.inverse();

    // Multiply Ainv by rhs vector to get [a, b]
    //   a = Ainv(0,0)*sumFX + Ainv(0,1)*sumF
    //   b = Ainv(1,0)*sumFX + Ainv(1,1)*sumF
    double a = Ainv(0, 0) * sumFX + Ainv(0, 1) * sumF;
    double b = Ainv(1, 0) * sumFX + Ainv(1, 1) * sumF;

    return {a, b};
}

// evaluate(): Returns y = a*x + b  at given x
double LeastSquares::evaluate(double x, double a, double b) const
{
    return a * x + b;
}

// rmsError(): Computes Root Mean Square error
//
// From notes (image 1):
//   rms = sqrt( (1/(N+1)) * Σ(yi - fi)² )
//   where fi = a*xi + b  (fitted value)
//         yi = actual data value
//
// Uses Matrix members:
//   rows        → N
//   data[i][0]  → xi
//   data[i][1]  → actual fi (yi)

double LeastSquares::rmsError(double a, double b) const
{
    double S = 0.0;
    int    N = rows;   // rows from Matrix

    for (int i = 0; i < N; ++i) {
        double xi      = data[i][0];          // xi from Matrix
        double yi      = data[i][1];          // actual y from Matrix
        double fi      = a * xi + b;          // fitted value
        double residual = yi - fi;
        S += residual * residual;             // S = Σ(yi - fi)²
    }

    // rms = sqrt( S / N )  — matches notes formula
    return sqrt(S / (double)N);
}