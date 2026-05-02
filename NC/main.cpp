#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <filesystem>
#include "include/Matrix.hpp"
#include "include/GaussElimination.hpp"
#include "include/GaussJacobi.hpp"
#include "include/LUCrout.hpp"
#include "include/LUDoolittle.hpp"
#include "include/LUCholesky.hpp"
#include "include/GaussSeidel.hpp"
#include "include/EigenValue.hpp"
#include "include/Lagrange.hpp"
#include "include/LeastSquares.hpp"

using namespace std;
ofstream fout;

// Output helper: writes to both cout and output file
void write(const string& s) {
    cout << s;
    fout << s;
}

// I/O helpers
bool readMatrix(const string& path, vector<vector<double>>& M, int& n) {
    ifstream fin(path);
    if (!fin) return false;
    int r, c;
    fin >> r >> c;
    n = r;
    M.assign(r, vector<double>(c));
    for (auto& row : M)
        for (auto& val : row)
            fin >> val;
    return true;
}

bool readVector(const string& path, vector<double>& b, int n) {
    ifstream fin(path);
    if (!fin) return false;
    b.resize(n);
    for (auto& val : b) fin >> val;
    return true;
}

// Read interpolation input file
bool readInterpolationInput(const string& path, vector<double>& xs,
                            vector<double>& ys, double& xQuery) {
    ifstream fin(path);
    if (!fin) return false;
    string line;
    int n = 0;
    while (getline(fin, line)) {
        if (!line.empty() && line[0] != '#') {
            n = stoi(line);
            break;
        }
    }
    if (n <= 0) return false;
    xs.resize(n);
    ys.resize(n);
    for (int i = 0; i < n; ++i)
        fin >> xs[i] >> ys[i];
    fin >> xQuery;
    return true;
}

// Read least squares input file into an (n x 2) Matrix
// Format: n, then n lines of "x f"
//   Matrix col-0 = xi,  col-1 = fi
// Returns an empty 0x0 Matrix on failure (check mat.isNull()).
Matrix readLeastSquaresInput(const string& path) {
    ifstream fin(path);
    if (!fin) return Matrix();           // empty sentinel
    string line;
    int n = 0;
    while (getline(fin, line)) {
        if (!line.empty() && line[0] != '#') {
            n = stoi(line);
            break;
        }
    }
    if (n <= 0) return Matrix();

    // Build n x 2 matrix; readFromFile reads row-major: x0 f0 x1 f1 ...
    Matrix mat(n, 2);
    mat.readFromFile(fin);               // uses Matrix's own file-reader
    return mat;
}

// Print helpers
void printMatrix(const vector<vector<double>>& M, int rows, const string& name) {
    write("\n" + name + "\n");
    for (int i = 0; i < rows; i++) {
        write("  ");
        for (int j = 0; j < (int)M[i].size(); j++) {
            char buf[32];
            snprintf(buf, sizeof(buf), "%12.6f", M[i][j]);
            write(string(buf));
        }
        write("\n");
    }
}

void printSolution(const vector<double>& x) {
    write("\nSolution:\n");
    for (int i = 0; i < (int)x.size(); i++) {
        char buf[64];
        snprintf(buf, sizeof(buf), "  x%d = %.8f\n", i + 1, x[i]);
        write(string(buf));
    }
}

void printVector(const vector<double>& v, const string& name) {
    write(name + " = [ ");
    for (double x : v) {
        char buf[32];
        snprintf(buf, sizeof(buf), "%.6f  ", x);
        write(string(buf));
    }
    write("]\n");
}

void printSeparator(const string& title = "") {
    if (!title.empty()) write("  " + title + "\n");
}

// Solve-system helpers
template<typename Solver>
void fillSolver(Solver& s, const vector<vector<double>>& A,
                const vector<double>& b, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) s(i, j) = A[i][j];
        s(i, n) = b[i];
    }
}

vector<double> solveSystem(const vector<vector<double>>& A, const vector<double>& b, int method, int n) {
    vector<double> x;
    switch (method) {
        case 1: { GaussElimination s(n,n+1); fillSolver(s,A,b,n); x=s.solve(); break; }
        case 2: { GaussJacobi      s(n,n+1); fillSolver(s,A,b,n); x=s.solve(); break; }
        case 3: { LUCrout          s(n,n+1); fillSolver(s,A,b,n); x=s.solve(); break; }
        case 4: { LUDoolittle      s(n,n+1); fillSolver(s,A,b,n); x=s.solve(); break; }
        case 5: { LUCholesky       s(n,n+1); fillSolver(s,A,b,n); x=s.solve(); break; }
        case 6: { GaussSeidel      s(n,n+1); fillSolver(s,A,b,n); x=s.solve(); break; }
    }
    return x;
}

void runSolve(const string& title, const vector<vector<double>>& A, const vector<double>& b, int n) {
    printSeparator(title);
    printMatrix(A, n, "Matrix A");
    write("\nChoose Solver:\n");
    write("  1. Gauss Elimination\n");
    write("  2. Gauss Jacobi\n");
    write("  3. LU Crout\n");
    write("  4. LU Doolittle\n");
    write("  5. LU Cholesky\n");
    write("  6. Gauss Seidel\n");
    write("Enter choice: ");
    int method; cin >> method;
    vector<double> x = solveSystem(A, b, method, n);
    printSolution(x);
}

void runEigenvalues(const vector<vector<double>>& Adata, int n) {
    Matrix mat(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            mat(i, j) = Adata[i][j];
    EigenValue ev(mat);

    printSeparator("EIGENVALUE ANALYSIS");
    printMatrix(Adata, n, "Matrix A  (n = " + to_string(n) + ")");

    write("\nSelect Eigenvalue Method:\n");
    write("  1. Power Method          (largest eigenvalue + eigenvector)\n");
    write("  2. Inverse Power Method  (smallest eigenvalue + eigenvector)\n");
    write("  3. QR Algorithm          (ALL eigenvalues, sorted by magnitude)\n");
    write("  4. Gershgorin Bounds     (eigenvalue region estimation)\n");
    write("  5. Condition Number      (cond = |λmax| / |λmin|)\n");
    write("  6. Run ALL methods\n");
    write("Enter choice: ");
    int choice; cin >> choice;
    write("\n");

    auto doPower = [&]() {
        printSeparator("Method 1: Power Iteration");
        write("Theory: v_{k+1} = A·v_k / ||A·v_k||  ->  converges to λmax\n\n");
        vector<double> vec;
        double lmax = ev.powerMethod(vec);
        char buf[64];
        snprintf(buf, sizeof(buf), "  λmax = %.10f\n", lmax);
        write(string(buf));
        write("  ");
        printVector(vec, "Eigenvector v");
    };
    auto doInverse = [&]() {
        printSeparator("Method 2: Inverse Power Iteration");
        write("Theory: solve A·v_{k+1} = v_k  ->  λmin = 1/μ\n\n");
        vector<double> vec;
        double lmin = ev.inversePowerMethod(vec);
        char buf[64];
        snprintf(buf, sizeof(buf), "  λmin = %.10f\n", lmin);
        write(string(buf));
        write("  ");
        printVector(vec, "Eigenvector v");
    };
    auto doQR = [&]() {
        printSeparator("Method 3: QR Algorithm");
        write("Theory: A_{k+1} = R_k·Q_k  (Wilkinson shift)\n\n");
        vector<double> eigs = ev.qrAlgorithm();
        write("  All eigenvalues (sorted descending by |λ|):\n");
        for (int i = 0; i < (int)eigs.size(); i++) {
            char buf[64];
            snprintf(buf, sizeof(buf), "    λ%-2d = %+.10f\n", i+1, eigs[i]);
            write(string(buf));
        }
        double prod = 1.0; for (double e : eigs) prod *= e;
        double traceSum = 0.0; for (double e : eigs) traceSum += e;
        char buf2[128];
        snprintf(buf2, sizeof(buf2), "\n  Product(λ) = %.6f\n  Sum(λ) = %.6f\n", prod, traceSum);
        write(string(buf2));
    };
    auto doGershgorin = [&]() {
        printSeparator("Method 4: Gershgorin Circle Theorem");
        auto disks = ev.gershgorinBounds();
        double globalLo = 1e30, globalHi = -1e30;
        for (auto& d : disks) {
            char buf[160];
            snprintf(buf, sizeof(buf),
                "  D_%d: center=%8.4f,  radius=%8.4f  =>  [%8.4f, %8.4f]\n",
                d.row+1, d.center, d.radius, d.lo(), d.hi());
            write(string(buf));
            if (d.lo() < globalLo) globalLo = d.lo();
            if (d.hi() > globalHi) globalHi = d.hi();
        }
        char buf2[128];
        snprintf(buf2, sizeof(buf2), "\n  Union: [%.4f, %.4f]\n", globalLo, globalHi);
        write(string(buf2));
    };
    auto doCond = [&]() {
        printSeparator("Condition Number");
        double cond = ev.conditionNumber();
        char buf[64];
        snprintf(buf, sizeof(buf), "  cond(A) = %.6f\n", cond);
        write(string(buf));
        if (cond > 1e6)       write("  [!] SEVERELY ill-conditioned\n");
        else if (cond > 1000) write("  [!] Ill-conditioned\n");
        else if (cond > 100)  write("  [~] Moderately conditioned\n");
        else                  write("  [OK] Well-conditioned\n");
    };

    switch (choice) {
        case 1: doPower();      break;
        case 2: doInverse();    break;
        case 3: doQR();         break;
        case 4: doGershgorin(); break;
        case 5: doCond();       break;
        case 6: doPower(); doInverse(); doQR(); doGershgorin(); doCond(); break;
        default: write("  [Error] Invalid choice.\n");
    }
}

// INTERPOLATION
void runInterpolation() {
    printSeparator("INTERPOLATION");
    vector<double> xs, ys;
    double xQuery;
    if (!readInterpolationInput("input/input_interpolation.txt", xs, ys, xQuery)) {
        write("[ERROR] Cannot open input/input_interpolation.txt\n");
        return;
    }
    int n = (int)xs.size();
    char buf[128];

    write("\nData Points:\n");
    write("   i        x            y\n");
    write("  ---  -----------  -----------\n");
    for (int i = 0; i < n; i++) {
        snprintf(buf, sizeof(buf), "  %3d  %11.6f  %11.6f\n", i, xs[i], ys[i]);
        write(string(buf));
    }
    snprintf(buf, sizeof(buf), "\nInterpolate at x = %.6f\n", xQuery);
    write(string(buf));

    write("\nChoose Interpolation Method:\n");
    write("  1. Lagrange Interpolation\n");
    write("Enter choice: ");
    int method; cin >> method;
    write("\n");

    switch (method) {
        case 1: {
            printSeparator("Lagrange Interpolation");
            write("Theory: P(x) = sum of y_i * L_i(x)\n");
            write("        L_i(x) = product_{j!=i} (x - x_j) / (x_i - x_j)\n\n");
            try {
                Lagrange lag(xs, ys);
                write("Basis polynomials at x_query:\n");
                for (int i = 0; i < n; i++) {
                    double L = lag.basisPoly(i, xQuery);
                    snprintf(buf, sizeof(buf),
                             "  L_%d(%.4f) = %12.8f   [y_%d = %.6f]\n",
                             i, xQuery, L, i, ys[i]);
                    write(string(buf));
                }
                double result = lag.interpolate(xQuery);
                write("\n");
                snprintf(buf, sizeof(buf), "  P(%.6f) = %.10f\n", xQuery, result);
                write(string(buf));
            } catch (const exception& e) {
                write(string("  [Error] ") + e.what() + "\n");
            }
            break;
        }
        default: write("  [Error] Invalid choice.\n");
    }
}

// ── LEAST SQUARES ─────────────────────────────────────────────────────────────
void runLeastSquares() {
    printSeparator("LEAST SQUARES METHOD");

    write("Theory:  Fit y = a*x + b  to N data points\n");
    write("         Minimize S = sum(fi - (a*xi + b))^2\n");
    write("         Normal equations:\n");
    write("           a*sum(xi^2) + b*sum(xi)  = sum(fi*xi)\n");
    write("           a*sum(xi)   + b*N        = sum(fi)\n\n");

    // ── Data flow: input_leastsquares.txt → Matrix(n x 2) → LeastSquares ─────
    // Matrix stores xi in col-0, fi in col-1 (inherited data[][]).
    // LeastSquares inherits rows/cols/data — no duplicate storage.
    Matrix dataMat = readLeastSquaresInput("input/input_leastsquares.txt");
    if (dataMat.isNull()) {
        write("[ERROR] Cannot open input/input_leastsquares.txt\n");
        write("        Expected format:\n");
        write("          n\n");
        write("          x_0  f_0\n");
        write("          x_1  f_1\n");
        write("          ...\n");
        return;
    }

    // Pass Matrix to LeastSquares constructor (inheritance)
    LeastSquares ls(dataMat);
    int N = ls.numPoints();           // uses inherited rows
    char buf[256];

    // Print data table using accessor methods xAt(i) / fAt(i)
    // which read from inherited data[i][0] and data[i][1]
    write("Data Points Table (read via Matrix inheritance):\n");
    write("   i     xi        fi        xi^2      fi*xi\n");
    write("  -----------------------------------------------------------\n");

    double sumX  = 0, sumX2 = 0, sumF = 0, sumFX = 0;
    for (int i = 0; i < N; i++) {
        double xi  = ls.xAt(i);       // inherited data[i][0]
        double fi  = ls.fAt(i);       // inherited data[i][1]
        double xi2 = xi * xi;
        double fxi = fi * xi;
        sumX  += xi;
        sumX2 += xi2;
        sumF  += fi;
        sumFX += fxi;
        snprintf(buf, sizeof(buf),
                 "  %3d  %8.4f  %8.4f  %9.4f  %9.4f\n",
                 i, xi, fi, xi2, fxi);
        write(string(buf));
    }
    snprintf(buf, sizeof(buf),
             "  Sum  %8.4f  %8.4f  %9.4f  %9.4f   N=%d\n",
             sumX, sumF, sumX2, sumFX, N);
    write(string(buf));

    // Normal equations display
    write("\nNormal Equations:\n");
    snprintf(buf, sizeof(buf),
             "  a * %.4f + b * %.4f = %.4f\n", sumX2, sumX, sumFX);
    write(string(buf));
    snprintf(buf, sizeof(buf),
             "  a * %.4f + b * %d   = %.4f\n", sumX, N, sumF);
    write(string(buf));

    // Matrix form display
    write("\nMatrix form  [A]{coeff} = {rhs}:\n");
    snprintf(buf, sizeof(buf),
             "  [ %9.4f  %9.4f ] [a]   [ %9.4f ]\n", sumX2, sumX, sumFX);
    write(string(buf));
    snprintf(buf, sizeof(buf),
             "  [ %9.4f  %9.4f ] [b] = [ %9.4f ]\n", sumX, (double)N, sumF);
    write(string(buf));

    // Solve and print results
    try {
        vector<double> coeff = ls.fit();   // returns {a, b}
        double a = coeff[0];
        double b = coeff[1];

        write("\n--- Solution ---\n");
        snprintf(buf, sizeof(buf), "  a = %.6f\n", a);
        write(string(buf));
        snprintf(buf, sizeof(buf), "  b = %.6f\n", b);
        write(string(buf));

        // Print fitted equation
        if (b >= 0)
            snprintf(buf, sizeof(buf), "\n  Fitted line:  y = %.4f * x + %.4f\n", a, b);
        else
            snprintf(buf, sizeof(buf), "\n  Fitted line:  y = %.4f * x - %.4f\n", a, -b);
        write(string(buf));

        // Print residuals table
        write("\nResiduals Table:\n");
        write("  -------------------------------------------------------\n");
        write("   i     xi        fi     fitted_yi   residual  residual^2\n");

        double S = 0.0;
        for (int i = 0; i < N; i++) {
            double xi     = ls.xAt(i);          // inherited data[i][0]
            double fi     = ls.fAt(i);          // inherited data[i][1]
            double yi_fit = ls.evaluate(xi, a, b);
            double resid  = fi - yi_fit;
            double resid2 = resid * resid;
            S += resid2;
            snprintf(buf, sizeof(buf),
                     "  %3d  %7.4f  %7.4f  %10.4f  %9.4f  %9.4f\n",
                     i, xi, fi, yi_fit, resid, resid2);
            write(string(buf));
        }
        snprintf(buf, sizeof(buf),
                 "  S = sum(residual^2) = %.6f\n", S);
        write(string(buf));

        // RMS error
        double rms = ls.rmsError(a, b);
        snprintf(buf, sizeof(buf),
                 "\n  RMS Error = sqrt(S / N) = sqrt(%.6f / %d)\n", S, N);
        write(string(buf));
        snprintf(buf, sizeof(buf),
                 "  RMS Error = %.6f\n", rms);
        write(string(buf));

    } catch (const exception& e) {
        write(string("  [Error] ") + e.what() + "\n");
    }
}

int main() {

    filesystem::create_directories("output");

    fout.open("output/output.txt");
    if (!fout.is_open()) {
        cerr << "[ERROR] Cannot open output/output.txt for writing.\n";
        return 1;
    }

    vector<vector<double>> A, Atilde;
    vector<double> b, btilde;
    int nA, nAt;

    bool hasMatrices =
        readMatrix("input/input_A.txt",      A,      nA)  &&
        readMatrix("input/input_Atilde.txt",  Atilde, nAt) &&
        readVector("input/input_b.txt",       b,      nA)  &&
        readVector("input/input_btilde.txt",  btilde, nAt);

    // TOP-LEVEL MENU
    cout << "Numerical Computation Toolbox\n";
    cout << "\nSelect Category:\n";
    cout << "  1. Matrix Operations  (Det / Inverse / SLE / Eigenvalues)\n";
    cout << "  2. Interpolation      (Lagrange)\n";
    cout << "  3. Least Squares      (Linear Fit  y = ax + b)\n";
    cout << "Enter choice: ";
    int category; cin >> category;

    if (category == 3) {
        runLeastSquares();
    }
    else if (category == 2) {
        runInterpolation();
    }
    else if (category == 1) {

        if (!hasMatrices) {
            cout << "[ERROR] Missing matrix input files!\n";
            fout.close();
            return 1;
        }

        cout << "\nSelect Matrix:\n";
        cout << "  1. Work with A\n";
        cout << "  2. Work with A~\n";
        cout << "  3. All 4 systems (solve/eigen for both)\n";
        cout << "Enter choice: ";
        int mode; cin >> mode;

        cout << "\nSelect Operation:\n";
        cout << "  1. Determinant\n";
        cout << "  2. Inverse\n";
        cout << "  3. Solve Linear System\n";
        cout << "  4. Eigenvalue Analysis\n";
        cout << "Enter choice: ";
        int op; cin >> op;

        if (mode == 1 || mode == 2) {
            vector<vector<double>> Ause = (mode == 1) ? A : Atilde;
            vector<double>         buse = (mode == 1) ? b : btilde;
            int n = (mode == 1) ? nA : nAt;

            Matrix mat(n, n);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    mat(i, j) = Ause[i][j];

            printMatrix(Ause, n, "Matrix");

            if (op == 1) {
                double det = mat.determinant();
                char buf[64];
                snprintf(buf, sizeof(buf), "\nDeterminant = %.8f\n", det);
                write(string(buf));
            }
            else if (op == 2) {
                try {
                    Matrix inv = mat.inverse();
                    write("\nInverse Matrix:\n");
                    for (int i = 0; i < n; i++) {
                        write("  ");
                        for (int j = 0; j < n; j++) {
                            char buf[32];
                            snprintf(buf, sizeof(buf), "%12.6f", inv(i, j));
                            write(string(buf));
                        }
                        write("\n");
                    }
                } catch (...) { write("\nMatrix is singular.\n"); }
            }
            else if (op == 3) { runSolve("Solving System", Ause, buse, n); }
            else if (op == 4) { runEigenvalues(Ause, n); }
            else              { write("\n[Error] Unknown operation.\n"); }
        }
        else if (mode == 3) {
            if (op == 3) {
                runSolve("1. A  x = b",  A,      b,      nA);
                runSolve("2. A~ x = b",  Atilde, b,      nAt);
                runSolve("3. A  x = b~", A,      btilde, nA);
                runSolve("4. A~ x = b~", Atilde, btilde, nAt);
            }
            else if (op == 4) {
                write("\n=== Eigenvalue Analysis: Matrix A ===\n");
                runEigenvalues(A, nA);
                write("\n=== Eigenvalue Analysis: Matrix A~ ===\n");
                runEigenvalues(Atilde, nAt);
            }
            else { write("\n[Error] Mode 3 only supports operations 3 and 4.\n"); }
        }
    }
    else {
        cout << "\n[Error] Invalid category choice.\n";
    }

    fout.close();
    cout << "\nOutput saved to output/output.txt\n";
    return 0;
}






