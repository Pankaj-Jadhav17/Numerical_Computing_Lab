#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "include/Matrix.hpp"
#include "include/GaussElimination.hpp"
#include "include/GaussJacobi.hpp"
#include "include/LUCrout.hpp"
#include "include/LUDoolittle.hpp"
#include "include/LUCholesky.hpp"
#include "include/GaussSeidel.hpp"
#include "include/EigenValue.hpp"
#include "include/Lagrange.hpp"

using namespace std;

ofstream fout;

//  Output helper: writes to both cout and output file
void write(const string& s) {
    cout << s;
    fout << s;
}

//  I/O helpers
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
// Format: n, then n lines of "x y", then one query value
bool readInterpolationInput(const string& path, vector<double>& xs, vector<double>& ys,double& xQuery) {
    ifstream fin(path);
    if (!fin) return false;

    // skip comment lines
    string line;
    int n = 0;
    while (getline(fin, line)) {
        if (!line.empty() && line[0] != '#') {
            // first non-comment line is n
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

//  Print helpers
void printMatrix(const vector<vector<double>>& M, int rows,
                 const string& name) {
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
    write("\n══════════════════════════════════════════════\n");
    if (!title.empty()) write("  " + title + "\n");
    write("══════════════════════════════════════════════\n");
}

//  Solve-system helpers (methods 1-6)
template<typename Solver>
void fillSolver(Solver& s, const vector<vector<double>>& A, const vector<double>& b, int n) {
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
        write("Theory: v_{k+1} = A·v_k / ||A·v_k||  →  converges to λmax\n\n");
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
        write("Theory: solve A·v_{k+1} = v_k  →  λmin = 1/μ\n");
        write("        Eigenvector of A⁻¹ for largest eigenvalue\n\n");
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
        write("Theory: A_{k+1} = R_k·Q_k  (Wilkinson shift applied)\n");
        write("        Diagonal of A_∞ gives all eigenvalues\n");
        write("        det(A) = product of all eigenvalues\n\n");
        vector<double> eigs = ev.qrAlgorithm();
        write("  All eigenvalues (sorted descending by |λ|):\n");
        for (int i = 0; i < (int)eigs.size(); i++) {
            char buf[64];
            snprintf(buf, sizeof(buf), "    λ%-2d = %+.10f\n", i+1, eigs[i]);
            write(string(buf));
        }
        double prod = 1.0;
        for (double e : eigs) prod *= e;
        char buf2[128];
        snprintf(buf2, sizeof(buf2),
            "\n  Product(λ) = %.6f  (≈ det(A) — Cayley-Hamilton)\n", prod);
        write(string(buf2));
        double traceSum = 0.0;
        for (double e : eigs) traceSum += e;
        snprintf(buf2, sizeof(buf2),
            "  Sum(λ)     = %.6f  (= trace of A)\n", traceSum);
        write(string(buf2));
    };

    auto doGershgorin = [&]() {
        printSeparator("Method 4: Gershgorin Circle Theorem");
        write("Theory: Each eigenvalue lies in at least one disk D_i\n");
        write("        D_i = { z : |z - a_ii| ≤ R_i }\n");
        write("        R_i = Σ_{j≠i} |a_ij|   (sum of off-diagonal |entries|)\n\n");
        auto disks = ev.gershgorinBounds();
        double globalLo =  1e30, globalHi = -1e30;
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
        snprintf(buf2, sizeof(buf2),
            "\n  Union (all eigenvalues ∈): [%.4f, %.4f]\n",
            globalLo, globalHi);
        write(string(buf2));
    };

    auto doCond = [&]() {
        printSeparator("Condition Number");
        write("Theory: cond(A) = ||A||·||A⁻¹||  ≥  |λmax| / |λmin|\n");
        write("        High cond(A) => ill-conditioned => errors amplified\n\n");
        double cond = ev.conditionNumber();
        char buf[64];
        snprintf(buf, sizeof(buf), "  cond(A) = %.6f\n", cond);
        write(string(buf));
        if (cond > 1e6)
            write("  [!] SEVERELY ill-conditioned (cond > 1e6)\n");
        else if (cond > 1000)
            write("  [!] Ill-conditioned: errors in b are amplified ~"
                  + to_string((int)cond) + "x in solution x\n");
        else if (cond > 100)
            write("  [~] Moderately conditioned\n");
        else
            write("  [OK] Well-conditioned matrix\n");
    };

    switch (choice) {
        case 1: doPower();      break;
        case 2: doInverse();    break;
        case 3: doQR();         break;
        case 4: doGershgorin(); break;
        case 5: doCond();       break;
        case 6:
            doPower();
            doInverse();
            doQR();
            doGershgorin();
            doCond();
            break;
        default:
            write("  [Error] Invalid eigenvalue method choice.\n");
    }
}

// ── INTERPOLATION ────────────────────────────────────────
void runInterpolation() {
    printSeparator("INTERPOLATION");

    // Read data from file
    vector<double> xs, ys;
    double xQuery;
    if (!readInterpolationInput("input/input_interpolation.txt", xs, ys, xQuery)) {
        write("[ERROR] Cannot open input/input_interpolation.txt\n");
        write(" Expected format:\n");
        write(" n\n");
        write(" x_0  y_0\n");
        write(" x_1  y_1\n");
        write(" ...\n");
        write(" x_query\n");
        return;
    }

    int n = (int)xs.size();
    char buf[128];

    // Print data table
    write("\nData Points:\n");
    write("   i        x            y\n");
    write("  ---  -----------  -----------\n");
    for (int i = 0; i < n; i++) {
        snprintf(buf, sizeof(buf), "  %3d  %11.6f  %11.6f\n", i, xs[i], ys[i]);
        write(string(buf));
    }
    snprintf(buf, sizeof(buf), "\nInterpolate at x = %.6f\n", xQuery);
    write(string(buf));

    // Method selection
    write("\nChoose Interpolation Method:\n");
    write("  1. Lagrange Interpolation\n");
    write("Enter choice: ");
    int method; cin >> method;
    write("\n");

    switch (method) {
        case 1: {
            printSeparator("Lagrange Interpolation");
            write("Theory: P(x) = Σ y_i · L_i(x)\n");
            write("        L_i(x) = Π_{j≠i} (x - x_j) / (x_i - x_j)\n\n");

            try {
                Lagrange lag(xs, ys);

                // Show each basis polynomial value
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
                snprintf(buf, sizeof(buf),
                         "  P(%.6f) = %.10f\n", xQuery, result);
                write(string(buf));
            } catch (const exception& e) {
                write(string("  [Error] ") + e.what() + "\n");
            }
            break;
        }
        default:
            write("  [Error] Invalid interpolation method choice.\n");
    }
}

int main() {

    vector<vector<double>> A, Atilde;
    vector<double> b, btilde;
    int nA, nAt;

    bool hasMatrices =
        readMatrix("input/input_A.txt",      A,      nA)  &&
        readMatrix("input/input_Atilde.txt",  Atilde, nAt) &&
        readVector("input/input_b.txt",       b,      nA)  &&
        readVector("input/input_btilde.txt",  btilde, nAt);

    fout.open("output/output.txt");

    // ── TOP-LEVEL OPERATION ──────────────────────────────────
    cout <<"Numerical Computation Toolbox \n";
    cout << "\nSelect Category:\n";
    cout << "  1. Matrix Operations  (Det / Inverse / SLE / Eigenvalues)\n";
    cout << "  2. Interpolation      (Lagrange)\n";
    cout << "Enter choice: ";
    int category; cin >> category;

    if (category == 2) {
        runInterpolation();
    }
    else if (category == 1) {

        if (!hasMatrices) {
            cout << "[ERROR] Missing matrix input files!\n";
            fout.close();
            return 1;
        }

        // ── SELECT MATRIX ────────────────────────────────────
        cout << "\nSelect Matrix:\n";
        cout << "  1. Work with A\n";
        cout << "  2. Work with A~\n";
        cout << "  3. All 4 systems (solve/eigen for both)\n";
        cout << "Enter choice: ";
        int mode; cin >> mode;

        // ── SELECT OPERATION ─────────────────────────────────
        cout << "\nSelect Operation:\n";
        cout << "  1. Determinant\n";
        cout << "  2. Inverse\n";
        cout << "  3. Solve Linear System\n";
        cout << "  4. ★ Eigenvalue Analysis\n";
        cout << "Enter choice: ";
        int op; cin >> op;

        // ── SINGLE MATRIX MODE ───────────────────────────────
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
                } catch (...) {
                    write("\nMatrix is singular — no inverse exists.\n");
                }
            }
            else if (op == 3) {
                runSolve("Solving System", Ause, buse, n);
            }
            else if (op == 4) {
                runEigenvalues(Ause, n);
            }
            else {
                write("\n[Error] Unknown operation.\n");
            }
        }

        // ── ALL SYSTEMS MODE ─────────────────────────────────
        else if (mode == 3) {
            if (op == 3) {
                runSolve("1. A  x = b",   A,      b,      nA);
                runSolve("2. A~ x = b",   Atilde, b,      nAt);
                runSolve("3. A  x = b~",  A,      btilde, nA);
                runSolve("4. A~ x = b~",  Atilde, btilde, nAt);
            }
            else if (op == 4) {
                write("\n=== Eigenvalue Analysis: Matrix A ===\n");
                runEigenvalues(A, nA);
                write("\n=== Eigenvalue Analysis: Matrix A~ ===\n");
                runEigenvalues(Atilde, nAt);
            }
            else {
                write("\n[Error] Mode 3 only supports operations 3 and 4.\n");
            }
        }
    }
    else {
        cout << "\n[Error] Invalid category choice.\n";
    }

    fout.close();
    cout << "\nOutput saved to output/output.txt\n";
    return 0;
}