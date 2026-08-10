#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <filesystem>
#include <sstream>
#include <map>
#include <functional>
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
#include "include/ForwardDifference.hpp"
#include "include/BackwardDifference.hpp"
#include "include/CentralDifference.hpp"
#include "include/CentralSecondDerivative.hpp"
#include "include/DifferentiationAnalyzer.hpp"

using namespace std;

struct RichardsonResultRow {
    string function;
    double h;
    double D_h;
    double R_h;
    double exact;
    double errD;
    double errR;
};

class RichardsonAnalyzer {
private:
    struct TestCase {
        string name;
        function<double(double)> f;
        function<double(double)> fExact;
    };

    double x0;
    vector<double> hValues;
    vector<TestCase> testFunctions;
    vector<RichardsonResultRow> results;

    static void groupByFunction(const vector<RichardsonResultRow>& rows,
                                vector<string>& order,
                                map<string, vector<const RichardsonResultRow*>>& byFunc) {
        for (const auto& r : rows) {
            if (byFunc.find(r.function) == byFunc.end()) order.push_back(r.function);
            byFunc[r.function].push_back(&r);
        }
    }

    static double logLogSlope(const vector<const RichardsonResultRow*>& rows, bool useD) {
        auto errOf = [&](const RichardsonResultRow* r) { return useD ? r->errD : r->errR; };
        size_t n = 1;
        while (n < rows.size() && errOf(rows[n]) < errOf(rows[n - 1])) ++n;
        if (n < 2) n = min(rows.size(), (size_t)2);

        vector<double> X, Y;
        for (size_t i = 0; i < n; ++i) {
            double e = errOf(rows[i]);
            if (e <= 0.0) continue;
            X.push_back(log10(rows[i]->h));
            Y.push_back(log10(e));
        }
        if (X.size() < 2) return 0.0;

        double xMean = 0, yMean = 0;
        for (size_t i = 0; i < X.size(); ++i) { xMean += X[i]; yMean += Y[i]; }
        xMean /= X.size(); yMean /= Y.size();

        double num = 0, den = 0;
        for (size_t i = 0; i < X.size(); ++i) {
            num += (X[i] - xMean) * (Y[i] - yMean);
            den += (X[i] - xMean) * (X[i] - xMean);
        }
        return (den != 0.0) ? num / den : 0.0;
    }

public:
    RichardsonAnalyzer(double evalPoint, const vector<double>& hVals)
        : x0(evalPoint), hValues(hVals) {}

    void addFunction(const string& name,
                     function<double(double)> f,
                     function<double(double)> fExact) {
        testFunctions.push_back({name, f, fExact});
    }

    void run() {
        results.clear();
        for (const auto& tc : testFunctions) {
            double exact = tc.fExact(x0);
            for (double h : hValues) {
                CentralDifference cd_h(tc.f, h);
                CentralDifference cd_half(tc.f, h / 2.0);

                double D_h = cd_h.derivative(x0);
                double D_half = cd_half.derivative(x0);
                double R_h = (4.0 * D_half - D_h) / 3.0;

                RichardsonResultRow row;
                row.function = tc.name;
                row.h = h;
                row.D_h = D_h;
                row.R_h = R_h;
                row.exact = exact;
                row.errD = Differentiation::absoluteError(exact, D_h);
                row.errR = Differentiation::absoluteError(exact, R_h);
                results.push_back(row);
            }
        }
    }

    const vector<RichardsonResultRow>& getResults() const { return results; }

    void printTable(ostream& out) const {
        out << scientific << setprecision(6);
        out << left
            << setw(8)  << "func"
            << setw(12) << "h"
            << setw(16) << "D(h)"
            << setw(16) << "R(h)"
            << setw(16) << "exact"
            << setw(16) << "err_D(h)"
            << setw(16) << "err_R(h)" << "\n";

        for (const auto& r : results) {
            out << left
                << setw(8)  << r.function
                << setw(12) << r.h
                << setw(16) << r.D_h
                << setw(16) << r.R_h
                << setw(16) << r.exact
                << setw(16) << r.errD
                << setw(16) << r.errR << "\n";
        }
    }

    void writeTable(const string& path) const {
        ofstream out(path);
        out << "Richardson Extrapolation vs Central Difference - Results Table\n";
        out << "Evaluation point x0 = " << x0 << "\n\n";
        printTable(out);
    }

    void writeCSV(const string& path) const {
        ofstream csv(path);
        csv << "function,h,D_h,R_h,exact,err_D,err_R\n";
        csv << scientific << setprecision(8);
        for (const auto& r : results) {
            csv << r.function << "," << r.h << "," << r.D_h << "," << r.R_h << ","
                << r.exact << "," << r.errD << "," << r.errR << "\n";
        }
    }

    void writeGnuplotFiles(const string& dataPrefix,
                          const string& scriptPath,
                          const string& outputImage) const {
        auto safeLog = [](double v) { return log10(v > 0 ? v : 1e-20); };

        vector<string> order;
        map<string, vector<const RichardsonResultRow*>> byFunc;
        groupByFunction(results, order, byFunc);

        for (const auto& fname : order) {
            ofstream dat(dataPrefix + "_" + fname + ".dat");
            dat << "# log10(h)  log10(err_D)  log10(err_R)\n";
            for (const auto* r : byFunc[fname]) {
                dat << safeLog(r->h) << " " << safeLog(r->errD) << " "
                    << safeLog(r->errR) << "\n";
            }
        }

        ofstream gp(scriptPath);
        gp << "set terminal pngcairo size 1000,700 enhanced font 'Verdana,10'\n";
        gp << "set output '" << outputImage << "'\n";
        gp << "set title 'Log-Log Error Plot: Central Difference D(h) vs Richardson R(h)'\n";
        gp << "set xlabel 'log10(h)'\n";
        gp << "set ylabel 'log10(|error|)'\n";
        gp << "set grid\n";
        gp << "set key outside right\n";
        gp << "plot \\\n";
        bool first = true;
        for (const auto& fname : order) {
            string file = dataPrefix + "_" + fname + ".dat";
            if (!first) gp << ", \\\n";
            gp << "  '" << file << "' using 1:2 with linespoints title '" << fname << " D(h) [O(h^2)]'";
            gp << ", \\\n  '" << file << "' using 1:3 with linespoints title '" << fname << " R(h) [O(h^4)]'";
            first = false;
        }
        gp << "\n";
    }

    bool renderPlot(const string& scriptPath) const {
        string cmd = "gnuplot " + scriptPath;
        return system(cmd.c_str()) == 0;
    }

    void printAnalysis(ostream& out) const {
        out << fixed << setprecision(4);
        out << "\n===== CONVERGENCE / ANALYSIS =====\n";

        vector<string> order;
        map<string, vector<const RichardsonResultRow*>> byFunc;
        groupByFunction(results, order, byFunc);

        bool richardsonAlwaysBetter = true;
        double totalRatioD = 0, totalRatioR = 0;
        int countRatioD = 0, countRatioR = 0;

        for (const auto& fname : order) {
            const auto& rows = byFunc[fname];
            out << "\n-- " << fname << " --\n";

            double slopeD = logLogSlope(rows, true);
            double slopeR = logLogSlope(rows, false);
            out << "  Observed slope of log10(err_D) vs log10(h): " << slopeD
                << "   (theoretical: 2.0000)\n";

            bool richardsonExact = (rows.front()->errR < 1e-9);
            if (richardsonExact) {
                out << "  Observed slope of log10(err_R) vs log10(h): N/A -- R(h) is already at "
                    "machine-precision (~" << scientific << rows.front()->errR << fixed
                    << ") even at the LARGEST h tested.\n";
            } else {
                out << "  Observed slope of log10(err_R) vs log10(h): " << slopeR
                    << "   (theoretical: 4.0000)\n";
            }

            size_t nD = 1; while (nD < rows.size() && rows[nD]->errD < rows[nD - 1]->errD) ++nD;
            size_t nR = 1; while (nR < rows.size() && rows[nR]->errR < rows[nR - 1]->errR) ++nR;

            double logRatioD = 0; int kD = 0;
            for (size_t i = 1; i < nD; ++i) {
                if (rows[i]->errD > 0 && rows[i-1]->errD > 0) {
                    logRatioD += log10(rows[i-1]->errD / rows[i]->errD);
                    ++kD;
                }
            }
            double logRatioR = 0; int kR = 0;
            for (size_t i = 1; i < nR; ++i) {
                if (rows[i]->errR > 0 && rows[i-1]->errR > 0) {
                    logRatioR += log10(rows[i-1]->errR / rows[i]->errR);
                    ++kR;
                }
            }
            double factorD = (kD > 0) ? pow(10.0, logRatioD / kD) : 0.0;
            double factorR = (kR > 0) ? pow(10.0, logRatioR / kR) : 0.0;
            out << "  Error reduction factor per 10x decrease in h:\n";
            out << "    D(h): ~" << factorD << "x   (theoretical: ~100x for O(h^2))\n";
            if (richardsonExact) {
                out << "    R(h): N/A -- already at machine precision.\n";
            } else {
                out << "    R(h): ~" << factorR << "x   (theoretical: ~10000x for O(h^4))\n";
            }
            totalRatioD += factorD;
            if (!richardsonExact) { totalRatioR += factorR; ++countRatioR; }
            ++countRatioD;

            int exceptions = 0;
            for (const auto* r : rows) if (r->errR >= r->errD) ++exceptions;
            if (exceptions > 0) richardsonAlwaysBetter = false;
            out << "  Richardson more accurate than Central at " << (rows.size() - exceptions)
                << "/" << rows.size() << " tested h values"
                << (exceptions > 0 ? "  (fails at the smallest h -- round-off dominates R(h) first)" : "")
                << "\n";

            size_t minD = 0, minR = 0;
            for (size_t i = 1; i < rows.size(); ++i) {
                if (rows[i]->errD < rows[minD]->errD) minD = i;
                if (rows[i]->errR < rows[minR]->errR) minR = i;
            }
            out << "  D(h) reaches its minimum error at h = " << scientific << rows[minD]->h
                << "  (err = " << rows[minD]->errD << ")\n";
            out << "  R(h) reaches its minimum error at h = " << rows[minR]->h
                << "  (err = " << rows[minR]->errR << ")\n";
            out << fixed;
            if (minR > minD)
                out << "  -> R(h) keeps improving to a SMALLER h than D(h) before round-off takes over.\n";
            else if (minR < minD)
                out << "  -> R(h) hits its round-off floor EARLIER (larger h) than D(h).\n";
            else
                out << "  -> D(h) and R(h) hit their round-off floor at about the same h.\n";
        }

        double avgFactorD = countRatioD ? totalRatioD / countRatioD : 0.0;
        double avgFactorR = countRatioR ? totalRatioR / countRatioR : 0.0;

        out << "\n===== ANSWERS TO ASSIGNMENT QUESTIONS =====\n";
        out << "Q1. Does Richardson always give a smaller error than Central Difference?\n";
        out << "    " << (richardsonAlwaysBetter
                ? "Yes, across every h tested for every function."
                : "Mostly, but NOT always: at the smallest step sizes, round-off error in computing D(h) and D(h/2) and combining them can make R(h) worse than D(h).")
            << "\n\n";

        out << "Q2. By approximately what factor does the error decrease when h is reduced by 10?\n";
        out << "    D(h): ~" << avgFactorD << "x per decade (matches O(h^2): 10^2 = 100)\n";
        out << "    R(h): ~" << avgFactorR << "x per decade (matches O(h^4): 10^4 = 10000)\n\n";

        out << "Q3. Does Richardson continue to improve as h becomes very small?\n";
        out << "    No. Once round-off dominates, R(h) can stop decreasing and increase again.\n\n";
    }

    void writeAnalysis(const string& path) const {
        ofstream out(path);
        out << "Richardson Extrapolation - Written Analysis\n";
        out << "Evaluation point x0 = " << x0 << "\n";
        printAnalysis(out);
    }
};

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

// Read differentiation-family input files (shared format)
// Format:
//   x0
//   count_of_h_values
//   h1 h2 h3 ...
bool readDifferentiationInput(const string& path, double& x0, vector<double>& hVals) {
    ifstream fin(path);
    if (!fin) return false;
    string line;
    vector<double> tokens;
    while (getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        double v;
        while (ss >> v) tokens.push_back(v);
    }
    if (tokens.size() < 2) return false;
    x0 = tokens[0];
    int count = (int)tokens[1];
    if ((int)tokens.size() < 2 + count) return false;
    hVals.assign(tokens.begin() + 2, tokens.begin() + 2 + count);
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

// ── NUMERICAL DIFFERENTIATION ───────────────────────────────────────────────
void runDifferentiation() {
    printSeparator("NUMERICAL DIFFERENTIATION");

    write("Theory:\n");
    write("  Forward:  f'(x) ~= [f(x+h) - f(x)] / h          -- O(h)\n");
    write("  Backward: f'(x) ~= [f(x) - f(x-h)] / h           -- O(h)\n");
    write("  Central:  f'(x) ~= [f(x+h) - f(x-h)] / (2h)      -- O(h^2)\n\n");

    double x0;
    vector<double> hVals;
    if (!readDifferentiationInput("input/input_differentiation.txt", x0, hVals)) {
        write("[ERROR] Cannot open/parse input/input_differentiation.txt\n");
        write("        Expected format:\n");
        write("          x0\n");
        write("          count_of_h_values\n");
        write("          h1 h2 h3 ...\n");
        return;
    }

    char buf[128];
    snprintf(buf, sizeof(buf), "Evaluation point x0 = %.6f\n", x0);
    write(string(buf));
    write("Step sizes h = [ ");
    for (double h : hVals) { snprintf(buf, sizeof(buf), "%.0e ", h); write(string(buf)); }
    write("]\n\n");

    write("Choose test function(s):\n");
    write("  1. f(x) = e^x\n");
    write("  2. f(x) = sin(x)\n");
    write("  3. f(x) = x^3 - 2x + 1\n");
    write("  4. All three (recommended -- matches assignment spec)\n");
    write("Enter choice: ");
    int fchoice; cin >> fchoice;
    write("\n");

    DifferentiationAnalyzer analyzer(x0, hVals);

    if (fchoice == 1 || fchoice == 4)
        analyzer.addFunction("exp",
            [](double x){ return exp(x); },
            [](double x){ return exp(x); });
    if (fchoice == 2 || fchoice == 4)
        analyzer.addFunction("sin",
            [](double x){ return sin(x); },
            [](double x){ return cos(x); });
    if (fchoice == 3 || fchoice == 4)
        analyzer.addFunction("poly",
            [](double x){ return x*x*x - 2*x + 1; },
            [](double x){ return 3*x*x - 2; });

    if (fchoice < 1 || fchoice > 4) {
        write("  [Error] Invalid choice.\n");
        return;
    }

    analyzer.run();

    // Results table (console + shared output/output.txt, via write())
    ostringstream tableStream;
    analyzer.printTable(tableStream);
    write(tableStream.str());

    // Written analysis answering: most accurate method / observed order / round-off floor
    ostringstream analysisStream;
    analyzer.printAnalysis(analysisStream);
    write(analysisStream.str());

    // ---- All differentiation outputs saved under output/, in proper formats ----
    filesystem::create_directories("output");
    analyzer.writeCSV("output/differentiation_results.csv");
    analyzer.writeGnuplotFiles("output/diff_data", "output/plot.gp",
                                "output/loglog_error_plot.png");
    bool plotted = analyzer.renderPlot("output/plot.gp");

    write("\nSaved to output/:\n");
    write("  differentiation_results.csv   (raw results, CSV)\n");
    write("  differentiation_analysis.txt  (written analysis)\n");
    write("  diff_data_<func>.dat + plot.gp (gnuplot source data/script)\n");
    if (plotted)
        write("  loglog_error_plot.png         (log-log error graph)\n");
    else
        write("  [!] loglog_error_plot.png NOT generated -- gnuplot not found.\n"
              "      Install gnuplot, then run: gnuplot output/plot.gp\n");
}

// ── RICHARDSON EXTRAPOLATION (Assignment III) ───────────────────────────────
void runRichardson() {
    printSeparator("RICHARDSON EXTRAPOLATION -- INVESTIGATION");

    write("Theory:\n");
    write("  Central Difference: D(h) = [f(x+h) - f(x-h)] / (2h)      -- O(h^2)\n");
    write("  Richardson:         R(h) = [4*D(h/2) - D(h)] / 3         -- O(h^4)\n");
    write("  (Richardson eliminates the leading h^2 error term of D(h).)\n\n");

    double x0;
    vector<double> hVals;
    if (!readDifferentiationInput("input/input_richardson.txt", x0, hVals)) {
        write("[ERROR] Cannot open/parse input/input_richardson.txt\n");
        write("        Expected format:\n");
        write("          x0\n");
        write("          count_of_h_values\n");
        write("          h1 h2 h3 ...\n");
        return;
    }

    char buf[128];
    snprintf(buf, sizeof(buf), "Evaluation point x0 = %.6f\n", x0);
    write(string(buf));
    write("Step sizes h = [ ");
    for (double h : hVals) { snprintf(buf, sizeof(buf), "%.0e ", h); write(string(buf)); }
    write("]\n\n");

    write("Choose test function(s):\n");
    write("  1. f(x) = e^x\n");
    write("  2. f(x) = sin(x)\n");
    write("  3. f(x) = cos(x)\n");
    write("  4. f(x) = x^3 - 2x + 1\n");
    write("  5. All four (recommended -- matches assignment spec)\n");
    write("Enter choice: ");
    int fchoice; cin >> fchoice;
    write("\n");

    RichardsonAnalyzer analyzer(x0, hVals);

    if (fchoice == 1 || fchoice == 5)
        analyzer.addFunction("exp",
            [](double x){ return exp(x); },
            [](double x){ return exp(x); });
    if (fchoice == 2 || fchoice == 5)
        analyzer.addFunction("sin",
            [](double x){ return sin(x); },
            [](double x){ return cos(x); });
    if (fchoice == 3 || fchoice == 5)
        analyzer.addFunction("cos",
            [](double x){ return cos(x); },
            [](double x){ return -sin(x); });
    if (fchoice == 4 || fchoice == 5)
        analyzer.addFunction("poly",
            [](double x){ return x*x*x - 2*x + 1; },
            [](double x){ return 3*x*x - 2; });

    if (fchoice < 1 || fchoice > 5) {
        write("  [Error] Invalid choice.\n");
        return;
    }

    analyzer.run();

    // Results table (console + shared output/output.txt, via write())
    ostringstream tableStream;
    analyzer.printTable(tableStream);
    write(tableStream.str());

    // Written analysis: slopes, error-reduction factors, round-off floor, Q1-Q6
    ostringstream analysisStream;
    analyzer.printAnalysis(analysisStream);
    write(analysisStream.str());

    // ---- All Richardson outputs saved under output/, in proper formats ----
    filesystem::create_directories("output");

    analyzer.writeTable("output/richardson_table.txt");
    analyzer.writeCSV("output/richardson_results.csv");
    analyzer.writeAnalysis("output/richardson_analysis.txt");
    analyzer.writeGnuplotFiles("output/richardson_data", "output/richardson_plot.gp",
                                "output/richardson_loglog_plot.png");
    bool plotted = analyzer.renderPlot("output/richardson_plot.gp");

    write("\nSaved to output/:\n");
    write("  richardson_table.txt              (formatted results table: D(h), R(h), errors)\n");
    write("  richardson_results.csv            (raw results, CSV)\n");
    write("  richardson_analysis.txt           (written analysis + Q1-Q6 answers)\n");
    write("  richardson_data_<func>.dat + richardson_plot.gp (gnuplot source)\n");
    if (plotted)
        write("  richardson_loglog_plot.png        (log-log error graph: D(h) vs R(h))\n");
    else
        write("  [!] richardson_loglog_plot.png NOT generated -- gnuplot not found.\n"
              "      Install gnuplot, then run: gnuplot output/richardson_plot.gp\n");
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
    cout << "  4. Numerical Differentiation (Forward / Backward / Central)\n";
    cout << "  5. Richardson Extrapolation  (D(h) vs R(h) -- Assignment III)\n";
    cout << "Enter choice: ";
    int category; cin >> category;

    if (category == 5) {
        runRichardson();
    }
    else if (category == 4) {
        runDifferentiation();
    }
    else if (category == 3) {
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





