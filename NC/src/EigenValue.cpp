#include "../include/EigenValue.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
using namespace std;

EigenValue::EigenValue(int n) : Matrix(n, n) {}

EigenValue::EigenValue(const Matrix& m) : Matrix(m) {
    if (rows != cols)
        throw invalid_argument("[EigenValue] Matrix must be square (n×n).");
}

vector<vector<double>> EigenValue::toRaw() const {
    vector<vector<double>> A(rows, vector<double>(cols));
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            A[i][j] = data[i][j];
    return A;
}

double EigenValue::vecNorm(const vector<double>& v) const {
    double s = 0.0;
    for (double x : v) s += x * x;
    return sqrt(s);
}

void EigenValue::vecNormalize(vector<double>& v) const {
    double n = vecNorm(v);
    if (n < 1e-15) return;   // zero vector — leave as is
    for (double& x : v) x /= n;
}

double EigenValue::vecDot(const vector<double>& a,
                          const vector<double>& b) const {
    double s = 0.0;
    for (int i = 0; i < (int)a.size(); i++) s += a[i] * b[i];
    return s;
}

vector<double> EigenValue::matvec(const vector<vector<double>>& A,
                                  const vector<double>& v) const {
    int n = (int)v.size();
    vector<double> w(n, 0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            w[i] += A[i][j] * v[j];
    return w;
}

void EigenValue::qrDecompose(const vector<vector<double>>& A,
                             vector<vector<double>>& Q,
                             vector<vector<double>>& R) const {
    int n = (int)A.size();
    Q.assign(n, vector<double>(n, 0.0));
    R.assign(n, vector<double>(n, 0.0));

    // Extract columns of A into cols[j]
    vector<vector<double>> cols(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            cols[j][i] = A[i][j];
    for (int j = 0; j < n; j++) {
        R[j][j] = vecNorm(cols[j]);
        if (R[j][j] < 1e-14) {
            cols[j].assign(n, 0.0);
            cols[j][j] = 1.0;
            R[j][j]    = 1.0;
        }
        // q_j = cols[j] / R[j][j]
        vector<double> qj = cols[j];
        for (double& x : qj) x /= R[j][j];
        for (int i = 0; i < n; i++) Q[i][j] = qj[i];

        // Project out q_j from remaining columns
        for (int k = j + 1; k < n; k++) {
            R[j][k] = vecDot(qj, cols[k]);
            for (int i = 0; i < n; i++) cols[k][i] -= R[j][k] * qj[i];
        }
    }
}

vector<double> EigenValue::gaussSolve(vector<vector<double>> M,
                                      vector<double> rhs) const {
    int n = (int)rhs.size();

    // Forward elimination with partial pivoting
    for (int col = 0; col < n; col++) {
        // Find pivot row
        int pivot = col;
        for (int row = col + 1; row < n; row++)
            if (fabs(M[row][col]) > fabs(M[pivot][col])) pivot = row;

        swap(M[col], M[pivot]);
        swap(rhs[col], rhs[pivot]);

        if (fabs(M[col][col]) < 1e-15)
            throw runtime_error("[gaussSolve] Singular matrix encountered.");

        for (int row = col + 1; row < n; row++) {
            double factor = M[row][col] / M[col][col];
            for (int k = col; k < n; k++)
                M[row][k] -= factor * M[col][k];
            rhs[row] -= factor * rhs[col];
        }
    }

    // Back substitution
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = rhs[i];
        for (int j = i + 1; j < n; j++)
            x[i] -= M[i][j] * x[j];
        x[i] /= M[i][i];
    }
    return x;
}

double EigenValue::powerMethod(vector<double>& eigenvec,
                               int maxIter,
                               double tol) const {
    int n = rows;
    auto A = toRaw();

    // Initial vector — ones, then normalise
    eigenvec.assign(n, 1.0);
    vecNormalize(eigenvec);

    double lambda     = 0.0;
    double lambda_old = 1e30;

    for (int iter = 0; iter < maxIter; iter++) {
        vector<double> w = matvec(A, eigenvec);

        lambda = vecDot(eigenvec, w);   // v normalised → denominator = 1

        double change = fabs(lambda - lambda_old);
        lambda_old = lambda;

        vecNormalize(w);
        eigenvec = w;

        if (change < tol) {
            cout << "  (Power Method converged in " << iter + 1
                 << " iterations)\n";
            return lambda;
        }
    }
    cout << "  [Warning] Power Method reached max iterations ("
         << maxIter << "). Result may be approximate.\n";
    return lambda;
}

double EigenValue::inversePowerMethod(vector<double>& eigenvec,
                                     int maxIter,
                                     double tol) const {
    int n = rows;
    auto A = toRaw();

    // Start with a random-ish vector (not all ones to avoid symmetry bias)
    eigenvec.resize(n);
    for (int i = 0; i < n; i++) eigenvec[i] = 1.0 + 0.1 * i;
    vecNormalize(eigenvec);

    vector<double> v_old = eigenvec;

    for (int iter = 0; iter < maxIter; iter++) {
        vector<double> w;
        try {
            w = gaussSolve(A, eigenvec);    // solve A·w = v  (= A⁻¹ v)
        } catch (...) {
            cout << "  [Warning] Matrix is singular; λmin ≈ 0.\n";
            return 0.0;
        }

        vecNormalize(w);

        // Convergence check: eigenvector stopped changing
        double diff = 0.0;
        for (int i = 0; i < n; i++)
            diff = max(diff, fabs(fabs(w[i]) - fabs(eigenvec[i])));

        eigenvec = w;

        if (diff < tol) {
            // Compute eigenvalue: λ = v' A v  (Rayleigh quotient of A)
            vector<double> Av = matvec(A, eigenvec);
            double lambda = vecDot(eigenvec, Av);
            cout << "  (Inverse Power Method converged in " << iter + 1
                 << " iterations)\n";
            return lambda;
        }
    }

    // Final Rayleigh quotient of A
    vector<double> Av = matvec(A, eigenvec);
    double lambda = vecDot(eigenvec, Av);
    cout << "  (Inverse Power Method: " << maxIter
         << " iters, returning Rayleigh estimate)\n";
    return lambda;
}

vector<double> EigenValue::qrAlgorithm(int maxIter, double tol) const {
    int n = rows;
    auto Ak = toRaw();

    for (int iter = 0; iter < maxIter; iter++) {
        double a = Ak[n-2][n-2], b = Ak[n-2][n-1];
        double c = Ak[n-1][n-2], d = Ak[n-1][n-1];
        double tr   = a + d;
        double disc = (a - d) * (a - d) + 4.0 * b * c;
        double shift;
        if (disc >= 0.0) {
            double sq = sqrt(disc);
            double e1 = (tr + sq) / 2.0;
            double e2 = (tr - sq) / 2.0;
            shift = (fabs(e1 - d) < fabs(e2 - d)) ? e1 : e2;
        } else {
            shift = d;
        }

        // ── Shift: A_k − σI ──────────────────────────────────
        for (int i = 0; i < n; i++) Ak[i][i] -= shift;

        // ── QR decomposition of Ak ────────────────────────────
        vector<vector<double>> Q, R;
        qrDecompose(Ak, Q, R);

        // ── A_{k+1} = R·Q ─────────────────────────────────────
        vector<vector<double>> Anew(n, vector<double>(n, 0.0));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int k = 0; k < n; k++)
                    Anew[i][j] += R[i][k] * Q[k][j];

        // ── Un-shift: A_{k+1} + σI ────────────────────────────
        for (int i = 0; i < n; i++) Anew[i][i] += shift;
        Ak = Anew;

        // ── Convergence check: sub-diagonal elements → 0 ─────
        double offDiag = 0.0;
        for (int i = 1; i < n; i++) offDiag += fabs(Ak[i][i-1]);
        if (offDiag < tol) {
            cout << "  (QR Algorithm converged in " << iter + 1
                 << " iterations)\n";
            break;
        }
    }

    // Diagonal of Ak ≈ eigenvalues
    vector<double> eigenvalues(n);
    for (int i = 0; i < n; i++) eigenvalues[i] = Ak[i][i];

    // Sort descending by magnitude
    sort(eigenvalues.begin(), eigenvalues.end(),
         [](double a, double b){ return fabs(a) > fabs(b); });

    return eigenvalues;
}

vector<EigenValue::GershgorinDisk> EigenValue::gershgorinBounds() const {
    int n = rows;
    vector<GershgorinDisk> disks(n);

    for (int i = 0; i < n; i++) {
        disks[i].row    = i;
        disks[i].center = data[i][i];
        disks[i].radius = 0.0;
        for (int j = 0; j < n; j++)
            if (j != i) disks[i].radius += fabs(data[i][j]);
    }
    return disks;
}

double EigenValue::conditionNumber() const {
    vector<double> eigs = qrAlgorithm();
    if (eigs.empty()) return 1.0;

    double lmax = 0.0, lmin = 1e30;
    for (double e : eigs) {
        double ae = fabs(e);
        if (ae > lmax) lmax = ae;
        if (ae < lmin) lmin = ae;
    }
    if (lmin < 1e-15) return 1e15;   // singular matrix
    return lmax / lmin;
}