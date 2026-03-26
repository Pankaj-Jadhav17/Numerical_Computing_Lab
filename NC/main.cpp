#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include "include/Matrix.hpp"
#include "include/GaussElimination.hpp"
#include "include/GaussJacobi.hpp"
#include "include/LUCrout.hpp"
#include "include/LUDoolittle.hpp"
#include "include/LUCholesky.hpp"
#include "include/GaussSeidel.hpp"

using namespace std;

ofstream fout;

// OUTPUT 
void write(const string& s) {
    cout << s;
    fout << s;
}

// INPUT 
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

// PRINT 
void printMatrix(const vector<vector<double>>& M, int n, const string& name) {
    write("\n" + name + "\n");
    for (int i = 0; i < n; i++) {
        write("  ");
        for (int j = 0; j < n; j++) {
            char buf[32];
            snprintf(buf, sizeof(buf), "%10.4f", M[i][j]);
            write(string(buf));
        }
        write("\n");
    }
}

void printSolution(const vector<double>& x) {
    write("\nSolution:\n");
    for (int i = 0; i < (int)x.size(); i++) {
        char buf[64];
        snprintf(buf, sizeof(buf), "  x%d = %.8f\n", i+1, x[i]);
        write(string(buf));
    }
}

//  SOLVER 
template<typename Solver>
void fillSolver(Solver& s,
                const vector<vector<double>>& A,
                const vector<double>& b, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            s(i,j) = A[i][j];
        s(i,n) = b[i];
    }
}

vector<double> solveSystem(const vector<vector<double>>& A,
                           const vector<double>& b,
                           int method, int n) {

    vector<double> x;

    switch (method) {
        case 1: {
            GaussElimination s(n,n+1);
            fillSolver(s,A,b,n);
            x = s.solve();
            break;
        }
        case 2: {
            GaussJacobi s(n,n+1);
            fillSolver(s,A,b,n);
            x = s.solve();
            break;
        }
        case 3: {
            LUCrout s(n,n+1);
            fillSolver(s,A,b,n);
            x = s.solve();
            break;
        }
        case 4: {
            LUDoolittle s(n,n+1);
            fillSolver(s,A,b,n);
            x = s.solve();
            break;
        }
        case 5: {
            LUCholesky s(n,n+1);
            fillSolver(s,A,b,n);
            x = s.solve();
            break;
        }
        case 6: {
            GaussSeidel s(n,n+1);
            fillSolver(s,A,b,n);
            x = s.solve();
            break;
        }
    }
    return x;
}

// RUN SYSTEM 
void runSolve(const string& title,
              const vector<vector<double>>& A,
              const vector<double>& b,
              int n) {

    write("\n==============================\n");
    write(title);
    write("\n==============================\n");

    printMatrix(A, n, "Matrix A");

    write("\nChoose Method:\n");
    write("1. Gauss Elimination\n");
    write("2. Gauss Jacobi\n");
    write("3. LU Crout\n");
    write("4. LU Doolittle\n");
    write("5. LU Cholesky\n");
    write("6. Gauss Seidel\n");

    int method;
    cin >> method;

    vector<double> x = solveSystem(A, b, method, n);
    printSolution(x);
}

//  MAIN 
int main() {

    vector<vector<double>> A, Atilde;
    vector<double> b, btilde;
    int nA, nAt;

    if (!readMatrix("input/input_A.txt", A, nA) ||
        !readMatrix("input/input_Atilde.txt", Atilde, nAt) ||
        !readVector("input/input_b.txt", b, nA) ||
        !readVector("input/input_btilde.txt", btilde, nAt)) {

        cout << "[ERROR] Missing input files!\n";
        return 1;
    }

    fout.open("output/output.txt");

    // SELECT MODE 
    cout << "\nSelect Mode:\n";
    cout << "1. Work with A\n";
    cout << "2. Work with A~\n";
    cout << "3. Run all 4 systems\n";
    cout << "Enter choice: ";

    int mode;
    cin >> mode;

    //  OPERATION 
    cout << "\nSelect Operation:\n";
    cout << "1. Determinant\n";
    cout << "2. Inverse\n";
    cout << "3. Solve System\n";
    cout << "Enter choice: ";

    int op;
    cin >> op;

    // CASE 1 & 2 
    if (mode == 1 || mode == 2) {

        vector<vector<double>> Ause = (mode == 1) ? A : Atilde;
        vector<double> buse = (mode == 1) ? b : btilde;
        int n = (mode == 1) ? nA : nAt;

        Matrix mat(n,n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                mat(i,j) = Ause[i][j];

        printMatrix(Ause, n, "Matrix");

        if (op == 1) {
            double det = mat.determinant();
            write("\nDeterminant = " + to_string(det) + "\n");
        }
        else if (op == 2) {
            try {
                Matrix inv = mat.inverse();
                write("\nInverse Matrix:\n");
                for (int i = 0; i < n; i++) {
                    write("  ");
                    for (int j = 0; j < n; j++) {
                        char buf[32];
                        snprintf(buf,sizeof(buf),"%10.4f",inv(i,j));
                        write(string(buf));
                    }
                    write("\n");
                }
            } catch (...) {
                write("\nMatrix is singular\n");
            }
        }
        else if (op == 3) {
            runSolve("Solving System", Ause, buse, n);
        }
    }

    // ALL 4 CASES 
    else if (mode == 3 && op == 3) {
        runSolve("1. A x = b", A, b, nA);
        runSolve("2. A~ x = b", Atilde, b, nAt);
        runSolve("3. A x = b~", A, btilde, nA);
        runSolve("4. A~ x = b~", Atilde, btilde, nAt);
    }

    fout.close();
    cout << "\nOutput saved to output/output.txt\n";

    return 0;
}















// #include <iostream>
// #include <fstream>
// #include <vector>
// #include <string>
// #include <cmath>

// #include "include/Matrix.hpp"
// #include "include/GaussElimination.hpp"
// #include "include/GaussJacobi.hpp"
// #include "include/LUCrout.hpp"
// #include "include/LUDoolittle.hpp"
// #include "include/LUCholesky.hpp"

// using namespace std;

// ofstream fout;

// // ================= OUTPUT =================
// void write(const string& s) {
//     cout << s;
//     fout << s;
// }

// // ================= INPUT =================
// bool readMatrix(const string& path, vector<vector<double>>& M, int& n) {
//     ifstream fin(path);
//     if (!fin) return false;

//     int r, c;
//     fin >> r >> c;
//     n = r;

//     M.assign(r, vector<double>(c));
//     for (auto& row : M)
//         for (auto& val : row)
//             fin >> val;

//     return true;
// }

// bool readVector(const string& path, vector<double>& b, int n) {
//     ifstream fin(path);
//     if (!fin) return false;

//     b.resize(n);
//     for (auto& val : b) fin >> val;

//     return true;
// }

// // ================= PRINT =================
// void printMatrix(const vector<vector<double>>& M, int n, const string& name) {
//     write("\n" + name + "\n");
//     for (int i = 0; i < n; i++) {
//         write("  ");
//         for (int j = 0; j < n; j++) {
//             char buf[32];
//             snprintf(buf, sizeof(buf), "%10.4f", M[i][j]);
//             write(string(buf));
//         }
//         write("\n");
//     }
// }

// void printSolution(const vector<double>& x) {
//     write("\nSolution:\n");
//     for (int i = 0; i < (int)x.size(); i++) {
//         char buf[64];
//         snprintf(buf, sizeof(buf), "  x%d = %.8f\n", i+1, x[i]);
//         write(string(buf));
//     }
// }

// // ================= SOLVER =================
// template<typename Solver>
// void fillSolver(Solver& s,
//                 const vector<vector<double>>& A,
//                 const vector<double>& b, int n) {
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++)
//             s(i,j) = A[i][j];
//         s(i,n) = b[i];
//     }
// }

// vector<double> solveSystem(const vector<vector<double>>& A,
//                            const vector<double>& b,
//                            int method, int n) {

//     vector<double> x;

//     switch (method) {
//         case 1: {
//             GaussElimination s(n,n+1);
//             fillSolver(s,A,b,n);
//             x = s.solve();
//             break;
//         }
//         case 2: {
//             GaussJacobi s(n,n+1);
//             fillSolver(s,A,b,n);
//             x = s.solve();
//             break;
//         }
//         case 3: {
//             LUCrout s(n,n+1);
//             fillSolver(s,A,b,n);
//             x = s.solve();
//             break;
//         }
//         case 4: {
//             LUDoolittle s(n,n+1);
//             fillSolver(s,A,b,n);
//             x = s.solve();
//             break;
//         }
//         case 5: {
//             LUCholesky s(n,n+1);
//             fillSolver(s,A,b,n);
//             x = s.solve();
//             break;
//         }
//     }
//     return x;
// }

// // ================= RUN SYSTEM =================
// void runSolve(const string& title,
//               const vector<vector<double>>& A,
//               const vector<double>& b,
//               int n) {

//     write("\n==============================\n");
//     write(title);
//     write("\n==============================\n");

//     printMatrix(A, n, "Matrix A");

//     write("\nChoose Method:\n");
//     write("1. Gauss Elimination\n");
//     write("2. Gauss Jacobi\n");
//     write("3. LU Crout\n");
//     write("4. LU Doolittle\n");
//     write("5. LU Cholesky\n");

//     int method;
//     cin >> method;

//     vector<double> x = solveSystem(A, b, method, n);
//     printSolution(x);
// }

// // ================= MAIN =================
// int main() {

//     vector<vector<double>> A, Atilde;
//     vector<double> b, btilde;
//     int nA, nAt;

//     if (!readMatrix("input/input_A.txt", A, nA) ||
//         !readMatrix("input/input_Atilde.txt", Atilde, nAt) ||
//         !readVector("input/input_b.txt", b, nA) ||
//         !readVector("input/input_btilde.txt", btilde, nAt)) {

//         cout << "[ERROR] Missing input files!\n";
//         return 1;
//     }

//     fout.open("output/output.txt");

//     // ================= SELECT MODE =================
//     cout << "\nSelect Mode:\n";
//     cout << "1. Work with A\n";
//     cout << "2. Work with A~\n";
//     cout << "3. Run all 4 systems\n";
//     cout << "Enter choice: ";

//     int mode;
//     cin >> mode;

//     // ================= OPERATION =================
//     cout << "\nSelect Operation:\n";
//     cout << "1. Determinant\n";
//     cout << "2. Inverse\n";
//     cout << "3. Solve System\n";
//     cout << "Enter choice: ";

//     int op;
//     cin >> op;

//     // ================= CASE 1 & 2 =================
//     if (mode == 1 || mode == 2) {

//         vector<vector<double>> Ause = (mode == 1) ? A : Atilde;
//         vector<double> buse = (mode == 1) ? b : btilde;
//         int n = (mode == 1) ? nA : nAt;

//         Matrix mat(n,n);
//         for (int i = 0; i < n; i++)
//             for (int j = 0; j < n; j++)
//                 mat(i,j) = Ause[i][j];

//         printMatrix(Ause, n, "Matrix");

//         if (op == 1) {
//             double det = mat.determinant();
//             write("\nDeterminant = " + to_string(det) + "\n");
//         }
//         else if (op == 2) {
//             try {
//                 Matrix inv = mat.inverse();
//                 write("\nInverse Matrix:\n");
//                 for (int i = 0; i < n; i++) {
//                     write("  ");
//                     for (int j = 0; j < n; j++) {
//                         char buf[32];
//                         snprintf(buf,sizeof(buf),"%10.4f",inv(i,j));
//                         write(string(buf));
//                     }
//                     write("\n");
//                 }
//             } catch (...) {
//                 write("\nMatrix is singular\n");
//             }
//         }
//         else if (op == 3) {
//             runSolve("Solving System", Ause, buse, n);
//         }
//     }

//     // ================= ALL 4 CASES =================
//     else if (mode == 3 && op == 3) {
//         runSolve("1. A x = b", A, b, nA);
//         runSolve("2. A~ x = b", Atilde, b, nAt);
//         runSolve("3. A x = b~", A, btilde, nA);
//         runSolve("4. A~ x = b~", Atilde, btilde, nAt);
//     }

//     fout.close();
//     cout << "\nOutput saved to output/output.txt\n";

//     return 0;
// }
