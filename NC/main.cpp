#include <iostream>
#include <fstream>
#include <vector>
#include "include/Matrix.hpp"
#include "include/GaussElimination.hpp"
#include "include/GaussJacobi.hpp"
#include "include/LUCrout.hpp"
#include "include/LUDoolittle.hpp"
#include "include/LUCholesky.hpp"

using namespace std;

int main()
{
    ifstream fin("input/input.txt");
    ofstream fout("output/output.txt");

    if(!fin)
    {
        cout << "Input file not found\n";
        return 0;
    }

    int r, c;
    fin >> r >> c;

    cout << "\nChoose Method:\n";
    cout << "1. Gauss Elimination\n";
    cout << "2. Gauss Jacobi\n";
    cout << "3. LU Crout\n";
    cout << "4. LU Doolittle\n";
    cout << "5. LU Cholesky\n";
    cout << "6. Determinant\n";        
    cout << "7. Inverse\n";            
    cout << "8. Transpose\n";          
    cout << "Enter choice: ";

    int choice;
    cin >> choice;

    vector<double> solution;
    switch(choice)
    {
        case 1:
        {
            GaussElimination solver(r, c);
            fin >> solver;
            fout << "Using Gauss Elimination\n";
            fout << "Input Matrix:\n" << solver;
            solution = solver.solve();
            break;
        }
        case 2:
        {
            GaussJacobi solver(r, c);
            fin >> solver;
            fout << "Using Gauss Jacobi\n";
            fout << "Input Matrix:\n" << solver;
            solution = solver.solve();
            break;
        }
        case 3:
        {
            LUCrout solver(r, c);
            fin >> solver;
            fout << "Using LU Crout\n";
            fout << "Input Matrix:\n" << solver;
            solution = solver.solve();
            break;
        }
        case 4:
        {
            LUDoolittle solver(r, c);
            fin >> solver;
            fout << "Using LU Doolittle\n";
            fout << "Input Matrix:\n" << solver;
            solution = solver.solve();
            break;
        }
        case 5:
        {
            LUCholesky solver(r, c);
            fin >> solver;
            fout << "Using LU Cholesky\n";
            fout << "Input Matrix:\n" << solver;
            solution = solver.solve();
            break;
        }
        case 6:
        {
            Matrix A(r, r);       // square matrix only
            fin >> A;

            double det = A.determinant();

            fout << "Operation: Determinant\n";
            fout << "Input Matrix:\n" << A;
            fout << "\nDeterminant = " << det << "\n";

            cout << "\nDeterminant = " << det << "\n";
            break;
        }
        case 7:
        {
            Matrix A(r, r);       // square matrix only
            fin >> A;

            fout << "Operation: Inverse\n";
            fout << "Input Matrix:\n" << A;

            try
            {
                Matrix inv = A.inverse();
                fout << "\nInverse Matrix:\n" << inv;
                cout << "\nInverse computed successfully.\n";
            }
            catch(exception& e)
            {
                fout << "\nError: " << e.what() << "\n";
                cout << "\nError: " << e.what() << "\n";
            }
            break;
        }
        case 8:
        {
            Matrix A(r, c);
            fin >> A;

            Matrix T = A.transpose();

            fout << "Operation: Transpose\n";
            fout << "Input Matrix:\n" << A;
            fout << "\nTranspose Matrix:\n" << T;

            cout << "\nTranspose computed successfully.\n";
            break;
        }
        default:
        {
            cout << "Invalid choice\n";
            return 0;
        }
    }
    // Print solution (only for cases 1-5)
    if(!solution.empty())
    {
        fout << "\nSolution:\n";
        for(size_t i = 0; i < solution.size(); i++)
        {
            fout << "x" << i+1 << " = " << solution[i] << endl;
        }
    }

    fin.close();
    fout.close();

    cout << "\nProgram executed successfully.\n";
    cout << "Check output/output.txt\n";

    return 0;
}












































































































// #include <iostream>
// #include <fstream>
// #include <iomanip>
// #include <string>
// #include <vector>
// #include <limits>
// #include <ctime>
// #include <chrono>
// #include <stdexcept>
// #include <cmath>

// #include "include/Matrix.hpp"
// #include "include/GaussElimination.hpp"
// #include "include/GaussJacobi.hpp"
// #include "include/LUCrout.hpp"
// #include "include/LUDoolittle.hpp"
// #include "include/LUCholesky.hpp"

// using namespace std;

// //  UI HELPERS
// // ================================================================

// void line(char c = '=', int n = 60) { cout << string(n, c) << "\n"; }

// void header(const string& t) {
//     cout << "\n"; line('=');
//     int p = max(0, (60 - (int)t.size()) / 2);
//     cout << string(p,' ') << t << "\n";
//     line('=');
// }

// void section(const string& t) {
//     cout << "\n"; line('-', 55);
//     cout << "  " << t << "\n";
//     line('-', 55);
// }

// // ── validated integer input ──────────────────────────────────────
// int getInt(const string& prompt, int lo, int hi) {
//     int v;
//     while (true) {
//         cout << prompt;
//         if (cin >> v && v >= lo && v <= hi) return v;
//         cout << "  Enter a number between " << lo << " and " << hi << ".\n";
//         cin.clear();
//         cin.ignore(numeric_limits<streamsize>::max(), '\n');
//     }
// }

// bool yesNo(const string& prompt) {
//     char c;
//     while (true) {
//         cout << prompt << " (y/n): ";
//         cin >> c; c = tolower(c);
//         if (c == 'y') return true;
//         if (c == 'n') return false;
//         cout << "  Enter y or n.\n";
//     }
// }

// // ── pretty print a Matrix object ────────────────────────────────
// void printMatrix(const Matrix& M, int r, int c, ostream& out = cout) {
//     out << "\n";
//     for (int i = 0; i < r; i++) {
//         out << "  |";
//         for (int j = 0; j < c; j++)
//             out << setw(12) << fixed << setprecision(4) << M(i, j);
//         out << "  |\n";
//     }
//     out << "\n";
// }

// // ================================================================
// //  ENTER A SQUARE MATRIX (n x n) FROM KEYBOARD OR FILE
// // ================================================================
// Matrix inputSquareMatrix(int& n, const string& label = "Matrix") {
//     section("INPUT: " + label);
//     cout << "\n  [1] Enter manually   [2] Load from file\n\n";
//     int ch = getInt("  Choice: ", 1, 2);

//     if (ch == 1) {
//         n = getInt("\n  Size (NxN, 1-8): ", 1, 8);
//         Matrix M(n, n);
//         cout << "\n  Enter " << n << "x" << n << " matrix row by row:\n\n";
//         for (int i = 0; i < n; i++) {
//             cout << "  Row " << i+1 << " [" << n << " values]: ";
//             for (int j = 0; j < n; j++) cin >> M(i, j);
//         }
//         return M;
//     } else {
//         cin.ignore(numeric_limits<streamsize>::max(), '\n');
//         cout << "\n  File path [Enter = input/input.txt]: ";
//         string path; getline(cin, path);
//         if (path.empty()) path = "input/input.txt";
//         ifstream fin(path);
//         if (!fin) {
//             cout << "  File not found. Switching to manual.\n";
//             ch = 1; // fallback — redo as manual
//             n = getInt("\n  Size (NxN, 1-8): ", 1, 8);
//             Matrix M(n, n);
//             cout << "\n  Enter " << n << "x" << n << " matrix:\n";
//             for (int i = 0; i < n; i++) {
//                 cout << "  Row " << i+1 << ": ";
//                 for (int j = 0; j < n; j++) cin >> M(i, j);
//             }
//             return M;
//         }
//         int r, c; fin >> r >> c;
//         n = r;
//         // Read all columns, but only keep the square (r x r) part
//         vector<vector<double>> tmp(r, vector<double>(c, 0.0));
//         for (int i = 0; i < r; i++)
//             for (int j = 0; j < c; j++)
//                 fin >> tmp[i][j];
//         fin.close();
//         int use_c = min(r, c);  // square part only
//         Matrix M(r, use_c);
//         for (int i = 0; i < r; i++)
//             for (int j = 0; j < use_c; j++)
//                 M(i, j) = tmp[i][j];
//         cout << "  Loaded " << r << "x" << use_c << " square matrix from " << path << "\n";
//         return M;
//     }
// }

// // ================================================================
// //  ENTER AUGMENTED MATRIX (n x n+1) FOR SLE SOLVING
// // ================================================================
// Matrix inputAugmented(int& n) {
//     section("INPUT: Augmented Matrix [A|b]");
//     cout << "\n  [1] Enter manually   [2] Load from file\n\n";
//     int ch = getInt("  Choice: ", 1, 2);

//     if (ch == 1) {
//         n = getInt("\n  Number of equations (1-8): ", 1, 8);
//         int c = n + 1;
//         Matrix M(n, c);
//         cout << "\n  Enter augmented matrix [" << n << "x" << c
//              << "] — last column is b.\n\n";
//         for (int i = 0; i < n; i++) {
//             cout << "  Row " << i+1 << " [" << c << " values]: ";
//             for (int j = 0; j < c; j++) cin >> M(i, j);
//         }
//         return M;
//     } else {
//         cin.ignore(numeric_limits<streamsize>::max(), '\n');
//         cout << "\n  File path [Enter = input/input.txt]: ";
//         string path; getline(cin, path);
//         if (path.empty()) path = "input/input.txt";
//         ifstream fin(path);
//         if (!fin) {
//             cout << "  File not found. Switching to manual.\n";
//             n = getInt("\n  Number of equations: ", 1, 8);
//             int c = n + 1;
//             Matrix M(n, c);
//             for (int i = 0; i < n; i++) {
//                 cout << "  Row " << i+1 << ": ";
//                 for (int j = 0; j < c; j++) cin >> M(i, j);
//             }
//             return M;
//         }
//         int r, c; fin >> r >> c;
//         n = r;
//         Matrix M(r, c);
//         for (int i = 0; i < r; i++)
//             for (int j = 0; j < c; j++)
//                 fin >> M(i, j);
//         fin.close();
//         cout << "  Loaded from " << path << "\n";
//         return M;
//     }
// }

// // ================================================================
// //  SAVE RESULT TO FILE (with timestamp)
// // ================================================================
// void saveToFile(const string& title, const string& body) {
//     if (!yesNo("\n  Save result to file?")) return;

//     cin.ignore(numeric_limits<streamsize>::max(), '\n');
//     cout << "  File path [Enter = output/output.txt]: ";
//     string path; getline(cin, path);
//     if (path.empty()) path = "output/output.txt";

//     time_t now = time(nullptr);
//     char ts[32];
//     strftime(ts, sizeof(ts), "%Y-%m-%d %H:%M:%S", localtime(&now));

//     ofstream fout(path, ios::app);
//     fout << "\n" << string(60,'=') << "\n";
//     fout << "  Timestamp : " << ts  << "\n";
//     fout << "  Operation : " << title << "\n";
//     fout << string(60,'-') << "\n";
//     fout << body;
//     fout << string(60,'=') << "\n";
//     fout.close();
//     cout << "  Saved to: " << path << "\n";
// }

// // helper: matrix to string for saving
// string matrixToStr(const Matrix& M, int r, int c) {
//     ostringstream ss;
//     for (int i = 0; i < r; i++) {
//         ss << "  |";
//         for (int j = 0; j < c; j++)
//             ss << setw(12) << fixed << setprecision(4) << M(i,j);
//         ss << "  |\n";
//     }
//     return ss.str();
// }

// // ================================================================
// //  OPERATIONS
// // ================================================================

// // ─── A. CHECK MATRIX PROPERTIES ────────────────────────────────
// void opProperties() {
//     header("CHECK MATRIX PROPERTIES");
//     int n;
//     Matrix M = inputSquareMatrix(n, "Matrix A");

//     section("Matrix A");
//     printMatrix(M, n, n);

//     auto yn = [](bool b) { return b ? "YES" : "NO"; };

//     section("PROPERTIES RESULT");
//     cout << "\n";
//     cout << "  Is Square            : " << yn(M.isSquare())            << "\n";
//     cout << "  Is Symmetric         : " << yn(M.isSymmetric())         << "\n";
//     cout << "  Is Identity          : " << yn(M.isIdentity())          << "\n";
//     cout << "  Is Null (Zero)       : " << yn(M.isNull())              << "\n";
//     cout << "  Is Diagonal          : " << yn(M.isDiagonal())          << "\n";
//     cout << "  Is Diagonally Dom.   : " << yn(M.isDiagonallyDominant()) << "\n";

//     // build save body
//     ostringstream body;
//     body << "\n  Input Matrix:\n" << matrixToStr(M, n, n) << "\n";
//     body << "  Is Square            : " << yn(M.isSquare())            << "\n";
//     body << "  Is Symmetric         : " << yn(M.isSymmetric())         << "\n";
//     body << "  Is Identity          : " << yn(M.isIdentity())          << "\n";
//     body << "  Is Null (Zero)       : " << yn(M.isNull())              << "\n";
//     body << "  Is Diagonal          : " << yn(M.isDiagonal())          << "\n";
//     body << "  Is Diagonally Dom.   : " << yn(M.isDiagonallyDominant()) << "\n";
//     saveToFile("Matrix Properties", body.str());
// }

// // ─── B. TRANSPOSE ───────────────────────────────────────────────
// void opTranspose() {
//     header("TRANSPOSE");
//     int n;
//     Matrix M = inputSquareMatrix(n, "Matrix A");
//     Matrix T = M.transpose();

//     section("Input Matrix A");        printMatrix(M, n, n);
//     section("Transpose of A  (A^T)"); printMatrix(T, n, n);

//     ostringstream body;
//     body << "\n  Input A:\n" << matrixToStr(M, n, n)
//          << "\n  Transpose A^T:\n" << matrixToStr(T, n, n);
//     saveToFile("Transpose", body.str());
// }

// // ─── C. DETERMINANT ─────────────────────────────────────────────
// void opDeterminant() {
//     header("DETERMINANT");
//     int n;
//     Matrix M = inputSquareMatrix(n, "Matrix A");

//     if (!M.isSquare()) {
//         cout << "\n  [!] Determinant requires a square matrix.\n";
//         return;
//     }

//     section("Input Matrix A"); printMatrix(M, n, n);

//     auto t0 = chrono::high_resolution_clock::now();
//     double det = M.determinant();
//     auto t1 = chrono::high_resolution_clock::now();
//     double ms = chrono::duration<double,milli>(t1-t0).count();

//     section("RESULT");
//     cout << "\n  Determinant of A  =  "
//          << fixed << setprecision(8) << det << "\n";
//     cout << "  Time: " << fixed << setprecision(4) << ms << " ms\n";

//     if (fabs(det) < 1e-12)
//         cout << "\n  [!] Matrix is SINGULAR (det = 0). Inverse does not exist.\n";
//     else
//         cout << "\n  [OK] Matrix is NON-SINGULAR. Inverse exists.\n";

//     ostringstream body;
//     body << "\n  Input A:\n" << matrixToStr(M, n, n)
//          << "\n  Determinant = " << fixed << setprecision(8) << det << "\n"
//          << "  Time: " << ms << " ms\n";
//     saveToFile("Determinant", body.str());
// }

// // ─── D. INVERSE ─────────────────────────────────────────────────
// void opInverse() {
//     header("INVERSE");
//     int n;
//     Matrix M = inputSquareMatrix(n, "Matrix A");

//     section("Input Matrix A"); printMatrix(M, n, n);

//     try {
//         auto t0 = chrono::high_resolution_clock::now();
//         Matrix inv = M.inverse();
//         auto t1 = chrono::high_resolution_clock::now();
//         double ms = chrono::duration<double,milli>(t1-t0).count();

//         section("Inverse of A  (A^-1)"); printMatrix(inv, n, n);
//         cout << "  Time: " << fixed << setprecision(4) << ms << " ms\n";

//         // Verify: A * A^-1 should be Identity
//         Matrix check = M * inv;
//         cout << "\n  Verification  A * A^-1:\n";
//         printMatrix(check, n, n);
//         cout << "  Is Identity? " << (check.isIdentity() ? "YES [OK]" : "NO [check error]") << "\n";

//         ostringstream body;
//         body << "\n  Input A:\n" << matrixToStr(M, n, n)
//              << "\n  Inverse A^-1:\n" << matrixToStr(inv, n, n)
//              << "  Time: " << ms << " ms\n";
//         saveToFile("Inverse", body.str());

//     } catch (const exception& e) {
//         cout << "\n  [!] " << e.what() << "\n";
//     }
// }

// // ─── E. MATRIX ADDITION ─────────────────────────────────────────
// void opAddition() {
//     header("MATRIX ADDITION  (A + B)");
//     int n1, n2;
//     Matrix A = inputSquareMatrix(n1, "Matrix A");
//     Matrix B = inputSquareMatrix(n2, "Matrix B");

//     section("Matrix A"); printMatrix(A, n1, n1);
//     section("Matrix B"); printMatrix(B, n2, n2);

//     try {
//         Matrix C = A + B;
//         section("Result:  A + B"); printMatrix(C, n1, n1);

//         ostringstream body;
//         body << "\n  A:\n" << matrixToStr(A,n1,n1)
//              << "\n  B:\n" << matrixToStr(B,n2,n2)
//              << "\n  A+B:\n" << matrixToStr(C,n1,n1);
//         saveToFile("Matrix Addition", body.str());
//     } catch (const exception& e) {
//         cout << "\n  [!] " << e.what() << "\n";
//     }
// }

// // ─── F. MATRIX SUBTRACTION ──────────────────────────────────────
// void opSubtraction() {
//     header("MATRIX SUBTRACTION  (A - B)");
//     int n1, n2;
//     Matrix A = inputSquareMatrix(n1, "Matrix A");
//     Matrix B = inputSquareMatrix(n2, "Matrix B");

//     section("Matrix A"); printMatrix(A, n1, n1);
//     section("Matrix B"); printMatrix(B, n2, n2);

//     try {
//         Matrix C = A - B;
//         section("Result:  A - B"); printMatrix(C, n1, n1);

//         ostringstream body;
//         body << "\n  A:\n" << matrixToStr(A,n1,n1)
//              << "\n  B:\n" << matrixToStr(B,n2,n2)
//              << "\n  A-B:\n" << matrixToStr(C,n1,n1);
//         saveToFile("Matrix Subtraction", body.str());
//     } catch (const exception& e) {
//         cout << "\n  [!] " << e.what() << "\n";
//     }
// }

// // ─── G. MATRIX MULTIPLICATION ───────────────────────────────────
// void opMultiplication() {
//     header("MATRIX MULTIPLICATION  (A x B)");

//     // For multiplication, rows of A can differ from cols of B
//     // So we ask for r,c separately
//     section("Enter Matrix A");
//     cout << "\n  [1] Keyboard  [2] File\n\n";
//     int ch = getInt("  Choice: ", 1, 2);
//     int rA, cA;

//     Matrix A(1,1);
//     if (ch == 1) {
//         rA = getInt("  Rows of A (1-8): ", 1, 8);
//         cA = getInt("  Cols of A (1-8): ", 1, 8);
//         A = Matrix(rA, cA);
//         cout << "\n  Enter A [" << rA << "x" << cA << "]:\n";
//         for (int i = 0; i < rA; i++) {
//             cout << "  Row " << i+1 << ": ";
//             for (int j = 0; j < cA; j++) cin >> A(i,j);
//         }
//     } else {
//         cin.ignore(numeric_limits<streamsize>::max(), '\n');
//         cout << "  File for A [Enter = input/input.txt]: ";
//         string p; getline(cin,p);
//         if (p.empty()) p = "input/input.txt";
//         ifstream fin(p);
//         if (!fin) { cout << "  File not found.\n"; return; }
//         fin >> rA >> cA;
//         A = Matrix(rA, cA);
//         for (int i = 0; i < rA; i++)
//             for (int j = 0; j < cA; j++) fin >> A(i,j);
//         fin.close();
//     }

//     section("Enter Matrix B");
//     cout << "\n  [1] Keyboard  [2] File\n\n";
//     ch = getInt("  Choice: ", 1, 2);
//     int rB, cB;
//     Matrix B(1,1);
//     if (ch == 1) {
//         rB = getInt("  Rows of B (must equal cols of A = " + to_string(cA) + "): ", cA, cA);
//         cB = getInt("  Cols of B (1-8): ", 1, 8);
//         B = Matrix(rB, cB);
//         cout << "\n  Enter B [" << rB << "x" << cB << "]:\n";
//         for (int i = 0; i < rB; i++) {
//             cout << "  Row " << i+1 << ": ";
//             for (int j = 0; j < cB; j++) cin >> B(i,j);
//         }
//     } else {
//         cin.ignore(numeric_limits<streamsize>::max(), '\n');
//         cout << "  File for B: ";
//         string p; getline(cin,p);
//         if (p.empty()) p = "input/input.txt";
//         ifstream fin(p);
//         if (!fin) { cout << "  File not found.\n"; return; }
//         fin >> rB >> cB;
//         B = Matrix(rB, cB);
//         for (int i = 0; i < rB; i++)
//             for (int j = 0; j < cB; j++) fin >> B(i,j);
//         fin.close();
//     }

//     section("Matrix A"); printMatrix(A, rA, cA);
//     section("Matrix B"); printMatrix(B, rB, cB);

//     try {
//         auto t0 = chrono::high_resolution_clock::now();
//         Matrix C = A * B;
//         auto t1 = chrono::high_resolution_clock::now();
//         double ms = chrono::duration<double,milli>(t1-t0).count();

//         section("Result:  A x B  [" + to_string(rA) + "x" + to_string(cB) + "]");
//         printMatrix(C, rA, cB);
//         cout << "  Time: " << fixed << setprecision(4) << ms << " ms\n";

//         ostringstream body;
//         body << "\n  A[" << rA << "x" << cA << "]:\n" << matrixToStr(A,rA,cA)
//              << "\n  B[" << rB << "x" << cB << "]:\n" << matrixToStr(B,rB,cB)
//              << "\n  A*B[" << rA << "x" << cB << "]:\n" << matrixToStr(C,rA,cB)
//              << "  Time: " << ms << " ms\n";
//         saveToFile("Matrix Multiplication", body.str());
//     } catch (const exception& e) {
//         cout << "\n  [!] " << e.what() << "\n";
//     }
// }

// // ─── H. CHECK TRANSPOSE RELATION ────────────────────────────────
// void opIsTranspose() {
//     header("CHECK: Is B the Transpose of A?");
//     int n1, n2;
//     Matrix A = inputSquareMatrix(n1, "Matrix A");
//     Matrix B = inputSquareMatrix(n2, "Matrix B");

//     section("Matrix A"); printMatrix(A, n1, n1);
//     section("Matrix B"); printMatrix(B, n2, n2);

//     bool result = A.isTranspose(B);
//     section("RESULT");
//     cout << "\n  Is B the transpose of A?  "
//          << (result ? "YES" : "NO") << "\n";

//     ostringstream body;
//     body << "\n  A:\n" << matrixToStr(A,n1,n1)
//          << "\n  B:\n" << matrixToStr(B,n2,n2)
//          << "\n  Is B the transpose of A? " << (result ? "YES" : "NO") << "\n";
//     saveToFile("Check Transpose Relation", body.str());
// }

// // ─── I. SOLVE SLE  (single method) ──────────────────────────────
// void opSolveSLE(int methodChoice) {
//     static const string names[] = {
//         "", "Gauss Elimination", "Gauss Jacobi",
//         "LU Crout", "LU Doolittle", "LU Cholesky"
//     };
//     header("SOLVE SLE  --  " + names[methodChoice]);

//     int n;
//     Matrix Aug = inputAugmented(n);
//     int c = n + 1;

//     section("Augmented Matrix [A|b]");
//     cout << "\n";
//     // column labels
//     cout << "  |";
//     for (int j = 0; j < n; j++)
//         cout << setw(10) << right << ("x"+to_string(j+1));
//     cout << "  |" << setw(10) << "b" << "  |\n";
//     cout << "  " << string(10*n+16, '-') << "\n";
//     for (int i = 0; i < n; i++) {
//         cout << "  |";
//         for (int j = 0; j < c; j++) {
//             if (j == c-1) cout << "  |";
//             cout << setw(10) << fixed << setprecision(4) << Aug(i,j);
//         }
//         cout << "  |\n";
//     }
//     cout << "\n";

//     // Build a fresh solver from Aug data
//     auto fill = [&](auto& solver) {
//         for (int i = 0; i < n; i++)
//             for (int j = 0; j < c; j++)
//                 solver(i,j) = Aug(i,j);
//     };

//     vector<double> sol;
//     auto t0 = chrono::high_resolution_clock::now();

//     try {
//         switch (methodChoice) {
//             case 1: { GaussElimination s(n,c); fill(s); sol = s.solve(); break; }
//             case 2: { GaussJacobi     s(n,c); fill(s); sol = s.solve(); break; }
//             case 3: { LUCrout         s(n,c); fill(s); sol = s.solve(); break; }
//             case 4: { LUDoolittle     s(n,c); fill(s); sol = s.solve(); break; }
//             case 5: { LUCholesky      s(n,c); fill(s); sol = s.solve(); break; }
//         }
//     } catch (const exception& e) {
//         cout << "\n  [!] " << e.what() << "\n"; return;
//     }

//     auto t1 = chrono::high_resolution_clock::now();
//     double ms = chrono::duration<double,milli>(t1-t0).count();

//     section("SOLUTION");
//     cout << "\n";
//     for (size_t i = 0; i < sol.size(); i++)
//         cout << "    x" << left << setw(4) << i+1
//              << "=  " << right << setw(16)
//              << fixed << setprecision(8) << sol[i] << "\n";
//     cout << "\n  Time: " << fixed << setprecision(4) << ms << " ms\n";

//     // save
//     ostringstream body;
//     body << "\n  Method: " << names[methodChoice] << "\n";
//     body << "\n  Augmented Matrix:\n";
//     for (int i = 0; i < n; i++) {
//         body << "  |";
//         for (int j = 0; j < c; j++) {
//             if (j==c-1) body << "  |";
//             body << setw(10) << fixed << setprecision(4) << Aug(i,j);
//         }
//         body << "  |\n";
//     }
//     body << "\n  Solution:\n";
//     for (size_t i = 0; i < sol.size(); i++)
//         body << "    x" << i+1 << " = " << fixed << setprecision(8) << sol[i] << "\n";
//     body << "  Time: " << ms << " ms\n";
//     saveToFile("Solve SLE - " + names[methodChoice], body.str());
// }

// // ─── J. SOLVE SLE — ALL METHODS (comparison table) ──────────────
// void opSolveAllMethods() {
//     header("SOLVE SLE  --  ALL 5 METHODS COMPARISON");

//     int n;
//     Matrix Aug = inputAugmented(n);
//     int c = n + 1;

//     section("Augmented Matrix [A|b]");
//     cout << "\n";
//     cout << "  |";
//     for (int j = 0; j < n; j++)
//         cout << setw(10) << right << ("x"+to_string(j+1));
//     cout << "  |" << setw(10) << "b" << "  |\n";
//     cout << "  " << string(10*n+16, '-') << "\n";
//     for (int i = 0; i < n; i++) {
//         cout << "  |";
//         for (int j = 0; j < c; j++) {
//             if (j==c-1) cout << "  |";
//             cout << setw(10) << fixed << setprecision(4) << Aug(i,j);
//         }
//         cout << "  |\n";
//     }

//     const int NM = 5;
//     string names[NM] = {"Gauss Elim","Gauss Jacobi","LU Crout","LU Doolittle","LU Cholesky"};
//     vector<vector<double>> sols(NM);
//     double times[NM];

//     auto fill = [&](auto& solver) {
//         for (int i = 0; i < n; i++)
//             for (int j = 0; j < c; j++)
//                 solver(i,j) = Aug(i,j);
//     };

//     section("Running all methods...");
//     for (int m = 0; m < NM; m++) {
//         cout << "  [" << m+1 << "/5] " << names[m] << "... ";
//         cout.flush();
//         auto t0 = chrono::high_resolution_clock::now();
//         try {
//             switch (m) {
//                 case 0: { GaussElimination s(n,c); fill(s); sols[m]=s.solve(); break; }
//                 case 1: { GaussJacobi     s(n,c); fill(s); sols[m]=s.solve(); break; }
//                 case 2: { LUCrout         s(n,c); fill(s); sols[m]=s.solve(); break; }
//                 case 3: { LUDoolittle     s(n,c); fill(s); sols[m]=s.solve(); break; }
//                 case 4: { LUCholesky      s(n,c); fill(s); sols[m]=s.solve(); break; }
//             }
//         } catch (...) { sols[m] = {}; }
//         auto t1 = chrono::high_resolution_clock::now();
//         times[m] = chrono::duration<double,milli>(t1-t0).count();
//         cout << "Done (" << fixed << setprecision(3) << times[m] << " ms)\n";
//     }

//     section("RESULTS TABLE");
//     cout << "\n  " << left << setw(7) << "Var";
//     for (int m = 0; m < NM; m++)
//         cout << setw(16) << right << names[m];
//     cout << "\n  " << string(7 + 16*NM, '-') << "\n";

//     for (int i = 0; i < n; i++) {
//         cout << "  " << left << setw(7) << ("x"+to_string(i+1));
//         for (int m = 0; m < NM; m++) {
//             if (!sols[m].empty())
//                 cout << setw(16) << right << fixed << setprecision(6) << sols[m][i];
//             else
//                 cout << setw(16) << right << "N/A";
//         }
//         cout << "\n";
//     }
//     cout << "\n  " << left << setw(7) << "ms";
//     for (int m = 0; m < NM; m++)
//         cout << setw(16) << right << fixed << setprecision(4) << times[m];
//     cout << "\n\n";

//     // save
//     ostringstream body;
//     body << "\n  Augmented Matrix:\n";
//     for (int i = 0; i < n; i++) {
//         body << "  |";
//         for (int j = 0; j < c; j++) {
//             if (j==c-1) body << "  |";
//             body << setw(10) << fixed << setprecision(4) << Aug(i,j);
//         }
//         body << "  |\n";
//     }
//     body << "\n  Results:\n  " << left << setw(7) << "Var";
//     for (int m = 0; m < NM; m++) body << setw(16) << right << names[m];
//     body << "\n  " << string(7+16*NM, '-') << "\n";
//     for (int i = 0; i < n; i++) {
//         body << "  " << left << setw(7) << ("x"+to_string(i+1));
//         for (int m = 0; m < NM; m++)
//             if (!sols[m].empty())
//                 body << setw(16) << right << fixed << setprecision(6) << sols[m][i];
//             else body << setw(16) << right << "N/A";
//         body << "\n";
//     }
//     body << "\n  Time (ms): ";
//     for (int m = 0; m < NM; m++) body << names[m] << "=" << fixed << setprecision(4) << times[m] << "  ";
//     body << "\n";
//     saveToFile("SLE - All Methods Comparison", body.str());
// }

// // ================================================================
// //  MAIN MENU
// // ================================================================
// void showMenu() {
//     cout << "\n";
//     line('=');
//     cout << "            MATRIX OPERATIONS MENU\n";
//     line('=');
//     cout << "\n";
//     cout << "  ── MATRIX PROPERTIES ──────────────────────\n";
//     cout << "   [1]  Check Properties   (Square/Symmetric/Diagonal/etc.)\n";
//     cout << "\n";
//     cout << "  ── SINGLE MATRIX OPERATIONS ────────────────\n";
//     cout << "   [2]  Transpose\n";
//     cout << "   [3]  Determinant\n";
//     cout << "   [4]  Inverse\n";
//     cout << "\n";
//     cout << "  ── TWO MATRIX OPERATIONS ───────────────────\n";
//     cout << "   [5]  Addition          \n";
//     cout << "   [6]  Subtraction       \n";
//     cout << "   [7]  Multiplication    \n";
//     cout << "   [8]  Check if Transpose Pair  (is B = A^T?)\n";
//     cout << "\n";
//     cout << "  ── SOLVE LINEAR EQUATIONS (Ax = b) ─────────\n";
//     cout << "   [9]  Gauss Elimination\n";
//     cout << "  [10]  Gauss Jacobi\n";
//     cout << "  [11]  LU Crout\n";
//     cout << "  [12]  LU Doolittle\n";
//     cout << "  [13]  LU Cholesky\n";
//     cout << "  [14]  Run ALL 5 Methods  (comparison table)\n";
//     cout << "\n";
//     line('-', 55);
//     cout << "   [0]  Exit\n";
//     line('=');
// }

// // ================================================================
// //  MAIN
// // ================================================================
// int main() {

//     header("NUMERICAL COMPUTING  --  MATRIX SOLVER");
//     cout << "\n  All matrix operations available at runtime.\n";
//     cout << "  Results can be saved to file after each operation.\n";

//     while (true) {
//         showMenu();
//         int choice = getInt("\n  Enter choice: ", 0, 14);

//         if (choice == 0) break;

//         switch (choice) {
//             case  1: opProperties();         break;
//             case  2: opTranspose();           break;
//             case  3: opDeterminant();         break;
//             case  4: opInverse();             break;
//             case  5: opAddition();            break;
//             case  6: opSubtraction();         break;
//             case  7: opMultiplication();      break;
//             case  8: opIsTranspose();         break;
//             case  9: opSolveSLE(1);           break;  // Gauss Elimination
//             case 10: opSolveSLE(2);           break;  // Gauss Jacobi
//             case 11: opSolveSLE(3);           break;  // LU Crout
//             case 12: opSolveSLE(4);           break;  // LU Doolittle
//             case 13: opSolveSLE(5);           break;  // LU Cholesky
//             case 14: opSolveAllMethods();     break;
//             default: cout << "  Invalid choice.\n";
//         }

//         cout << "\n";
//         line('-', 55);
//         cout << "  Operation complete. Returning to main menu...\n";
//     }

//     line('=');
//     cout << "  Thank you for using NC Matrix Solver!\n";
//     line('=');
//     cout << "\n";
//     return 0;
// }





















































// // #include <iostream>
// // #include <fstream>
// // #include "include/Matrix.hpp"

// // using namespace std;

// // int main() {

// //     ifstream fin("input/input.txt");
// //     ofstream fout("output/output.txt");

// //     if(!fin) {
// //         cout << "Input file not found\n";
// //         return 0;
// //     }

// //     int r,c;

// //     fin >> r >> c;

// //     Matrix A(r,c);

// //     fin >> A;

// //     fout << "Input Matrix:\n";
// //     fout << A;

// //     Matrix T = A.transpose();

// //     fout << "\nTranspose Matrix:\n";
// //     fout << T;

// //     if(A.isSquare()) {
// //         fout << "\nDeterminant = " << A.determinant() << endl;
// //     }

// //     fin.close();
// //     fout.close();

// //     cout << "Program executed successfully.\n";
// //     cout << "Check output.txt for results.\n";

// //     return 0;
// // }