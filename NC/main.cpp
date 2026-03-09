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

    int r,c;
    fin >> r >> c;

    cout << "\nChoose Method:\n";
    cout << "1. Gauss Elimination\n";
    cout << "2. Gauss Jacobi\n";
    cout << "3. LU Crout\n";
    cout << "4. LU Doolittle\n";
    cout << "5. LU Cholesky\n";
    cout << "Enter choice: ";

    int choice;
    cin >> choice;

    vector<double> solution;

    switch(choice)
    {
        case 1:
        {
            GaussElimination solver(r,c);
            fin >> solver;

            fout << "Using Gauss Elimination\n";
            fout << "Input Matrix:\n" << solver;

            solution = solver.solve();
            break;
        }

        case 2:
        {
            GaussJacobi solver(r,c);
            fin >> solver;

            fout << "Using Gauss Jacobi\n";
            fout << "Input Matrix:\n" << solver;

            solution = solver.solve();
            break;
        }

        case 3:
        {
            LUCrout solver(r,c);
            fin >> solver;

            fout << "Using LU Crout\n";
            fout << "Input Matrix:\n" << solver;

            solution = solver.solve();
            break;
        }

        case 4:
        {
            LUDoolittle solver(r,c);
            fin >> solver;

            fout << "Using LU Doolittle\n";
            fout << "Input Matrix:\n" << solver;

            solution = solver.solve();
            break;
        }

        case 5:
        {
            LUCholesky solver(r,c);
            fin >> solver;

            fout << "Using LU Cholesky\n";
            fout << "Input Matrix:\n" << solver;

            solution = solver.solve();
            break;
        }

        default:
        {
            cout << "Invalid choice\n";
            return 0;
        }
    }

    fout << "\nSolution:\n";

    for(size_t i = 0; i < solution.size(); i++)
    {
        fout << "x" << i+1 << " = " << solution[i] << endl;
    }

    fin.close();
    fout.close();

    cout << "\nProgram executed successfully.\n";
    cout << "Check output/output.txt\n";

    return 0;
}




















// #include <iostream>
// #include <fstream>
// #include "Matrix.hpp"
// #include "GaussElimination.hpp"
// #include "GaussJacobi.hpp"

// using namespace std;

// int main()
// {
//     ifstream fin("input/input.txt");
//     ofstream fout("output/output.txt");

//     int r,c;

//     fin >> r >> c;

//     Matrix A(r,c);

//     fin >> A;

//     fout << "Input Matrix:\n";
//     fout << A;

//     if(A.isDiagonallyDominant())
//     {
//         fout << "\nUsing Jacobi Method\n";

//         GaussJacobi solver(A);
//         solver.solve(fout);
//     }
//     else
//     {
//         fout << "\nUsing Gauss Elimination\n";

//         GaussElimination solver(A);
//         solver.solve(fout);
//     }

//     fin.close();
//     fout.close();

//     return 0;
// }












// #include <iostream>
// #include <fstream>
// #include "include/Matrix.hpp"

// using namespace std;

// int main() {

//     ifstream fin("input/input.txt");
//     ofstream fout("output/output.txt");

//     if(!fin) {
//         cout << "Input file not found\n";
//         return 0;
//     }

//     int r,c;

//     fin >> r >> c;

//     Matrix A(r,c);

//     fin >> A;

//     fout << "Input Matrix:\n";
//     fout << A;

//     Matrix T = A.transpose();

//     fout << "\nTranspose Matrix:\n";
//     fout << T;

//     if(A.isSquare()) {
//         fout << "\nDeterminant = " << A.determinant() << endl;
//     }

//     fin.close();
//     fout.close();

//     cout << "Program executed successfully.\n";
//     cout << "Check output.txt for results.\n";

//     return 0;
// }