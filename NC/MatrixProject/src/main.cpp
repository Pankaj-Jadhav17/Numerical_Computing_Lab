#include "../include/gaussian.hpp"
#include <iostream>
#include <string>
#include <iomanip>
#include <chrono>
using namespace std;

extern Matrix readFromFile(const string &);
extern void writeToFile(const string &, vector<double> &);

int main()
{
    try
    {
        //Choose input file 
        cout << "\nGaussian Elimination Solver \n\n";
        cout << "Select system to solve:\n";
        cout << "  1. system_3x3.txt\n";
        cout << "  2. system_4x4.txt\n";
        cout << "  3. system_5x5.txt\n";
        cout << "  4. system_6x6.txt\n";
        cout << "  5. system_20x20.txt\n";
        cout << "  6. system_100x100.txt\n";
        cout << "  7. system_200x200.txt\n";
        cout << "  8. system_500x500.txt\n";
        cout << "Enter choice (1-8): ";

        int fileChoice;
        cin >> fileChoice;

        string inputFile;
        switch (fileChoice)
        {
            case 1: inputFile = "input/system_3x3.txt";     break;
            case 2: inputFile = "input/system_4x4.txt";     break;
            case 3: inputFile = "input/system_5x5.txt";     break;
            case 4: inputFile = "input/system_6x6.txt";     break;
            case 5: inputFile = "input/system_20x20.txt";   break;
            case 6: inputFile = "input/system_100x100.txt"; break;
            case 7: inputFile = "input/system_200x200.txt"; break;
            case 8: inputFile = "input/system_500x500.txt"; break;
            default:
                cout << "Invalid choice. Using system_4x4.txt\n";
                inputFile = "input/system_4x4.txt";
        }

        //Choose method
        cout << "\nSelect method:\n";
        cout << "  1. Without Pivoting\n";
        cout << "  2. Partial Pivoting (Recommended)\n";
        cout << "Enter choice (1 or 2): ";

        int methodChoice;
        cin >> methodChoice;

        bool usePivot = (methodChoice == 2);
        cout << "\nUsing  : " << (usePivot ? "Partial Pivoting" : "No Pivoting") << "\n";
        cout << "Reading: " << inputFile << "\n";

        //Read 
        Matrix A = readFromFile(inputFile);
        int n = A.getRows();
        cout << "Size   : " << n << " x " << n << "\n\n";

        //Solve (with timing) 
        Solver *solver = new GaussianElimination(usePivot);

        auto start = chrono::high_resolution_clock::now();
        vector<double> solution = solver->solve(A);
        auto end   = chrono::high_resolution_clock::now();

        delete solver;

        double ms = chrono::duration<double, milli>(end - start).count();
        cout << "Time   : " << fixed << setprecision(3) << ms << " ms\n\n";

        //Print solution 
        cout << "Solution\n";
        cout << fixed << setprecision(6);

        // For large systems only print first 10 and last 5 values
        if (n <= 20)
        {
            for (int i = 0; i < n; i++)
                cout << "  x" << i + 1 << " = " << solution[i] << "\n";
        }
        else
        {
            cout << "  (Showing first 10 and last 5 — full list in output file)\n";
            for (int i = 0; i < 10; i++)
                cout << "  x" << i + 1 << " = " << solution[i] << "\n";
            cout << "  ...\n";
            for (int i = n - 5; i < n; i++)
                cout << "  x" << i + 1 << " = " << solution[i] << "\n";
        }
        //Save to file 
        writeToFile("output/result_all.txt", solution);
        cout << "\nFull solution saved to output/result_all.txt\n";
    }
    catch (const exception &e)
    {
        cout << "\n[ERROR] " << e.what() << "\n";
        return 1;
    }

    return 0;
}








// ofstream perf("output/performance.txt", ios::app);
// perf << n << " " << ms << endl;
// perf.close();








// #include "../include/gaussian.hpp"
// #include <iostream>
// #include <string>
// #include <iomanip>

// extern Matrix readFromFile(const string &);
// extern void writeToFile(const string &, vector<double> &);

// int main()
// {
//     try
//     {
//         //Choose input file
//         cout << "\n Gaussian Elimination Solver\n\n";
//         cout << "Select system to solve:\n";
//         cout << "  1. system_3x3.txt\n";
//         cout << "  2. system_4x4.txt\n";
//         cout << "  3. system_5x5.txt\n";
//         cout << "  4. system_6x6.txt\n";
//         cout << "  5. system_20x20.txt\n";
//         cout << "Enter choice (1-5): ";

//         int fileChoice;
//         cin >> fileChoice;

//         string inputFile;
//         switch (fileChoice)
//         {
//             case 1: inputFile = "input/system_3x3.txt";   break;
//             case 2: inputFile = "input/system_4x4.txt";   break;
//             case 3: inputFile = "input/system_5x5.txt";   break;
//             case 4: inputFile = "input/system_6x6.txt";   break;
//             case 5: inputFile = "input/system_20x20.txt"; break;
//             default:
//                 cout << "Invalid choice. Using system_4x4.txt\n";
//                 inputFile = "input/system_4x4.txt";
//         }

//         //Choose method
//         cout << "\nSelect method:\n";
//         cout << "  1. Without Pivoting\n";
//         cout << "  2. Partial Pivoting (Recommended)\n";
//         cout << "Enter choice (1 or 2): ";

//         int methodChoice;
//         cin >> methodChoice;

//         bool usePivot = (methodChoice == 2);
//         cout << "\nUsing: " << (usePivot ? "Partial Pivoting" : "No Pivoting") << "\n";
//         cout << "Reading: " << inputFile << "\n";

//         // Read, Solve, Write 
//         Matrix A = readFromFile(inputFile);

//         int n = A.getRows();
//         cout << "Matrix size: " << n << "x" << n << "\n\n";

//         Solver *solver = new GaussianElimination(usePivot);
//         vector<double> solution = solver->solve(A);
//         delete solver;

//         // Print to screen
//         cout << " Solution \n";
//         cout << fixed << setprecision(6);
//         for (int i = 0; i < (int)solution.size(); i++)
//             cout << "  x" << i + 1 << " = " << solution[i] << "\n";

//         // Write to file
//         writeToFile("output/result_all.txt", solution);
//         cout << "\nSaved to output/result_all.txt\n";
//     }
//     catch (const exception &e)
//     {
//         cout << "\n[ERROR] " << e.what() << "\n";
//         return 1;
//     }

//     return 0;
// }










// #include "../include/gaussian.hpp"
// #include <iostream>

// extern Matrix readFromFile(string);
// extern void writeToFile(string, vector<double>&);

// int main()
// {
//     try
//     {
//         int choice;
//         cout << "1. Without Pivoting\n2. Partial Pivoting\n";
//         cin >> choice;

//         Matrix A = readFromFile("../input/system_4x4.txt");

//         Solver* solver;

//         if (choice == 1)
//             solver = new GaussianElimination(false);
//         else
//             solver = new GaussianElimination(true);

//         vector<double> solution = solver->solve(A);

//         writeToFile("../output/result_all.txt", solution);

//         delete solver;

//         cout << "Solved! Check output file.\n";
//     }
//     catch (exception &e)
//     {
//         cout << "Error: " << e.what();
//     }
// }









// #include "include/matrix.hpp"
// int main()
// {
//     try
//     {
//         int choice;
//         cout << "1. Addition\n2. Subtraction\n3. Gaussian Elimination\n";
//         cin >> choice;

//         if (choice == 1 || choice == 2)
//         {
//             int r, c;
//             cin >> r >> c;

//             MatrixBase *A = new Matrix(r, c);
//             MatrixBase *B = new Matrix(r, c);

//             A->read();
//             B->read();

//             MatrixBase *C;

//             if (choice == 1)
//                 C = A->add(*B);
//             else
//                 C = A->subtract(*B);

//             C->display();

//             delete A;
//             delete B;
//             delete C;
//         }
//         else if (choice == 3)
//         {
//             int n;
//             cin >> n;

//             MatrixBase *aug = new Matrix(n, n + 1);
//             aug->read();

//             aug->gaussianElimination();
//             aug->backSubstitution();

//             delete aug;
//         }
//     }
//     catch (exception &e)
//     {
//         cout << "Error: " << e.what();
//     }

//     return 0;
// }