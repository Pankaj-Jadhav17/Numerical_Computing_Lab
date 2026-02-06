#include "BM.hpp"
#include <iostream>
#include <iomanip>

using namespace std;
void header() {
    cout << "\nBisection Method for f(x) = 4x^3 - 3x\n";
    cout << "Roots: x = 0, ±0.866025\n\n";
}
void chooseInterval(double& a, double& b) {
    int c;
    cout << "1) [-0.5, 0.5]\n"
         << "2) [0.5, 1.0]\n"
         << "3) [-1.0, -0.5]\n"
         << "4) Custom\n"
         << "Choice: ";
    cin >> c;

    if (c == 1) a = -0.5, b = 0.5;
    else if (c == 2) a = 0.5, b = 1.0;
    else if (c == 3) a = -1.0, b = -0.5;
    else {
        cout << "Enter a b: ";
        cin >> a >> b;
    }
}

int main() {
    header();

    double a, b;
    int k;

    chooseInterval(a, b);

    cout << "Enter k for tolerance (10^-k): ";
    cin >> k;

    BisectionMethod bm(a, b, k);

    if (!bm.validateInterval()) {
        cout << "Invalid interval\n";
        return 1;
    }

    double root = bm.findRoot();

    char ch;
    cout << "Show iteration table? (y/n): ";
    cin >> ch;

    if (ch == 'y' || ch == 'Y')
        bm.displayIterationTable();

    cout << fixed << setprecision(10);
    cout << "\nRoot found: " << root << endl;
    cout << "Iterations: " << bm.getActualIterations() << endl;
    cout << "Tolerance:  " << bm.getTolerance() << endl;
    return 0;
}





























// #include <iostream>
// #include <cmath>
// #include "BM.hpp"


// using namespace std;

// // Example functions to solve
// // 1. f(x) = x^2 - 2  (find sqrt(2))
// double function1(double x) {
//     return x * x - 2;
// }

// // 2. f(x) = x^3 - x - 2  (cubic equation)
// double function2(double x) {
//     return x * x * x - x - 2;
// }

// // 3. f(x) = cos(x) - x  (transcendental equation)
// double function3(double x) {
//     return cos(x) - x;
// }

// // 4. f(x) = e^x - 3x  (exponential equation)
// double function4(double x) {
//     return exp(x) - 3 * x;
// }

// // 5. f(x) = x^3 - 6x^2 + 11x - 6  (polynomial with known roots: 1, 2, 3)
// double function5(double x) {
//     return x * x * x - 6 * x * x + 11 * x - 6;
// }

// // Display menu
// void displayMenu() {
//     cout << "\n" << string(50, '=') << endl;
//     cout << "BISECTION METHOD - EQUATION SOLVER" << endl;
//     cout << string(50, '=') << endl;
//     cout << "Choose an equation to solve:" << endl;
//     cout << "1. f(x) = 4*x^3 - 3*x" << endl;
//     // cout << "1. f(x) = x^2 - 2          (find sqrt(2))" << endl;
//     // cout << "2. f(x) = x^3 - x - 2      (cubic)" << endl;
//     // cout << "3. f(x) = cos(x) - x       (transcendental)" << endl;
//     // cout << "4. f(x) = e^x - 3x         (exponential)" << endl;
//     // cout << "5. f(x) = x^3 - 6x^2 + 11x - 6  (polynomial)" << endl;
//     cout << "6. Custom function (modify code)" << endl;
//     cout << "0. Exit" << endl;
//     cout << string(50, '-') << endl;
// }

// // Get function based on choice
// function<double(double)> getFunction(int choice) {
//     switch (choice) {
//         case 1: return function1;
//         case 2: return function2;
//         case 3: return function3;
//         case 4: return function4;
//         case 5: return function5;
//         default: return function1;
//     }
// }

// // Get suggested intervals for each function
// void getSuggestedInterval(int choice, double& a, double& b) {
//     switch (choice) {
//         case 1: a = 1.0; b = 2.0; break;      // sqrt(2) is between 1 and 2
//         case 2: a = 1.0; b = 2.0; break;      // root is between 1 and 2
//         case 3: a = 0.0; b = 1.0; break;      // root is between 0 and 1
//         case 4: a = 0.0; b = 1.0; break;      // one root between 0 and 1
//         case 5: a = 0.5; b = 1.5; break;      // root at x=1
//         default: a = 1.0; b = 2.0; break;
//     }
// }

// int main() {
//     while (true) {
//         displayMenu();
        
//         int choice;
//         cout << "Enter your choice: ";
//         cin >> choice;
        
//         if (choice == 0) {
//             cout << "\nThank you for using Bisection Method Solver!" << endl;
//             break;
//         }
        
//         if (choice < 1 || choice > 6) {
//             cout << "Invalid choice! Please try again." << endl;
//             continue;
//         }
        
//         // Get the function
//         auto func = getFunction(choice);
        
//         // Get interval
//         double a, b;
//         char useDefault;
        
//         getSuggestedInterval(choice, a, b);
//         cout << "\nSuggested interval: [" << a << ", " << b << "]" << endl;
//         cout << "Use default interval? (y/n): ";
//         cin >> useDefault;
        
//         if (useDefault == 'n' || useDefault == 'N') {
//             cout << "Enter left endpoint (a): ";
//             cin >> a;
//             cout << "Enter right endpoint (b): ";
//             cin >> b;
//         }
        
//         // Get tolerance
//         int k;
//         cout << "Enter k for tolerance (10^-k): ";
//         cin >> k;
//         double tolerance = pow(10, -k);
        
//         // Create BisectionMethod object
//         BisectionMethod bm(func, a, b, tolerance);
        
//         // Display theoretical iterations
//         cout << "\nTheoretical minimum iterations: " 
//              << bm.getTheoreticalIterations() << endl;
        
//         // Find root
//         cout << "\nSolving..." << endl;
//         double root = bm.findRoot();
        
//         if (!isnan(root)) {
//             // Display iteration table
//             char showTable;
//             cout << "\nShow iteration table? (y/n): ";
//             cin >> showTable;
            
//             if (showTable == 'y' || showTable == 'Y') {
//                 bm.displayIterationTable();
//             }
            
//             // Display results
//             bm.displayResults(root, bm.getHistory().size());
            
//             // Special case: show actual sqrt(2) for function 1
//             if (choice == 1) {
//                 cout << "\nActual sqrt(2) = " << sqrt(2) << endl;
//                 cout << "Absolute error = " << abs(root - sqrt(2)) << endl;
//             }
//         }
        
//         cout << "\nPress Enter to continue...";
//         cin.ignore();
//         cin.get();
//     }
    
//     return 0;
// }