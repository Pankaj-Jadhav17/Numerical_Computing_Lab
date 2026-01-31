// f(x)=x^2−2
// Initial interval [1,2]
// Stopping condition
// |mcurrent - mprevious|<10^-k
// Also computes theoretical number of iterations using
// n≥log⁡(b−a)−log⁡(tolerance)/log⁡2

#include <iostream>
#include <cmath>
using namespace std;

// Step 1: Define the function p(x) = x^2 - 2
double p(double x) {
    return x * x - 2;
}

int main() {
    double a = 1.0, b = 2.0;   // Initial interval [1, 2]
    int k;                     // tolerance power
    double tolerance;

    cout << "Enter k (tolerance = 10^-k): ";
    cin >> k;

    tolerance = pow(10, -k);

    // Step 2: Theoretical number of iterations
    int n_theory = ceil((log(b - a) - log(tolerance)) / log(2));

    cout << "\nTheoretical minimum iterations = " << n_theory << endl;

    // Step 3: Bisection iterations
    double m_prev = a;
    double m_curr;
    int iterations = 0;

    cout << "\nIteration details:\n";
    cout << "Iter\tMidpoint\tError\n";

    while (true) {
        m_curr = (a + b) / 2.0;
        iterations++;

        double error = fabs(m_curr - m_prev);

        cout << iterations << "\t"
             << m_curr << "\t"
             << error << endl;

        // Step 4: Stopping condition
        if (error < tolerance)
            break;

        // Update interval
        if (p(a) * p(m_curr) < 0)
            b = m_curr;
        else
            a = m_curr;

        m_prev = m_curr;
    }

    cout << "\nApproximate root = " << m_curr << endl;
    cout << "Actual sqrt(2)   = " << sqrt(2) << endl;
    cout << "Total iterations = " << iterations << endl;

    return 0;
}



//What This Program Does (Concept Mapping)
// | Concept                  | Code Implementation                 |
// | ------------------------ | ----------------------------------- |
// | ( p(x) = x^2 - 2 )       | `double p(double x)`                |
// | Initial interval ([1,2]) | `double a = 1.0, b = 2.0;`          |
// | Midpoint formula         | `(a + b) / 2.0`                     |
// | Sign check               | `p(a) * p(m_curr) < 0`              |
// | Stopping condition       | `fabs(m_curr - m_prev) < tolerance` |
// | Tolerance (10^{-k})      | `pow(10, -k)`                       |
// | Iteration formula        | `ceil((log(b-a)-log(tol))/log(2))`  |

// #include <iostream>
// #include <cmath>
// using namespace std;

// double p(double x) {
//     return x * x - 2;
// }

// int main() {
//     double a = 1.0, b = 2.0;
//     int k;
//     cout << "Enter k (tolerance = 10^-k): ";
//     cin >> k;

//     double tolerance = pow(10, -k);

//     int n_theory = ceil((log(b - a) - log(tolerance)) / log(2));
//     cout << "\nTheoretical minimum iterations = " << n_theory << endl;

//     cout << "\nIter\tMidpoint\tError\n";

//     for (int i = 1; i <= n_theory; i++) {
//         double m = (a + b) / 2.0;
//         double error = (b - a) / 2.0;

//         cout << i << "\t" << m << "\t" << error << endl;

//         if (p(a) * p(m) < 0)
//             b = m;
//         else
//             a = m;
//     }

//     cout << "\nApproximate root = " << (a + b) / 2.0 << endl;
//     cout << "Actual sqrt(2)   = " << sqrt(2) << endl;

//     return 0;
// }
