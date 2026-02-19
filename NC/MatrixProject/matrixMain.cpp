#include "include/matrix.hpp"

int main()
{
    try
    {
        MatrixBase *A = new Matrix(2, 2);
        MatrixBase *B = new Matrix(2, 2);

        cout << "Enter Matrix A:\n";
        A->read();

        cout << "Enter Matrix B:\n";
        B->read();

        MatrixBase *C = A->add(*B);
        cout << "\nA + B:\n";
        C->display();

        // Solve equations
        int n;
        cout << "\nEnter number of variables: ";
        cin >> n;

        MatrixBase *aug = new Matrix(n, n + 1);
        aug->read();

        aug->gaussianElimination();
        aug->backSubstitution();

        delete A;
        delete B;
        delete C;
        delete aug;
    }
    catch (exception &e)
    {
        cout << "Error: " << e.what();
    }

    return 0;
}
