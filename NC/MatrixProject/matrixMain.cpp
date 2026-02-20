#include "include/matrix.hpp"
int main()
{
    try
    {
        int choice;
        cout << "1. Addition\n2. Subtraction\n3. Gaussian Elimination\n";
        cin >> choice;

        if (choice == 1 || choice == 2)
        {
            int r, c;
            cin >> r >> c;

            MatrixBase *A = new Matrix(r, c);
            MatrixBase *B = new Matrix(r, c);

            A->read();
            B->read();

            MatrixBase *C;

            if (choice == 1)
                C = A->add(*B);
            else
                C = A->subtract(*B);

            C->display();

            delete A;
            delete B;
            delete C;
        }
        else if (choice == 3)
        {
            int n;
            cin >> n;

            MatrixBase *aug = new Matrix(n, n + 1);
            aug->read();

            aug->gaussianElimination();
            aug->backSubstitution();

            delete aug;
        }
    }
    catch (exception &e)
    {
        cout << "Error: " << e.what();
    }

    return 0;
}