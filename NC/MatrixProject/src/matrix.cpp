#include "../include/matrix.hpp"

// Default constructor
Matrix::Matrix()
{
    this->rows = 0;
    this->cols = 0;
}

// Parameterized constructor
Matrix::Matrix(int r, int c)
{
    if (r <= 0 || c <= 0)
        throw invalid_argument("Invalid matrix size");

    this->rows = r;
    this->cols = c;
    mat.resize(rows, vector<double>(cols, 0));
}

// Virtual constructor
MatrixBase* Matrix::create(int r, int c)
{
    return new Matrix(r, c);
}

// Read matrix
void Matrix::read()
{
    cout << "Enter elements:\n";
    for (int i = 0; i < this->rows; i++)
        for (int j = 0; j < this->cols; j++)
            cin >> mat[i][j];
}

// Display matrix
void Matrix::display()
{
    for (auto &row : mat)
    {
        for (auto val : row)
            cout << val << " ";
        cout << endl;
    }
}

// Addition (method overriding)
MatrixBase* Matrix::add(MatrixBase &m)
{
    Matrix &other = dynamic_cast<Matrix&>(m);

    if (this->rows != other.rows || this->cols != other.cols)
        throw logic_error("Addition not possible");

    Matrix *result = new Matrix(this->rows, this->cols);

    for (int i = 0; i < this->rows; i++)
        for (int j = 0; j < this->cols; j++)
            result->mat[i][j] = this->mat[i][j] + other.mat[i][j];

    return result;
}

// Subtraction
MatrixBase* Matrix::subtract(MatrixBase &m)
{
    Matrix &other = dynamic_cast<Matrix&>(m);

    Matrix *result = new Matrix(this->rows, this->cols);

    for (int i = 0; i < this->rows; i++)
        for (int j = 0; j < this->cols; j++)
            result->mat[i][j] = this->mat[i][j] - other.mat[i][j];

    return result;
}

// Gaussian elimination
void Matrix::gaussianElimination()
{
    for (int k = 0; k < this->rows; k++)
    {
        if (mat[k][k] == 0)
            throw runtime_error("Zero pivot");

        double pivot = mat[k][k];

        for (int j = k; j < this->cols; j++)
            mat[k][j] /= pivot;

        for (int i = k + 1; i < this->rows; i++)
        {
            double factor = mat[i][k];

            for (int j = k; j < this->cols; j++)
                mat[i][j] -= factor * mat[k][j];
        }
    }

    cout << "\nRow Echelon Form:\n";
    display();
}

// Back substitution
void Matrix::backSubstitution()
{
    if (this->cols != this->rows + 1)
        throw logic_error("Need augmented matrix");

    vector<double> x(this->rows);

    for (int i = this->rows - 1; i >= 0; i--)
    {
        x[i] = mat[i][this->cols - 1];

        for (int j = i + 1; j < this->rows; j++)
            x[i] -= mat[i][j] * x[j];
    }

    cout << "\nSolutions:\n";
    for (int i = 0; i < this->rows; i++)
        cout << "x" << i + 1 << " = " << x[i] << endl;
}
