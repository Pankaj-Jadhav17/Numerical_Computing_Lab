#include "../include/matrix.hpp"

// Parent constructors
MatrixBase::MatrixBase() : rows(0), cols(0) {}
MatrixBase::MatrixBase(int r, int c) : rows(r), cols(c) {}
MatrixBase::~MatrixBase() {}

//Child constructors
Matrix::Matrix() : MatrixBase() {}

Matrix::Matrix(int r, int c) : MatrixBase(r, c)
{
    if (r <= 0 || c <= 0)
        throw invalid_argument("Invalid matrix size");

    this->mat.resize(rows, vector<double>(cols, 0));
}

// Copy constructor
Matrix::Matrix(const Matrix &other)
{
    this->rows = other.rows;
    this->cols = other.cols;
    this->mat = other.mat;
}

// Destructor
Matrix::~Matrix()
{
    cout << "Matrix destroyed\n";
}

// Read
void Matrix::read()
{
    cout << "Enter elements:\n";
    for (int i = 0; i < this->rows; i++)
        for (int j = 0; j < this->cols; j++)
            cin >> this->mat[i][j];
}

// Display
void Matrix::display() const
{
    for (auto &row : mat)
    {
        for (auto val : row)
            cout << val << " ";
        cout << endl;
    }
}

// Addition
MatrixBase* Matrix::add(const MatrixBase &m) const
{
    const Matrix &other = dynamic_cast<const Matrix&>(m);

    if (this->rows != other.rows || this->cols != other.cols)
        throw logic_error("Addition not possible");

    Matrix *result = new Matrix(this->rows, this->cols);

    for (int i = 0; i < this->rows; i++)
        for (int j = 0; j < this->cols; j++)
            result->mat[i][j] = this->mat[i][j] + other.mat[i][j];

    return result;
}

// Subtraction
MatrixBase* Matrix::subtract(const MatrixBase &m) const
{
    const Matrix &other = dynamic_cast<const Matrix&>(m);

    Matrix *result = new Matrix(this->rows, this->cols);

    for (int i = 0; i < this->rows; i++)
        for (int j = 0; j < this->cols; j++)
            result->mat[i][j] = this->mat[i][j] - other.mat[i][j];

    return result;
}

// Gaussian Elimination
void Matrix::gaussianElimination()
{
    for (int k = 0; k < this->rows; k++)
    {
        if (mat[k][k] == 0)
            throw runtime_error("Zero pivot");

        for (int i = k + 1; i < this->rows; i++)
        {
            double factor = mat[i][k] / mat[k][k];

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
    if (cols != rows + 1)
        throw logic_error("Need augmented matrix");

    vector<double> x(rows);

    for (int i = rows - 1; i >= 0; i--)
    {
        x[i] = mat[i][cols - 1];

        for (int j = i + 1; j < rows; j++)
            x[i] -= mat[i][j] * x[j];

        x[i] /= mat[i][i];
    }

    cout << "\nSolutions:\n";
    for (int i = 0; i < rows; i++)
        cout << "x" << i + 1 << " = " << x[i] << endl;
}




























































// #include "../include/matrix.hpp"

// // Default constructor
// Matrix::Matrix()
// {
//     this->rows = 0;
//     this->cols = 0;
// }

// // Parameterized constructor
// Matrix::Matrix(int r, int c)
// {
//     if (r <= 0 || c <= 0)
//         throw invalid_argument("Invalid matrix size");

//     this->rows = r;
//     this->cols = c;
//     mat.resize(rows, vector<double>(cols, 0));
// }

// // Virtual constructor
// MatrixBase* Matrix::create(int r, int c)
// {
//     return new Matrix(r, c);
// }

// // Read matrix
// void Matrix::read()
// {
//     cout << "Enter elements:\n";
//     for (int i = 0; i < this->rows; i++)
//         for (int j = 0; j < this->cols; j++)
//             cin >> mat[i][j];
// }

// // Display matrix
// void Matrix::display()
// {
//     for (auto &row : mat)
//     {
//         for (auto val : row)
//             cout << val << " ";
//         cout << endl;
//     }
// }

// // Addition (method overriding)
// MatrixBase* Matrix::add(MatrixBase &m)
// {
//     Matrix &other = dynamic_cast<Matrix&>(m);

//     if (this->rows != other.rows || this->cols != other.cols)
//         throw logic_error("Addition not possible");

//     Matrix *result = new Matrix(this->rows, this->cols);

//     for (int i = 0; i < this->rows; i++)
//         for (int j = 0; j < this->cols; j++)
//             result->mat[i][j] = this->mat[i][j] + other.mat[i][j];

//     return result;
// }

// // Subtraction
// MatrixBase* Matrix::subtract(MatrixBase &m)
// {
//     Matrix &other = dynamic_cast<Matrix&>(m);

//     Matrix *result = new Matrix(this->rows, this->cols);

//     for (int i = 0; i < this->rows; i++)
//         for (int j = 0; j < this->cols; j++)
//             result->mat[i][j] = this->mat[i][j] - other.mat[i][j];

//     return result;
// }

// // Gaussian elimination
// void Matrix::gaussianElimination()
// {
//     for (int k = 0; k < this->rows; k++)
//     {
//         if (mat[k][k] == 0)
//             throw runtime_error("Zero pivot");

//         double pivot = mat[k][k];

//         for (int j = k; j < this->cols; j++)
//             mat[k][j] /= pivot;

//         for (int i = k + 1; i < this->rows; i++)
//         {
//             double factor = mat[i][k];

//             for (int j = k; j < this->cols; j++)
//                 mat[i][j] -= factor * mat[k][j];
//         }
//     }

//     cout << "\nRow Echelon Form:\n";
//     display();
// }

// // Back substitution
// void Matrix::backSubstitution()
// {
//     if (this->cols != this->rows + 1)
//         throw logic_error("Need augmented matrix");

//     vector<double> x(this->rows);

//     for (int i = this->rows - 1; i >= 0; i--)
//     {
//         x[i] = mat[i][this->cols - 1];

//         for (int j = i + 1; j < this->rows; j++)
//             x[i] -= mat[i][j] * x[j];
//     }

//     cout << "\nSolutions:\n";
//     for (int i = 0; i < this->rows; i++)
//         cout << "x" << i + 1 << " = " << x[i] << endl;
// }
