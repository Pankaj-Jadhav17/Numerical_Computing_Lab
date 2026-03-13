#include "../include/Matrix.hpp"
#include <cmath>
#include <stdexcept>

using namespace std;

//  CONSTRUCTORS
Matrix::Matrix() {
    rows = 0;
    cols = 0;
}

Matrix::Matrix(int r, int c) {
    rows = r;
    cols = c;
    data.resize(r, vector<double>(c, 0.0));
}

Matrix::Matrix(const Matrix& m) {
    rows = m.rows;
    cols = m.cols;
    data = m.data;
}

//  FILE I/O
void Matrix::readFromFile(ifstream& fin) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            fin >> data[i][j];
}

void Matrix::displayToFile(ofstream& fout) const {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++)
            fout << data[i][j] << " ";
        fout << endl;
    }
}

//  ELEMENT ACCESS
double& Matrix::operator()(int i, int j) {
    return data[i][j];
}

double Matrix::operator()(int i, int j) const {
    return data[i][j];
}

//  STREAM OPERATORS
istream& operator>>(istream& in, Matrix& m) {
    for (int i = 0; i < m.rows; i++)
        for (int j = 0; j < m.cols; j++)
            in >> m.data[i][j];
    return in;
}

ostream& operator<<(ostream& out, const Matrix& m) {
    for (int i = 0; i < m.rows; i++) {
        for (int j = 0; j < m.cols; j++) {
            out << m.data[i][j];
            if (j < m.cols - 1) out << " ";
        }
        out << endl;
    }
    return out;
}

//  ARITHMETIC OPERATORS
Matrix Matrix::operator+(const Matrix& m) const {
    if (rows != m.rows || cols != m.cols)
        throw invalid_argument("Matrix dimensions must match for addition.");
    Matrix result(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.data[i][j] = data[i][j] + m.data[i][j];
    return result;
}

Matrix Matrix::operator-(const Matrix& m) const {
    if (rows != m.rows || cols != m.cols)
        throw invalid_argument("Matrix dimensions must match for subtraction.");
    Matrix result(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.data[i][j] = data[i][j] - m.data[i][j];
    return result;
}

// A(m x n) * B(n x p) = C(m x p)
Matrix Matrix::operator*(const Matrix& m) const {
    if (cols != m.rows)
        throw invalid_argument("Matrix dimensions incompatible for multiplication (A cols must equal B rows).");
    Matrix result(rows, m.cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < m.cols; j++)
            for (int k = 0; k < cols; k++)
                result.data[i][j] += data[i][k] * m.data[k][j];
    return result;
}

bool Matrix::operator==(const Matrix& m) const {
    if (rows != m.rows || cols != m.cols) return false;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (fabs(data[i][j] - m.data[i][j]) > 1e-9) return false;
    return true;
}

//  PROPERTY CHECKS
bool Matrix::isSquare() const {
    return rows == cols;
}

bool Matrix::isSymmetric() const {
    if (!isSquare()) return false;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (fabs(data[i][j] - data[j][i]) > 1e-9) return false;
    return true;
}

bool Matrix::isIdentity() const {
    if (!isSquare()) return false;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            if (fabs(data[i][j] - expected) > 1e-9) return false;
        }
    return true;
}

bool Matrix::isNull() const {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (fabs(data[i][j]) > 1e-9) return false;
    return true;
}

bool Matrix::isDiagonal() const {
    if (!isSquare()) return false;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (i != j && fabs(data[i][j]) > 1e-9) return false;
    return true;
}

bool Matrix::isDiagonallyDominant() const {
    if (!isSquare()) return false;
    for (int i = 0; i < rows; i++) {
        double diag = fabs(data[i][i]);
        double sum  = 0.0;
        for (int j = 0; j < cols; j++)
            if (j != i) sum += fabs(data[i][j]);
        if (diag < sum) return false;
    }
    return true;
}

bool Matrix::isTranspose(const Matrix& m) const {
    if (rows != m.cols || cols != m.rows) return false;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (fabs(data[i][j] - m.data[j][i]) > 1e-9) return false;
    return true;
}

//  TRANSFORMATIONS
void Matrix::makeDiagonallyDominant() {
    for (int col = 0; col < rows; col++) {
        int maxRow = col;
        for (int row = col + 1; row < rows; row++)
            if (fabs(data[row][col]) > fabs(data[maxRow][col]))
                maxRow = row;
        if (maxRow != col)
            swap(data[col], data[maxRow]);
    }
}

Matrix Matrix::transpose() const {
    Matrix t(cols, rows);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            t.data[j][i] = data[i][j];
    return t;
}

//  DETERMINANT  (recursive cofactor expansion)
double Matrix::determinant() const {
    if (rows != cols) return 0.0;
    if (rows == 1) return data[0][0];
    if (rows == 2)
        return data[0][0]*data[1][1] - data[0][1]*data[1][0];

    double det = 0.0;
    for (int p = 0; p < cols; p++) {
        Matrix sub(rows-1, cols-1);
        for (int i = 1; i < rows; i++) {
            int ci = 0;
            for (int j = 0; j < cols; j++) {
                if (j == p) continue;
                sub.data[i-1][ci++] = data[i][j];
            }
        }
        det += data[0][p] * sub.determinant() * ((p % 2 == 0) ? 1 : -1);
    }
    return det;
}

// INVERSE  (Gauss-Jordan on augmented [A | I])
Matrix Matrix::inverse() const {
    if (!isSquare())
        throw invalid_argument("Inverse requires a square matrix.");

    double det = determinant();
    if (fabs(det) < 1e-12)
        throw invalid_argument("Matrix is singular (determinant = 0). Inverse does not exist.");

    int n = rows;

    // Build augmented [A | I]
    vector<vector<double>> aug(n, vector<double>(2*n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            aug[i][j] = data[i][j];
        aug[i][n + i] = 1.0;   // identity part
    }

    // Forward elimination with partial pivoting
    for (int col = 0; col < n; col++) {
        // Find pivot row
        int pivot = col;
        for (int row = col+1; row < n; row++)
            if (fabs(aug[row][col]) > fabs(aug[pivot][col]))
                pivot = row;
        swap(aug[col], aug[pivot]);

        double pivotVal = aug[col][col];
        // Normalize pivot row
        for (int j = 0; j < 2*n; j++)
            aug[col][j] /= pivotVal;

        // Eliminate all other rows
        for (int row = 0; row < n; row++) {
            if (row == col) continue;
            double factor = aug[row][col];
            for (int j = 0; j < 2*n; j++)
                aug[row][j] -= factor * aug[col][j];
        }
    }

    // Extract right half as inverse
    Matrix inv(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv.data[i][j] = aug[i][n + j];

    return inv;
}











































































// #include "../include/Matrix.hpp"

// using namespace std;

// Matrix::Matrix() {
//     rows = 0;
//     cols = 0;
// }

// Matrix::Matrix(int r,int c) {
//     rows = r;
//     cols = c;
//     data.resize(r, vector<double>(c,0));
// }
// //
// Matrix::Matrix(const Matrix &m) {
//     rows = m.rows;
//     cols = m.cols;
//     data = m.data;
// }

// void Matrix::readFromFile(ifstream &fin) 
// {
//     for(int i=0;i<rows;i++)
//         for(int j=0;j<cols;j++)
//             fin >> data[i][j];
// }

// void Matrix::displayToFile(ofstream &fout) const {

//     for(int i=0;i<rows;i++) {
//         for(int j=0;j<cols;j++) {
//             fout << data[i][j] << " ";
//         }
//         fout << endl;
//     }
// }

// double& Matrix::operator()(int i,int j) {
//     return data[i][j];
// }

// double Matrix::operator()(int i,int j) const {
//     return data[i][j];
// }

// istream& operator>>(istream &in, Matrix &m) {

//     for(int i=0;i<m.rows;i++)
//         for(int j=0;j<m.cols;j++)
//             in >> m.data[i][j];

//     return in;
// }

// ostream& operator<<(ostream &out, const Matrix &m)
// {
//     for(int i = 0; i < m.rows; i++)
//     {
//         for(int j = 0; j < m.cols; j++)
//         {
//             out << m.data[i][j];

//             if(j < m.cols - 1)
//                 out << " ";
//         }
//         out << endl;
//     }

//     return out;
// }

// Matrix Matrix::transpose() const {

//     Matrix t(cols,rows);

//     for(int i=0;i<rows;i++)
//         for(int j=0;j<cols;j++)
//             t.data[j][i] = data[i][j];

//     return t;
// }

// bool Matrix::isSquare() const {
//     return rows == cols;
// }

// double Matrix::determinant() const {

//     if(rows != cols) return 0;

//     if(rows == 1)
//         return data[0][0];

//     if(rows == 2)
//         return data[0][0]*data[1][1] - data[0][1]*data[1][0];

//     double det = 0;

//     for(int p=0;p<cols;p++) {

//         Matrix sub(rows-1, cols-1);

//         for(int i=1;i<rows;i++) {

//             int colIndex = 0;

//             for(int j=0;j<cols;j++) {

//                 if(j == p) continue;

//                 sub.data[i-1][colIndex] = data[i][j];
//                 colIndex++;
//             }
//         }

//         det += data[0][p] * sub.determinant() * ((p%2==0)?1:-1);
//     }

//     return det;
// }