#include "../include/Matrix.hpp"

using namespace std;

Matrix::Matrix() {
    rows = 0;
    cols = 0;
}

Matrix::Matrix(int r,int c) {
    rows = r;
    cols = c;
    data.resize(r, vector<double>(c,0));
}

Matrix::Matrix(const Matrix &m) {
    rows = m.rows;
    cols = m.cols;
    data = m.data;
}

void Matrix::readFromFile(ifstream &fin) {

    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
            fin >> data[i][j];
}

void Matrix::displayToFile(ofstream &fout) const {

    for(int i=0;i<rows;i++) {
        for(int j=0;j<cols;j++) {
            fout << data[i][j] << " ";
        }
        fout << endl;
    }
}

double& Matrix::operator()(int i,int j) {
    return data[i][j];
}

double Matrix::operator()(int i,int j) const {
    return data[i][j];
}

istream& operator>>(istream &in, Matrix &m) {

    for(int i=0;i<m.rows;i++)
        for(int j=0;j<m.cols;j++)
            in >> m.data[i][j];

    return in;
}

ostream& operator<<(ostream &out, const Matrix &m)
{
    for(int i = 0; i < m.rows; i++)
    {
        for(int j = 0; j < m.cols; j++)
        {
            out << m.data[i][j];

            if(j < m.cols - 1)
                out << " ";
        }
        out << endl;
    }

    return out;
}

Matrix Matrix::transpose() const {

    Matrix t(cols,rows);

    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
            t.data[j][i] = data[i][j];

    return t;
}

bool Matrix::isSquare() const {
    return rows == cols;
}

double Matrix::determinant() const {

    if(rows != cols) return 0;

    if(rows == 1)
        return data[0][0];

    if(rows == 2)
        return data[0][0]*data[1][1] - data[0][1]*data[1][0];

    double det = 0;

    for(int p=0;p<cols;p++) {

        Matrix sub(rows-1, cols-1);

        for(int i=1;i<rows;i++) {

            int colIndex = 0;

            for(int j=0;j<cols;j++) {

                if(j == p) continue;

                sub.data[i-1][colIndex] = data[i][j];
                colIndex++;
            }
        }

        det += data[0][p] * sub.determinant() * ((p%2==0)?1:-1);
    }

    return det;
}