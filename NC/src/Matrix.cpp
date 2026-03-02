#include "../include/Matrix.hpp"
#include <algorithm>
#include <stdexcept>  // for std::out_of_range
#include <fstream>   // for file I/O

using namespace std;           // bring std names into the cpp file

Matrix::Matrix(int r, int c)
{
    rows = r;
    cols = c;
    data.resize(r, vector<double>(c, 0.0));
}

int Matrix::getRows() const { return rows; }
int Matrix::getCols() const { return cols; }

double Matrix::get(int i, int j) const
{
    if (i < 0 || i >= rows || j < 0 || j >= cols)
        throw out_of_range("Matrix index out of range");
    return data[i][j];
}

void Matrix::set(int i, int j, double value)
{
    if (i < 0 || i >= rows || j < 0 || j >= cols)
        throw out_of_range("Matrix index out of range");
    data[i][j] = value;
}

void Matrix::swapRows(int r1, int r2)
{
    if (r1 < 0 || r1 >= rows || r2 < 0 || r2 >= rows)
        throw out_of_range("Row index out of range");
    swap(data[r1], data[r2]);
}

// read the full matrix (including augmented column) from a text file
void Matrix::readFromFile(const string &file)
{
    ifstream fin(file);
    if (!fin)
        throw runtime_error("Unable to open file " + file);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            if (!(fin >> data[i][j]))
                throw runtime_error("Insufficient data in file " + file);
        }
    }
}

// write solution vector to file (one value per line)
void Matrix::writeSolution(const string &file, const vector<double> &x)
{
    ofstream fout(file);
    if (!fout)
        throw runtime_error("Unable to open file " + file);

    for (double v : x)
        fout << v << "\n";
}







