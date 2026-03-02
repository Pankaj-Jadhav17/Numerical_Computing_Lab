#include "../include/matrix.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

Matrix readFromFile(const string &file)
{
    ifstream fin(file);
    if (!fin.is_open())
        throw runtime_error("Cannot open input file: " + file);

    string line;
    int n = 0;

    // Skip comment lines (starting with '#') and read n
    while (getline(fin, line))
    {
        if (line.empty() || line[0] == '#')
            continue;

        istringstream iss(line);
        if (iss >> n)
            break;
    }

    if (n <= 0)
        throw runtime_error("Invalid or missing matrix size in file: " + file);

    Matrix A(n, n + 1); 
    int row = 0;

    // Read n rows, skipping any comment lines in between
    while (getline(fin, line) && row < n)
    {
        if (line.empty() || line[0] == '#')
            continue;

        istringstream iss(line);
        for (int j = 0; j < n + 1; j++)
        {
            double val;
            if (!(iss >> val))
                throw runtime_error("Not enough values in row " + to_string(row + 1)
                                    + " of file: " + file);
            A.set(row, j, val);
        }
        row++;
    }

    if (row < n)
        throw runtime_error("File ended before all rows were read: " + file);

    return A;
}

//Write solution to file 
void writeToFile(const string &file, vector<double> &x)
{
    ofstream fout(file);
    if (!fout.is_open())
        throw runtime_error("Cannot open output file: " + file);

    for (int i = 0; i < (int)x.size(); i++)
        fout << "x" << i + 1 << " = " << x[i] << "\n";
}



























// #include "../include/matrix.hpp"
// #include <fstream>

// Matrix readFromFile(string file)
// {
//     ifstream fin(file);

//     int n;
//     fin >> n;

//     Matrix A(n, n + 1);

//     for (int i = 0; i < n; i++)
//         for (int j = 0; j < n + 1; j++)
//         {
//             double val;
//             fin >> val;
//             A.set(i, j, val);
//         }

//     return A;
// }

// void writeToFile(string file, vector<double> &x)
// {
//     ofstream fout(file);

//     for (int i = 0; i < x.size(); i++)
//         fout << "x" << i + 1 << " = " << x[i] << endl;
// }