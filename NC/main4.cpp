#include <iostream>
#include "include/GE.hpp"

using namespace std;

int main()
{
    int n = 4;

    SLE *solver = new GE(n, true);

    solver->readFromFile("input/system_4x4.txt");

    vector<double> solution = solver->solve();

    solver->writeSolution("output/result.txt", solution);

    delete solver;
}