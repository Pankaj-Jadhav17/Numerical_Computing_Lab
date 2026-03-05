#include <iostream>
#include <fstream>
#include <memory>   // for std::unique_ptr 
#include "include/GaussElimination.hpp"
#include "include/LU.hpp"

int main() {
    std::ifstream inFile("input.txt");
    std::ofstream outFile("output.txt");

    if (!inFile) {
        std::cerr << "Error: Could not open input.txt\n";
        return 1;
    }
    if (!outFile) {
        std::cerr << "Error: Could not open output.txt\n";
        return 1;
    }

    int n;
    inFile >> n;
    if (n <= 0) {
        std::cerr << "Error: Matrix size must be a positive integer.\n";
        return 1;
    }

    int method;
    inFile >> method;

    std::unique_ptr<SLE> solver;

    if (method == 1) {
        int pivotInput;
        inFile >> pivotInput;
        bool pivot = (pivotInput != 0);
        solver = std::make_unique<GaussElimination>(n, n + 1, pivot);
    }
    else if (method == 2) {
        solver = std::make_unique<LU>(n, n + 1);
    }
    else {
        std::cerr << "Error: Invalid method. Use 1 (Gauss) or 2 (LU).\n";
        return 1;
    }

    solver->readFromFile(inFile);
    solver->solve(outFile);

    std::cout << "Done! Results written to output.txt\n";

    // No need to call delete — unique_ptr cleans up automatically
    return 0;
}