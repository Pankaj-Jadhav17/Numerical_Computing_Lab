#include "../include/Interpolation.hpp"
#include <stdexcept>
using namespace std;

Interpolation::Interpolation(const vector<double>& xs, const vector<double>& ys)
    : Matrix((int)xs.size(), 2)   // rows = n, cols = 2 — Matrix allocates data[][]
{
    // Validate input
    if (xs.size() != ys.size() || xs.empty())
        throw invalid_argument("Interpolation: x and y must be non-empty and same length.");

    // Fill Matrix data[][] directly
    // rows is now set to xs.size() by Matrix constructor
    for (int i = 0; i < rows; ++i) {
        data[i][0] = xs[i];   // col 0 = x value
        data[i][1] = ys[i];   // col 1 = y value
    }
}