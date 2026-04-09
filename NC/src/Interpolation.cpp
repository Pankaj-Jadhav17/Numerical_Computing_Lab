#include "../include/Interpolation.hpp"
#include <stdexcept>
using namespace std;            

Interpolation::Interpolation(const vector<double>& xs, const vector<double>& ys)
    : Matrix((int)xs.size(), 2),
      n((int)xs.size()),
      xData(xs),
      yData(ys)
{
    if (xs.size() != ys.size() || xs.empty())
        throw invalid_argument("Interpolation: x and y must be non-empty and same length.");

    for (int i = 0; i < n; ++i) {
        data[i][0] = xs[i];
        data[i][1] = ys[i];
    }
}