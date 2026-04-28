#include "../include/CurveFitting.hpp"
#include <stdexcept>
using namespace std;

CurveFitting::CurveFitting(const Matrix& m) : Matrix(m)
{
    if (cols != 2)
        throw invalid_argument(
            "CurveFitting: Matrix must have exactly 2 columns "
            "(col 0 = x values, col 1 = f values).");

    if (rows < 2)
        throw invalid_argument(
            "CurveFitting: Need at least 2 data points (rows).");
}

// Out-of-line virtual destructor — anchors the vtable in this translation unit.
CurveFitting::~CurveFitting() = default;