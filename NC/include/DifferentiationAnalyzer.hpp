#ifndef DIFFERENTIATION_ANALYZER_HPP
#define DIFFERENTIATION_ANALYZER_HPP

#include <vector>
#include <string>
#include <functional>
#include <ostream>
using namespace std;

// One row of the comparison table for a (function, h) pair.
struct DiffResultRow {
    string function;
    double h;
    double forward, backward, central;
    double exact;
    double errForward, errBackward, errCentral;
};

// Registers test functions, runs Forward/Backward/Central difference
// (via the Differentiation hierarchy) across a list of step sizes h,
// and reports results as a table / CSV / gnuplot log-log plot.
//
// This class does NOT inherit Matrix or Differentiation — it is a
// driver/orchestrator, analogous to the free functions in main.cpp
// (runSolve, runEigenvalues, runLeastSquares) but packaged as a
// reusable class since it holds state (results) across multiple calls.
class DifferentiationAnalyzer {
private:
    struct TestCase {
        string name;
        function<double(double)> f;
        function<double(double)> fExact;
    };

    double x0;                       // evaluation point
    vector<double> hValues;          // step sizes to test
    vector<TestCase> testFunctions;  // registered functions
    vector<DiffResultRow> results;   // filled in by run()

public:
    DifferentiationAnalyzer(double evalPoint, const vector<double>& hVals);

    // Register a function to test: name, f(x), and its EXACT derivative f'(x)
    void addFunction(const string& name,
                      function<double(double)> f,
                      function<double(double)> fExact);

    // Runs Forward/Backward/Central difference for every function & every h
    void run();

    const vector<DiffResultRow>& getResults() const { return results; }

    // ---- Output ----
    void printTable(ostream& out) const;
    void writeCSV(const string& path) const;

    // Writes one .dat file per function (log10(h), log10(err_fwd/bwd/ctr))
    // plus a plot.gp gnuplot script, matching data file naming "prefix_<fn>.dat"
    void writeGnuplotFiles(const string& dataPrefix,
                            const string& scriptPath,
                            const string& outputImage) const;

    // Runs "gnuplot <scriptPath>"; returns true on success (false if gnuplot missing)
    bool renderPlot(const string& scriptPath) const;

    // ---- Analysis (answers the assignment's three questions from real data) ----
    // Observed convergence order between two consecutive h values:
    // order = log(err1/err2) / log(h1/h2)
    static double observedOrder(double h1, double err1, double h2, double err2);

    // Prints a written analysis: most accurate method, observed convergence
    // rate per method, and the h beyond which error stops decreasing (round-off floor).
    void printAnalysis(ostream& out) const;
};

#endif
