#ifndef RICHARDSON_ANALYZER_HPP
#define RICHARDSON_ANALYZER_HPP

#include <vector>
#include <string>
#include <functional>
#include <ostream>
using namespace std;

// One row of the comparison table for a (function, h) pair.
struct RichardsonResultRow {
    string function;
    double h;
    double D_h;   // Central Difference approximation
    double R_h;   // Richardson Extrapolation approximation
    double exact;
    double errD;  // |exact - D_h|
    double errR;  // |exact - R_h|
};

// Registers test functions, runs Central Difference (D(h)) and
// Richardson Extrapolation (R(h)) across a list of step sizes h,
// and reports results as a table / CSV / gnuplot log-log plot,
// plus a computed convergence-order + written analysis.
//
// Standalone driver class (does not inherit Matrix or Differentiation),
// analogous to DifferentiationAnalyzer from Assignment II.
class RichardsonAnalyzer {
private:
    struct TestCase {
        string name;
        function<double(double)> f;
        function<double(double)> fExact;
    };

    double x0;
    vector<double> hValues;
    vector<TestCase> testFunctions;
    vector<RichardsonResultRow> results;

public:
    RichardsonAnalyzer(double evalPoint, const vector<double>& hVals);

    void addFunction(const string& name,
                      function<double(double)> f,
                      function<double(double)> fExact);

    // Runs Central Difference + Richardson Extrapolation for every function & h
    void run();

    const vector<RichardsonResultRow>& getResults() const { return results; }

    // ---- Output ----
    void printTable(ostream& out) const;
    void writeTable(const string& path) const;
    void writeCSV(const string& path) const;
    void writeGnuplotFiles(const string& dataPrefix,
                            const string& scriptPath,
                            const string& outputImage) const;
    bool renderPlot(const string& scriptPath) const;

    // ---- Convergence analysis ----
    // Least-squares slope of log10(h) vs log10(error), using only the
    // "well-behaved" (pre-round-off) points for a given function/column.
    // colSelector picks errD or errR from a row.
    static double logLogSlope(const vector<const RichardsonResultRow*>& rows, bool useD);

    // Prints the full written analysis, answering Q1-Q6 from the
    // assignment sheet directly from the computed data.
    void printAnalysis(ostream& out) const;
    void writeAnalysis(const string& path) const;
};

#endif
