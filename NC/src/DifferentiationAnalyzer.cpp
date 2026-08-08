#include "../include/DifferentiationAnalyzer.hpp"
#include "../include/ForwardDifference.hpp"
#include "../include/BackwardDifference.hpp"
#include "../include/CentralDifference.hpp"

#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <map>
#include <algorithm>

DifferentiationAnalyzer::DifferentiationAnalyzer(double evalPoint, const vector<double>& hVals)
    : x0(evalPoint), hValues(hVals) {}

void DifferentiationAnalyzer::addFunction(const string& name,
                                           function<double(double)> f,
                                           function<double(double)> fExact) {
    testFunctions.push_back({name, f, fExact});
}

void DifferentiationAnalyzer::run() {
    results.clear();
    for (const auto& tc : testFunctions) {
        double exact = tc.fExact(x0);
        for (double h : hValues) {
            ForwardDifference  fwd(tc.f, h);
            BackwardDifference bwd(tc.f, h);
            CentralDifference  ctr(tc.f, h);

            double fd = fwd.derivative(x0);
            double bd = bwd.derivative(x0);
            double cd = ctr.derivative(x0);

            DiffResultRow row;
            row.function    = tc.name;
            row.h           = h;
            row.forward     = fd;
            row.backward    = bd;
            row.central     = cd;
            row.exact       = exact;
            row.errForward  = Differentiation::absoluteError(exact, fd);
            row.errBackward = Differentiation::absoluteError(exact, bd);
            row.errCentral  = Differentiation::absoluteError(exact, cd);

            results.push_back(row);
        }
    }
}

void DifferentiationAnalyzer::printTable(ostream& out) const {
    out << scientific << setprecision(6);
    out << left
        << setw(8)  << "func"
        << setw(12) << "h"
        << setw(16) << "forward"
        << setw(16) << "backward"
        << setw(16) << "central"
        << setw(16) << "err_fwd"
        << setw(16) << "err_bwd"
        << setw(16) << "err_ctr" << "\n";

    for (const auto& r : results) {
        out << left
            << setw(8)  << r.function
            << setw(12) << r.h
            << setw(16) << r.forward
            << setw(16) << r.backward
            << setw(16) << r.central
            << setw(16) << r.errForward
            << setw(16) << r.errBackward
            << setw(16) << r.errCentral << "\n";
    }
}

void DifferentiationAnalyzer::writeCSV(const string& path) const {
    ofstream csv(path);
    csv << "function,h,forward,backward,central,exact,"
           "err_forward,err_backward,err_central\n";
    csv << scientific << setprecision(8);
    for (const auto& r : results) {
        csv << r.function << "," << r.h << "," << r.forward << "," << r.backward << ","
            << r.central << "," << r.exact << "," << r.errForward << ","
            << r.errBackward << "," << r.errCentral << "\n";
    }
}

void DifferentiationAnalyzer::writeGnuplotFiles(const string& dataPrefix,
                                                 const string& scriptPath,
                                                 const string& outputImage) const {
    auto safeLog = [](double v) { return log10(v > 0 ? v : 1e-20); };

    // group results by function name, preserving first-seen order
    vector<string> order;
    map<string, vector<const DiffResultRow*>> byFunc;
    for (const auto& r : results) {
        if (byFunc.find(r.function) == byFunc.end()) order.push_back(r.function);
        byFunc[r.function].push_back(&r);
    }

    for (const auto& fname : order) {
        ofstream dat(dataPrefix + "_" + fname + ".dat");
        dat << "# log10(h)  err_forward  err_backward  err_central\n";
        for (const auto* r : byFunc[fname]) {
            dat << safeLog(r->h) << " " << safeLog(r->errForward) << " "
                << safeLog(r->errBackward) << " " << safeLog(r->errCentral) << "\n";
        }
    }

    ofstream gp(scriptPath);
    gp << "set terminal pngcairo size 900,650 enhanced font 'Verdana,10'\n";
    gp << "set output '" << outputImage << "'\n";
    gp << "set title 'Log-Log Error Plot: Forward/Backward/Central Difference'\n";
    gp << "set xlabel 'log10(h)'\n";
    gp << "set ylabel 'log10(|error|)'\n";
    gp << "set grid\n";
    gp << "set key outside right\n";
    gp << "plot \\\n";
    bool first = true;
    for (const auto& fname : order) {
        string file = dataPrefix + "_" + fname + ".dat";
        if (!first) gp << ", \\\n";
        gp << "  '" << file << "' using 1:2 with linespoints title '" << fname << " forward'";
        gp << ", \\\n  '" << file << "' using 1:3 with linespoints title '" << fname << " backward'";
        gp << ", \\\n  '" << file << "' using 1:4 with linespoints title '" << fname << " central'";
        first = false;
    }
    gp << "\n";
}

bool DifferentiationAnalyzer::renderPlot(const string& scriptPath) const {
    string cmd = "gnuplot " + scriptPath;
    return system(cmd.c_str()) == 0;
}

double DifferentiationAnalyzer::observedOrder(double h1, double err1, double h2, double err2) {
    if (err1 <= 0.0 || err2 <= 0.0 || h1 == h2) return 0.0;
    return log(err1 / err2) / log(h1 / h2);
}

void DifferentiationAnalyzer::printAnalysis(ostream& out) const {
    out << "\n===== ANALYSIS =====\n";

    // group by function, preserving order, and keep rows sorted by h descending
    // (results were generated in that order already, so just iterate)
    vector<string> order;
    map<string, vector<const DiffResultRow*>> byFunc;
    for (const auto& r : results) {
        if (byFunc.find(r.function) == byFunc.end()) order.push_back(r.function);
        byFunc[r.function].push_back(&r);
    }

    for (const auto& fname : order) {
        const auto& rows = byFunc[fname];
        out << "\n-- " << fname << " --\n";

        // Which method is most accurate at the smallest well-behaved h (not the floor)?
        // Use the middle h value as a representative, well-conditioned point.
        const DiffResultRow* mid = rows[rows.size() / 2];
        out << "  At h = " << mid->h << ":  err_fwd=" << mid->errForward
            << "  err_bwd=" << mid->errBackward
            << "  err_ctr=" << mid->errCentral << "\n";
        double best = min({mid->errForward, mid->errBackward, mid->errCentral});
        string bestName = (best == mid->errCentral) ? "Central" :
                           (best == mid->errForward) ? "Forward" : "Backward";
        out << "  -> Most accurate at this h: " << bestName << " difference\n";

        // Observed convergence order between consecutive h values (first pair,
        // before round-off takes over)
        if (rows.size() >= 2) {
            double p_fwd = observedOrder(rows[0]->h, rows[0]->errForward,
                                          rows[1]->h, rows[1]->errForward);
            double p_bwd = observedOrder(rows[0]->h, rows[0]->errBackward,
                                          rows[1]->h, rows[1]->errBackward);
            double p_ctr = observedOrder(rows[0]->h, rows[0]->errCentral,
                                          rows[1]->h, rows[1]->errCentral);
            out << "  Observed order (between h=" << rows[0]->h << " and h=" << rows[1]->h << "):\n";
            out << "    Forward  ~ O(h^" << p_fwd << ")  (expected ~1)\n";
            out << "    Backward ~ O(h^" << p_bwd << ")  (expected ~1)\n";
            out << "    Central  ~ O(h^" << p_ctr << ")  (expected ~2)\n";
        }

        // Find h beyond which central-difference error stops decreasing (round-off floor)
        for (size_t i = 1; i < rows.size(); ++i) {
            if (rows[i]->errCentral >= rows[i - 1]->errCentral) {
                out << "  Error stops decreasing (round-off floor) beyond h = "
                    << rows[i - 1]->h << " (central difference)\n";
                break;
            }
        }
    }
}
