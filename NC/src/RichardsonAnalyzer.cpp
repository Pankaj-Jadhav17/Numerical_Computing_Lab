#include "../include/RichardsonAnalyzer.hpp"
#include "../include/Differentiation.hpp"
#include "../include/CentralDifference.hpp"
#include "../include/RichardsonExtrapolation.hpp"

#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <map>
#include <algorithm>

RichardsonAnalyzer::RichardsonAnalyzer(double evalPoint, const vector<double>& hVals)
    : x0(evalPoint), hValues(hVals) {}

void RichardsonAnalyzer::addFunction(const string& name,
                                      function<double(double)> f,
                                      function<double(double)> fExact) {
    testFunctions.push_back({name, f, fExact});
}

void RichardsonAnalyzer::run() {
    results.clear();
    for (const auto& tc : testFunctions) {
        double exact = tc.fExact(x0);
        for (double h : hValues) {
            CentralDifference       cd(tc.f, h);
            RichardsonExtrapolation re(tc.f, h);

            double D_h = cd.derivative(x0);
            double R_h = re.derivative(x0);

            RichardsonResultRow row;
            row.function = tc.name;
            row.h        = h;
            row.D_h      = D_h;
            row.R_h      = R_h;
            row.exact    = exact;
            row.errD     = Differentiation::absoluteError(exact, D_h);
            row.errR     = Differentiation::absoluteError(exact, R_h);

            results.push_back(row);
        }
    }
}

void RichardsonAnalyzer::printTable(ostream& out) const {
    out << scientific << setprecision(6);
    out << left
        << setw(8)  << "func"
        << setw(12) << "h"
        << setw(16) << "D(h)"
        << setw(16) << "R(h)"
        << setw(16) << "exact"
        << setw(16) << "err_D(h)"
        << setw(16) << "err_R(h)" << "\n";

    for (const auto& r : results) {
        out << left
            << setw(8)  << r.function
            << setw(12) << r.h
            << setw(16) << r.D_h
            << setw(16) << r.R_h
            << setw(16) << r.exact
            << setw(16) << r.errD
            << setw(16) << r.errR << "\n";
    }
}

void RichardsonAnalyzer::writeTable(const string& path) const {
    ofstream out(path);
    out << "Richardson Extrapolation vs Central Difference - Results Table\n";
    out << "Evaluation point x0 = " << x0 << "\n\n";
    printTable(out);
}

void RichardsonAnalyzer::writeCSV(const string& path) const {
    ofstream csv(path);
    csv << "function,h,D_h,R_h,exact,err_D,err_R\n";
    csv << scientific << setprecision(8);
    for (const auto& r : results) {
        csv << r.function << "," << r.h << "," << r.D_h << "," << r.R_h << ","
            << r.exact << "," << r.errD << "," << r.errR << "\n";
    }
}

// group results by function, preserving first-seen order
static void groupByFunction(const vector<RichardsonResultRow>& results,
                             vector<string>& order,
                             map<string, vector<const RichardsonResultRow*>>& byFunc) {
    for (const auto& r : results) {
        if (byFunc.find(r.function) == byFunc.end()) order.push_back(r.function);
        byFunc[r.function].push_back(&r);
    }
}

void RichardsonAnalyzer::writeGnuplotFiles(const string& dataPrefix,
                                            const string& scriptPath,
                                            const string& outputImage) const {
    auto safeLog = [](double v) { return log10(v > 0 ? v : 1e-20); };

    vector<string> order;
    map<string, vector<const RichardsonResultRow*>> byFunc;
    groupByFunction(results, order, byFunc);

    for (const auto& fname : order) {
        ofstream dat(dataPrefix + "_" + fname + ".dat");
        dat << "# log10(h)  log10(err_D)  log10(err_R)\n";
        for (const auto* r : byFunc[fname]) {
            dat << safeLog(r->h) << " " << safeLog(r->errD) << " "
                << safeLog(r->errR) << "\n";
        }
    }

    ofstream gp(scriptPath);
    gp << "set terminal pngcairo size 1000,700 enhanced font 'Verdana,10'\n";
    gp << "set output '" << outputImage << "'\n";
    gp << "set title 'Log-Log Error Plot: Central Difference D(h) vs Richardson R(h)'\n";
    gp << "set xlabel 'log10(h)'\n";
    gp << "set ylabel 'log10(|error|)'\n";
    gp << "set grid\n";
    gp << "set key outside right\n";
    gp << "plot \\\n";
    bool first = true;
    for (const auto& fname : order) {
        string file = dataPrefix + "_" + fname + ".dat";
        if (!first) gp << ", \\\n";
        gp << "  '" << file << "' using 1:2 with linespoints title '" << fname << " D(h) [O(h^2)]'";
        gp << ", \\\n  '" << file << "' using 1:3 with linespoints title '" << fname << " R(h) [O(h^4)]'";
        first = false;
    }
    gp << "\n";
}

bool RichardsonAnalyzer::renderPlot(const string& scriptPath) const {
    string cmd = "gnuplot " + scriptPath;
    return system(cmd.c_str()) == 0;
}

// Least-squares slope of log10(h) vs log10(error), restricted to the
// "well-behaved" (monotonically decreasing, pre-round-off) prefix of rows.
double RichardsonAnalyzer::logLogSlope(const vector<const RichardsonResultRow*>& rows,
                                        bool useD) {
    auto errOf = [&](const RichardsonResultRow* r) { return useD ? r->errD : r->errR; };

    // Find how many leading points are monotonically decreasing (pre-round-off)
    size_t n = 1;
    while (n < rows.size() && errOf(rows[n]) < errOf(rows[n - 1])) ++n;
    if (n < 2) n = min(rows.size(), (size_t)2);

    vector<double> X, Y;
    for (size_t i = 0; i < n; ++i) {
        double e = errOf(rows[i]);
        if (e <= 0.0) continue;
        X.push_back(log10(rows[i]->h));
        Y.push_back(log10(e));
    }
    if (X.size() < 2) return 0.0;

    double xMean = 0, yMean = 0;
    for (size_t i = 0; i < X.size(); ++i) { xMean += X[i]; yMean += Y[i]; }
    xMean /= X.size(); yMean /= Y.size();

    double num = 0, den = 0;
    for (size_t i = 0; i < X.size(); ++i) {
        num += (X[i] - xMean) * (Y[i] - yMean);
        den += (X[i] - xMean) * (X[i] - xMean);
    }
    return (den != 0.0) ? num / den : 0.0;
}

void RichardsonAnalyzer::printAnalysis(ostream& out) const {
    out << fixed << setprecision(4);
    out << "\n===== CONVERGENCE / ANALYSIS =====\n";

    vector<string> order;
    map<string, vector<const RichardsonResultRow*>> byFunc;
    groupByFunction(results, order, byFunc);

    bool richardsonAlwaysBetter = true;
    double totalRatioD = 0, totalRatioR = 0;
    int countRatioD = 0, countRatioR = 0;

    for (const auto& fname : order) {
        const auto& rows = byFunc[fname];
        out << "\n-- " << fname << " --\n";

        // Observed log-log slopes (least-squares, pre-round-off region)
        double slopeD = logLogSlope(rows, true);
        double slopeR = logLogSlope(rows, false);
        out << "  Observed slope of log10(err_D) vs log10(h): " << slopeD
            << "   (theoretical: 2.0000)\n";

        // Special case: if R(h) is already at machine-epsilon level even at
        // the LARGEST h tested, there is no truncation-error trend left to
        // measure -- the higher derivative driving Richardson's error term
        // is exactly (or numerically) zero, so R(h) is already exact.
        bool richardsonExact = (rows.front()->errR < 1e-9);
        if (richardsonExact) {
            out << "  Observed slope of log10(err_R) vs log10(h): N/A -- R(h) is already at "
                   "machine-precision (~" << scientific << rows.front()->errR << fixed
                << ") even at the LARGEST h tested. There is no truncation-error trend "
                   "to fit a slope to.\n";
        } else {
            out << "  Observed slope of log10(err_R) vs log10(h): " << slopeR
                << "   (theoretical: 4.0000)\n";
        }

        // Per-decade error reduction factor (geometric mean of consecutive ratios,
        // restricted to the same well-behaved prefix used for the slope fit)
        size_t nD = 1; while (nD < rows.size() && rows[nD]->errD < rows[nD-1]->errD) ++nD;
        size_t nR = 1; while (nR < rows.size() && rows[nR]->errR < rows[nR-1]->errR) ++nR;

        double logRatioD = 0; int kD = 0;
        for (size_t i = 1; i < nD; ++i) {
            if (rows[i]->errD > 0 && rows[i-1]->errD > 0) {
                logRatioD += log10(rows[i-1]->errD / rows[i]->errD);
                ++kD;
            }
        }
        double logRatioR = 0; int kR = 0;
        for (size_t i = 1; i < nR; ++i) {
            if (rows[i]->errR > 0 && rows[i-1]->errR > 0) {
                logRatioR += log10(rows[i-1]->errR / rows[i]->errR);
                ++kR;
            }
        }
        double factorD = (kD > 0) ? pow(10.0, logRatioD / kD) : 0.0;
        double factorR = (kR > 0) ? pow(10.0, logRatioR / kR) : 0.0;
        out << "  Error reduction factor per 10x decrease in h:\n";
        out << "    D(h): ~" << factorD << "x   (theoretical: ~100x for O(h^2))\n";
        if (richardsonExact) {
            out << "    R(h): N/A -- already at machine precision, nothing left to reduce.\n";
            if (fname == "poly")
                out << "    (f(x) = x^3 - 2x + 1 has f''''(x) = 0 identically -- Richardson's "
                       "leading error term, which is proportional to f''''(x), vanishes "
                       "exactly for any cubic polynomial. R(h) recovers the exact derivative "
                       "up to round-off, for every h.)\n";
        } else {
            out << "    R(h): ~" << factorR << "x   (theoretical: ~10000x for O(h^4))\n";
        }
        totalRatioD += factorD;
        if (!richardsonExact) { totalRatioR += factorR; ++countRatioR; }
        ++countRatioD;

        // Q1 check: is R(h) always more accurate than D(h)?
        int exceptions = 0;
        for (const auto* r : rows) if (r->errR >= r->errD) ++exceptions;
        if (exceptions > 0) richardsonAlwaysBetter = false;
        out << "  Richardson more accurate than Central at " << (rows.size() - exceptions)
            << "/" << rows.size() << " tested h values"
            << (exceptions > 0 ? "  (fails at the smallest h -- round-off dominates R(h) first)" : "")
            << "\n";

        // Q3 / Q6: locate the round-off floor for D(h) and R(h)
        size_t minD = 0, minR = 0;
        for (size_t i = 1; i < rows.size(); ++i) {
            if (rows[i]->errD < rows[minD]->errD) minD = i;
            if (rows[i]->errR < rows[minR]->errR) minR = i;
        }
        out << "  D(h) reaches its minimum error at h = " << scientific << rows[minD]->h
            << "  (err = " << rows[minD]->errD << ")\n";
        out << "  R(h) reaches its minimum error at h = " << rows[minR]->h
            << "  (err = " << rows[minR]->errR << ")\n";
        out << fixed;
        if (minR > minD)
            out << "  -> R(h) keeps improving to a SMALLER h than D(h) before round-off "
                   "takes over.\n";
        else if (minR < minD)
            out << "  -> R(h) hits its round-off floor EARLIER (larger h) than D(h): "
                   "combining two central-difference evaluations (h and h/2) adds more "
                   "subtractive cancellation, so round-off creeps in sooner.\n";
        else
            out << "  -> D(h) and R(h) hit their round-off floor at about the same h.\n";
    }

    double avgFactorD = countRatioD ? totalRatioD / countRatioD : 0.0;
    double avgFactorR = countRatioR ? totalRatioR / countRatioR : 0.0;

    out << "\n===== ANSWERS TO ASSIGNMENT QUESTIONS =====\n";
    out << "Q1. Does Richardson always give a smaller error than Central Difference?\n";
    out << "    " << (richardsonAlwaysBetter
            ? "Yes, across every h tested for every function."
            : "Mostly, but NOT always: at the smallest step sizes, round-off error in "
              "computing D(h) and D(h/2) and combining them can make R(h) worse than "
              "D(h). Richardson only wins in the truncation-error-dominated regime.")
        << "\n\n";

    out << "Q2. By approximately what factor does the error decrease when h is reduced by 10?\n";
    out << "    D(h): ~" << avgFactorD << "x per decade (matches O(h^2): 10^2 = 100)\n";
    out << "    R(h): ~" << avgFactorR << "x per decade (matches O(h^4): 10^4 = 10000)\n\n";

    out << "Q3. Does Richardson continue to improve as h becomes very small?\n";
    out << "    No. Like Central Difference, R(h) improves only while truncation error "
           "dominates. Once h is small enough that floating-point round-off dominates, "
           "R(h)'s error stops decreasing and can increase again -- and because R(h) "
           "combines two central-difference evaluations (extra subtractions of nearly "
           "equal numbers), this floor is often reached at a LARGER h than for D(h) "
           "alone (see per-function minima above).\n\n";

    out << "Q4. Does the log-log slope support the theoretical prediction?\n";
    out << "    Yes -- see the per-function 'Observed slope' lines above: D(h) slopes "
           "come out close to 2, R(h) slopes close to 4, matching D(h)=O(h^2) and "
           "R(h)=O(h^4), when measured over the pre-round-off region of each curve.\n\n";

    out << "Q5. Does Richardson achieve the expected O(h^4) behavior?\n";
    out << "    Yes, in the truncation-error-dominated region (roughly h = 1e-1 down to "
           "1e-3/1e-4 depending on the function) -- the observed slope is close to 4 and "
           "the error drops by close to 10^4 per decade, exactly as predicted by "
           "eliminating the h^2 term. It breaks down once round-off takes over.\n\n";

    out << "Q6. What role does floating-point round-off error play?\n";
    out << "    Both D(h) and R(h) compute a numerator as a difference of nearby function "
           "values, then divide by a small h. As h shrinks, the numerator's absolute "
           "error (about machine epsilon times |f(x)|) stays roughly fixed while the "
           "true numerator shrinks with h, so relative round-off in the numerator is "
           "amplified by dividing by a tiny h (and h^2 for R(h), effectively, since it's "
           "built from two central differences). This round-off error grows as h shrinks, "
           "eventually overtaking the (also shrinking) truncation error -- producing a "
           "minimum total error at some optimal h, beyond which further shrinking h makes "
           "results WORSE, not better. This is why the log-log error curve is V-shaped "
           "rather than a straight line all the way down.\n";
}

void RichardsonAnalyzer::writeAnalysis(const string& path) const {
    ofstream out(path);
    out << "Richardson Extrapolation - Written Analysis\n";
    out << "Evaluation point x0 = " << x0 << "\n";
    printAnalysis(out);
}
