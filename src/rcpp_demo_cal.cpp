#include <Rcpp.h>
#include <RcppEigen.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>
#include <cmath>

#include <thread>
#include <mutex>
/*
#define _DEBUG_
*/

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

class SolveM
{
public:
    Eigen::MatrixXd pvals;
    Eigen::MatrixXd abs_diffs;
    SolveM(const std::string &filename) : filename(filename)
    {
        read_data();
        precompute_design_matrix();
    }

    void compute()
    {
        size_t num_threads = std::thread::hardware_concurrency();
        std::vector<std::thread> threads;
        std::mutex mtx;

        for (size_t t = 0; t < num_threads; ++t)
        {
            threads.emplace_back([this, t, num_threads, &mtx]()
                                 {
                for (size_t i = t; i < unique_celltypes.size(); i += num_threads) {
                    const auto& celltype = unique_celltypes[i];
                    
                    Eigen::VectorXd local_pvals(df.size());
                    Eigen::VectorXd local_abs_diffs(df.size());

                    for (size_t k = 0; k < df.size(); ++k) {
                        local_pvals(k) = pval_f(df[k].values, i);
                        local_abs_diffs(k) = abs_diffs_r(df[k].values, i);
                    }

                    std::lock_guard<std::mutex> lock(mtx);
                    pvals.col(i) = local_pvals;
                    abs_diffs.col(i) = local_abs_diffs;
                } });
        }

        for (auto &thread : threads)
        {
            thread.join();
        }
    }

    void print_results()
    {
        std::cout << "P-values:\n"
                  << pvals << "\n\n";
        std::cout << "Absolute differences:\n"
                  << abs_diffs << "\n";
    }

private:
    std::vector<Eigen::MatrixXd> design_matrices;
    std::vector<std::vector<bool>> is_celltype;

    struct DataFrameRow
    {
        std::vector<double> values;
    };

    using DataFrame = std::vector<DataFrameRow>;

    std::string filename;
    DataFrame df;
    std::vector<std::string> celltypes;
    std::vector<std::string> unique_celltypes;
    std::map<std::string, size_t> celltype_index;

    void read_data()
    {
        std::ifstream file(filename);
        std::string line;

        // Read header
        std::getline(file, line);
        celltypes = get_celltypes(line);

        // Read data
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            std::string value;
            DataFrameRow row;

            // Skip first two columns
            std::getline(iss, value, '\t');
            std::getline(iss, value, '\t');

            while (std::getline(iss, value, '\t'))
            {
                row.values.push_back(std::stod(value));
            }
            df.push_back(row);
        }

        // Get unique celltypes
        std::set<std::string> unique_set(celltypes.begin(), celltypes.end());
        unique_celltypes = std::vector<std::string>(unique_set.begin(), unique_set.end());

        // Create celltype index map
        for (size_t i = 0; i < unique_celltypes.size(); ++i)
        {
            celltype_index[unique_celltypes[i]] = i;
        }

        // Initialize pvals and abs_diffs matrices
        pvals = Eigen::MatrixXd::Constant(df.size(), unique_celltypes.size(), NAN);
        abs_diffs = Eigen::MatrixXd::Constant(df.size(), unique_celltypes.size(), NAN);
    }

    std::vector<std::string> get_celltypes(const std::string &header_line)
    {
        std::istringstream iss(header_line);
        std::string cell;
        std::vector<std::string> celltypes;

        // Skip first two columns
        std::getline(iss, cell, '\t');
        std::getline(iss, cell, '\t');

        while (std::getline(iss, cell, '\t'))
        {
            size_t first_dot = cell.find('.');
            size_t second_dot = cell.find('.', first_dot + 1);
            if (first_dot != std::string::npos && second_dot != std::string::npos)
            {
                celltypes.push_back(cell.substr(first_dot + 1, second_dot - first_dot - 1));
            }
        }

        return celltypes;
    }

    void precompute_design_matrix()
    {
        design_matrices.resize(unique_celltypes.size());
        is_celltype.resize(unique_celltypes.size(), std::vector<bool>(celltypes.size()));

        for (size_t i = 0; i < unique_celltypes.size(); ++i)
        {
            const auto &celltype = unique_celltypes[i];

            // Create is_celltype vector
            for (size_t j = 0; j < celltypes.size(); ++j)
            {
                is_celltype[i][j] = (celltypes[j] == celltype);
            }

            // Create design matrix
            Eigen::MatrixXd X = Eigen::MatrixXd::Ones(celltypes.size(), unique_celltypes.size());
            for (size_t j = 0; j < unique_celltypes.size(); ++j)
            {
                if (j != i)
                { // Skip the current celltype (it's the baseline)
                    for (size_t k = 0; k < celltypes.size(); ++k)
                    {
                        X(k, j) = (celltypes[k] == unique_celltypes[j]) ? 1.0 : 0.0;
                    }
                }
            }
            design_matrices[i] = X;
        }
    }

    double pval_f(const std::vector<double> &x, size_t celltype_index)
    {
        Eigen::VectorXd y = Eigen::Map<const Eigen::VectorXd>(x.data(), x.size());

        // Use precomputed design matrix
        const Eigen::MatrixXd &X = design_matrices[celltype_index];

        // Fit linear model
        Eigen::VectorXd coeffs = (X.transpose() * X).ldlt().solve(X.transpose() * y);

        // Calculate residuals
        Eigen::VectorXd residuals = y - X * coeffs;

        // Calculate degrees of freedom and MSE
        int df = x.size() - coeffs.size();
        double mse = residuals.squaredNorm() / df;

        // Calculate standard errors
        Eigen::VectorXd se = ((X.transpose() * X).inverse().diagonal() * mse).cwiseSqrt();

        // Calculate t-statistics
        Eigen::VectorXd t_stats = coeffs.cwiseQuotient(se);
        Eigen::VectorXd t_stats_abs = t_stats.cwiseAbs();

        // Return the maximum p-value
        return t_stats_abs.minCoeff();
    }

    double pval_f2(const std::vector<double> &x, const std::vector<std::string> &tmp)
    {
        Eigen::VectorXd y = Eigen::Map<const Eigen::VectorXd>(x.data(), x.size());

        // Get unique categories (excluding 'A')
        std::set<std::string> categories;
        for (const auto &category : tmp)
        {
            if (category != "A")
            {
                categories.insert(category);
            }
        }

        // Create design matrix
        Eigen::MatrixXd X = Eigen::MatrixXd::Ones(x.size(), categories.size() + 1);
        int col = 1;
        for (const auto &category : categories)
        {
            for (size_t i = 0; i < tmp.size(); ++i)
            {
                X(i, col) = (tmp[i] == category) ? 1.0 : 0.0;
            }
            ++col;
        }
#ifdef _DEBUG_
        std::cout << "\033[0;31mCAT: \n"
                  << std::endl;
        for (auto categories_i : categories)
        {
            std::cout << categories_i << "\t";
        }
        std::cout << "\nX: \n"
                  << X << "\033[0m" << std::endl;
#endif

        // Fit linear model
        Eigen::VectorXd coeffs = (X.transpose() * X).ldlt().solve(X.transpose() * y);

        // Calculate residuals
        Eigen::VectorXd residuals = y - X * coeffs;

        // Calculate degrees of freedom and MSE
        int df = x.size() - coeffs.size();
        double mse = residuals.squaredNorm() / df;

        // Calculate standard errors
        Eigen::VectorXd se = ((X.transpose() * X).inverse().diagonal() * mse).cwiseSqrt();

        // Calculate t-statistics
        Eigen::VectorXd t_stats = coeffs.cwiseQuotient(se);
        Eigen::VectorXd t_stats_abs = t_stats.cwiseAbs();

        // Calculate p-values
        /*
        Eigen::VectorXd p_values(t_stats.size());
        for (int i = 0; i < t_stats.size(); ++i)
        {
            double t = std::abs(t_stats(i));
            // This is an approximation of the p-value for the t-distribution
            p_values(i) = (1 - std::erf(t / std::sqrt(2)));
        }
        */

        // Return the maximum p-value
        return t_stats_abs.minCoeff();
    }

    double abs_diffs_r(const std::vector<double> &x, size_t celltype_index)
    {
        const auto &is_a = is_celltype[celltype_index];

        double sum_a = 0.0, sum_not_a = 0.0;
        size_t count_a = 0, count_not_a = 0;
        double min_not_a = std::numeric_limits<double>::max();
        double max_not_a = std::numeric_limits<double>::lowest();

        for (size_t i = 0; i < x.size(); ++i)
        {
            if (is_a[i])
            {
                sum_a += x[i];
                ++count_a;
            }
            else
            {
                sum_not_a += x[i];
                ++count_not_a;
                min_not_a = std::min(min_not_a, x[i]);
                max_not_a = std::max(max_not_a, x[i]);
            }
        }

        if (count_a == 0 || count_not_a == 0)
        {
            return 0.0; // Return 0 if all samples are in one group
        }

        double mean_a = sum_a / count_a;
        double mean_not_a = sum_not_a / count_not_a;

        if (mean_a < mean_not_a)
        {
            return min_not_a - mean_a;
        }
        else
        {
            return mean_a - max_not_a;
        }
    }
};

// [[Rcpp::export]]
List getresult(std::string filename)
{
    SolveM solver(filename);
    solver.compute();
    return List::create(Named("tstats") = wrap(solver.pvals), Named("abs_diffs") = wrap(solver.abs_diffs));
}
