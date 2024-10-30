#ifndef UTILS_H
#define UTILS_H

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <numeric>
#include <cmath>
#include <variant> // For std::variant
#include <type_traits> // For std::holds_alternative

#include <stdexcept>


// Small constant to avoid log(0)
const double EPSILON = 1e-12;


// Define a variant that can hold either double or string values in the vectors
using ColumnData = std::variant<std::vector<double>, std::vector<std::string>>;

// Function to open the file and write the header
std::ofstream get_of_perturbation(const std::string& filename, const std::string& sep = ",");

// Function to open the file and write the header
std::ofstream get_of_trajectories(const std::string& filename, size_t AbsAbun_size, const std::string& sep = ",");

// Function to write time and absolute abundances (AbsAbun) to the file
void write_of_trajectories(std::ofstream& file, double time, const std::vector<double>& AbsAbun, const std::string& sep = ",");

// Function to print the interaction matrix
void print_interaction_matrix(const std::vector<std::vector<double>>& eps, int num_species, const std::string& sep = ",");

// Function to calculate the total absolute abundance
double calc_TotalAbun(const std::vector<double>& absoluteAbundances);

// Function to calculate the relative abundances
std::vector<double> calc_RelAbun(const std::vector<double>& absoluteAbundances);

// Function to calculate Shannon index based on the relative abundances
double calc_ShannonIndex(const std::vector<double>& relativeAbundances);

// Helper function to remove leading and trailing quotes
std::string removeQuotes(const std::string& str);


std::map<std::string, std::vector<double>> readCSVWithHeaders(const std::string& fileName);

// Function to read a CSV file with headers and save each column independently
std::map<std::string, std::vector<double>> readCSVWithHeaders01(const std::string& fileName);

// Function to read a CSV file with headers and save each column independently
std::map<std::string, ColumnData> readCSVWithHeaders02(const std::string& fileName);

// Function to read a vector of vectors of doubles from a CSV file
std::vector<std::vector<double>> readMatrixFromCSVFile(const std::string& fileName);

// Function to calculate the Lotka-Volterra equations for any number of species
void lotka_volterra(const std::vector<double>& AbsAbun, std::vector<double>& dAbsAbun_dt, 
                    const std::vector<double>& alpha, 
                    const std::vector<std::vector<double>>& eps);

// Runge-Kutta 4th order method for solving ODEs
void rk4_step_lotka_volterra(std::vector<double>& AbsAbun, double t, double dt, 
                              const std::vector<double>& alpha, 
                              const std::vector<std::vector<double>>& eps);

// Function to calculate the Lotka-Volterra equations for any number of species with a linear perturbation proportional to the species abundance
void lotka_volterra_perturb_linear(const std::vector<double>& AbsAbun, std::vector<double>& dAbsAbun_dt, 
            const std::vector<double>& alpha, 
            const std::vector<std::vector<double>>& eps, 
            const std::vector<double>& gamma);

// Runge-Kutta 4th order method for solving ODEs for linear perturbation
void rk4_step_lotka_volterra_perturb_linear(std::vector<double>& AbsAbun, double t, double dt, 
            const std::vector<double>& alpha,
            const std::vector<std::vector<double>>& eps,
            const std::vector<double>& gamma);


// Function to calculate the relative abundances
double calc_bray_curtis_dissimilarity(const std::vector<double>& x, const std::vector<double>& y);

// Function to calculate the antibiotic concentration
double antibiotic_concentration(double t);

// Function to calculate the Lotka-Volterra equations for any number of species with a linear perturbation proportional to the species abundance
void lotka_volterra_Ab_pertur_001(const std::vector<double>& AbsAbun, std::vector<double>& dAbsAbun_dt, 
            const std::vector<double>& alpha, 
            const std::vector<std::vector<double>>& eps, 
            const std::vector<double>& gamma, 
            double t);


// Runge-Kutta 4th order method for solving ODEs for linear perturbation
void rk4_step_lotka_volterra_Ab_pertur_001(std::vector<double>& AbsAbun, double t, double dt, 
            const std::vector<double>& alpha,
            const std::vector<std::vector<double>>& eps,
            const std::vector<double>& gamma);


// Kullback-Leibler divergence with log base 2
double calc_KLDivergence(const std::vector<double>& P, const std::vector<double>& Q);


// Jensen-Shannon Divergence
double calc_JensenShannonDivergence(const std::vector<double>& P, const std::vector<double>& Q);



// Function to calculate the Lotka-Volterra equations for any number of species with a linear perturbation proportional to the species abundance
void lotka_volterra_w_perturbation(const std::vector<double>& AbsAbun, std::vector<double>& dAbsAbun_dt, 
            const std::vector<double>& alpha, 
            const std::vector<std::vector<double>>& eps, 
            const std::vector<double>& gamma, 
            double perturb_factor);


// Runge-Kutta 4th order method for solving ODEs for linear perturbation
void rk4_step_lotka_volterra_w_perturbation(std::vector<double>& AbsAbun, double t, double dt, 
            const std::vector<double>& alpha,
            const std::vector<std::vector<double>>& eps,
            const std::vector<double>& gamma,
            double pertub_factor);



#endif // UTILS_H
