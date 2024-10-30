#include "utils.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <iomanip>
#include <string>
#include <sstream>
#include <numeric>
#include <stdexcept>


// Function to open the file and write the header
std::ofstream get_of_perturbation(const std::string& filename, const std::string& sep) {
    // Open the file
    std::ofstream file(filename);

    // Check if file is open
    if (!file) {
        std::cerr << "Error: Unable to open file: " << filename << std::endl;
        return file;  // Return the file object even if there's an error, so it can be checked later
    }

    // Write the header
    file << "Time" 
        << sep << "ABfactor"
        << std::endl;
    return file;
}


// Function to open the file and write the header
std::ofstream get_of_trajectories(const std::string& filename, size_t AbsAbun_size, const std::string& sep) {
    // Open the file
    std::ofstream file(filename);

    // Check if file is open
    if (!file) {
        std::cerr << "Error: Unable to open file: " << filename << std::endl;
        return file;  // Return the file object even if there's an error, so it can be checked later
    }

    // Write the header
    file << "Time";
    for (size_t i = 1; i <= AbsAbun_size; ++i) {
        file << sep << "Species." << i;
    }
    file << "\n";

    return file;
}

// Function to write time and absolute abundances (AbsAbun) to the file
void write_of_trajectories(std::ofstream& file, double time, const std::vector<double>& AbsAbun, const std::string& sep) {
    // Write time first
    file << time;

    // Write each element of AbsAbun separated by the defined separator
    for (const auto& abundance : AbsAbun) {
        file << sep << abundance;
    }
    file << "\n";
}


// Print the interaction matrix 
// TODO: GENERALIZE TO REDIRECT TO A FILE
void print_interaction_matrix(const std::vector<std::vector<double>>& eps, int num_species, const std::string& sep){
    for (int i = 0; i < num_species; i++) {
        for (int j = 0; j < num_species-1; j++) {
            std::cout<< eps[i][j] << sep;
            // std::cout<< dAbsAbun_dt[i] << std::endl;
        }
        std::cout << eps[i][num_species-1] << std::endl;
    }
}


// // Print each column to verify the result
    // TODO: MAKE AN OUTPUT FUNCTION
    // for (const auto& [header, column] : dataColumns) {
    //     std::cout << header << ": ";
    //     for (double value : column) {
    //         std::cout << value << " ";
    //     }
    //     std::cout << std::endl;
    // }

double calc_TotalAbun(const std::vector<double>& absoluteAbundances){
    double totalAbundance = std::accumulate(absoluteAbundances.begin(), absoluteAbundances.end(), 0.0);
    return totalAbundance;
}

// Function to calculate the relative abundances
std::vector<double> calc_RelAbun(const std::vector<double>& absoluteAbundances) {
    // Sum of all absolute abundances
    double totalAbundance = calc_TotalAbun(absoluteAbundances);
    // double totalAbundance = std::accumulate(absoluteAbundances.begin(), absoluteAbundances.end(), 0.0);
    
    // Vector to store relative abundances
    std::vector<double> relativeAbundances(absoluteAbundances.size());
    
    // Calculate relative abundance for each species
    for (size_t i = 0; i < absoluteAbundances.size(); ++i) {
        relativeAbundances[i] = absoluteAbundances[i] / totalAbundance;
    }

    return relativeAbundances;
}

// Function to calculate Shannon index based on the relative abundances
double calc_ShannonIndex(const std::vector<double>& relativeAbundances) {
    double shannonIndex = 0.0;
    
    for (double relAbun : relativeAbundances) {
        if (relAbun > 0) { // To avoid log(0) and NaNs
            shannonIndex += -relAbun * std::log2(relAbun);
        }
    }
    
    return shannonIndex;
}

// Function to read a CSV file with headers and save each column independently
// FIXME: ADD csvsep as default ',' 
std::map<std::string, std::vector<double>> readCSVWithHeaders01(const std::string& fileName) {
    std::ifstream file(fileName);
    std::map<std::string, std::vector<double>> dataColumns;
    std::vector<std::string> headers;

    // Check if the file was opened successfully
    if (!file) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return dataColumns; // Return empty map if file can't be opened
    }

    std::string line;

    // Read the header line
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string header;

        // Read each header, store it in headers list, and initialize an empty vector in dataColumns
        while (std::getline(ss, header, ',')) {
            headers.push_back(header);
            dataColumns[header] = std::vector<double>();  // Create an empty vector for each column
        }
    }

    // Read the remaining lines (the data)
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        
        // Iterate using the headers list to ensure column order
        for (const auto& header : headers) {
            if (std::getline(ss, value, ',')) {
                dataColumns[header].push_back(std::stod(value));  // Convert string to double and add to column
            }
        }
    }

    // Close the file
    file.close();

    return dataColumns;
}



// Helper function to remove leading and trailing quotes
std::string removeQuotes(const std::string& str) {
    if (str.size() > 1 && (str.front() == '"' || str.front() == '\'') && str.front() == str.back()) {
        return str.substr(1, str.size() - 2);  // Remove both first and last character if they are matching quotes
    }
    return str;
}

std::map<std::string, std::vector<double>> readCSVWithHeaders(const std::string& fileName) {
    std::ifstream file(fileName);
    std::map<std::string, std::vector<double>> dataColumns;
    std::vector<std::string> headers;

    // Check if the file was opened successfully
    if (!file) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return dataColumns; // Return empty map if file can't be opened
    }

    std::string line;

    // Read the header line
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string header;

        // Read each header, remove any surrounding quotes, and initialize an empty vector in dataColumns
        while (std::getline(ss, header, ',')) {
            std::string cleanedHeader = removeQuotes(header);
            headers.push_back(cleanedHeader);
            dataColumns[cleanedHeader] = std::vector<double>();  // Create an empty vector for each column
        }
    }

    // Read the remaining lines (the data)
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        
        // Iterate using the headers list to ensure column order
        for (const auto& header : headers) {
            if (std::getline(ss, value, ',')) {
                try {
                    dataColumns[header].push_back(std::stod(value));  // Convert string to double and add to column
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Non-numeric data found in column '" << header << "': " << value << std::endl;
                    dataColumns[header].push_back(0.0);  // Optional: Handle non-numeric data
                }
            }
        }
    }

    // Close the file
    file.close();

    return dataColumns;
}




// FIXME: ADD csvsep as default ',' 
std::map<std::string, ColumnData> readCSVWithHeaders02(const std::string& fileName) {
    std::ifstream file(fileName);
    std::map<std::string, ColumnData> dataColumns;

    // Check if the file was opened successfully
    if (!file) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return dataColumns; // Return empty map if file can't be opened
    }

    std::string line;

    // Read the header line
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string header;

        // Read each header and initialize an empty vector for it
        while (std::getline(ss, header, ',')) {
            // Assume all columns are strings at first
            dataColumns[header] = std::vector<std::string>();
        }
    }

    // Read the remaining lines (the data)
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        // Iterate over the map to access headers and fill columns
        auto it = dataColumns.begin();
        while (std::getline(ss, value, ',')) {
            // Check if this column already contains double data
            if (std::holds_alternative<std::vector<double>>(it->second)) {
                try {
                    // Try converting the value to double and adding to the column
                    std::get<std::vector<double>>(it->second).push_back(std::stod(value));
                } catch (const std::invalid_argument& e) {
                    // If conversion fails, print a message and exit
                    std::cerr << "Error: Non-numeric value '" << value << "' in numeric column '" << it->first << "'\n";
                    return dataColumns;
                }
            } else {
                // Add string values to the string vector
                std::get<std::vector<std::string>>(it->second).push_back(value);
            }

            ++it; // Move to the next column
            if (it == dataColumns.end()) {
                break;
            }
        }
    }

    // Close the file
    file.close();

    return dataColumns;
}

// Function to read a CSV file with headers and save each column independently with string as columns too.
// Function to read a vector of vectors of doubles from a CSV file
std::vector<std::vector<double>> readMatrixFromCSVFile(const std::string& fileName) {
    std::ifstream file(fileName);
    std::vector<std::vector<double>> vecOfVecs;

    // Check if the file was opened successfully
    if (!file) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return vecOfVecs; // Return empty vector if file can't be opened
    }

    std::string line;
    // Read file line by line
    while (std::getline(file, line)) {
        std::vector<double> innerVec;  // Vector to store a row
        std::stringstream ss(line);
        std::string value;

        // Extract numbers from the line, split by commas
        while (std::getline(ss, value, ',')) {
            innerVec.push_back(std::stod(value));  // Convert string to double
        }

        // Add the inner vector (row) to the outer vector (2D matrix)
        vecOfVecs.push_back(innerVec);
    }

    // Close the file
    file.close();

    return vecOfVecs;
}


// Function to calculate the Lotka-Volterra equations for any number of species
void lotka_volterra(const std::vector<double>& AbsAbun, std::vector<double>& dAbsAbun_dt, 
                    const std::vector<double>& alpha, const std::vector<std::vector<double>>& eps) {
    int num_species = AbsAbun.size();
    
    for (int i = 0; i < num_species; i++) {
        dAbsAbun_dt[i] = alpha[i] * AbsAbun[i]; // logistic growth
        for (int j = 0; j < num_species; j++) {
            dAbsAbun_dt[i] += eps[i][j] * AbsAbun[i] * AbsAbun[j]; // Interaction term
            // std::cout<< dAbsAbun_dt[i] << std::endl;
        }
    }
}


// Runge-Kutta 4th order method for solving ODEs
void rk4_step_lotka_volterra(std::vector<double>& AbsAbun, double t, double dt, 
            const std::vector<double>& alpha,
            const std::vector<std::vector<double>>& eps) {
    int n = AbsAbun.size();
    std::vector<double> k1(n), k2(n), k3(n), k4(n), AbsAbun_temp(n);

    // k1 = f(AbsAbun, t)
    lotka_volterra(AbsAbun, k1, alpha, eps);

    // k2 = f(AbsAbun + dt/2 * k1, t + dt/2)
    for (int i = 0; i < n; i++) AbsAbun_temp[i] = AbsAbun[i] + dt / 2.0 * k1[i];
    lotka_volterra(AbsAbun_temp, k2, alpha, eps);

    // k3 = f(AbsAbun + dt/2 * k2, t + dt/2)
    for (int i = 0; i < n; i++) AbsAbun_temp[i] = AbsAbun[i] + dt / 2.0 * k2[i];
    lotka_volterra(AbsAbun_temp, k3, alpha, eps);

    // k4 = f(AbsAbun + dt * k3, t + dt)
    for (int i = 0; i < n; i++) AbsAbun_temp[i] = AbsAbun[i] + dt * k3[i];
    lotka_volterra(AbsAbun_temp, k4, alpha, eps);

    // Update AbsAbun using the RK4 formula
    for (int i = 0; i < n; i++) {
        AbsAbun[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}



// Function to calculate the Lotka-Volterra equations for any number of species with a linear perturbation proportional to the species abundance
void lotka_volterra_perturb_linear(const std::vector<double>& AbsAbun, std::vector<double>& dAbsAbun_dt, 
                    const std::vector<double>& alpha, 
                    const std::vector<std::vector<double>>& eps, 
                    const std::vector<double>& gamma) {
    int num_species = AbsAbun.size();
    
    for (int i = 0; i < num_species; i++) {
        dAbsAbun_dt[i] = alpha[i] * AbsAbun[i] - gamma[i] * AbsAbun[i]; // logistic growth
        for (int j = 0; j < num_species; j++) {
            dAbsAbun_dt[i] += eps[i][j] * AbsAbun[i] * AbsAbun[j]; // Interaction term
            // std::cout<< dAbsAbun_dt[i] << std::endl;
        }
    }
}


// Runge-Kutta 4th order method for solving ODEs for linear perturbation
void rk4_step_lotka_volterra_perturb_linear(std::vector<double>& AbsAbun, double t, double dt, 
            const std::vector<double>& alpha,
            const std::vector<std::vector<double>>& eps,
            const std::vector<double>& gamma) {
    int n = AbsAbun.size();
    std::vector<double> k1(n), k2(n), k3(n), k4(n), AbsAbun_temp(n);

    // k1 = f(AbsAbun, t)
    lotka_volterra_perturb_linear(AbsAbun, k1, alpha, eps, gamma);

    // k2 = f(AbsAbun + dt/2 * k1, t + dt/2)
    for (int i = 0; i < n; i++) AbsAbun_temp[i] = AbsAbun[i] + dt / 2.0 * k1[i];
    lotka_volterra_perturb_linear(AbsAbun_temp, k2, alpha, eps, gamma);

    // k3 = f(AbsAbun + dt/2 * k2, t + dt/2)
    for (int i = 0; i < n; i++) AbsAbun_temp[i] = AbsAbun[i] + dt / 2.0 * k2[i];
    lotka_volterra_perturb_linear(AbsAbun_temp, k3, alpha, eps, gamma);

    // k4 = f(AbsAbun + dt * k3, t + dt)
    for (int i = 0; i < n; i++) AbsAbun_temp[i] = AbsAbun[i] + dt * k3[i];
    lotka_volterra_perturb_linear(AbsAbun_temp, k4, alpha, eps, gamma);

    // Update AbsAbun using the RK4 formula
    for (int i = 0; i < n; i++) {
        AbsAbun[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}


// Function to calculate the relative abundances
double calc_bray_curtis_dissimilarity(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Vectors must be of the same length.");
    }
    
    double numerator = 0.0;
    double denominator = 0.0;
    
    for (size_t i = 0; i < x.size(); ++i) {
        numerator += std::abs(x[i] - y[i]);
        denominator += x[i] + y[i];
    }
    
    if (denominator == 0.0) {
        throw std::invalid_argument("Denominator is zero, cannot calculate Bray-Curtis dissimilarity.");
    }
    
    return numerator / denominator;
}

// Function to calculate the antibiotic concentration
// FIXME: Generalize parameters
double antibiotic_concentration(double t){
    double ka = 0.99;
    double ke = 0.21;
    double Tend = 5;
    double Ct = 0.0;

    if(t < 0){
        Ct = 0.0;
    }if (t < Tend){
        Ct = (1.0/(1-exp(-ka*Tend)))*(1-exp(-ka*t));
    }else{
        Ct = exp(-ke*(t-Tend));
    }

    if (Ct < 1e-10) Ct = 0.0;
    return(Ct);
}


// Function to calculate the Lotka-Volterra equations for any number of species with a linear perturbation proportional to the species abundance
void lotka_volterra_Ab_pertur_001(const std::vector<double>& AbsAbun, std::vector<double>& dAbsAbun_dt, 
            const std::vector<double>& alpha, 
            const std::vector<std::vector<double>>& eps, 
            const std::vector<double>& gamma, 
            double t) {
    int num_species = AbsAbun.size();
    
    for (int i = 0; i < num_species; i++) {
        dAbsAbun_dt[i] = alpha[i] * AbsAbun[i] - gamma[i] * antibiotic_concentration(t) * AbsAbun[i]; // logistic growth
        for (int j = 0; j < num_species; j++) {
            dAbsAbun_dt[i] += eps[i][j] * AbsAbun[i] * AbsAbun[j]; // Interaction term
            // std::cout<< dAbsAbun_dt[i] << std::endl;
        }
    }
}


// Runge-Kutta 4th order method for solving ODEs for linear perturbation
void rk4_step_lotka_volterra_Ab_pertur_001(std::vector<double>& AbsAbun, double t, double dt, 
            const std::vector<double>& alpha,
            const std::vector<std::vector<double>>& eps,
            const std::vector<double>& gamma) {
    int n = AbsAbun.size();
    std::vector<double> k1(n), k2(n), k3(n), k4(n), AbsAbun_temp(n);

    // k1 = f(AbsAbun, t)
    lotka_volterra_Ab_pertur_001(AbsAbun, k1, alpha, eps, gamma, t);

    // k2 = f(AbsAbun + dt/2 * k1, t + dt/2)
    for (int i = 0; i < n; i++) AbsAbun_temp[i] = AbsAbun[i] + dt / 2.0 * k1[i];
    lotka_volterra_Ab_pertur_001(AbsAbun_temp, k2, alpha, eps, gamma, t);

    // k3 = f(AbsAbun + dt/2 * k2, t + dt/2)
    for (int i = 0; i < n; i++) AbsAbun_temp[i] = AbsAbun[i] + dt / 2.0 * k2[i];
    lotka_volterra_Ab_pertur_001(AbsAbun_temp, k3, alpha, eps, gamma, t);

    // k4 = f(AbsAbun + dt * k3, t + dt)
    for (int i = 0; i < n; i++) AbsAbun_temp[i] = AbsAbun[i] + dt * k3[i];
    lotka_volterra_Ab_pertur_001(AbsAbun_temp, k4, alpha, eps, gamma, t);

    // Update AbsAbun using the RK4 formula
    for (int i = 0; i < n; i++) {
        AbsAbun[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}



// // Small constant to avoid log(0)
// const double EPSILON = 1e-12;

// Kullback-Leibler divergence with log base 2
double calc_KLDivergence(const std::vector<double>& P, const std::vector<double>& Q) {
    double klDiv = 0.0;
    for (size_t i = 0; i < P.size(); ++i) {
        if (P[i] > 0) {
            klDiv += P[i] * std::log2(P[i] / (Q[i] + EPSILON) + EPSILON);
        }
    }
    return klDiv;
}

// Jensen-Shannon Divergence
double calc_JensenShannonDivergence(const std::vector<double>& P, const std::vector<double>& Q) {
    if (P.size() != Q.size()) {
        throw std::invalid_argument("Vectors P and Q must be of the same length.");
    }
    
    // Calculate the average distribution M = 0.5 * (P + Q)
    std::vector<double> M(P.size());
    for (size_t i = 0; i < P.size(); ++i) {
        M[i] = 0.5 * (P[i] + Q[i]);
    }
    
    // Calculate JS Divergence: 0.5 * (KL(P || M) + KL(Q || M))
    return 0.5 * (calc_KLDivergence(P, M) + calc_KLDivergence(Q, M));
}




// Function to calculate the Lotka-Volterra equations for any number of species with a linear perturbation proportional to the species abundance
void lotka_volterra_w_perturbation(const std::vector<double>& AbsAbun, std::vector<double>& dAbsAbun_dt, 
            const std::vector<double>& alpha, 
            const std::vector<std::vector<double>>& eps, 
            const std::vector<double>& gamma, 
            double perturb_factor) {
    int num_species = AbsAbun.size();
    
    for (int i = 0; i < num_species; i++) {
        dAbsAbun_dt[i] = alpha[i] * AbsAbun[i] - gamma[i] * perturb_factor * AbsAbun[i]; // logistic growth
        for (int j = 0; j < num_species; j++) {
            dAbsAbun_dt[i] += eps[i][j] * AbsAbun[i] * AbsAbun[j]; // Interaction term
            // std::cout<< dAbsAbun_dt[i] << std::endl;
        }
    }
}


// Runge-Kutta 4th order method for solving ODEs for linear perturbation
void rk4_step_lotka_volterra_w_perturbation(std::vector<double>& AbsAbun, double t, double dt, 
            const std::vector<double>& alpha,
            const std::vector<std::vector<double>>& eps,
            const std::vector<double>& gamma,
            double pertub_factor) {
    int n = AbsAbun.size();
    std::vector<double> k1(n), k2(n), k3(n), k4(n), AbsAbun_temp(n);

    // k1 = f(AbsAbun, t)
    lotka_volterra_w_perturbation(AbsAbun, k1, alpha, eps, gamma, pertub_factor);

    // k2 = f(AbsAbun + dt/2 * k1, t + dt/2)
    for (int i = 0; i < n; i++) AbsAbun_temp[i] = AbsAbun[i] + dt / 2.0 * k1[i];
    lotka_volterra_w_perturbation(AbsAbun_temp, k2, alpha, eps, gamma, pertub_factor);

    // k3 = f(AbsAbun + dt/2 * k2, t + dt/2)
    for (int i = 0; i < n; i++) AbsAbun_temp[i] = AbsAbun[i] + dt / 2.0 * k2[i];
    lotka_volterra_w_perturbation(AbsAbun_temp, k3, alpha, eps, gamma, pertub_factor);

    // k4 = f(AbsAbun + dt * k3, t + dt)
    for (int i = 0; i < n; i++) AbsAbun_temp[i] = AbsAbun[i] + dt * k3[i];
    lotka_volterra_w_perturbation(AbsAbun_temp, k4, alpha, eps, gamma, pertub_factor);

    // Update AbsAbun using the RK4 formula
    for (int i = 0; i < n; i++) {
        AbsAbun[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}






