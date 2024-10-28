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

#include "utils.h"


// Function to calculate the relative abundances
// FIXME: look for the right definition of 
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



int main() {
    
    // TODO: 
    std::string csvsep = ",";

    // Inputs
    std::string fln_species = "species_N_7.csv";
    std::string fln_interactions = "interaction_matrix_N_7.mat";

    // Outputs
    std::string fln_diversity = "df_diversity_N_7.csv";
    std::string fln_trajectories = "df_trajectories_N_7.csv";
    double t_perturbation = 10;


    // Loading data
    // Call the function to read the CSV and store each column independently
    std::map<std::string, std::vector<double>> dataColumns = readCSVWithHeaders(fln_species);
    // std::map<std::string, ColumnData> dataColumns = readCSVWithHeaders02(fln_species);

    // Interaction matrix (eps, NxN matrix)
    std::vector< std::vector<double> > eps = readMatrixFromCSVFile(fln_interactions);

    // Initial populations of the species (Absolute Abundance)
    std::vector<double> AbsAbun = dataColumns["AbsAbun"];
    std::vector<double> AbsAbun_0 = AbsAbun; // Reference to make calculations from initial distribution as reference.

    double dbl_total_AbsAbun = std::accumulate(AbsAbun.begin(), AbsAbun.end(), 0.0);
    std::cout << "dbl_total_AbsAbun: " << dbl_total_AbsAbun << std::endl;

    std::vector<double> RelAbun = std::vector<double>(AbsAbun.size()); //{0.36807, 0.31023, 0.3561, 0.54006, 0.70898, 0.47064, 0.2297, 0.83005, 0.39181, 0.29075, 0.32367};
    for (size_t i =0 ; i<RelAbun.size(); i++) RelAbun[i] = AbsAbun[i]/dbl_total_AbsAbun;

    // Growth rates for each species (alpha)
    std::vector<double> alpha = dataColumns["alpha"];

    std::vector<double> gamma = dataColumns["gamma"];
    for (int i=0; i < AbsAbun.size(); i++){
        std::cout << i << ": " << AbsAbun[i] << csvsep << alpha[i] << csvsep << gamma[i] << std::endl;
    }

    // Input validations:
    // 1. The length of AbsAbun, alpha and all extra columns have the same number of and all bigger than 2
    // 2. Then the interaction matrix has dimensions of num_species x num_species

    int num_species = AbsAbun.size();; // You can change this to any number of species

    // Open the file to be use to write the trajectories for each species time evolution.
    std::ofstream of_trajectories = get_of_trajectories(fln_trajectories, num_species);
    // Check if the file was opened successfully
    if (!of_trajectories.is_open()) {
        return 1;  // Exit if the file couldn't be opened
    }

    // Simulation parameters
    double t = 0.0;       // Initial time
    double t_end = 450.0; // End time
    double dt = 0.01;      // Time step
    int steps = static_cast<int>(t_end / dt);

    std::ofstream of_diversity(fln_diversity);
    of_diversity << "time" << csvsep << "total_AbsAbun" << csvsep << "Shannon" << csvsep << "BrayCurtis" << "\n";
    
    double dbl_shannon = 0.0;
    double dbl_bray_curtis = 0.0;
    // Time evolution
    for (int step = 0; step <= steps; ++step) {
        if (step % 10 == 0) { // Print every 10 steps for readability
            write_of_trajectories(of_trajectories, t, AbsAbun);
            
            // of_trajectories << t << csvsep;
            // for (int i = 0; i < num_species; i++) {
            //     of_trajectories << AbsAbun[i] << csvsep;
            //     dbl_total_AbsAbun += AbsAbun[i];
            // }
            // of_trajectories << "\n";

            dbl_total_AbsAbun = accumulate(AbsAbun.begin(), AbsAbun.end(), 0.0);
            RelAbun = calc_RelAbun(AbsAbun);
            dbl_shannon = calc_ShannonIndex(RelAbun);

            for (int i=0; i < AbsAbun_0.size(); i++){
                std::cout << AbsAbun_0[i] << csvsep;
            }
            std::cout <<  std::endl;

            try {
                dbl_bray_curtis = calc_bray_curtis_dissimilarity(AbsAbun, AbsAbun_0);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Error: " << e.what() << std::endl;
            }

            of_diversity << t << csvsep << dbl_total_AbsAbun<< csvsep << dbl_shannon << csvsep << dbl_bray_curtis << "\n";

        }
        // Calculate Shannon
        if (t> 0 & t < t_perturbation){
            rk4_step_lotka_volterra_perturb_linear(AbsAbun, t, dt, alpha, eps, gamma); // Advance the system by one time step
        }else{
            rk4_step_lotka_volterra(AbsAbun, t, dt, alpha, eps); // Advance the system by one time step
        }
        
        t += dt;
    }

    of_trajectories.close();
    of_diversity.close();

    print_interaction_matrix(eps, num_species);

    return 0;
}
