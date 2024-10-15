#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <iomanip>
#include <string>
#include <sstream>
#include <numeric>




// Function to open the file and write the header
std::ofstream get_of_trajectories(const std::string& filename, size_t AbsAbun_size, const std::string& sep = ",") {
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
void write_of_trajectories(std::ofstream& file, double time, const std::vector<double>& AbsAbun, const std::string& sep = ",") {
    // Write time first
    file << time;

    // Write each element of AbsAbun separated by the defined separator
    for (const auto& abundance : AbsAbun) {
        file << sep << abundance;
    }
    file << "\n";
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
std::map<std::string, std::vector<double>> readCSVWithHeaders(const std::string& fileName) {
    std::ifstream file(fileName);
    std::map<std::string, std::vector<double>> dataColumns;

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
            dataColumns[header] = std::vector<double>();  // Create an empty vector for each column
        }
    }

    // Read the remaining lines (the data)
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        size_t columnIndex = 0;

        // Iterate over the map to access headers and fill columns
        for (auto& [header, column] : dataColumns) {
            if (std::getline(ss, value, ',')) {
                column.push_back(std::stod(value));  // Convert string to double and add to column
            }
        }
    }

    // Close the file
    file.close();

    return dataColumns;
}


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
void rk4_step(std::vector<double>& AbsAbun, double t, double dt, const std::vector<double>& alpha,
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

int main() {
    
    // TODO: 
    std::string csvsep = ",";

    // Inputs
    std::string fln_species = "species_N_7.csv";
    std::string fln_interactions = "interaction_matrix_N_7.mat";

    // Outputs
    std::string fln_diversity = "df_diversity_N_7.csv";
    std::string fln_trajectories = "df_trajectories_N_7.csv";


    // Loading data
    // Call the function to read the CSV and store each column independently
    std::map<std::string, std::vector<double>> dataColumns = readCSVWithHeaders(fln_species);

    // Interaction matrix (eps, NxN matrix)
    std::vector< std::vector<double> > eps = readMatrixFromCSVFile(fln_interactions);

    // Initial populations of the species (Absolute Abundance)
    std::vector<double> AbsAbun = dataColumns["AbsAbun"];

    double dbl_total_AbsAbun = std::accumulate(AbsAbun.begin(), AbsAbun.end(), 0);
    std::vector<double> RelAbun = std::vector<double>(AbsAbun.size()); //{0.36807, 0.31023, 0.3561, 0.54006, 0.70898, 0.47064, 0.2297, 0.83005, 0.39181, 0.29075, 0.32367};
    for (size_t i =0 ; i<RelAbun.size(); i++) RelAbun[i] = AbsAbun[i]/dbl_total_AbsAbun;

    // Growth rates for each species (alpha)
    std::vector<double> alpha = dataColumns["alpha"];

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

    for (int i = 0; i < AbsAbun.size(); i++) {         // Loop through rows
        for (int j = 0; j < AbsAbun.size(); j++) {     // Loop through columns
            std::cout << eps[i][j] << " ";  // Print each element
        }
        std::cout << std::endl;  // New line after each row
    }



    // Simulation parameters
    double t = 0.0;       // Initial time
    double t_end = 45.0; // End time
    double dt = 0.01;      // Time step
    int steps = static_cast<int>(t_end / dt);

    std::ofstream of_diversity(fln_diversity);
    of_diversity << "time" << csvsep << "total_AbsAbun" << csvsep << "Shannon" << "\n";
    // ofile.open("DF_DIVERSITY.CSV");
    // Output formatting

    // std::cout << std::fixed << std::setprecision(2);
    // std::cout << "Time\t";
    // for (int i = 0; i < num_species; i++) std::cout << "Species " << i + 1 << "\t";
    // std::cout << "\n";

    
    double dbl_shannon = 0.0;
    // Time evolution
    for (int step = 0; step <= steps; ++step) {
        // if (step % 10 == 0) { // Print every 10 steps for readability
            write_of_trajectories(of_trajectories, t, AbsAbun);
            
            // of_trajectories << t << csvsep;
            // for (int i = 0; i < num_species; i++) {
            //     of_trajectories << AbsAbun[i] << csvsep;
            //     dbl_total_AbsAbun += AbsAbun[i];
            // }
            // of_trajectories << "\n";

            dbl_total_AbsAbun = accumulate(AbsAbun.begin(), AbsAbun.end(), 0);
            RelAbun = calc_RelAbun(AbsAbun);
            dbl_shannon = calc_ShannonIndex(RelAbun);

            of_diversity << t << csvsep << dbl_total_AbsAbun<< csvsep << dbl_shannon << "\n";

        // }
        // Calculate Shannon

        rk4_step(AbsAbun, t, dt, alpha, eps); // Advance the system by one time step
        t += dt;
    }

    of_trajectories.close();
    of_diversity.close();

    return 0;
}
