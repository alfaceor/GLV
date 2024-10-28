#include <gtest/gtest.h>
#include <vector>
#include <stdexcept>
#include <cmath>

// Function under test
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
    // std::cout << numerator << " " << denominator << std::endl;
    
    return numerator / denominator;
}

// Test suite
class BrayCurtisDissimilarityTest : public ::testing::Test {
protected:
    // Utility for comparing floating-point numbers
    static constexpr double tolerance = 1e-6;
};

// Test for equal vectors (expecting 0.0 dissimilarity)
TEST_F(BrayCurtisDissimilarityTest, EqualVectors) {
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> y = {1.0, 2.0, 3.0};
    
    double result = calc_bray_curtis_dissimilarity(x, y);
    EXPECT_NEAR(result, 0.0, tolerance);
}

// Test for different vectors (positive dissimilarity)
TEST_F(BrayCurtisDissimilarityTest, DifferentVectors) {
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> y = {2.0, 3.0, 4.0};
    
    double result = calc_bray_curtis_dissimilarity(x, y);
    // EXPECT_NEAR(result, 0.166667, tolerance);
    EXPECT_NEAR(result, 0.2, tolerance);
}

// Test for vectors with mismatched lengths (expecting an exception)
TEST_F(BrayCurtisDissimilarityTest, MismatchedLengthVectors) {
    std::vector<double> x = {1.0, 2.0};
    std::vector<double> y = {1.0, 2.0, 3.0};
    
    EXPECT_THROW(calc_bray_curtis_dissimilarity(x, y), std::invalid_argument);
}

// Test for vectors where the denominator is zero (expecting an exception)
TEST_F(BrayCurtisDissimilarityTest, ZeroDenominator) {
    std::vector<double> x = {0.0, 0.0, 0.0};
    std::vector<double> y = {0.0, 0.0, 0.0};
    
    EXPECT_THROW(calc_bray_curtis_dissimilarity(x, y), std::invalid_argument);
}

// Test for large vectors with small floating-point differences
TEST_F(BrayCurtisDissimilarityTest, LargeVectorsWithSmallDifferences) {
    std::vector<double> x(1000, 1.0);  // Vector of 1000 elements, all 1.0
    std::vector<double> y(1000, 1.001); // Vector of 1000 elements, all 1.001
    
    double result = calc_bray_curtis_dissimilarity(x, y);
    EXPECT_NEAR(result, 0.0005, tolerance);
}

// Main function for running tests
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

