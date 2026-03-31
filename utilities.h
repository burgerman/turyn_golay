//
// Created by Wilfried Wu on 2024-05-29.
//

#ifndef TURYN_GOLAY_UTILITIES_H
#define TURYN_GOLAY_UTILITIES_H

#include <vector>
#include <string>
#include <unordered_map>
#include <climits>
#include <array>

/**
 * @brief Structure to hold precomputed PAF contribution vector.
 * For length m, we care about shifts 1 to m/2.
 */
struct PAFVector {
    std::vector<int> values;
    
    bool operator==(const PAFVector& other) const {
        return values == other.values;
    }
};

/**
 * @brief Hash function for PAFVector to enable use in unordered_map.
 */
struct PAFVectorHash {
    std::size_t operator()(const PAFVector& v) const {
        std::size_t seed = 0;
        for (int val : v.values) {
            seed ^= std::hash<int>{}(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

/**
 * @brief Structure representing a sequence candidate with its precomputed properties.
 */
struct SequenceCandidate {
    std::vector<int> seq;
    int squaredSum;
    PAFVector pafContrib;
};

/**
 * @brief Filters an array of integers by a modulo condition and stores the result.
 * @param arr Input array.
 * @param len Length of the input array.
 * @param mod Modulo value to check against.
 * @param arr_res Output array for filtered elements.
 * @param letter Identifier for the sequence being processed.
 * @return True if filtering was successful.
 */
bool filter_array_element_by_mod(int arr[], int len, int mod, int arr_res[], char letter);

/**
 * @brief Filters a vector of integers by a modulo condition and returns the result.
 * @param sequence Input vector.
 * @param len Expected length of the sequence.
 * @param mod Modulo value to check against.
 * @param letter Identifier for the sequence being processed.
 * @return Filtered vector of integers.
 */
std::vector<int> filter_array_element_by_mod(std::vector<int>& sequence, int len, int mod, char letter);

/**
 * @brief Computes the Periodic Autocorrelation Function (PAF) decomposition for an array.
 * @param s Shift value.
 * @param m Sequence length.
 * @param arr_res Input array.
 * @param letter Identifier for the sequence.
 * @return The computed decomposition value.
 */
int naf_polynomial_decomposition (int s, int m, int arr_res[], char letter);

/**
 * @brief Computes the Periodic Autocorrelation Function (PAF) decomposition for a vector.
 * @param s Shift value.
 * @param m Sequence length.
 * @param sequence Input vector.
 * @return The computed decomposition value.
 */
int naf_polynomial_decomposition (int s, int m, std::vector<int>& sequence);

/**
 * @brief Verifies condition 3 for sequences K and R in Step 2 of the algorithm.
 * @param n Parameter n.
 * @param m Parameter m.
 * @param sequence1 Pointer to the first sequence.
 * @param sequence2 Pointer to the second sequence.
 * @return True if the condition is satisfied.
 */
bool step2_cond3_k_r(int n, int m, int sequence1[], int sequence2[]);

/**
 * @brief Prints a sequence represented by an array.
 * @param arr Input array.
 * @param len Length of the array.
 * @param letter Identifier for the sequence.
 */
void print_sequence(int arr[], int len, char letter);

/**
 * @brief Prints a sequence represented by a vector.
 * @param sequence Input vector.
 * @param len Length of the vector.
 * @param letter Identifier for the sequence.
 */
void print_sequence(std::vector<int>& sequence, int len, char letter);

// New functions moved from main.cpp

/**
 * @brief Helper function for nsoks (n squares of k sums) with memoization.
 * @param n Target sum.
 * @param k Number of squares.
 * @param maxsq Maximum square root value to consider for the next element.
 * @return A vector of all possible combinations of square roots whose squares sum to n.
 */
std::vector<std::vector<int>> nsoksAUX(int n, int k, int maxsq = INT_MAX);

/**
 * @brief Finds all combinations of k squares that sum up to nn.
 * @param nn Target sum.
 * @param kk Number of squares.
 * @return A vector of vectors containing the square roots of the combinations.
 */
std::vector<std::vector<int>> nsoks(int nn, int kk);

/**
 * @brief Computes the Non-periodic Autocorrelation Function (NAF) for an array at a given shift.
 * @param arr Input array.
 * @param n Length of the array.
 * @param s Shift value.
 * @return The computed NAF value.
 */
int npaf_local(const int arr[], int n, int s);

/**
 * @brief Verifies if four sequences (A, B, C, D) satisfy the Turyn-Golay property.
 * @param n Length parameter.
 * @param A Sequence A.
 * @param B Sequence B.
 * @param C Sequence C.
 * @param D Sequence D.
 * @return True if the total NAF sum is 0 for all shifts.
 */
bool sequence_verify(int n, const std::vector<int>& A, const std::vector<int>& B, const std::vector<int>& C, const std::vector<int>& D);

/**
 * @brief Removes duplicate elements from a vector.
 * @param vec Vector to process.
 */
void removeDuplicates(std::vector<int>& vec);

/**
 * @brief Recursively finds all sequences with elements {-1, 0, 1} that sum to a target value.
 * @param m Sequence length.
 * @param targetSum Target sum of elements.
 * @param current_j Current index being processed.
 * @param current_seq Current sequence being built.
 * @param currentSum Current sum of the sequence.
 * @param solutions Vector to store valid solutions.
 * @param lastZero If true, the last element is fixed to 0.
 */
void find_sequences(int m, int targetSum, int current_j, std::vector<int>& current_seq, int currentSum, std::vector<std::vector<int>>& solutions, bool lastZero = false);

/**
 * @brief Verifies condition 3 for K and R sequences using vector input.
 * @param n Parameter n.
 * @param m Parameter m.
 * @param sequence1 Vector sequence 1.
 * @param sequence2 Vector sequence 2.
 * @return True if the condition is satisfied.
 */
bool step2Condition3_k_r(int n, int m, const std::vector<int>& sequence1, const std::vector<int>& sequence2);

/**
 * @brief Verifies condition 3 for P and Q sequences using vector input.
 * @param n Parameter n.
 * @param m Parameter m.
 * @param sequence1 Vector sequence 1.
 * @param sequence2 Vector sequence 2.
 * @return True if the condition is satisfied.
 */
bool step2Condition3_p_q(int n, int m, const std::vector<int>& sequence1, const std::vector<int>& sequence2);

/**
 * @brief Calculates the sum of squares of elements in a sequence.
 * @param sequence Input sequence.
 * @return Sum of squares.
 */
int getSquaredSum(const std::vector<int>& sequence);

/**
 * @brief Verifies condition 5 for a set of four sequences K, R, P, Q.
 * @param m Sequence length.
 * @param sequenceK Sequence K.
 * @param sequenceR Sequence R.
 * @param sequenceP Sequence P.
 * @param sequenceQ Sequence Q.
 * @return True if condition 5 is satisfied for all shifts.
 */
bool step2_condition5(int m, const std::vector<int>& sequenceK, const std::vector<int>& sequenceR,
                       const std::vector<int>& sequenceP, const std::vector<int>& sequenceQ);

/**
 * @brief Orchestrates the finding and verification of solutions for a given quadruple.
 * @param n Parameter n.
 * @param m Parameter m.
 * @param k11 Target sum for sequence K.
 * @param r11 Target sum for sequence R.
 * @param p11 Target sum for sequence P.
 * @param q11 Target sum for sequence Q.
 * @param rank MPI rank of the current process.
 * @param size Total number of MPI processes.
 * @return True if a solution is found and verified.
 */
bool getAllSolutions(int n, int m, int k11, int r11, int p11, int q11, int rank, int size);

/**
 * @brief High-level function to find quadruples of target sums and search for solutions.
 * @param n Parameter n.
 * @param m Parameter m.
 * @param rank MPI rank.
 * @param size MPI size.
 */
void findQuadruple(int n, int m, int rank, int size);

#endif //TURYN_GOLAY_UTILITIES_H
