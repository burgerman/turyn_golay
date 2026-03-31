//
// Created by Wilfried Wu on 2024-05-29.
//

#include <cstdio>
#include <vector>
#include <iostream>
#include <cmath>
#include <numeric>
#include <unordered_set>
#include <sstream>
#include <unordered_map>
#include <set>
#include <tuple>
#include <climits>
#include <algorithm>
#include <mpi.h>
#include "utilities.h"

// Threshold for total combinations to prevent excessive computation on a single node
const long long THRESHOLD_V = 1000000000000LL;

// Structure to serve as a key for memoization in nsoksAUX
struct MemoKey {
    int n, k, maxsq;
    bool operator==(const MemoKey& other) const {
        return n == other.n && k == other.k && maxsq == other.maxsq;
    }
};

// Hash function for MemoKey to be used in unordered_map
struct MemoKeyHash {
    std::size_t operator()(const MemoKey& k) const {
        return (std::size_t(k.n) << 16) ^ (std::size_t(k.k) << 8) ^ std::size_t(k.maxsq);
    }
};

// Global memoization table for nsoks calculations
std::unordered_map<MemoKey, std::vector<std::vector<int>>, MemoKeyHash> memo;

/**
 * Computes the Periodic Autocorrelation Function (PAF) contribution for a sequence at shift s.
 */
int naf_polynomial_decomposition (int s, int m, std::vector<int>& sequence) {
    int i;
    int res = 0;
    for(i=1;i<=m-s;i++) {
        res += sequence[i-1]*sequence[i-1+s];
    }
    return res;
}

/**
 * Prints the sequence in a readable format: Sequence [letter]: [e1, e2, ..., en]
 */
void print_sequence(std::vector<int>& sequence, int len, char letter) {
    int i;
    printf("Sequence %c: [", letter);
    for (i=0; i<len;i++) {
        if(i!= len-1) {
            printf("%d, ", sequence[i]);
        } else {
            printf("%d]\n", sequence[i]);
        }
    }
}

/**
 * Recursive function with memoization to find all sets of k integers whose squares sum to n.
 * Enforces non-increasing order using maxsq to avoid permutations.
 */
std::vector<std::vector<int>> nsoksAUX(int n, int k, int maxsq) {
    MemoKey key = {n, k, maxsq};
    auto it = memo.find(key);
    if (it != memo.end()) return it->second;

    std::vector<std::vector<int>> seqs;
    if (k == 1) {
        int root = static_cast<int>(std::sqrt(n));
        if (root * root == n && root <= maxsq) seqs.push_back({root});
        return memo[key] = seqs;
    }

    // Pruning: minTs ensures we can reach the target sum with remaining k elements
    int maxTs = std::min(static_cast<int>(std::sqrt(n - k + 1)), maxsq);
    int minTs = static_cast<int>(std::sqrt(n / k));
    if (minTs * minTs * k < n) minTs++;

    for (int trialSquare = minTs; trialSquare <= maxTs; ++trialSquare) {
        std::vector<std::vector<int>> val = nsoksAUX(n - trialSquare * trialSquare, k - 1, trialSquare);
        for (auto& subSeq : val) {
            std::vector<int> seq = subSeq;
            seq.push_back(trialSquare);
            seqs.push_back(std::move(seq));
        }
    }
    return memo[key] = seqs;
}

/**
 * Wrapper for nsoksAUX that handles cases where fewer than kk squares are used by padding with zeros.
 */
std::vector<std::vector<int>> nsoks(int nn, int kk) {
    std::vector<std::vector<int>> res;
    for (int kkk = 1; kkk <= kk; ++kkk) {
        std::vector<std::vector<int>> aux_res = nsoksAUX(nn, kkk);
        for (auto& seq : aux_res) {
            while (seq.size() < static_cast<std::size_t>(kk)) seq.insert(seq.begin(), 0);
            res.push_back(std::move(seq));
        }
    }
    return res;
}

/**
 * Local NPAF calculation for a raw array.
 */
int npaf_local(const int arr[], int n, int s) {
    int sum = 0;
    for (int i = 0; i < n - s; i++) sum += arr[i] * arr[i + s];
    return sum;
}

/**
 * Verifies if the sum of NPAF of four sequences is zero for all shifts s from 1 to n.
 * This is the defining property for the sequences being searched.
 */
bool sequence_verify(int n, const std::vector<int>& A, const std::vector<int>& B, const std::vector<int>& C, const std::vector<int>& D) {
    for (int s = 1; s <= n; ++s) {
        int sum = npaf_local(A.data(), n + 1, s) + 
                  npaf_local(B.data(), n + 1, s) + 
                  npaf_local(C.data(), n, s) + 
                  npaf_local(D.data(), n, s);
        if (sum != 0) return false;
    }
    return true;
}

/**
 * Utility to sort and remove duplicates from a vector.
 */
void removeDuplicates(std::vector<int>& vec) {
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

/**
 * Generates all sequences of length m with specified target sum using recursion.
 * elements are restricted to {-1, 1} unless lastZero is true, where the last element is 0.
 */
void find_sequences(int m, int targetSum, int current_j, std::vector<int>& current_seq, int currentSum, std::vector<std::vector<int>>& solutions, bool lastZero) {
    if (current_j > m) {
        if (currentSum == targetSum) solutions.push_back(current_seq);
        return;
    }
    if (lastZero && current_j == m) {
        current_seq[current_j - 1] = 0;
        find_sequences(m, targetSum, current_j + 1, current_seq, currentSum, solutions, lastZero);
    } else {
        current_seq[current_j - 1] = 1;
        find_sequences(m, targetSum, current_j + 1, current_seq, currentSum + 1, solutions, lastZero);
        current_seq[current_j - 1] = -1;
        find_sequences(m, targetSum, current_j + 1, current_seq, currentSum - 1, solutions, lastZero);
    }
}

/**
 * Verifies specific modulo and sum constraints for K and R sequences as per Turyn's construction.
 */
bool step2Condition3_k_r(int n, int m, const std::vector<int>& sequence1, const std::vector<int>& sequence2) {
    bool isMultipleOfM = (n % m == 0);
    for (int i = 1; i <= m; ++i) {
        int maxV = (n + 1 - i) / m + 1;
        if (std::abs(sequence1[i - 1]) <= maxV && std::abs(sequence2[i - 1]) <= maxV) {
            if (maxV % 2 == 0) {
                if ((sequence1[i - 1] % 2 != 0) || (sequence2[i - 1] % 2 != 0)) return false;
            } else {
                if ((sequence1[i - 1] % 2 == 0) || (sequence2[i - 1] % 2 == 0)) return false;
            }
        }
    }
    for (int j = 2; j <= m; ++j) {
        int index = n + 2 - j;
        int sum = sequence1[j - 1] + sequence2[j - 1] + sequence1[index - 1] + sequence2[index - 1];
        if (j == (n + 1) && index == 1) {
            if (isMultipleOfM) { if ((sum % 4) != 0) return false; }
            else { if ((sum % 4) != 2) return false; }
        } else if ((sum % 4) != 0) return false;
    }
    return true;
}

/**
 * Verifies specific modulo and sum constraints for P and Q sequences.
 */
bool step2Condition3_p_q(int n, int m, const std::vector<int>& sequence1, const std::vector<int>& sequence2) {
    for (int i = 1; i < m; ++i) {
        int maxV = (n - i) / m + 1;
        if (std::abs(sequence1[i - 1]) <= maxV && std::abs(sequence2[i - 1]) <= maxV) {
            if (maxV % 2 == 0) {
                if ((sequence1[i - 1] % 2 != 0) || (sequence2[i - 1] % 2 != 0)) return false;
            } else {
                if ((sequence1[i - 1] % 2 == 0) || (sequence2[i - 1] % 2 == 0)) return false;
            }
        }
    }
    for (int j = 1; j < m; ++j) {
        int index = n + 1 - j;
        int sum1 = sequence1[j - 1] + sequence2[j - 1] + sequence1[index - 1] + sequence2[index - 1];
        if (sum1 % 4 != 0) return false;
    }
    return true;
}

/**
 * Simple utility to get sum of squares of elements.
 */
int getSquaredSum(const std::vector<int>& sequence) {
    int sum = 0;
    for (int x : sequence) sum += x * x;
    return sum;
}

/**
 * Validates condition 5: The total sum of PAF contributions from all four sequences must be zero.
 */
bool step2_condition5(int m, const std::vector<int>& sequenceK, const std::vector<int>& sequenceR,
                       const std::vector<int>& sequenceP, const std::vector<int>& sequenceQ) {
    for (int i = 1; i <= m / 2; ++i) {
        int sum = 0;
        sum += naf_polynomial_decomposition(i, m, const_cast<std::vector<int>&>(sequenceK));
        sum += naf_polynomial_decomposition(m - i, m, const_cast<std::vector<int>&>(sequenceK));
        sum += naf_polynomial_decomposition(i, m, const_cast<std::vector<int>&>(sequenceR));
        sum += naf_polynomial_decomposition(m - i, m, const_cast<std::vector<int>&>(sequenceR));
        sum += naf_polynomial_decomposition(i, m, const_cast<std::vector<int>&>(sequenceP));
        sum += naf_polynomial_decomposition(m - i, m, const_cast<std::vector<int>&>(sequenceP));
        sum += naf_polynomial_decomposition(i, m, const_cast<std::vector<int>&>(sequenceQ));
        sum += naf_polynomial_decomposition(m - i, m, const_cast<std::vector<int>&>(sequenceQ));
        if (sum != 0) return false;
    }
    return true;
}

/**
 * Core search logic:
 * 1. Generates candidate sequences for K, R, P, Q based on target sums.
 * 2. Prunes combinations using condition 3.
 * 3. Iterates through remaining combinations in parallel (MPI).
 * 4. Verifies candidate sets using condition 5 and final NPAF check.
 */
bool getAllSolutions(int n, int m, int k11, int r11, int p11, int q11, int rank, int size) {
    std::vector<int> seq(m);
    std::vector<std::vector<int>> k_sols, r_sols, p_sols, q_sols;

    find_sequences(m, k11, 1, seq, 0, k_sols);
    if (k11 != r11) find_sequences(m, r11, 1, seq, 0, r_sols); else r_sols = k_sols;
    find_sequences(m, p11, 1, seq, 0, p_sols, true);
    if (p11 != q11) find_sequences(m, q11, 1, seq, 0, q_sols, true); else q_sols = p_sols;

    std::vector<std::pair<int, int>> valid_kr;
    for (int i = 0; i < (int)k_sols.size(); ++i)
        for (int j = 0; j < (int)r_sols.size(); ++j)
            if (step2Condition3_k_r(n, m, k_sols[i], r_sols[j])) valid_kr.push_back({i, j});

    std::vector<std::pair<int, int>> valid_pq;
    for (int i = 0; i < (int)p_sols.size(); ++i)
        for (int j = 0; j < (int)q_sols.size(); ++j)
            if (step2Condition3_p_q(n, m, p_sols[i], q_sols[j])) valid_pq.push_back({i, j});

    long long combinations = static_cast<long long>(valid_kr.size()) * valid_pq.size();
    if (rank == 0) std::cout << "Pairs KR: " << valid_kr.size() << " PQ: " << valid_pq.size() << " -> Combinations: " << combinations << std::endl;
    if (combinations > THRESHOLD_V) return false;

    int cond = 4 * n + 2;
    int found_local = 0;
    int found_global = 0;
    
    // Dynamic check interval: scale with n to reduce MPI synchronization overhead 
    // as the search space and work per iteration grow.
    long long check_interval = std::max(1000LL, static_cast<long long>(n) * 500);
    long long counter = 0;

    // Distribute work among MPI ranks
    for (std::size_t idx = rank; idx < valid_kr.size(); idx += size) {
        auto pair_kr = valid_kr[idx];
        const auto &sk = k_sols[pair_kr.first], &sr = r_sols[pair_kr.second];
        int sk_sum = getSquaredSum(sk), sr_sum = getSquaredSum(sr);

        for (const auto& pair_pq : valid_pq) {
            // Periodic check for early exit if another rank found a solution
            if (++counter >= check_interval) {
                counter = 0;
                MPI_Allreduce(&found_local, &found_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                if (found_global) goto exit_loops;
            }

            const auto &sp = p_sols[pair_pq.first], &sq = q_sols[pair_pq.second];
            if (sk_sum + sr_sum + getSquaredSum(sp) + getSquaredSum(sq) == cond) {
                if (step2_condition5(m, sk, sr, sp, sq)) {
                    if (sequence_verify(n, sk, sr, sp, sq)) {
                        std::cout << "Rank " << rank << " found verified solution!" << std::endl;
                        print_sequence(const_cast<std::vector<int>&>(sk), m, 'K');
                        print_sequence(const_cast<std::vector<int>&>(sr), m, 'R');
                        print_sequence(const_cast<std::vector<int>&>(sp), m, 'P');
                        print_sequence(const_cast<std::vector<int>&>(sq), m, 'Q');
                        found_local = 1;
                        MPI_Allreduce(&found_local, &found_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                        goto exit_loops;
                    }
                }
            }
        }
    }

exit_loops:
    if (found_global == 0) MPI_Allreduce(&found_local, &found_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    return found_global != 0;
}

/**
 * Entry point for the algorithm search. Finds quadruples of target sums whose squares sum to 4n+2.
 * Based on n's parity, assigns these sums to K,R and P,Q sequences.
 */
void findQuadruple(int n, int m, int rank, int size) {
    int quadruple_sum = 4 * n + 2;
    if (rank == 0) std::cout << "Target 4*n+2 = " << quadruple_sum << std::endl;
    std::vector<std::vector<int>> results = nsoks(quadruple_sum, 4);
    for (const auto& seq : results) {
        std::vector<int> evens, odds;
        for (int e : seq) if (e % 2 == 0) evens.push_back(e); else odds.push_back(e);
        
        // Construction requires exactly 2 even and 2 odd square roots
        if (evens.size() != 2 || odds.size() != 2) continue;
        
        std::sort(evens.rbegin(), evens.rend());
        std::sort(odds.rbegin(), odds.rend());
        
        int k11, r11, p11, q11;
        if (n % 2 == 0) { 
            k11 = odds[0]; r11 = odds[1]; p11 = evens[0]; q11 = evens[1]; 
        } else { 
            k11 = evens[0]; r11 = evens[1]; p11 = odds[0]; q11 = odds[1]; 
        }
        
        if (rank == 0) std::cout << "NSOKS Quadruple: " << k11 << " " << r11 << " " << p11 << " " << q11 << std::endl;
        if (getAllSolutions(n, m, k11, r11, p11, q11, rank, size)) break;
    }
}
