#include <iostream>
#include <cmath>
#include <numeric>
#include <list>
#include <bitset>
#include <unordered_set>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <set>
#include <tuple>
#include <climits>
#include <algorithm>
#include "utilities.h"
#include <omp.h>

int THRESHOLD_V = 1000000000;
using VectorTuple = std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>>;
int npaf(const int arr[], int n, int s) {
    int i, npaf;
    npaf = 0;
    for (i=0; i<n-s; i++) {
        npaf += arr[i]*arr[i+s];
    }
    return npaf;
}

std::string generateMemoKey(int n, int k, int maxsq) {
//    return std::to_string(n) + "," + std::to_string(k) + "," + std::to_string(maxsq);
    std::ostringstream oss;
    oss << n << "," << k << "," << maxsq;
    return oss.str();
}

int iquo(int n, int k, int &r) {
    r = n%k;
    return n/k;
}

const int DEFAULT_MAX_SQ = INT_MAX;
std::unordered_map<std::string, std::list<std::list<int>>> memo;

std::list< std::list<int>> nsoksAUX(int n, int k, int maxsq = DEFAULT_MAX_SQ) {
    std::string memoKey = generateMemoKey(n, k, maxsq);
    if(memo.count(memoKey)) {
        return memo[memoKey];
    }
    int minTs,trialSquare, maxTs, r;
    maxTs = static_cast<int>(std::sqrt(n-k+1));
    std::list<std::list<int>> seqs;
    std::list<int> maxTs_list;
    if (k==1) {
        if(maxTs*maxTs == n) {
            maxTs_list.push_back(maxTs);
            seqs.push_back(maxTs_list);
            return seqs;
        } else {
            return seqs;
        }
    }
    if (maxTs*maxTs > (n-k+1)) {
        maxTs--;
    }
    maxTs = std::min(maxTs, maxsq);
    int quo = iquo(n, k, r);
    int has_remainder = (r!=0) ? 1:0;
    minTs = static_cast<int>(std::sqrt(quo + has_remainder));
    if(k * minTs * minTs < n) {
        minTs++;
    }
    if(minTs>maxTs) {
        return seqs;
    }
    for (trialSquare = minTs; trialSquare <= maxTs; trialSquare++) {
        std::list< std::list<int>> val = nsoksAUX(n-trialSquare*trialSquare,k-1,trialSquare);
        for( std::list<int>& subSeq : val) {
            std::list<int> seq;
            seq.insert(seq.end(), subSeq.begin(), subSeq.end());
            seq.push_back(trialSquare);
            seqs.push_back(seq);
        }
    }
    memo[memoKey] = seqs;
    return seqs;
}

std::list<std::list<int>> nsoks(int nn, int kk) {
    int kkk, ii, num_of_zeros;
    std::list<std::list<int>> res;
    std::list<std::list<int>> aux_res;
    for(kkk=1; kkk<=kk; kkk++) {
        aux_res = nsoksAUX(nn, kkk);
        res.insert(res.end(), aux_res.begin(), aux_res.end());
    }
    /** if an element of the result has less than kk elements,
     * then replace it with an element with the right number of zeros
     * **/
    for (auto i=res.begin(); i!=res.end(); i++) {
        std::list<int> &seq = *i;
        if(seq.size()!=kk) {
            num_of_zeros = static_cast<int>(kk - seq.size());
            for(ii=1; ii<=num_of_zeros; ii++){
                seq.push_front(0);
            }
        }
    }
    return res;
}

void removeDuplicates(std::vector<int>& vec) {
    // Use an unordered_set to track unique elements
    std::unordered_set<int> uniqueElements;
    std::vector<int> result;

    for (const int& num : vec) {
        // If the element is not in the set, add it to the set and result vector
        if (uniqueElements.insert(num).second) {
            result.push_back(num);
        }
    }
    // Swap the original vector with the result vector
    vec.swap(result);
}

void get_Solution_to_k11_r11(int n, int m, int k11, int r11, int current_j, int currentSumK, int currentSumR,
                             std::vector<int>& sequence_k, std::vector<std::vector<int>>& k_solutions, std::vector<int>& sequence_r,
                             std::vector<std::vector<int>>& r_solutions) {
    if(current_j > m) {
        if(currentSumK == k11) {
            k_solutions.push_back(sequence_k);
        }
        if(currentSumR == r11) {
            r_solutions.push_back(sequence_r);
        }
        return;
    }
    int maxV = static_cast<int>(std::floor((n + 1 - current_j) / m) + 1);
    bool isOdd = (maxV % 2 != 0);
    int v1, v2;
    if(m>=n+1) {
        for(v1 = -1; v1<=1; v1++) {
            for(v2 = -1; v2<=1; v2++) {
                if(v1==0 || v2 ==0) {
                    continue;
                }
                else {
                    if((isOdd && (v1%2 != 0 && v2%2 !=0)) || (!isOdd && (v1%2 == 0 && v2%2 ==0))) {
                        if( (v1+currentSumK <= k11 && v2+currentSumR <= r11) && (std::abs(v1)<=maxV && std::abs(v2)<=maxV)){
                            sequence_k[current_j-1] =  v1;
                            sequence_r[current_j-1] =  v2;
                            get_Solution_to_k11_r11(n, m, k11, r11, current_j+1, currentSumK+v1,
                                                    currentSumR+v2, sequence_k, k_solutions, sequence_r, r_solutions);
                        }
                    }
                }
            }
        }
    } else {
        for(v1 = -maxV; v1<=maxV; v1++) {
            for(v2 = -maxV; v2<=maxV; v2++) {
                if((isOdd && (v1%2 != 0 && v2%2 !=0)) || (!isOdd && (v1%2 == 0 && v2%2 ==0))) {
                    if((v1+currentSumK <= k11 && v2+currentSumR <= r11) && (std::abs(v1)<=maxV && std::abs(v2)<=maxV)){
                        sequence_k[current_j-1] =  v1;
                        sequence_r[current_j-1] =  v2;
                        get_Solution_to_k11_r11(n, m, k11, r11, current_j+1, currentSumK+v1,
                                                currentSumR+v2, sequence_k, k_solutions, sequence_r, r_solutions);
                    }
                }
            }
        }
    }
}


void get_Solution_to_p11_q11(int n, int m, int p11, int q11, int current_j, int currentSumP, int currentSumQ,
                             std::vector<int>& sequence_p, std::vector<std::vector<int>>& p_solutions, std::vector<int>& sequence_q,
                             std::vector<std::vector<int>>& q_solutions) {
    if(current_j > m) {
        if(currentSumP == p11) {
                p_solutions.push_back(sequence_p);
        }
        if(currentSumQ == q11) {
            q_solutions.push_back(sequence_q);
        }
        return;
    }
    int maxV = static_cast<int>(std::floor((n - current_j) / m) + 1);
    bool isOdd = (maxV % 2 != 0);
    int v1, v2;
    if(m>=n+1) {
        if(current_j == m) {
            sequence_p[current_j-1] = 0;
            sequence_q[current_j-1] =  0;
            get_Solution_to_p11_q11(n, m, p11, q11, current_j+1, currentSumP,
                                    currentSumQ, sequence_p, p_solutions, sequence_q, q_solutions);
        } else {
            for(v1 = -1; v1<=1; v1++) {
                for(v2 = -1; v2<=1; v2++) {
                    if(v1==0 || v2 == 0) {
                        continue;
                    }
                    else {
                        if((isOdd && (v1%2 != 0 && v2%2 !=0)) || (!isOdd && (v1%2 == 0 && v2%2 ==0))) {
                            if( (v1+currentSumP <= p11 && v2+currentSumQ <= q11) && (std::abs(v1)<=maxV && std::abs(v2)<=maxV)){
                                sequence_p[current_j-1] =  v1;
                                sequence_q[current_j-1] =  v2;
                                get_Solution_to_p11_q11(n, m, p11, q11, current_j+1, currentSumP+v1,
                                                        currentSumQ+v2, sequence_p, p_solutions, sequence_q, q_solutions);

                            }
                        }
                    }

                }
            }
        }
    } else {
        for(v1 = -maxV; v1<=maxV; v1++) {
            for(v2 = -maxV; v2<=maxV; v2++) {
                if((isOdd && (v1%2 != 0 && v2%2 !=0)) || (!isOdd && (v1%2 == 0 && v2%2 ==0))) {
                    if((v1+currentSumP <= p11 && v2+currentSumQ <= q11) && (std::abs(v1)<=maxV && std::abs(v2)<=maxV)){
                        sequence_p[current_j-1] =  v1;
                        sequence_q[current_j-1] =  v2;
                        get_Solution_to_p11_q11(n, m, p11, q11, current_j+1, currentSumP+v1,
                                                currentSumQ+v2, sequence_p, p_solutions, sequence_q, q_solutions);
                    }
                }
            }
        }
    }
}

void getSolutionToK11(int n, int m, int seq_11, int current_j, std::vector<int>& sequence_k, int currentSum, std::vector<std::vector<int>>& k_solutions) {
    if(current_j > m) {
        if(currentSum == seq_11) {
            k_solutions.push_back(sequence_k);
        }
        return;
    }
    int maxV = static_cast<int>(std::floor((n + 1 - current_j) / m) + 1);
    bool isOdd = (maxV % 2 != 0);
    int k;
    if(m>=n+1) {
        for(k = -1; k<=1; k++) {
            if(k==0) {
                continue;
            }
            if( (isOdd && (k%2 != 0)) || (!isOdd && (k%2 == 0))) {
                if( (k+currentSum <= seq_11) && (std::abs(k)<=maxV)){
                    sequence_k[current_j-1] =  k;
                    getSolutionToK11(n, m, seq_11, current_j+1, sequence_k, currentSum+ k, k_solutions);
                }
            }
        }
    } else {
        for(k = -maxV; k<=maxV; k++) {
            if( (isOdd && (k%2 != 0)) || (!isOdd && (k%2 == 0))) {
                if( (k+currentSum <= seq_11) && (std::abs(k)<=maxV)){
                    sequence_k[current_j-1] =  k;
                    getSolutionToK11(n, m, seq_11, current_j+1, sequence_k, currentSum+ k, k_solutions);
                }
            }
        }
    }

}


void getSolutionToR11(int n, int m, int seq_11, int current_j, std::vector<int>& sequence_r, int currentSum, std::vector<std::vector<int>>& r_solutions) {
    if(current_j > m) {
        if(currentSum == seq_11) {
            r_solutions.push_back(sequence_r);
        }
        return;
    }
    int maxV = static_cast<int>(std::floor((n + 1 - current_j) / m) + 1);
    bool isOdd = (maxV % 2 != 0);
    int k;
    if(m>=n+1) {
        for(k = -1; k<=1; k++) {
            if(k==0) {
                continue;
            }
            if( (isOdd && (k%2 != 0)) || (!isOdd && (k%2 == 0))) {
                if( (k+currentSum <= seq_11) && (std::abs(k)<=maxV)){
                    sequence_r[current_j-1] =  k;
                    getSolutionToR11(n, m, seq_11, current_j+1, sequence_r, currentSum+ k, r_solutions);
                }
            }
        }
    } else {
        for(k = -maxV; k<=maxV; k++) {
            if( (isOdd && (k%2 != 0)) || (!isOdd && (k%2 == 0))) {
                if( (k+currentSum <= seq_11) && (std::abs(k)<=maxV)){
                    sequence_r[current_j-1] =  k;
                    getSolutionToR11(n, m, seq_11, current_j+1, sequence_r, currentSum+ k, r_solutions);
                }
            }
        }
    }

}

void getSolutionToP11(int n, int m, int seq_11, int current_j, std::vector<int>& sequence_p, int currentSum, std::vector<std::vector<int>>& p_solutions) {
    if(current_j > m) {
        if(currentSum == seq_11) {
            p_solutions.push_back(sequence_p);
        }
        return;
    }
    int maxV = static_cast<int>(std::floor((n - current_j) / m) + 1);
    bool isOdd = (maxV % 2 != 0);
    int k;
    if(m>=n+1) {
        if(current_j == m) {
            sequence_p[current_j-1] =  0;
            if(currentSum == seq_11) {
                p_solutions.push_back(sequence_p);
            }
            return;
        }
        for(k = -1; k<=1; k++) {
            if(k==0) {
                continue;
            }
            if( (isOdd && (k%2 != 0)) || (!isOdd && (k%2 == 0))) {
                if((k+currentSum <= seq_11) &&(std::abs(k)<=maxV)){
                    sequence_p[current_j-1] =  k;
                    getSolutionToP11(n, m, seq_11, current_j+1, sequence_p, currentSum+ k, p_solutions);
                }
            }
        }
    } else {
        for(k = -maxV; k<=maxV; k++) {
            if( (isOdd && (k%2 != 0)) || (!isOdd && (k%2 == 0))) {
                if((k+currentSum <= seq_11) &&(std::abs(k)<=maxV)){
                    sequence_p[current_j-1] =  k;
                    getSolutionToP11(n, m, seq_11, current_j+1, sequence_p, currentSum+ k, p_solutions);
                }
            }
        }
    }
}

void getSolutionToQ11(int n, int m, int seq_11, int current_j, std::vector<int>& sequence_q, int currentSum, std::vector<std::vector<int>>& q_solutions) {
    if(current_j > m) {
        if(currentSum == seq_11) {
            q_solutions.push_back(sequence_q);
        }
        return;
    }
    int maxV = static_cast<int>(std::floor((n - current_j) / m) + 1);
    bool isOdd = (maxV % 2 != 0);
    int k;
    if(m>=n+1) {
        if(current_j == m) {
            sequence_q[current_j-1] =  0;
            if(currentSum == seq_11) {
                q_solutions.push_back(sequence_q);
            }
            return;
        }
        for(k = -1; k<=1; k++) {
            if(k==0) {
                continue;
            }
            if( (isOdd && (k%2 != 0)) || (!isOdd && (k%2 == 0))) {
                if((k+currentSum <= seq_11) &&(std::abs(k)<=maxV)){
                    sequence_q[current_j-1] =  k;
                    getSolutionToQ11(n, m, seq_11, current_j+1, sequence_q, currentSum+ k, q_solutions);
                }
            }
        }

    } else {
        for(k = -maxV; k<=maxV; k++) {
            if( (isOdd && (k%2 != 0)) || (!isOdd && (k%2 == 0))) {
                if((k+currentSum <= seq_11) &&(std::abs(k)<=maxV)){
                    sequence_q[current_j-1] =  k;
                    getSolutionToQ11(n, m, seq_11, current_j+1, sequence_q, currentSum+ k, q_solutions);
                }
            }
        }
    }
}

bool step2Condition3_k_r(int n, int m, std::vector<int>& sequence1, std::vector<int>& sequence2) {
    int i, j, index, sum;
    bool isMultipleOfM = (n%m == 0);
    int maxV;
    for(i=1; i<=m; i++) {
        maxV = static_cast<int>(std::floor((n + 1 - i) / m) + 1);
        if(std::abs(sequence1[i-1])<=maxV && std::abs(sequence2[i-1])<=maxV) {
            if(maxV%2==0) {
                if((sequence1[i-1]%2!=0)||(sequence2[i-1]%2!=0)) {
                    return false;
                }
            } else {
                if((sequence1[i-1]%2==0)||(sequence2[i-1]%2==0)) {
                    return false;
                }
            }
        }
    }

    for(j=2; j<=m; j++) {
        index = n+2-j;
        sum = sequence1[j-1]+sequence2[j-1]+sequence1[index-1]+sequence2[index-1];
//        printf("k_%d%d + r_%d%d + k_%d%d + r_%d%d = %d\n",j,m,j,m,index,m,index,m,sum);
        if(j==(n+1) && index==1) {
            if (isMultipleOfM) {
                if((sum % 4) != 0){
                    return false;
                }
            } else {
                if((sum % 4) != 2){
                    return false;
                }
            }
        } else {
            if((sum%4)!=0) {
                return false;
            }
        }
    }
    return true;
}

//bool step2Condition3_p_q(int n, int m, std::vector<int> sequence_p, std::vector<int> sequence_q) {
//    int i, j, index, sum1;
//
//    for(j=1; j<=m; j++) {
//        index = n+1-j;
//        sum1 = sequence_p[j-1]+sequence_q[j-1]+sequence_p[index+1]+sequence_q[index+1];
//        if(sum1%4!=0) {
//            return false;
//        }
//    }
//    return true;
//}

bool step2Condition3_p_q(int n, int m, std::vector<int>& sequence1, std::vector<int>& sequence2) {
    int i, j, index, sum1;
    int maxV;
    for(i=1; i<m; i++) {
        maxV = static_cast<int>(std::floor((n - i) / m) + 1);
        if(std::abs(sequence1[i-1])<=maxV && std::abs(sequence2[i-1])<=maxV)  {
            if(maxV%2==0) {
                if((sequence1[i-1]%2!=0)||(sequence2[i-1]%2!=0)) {
                    return false;
                }
            } else {
                if((sequence1[i-1]%2==0)||(sequence2[i-1]%2==0)) {
                    return false;
                }
            }
        }
    }

    for(j=1; j<m; j++) {
        index = n+1-j;
        sum1 = sequence1[j-1]+sequence2[j-1]+sequence1[index-1]+sequence2[index-1];
        if(sum1%4!=0) {
            return false;
        }
    }
    return true;
}

int getSquaredSum(int m, std::vector<int>& sequence) {
    int sum =0;
    int i;
    for(i=0; i<m; i++) {
        sum += sequence[i] * sequence[i];
    }
    return sum;
}

//int getSquaredSum(int m, std::vector<int>& sequence) {
    // Ensure m is not larger than the sequence size
//    m = std::min(m, static_cast<int>(sequence.size()));
//
//    return std::accumulate(sequence.begin(), sequence.begin() + m, 0,
//                           [](int sum, int value) {
//                               return sum + value * value;
//                           });
//}

bool step2_condition5 (int m, std::vector<int>& sequenceK, std::vector<int>& sequenceR,
                       std::vector<int>& sequenceP, std::vector<int>& sequenceQ) {
    int i, sum;
    for (i=1; i<=m/2; i++) {
        sum = 0;
        sum += naf_polynomial_decomposition (i, m, sequenceK);
        sum += naf_polynomial_decomposition (m-i, m, sequenceK);
        sum += naf_polynomial_decomposition (i, m, sequenceR);
        sum += naf_polynomial_decomposition (m-i, m, sequenceR);
        sum += naf_polynomial_decomposition (i, m, sequenceP);
        sum += naf_polynomial_decomposition (m-i, m, sequenceP);
        sum += naf_polynomial_decomposition (i, m, sequenceQ);
        sum += naf_polynomial_decomposition (m-i, m, sequenceQ);
        if (sum != 0) {
            return false;
        }
    }
    return true;
}

bool step2_condition4 (int n, int m, std::vector<int>& filtered_k_solutions, std::vector<int>& filtered_r_solutions,
                              std::vector<int>& filtered_p_solutions, std::vector<int>& filtered_q_solutions,
                              std::vector<std::vector<int>>& k_solutions,  std::vector<std::vector<int>>& r_solutions,
                              std::vector<std::vector<int>>& p_solutions,  std::vector<std::vector<int>>& q_solutions) {

    int i, j, k, l, result;
    int cond = 4*n+2;
//    std::set<VectorTuple> uniqueCombinations;
    for (i=0; i<filtered_k_solutions.size(); i++) {
        for(j=0; j<filtered_r_solutions.size(); j++) {
            for(k=0; k<filtered_p_solutions.size(); k++) {
                for(l=0; l<filtered_q_solutions.size(); l++) {
                    auto& sequenceK = k_solutions[filtered_k_solutions[i]];
                    auto& sequenceR = r_solutions[filtered_r_solutions[j]];
                    auto& sequenceP = p_solutions[filtered_p_solutions[k]];
                    auto& sequenceQ = q_solutions[filtered_q_solutions[l]];
//                    std::sort(sequenceK.begin(), sequenceK.end());
//                    std::sort(sequenceR.begin(), sequenceR.end());
//                    std::sort(sequenceP.begin(), sequenceP.end());
//                    std::sort(sequenceQ.begin(), sequenceQ.end());
//                    if(!uniqueCombinations.insert(std::make_tuple(sequenceK, sequenceR, sequenceP, sequenceQ)).second) {
//                        continue;
//                    }

                    result = getSquaredSum(m, sequenceK) +
                             getSquaredSum(m, sequenceR) +
                             getSquaredSum(m, sequenceP) +
                             getSquaredSum(m, sequenceQ);
                    if (result == cond) {
                        bool cond5 = step2_condition5(m, sequenceK,sequenceR,sequenceP,sequenceQ);
                        if (cond5) {
                                if(m==n+1) {
                                    print_sequence(sequenceK, m, 'K');
                                    print_sequence(sequenceR, m, 'R');
                                    print_sequence(sequenceP, m, 'P');
                                    print_sequence(sequenceQ, m, 'Q');
                                    printf("\n");
                                }
                                return true;
                            }
                        }
                    }
                }
            }
        }
    return false;
}

bool getAllSolutions(int n, int m, int k11, int r11, int p11, int q11) {
    printf("k11: %d\n", k11);
    printf("r11: %d\n", r11);
    printf("p11: %d\n", p11);
    printf("q11: %d\n", q11);
    std::vector<int> sequence_k(m);
    std::vector<int> sequence_r(m);
    std::vector<int> sequence_p(m);
    std::vector<int> sequence_q(m);
    std::vector<std::vector<int>> k_solutions;
    std::vector<std::vector<int>> r_solutions;
    std::vector<std::vector<int>> p_solutions;
    std::vector<std::vector<int>> q_solutions;

    std::vector<int> filtered_k_solutions;
    std::vector<int> filtered_r_solutions;
    std::vector<int> filtered_p_solutions;
    std::vector<int> filtered_q_solutions;
    if(k11!=r11) {
        get_Solution_to_k11_r11(n, m, k11, r11, 1, 0,0,
                                sequence_k, k_solutions, sequence_r, r_solutions);
    } else {
        getSolutionToK11(n, m, k11, 1, sequence_k, 0, k_solutions);
        r_solutions = k_solutions;
    }
    if(p11!=q11) {
        get_Solution_to_p11_q11(n,m,p11, q11, 1, 0, 0,
                                sequence_p, p_solutions, sequence_q, q_solutions);
    }
    else {
        getSolutionToP11(n, m, p11, 1, sequence_p, 0, p_solutions);
        q_solutions = p_solutions;
    }
    int i, j;
    bool cond_k_r, cond_p_q;
    for(i=0; i<k_solutions.size(); i++) {
        for(j=0; j<r_solutions.size(); j++) {
            cond_k_r = step2Condition3_k_r(n,m, k_solutions[i], r_solutions[j]);
            if(cond_k_r) {
                filtered_k_solutions.push_back(i);
                filtered_r_solutions.push_back(j);
            }
        }
    }
    for(i=0; i<p_solutions.size(); i++) {
        for(j=0; j<q_solutions.size(); j++) {
            cond_p_q = step2Condition3_p_q(n,m, p_solutions[i], q_solutions[j]);
            if(cond_p_q) {
                filtered_p_solutions.push_back(i);
                filtered_q_solutions.push_back(j);
            }
        }
    }
    removeDuplicates(filtered_k_solutions);
    if(k11==r11) {
        filtered_r_solutions = filtered_k_solutions;
    }
    else {
        removeDuplicates(filtered_r_solutions);
    }
    removeDuplicates(filtered_p_solutions);
    if(p11 == q11) {
        filtered_q_solutions = filtered_p_solutions;
    }
    else {
        removeDuplicates(filtered_q_solutions);
    }

    if(filtered_k_solutions.size()*filtered_r_solutions.size()*filtered_p_solutions.size()* filtered_q_solutions.size()>THRESHOLD_V) {
        return false;
    }

    printf("k11 solutions found: %zu \n", filtered_k_solutions.size());
    printf("r11 solutions found: %zu \n", filtered_r_solutions.size());
    printf("p11 solutions found: %zu \n", filtered_p_solutions.size());
    printf("q11 solutions found: %zu \n", filtered_q_solutions.size());
    return step2_condition4(n, m, filtered_k_solutions, filtered_r_solutions, filtered_p_solutions, filtered_q_solutions,
                     k_solutions, r_solutions, p_solutions, q_solutions);
}

void findQuadruple(int n, int m) {
    int i, k11, r11, p11, q11;
    int quadruple_sum = 4*n+2;
    printf("4*n+2 = %d\n", quadruple_sum);
    int sums[4];
    // check 4n+2
    std::list< std::list<int>> result = nsoks(quadruple_sum, 4);
    for(auto itr = result.begin(); itr != result.end(); itr++) {
        i=0;
        std::vector<int> evens;
        std::vector<int> odds;
        std::list<int>& seq = *itr;
        printf("NSOKS: ");
        for(int& e: seq) {
            sums[i++] = e;
            printf("%d ", e);
            if(e%2==0) {
                evens.push_back(e);
            } else {
                odds.push_back(e);
            }
        }
        printf("\n");
        if(evens.size() != 2 || odds.size()!=2) {
            continue;
        }
        if(evens[0] < evens[1]) {
            evens[0] = evens[0] ^ evens[1];
            evens[1] = evens[0] ^ evens[1];
            evens[0] = evens[0] ^ evens[1];
        }
        if(odds[0] < odds[1]) {
            odds[0] = odds[0] ^ odds[1];
            odds[1] = odds[0] ^ odds[1];
            odds[0] = odds[0] ^ odds[1];
        }
        if(n%2==0) {
            printf("n is a even\n");
            k11 = odds[0];
            r11 = odds[1];
            p11 = evens[0];
            q11 = evens[1];
        } else {
            printf("n is an odd\n");
            k11 = evens[0];
            r11 = evens[1];
            p11 = odds[0];
            q11 = odds[1];
        }

        int j;
        for(j=2; j<=m; j++) {
            printf("current m=%d \n", j);
            if(getAllSolutions(n, j, k11, r11, p11, q11)){
                printf("--------------------------------------------------------------------------\n");
                printf("\n");
                continue;
            } else {
                break;
            }
        }
    }
}

int main(int argc, char *argv[]) {
    int i,j;
    int n = 9;

    for(i=n; i<=n; i++) {
        int m = n+1;
        printf("n=%d\n", i);
        findQuadruple(i, m);
        printf("--------------------------------------------------------------------------\n");
    }
    return 0;
}