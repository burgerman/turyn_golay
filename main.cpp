#include <iostream>
#include <cmath>
#include <numeric>
#include <list>
#include <bitset>
#include <unordered_set>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <climits>
#include <algorithm>
#include "utilities.h"
#include <omp.h>

void construct_quadruple (const int a[], const int b[], int** quadruple, int n) {
    int i;
    for(i=0; i<4; i++) {
        quadruple[i] = (int *)malloc((n+1) * sizeof(int));
    }
    quadruple[0][0] = quadruple[2][0] =1;
    quadruple[1][0] = quadruple[3][0] = -1;
    for(i=1; i<n+1; i++) {
        quadruple[0][i] = quadruple[1][i] = a[i-1];
        quadruple[2][i] = quadruple[3][i] = b[i-1];
    }
}

int npaf(const int arr[], int n, int s) {
    int i, npaf;
    npaf = 0;
    for (i=0; i<n-s; i++) {
        npaf += arr[i]*arr[i+s];
    }
    return npaf;
}

int sum_arr(int arr[], int size) {
    return std::accumulate(arr, arr+size, 0);
}

void reversal(const int arr[], int* reversed_arr, int len) {
    int start = 0;
    int end = len - 1;
    while (start < end) {
        reversed_arr[start] = arr[end];
        reversed_arr[end] = arr[start];
        start++;
        end--;
    }
}

void addition (const int a[], const int b[], int* added_arr, int len) {
    int i;
    for (i=0; i<len; i++) {
        added_arr[i] = a[i] + b[i];
    }
}

void subtraction (const int a[], const int b[], int* sub_arr, int len) {
    int i;
    for (i=0; i<len; i++) {
        sub_arr[i] = a[i] - b[i];
    }
}

void tensor_product (const int a[], const int b[], int* tp_arr, int len_m, int len_n ) {
    int i, j;
    for (i = 0; i < len_m; i++) {
        for(j=0; j<len_n; j++) {
            tp_arr[i*len_n+j] = a[i]*b[j];
        }
    }
}

void seq_combine (const int a[], const int b[], int* combined_seq, int len) {
    int i, j;
    int idx = 0;
    for (i=0; i<len; i++) {
        combined_seq[idx++] = a[i];
        combined_seq[idx++] = b[i];
    }
    combined_seq[idx] = a[len];
}

bool diophantine_sum_func (int** quadruple, int n) {
    int x, y, z, w, seq_sum, condition;
    x = sum_arr(quadruple[0], n);
    y = sum_arr(quadruple[1], n);
    z = sum_arr(quadruple[2], n);
    w = sum_arr(quadruple[3], n);
    seq_sum = x*x + y*y + z*z + w*w;
    condition = 4*n;
    if(seq_sum != condition) {
        return false;
    }
    else {
        return true;
    }
}

bool turyn_quadruples_func (int** quadruple, int n) {
    int i, seq_sum;
    for (i=1; i<=n; i++) {
        seq_sum = npaf(quadruple[0], n+1, i)+ npaf(quadruple[1], n+1, i) + npaf(quadruple[2], n, i) + npaf(quadruple[3], n, i);
        if(seq_sum!=0) {
            return false;
        }
    }
    return true;
}

void print_bs (int** quadruple, int n) {
    int i;
    for (i=0; i<=n; i++) {
        printf("%d", quadruple[0][i]);
    }
    printf("\n");
    for (i=0; i<=n; i++) {
        printf("%d", quadruple[1][i]);
    }
    printf("\n");
    for (i=0; i<n; i++) {
        printf("%d", quadruple[2][i]);
    }
    printf("\n");
    for (i=0; i<n; i++) {
        printf("%d", quadruple[3][i]);
    }
    printf("\n");
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

const int N = 5;
const int N_PLUS_ONE = N+1;

std::unordered_map<std::string, std::unordered_set<std::bitset<N_PLUS_ONE>>> seq_n_plus_one_map;
std::unordered_map<std::string, std::unordered_set<std::bitset<N>>> seq_n_map;

void seqFinder(int n, int sum, int solution, std::bitset<N> seq, std::unordered_set<std::bitset<N>> &solution_seqs) {
    std::string key = std::to_string(n) + "_" + std::to_string(solution);
    if(seq_n_map.count(key)) {
        solution_seqs = seq_n_map[key];
        return;
    }

    if(n == 0) {
        if(sum == solution) {
            solution_seqs.insert(seq);
        }
        return;
    }
    seq[n-1] = false;
    seqFinder(n-1, sum-1, solution, seq, solution_seqs);
    seq[n-1] = true;
    seqFinder(n-1, sum+1, solution, seq, solution_seqs);

    if(n==N) {
        seq_n_map[key] = solution_seqs;
    }

}

void seqPlusOneFinder(int n, int sum, int solution, std::bitset<N_PLUS_ONE> &seq, std::unordered_set<std::bitset<N_PLUS_ONE>> &solution_seqs) {
    if(n == N_PLUS_ONE) {
        std::string key = std::to_string(seq.size()) + "_" + std::to_string(solution);
        if(seq_n_plus_one_map.count(key)) {
            solution_seqs = seq_n_plus_one_map[key];
            return;
        } else {
            seq_n_plus_one_map[key] = solution_seqs;
        }
    }
    if(n == 0) {
        if(sum == solution) {
            solution_seqs.insert(seq);
        }
    } else {
        seq[n-1] = false;
        seqPlusOneFinder(n-1, sum-1, solution, seq, solution_seqs);
        seq[n-1] = true;
        seqPlusOneFinder(n-1, sum+1, solution, seq, solution_seqs);
    }
}

void get_Solution_to_k11_r11(int n, int m, int k11, int r11, int current_j, int currentSumK, int currentSumR,
                             std::vector<int>& sequence_k, std::vector<std::vector<int>>& k_solutions, std::vector<int>& sequence_r,
                             std::vector<std::vector<int>>& r_solutions) {
    if(current_j > m) {
        if(currentSumK == k11) {
            #pragma omp critical
            {
                k_solutions.push_back(sequence_k);
            }

        }
        if(currentSumR == r11) {
            #pragma omp critical
            {
                r_solutions.push_back(sequence_r);
            }
        }
        return;
    }
    int maxV = static_cast<int>(std::floor((n + 1 - current_j) / m) + 1);
    bool isOdd = (maxV % 2 != 0);
    int v1, v2;
    if(m>=n+1) {
        #pragma omp parallel for collapse(2) schedule(dynamic) private(v1, v2) shared(sequence_k, sequence_r)
        for(v1 = -1; v1<=1; v1++) {
            for(v2 = -1; v2<=1; v2++) {
                if((isOdd && (v1%2 != 0 && v2%2 !=0)) || (!isOdd && (v1%2 == 0 && v2%2 ==0))) {
                    if( (v1+currentSumK <= k11 && v2+currentSumR <= r11) && (std::abs(v1)<=maxV && std::abs(v2)<=maxV)){
                        #pragma omp task firstprivate(v1, v2)
                        {
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
        #pragma omp parallel for collapse(2) schedule(dynamic) private(v1, v2) shared(sequence_k, sequence_r)
        for(v1 = -maxV; v1<=maxV; v1++) {
            for(v2 = -maxV; v2<=maxV; v2++) {
                if((isOdd && (v1%2 != 0 && v2%2 !=0)) || (!isOdd && (v1%2 == 0 && v2%2 ==0))) {
                    if((v1+currentSumK <= k11 && v2+currentSumR <= r11) && (std::abs(v1)<=maxV && std::abs(v2)<=maxV)){
                        #pragma omp task firstprivate(v1, v2)
                        {
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
}


void get_Solution_to_p11_q11(int n, int m, int p11, int q11, int current_j, int currentSumP, int currentSumQ,
                             std::vector<int>& sequence_p, std::vector<std::vector<int>>& p_solutions, std::vector<int>& sequence_q,
                             std::vector<std::vector<int>>& q_solutions) {
    if(current_j > m) {
        if(currentSumP == p11) {
            #pragma omp critical
            {
                p_solutions.push_back(sequence_p);
            }

        }
        if(currentSumQ == q11) {
            #pragma omp critical
            {
                q_solutions.push_back(sequence_q);
            }
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
            #pragma omp parallel for collapse(2) schedule(dynamic) private(v1, v2) shared(sequence_p, sequence_q)
            for(v1 = -1; v1<=1; v1++) {
                for(v2 = -1; v2<=1; v2++) {
                    if((isOdd && (v1%2 != 0 && v2%2 !=0)) || (!isOdd && (v1%2 == 0 && v2%2 ==0))) {
                        if( (v1+currentSumP <= p11 && v2+currentSumQ <= q11) && (std::abs(v1)<=maxV && std::abs(v2)<=maxV)){
                            #pragma omp task firstprivate(v1, v2)
                            {
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
        #pragma omp parallel for collapse(2) schedule(dynamic) private(v1, v2) shared(sequence_p, sequence_q)
        for(v1 = -maxV; v1<=maxV; v1++) {
            for(v2 = -maxV; v2<=maxV; v2++) {
                if((isOdd && (v1%2 != 0 && v2%2 !=0)) || (!isOdd && (v1%2 == 0 && v2%2 ==0))) {
                    if((v1+currentSumP <= p11 && v2+currentSumQ <= q11) && (std::abs(v1)<=maxV && std::abs(v2)<=maxV)){
                        #pragma omp task firstprivate(v1, v2)
                        {
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
    if(m==n+1) {
        for(k = -1; k<=1; k++) {
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
    if(m==n+1) {
        for(k = -1; k<=1; k++) {
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
    if(m==n+1) {
        if(current_j == m) {
            sequence_p[current_j-1] =  0;
            if(currentSum == seq_11) {
                p_solutions.push_back(sequence_p);
            }
            return;
        }
        for(k = -1; k<=1; k++) {
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
    if(m==n+1) {
        if(current_j == m) {
            sequence_q[current_j-1] =  0;
            if(currentSum == seq_11) {
                q_solutions.push_back(sequence_q);
            }
            return;
        }
        for(k = -1; k<=1; k++) {
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

//int getSquaredSum(int m, std::vector<int>& sequence) {
//    int sum =0;
//    int i;
//    for(i=0; i<m; i++) {
//        sum += sequence[i] * sequence[i];
//    }
//    return sum;
//}

int getSquaredSum(int m, std::vector<int>& sequence) {
    // Ensure m is not larger than the sequence size
    m = std::min(m, static_cast<int>(sequence.size()));

    return std::accumulate(sequence.begin(), sequence.begin() + m, 0,
                           [](int sum, int value) {
                               return sum + value * value;
                           });
}

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
    bool found = false;
    int cond = 4*n+2;
    #pragma omp parallel for collapse(4) schedule(dynamic) num_threads(8) private(i,j,k,l,result) shared(found)
    for (i=0; i<filtered_k_solutions.size(); i++) {
        for(j=0; j<filtered_r_solutions.size(); j++) {
            for(k=0; k<filtered_p_solutions.size(); k++) {
                for(l=0; l<filtered_q_solutions.size(); l++) {
                    auto& sequenceK = k_solutions[filtered_k_solutions[i]];
                    auto& sequenceR = r_solutions[filtered_r_solutions[j]];
                    auto& sequenceP = p_solutions[filtered_p_solutions[k]];
                    auto& sequenceQ = q_solutions[filtered_q_solutions[l]];
                    result = getSquaredSum(m, sequenceK) +
                    getSquaredSum(m, sequenceR) +
                    getSquaredSum(m, sequenceP) +
                    getSquaredSum(m, sequenceQ);
                    if (result == cond) {
                        bool cond5 = step2_condition5(m, sequenceK,sequenceR,sequenceP,sequenceQ);
                        if (cond5) {
                            #pragma omp critical
                            {
                                if(!found) {
                                    print_sequence(sequenceK, m, 'K');
                                    print_sequence(sequenceR, m, 'R');
                                    print_sequence(sequenceP, m, 'P');
                                    print_sequence(sequenceQ, m, 'Q');
                                    printf("\n");
                                    found = true;
                                    return found;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return found;
}

bool getAllSolutions(int n, int m, int k11, int r11, int p11, int q11) {
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
    get_Solution_to_k11_r11(n, m, k11, r11, 1, 0,0,
                            sequence_k, k_solutions, sequence_r, r_solutions);
    get_Solution_to_p11_q11(n,m,p11, q11, 1, 0, 0,
                            sequence_p, p_solutions, sequence_q, q_solutions);
//    getSolutionToK11(n, m, k11, 1, sequence_k, 0, k_solutions);
//    getSolutionToR11(n, m, r11, 1, sequence_r, 0, r_solutions);
//    getSolutionToP11(n, m, p11, 1, sequence_p, 0, p_solutions);
//    getSolutionToQ11(n, m, q11, 1, sequence_q, 0, q_solutions);
    int i, j;
    bool cond_k_r, cond_p_q, res;
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
    removeDuplicates(filtered_r_solutions);
    removeDuplicates(filtered_p_solutions);
    removeDuplicates(filtered_q_solutions);
    printf("K solutions found: %zu \n", filtered_k_solutions.size());
    printf("R solutions found: %zu \n", filtered_r_solutions.size());
    printf("P solutions found: %zu \n", filtered_p_solutions.size());
    printf("Q solutions found: %zu \n", filtered_q_solutions.size());
    res = step2_condition4(n, m, filtered_k_solutions, filtered_r_solutions, filtered_p_solutions, filtered_q_solutions,
                     k_solutions, r_solutions, p_solutions, q_solutions);
    if(res) {
        return true;
    }
    else {
        return false;
    }
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
        printf("number of evens: %zu\n", evens.size());
        printf("number of odds : %zu\n", odds.size());
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

        printf("k11: %d\n", k11);
        printf("r11: %d\n", r11);
        printf("p11: %d\n", p11);
        printf("q11: %d\n", q11);

        int j;
        for(j=n+1; j<=m; j++) {
            if(getAllSolutions(n, j, k11, r11, p11, q11)){
                printf("--------------------------------------------------------------------------\n");
                printf("\n");
                break;
            }
        }
    }
}

int main(int argc, char *argv[]) {
    int i,j;
    int n = 9;
//    int A[] = {-1,1,-1,1,-1,-1,1,1,1,-1,-1,-1,1,1,1,1,1,-1,-1,-1};
//    int B[] ={1,-1,-1,1,-1,-1,1,-1,-1,1,1,1,-1,1,1,1,1,-1,1,-1};
//    int C[] ={1,-1,1,1,1,1,1,-1,1,1,1,-1,-1,1,1,-1,1,-1,-1};
//    int D[] ={1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,1,1,-1,1,1,-1};
    int A[] = {1, 1, 1, 1, 1, -1};
    int B[] = {1, 1, -1, -1, 1, 1};
    int C[] = {1, 1, -1, 1, -1};
    int D[] = {1, -1, 1, -1, -1};
//    int A[] = {1,1,1,1,1,1,-1,1,-1,1,-1,1};
//    int B[] = {1,1,1,1,-1,-1,1,-1,-1,1,-1,-1};
//    int C[] = {1,1,1,-1,-1,-1,1,1,-1,1,1};
//    int D[] = {1,-1,1,-1,-1,1,1,1,-1,-1,1};
    int mod_start = 2;
    int mod_end = n+1;
    int len_A = sizeof(A)/ sizeof(int);
    int len_B = sizeof(B)/ sizeof(int);
    int len_C = sizeof(C)/ sizeof(int);
    int len_D = sizeof(D)/ sizeof(int);
//    int arr_res_A[m];
//    int arr_res_B[m];
//    int arr_res_C[m];
//    int arr_res_D[m];
    int nk, nr, np, nq, sum;
    bool k_res, r_res, p_res, q_res;
    for(i=n; i<=9; i++) {
        int m = n+1;
        printf("n=%d\n", i);
        findQuadruple(i, m);
        printf("--------------------------------------------------------------------------\n");
    }

//    filter_array_element_by_mod(A, len_A, m, arr_res_A, 'k');
//    filter_array_element_by_mod(B, len_B, m, arr_res_B, 'r');
//    bool isCon3 = step2_cond3_k_r(n, m, arr_res_A, arr_res_B);
//    if(isCon3) {
//        printf("Condition 3 matches");
//    }



//    for (i=1; i<=m; i++) {
//        printf("When m is: %d \n", i);
//        printf("Sequence A: \n");
//        k_res = filter_array_element_by_mod(A, len_A, i, arr_res_A, 'k');
//        if(!k_res){
//            printf("Seq A doesn't satisfy the condition in step 2");
//            break;
//        } else if (!(0<=arr_res_A[1] && arr_res_A[1]<=arr_res_A[0]) && m==2) {
//            printf("Seq A doesn't satisfy the condition in step 1");
//            break;
//        }
//        nk = naf_polynomial_decomposition (0, i, arr_res_A, 'k');
//        printf("N_K(%d)=%d\n", 0, nk);
//        printf("Sequence B: \n");
//        r_res = filter_array_element_by_mod(B, len_B, i, arr_res_B, 'r');
//        if(!r_res){
//            printf("Seq B doesn't satisfy the condition in step 2");
//            break;
//        } else if (arr_res_B[1] > arr_res_A[0] && m==2) {
//            printf("Seq B doesn't satisfy the condition in step 1");
//            break;
//        }
//        nr = naf_polynomial_decomposition (0, i, arr_res_B, 'r');
//        printf("N_R(%d)=%d\n", 0, nr);
//        printf("Sequence C: \n");
//        p_res = filter_array_element_by_mod(C, len_C, i, arr_res_C, 'p');
//        if(!p_res){
//            printf("Seq C doesn't satisfy the condition in step 2");
//            break;
//        } else if (arr_res_C[1]>arr_res_C[0] && m==2) {
//            printf("Seq C doesn't satisfy the condition in step 1");
//            break;
//        }
//        np = naf_polynomial_decomposition (0, i, arr_res_C, 'p');
//        printf("N_P(%d)=%d\n", 0, np);
//        printf("Sequence D: \n");
//        q_res = filter_array_element_by_mod(D, len_D, i, arr_res_D, 'q');
//        if(!q_res){
//            printf("Seq D doesn't satisfy the condition in step 2");
//            break;
//        } else if (arr_res_D[1] > arr_res_D[0] && m==2) {
//            printf("Seq D doesn't satisfy the condition in step 1");
//            break;
//        }
//        nq = naf_polynomial_decomposition (0, i, arr_res_D, 'q');
//        printf("N_Q(%d)=%d\n", 0, nq);
//        sum = nk+nr+np+nq;
//        printf("Sum=%d\n", sum);
//        print_sequence(arr_res_A, i, 'K');
//        print_sequence(arr_res_B, i, 'R');
//        print_sequence(arr_res_C, i, 'P');
//        print_sequence(arr_res_D, i, 'Q');
//        printf("\n");
//        printf("\n");
//    }

//    int sum_part1=0;
//    for (i=1; i<=m/2; i++) {
//        printf("When s is: %d \n", i);
//        printf("Sequence A: \n");
//        k_res = filter_array_element_by_mod(A, len_A, m, arr_res_A, 'k');
//        if(!k_res){
//            printf("Seq A doesn't satisfy the condition in step 2\n");
//            break;
//        } else if (!(0<=arr_res_A[1] && arr_res_A[1]<=arr_res_A[0])) {
//            printf("Seq A doesn't satisfy the condition in step 1\n");
//            break;
//        }
//        nk = naf_polynomial_decomposition (i, m, arr_res_A, 'k');
//        printf("N_K(%d)=%d\n", 0, nk);
//        printf("Sequence B: \n");
//        r_res = filter_array_element_by_mod(B, len_B, m, arr_res_B, 'r');
//        if(!r_res){
//            printf("Seq B doesn't satisfy the condition in step 2\n");
//            break;
//        } else if (arr_res_B[1] > arr_res_A[0]) {
//            printf("Seq B doesn't satisfy the condition in step 1\n");
//            break;
//        }
//        nr = naf_polynomial_decomposition (i, m, arr_res_B, 'r');
//        printf("N_R(%d)=%d\n", 0, nr);
//        printf("Sequence C: \n");
//        p_res = filter_array_element_by_mod(C, len_C, m, arr_res_C, 'p');
//        if(!p_res){
//            printf("Seq C doesn't satisfy the condition in step 2\n");
//            break;
//        } else if (!(0<=arr_res_C[1] && arr_res_C[1]<=arr_res_C[0])) {
//            printf("Seq C doesn't satisfy the condition in step 1\n");
//            break;
//        }
//        np = naf_polynomial_decomposition (i, m, arr_res_C, 'p');
//        printf("N_P(%d)=%d\n", 0, np);
//        printf("Sequence D: \n");
//        q_res = filter_array_element_by_mod(D, len_D, m, arr_res_D, 'q');
//        if(!q_res){
//            printf("Seq D doesn't satisfy the condition in step 2\n");
//            break;
//        } else if (arr_res_D[1] > arr_res_D[0]) {
//            printf("Seq D doesn't satisfy the condition in step 1\n");
//            break;
//        }
//        nq = naf_polynomial_decomposition (i, m, arr_res_D, 'q');
//        printf("N_Q(%d)=%d\n", 0, nq);
//        sum = nk+nr+np+nq;
//        printf("Sum=%d\n", sum);
//        print_sequence(arr_res_A, m, 'K');
//        print_sequence(arr_res_B, m, 'R');
//        print_sequence(arr_res_C, m, 'P');
//        print_sequence(arr_res_D, m, 'Q');
//        printf("\n");
//        printf("\n");
//
//        sum_part1 += sum;
//    }
//    printf("sum_part1: %d \n", sum_part1);
//
//    int sum_part2 = 0;
//    for (i=1; i<=m/2; i++) {
//        printf("When s is: %d \n", m-i);
//        printf("Sequence A: \n");
//        k_res = filter_array_element_by_mod(A, len_A, m, arr_res_A, 'k');
//        if(!k_res){
//            printf("Seq A doesn't satisfy the condition in step 2\n");
//            break;
//        } else if (!(0<=arr_res_A[1] && arr_res_A[1]<=arr_res_A[0])) {
//            printf("Seq A doesn't satisfy the condition in step 1\n");
//            break;
//        }
//        nk = naf_polynomial_decomposition (m-i, m, arr_res_A, 'k');
//        printf("N_K(%d)=%d\n", 0, nk);
//        printf("Sequence B: \n");
//        r_res = filter_array_element_by_mod(B, len_B, m, arr_res_B, 'r');
//        if(!r_res){
//            printf("Seq B doesn't satisfy the condition in step 2\n");
//            break;
//        } else if (arr_res_B[1] > arr_res_A[0]) {
//            printf("Seq B doesn't satisfy the condition in step 1\n");
//            break;
//        }
//        nr = naf_polynomial_decomposition (m-i, m, arr_res_B, 'r');
//        printf("N_R(%d)=%d\n", 0, nr);
//        printf("Sequence C: \n");
//        p_res = filter_array_element_by_mod(C, len_C, m, arr_res_C, 'p');
//        if(!p_res){
//            printf("Seq C doesn't satisfy the condition in step 2\n");
//            break;
//        } else if (!(0<=arr_res_C[1] && arr_res_C[1]<=arr_res_C[0])) {
//            printf("Seq C doesn't satisfy the condition in step 1\n");
//            break;
//        }
//        np = naf_polynomial_decomposition (m-i, m, arr_res_C, 'p');
//        printf("N_P(%d)=%d\n", 0, np);
//        printf("Sequence D: \n");
//        q_res = filter_array_element_by_mod(D, len_D, m, arr_res_D, 'q');
//        if(!q_res){
//            printf("Seq D doesn't satisfy the condition in step 2\n");
//            break;
//        } else if (arr_res_D[1] > arr_res_D[0]) {
//            printf("Seq D doesn't satisfy the condition in step 1\n");
//            break;
//        }
//        nq = naf_polynomial_decomposition (m-i, m, arr_res_D, 'q');
//        printf("N_Q(%d)=%d\n", 0, nq);
//        sum = nk+nr+np+nq;
//        printf("Sum=%d\n", sum);
//        print_sequence(arr_res_A, m, 'K');
//        print_sequence(arr_res_B, m, 'R');
//        print_sequence(arr_res_C, m, 'P');
//        print_sequence(arr_res_D, m, 'Q');
//        printf("\n");
//        printf("\n");
//        sum_part2 += sum;
//    }
//    printf("sum_part2: %d \n", sum_part2);
//    if(sum_part1+sum_part2==0) {
//        printf("step v matched\n");
//    } else {
//        printf("step v failed to match\n");
//    }


//    std::unordered_set<std::bitset<N_PLUS_ONE>> solution_seqs_A;
//    std::unordered_set<std::bitset<N_PLUS_ONE>> solution_seqs_B;
//    std::unordered_set<std::bitset<N>> solution_seqs_C;
//    std::unordered_set<std::bitset<N>> solution_seqs_D;
//    std::bitset<N_PLUS_ONE> seq_plus_one_a;
//    std::bitset<N_PLUS_ONE> seq_plus_one_b;
//    std::bitset<N> seq_c;
//    std::bitset<N> seq_d;

//    seqPlusOneFinder(N_PLUS_ONE, 0, 2, seq_plus_one_a, solution_seqs_A);
//    std::bitset<N_PLUS_ONE> query("10010010011101111010");
//    if(solution_seqs_A.count(query)) {
//        printf("Found sequence in the sequence set");
//    }


//    int** quadruple = (int **)malloc(4 * sizeof(int *));
//    for (auto itr = result.begin(); itr != result.end(); itr++) {
//        i=0;
//        std::list<int>& seq = *itr;
//        for(int& e: seq) {
//            sums[i++] = e;
//        }
//
//        seqPlusOneFinder(N_PLUS_ONE, 0, sums[0], seq_plus_one_a, solution_seqs_A);
//        if(sums[0] != sums[1]) {
//            seqPlusOneFinder(N_PLUS_ONE, 0, sums[1], seq_plus_one_b, solution_seqs_B);
//        } else {
//            solution_seqs_B = solution_seqs_A;
//        }
//
//        seqFinder(N, 0, sums[2], seq_c, solution_seqs_C);
//        if(sums[2] != sums[3]) {
//            seqFinder(N, 0, sums[3], seq_d, solution_seqs_D);
//        } else {
//            solution_seqs_D = solution_seqs_C;
//        }
//        for (const auto& bitset_a : solution_seqs_A) {
//            for (const auto& bitset_b : solution_seqs_B) {
//                if(bitset_b!=bitset_a) {
//                    for (const auto& bitset_c : solution_seqs_C) {
//                        for (const auto& bitset_d : solution_seqs_D) {
//                            if(bitset_d!=bitset_c) {
//                                quadruple[0] = (int*)malloc(N_PLUS_ONE * sizeof(int));
//                                quadruple[1] = (int*)malloc(N_PLUS_ONE * sizeof(int));
//                                quadruple[2] = (int*)malloc(N * sizeof(int));
//                                quadruple[3] = (int*)malloc(N * sizeof(int));
//                                for(j=0; j<N_PLUS_ONE; j++) {
//                                    quadruple[0][j] = bitset_a[j]?1:-1;
//                                    quadruple[1][j] = bitset_b[j]?1:-1;
//                                }
//                                for(j=0; j<N; j++) {
//                                    quadruple[2][j] = bitset_c[j]?1:-1;
//                                    quadruple[3][j] = bitset_d[j]?1:-1;
//                                }
//                                if(turyn_quadruples_func(quadruple, N)){
//                                    print_bs(quadruple, N);
//                                }
//                                for(j=0; j<4; j++) {
//                                    free(quadruple[j]);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    free(quadruple);
    return 0;
}