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
            num_of_zeros = kk - seq.size();
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

bool condition2_element_verify_k_r(int n, int m, std::vector<int>& seq, std::vector<int>& seq2) {
    int i, current_j, max_v_j;
    bool isOdd_j;
    for(i=0; i<seq.size(); i++) {
        current_j = i+1;
        max_v_j = std::floor((n+1-current_j)/m) +1;
        isOdd_j = (max_v_j % 2 != 0);
        if(isOdd_j) {
            if((seq[i]%2==0) || (seq2[i]%2==0)) {
                return false;
            }
        } else {
            if((seq[i]%2!=0) || (seq2[i]%2!=0)) {
                return false;
            }
        }
    }
    return true;
}

bool condition2_element_verify_p_q(int n, int m, std::vector<int>& seq, std::vector<int>& seq2) {
    int i, current_j, max_v_j;
    bool isOdd_j;
    for(i=0; i<seq.size(); i++) {
        current_j = i+1;
        max_v_j = std::floor((n-current_j)/m)+1;
        isOdd_j = (max_v_j % 2 != 0);
        if(isOdd_j) {
            if((seq[i]%2==0) || (seq2[i]%2==0)) {
                return false;
            }
        } else {
            if((seq[i]%2!=0) || (seq2[i]%2!=0)) {
                return false;
            }
        }

    }
    return true;
}

std::vector<std::pair<int, int>> step3_pair_find_k_r(int n, int m, int element, int current_j) {
    std::vector<std::pair<int, int>> pairs;
    int max_v_j = std::floor((n+1-current_j)/m) +1;
    int max_v_j_plus_m = std::floor((n+1-current_j+m/2)/m) +1;
    int i, j;
    if(m==n+1) {
        for (i=-1; i<=1; i++) {
            for(j=-1; j<=1; j++) {
                if(i+j == element) {
                    pairs.emplace_back(i, j);
                }
            }
        }
        return pairs;
    }
    for (i=-max_v_j; i<=max_v_j; i++) {
        for(j=-max_v_j_plus_m; j<=max_v_j_plus_m; j++) {
            if(i+j == element) {
                pairs.emplace_back(i, j);
            }
        }
    }

//    if(!pairs.empty() && pairs.size()>1) {
//        std::sort(pairs.begin(), pairs.end());
//        auto last = std::unique(pairs.begin(), pairs.end());
//        pairs.erase(last, pairs.end());
//    }
    return pairs;
}

std::vector<std::pair<int, int>> step3_pair_find_p_q(int n, int m, int element, int current_j) {
    std::vector<std::pair<int, int>> pairs;
    int max_v_j = std::floor((n-current_j)/m) +1;
    int max_v_j_plus_m = std::floor((n-current_j+m/2)/m) +1;
    int i, j;
    if(m==n+1) {
        for (i=-1; i<=1; i++) {
            for(j=-1; j<=1; j++) {
                if(i+j == element) {
                    pairs.emplace_back(i, j);
                }
            }
        }
        return pairs;
    }
    for (i=-max_v_j; i<=max_v_j; i++) {
        for(j=-max_v_j_plus_m; j<=max_v_j_plus_m; j++) {
            if(i+j == element) {
                pairs.emplace_back(i, j);
            }
        }
    }
//    if(!pairs.empty() && pairs.size()>1) {
//        std::sort(pairs.begin(), pairs.end());
//        auto last = std::unique(pairs.begin(), pairs.end());
//        pairs.erase(last, pairs.end());
//    }
    return pairs;
}


void getSolutionToK11(int n, int m, int seq_11, int current_j, std::vector<int>& sequence_k, int currentSum, std::vector<std::vector<int>>& k_solutions) {
    if(current_j > m) {
        if(currentSum == seq_11) {
            k_solutions.push_back(sequence_k);
        }
        return;
    }
    int maxV = std::floor((n + 1 - current_j) / m) + 1;
    bool isOdd = (maxV % 2 != 0);
    int k;
    std::vector<std::pair<int, int>> pairs;
    if(m>=n+1) {
        for(k =-1; k<=1; k++) {
            if( (isOdd && (k%2 != 0)) || (!isOdd && (k%2 == 0))) {
                if(k+currentSum <= seq_11){
                    pairs = step3_pair_find_k_r(n, m+m, k, current_j);
                    for(std::pair<int, int> & p : pairs) {
                        sequence_k[current_j-1] = p.first;
                        sequence_k[current_j+m-1] = p.second;
                        getSolutionToK11(n, m, seq_11, current_j+1, sequence_k, currentSum + k, k_solutions);
                    }
                }
            }
        }
    } else {
        for(k = -maxV; k<=maxV; k++) {
            if( (isOdd && (k%2 != 0)) || (!isOdd && (k%2 == 0))) {
                if(k+currentSum <= seq_11){
                    pairs = step3_pair_find_k_r(n, m+m, k, current_j);
                    for(std::pair<int, int> & p : pairs) {
                        sequence_k[current_j-1] = p.first;
                        sequence_k[current_j+m-1] = p.second;
                        getSolutionToK11(n, m, seq_11, current_j+1, sequence_k, currentSum + k, k_solutions);
                    }
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
    int maxV = std::floor((n + 1 - current_j) / m) + 1;
    bool isOdd = (maxV % 2 != 0);
    int k;
    std::vector<std::pair<int, int>> pairs;
    if(m>=n+1) {
        for(k = -1; k<=1; k++) {
            if( (isOdd && (k%2 != 0)) || (!isOdd && (k%2 == 0))) {
                if(k+currentSum <= seq_11){
                    pairs = step3_pair_find_k_r(n, m+m, k, current_j);
                    for(std::pair<int, int> & p : pairs) {
                        sequence_r[current_j-1] = p.first;
                        sequence_r[current_j+m-1] = p.second;
                        getSolutionToR11(n, m, seq_11, current_j+1, sequence_r, currentSum + k, r_solutions);
                    }
                }
            }
        }
    }
    else {
        for(k = -maxV; k<=maxV; k++) {
            if( (isOdd && (k%2 != 0)) || (!isOdd && (k%2 == 0))) {
                if( k+currentSum <= seq_11){
                    pairs = step3_pair_find_k_r(n, m+m, k, current_j);
                    for(std::pair<int, int> & p : pairs) {
                        sequence_r[current_j-1] = p.first;
                        sequence_r[current_j+m-1] = p.second;
                        getSolutionToR11(n, m, seq_11, current_j+1, sequence_r, currentSum + k, r_solutions);
                    }
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
    int maxV = std::floor((n - current_j) / m) + 1;
    bool isOdd = (maxV % 2 != 0);
    int k;
    std::vector<std::pair<int, int>> pairs;
    for(k = -maxV; k<=maxV; k++) {
        if( (isOdd && (k%2 != 0)) || (!isOdd && (k%2 == 0))) {
            if(k+currentSum <= seq_11){
                pairs = step3_pair_find_p_q(n, m+m, k, current_j);
                for(std::pair<int, int> & p : pairs) {
                    sequence_p[current_j-1] = p.first;
                    sequence_p[current_j+m-1] = p.second;
                    getSolutionToP11(n, m, seq_11, current_j+1, sequence_p, currentSum + k, p_solutions);
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
    int maxV = std::floor((n - current_j) / m) + 1;
    bool isOdd = (maxV % 2 != 0);
    int k;
    std::vector<std::pair<int, int>> pairs;

    for(k = -maxV; k<=maxV; k++) {
        if( (isOdd && (k%2 != 0)) || (!isOdd && (k%2 == 0))) {
            if(k+currentSum <= seq_11){
                pairs = step3_pair_find_p_q(n, m+m, k, current_j);
                for(std::pair<int, int> & p : pairs) {
                    sequence_q[current_j-1] = p.first;
                    sequence_q[current_j+m-1] = p.second;
                    getSolutionToQ11(n, m, seq_11, current_j+1, sequence_q, currentSum + k, q_solutions);
                }
            }
        }
    }

}

bool condition3_k_r(int n, int m, std::vector<int>& sequence1, std::vector<int>& sequence2) {
    int j, index, sum;
    bool isMultipleOfM = (n%m == 0);

    for(j=2; j<=m/2; j++) {
        index = n+2-j;
        sum = sequence1[j-1]+sequence1[(j+m/2)-1]+sequence2[j-1]+sequence2[(j+m/2)-1]
              +sequence1[index-1]+sequence1[(index+m/2)-1]+sequence2[index-1]+sequence2[(index+m/2)-1];
//        printf("k_%d,%d + k_%d,%d + r_%d,%d + r_%d,%d + k_%d,%d + k_%d,%d + r_%d,%d + r_%d,%d= %d\n",j,m, (j+m/2),m, j,m,(j+m/2), m, index,m, (index+m/2),m, index,m, (index+m/2), m, sum);
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

bool condition3_p_q(int n, int m, std::vector<int>& sequence1, std::vector<int>& sequence2) {
    int j, index, sum;
    for(j=1; j<=m/2; j++) {
        index = n+1-j;
        sum = sequence1[j-1]+sequence1[(j+m/2)-1]+sequence2[j-1]+sequence2[(j+m/2)-1]
              +sequence1[index-1]+sequence1[(index+m/2)-1]+sequence2[index-1]+sequence2[(index+m/2)-1];
        if(sum%4!=0) {
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

bool step2_condition5 (int n, int m, std::vector<int>& sequenceK, std::vector<int>& sequenceR,
                       std::vector<int>& sequenceP, std::vector<int>& sequenceQ) {
    if(m>=n+1) {
        int j;
        int sum2 =0;
        for(j=1; j<m; j++) {
            sum2+= naf_polynomial_decomposition (j, m, sequenceK, 'k');
            sum2+= naf_polynomial_decomposition (j, m, sequenceR, 'r');
            sum2+= naf_polynomial_decomposition (j, m, sequenceP, 'p');
            sum2+= naf_polynomial_decomposition (j, m, sequenceQ, 'q');
            if (sum2 != 0) {
                return false;
            }
        }
        return true;
    }
    else {
        int i, sum;
        for (i=1; i<=m/2; i++) {
            sum = 0;
            sum += naf_polynomial_decomposition (i, m, sequenceK, 'k');
            sum += naf_polynomial_decomposition (m-i, m, sequenceK, 'k');
            sum += naf_polynomial_decomposition (i, m, sequenceR, 'r');
            sum += naf_polynomial_decomposition (m-i, m, sequenceR, 'r');
            sum += naf_polynomial_decomposition (i, m, sequenceP, 'p');
            sum += naf_polynomial_decomposition (m-i, m, sequenceP, 'p');
            sum += naf_polynomial_decomposition (i, m, sequenceQ, 'q');
            sum += naf_polynomial_decomposition (m-i, m, sequenceQ, 'q');
            if (sum != 0) {
                return false;
            }
        }
        return true;
    }
}

void step2_condition4 (int n, int m, std::vector<int>& filtered_k_solutions, std::vector<int>& filtered_r_solutions,
                       std::vector<int>& filtered_p_solutions, std::vector<int>& filtered_q_solutions,
                       std::vector<std::vector<int>>& k_solutions,  std::vector<std::vector<int>>& r_solutions,
                       std::vector<std::vector<int>>& p_solutions,  std::vector<std::vector<int>>& q_solutions) {

    int i, j, k, l, result;
    int cond = 4*n+2;
    for (i=0; i<filtered_k_solutions.size(); i++) {
        for(j=0; j<filtered_r_solutions.size(); j++) {
            for(k=0; k<filtered_p_solutions.size(); k++) {
                for(l=0; l<filtered_q_solutions.size(); l++) {
                    result = getSquaredSum(m, k_solutions[filtered_k_solutions[i]]) +
                             getSquaredSum(m, r_solutions[filtered_r_solutions[j]]) +
                             getSquaredSum(m, p_solutions[filtered_p_solutions[k]]) +
                             getSquaredSum(m, q_solutions[filtered_q_solutions[l]]);
                    if (result == cond) {
                        bool cond5 = step2_condition5(n, m, k_solutions[filtered_k_solutions[i]],
                                                      r_solutions[filtered_r_solutions[j]],
                                                      p_solutions[filtered_p_solutions[k]],
                                                      q_solutions[filtered_q_solutions[l]]);
                        if (cond5) {
                            printf("2m = %d\n", m);
                            print_sequence(k_solutions[filtered_k_solutions[i]], m, 'K');
                            print_sequence(r_solutions[filtered_r_solutions[j]], m, 'R');
                            print_sequence(p_solutions[filtered_p_solutions[k]], m, 'P');
                            print_sequence(q_solutions[filtered_q_solutions[l]], m, 'Q');
                            printf("\n");
                        }
                    }
                }
            }
        }
    }


}

void getAllSolutions(int n, int m, int k11, int r11, int p11, int q11) {
    if(m>(n+1)/2) return;
    std::vector<int> sequence_k(2*m);
    std::vector<int> sequence_r(2*m);
    std::vector<int> sequence_p(2*m);
    std::vector<int> sequence_q(2*m);
    std::vector<std::vector<int>> k_solutions;
    std::vector<std::vector<int>> r_solutions;
    std::vector<std::vector<int>> p_solutions;
    std::vector<std::vector<int>> q_solutions;

    std::vector<int> filtered_k_solutions;
    std::vector<int> filtered_r_solutions;
    std::vector<int> filtered_p_solutions;
    std::vector<int> filtered_q_solutions;

    getSolutionToK11(n, m, k11, 1, sequence_k, 0, k_solutions);
    getSolutionToR11(n, m, r11, 1, sequence_r, 0, r_solutions);
    getSolutionToP11(n, m, p11, 1, sequence_p, 0, p_solutions);
    getSolutionToQ11(n, m, q11, 1, sequence_q, 0, q_solutions);
    int i, j;
    bool cond2_k_r, cond2_p_q, cond3_k_r, cond3_p_q;
    for(i=0; i<k_solutions.size(); i++) {
        for(j=0; j<r_solutions.size(); j++) {
            cond2_k_r = condition2_element_verify_k_r(n, m+m, k_solutions[i], r_solutions[j]);
            cond3_k_r = condition3_k_r(n,m+m, k_solutions[i], r_solutions[j]);
            if(cond2_k_r && cond3_k_r) {
                filtered_k_solutions.push_back(i);
                filtered_r_solutions.push_back(j);
            }
        }
    }
    for(i=0; i<p_solutions.size(); i++) {
        for(j=0; j<q_solutions.size(); j++) {
            cond2_p_q = condition2_element_verify_p_q(n, m+m, p_solutions[i], q_solutions[j]);
            cond3_p_q = condition3_p_q(n,m+m, p_solutions[i], q_solutions[j]);
            if(cond2_p_q && cond3_p_q) {
                filtered_p_solutions.push_back(i);
                filtered_q_solutions.push_back(j);
            }
        }
    }
    removeDuplicates(filtered_k_solutions);
    removeDuplicates(filtered_r_solutions);
    removeDuplicates(filtered_p_solutions);
    removeDuplicates(filtered_q_solutions);
    printf("K solutions found: %d \n", filtered_k_solutions.size());
    printf("R solutions found: %d \n", filtered_r_solutions.size());
    printf("P solutions found: %d \n", filtered_p_solutions.size());
    printf("Q solutions found: %d \n", filtered_q_solutions.size());
    step2_condition4(n, m+m, filtered_k_solutions, filtered_r_solutions, filtered_p_solutions, filtered_q_solutions,
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
        printf("number of evens: %d\n", evens.size());
        printf("number of odds : %d\n", odds.size());
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
        for(j=2; j<=m; j++) {
            getAllSolutions(n, j, k11, r11, p11, q11);
            printf("--------------------------------------------------------------------------\n");
            printf("\n");
        }
    }
}

int main(int argc, char *argv[]) {
    int n = 5;
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
    int m = n+1;
    findQuadruple(n, m);
    return 0;
}