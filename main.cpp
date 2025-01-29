#include <iostream>
#include <cmath>
#include <list>
#include <bitset>
#include <unordered_set>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <climits>
#include <algorithm>
#include <omp.h>
#include <mpi.h>
#include <numeric>
#include <set>


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
const int THRESHOLD = 600;
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

int naf_polynomial_decomposition (int s, int m, int* sequence) {
    int i;
    int res = 0;
    for(i=1;i<=m-s;i++) {
        res += sequence[i-1]*sequence[i-1+s];
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

void find_k11(int n, int m, int k11, int current_j, std::vector<int>& sequence_k, int currentSum, std::vector<std::vector<int> >& k_solutions) {
    if(k_solutions.size()>THRESHOLD) return;
    if(current_j > m) {
        if(currentSum == k11) {
            k_solutions.push_back(sequence_k);
        }
        return;
    }
    sequence_k[current_j-1] = 1;
    find_k11(n,m,k11, current_j+1, sequence_k, currentSum+1, k_solutions);
    sequence_k[current_j-1] = -1;
    find_k11(n,m,k11, current_j+1, sequence_k, currentSum-1, k_solutions);
}

void find_r11(int n, int m, int r11, int current_j, std::vector<int>& sequence_r, int currentSum, std::vector<std::vector<int> >& r_solutions) {
    if(r_solutions.size()>THRESHOLD) return;
    if(current_j > m) {
        if(currentSum == r11) {
            r_solutions.push_back(sequence_r);
        }
        return;
    }
    sequence_r[current_j-1] = 1;
    find_r11(n,m,r11, current_j+1, sequence_r, currentSum+1, r_solutions);
    sequence_r[current_j-1] = -1;
    find_r11(n,m,r11, current_j+1, sequence_r, currentSum-1, r_solutions);
}

void find_p11(int n, int m, int p11, int current_j, std::vector<int>& sequence_p, int currentSum, std::vector<std::vector<int> >& p_solutions) {
    if(p_solutions.size()>THRESHOLD) return;
    if(current_j == m) {
        if(currentSum == p11) {
            sequence_p[current_j-1] =0;
            p_solutions.push_back(sequence_p);
        }
        return;
    }
    sequence_p[current_j-1] = 1;
    find_p11(n,m,p11, current_j+1, sequence_p, currentSum+1, p_solutions);
    sequence_p[current_j-1] = -1;
    find_p11(n,m,p11, current_j+1, sequence_p, currentSum-1, p_solutions);
}

void find_q11(int n, int m, int q11, int current_j, std::vector<int>& sequence_q, int currentSum, std::vector<std::vector<int> >& q_solutions) {
    if(q_solutions.size()>THRESHOLD) return;
    if(current_j == m) {
        if(currentSum == q11) {
            sequence_q[current_j-1] =0;
            q_solutions.push_back(sequence_q);
        }
        return;
    }
    sequence_q[current_j-1] = 1;
    find_q11(n,m,q11, current_j+1, sequence_q, currentSum+1, q_solutions);
    sequence_q[current_j-1] = -1;
    find_q11(n,m,q11, current_j+1, sequence_q, currentSum-1, q_solutions);
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

void get_subarr(int start, int m, std::vector<int>& original_arr, int* sub_arr) {
    int i, j;
    int end = start+m;
    j=0;
    for (i=start; i<end; i++) {
        sub_arr[j] = original_arr[i];
        j++;
    }
}

int getSquaredSum(int m, int* sequence) {
    int sum =0;
    int i;
    for(i=0; i<m; i++) {
        sum += sequence[i] * sequence[i];
    }
    return sum;
}

//int getSquaredSum(int m, std::vector<int>& sequence) {
//     Ensure m is not larger than the sequence size
//    m = std::min(m, static_cast<int>(sequence.size()));
//
//    return std::accumulate(sequence.begin(), sequence.begin() + m, 0,
//                           [](int sum, int value) {
//                               return sum + value * value;
//                           });
//}

void print_sequence(int* sequence, int len, char letter) {
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

int step2_condition5 (int m, int* sequenceK, int* sequenceR, int* sequenceP, int* sequenceQ) {
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
            return 0;
        }
    }
    return 1;
}

void step2_condition4 (int n, int m, int &res, std::vector<int>& k11_solutions, std::vector<int>& r11_solutions, std::vector<int>& p11_solutions, std::vector<int>& q11_solutions, int** seq_k_r_p_q) {
    int result;
    int t;
    int cond = 4*n+2;
    int k11_size = (int)k11_solutions.size();
    int r11_size = (int)r11_solutions.size();
    int p11_size = (int)p11_solutions.size();
    int q11_size = (int)q11_solutions.size();

    for (int i=0; i<k11_size; i+=m) {
//        get_subarr(i, m, k11_solutions, seq_k_r_p_q[0]);
        for (t=0; t<m; t++) {
            seq_k_r_p_q[0][t] = k11_solutions[i+t];
        }

        printf("sequenceK: %d\n", seq_k_r_p_q[0][0]);
        for (int j=0; j<r11_size; j+=m) {
//            get_subarr(j, m, r11_solutions, seq_k_r_p_q[1]);
            for (t=0; t<m; t++) {
                seq_k_r_p_q[1][t] = r11_solutions[j+t];
            }
            printf("sequenceR: %d\n", seq_k_r_p_q[1][0]);
            for(int k=0;k<p11_size; k+=m){
//                get_subarr(k, m, p11_solutions, seq_k_r_p_q[2]);
                for (t=0; t<m; t++) {
                    seq_k_r_p_q[2][t] = r11_solutions[k+t];
                }
                printf("sequenceP: %d\n", seq_k_r_p_q[2][0]);
                for(int l=0;l<q11_size;l+=m) {
//                    get_subarr(l, m, q11_solutions, seq_k_r_p_q[3]);
                    for (t=0; t<m; t++) {
                        seq_k_r_p_q[3][t] = r11_solutions[l+t];
                    }
                    printf("sequenceQ: %d\n", seq_k_r_p_q[3][0]);
                    result = getSquaredSum(m, seq_k_r_p_q[0]) +
                             getSquaredSum(m, seq_k_r_p_q[1]) +
                             getSquaredSum(m, seq_k_r_p_q[2]) +
                             getSquaredSum(m, seq_k_r_p_q[3]);
                    if (result == cond) {
                        int cond5 = step2_condition5(m, seq_k_r_p_q[0],seq_k_r_p_q[1],seq_k_r_p_q[2],seq_k_r_p_q[3]);
                        if (cond5==1) {
                            print_sequence(seq_k_r_p_q[0], m, 'K');
                            print_sequence(seq_k_r_p_q[1], m, 'R');
                            print_sequence(seq_k_r_p_q[2], m, 'P');
                            print_sequence(seq_k_r_p_q[3], m, 'Q');
                            printf("\n");
                            res = 1;
                            return;
                        }
                    }
                }
            }
        }
    }
}

void flatten_solutions(std::vector<int>& filtered_indices, std::vector<std::vector<int>>& solutions, std::vector<int>& sent_solutions) {
    int i, j;
    for (i=0; i<filtered_indices.size(); i++) {
        auto seq = solutions[filtered_indices[i]];
        for(j=0; j<seq.size(); j++) {
            sent_solutions.push_back(seq[j]);
        }
    }
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int i, j, k, l, t, c, k11, r11, p11, q11, my_rank, p, p_rank;
    int sum_result, naf_sum, cond5, npaf_sum, step4_cond;
    int s1, s2;
    int ** seq_k_r_p_q;
    int res, global_found;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    if (p < 2) {
        printf("This program requires at least 2 processes.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int n = 19;
    int m = n+1;
    printf("n=%d\n", n);
    int quadruple_sum = 4*n+2;
    printf("4*n+2 = %d\n", quadruple_sum);
    int sums[4];
    // check 4n+2
    std::vector<int> received_k11s;
    std::vector<int> received_r11s;
    std::vector<int> received_p11s;
    std::vector<int> received_q11s;
    int k11_size, r11_size, p11_size, q11_size;
    std::list< std::list<int>> result = nsoks(quadruple_sum, 4);
    for(auto itr = result.begin(); itr != result.end(); itr++) {
        res = 0;
        global_found = 0;

        if(my_rank == 0) {
            j=0;
            std::vector<int> evens;
            std::vector<int> odds;
            std::list<int>& seq = *itr;
            printf("NSOKS: ");
            for(int& e: seq) {
                sums[j++] = e;
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
        }
        if(my_rank==0) {
            std::vector<int> send_solutions;
            printf("current m=%d \n", m);
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

            find_k11(n, m, k11, 1, sequence_k, 0, k_solutions);
            if(k11!=r11) {
                find_r11(n, m, r11, 1, sequence_r, 0, r_solutions);
            } else {
                r_solutions = k_solutions;
            }
            find_p11(n, m, p11, 1, sequence_p, 0, p_solutions);
            if(p11!=q11) {
                find_q11(n, m, q11, 1, sequence_q, 0, q_solutions);
            }
            else {
                q_solutions = p_solutions;
            }

            bool cond_k_r, cond_p_q;
            for(l=0; l<k_solutions.size(); l++) {
                for(i=0; i<r_solutions.size(); i++) {
                    cond_k_r = step2Condition3_k_r(n,m, k_solutions[l], r_solutions[i]);
                    if(cond_k_r) {
                        filtered_k_solutions.push_back(l);
                        filtered_r_solutions.push_back(i);
                    }
                }
            }
            for(l=0; l<p_solutions.size(); l++) {
                for(i=0; i<q_solutions.size(); i++) {
                    cond_p_q = step2Condition3_p_q(n,m, p_solutions[l], q_solutions[i]);
                    if(cond_p_q) {
                        filtered_p_solutions.push_back(l);
                        filtered_q_solutions.push_back(i);
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
            printf("k11 solutions found: %zu \n", filtered_k_solutions.size());
            printf("r11 solutions found: %zu \n", filtered_r_solutions.size());
            printf("p11 solutions found: %zu \n", filtered_p_solutions.size());
            printf("q11 solutions found: %zu \n", filtered_q_solutions.size());


            k11_size = (int)filtered_k_solutions.size()*m;
            r11_size = (int)filtered_r_solutions.size()*m;
            p11_size = (int)filtered_p_solutions.size()*m;
            q11_size = (int)filtered_q_solutions.size()*m;
            std::vector<int> sent_k11s;
            std::vector<int> sent_r11s;
            std::vector<int> sent_p11s;
            std::vector<int> sent_q11s;
            int sent_q_size = (int)filtered_q_solutions.size()/(p-1);
            sent_q_size = sent_q_size*m;
            flatten_solutions(filtered_k_solutions, k_solutions, sent_k11s);
            flatten_solutions(filtered_r_solutions, r_solutions, sent_r11s);
            flatten_solutions(filtered_p_solutions, p_solutions, sent_p11s);
            flatten_solutions(filtered_q_solutions, q_solutions, sent_q11s);
            int offset;
            int residual_size;
            for (p_rank=1; p_rank<p;p_rank++) {
                MPI_Send(&k11_size, 1, MPI_INT, p_rank, 0, MPI_COMM_WORLD);
                MPI_Send(&r11_size, 1, MPI_INT, p_rank, 0, MPI_COMM_WORLD);
                MPI_Send(&p11_size, 1, MPI_INT, p_rank, 0, MPI_COMM_WORLD);
                MPI_Send(sent_k11s.data(), k11_size, MPI_INT, p_rank, 0, MPI_COMM_WORLD);
                MPI_Send(sent_r11s.data(), r11_size, MPI_INT, p_rank, 0, MPI_COMM_WORLD);
                MPI_Send(sent_p11s.data(), p11_size, MPI_INT, p_rank, 0, MPI_COMM_WORLD);
                offset = p_rank-1;
                if(p_rank!=p-1) {
                    MPI_Send(&sent_q_size, 1, MPI_INT, p_rank, 0, MPI_COMM_WORLD);
                    MPI_Send(sent_q11s.data()+offset*sent_q_size, sent_q_size, MPI_INT, p_rank, 0, MPI_COMM_WORLD);
                } else {
                    residual_size = q11_size - offset*sent_q_size;
                    printf("residual_size: %d \n", residual_size);
                    MPI_Send(&residual_size, 1, MPI_INT, p_rank, 0, MPI_COMM_WORLD);
                    MPI_Send(sent_q11s.data()+offset*sent_q_size, residual_size, MPI_INT, p_rank, 0, MPI_COMM_WORLD);
                }
            }
            std::vector<std::vector<int>>().swap(k_solutions);
            std::vector<std::vector<int>>().swap(r_solutions);
            std::vector<std::vector<int>>().swap(p_solutions);
            std::vector<std::vector<int>>().swap(q_solutions);
            std::vector<int>().swap(filtered_k_solutions);
            std::vector<int>().swap(filtered_r_solutions);
            std::vector<int>().swap(filtered_p_solutions);
            std::vector<int>().swap(filtered_q_solutions);
            std::vector<int>().swap(sent_k11s);
            std::vector<int>().swap(sent_r11s);
            std::vector<int>().swap(sent_p11s);
            std::vector<int>().swap(sent_q11s);
        } else {
            MPI_Recv(&k11_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("k11 solutions received: %d \n", k11_size);
            MPI_Recv(&r11_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("r11 solutions received: %d \n", r11_size);
            MPI_Recv(&p11_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("p11 solutions received: %d \n", p11_size);
            received_k11s.reserve(k11_size);
            received_r11s.reserve(r11_size);
            received_p11s.reserve(p11_size);
            MPI_Recv(received_k11s.data(), k11_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(received_r11s.data(), r11_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(received_p11s.data(), p11_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&q11_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("q11 solutions received: %d \n", q11_size);
            received_q11s.reserve(q11_size);
            MPI_Recv(received_q11s.data(), q11_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            seq_k_r_p_q = (int **)malloc(4 * sizeof(int *));
            for(t=0; t<4; t++) {
                seq_k_r_p_q[t] = (int *)malloc(m * sizeof(int));
            }
            for (i=0; i<k11_size; i+=m) {
                if(res==1){
                    break;
                }
                for (t=0; t<m; t++) {
                    seq_k_r_p_q[0][t] = received_k11s[i+t];
                }

                for (j=0; j<r11_size; j+=m) {
                    if(res==1){
                        break;
                    }
                    for (t=0; t<m; t++) {
                        seq_k_r_p_q[1][t] = received_r11s[j+t];
                    }

                    for(k=0;k<p11_size; k+=m){
                        if(res==1){
                            break;
                        }
                        for (t=0; t<m; t++) {
                            seq_k_r_p_q[2][t] = received_p11s[k+t];
                        }

                        for(l=0;l<q11_size;l+=m) {

                            for (t=0; t<m; t++) {
                                seq_k_r_p_q[3][t] = received_q11s[l+t];
                            }

                            step4_cond = 1;
                            for(t=1; t<m; t++) {
                                s1 = t;
                                npaf_sum = 0;
                                for(c=1;c<=m-s1;c++) {
                                    npaf_sum += seq_k_r_p_q[0][c-1]*seq_k_r_p_q[0][c-1+s1];
                                    npaf_sum += seq_k_r_p_q[1][c-1]*seq_k_r_p_q[1][c-1+s1];
                                    npaf_sum += seq_k_r_p_q[2][c-1]*seq_k_r_p_q[2][c-1+s1];
                                    npaf_sum += seq_k_r_p_q[3][c-1]*seq_k_r_p_q[3][c-1+s1];
                                }
                                if(npaf_sum!=0) {
                                    step4_cond =0;
                                    break;
                                }
                            }
                            if(step4_cond==1) {
                                //condition 4
                                sum_result = 0;
                                for(t=0; t<m; t++) {
                                    sum_result += seq_k_r_p_q[0][t] * seq_k_r_p_q[0][t];
                                    sum_result += seq_k_r_p_q[1][t] * seq_k_r_p_q[1][t];
                                    sum_result += seq_k_r_p_q[2][t] * seq_k_r_p_q[2][t];
                                    sum_result += seq_k_r_p_q[3][t] * seq_k_r_p_q[3][t];
                                }

                                if (sum_result == quadruple_sum) {
                                    cond5=1;
                                    for (t=1; t<=m/2; t++) {
                                        naf_sum = 0;
                                        s1 = t;
                                        s2 = m-t;
                                        for(c=1;c<=m-s1;c++) {
                                            naf_sum += seq_k_r_p_q[0][c-1]*seq_k_r_p_q[0][c-1+s1];
                                            naf_sum += seq_k_r_p_q[1][c-1]*seq_k_r_p_q[1][c-1+s1];
                                            naf_sum += seq_k_r_p_q[2][c-1]*seq_k_r_p_q[2][c-1+s1];
                                            naf_sum += seq_k_r_p_q[3][c-1]*seq_k_r_p_q[3][c-1+s1];
                                        }
                                        for(c=1;c<=m-s2;c++) {
                                            naf_sum += seq_k_r_p_q[0][c-1]*seq_k_r_p_q[0][c-1+s2];
                                            naf_sum += seq_k_r_p_q[1][c-1]*seq_k_r_p_q[1][c-1+s2];
                                            naf_sum += seq_k_r_p_q[2][c-1]*seq_k_r_p_q[2][c-1+s2];
                                            naf_sum += seq_k_r_p_q[3][c-1]*seq_k_r_p_q[3][c-1+s2];
                                        }
                                        if (naf_sum != 0) {
                                            cond5 = 0;
                                            break;
                                        }
                                    }
                                    if (cond5==1) {
                                        if(step4_cond==1) {
                                            printf("process %d found a solution as below:\n", my_rank);
                                            printf("Sequence k11: [");
                                            for (t=0; t<m;t++) {
                                                if(t!= m-1) {
                                                    printf("%d, ", seq_k_r_p_q[0][t]);
                                                } else {
                                                    printf("%d]\n", seq_k_r_p_q[0][t]);
                                                }
                                            }
                                            printf("Sequence r11: [");
                                            for (t=0; t<m;t++) {
                                                if(t!= m-1) {
                                                    printf("%d, ", seq_k_r_p_q[1][t]);
                                                } else {
                                                    printf("%d]\n", seq_k_r_p_q[1][t]);
                                                }
                                            }
                                            printf("Sequence p11: [");
                                            for (t=0; t<m;t++) {
                                                if(t!= m-1) {
                                                    printf("%d, ", seq_k_r_p_q[2][t]);
                                                } else {
                                                    printf("%d]\n", seq_k_r_p_q[2][t]);
                                                }
                                            }
                                            printf("Sequence q11: [");
                                            for (t=0; t<m;t++) {
                                                if(t!= m-1) {
                                                    printf("%d, ", seq_k_r_p_q[3][t]);
                                                } else {
                                                    printf("%d]\n", seq_k_r_p_q[3][t]);
                                                }
                                            }
                                            printf("\n");
                                            res = 1;
                                            break;
                                        }
                                    } else {
                                        continue;
                                    }
                                } else{
                                    continue;
                                }
                            }
                        }
                    }
                }
            }
            for(t=0; t<4; t++) {
                free(seq_k_r_p_q[t]);
            }
            free(seq_k_r_p_q);
        }
        if(my_rank!=0) {
            MPI_Send(&res, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        } else {
            for(p_rank=1; p_rank<p;p_rank++) {
                MPI_Recv(&global_found, 1, MPI_INT, p_rank, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                res |= global_found;
            }
            global_found = res;
        }

        if(my_rank!=0) {
            MPI_Recv(&global_found, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        } else {
            for(p_rank=1; p_rank<p;p_rank++) {
                MPI_Send(&res, 1, MPI_INT, p_rank, 0, MPI_COMM_WORLD);
            }
        }
        if(global_found==1) break;
    }
    MPI_Finalize();
    return 0;
}