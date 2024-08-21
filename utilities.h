//
// Created by Wilfried Wu on 2024-05-29.
//

#ifndef TURYN_GOLAY_UTILITIES_H
#define TURYN_GOLAY_UTILITIES_H
bool filter_array_element_by_mod(int arr[], int len, int mod, int arr_res[], char letter);
std::vector<int> filter_array_element_by_mod(std::vector<int>& sequence, int len, int mod, char letter);
int naf_polynomial_decomposition (int s, int m, int arr_res[], char letter);
int naf_polynomial_decomposition (int s, int m, std::vector<int>& sequence);
bool step2_cond3_k_r(int n, int m, int sequence1[], int sequence2[]);
void print_sequence(int arr[], int len, char letter);
void print_sequence(std::vector<int>& sequence, int len, char letter);
#endif //TURYN_GOLAY_UTILITIES_H
