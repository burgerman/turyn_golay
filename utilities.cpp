//
// Created by Wilfried Wu on 2024-05-29.
//

#include <cstdio>
#include <cctype>
#include <cstdlib>
#include <vector>
#include <cmath>
#include "utilities.h"

bool filter_array_element_by_mod(int arr[], int len, int mod, int arr_res[], char letter) {
    int i, j;
    int res;
    double benchmark;
    bool output = true;
    for(i=1; i<=mod; i++){
        benchmark = (len-i)/mod + 1;
        res = 0;
        printf("%c_{%d,%d}: ", letter, i, mod);
        if(i!=mod) {
            for(j=1;j<=len;j++) {
                if((j%mod)==i) {
                    res += arr[j-1];
                    printf("%d ", arr[j-1]);
                }
            }
        } else {
            for(j=1;j<=len;j++) {
                if((j%mod)==0) {
                    res += arr[j-1];
                    printf("%d ", arr[j-1]);
                }
            }
        }
        printf("\n");
        printf("The sum of %c_{%d,%d} is: %d\n", letter, i, mod, res);
        arr_res[i-1] = res;
        if(std::abs(res) > benchmark) {
            printf("%c_{%d,%d} is greater than the benchmark\n", letter, i, mod);
            output = false;
        }
    }
    return output;
}


std::vector<int> filter_array_element_by_mod(std::vector<int>& sequence, int len, int mod, char letter) {
    int i, j;
    int res;
    double benchmark;
    std::vector<int> sequence_res(mod);
    for(i=1; i<=mod; i++){
        benchmark = (len-i)/mod + 1;
        res = 0;

        if(i!=mod) {
            for(j=1;j<=len;j++) {
                if((j%mod)==i) {
                    res += sequence[j-1];

                }
            }
        } else {
            for(j=1;j<=len;j++) {
                if((j%mod)==0) {
                    res += sequence[j-1];

                }
            }
        }

//        printf("%c_{%d,%d} = %d\n", letter, i, mod, res);
        sequence_res[i-1] = res;
        if(std::abs(res) > benchmark) {
            printf("%c_{%d,%d} is greater than the benchmark\n", letter, i, mod);
        }
    }
    return sequence_res;
}

int naf_polynomial_decomposition (int s, int m, int arr_res[], char letter) {
    int i;
    int res = 0;
    printf("N_%c(%d) = ", toupper(letter), s);
    for(i=1;i<=m-s;i++) {
        if (i!=m-s) {
            printf("%c_{%d,%d} * ", letter,i,m);
            printf("%c_{%d,%d} + ", letter,i+s,m);
        }
        else {
            printf("%c_{%d,%d} * ", letter,i,m);
            printf("%c_{%d,%d}", letter,i+s,m);
        }
        res += arr_res[i-1]*arr_res[i-1+s];
    }
    printf("\n");
    return res;
}

int naf_polynomial_decomposition (int s, int m, std::vector<int>& sequence, char letter) {
    int i;
    int res = 0;
//    printf("N_%c(%d) = ", toupper(letter), s);
    for(i=1;i<=m-s;i++) {
//        if (i!=m-s) {
//            printf("%c_{%d,%d} * ", letter,i,m);
//            printf("%c_{%d,%d} + ", letter,i+s,m);
//        }
//        else {
//            printf("%c_{%d,%d} * ", letter,i,m);
//            printf("%c_{%d,%d}", letter,i+s,m);
//        }
        res += sequence[i-1]*sequence[i-1+s];
    }
//    printf("%d", res);
//    printf("\n");
    return res;
}

void print_sequence(int arr[], int len, char letter) {
    int i;
    printf("Sequence %c: [", letter);
    for (i=0; i<len;i++) {
        if(i!= len-1) {
            printf("%d, ", arr[i]);
        } else {
            printf("%d]\n", arr[i]);
        }
    }
}

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