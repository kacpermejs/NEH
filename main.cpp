#include <iostream>
#include <cstdlib>
#include <fstream>
#include <algorithm>


int main()
{
    // N - liczba zadan, M - liczba maszyn, P - czasy na poszczegolnych maszynach (grupowane M)
    // T - czasy zakonczenia poszczegolnych zadan, X - kolejnosc
    int N, M, * P, * T, * X;
    int count;
    std::ifstream data("C:/Users/Patryk/Desktop/Studia/SterowanieProcesamiDyskretnymi/Zad3/NEH/data.neh.txt");
    std::string find = "data.";
    find.append("000:");
    std::string tmp;
    while (tmp != find) {
        data >> tmp;
    }
    data >> N;
    data >> M;

    P = (int*)malloc(M * N * sizeof(int));
    T = (int*)malloc(M * sizeof(int));
    X = (int*)malloc(N * sizeof(int));

    for (int i = 0; i < M * N; i++) {
        data >> P[i];
        std::cout << P[i] << " ";
    }
    std::cout << "\n\n";

    // Kolejnosc: 1 4 3 2
    X[0] = 0;
    X[1] = 3;
    X[2] = 2;
    X[3] = 1;

    for (int i = 0; i < M; i++) {
        T[i] = 0;
    }

    for (int n = 0; n < N; n++) {
        for (int m = 0; m < M; m++) {
            T[m + 1] = std::max(T[m + 1], T[m]) + P[(m)+X[n] * M];
        }
    }
    std::cout << T[M] << std::endl;
}