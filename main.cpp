#include <iostream>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <vector>

using namespace std;

int Cmax(int* T, int* P, int* X, int M, int N);

bool sortbysecdesc(const pair<int, int>& a,
    const pair<int, int>& b)
{
    if (a.second == b.second)
        return a.first < b.first;
    else
        return a.second > b.second;
}

int main()
{
    // N - liczba zadan, M - liczba maszyn, P - czasy na poszczegolnych maszynach (grupowane M)
    // T - czasy zakonczenia poszczegolnych zadan, X - kolejnosc
    int N, M, * P, * T, * X;
    //int* W;
    vector<pair<int, int>> W;

    bool *ordered;
    int count;
    std::ifstream data("C:/Users/kacpe/source/repos/NEH/neh.data.txt");
    std::string find = "data.";
    find.append("063:");
    std::string tmp;
    while (tmp != find) {
        data >> tmp;
    }
    data >> N;
    data >> M;
    
    P = new int[M*N];
    T = new int[M+1];

    X = (int*)malloc(N * sizeof(int));
    int* temp = new int[N];

    //W = new int[N];
    W.resize(N);
    ordered = new bool[N];

    for (int i = 0; i < M * N; i++) {
        data >> P[i];
        //cout << P[i] << " ";
    }
    //cout << "\n\n";
    
    //algorytm
    auto start = chrono::steady_clock::now();
    //Wagi
    cout << "Wagi:" << endl;
    for (int i = 0; i < N; i++)
    {
        //W[i] = 0;
        W[i].second = 0; //waga
        W[i].first = i; //index
        ordered[i] = false;
        for (int j = 0; j < M; j++)
        {
            //W[i] += P[i * M + j];
            W[i].second += P[i * M + j];
            
        }
        //cout << W[i].second << " ";
    }
    sort(W.begin(), W.end(), sortbysecdesc);
    /*
    for (int j = 0; j < N; j++)
    {
        cout << "Id: " << W[j].first << " W= " << W[j].second << endl;
    }*/
    
    //cout << endl;
    int highest = 0;
    int index = -1;
    int* newX;

    
    for (int i = 0; i < N; i++)
    {
        
        newX = (int*)malloc((i + 1) * sizeof(int));
        //szukamy najwy¿szej wagi z nieuszeregowanych
        /*
        highest = 0;
        for (int j = 0; j < N; j++)
        {
            if (W[j].second > highest && !ordered[j])
            {
                highest = W[j].second;
                index = j;
            }
                
        }
        */
        //index zadania z najwiêksz¹ wag¹ z pozosta³ych
        index = i;// W[i].first;

        cout << "Id: " << W[index].first << " W= " << W[index].second << endl;
        ordered[index] = true;
        index = W[index].first;
        

        

        int bestIndex = 0;
        int c;
        int best = INT32_MAX;
        
        int skip = 0;

        for (int k = 0; k < i + 1; k++)//2
        {
            //Wstawiamy do uszeregowania
            newX[k] = index;
            skip = 0;
            for (int l = 0; l < i; l++)//
            {
                if (l == k)
                    skip = 1;
                newX[l + skip] = X[l];
            }
            /*cout << "newX: ";
            for (int z = 0; z < i+1; z++)
            {
                cout << newX[z] << " ";
            }
            cout << endl;*/
            //sprawdzamy czy jest lepiej
            c = Cmax(T, P, newX, M, i + 1);
            if (c < best)
            {
                best = c;
                bestIndex = k;
            }
            
            
        }

        //stare=========================================================================
        //for (int k = 0; k < i + 1; k++)//2
        //{
        //    //Wstawiamy do uszeregowania
        //    newX[k] = index;
        //    skip = 0;
        //    for (int l = 0; l < i; l++)//
        //    {
        //        if (l == k)
        //            skip = 1;
        //        newX[l + skip] = X[l];
        //    }
        //    /*cout << "newX: ";
        //    for (int z = 0; z < i+1; z++)
        //    {
        //        cout << newX[z] << " ";
        //    }
        //    cout << endl;*/
        //    //sprawdzamy czy jest lepiej
        //    c = Cmax(T, P, newX, M, i + 1);
        //    if (c < best)
        //    {
        //        best = c;
        //        bestIndex = k;
        //    }
        //    
        //    
        //}
        //======================================================================
        
        //wynik czêœciowy
        
        for (int j = 0; j < i+1; j++)
        {
            temp[j] = X[j];
        }
        X[bestIndex] = index;
        skip = 0;
        for (int j = 0; j < i; j++)//
        {
            if (j == bestIndex)
                skip = 1;
            X[j + skip] = temp[j];
        }

        free(newX);
        
    }


    /*
    // Kolejnosc: 1 4 3 2
    X[0] = 0;
    X[1] = 3;
    X[2] = 2;
    X[3] = 1;
    */
    auto stop = chrono::steady_clock::now();
    chrono::duration<double> elapsed = stop - start;
    cout << "Czas: " << elapsed.count() << endl;
    cout << "Kolejnosc: ";
    for (int i = 0; i < N; i++)
    {
        cout << X[i] + 1 << " ";
    }
    cout << endl << Cmax(T,P,X,M,N) << endl;
    delete[] P;
    delete[] T;
    free(X);
    delete[] temp;
    //delete[] W;
    W.clear();
    //delete[] ordered;
    

}


int Cmax(int *T, int *P, int *X, int M, int N)
{
    for (int i = 0; i <= M; i++) {
        T[i] = 0;
    }

    for (int n = 0; n < N; n++) {
        for (int m = 0; m < M; m++) {
            T[m + 1] = max(T[m + 1], T[m]) + P[(m)+X[n] * M];
        }
    }
    return T[M];

}
