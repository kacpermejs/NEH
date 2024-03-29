#include <iostream>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <vector>

using namespace std;

int Cmax(int* T, int* P, int* X, int M, int N);
int Cmax1(int* T, int* P, const vector<int> &X, int M, int N)
{
    //zerowanie 
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
int Cmax2(vector<vector<int>>& Tr, vector<vector<int>>& Tl, int* P, const vector<int> &X, int M, int N, int Elem, int Place)
{

    P[Elem];//nowy element
    vector<int> R;
    vector<int> L;
    R.resize(M);
    L.resize(M);
    for (int i = 0; i < M; i++)
    {
        if (Place == 0 && i == 0)
        {
            R[i] = P[Elem * M + i];
        }
        else if (Place == 0 && i > 0)
        {
            R[i] = R[i - 1] + P[Elem * M + i];
        }
        else if (Place > 0 && i == 0)
        {
            R[i] = Tr[i][Place-1] + P[Elem * M + i];
        }
        else
        {
            R[i] = std::max(R[i - 1] + P[Elem * M + i], Tr[i][Place-1] + P[Elem * M + i]);
        }
    }

    
    if (Place == N - 1)
        return R[M - 1];

    int highest = 0;
    for (int i = 0; i < M; i++)
    {
        
        highest = max(R[i] + Tl[i][Place + 1], highest);
    }
    return highest;
}

bool sortbysecdesc(const pair<int, int>& a,
    const pair<int, int>& b)
{
    if (a.second == b.second)
        return a.first < b.first;
    else
        return a.second > b.second;
}

void algorytm0(int* T, int* P, int* X, int M, int N,  bool* ordered, vector<pair<int, int>> W, int* temp)
{
    //Wagi
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
        //szukamy najwy�szej wagi z nieuszeregowanych
        
        highest = 0;
        for (int j = 0; j < N; j++)
        {
            if (W[j].second > highest && !ordered[j])
            {
                highest = W[j].second;
                index = j;
            }

        }
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

            //sprawdzamy czy jest lepiej
            c = Cmax(T, P, newX, M, i + 1);
            if (c < best)
            {
                best = c;
                bestIndex = k;
            }


        }

        //wynik cz�ciowy

        for (int j = 0; j < i + 1; j++)
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
}

void algorytm1(int* T, int* P, vector<int> &X, int M, int N, bool* ordered, vector<pair<int, int>> W, int* temp)//poprawki
{

    //Wagi
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
    

    //cout << endl;
    int highest = 0;
    int index = -1;
    //int* newX;
    //newX = (int*)malloc(N * sizeof(int));
    //wersja z std::vector
    

    for (int i = 0; i < N; i++)
    {

        
        //index zadania z najwi�ksz� wag� z pozosta�ych
        index = i;// W[i].first;

        //cout << "Id: " << W[index].first << " W= " << W[index].second << endl;
        ordered[index] = true;
        index = W[index].first;

        int bestIndex = 0;
        int c;
        int best = INT32_MAX;

        for (int k = 0; k < i + 1; k++)//2
        {
            //wstawiamy na nowe miejsce
            X.insert(X.begin() + k, index);

            //sprawdzamy czy jest lepiej
            c = Cmax1(T, P, X, M, i + 1);
            if (c < best)
            {
                best = c;
                bestIndex = k;
            }
            X.erase(X.begin() + k);



        }

        //wynik cz�ciowy
        X.insert(X.begin() + bestIndex, index);

    }
    
}

void algorytm2(vector<vector<int>> Tr, vector<vector<int>> Tl, int* P, vector<int> &X, int M, int N, vector<pair<int, int>> W, int* temp)//qneh
{
    //Wagi
    for (int i = 0; i < N; i++)
    {
        //W[i] = 0;
        W[i].second = 0; //waga
        W[i].first = i; //index
        
        for (int j = 0; j < M; j++)
        {
            //W[i] += P[i * M + j];
            W[i].second += P[i * M + j];

        }
        //cout << W[i].second << " ";
    }
    sort(W.begin(), W.end(), sortbysecdesc);


    //cout << endl;
    int highest = 0;
    int index = -1;
    


    for (int l = 0; l < N; l++)
    {

        //index zadania z najwi�ksz� wag� z pozosta�ych
        index = l;// W[i].first;

        
        index = W[index].first;

        int bestIndex = 0;
        int c;
        int best = INT32_MAX;

        for (int k = 0; k < l + 1; k++)//2
        {
            //wstawiamy na nowe miejsce
            //X.insert(X.begin() + k, index);
            
            //sprawdzamy Cmax na k-tym miejscu

            //sprawdzamy czy jest lepiej
            c = Cmax2(Tr, Tl, P, X, M, l + 1, index, k);
            if (c < best)
            {
                best = c;
                bestIndex = k;
            }
            //X.erase(X.begin() + k);



        }

        //wynik cz�ciowy
        X.insert(X.begin() + bestIndex, index);
        vector<int> R;
        R.resize(M);
        
        for (int i = 0; i < M; i++)
        {
            Tr[i].insert(Tr.begin() + bestIndex, 0);
            Tl[i].insert(Tl.begin() + bestIndex, 0);
        }
        
        
        ///TODO naprawianie grafu
        for (int i = bestIndex; i < l+1; i++)//po zadaniach
        {
            for (int j = 0; j < M; j++)//po maszynach
            {
                if (i == 0 && j == 0)
                {
                    Tr[j][i] = P[X[i] * M + j];
                }
                else if (j > 0 && i == 0)
                {
                    Tr[j][i] = Tr[j - 1][i] + P[X[i] * M + j];//[j-1][i]
                }
                else if (i > 0 && j == 0)
                {
                    Tr[j][i] = Tr[j][i - 1] + P[X[i] * M + j];//[j][i-1]
                }
                else
                {
                    Tr[j][i] = std::max(Tr[j - 1][i] + P[X[i] * M + j], Tr[j][i - 1] + P[X[i] * M + j]);
                }
            }
        }


        for (int i = bestIndex; i >= 0; i--)
        {
            for (int j = M - 1; j >= 0; j--)
            {
                if (i == N - 1 && j == M - 1)
                {
                    Tl[j][i] = P[X[i] * M + j];
                }
                else if (j < M - 1 && i == N - 1)
                {
                    Tl[j][i] = Tl[j + 1][i] + P[X[i] * M + j];
                }
                else if (i < N - 1 && j == M - 1)
                {
                    Tl[j][i] = Tl[j][i + 1] + P[X[i] * M + j];
                }
                else
                {
                    Tl[j][i] = std::max(Tl[j + 1][i] + P[X[i] * M + j], Tl[j][i + 1] + P[X[i] * M + j]);
                }
            }
        }

    }
}

/*


*/

int main()
{
    // N - liczba zadan,
    // M - liczba maszyn, 
    // P - czasy na poszczegolnych maszynach (grupowane M)
    // T - czasy zakonczenia dla poszczeg�lnych maszyn (nadpisywana) , 
    // X - kolejnosc
    int N, M, * P, * T;
    int* Xold;
    //int* W;
    vector<pair<int, int>> W;
    vector<int> X;

    

    bool *ordered;
    int count;
    std::ifstream data("C:/Users/kacpe/source/repos/NEH/neh.data.txt");
    std::string find = "data.";
    find.append("000:");
    std::string tmp;
    while (tmp != find) {
        data >> tmp;
    }
    data >> N;
    data >> M;

    int* temp = new int[N];

    P = new int[M*N];
    T = new int[M+1];
    
    Xold = (int*)malloc(N * sizeof(int));
    X.resize(N);

    //qneh
    vector<vector<int>> TGraphR;
    TGraphR.resize(M);
    for (int i = 0; i < M; i++)
    {
        TGraphR[i].resize(N);
    }
    //int** TGraphR = new int* [N];
    ////Allocating the column space in heap dynamically
    //for (int i = 0; i < N; i++) {
    //    TGraphR[i] = new int[M];
    //}
    //qneh
    vector<vector<int>> TGraphL;
    TGraphL.resize(M);
    for (int i = 0; i < M; i++)
    {
        TGraphL[i].resize(N);
    }
    /*
    int** TGraphL = new int* [N];
    //Allocating the column space in heap dynamically
    for (int i = 0; i < N; i++) {
        TGraphL[i] = new int[M];
    }*/

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
    //algorytm2(T, P, Xold, M, N, ordered, W, temp);
    //algorytm1(T, P, X, M, N, ordered, W, temp);
    algorytm2(TGraphR, TGraphL, P, X, M, N, W, temp);



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
    //cout << endl << Cmax1(T,P,X,M,N) << endl;
    //cout << endl << Cmax(T,P,Xold,M,N) << endl;
    cout << endl << TGraphL[0][0] << endl;
    delete[] P;
    delete[] T;
    free(Xold);
    X.clear();
    delete[] temp;
    //delete[] W;
    W.clear();
    delete[] ordered;
    //delete[] TGraphR;
    //delete[] TGraphL;

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
