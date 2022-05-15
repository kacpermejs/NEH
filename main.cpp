#include <iostream>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <vector>

using namespace std;

struct Job
{
    int index;//numer zadania
    Job* Previous;
    Job* Next;
    vector<int>  R;
    vector<int>  L;
    

    Job(int m, int i): Previous(nullptr), Next(nullptr), index(i)
    {
        R.resize(m);
        L.resize(m);
    }

    Job(int m, Job* p, Job* n, int i) : Previous(p), Next(n), index(i)
    {
        R.resize(m);
        L.resize(m);
    }
    
};

class SchedulingGraphList
{
public:
    Job* Head;
    Job* Tail;
    int machines;
    int JobsOrdered;

    SchedulingGraphList(int M): machines(M), JobsOrdered(0), Head(nullptr), Tail(nullptr)
    {

    }
    

    void InsertAfter(Job* prevNode, int index)
    {
        if (prevNode == nullptr)//insert at head
        {
            //allocate and set pointers and job index
            Job* newNode = new Job(machines, nullptr, Head, index);

            if (Head == nullptr)//empty
            {
                //fix Tail
                Tail = newNode;
            }
            else
            {
                //fix previous
                Head->Previous = newNode;
            }
            Head = newNode;


        }
        else
        {
            //allocate and set pointers and job index
            Job* newNode = new Job(machines, prevNode, prevNode->Next, index);

            //reset prevNode next pointer
            prevNode->Next = newNode;

            //fix backward pointer
            if (newNode->Next != nullptr)//isn't last
                newNode->Next->Previous = newNode;
            else//is last
                Tail = newNode;
            //increment counter
            JobsOrdered++;
        }
    }
    

    
};

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
int Cmax2(SchedulingGraphList& Graph, int* P, const vector<int> &X, int M, int N, int Elem, Job* Place)
{

    /*vector<int> R;
    vector<int> L;
    R.resize(M);
    L.resize(M);*/
    Job newJob(M, Elem);
    for (int i = 0; i < M; i++)
    {
        if (Place == nullptr && i == 0)
        {
            newJob.R[i] = P[Elem * M + i];
        }
        else if (Place == nullptr && i > 0)
        {
            newJob.R[i] = newJob.R[i-1] + P[Elem * M + i];
        }
        else if (Place != nullptr && i == 0)
        {
            newJob.R[i] = Place->R[i] + P[Elem * M + i]; //Tr[i][Place-1] + P[Elem * M + i];
        }
        else
        {
            int a1 = newJob.R[i - 1] + P[Elem * M + i];
            int a2 = Place->R[i] + P[Elem * M + i];
            newJob.R[i] = std::max(a1 , a2);
        }
    }

    
    if (Place == Graph.Tail)
    {
        return newJob.R[M - 1];
    }
    
    int highest = 0;
    if (Place == nullptr)
    {
        for (int i = 0; i < M; i++)
        {
            highest = max(newJob.R[i] + Graph.Head->L[i], highest);
        }
    }
    else
    {
        for (int i = 0; i < M; i++)
        {
            highest = max(newJob.R[i] + Place->Next->L[i], highest);
        }
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
        //szukamy najwy¿szej wagi z nieuszeregowanych
        
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

        //wynik czêœciowy

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

        
        //index zadania z najwiêksz¹ wag¹ z pozosta³ych
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

        //wynik czêœciowy
        X.insert(X.begin() + bestIndex, index);

    }
    
}

void algorytm2(SchedulingGraphList& Graph, int* P, vector<int> &X, int M, int N, vector<pair<int, int>> W, int* temp)//qneh
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

        //index zadania z najwiêksz¹ wag¹ z pozosta³ych
        index = l;// W[i].first;

        
        index = W[index].first;

        int bestIndex = 0;
        int c;
        int best = INT32_MAX;
        Job* iterator = nullptr;
        Job* bestIterator = nullptr;

        for (int k = 0; k < l + 1; k++)//2
        {
            //wstawiamy na nowe miejsce
            //X.insert(X.begin() + k, index);
            
            //sprawdzamy Cmax na k-tym miejscu
            

            //sprawdzamy czy jest lepiej
            c = Cmax2(Graph, P, X, M, l + 1, index, iterator);
            if (c < best)
            {
                best = c;
                bestIterator = iterator;
                bestIndex = k;
            }
            //X.erase(X.begin() + k);
            ///TODO iterator increment
            if (Graph.Head != nullptr)
            {
                if (k == 0)
                {
                    iterator = Graph.Head;
                }
                else
                {
                    iterator = iterator->Next;
                }
            }
            


        }

        //wynik czêœciowy
        

        //wstawianie elementu na w³aœciw miejsce
        Graph.InsertAfter(bestIterator, index);
        if (bestIterator != nullptr)//wybrano inne ni¿ pierwsze miejsce
            iterator = bestIterator->Next;
        else//pierwsze miejsce
            iterator = Graph.Head;
        ///TODO naprawianie grafu
        for (int i = bestIndex; i < l+1; i++)//po zadaniach
        {
            for (int j = 0; j < M; j++)//po maszynach
            {
                if (i == 0 && j == 0)
                {
                    iterator->R[j] = P[iterator->index * M + j];
                }
                else if (j > 0 && i == 0)
                {
                    iterator->R[j] = iterator->R[j - 1] + P[iterator->index * M + j];//[j-1][i]
                }
                else if (i > 0 && j == 0)
                {
                    iterator->R[j] = iterator->Previous->R[j] + P[iterator->index * M + j];//[j][i-1]
                }
                else
                {
                    iterator->R[j] = std::max(iterator->R[j - 1] + P[iterator->index * M + j], iterator->Previous->R[j] + P[iterator->index * M + j]);
                }
            }
            iterator = iterator->Next;
        }

        if (bestIterator != nullptr)//wybrano inne ni¿ ostatnie miejsce
            iterator = bestIterator->Next;
        else//ostatnie miejsce
            iterator = Graph.Head;
        for (int i = bestIndex; i >= 0; i--)
        {
            for (int j = M - 1; j >= 0; j--)
            {
                if (i == l && j == M - 1)
                {
                    iterator->L[j] = P[iterator->index * M + j];
                }
                else if (i == l && j < M - 1)
                {
                    iterator->L[j] = iterator->L[j+1] + P[iterator->index * M + j];
                }
                else if (i < l && j == M - 1)
                {
                    iterator->L[j] = iterator->Next->L[j] + P[iterator->index * M + j];
                }
                else
                {
                    iterator->L[j] = std::max(iterator->L[j+1] + P[iterator->index * M + j], iterator->Next->L[j] + P[iterator->index * M + j]);
                }
            }
            iterator = iterator->Previous;
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
    // T - czasy zakonczenia dla poszczególnych maszyn (nadpisywana) , 
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
    
    vector<vector<int>> TGraphL;
    TGraphL.resize(M);
    for (int i = 0; i < M; i++)
    {
        TGraphL[i].resize(N);
    }
    SchedulingGraphList Graph = SchedulingGraphList(M);

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
    algorytm2(Graph, P, X, M, N, W, temp);

    auto stop = chrono::steady_clock::now();
    chrono::duration<double> elapsed = stop - start;
    cout << "Czas: " << elapsed.count() << endl;
    cout << "Kolejnosc: ";

    ///TODO kolejnoœæ z grafu
    Job* it=Graph.Head;
    for (int i = 0; i < N; i++)
    {
        //do qneh ===========
        X[i] = it->index;
        it = it->Next;
        //===================
        cout << X[i] + 1 << " ";
    }
    //cout << endl << Cmax1(T,P,X,M,N) << endl;
    //cout << endl << Cmax(T,P,Xold,M,N) << endl;
    //cout << endl << TGraphL[0][0] << endl;
    cout << endl << Graph.Head->L[0] << endl;
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
