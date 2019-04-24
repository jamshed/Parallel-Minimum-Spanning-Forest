#include<cstdio>
#include<cmath>
#include<random>

using namespace std;

#define N 10001
#define M 100001
#define GRAIN_SIZE 1
#define HEAD 0
#define TAIL 1

struct Edge
{
    int u, v;
    float w;

    void set_edge(int e_u, int e_v, float e_w)
    {
        u = e_u, v = e_v, w = e_w;
    }

    void print()
    {
        printf("(%d, %d): %f\n", u, v, w);
    }
} Edges[2 * M];


// Generate a thread-safe random number integer in the bound [low, high]
// Ref: https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers8
int random_int(const int &low, const int &high)
{
    static thread_local mt19937 generator;
    uniform_int_distribution<int> distribution(low, high);

    return distribution(generator);
}


// Generate a thread-safe random number float in the bound [low, high]

float random_float(const float &low, const float &high)
{
    static thread_local mt19937 generator;
    uniform_real_distribution<float> distribution(low, high);

    return distribution(generator);
}


// Takes the graph as input in an edge-list data structure

void input_graph(Edge *E, int &n, int &m)
{
    scanf("%d %d", &n, &m);

    for(int i = 1; i <= m; ++i)
        scanf("%d %d %f", &E[i].u, &E[i].v, &E[i].w);
}


// Generate a random m-length edge-list, with vertices in [1, n], and real-valued edge weights in [low, high]

void generate_random_edge_list(Edge *E, int n, int &m, float low, float high)
{
    for(int i = 1; i <= m; ++i)
        E[i].u = random_int(1, n), E[i].v = random_int(1, n), E[i].w = random_float(low, high);
}


// Print the m-length edge-list E, preceded with 'message'

void print_edge_list(Edge *E, int m, const char *message)
{
    puts(message);

    for(int i = 1; i <= m; ++i)
        E[i].print();

    printf("\n~~~END~~~\n\n");
}


// Checks if the m-length edge-list E is correctly sorted in non-decreasing order

bool is_non_decreasing(Edge *E, int m)
{
    for(int i = 2; i <= m; ++i)
        if(E[i].w < E[i - 1].w)
        {
            printf("Wrong order --> (%f, %f)\n\n", E[i - 1].w, E[i].w);
            return false;
        }

    return true;
}



// Insertion sort for the edge segment E[q : r]

void insertion_sort(Edge *E, int q, int r)
{
    for(int j = q + 1; j <= r; ++j)
    {
        Edge e = E[j];
        float key = e.w;
        int i = j - 1;

        while(i >= q && E[i].w > key)
        {
            E[i + 1] = E[i];
            i--;
        }

        E[i + 1] = e;
    }
}


// Computes the prefix sum of the n-length integer array X into the array S

void parallel_prefix_sum(int *X, int n, int *S)
{
    if(n == 1)
        S[1] = X[1];
    else
    {
        int *Y = new int[n / 2 + 1], *Z = new int[n / 2 + 1];

        for(int i = 1; i <= n / 2; ++i)
            Y[i] = X[2*i - 1] + X[2 * i];

        parallel_prefix_sum(Y, n / 2, Z);

        for(int i = 1; i <= n; ++i)
            if(i == 1)
                S[i] = X[1];
            else if(i % 2 == 0)
                S[i] = Z[i / 2];
            else
                S[i] = Z[(i - 1) / 2] + X[i];

        delete Y, delete Z;
    }
}


// Rearrange the edges of E[q : r] and return an index k in [q, r] such that all edges in E[q : k - 1] have smaller weights than E[k],
// and all edges in E[k + 1 : r] have higher weights than E[k], and E[k] is the pivot edge, which was situated at 'pivotIdx' earlier

int parallel_partition(Edge *E, int q, int r, int pivotIdx)
{
    int n = r - q + 1;

    if(n == 1)
        return q;

    Edge pivot = E[pivotIdx];
    Edge *B = new Edge[n + 1];
    int *LT = new int[n + 1], *GT = new int[n + 1];

    for(int i = 1; i <= n; ++i)
    {
        B[i] = E[q + i - 1];

        if(B[i].w < pivot.w)
            LT[i] = 1, GT[i] = 0;
        else if(B[i].w > pivot.w)
            LT[i] = 0, GT[i] = 1;
        else if(q + i - 1 < pivotIdx)
            LT[i] = 1, GT[i] = 0;
        else if(q + i - 1 > pivotIdx)
            LT[i] = 0, GT[i] = 1;
        else
            LT[i] = GT[i] = 0;
    }

    parallel_prefix_sum(LT, n, LT),
    parallel_prefix_sum(GT, n, GT);

    int k = q + LT[n];
    E[k] = pivot;

    for(int i = 1; i <= n; ++i)
        if(B[i].w < pivot.w)
            E[q + LT[i] - 1] = B[i];
        else if(B[i].w > pivot.w)
            E[k + GT[i]] = B[i];
        else if(q + i - 1 < pivotIdx)
            E[q + LT[i] - 1] = B[i];
        else if(q + i - 1 > pivotIdx)
            E[k + GT[i]] = B[i];


    delete B, delete LT, delete GT;

    return k;
}


// Randomized quicksort for the edge-list E[q : r]

void parallel_randomized_quicksort(Edge *E, int q, int r)
{
    int n = r - q + 1;

    if(n <= GRAIN_SIZE)
        insertion_sort(E, q, r);
    else
    {
        int randomIdx = random_int(q, r);
        int k = parallel_partition(E, q, r, randomIdx);

        parallel_randomized_quicksort(E, q, k - 1);
        parallel_randomized_quicksort(E, k + 1, r);
    }
}


// Replace each undirected edge (u, v, w) in the m-length edge-list E with
// two directed edges: (u, v, w) and (v, u, w)

void parallel_get_directed_edge_list(Edge *E, int &m)
{
    Edge *B = new Edge[m + 1];

    for(int i = 1; i <= m; ++i)
        B[i] = E[i];

    for(int i = 1; i <= m; ++i)
        E[2 * i - 1].set_edge(B[i].u, B[i].v, B[i].w),
        E[2 * i].set_edge(B[i].v, B[i].u, B[i].w);

    m *= 2;

    delete B;
}


// Radix sort in non-decreasing order for an n-length long integer array A of b-bit long integers

void parallel_radix_sort(long long *A, int n, int b)
{
    int *F0 = new int[n + 1], *F1 = new int[n + 1], *S0 = new int[n + 1], *S1 = new int[n + 1];
    long long *B = new long long[n + 1];


    for(int k = 0; k < b; ++k)
    {
        for(int i = 1; i <= n; ++i)
        {
            F1[i] = (A[i] >> k) & 1;
            F0[i] = 1 - F1[i];
        }

        parallel_prefix_sum(F0, n, S0);
        parallel_prefix_sum(F1, n, S1);

        for(int i = 1; i <= n; ++i)
            if(!F1[i])
                B[S0[i]] = A[i];
            else
                B[S0[n] + S1[i]] = A[i];

        for(int i = 1; i <= n; ++i)
            A[i] = B[i];
    }


    //delete F0, delete F1, delete S0, delete S1, delete B;
}


// For n vertices (in [1, n]) and m-length edge-list E, for each vertex u in [1, n],
// R[u] is set to the smallest index i such that E[i].u = u

void parallel_simulate_priority_CW_radix_sort(Edge *E, int n, int m, int *R)
{
    long long *A = new long long[m + 1];
    int k = (int)ceil(log2(m)) + 1, u, j;

    for(int i = 1; i <= m; ++i)
        A[i] = (E[i].u != E[i].v ? (((long long)E[i].u << k) | i) : 0);

    parallel_radix_sort(A, m, k + (int)ceil(log2(n)));

    for(int i = 1; i <= m; ++i)
        if(A[i])
        {
            u = (A[i] >> k);
            j = A[i] - (u << k);

            if(i == 1 || u != (A[i - 1] >> k))
                R[u] = j;
        }

    delete A;
}


inline bool coin_toss()
{
    return random_int(0, 1);
}


// For n vertices in [1, n] and m-length edge-list E, the minimum spanning forest for the graph is computed
// at the Boolean array MSF; MSF[i] = true if and only if the i'th edge is in the computed minimum spanning forest.
// Preconditions: MSF[1 : m] are set to false; for every undirected edge (u, v), both (u, v) and (v, u) are in E.

void parallel_randomized_MSF_priority_CW(Edge *E, int n, int m, bool *MSF)
{
    int *L = new int[n + 1], *R = new int[n + 1], u, v;
    bool *C = new bool[n + 1];
    Edge *B = new Edge[m + 1];

    parallel_randomized_quicksort(E, 1, m);

    for(int i = 1; i <= m; ++i)
        B[i] = E[i];

    for(int i = 1; i <= n; ++i)
        L[i] = i;

    bool edgesRemain = (m > 0);

    //printf("MSF routine with n = %d, m = %d\n", n, m);
    //print_edge_list(Edges, m, "Sorted");

    while(edgesRemain)
    {
        for(int i = 1; i <= n; ++i)
            C[i] = coin_toss();

        //print_edge_list(E, m, "Intermediate:");

        parallel_simulate_priority_CW_radix_sort(E, n, m, R);

        for(int i = 1; i <= m; ++i)
        {
            u = E[i].u, v = E[i].v;
            if(C[u] == TAIL && C[v] == HEAD && R[u] == i)
                L[u] = v, MSF[i] = true;//, printf("Weight %f included in MSF\n", E[i].w);
        }

        for(int i = 1; i <= m; ++i)
            E[i].u = L[E[i].u], E[i].v = L[E[i].v];

        edgesRemain = false;
        for(int i = 1; i <= m; ++i)
            if(E[i].u != E[i].v)
                edgesRemain = true;//, printf("(%d, %d) remain\n", E[i].u, E[i].v);
    }

    for(int i = 1; i <= m; ++i)
        E[i] = B[i];

    delete L, delete R, delete C, delete B;
}


void MSF(Edge *E, int n, int m)
{
    parallel_get_directed_edge_list(E, m);

    bool *MSF = new bool[m + 1];
    for(int i = 1; i <= m; ++i)
        MSF[i] = false;

    parallel_randomized_MSF_priority_CW(Edges, n, m, MSF);

    int edgeCount = 0;
    float sumCost = 0;

    for(int i = 1; i <= m; ++i)
        if(MSF[i])
            edgeCount++, sumCost += E[i].w;

    printf("Edge count in MSF = %d\nMSF Cost = %f\n", edgeCount, sumCost);

    printf("\nMSF edges:\n");
    for(int i = 1; i <= m; ++i)
        if(MSF[i])
            E[i].print();

    delete MSF;
}


void test()
{
    /*
    int n = 1000, m = 100000;
    float low = 0, high = 1000000;

    generate_random_edge_list(Edges, n, m, low, high);
    //print_edge_list(Edges, m, "Initial edge list:");

    parallel_randomized_quicksort(Edges, 1, m);
    //print_edge_list(Edges, m, "Sorted edge list:");

    if(!is_non_decreasing(Edges, m))
        puts("Incorrect sorting\n");
    else
        puts("Correct sorting\n");
    */

    int n, m;

    input_graph(Edges, n, m);
    print_edge_list(Edges, m, "Initial");

    /*parallel_randomized_quicksort(Edges, 1, m);
    print_edge_list(Edges, m, "Sorted");

    int *R = new int[n + 1];
    parallel_simulate_priority_CW_radix_sort(Edges, n, m, R);

    printf("\nRank:\n");
    for(int i = 1; i <= n; ++i)
        printf("R[%d] = %d\n", i, R[i]);

    delete R;
        */

    MSF(Edges, n, m);
}

int main()
{
    freopen("input.txt", "r", stdin);
    //freopen("output.txt", "w", stdout);

    //int n, m;

    //input_graph(n, m, Edges);

    //for(int i = 0; i < 10; ++i)
        test();

    return 0;
}
