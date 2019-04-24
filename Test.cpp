#include<cstdio>
#include<random>

using namespace std;

#define N 10000001
#define M 100000001
#define GRAIN_SIZE 32

struct Edge
{
    int u, v;
    float w;

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

void input_graph(int &n, int &m, Edge *E)
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

    printf("\nEND\n\n");
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

    //int lower = 0, higher = 0;
    for(int i = 1; i <= n; ++i)
    {
        B[i] = E[q + i - 1];
        //lower += (LT[i] = (B[i].w < pivot.w));
        //higher += (GT[i] = (B[i].w > pivot.w));

        if(B[i].w < pivot.w)
            LT[i] = 1, GT[i] = 0;//, lower++;
        else if(B[i].w > pivot.w)
            LT[i] = 0, GT[i] = 1;//, higher++;
        else if(q + i - 1 < pivotIdx)
            LT[i] = 1, GT[i] = 0;//, lower++;
        else if(q + i - 1 > pivotIdx)
            LT[i] = 0, GT[i] = 1;//, higher++;
        else
            LT[i] = GT[i] = 0;
    }

    parallel_prefix_sum(LT, n, LT),
    parallel_prefix_sum(GT, n, GT);

    /*if(LT[n] + GT[n] != n - 1)
    {
        puts("ERR in prefix sums");
        printf("prefix sums for %d elements: lt_n = %d, gt_n = %d; manually l = %d, r = %d\n", n, LT[n], GT[n], lower, higher);
        printf("lt_n + rt_n = %d, l + h = %d\n\n", LT[n] + GT[n], lower + higher);
    }*/

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


void test()
{
    int n = 100000, m = 10000000;
    float low = 0, high = 1000000;

    generate_random_edge_list(Edges, n, m, low, high);
    //print_edge_list(Edges, m, "Initial edge list:");

    parallel_randomized_quicksort(Edges, 1, m);
    //print_edge_list(Edges, m, "Sorted edge list:");

    if(!is_non_decreasing(Edges, m))
        puts("Incorrect sorting\n");
    else
        puts("Correct sorting\n");
}

int main()
{
    //freopen("input.txt", "r", stdin);
    //freopen("output.txt", "w", stdout);

    //int n, m;

    //input_graph(n, m, Edges);

    test();

    return 0;
}
