#include<cstdio>
#include<cmath>
#include<random>

using namespace std;

#define N 4000006
#define M 34000006
#define GRAIN_SIZE 1
#define HEAD 0
#define TAIL 1

const int P = 192; // #processing elements / cores

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


// Generates a thread-safe random integer in the bound [low, high]
// Ref: https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers8
int random_int(const int &low, const int &high)
{
    static thread_local mt19937 generator;
    uniform_int_distribution<int> distribution(low, high);

    return distribution(generator);
}


// Generates a thread-safe random float in the bound [low, high]
// Ref: https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers8

float random_float(const float &low, const float &high)
{
    static thread_local mt19937 generator;
    uniform_real_distribution<float> distribution(low, high);

    return distribution(generator);
}


// Simulate a coin toss and return 0 or 1 (HEAD or TAIL)

inline bool coin_toss()
{
    return random_int(0, 1);
}


// Generate a random m-length edge-list, with vertices in [1, n], and real-valued edge weights in [low, high]

void generate_random_edge_list(Edge *E, int n, int &m, float low, float high)
{
    for(int i = 1; i <= m; ++i)
        E[i].u = random_int(1, n), E[i].v = random_int(1, n), E[i].w = random_float(low, high);
}




// Take a weighted graph as input in an edge-list E;
// n is the number of vertices and m is the number of edges.

void input_graph(Edge *E, int &n, int &m)
{
    scanf("%d %d", &n, &m);

    for(int i = 1; i <= m; ++i)
        scanf("%d %d %f", &E[i].u, &E[i].v, &E[i].w);
}


// Print the m-length edge-list E, preceded with a 'message'

void print_edge_list(Edge *E, int m, const char *message)
{
    puts(message);

    for(int i = 1; i <= m; ++i)
        E[i].print();

    printf("\n---END---\n\n");
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
            Y[i] = X[2 * i - 1] + X[2 * i];

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


// Rearrange the edges of E[q : r] and return an index k in [q, r] such that
// all edges in E[q : k - 1] have smaller weights than E[k], and all edges in E[k + 1 : r] have higher weights than E[k],
// and E[k] is the pivot edge, which was situated at 'pivotIdx' earlier

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


// Randomized Quicksort for the edge-list E[q : r]

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
// two directed edges: (u, v, w) and (v, u, w).

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


    delete F0, delete F1, delete S0, delete S1, delete B;
}


// For n vertices (in [1, n]) and m-length edge-list E, for each vertex u in [1, n],
// R[u] is set to the smallest index i such that E[i].u = u and E[i] is not a loop.

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


// For n vertices in [1, n] and m-length edge-list E, the minimum spanning forest for the graph is computed
// at the Boolean array MSF; MSF[i] = true if and only if the i'th edge is in the computed minimum spanning forest.
// Preconditions: MSF[1 : m] are set to false; for every undirected edge (u, v), both (u, v) and (v, u) are in E.

void parallel_randomized_MSF_priority_CW_radix_sort(Edge *E, int n, int m, bool *MSF)
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

    while(edgesRemain)
    {
        for(int i = 1; i <= n; ++i)
            C[i] = coin_toss();

        parallel_simulate_priority_CW_radix_sort(E, n, m, R);

        for(int i = 1; i <= m; ++i)
        {
            u = E[i].u, v = E[i].v;
            if(C[u] == TAIL && C[v] == HEAD && R[u] == i)
                L[u] = v, MSF[i] = true;//, printf("Weight %f included in MSF\n", E[i].w);
        }

        edgesRemain = false;
        for(int i = 1; i <= m; ++i)
        {
            E[i].u = L[E[i].u], E[i].v = L[E[i].v];
            if(E[i].u != E[i].v)
                edgesRemain = true;//, printf("(%d, %d) remain\n", E[i].u, E[i].v);
        }
    }

    for(int i = 1; i <= m; ++i)
        E[i] = B[i];

    delete L, delete R, delete C, delete B;
}



// Ranking d-bit integer keys in n-length integer array S, using Counting sort (stable);
// R[i] provides the rank of S[i] when keys in S are sorted in non-decreasing order

void parallel_counting_rank(int *S, int n, int d, int *R)
{
    int **f = new int*[1 << d], **r1 = new int*[1 << d], *j_s = new int[P + 1], *j_e = new int[P + 1], *ofs = new int[P + 1];

    for(int i = 0; i < (1 << d); ++i)
        f[i] = new int[P + 1], r1[i] = new int[P + 1];

    for(int i = 1; i <= P; ++i)
    {
        for(int j = 0; j < (1 << d); ++j)
            f[j][i] = 0;

        j_s[i] = (i - 1) * (n / P) + 1, j_e[i] = (i < P ? i * (n / P) : n);
        for(int j = j_s[i]; j <= j_e[i]; ++j)
            f[S[j]][i]++;
    }

    for(int j = 0; j < (1 << d); ++j)
        parallel_prefix_sum(f[j], P, f[j]);

    for(int i = 1; i <= P; ++i)
    {
        ofs[i] = 1;
        for(int j = 0; j < (1 << d); ++j)
        {
            r1[j][i] = (i == 1 ? ofs[i] : ofs[i] + f[j][i - 1]);
            ofs[i] = ofs[i] + f[j][P];
        }

        for(int j = j_s[i]; j <= j_e[i]; ++j)
        {
            R[j] = r1[S[j]][i];
            r1[S[j]][i]++;
        }
    }


    for(int i = 0; i < (1 << d); ++i)
        delete f[i], delete r1[i];

    delete f, delete r1, delete j_s, delete j_e, delete ofs;
}


// Extract the bit segment from the s'th to the e'th position of the unsigned long integer n;
// [assuming that the extracted bit segment fits into an integer]

inline int extract_bit_segment(unsigned long long n, int s, int e)
{
    const int sz = 8 * sizeof(unsigned long long);

    return (n << (sz - 1 - e)) >> (sz - 1 - e + s);
}


// For an n-length b-bit long integer array A, stable sort it in non-decreasing order
// using Radix sort with ranking using Counting sort

void parallel_radix_sort_counting_rank(long long *A, int n, int b)
{
    int *S = new int[n + 1], *R = new int[n + 1];
    long long *B = new long long[n + 1];

    // printf("n = %d, P = %d\n", n, P);

    int d = (int)ceil(log2(n / (P * log2(n)))), q;

    for(int k = 0; k < b; k += d)
    {
        q = (k + d <= b ? d : b - k);

        for(int i = 1; i <= n; ++i)
            S[i] = extract_bit_segment((unsigned long long)A[i], k, k + q - 1);
            /*if(S[i] < 0)
            {
                printf("exiting as S[%d] = %d is invalid\n", i, S[i]);
                exit(1);
            }*/

        parallel_counting_rank(S, n, q, R);

        for(int i = 1; i <= n; ++i)
            B[R[i]] = A[i];

        for(int i = 1; i <= n; ++i)
            A[i] = B[i];
    }

    delete S, delete R, delete B;
}


// For n vertices (in [1, n]) and m-length edge-list E, for each vertex u in [1, n],
// R[u] is set to the smallest index i such that E[i].u = u and E[i] is not a loop.

void parallel_simulate_priority_CW_radix_sort_counting_rank(Edge *E, int n, int m, int *R)
{
    long long *A = new long long[m + 1];
    int k = (int)ceil(log2(m)) + 1, u, j;

    for(int i = 1; i <= m; ++i)
        A[i] = (E[i].u != E[i].v ? (((long long)E[i].u << k) | i) : 0);
        /*if(A[i] < 0)
        {
            printf("exiting as A[i] = %lld\n", A[i]);
            exit(1);
        }*/

    parallel_radix_sort_counting_rank(A, m, k + (int)ceil(log2(n)));

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


// For n vertices in [1, n] and m-length edge-list E, the minimum spanning forest for the graph is computed
// at the Boolean array MSF; MSF[i] = true if and only if the i'th edge is in the computed minimum spanning forest.
// Preconditions: MSF[1 : m] are set to false; for every undirected edge (u, v), both (u, v) and (v, u) are in E.

void parallel_randomized_MSF_priority_CW_radix_sort_counting_rank(Edge *E, int n, int m, bool *MSF)
{
    int *L = new int[n + 1], *R = new int[n + 1], u, v;
    bool *C = new bool[n + 1];
    Edge *B = new Edge[m + 1];

    parallel_randomized_quicksort(E, 1, m);

    for(int i = 1; i <= m; ++i)
        B[i] = E[i];
        /*if(E[i].u < 1 || E[i].u > n || E[i].v < 1 || E[i].v > n)
        {
            printf("Invalid edge after sorting (%d, %d)\n", E[i].u, E[i].v);
            printf("Exiting due to invalid edge.\n");
            exit(1);
        }*/

    for(int i = 1; i <= n; ++i)
        L[i] = i;

    bool edgesRemain = (m > 0);

    while(edgesRemain)
    {
        for(int i = 1; i <= n; ++i)
            C[i] = coin_toss();

        parallel_simulate_priority_CW_radix_sort_counting_rank(E, n, m, R);

        /*for(int i = 1; i <= m; ++i)
            if(E[i].u != E[i].v)
            {
                if(E[i].u < 1 || E[i].u > n || E[i].v < 1 || E[i].v > n)
                {
                    printf("Invalid edge after priority CW (%d, %d)\n", u, v);
                    printf("Exiting due to invalid edge.\n");
                    exit(1);
                }
                else if(E[i].u != E[i].v)
                {
                    int r = R[E[i].u];
                    if(r < 1 || r > m)
                    {
                        printf("Exiting due to invalid rank %d for vertex %d. Valid edge (%d, %d) (edge %d) exists.\n",
                               r, E[i].u, E[i].u, E[i].v, i);
                        exit(1);
                    }
                }
            }*/

        for(int i = 1; i <= m; ++i)
        {
            u = E[i].u, v = E[i].v;

            if(C[u] == TAIL && C[v] == HEAD && R[u] == i)
                L[u] = v, MSF[i] = true;
        }

        /*for(int i = 1; i <= n; ++i)
            if(L[i] < 1 || L[i] > n)
            {
                printf("Label corrupted for vertex %d, (label = %d)\n", i, L[i]);
                break;
            }*/

        edgesRemain = false;
        for(int i = 1; i <= m; ++i)
        {
            E[i].u = L[E[i].u], E[i].v = L[E[i].v];
            if(E[i].u != E[i].v)
                edgesRemain = true;//, printf("(%d, %d) remain\n", E[i].u, E[i].v);
        }
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

    //parallel_randomized_MSF_priority_CW_radix_sort(E, n, m, MSF);

    parallel_randomized_MSF_priority_CW_radix_sort_counting_rank(E, n, m, MSF);

    int edgeCount = 0;
    float sumCost = 0;

    for(int i = 1; i <= m; ++i)
        if(MSF[i])
            edgeCount++, sumCost += E[i].w;

    printf("Edge count in MSF = %d\nMSF Cost = %f\n", edgeCount, sumCost);

    /*
    printf("\nMSF edges:\n");
    for(int i = 1; i <= m; ++i)
        if(MSF[i])
            E[i].print();
    */

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
    Edge *E = Edges;

    //input_graph(Edges, n, m);
    //print_edge_list(Edges, m, "Initial");

    n = 10000, m = 1000000;
    generate_random_edge_list(E, n, m, 0, 1000);

    for(int i = 1; i <= m; ++i)
        if(E[i].u < 1 || E[i].u > n || E[i].v < 1 || E[i].v > n)
        {
            puts("Error in random edge-list generation.\n");
            exit(1);
        }

    /*parallel_randomized_quicksort(Edges, 1, m);
    print_edge_list(Edges, m, "Sorted");

    int *R = new int[n + 1];
    parallel_simulate_priority_CW_radix_sort(Edges, n, m, R);

    printf("\nRank:\n");
    for(int i = 1; i <= n; ++i)
        printf("R[%d] = %d\n", i, R[i]);

    delete R;
        */

    MSF(E, n, m);
}

int main()
{
    freopen("input.txt", "r", stdin);
    //freopen("output.txt", "w", stdout);

    //int n, m;

    //input_graph(n, m, Edges);

    //for(int i = 0; i < 10; ++i)
        test();

        //printf("extracted %d\n", extract_bit_segment(7, 1, 3));

    return 0;
}
