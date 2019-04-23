#include<cstdio>
#include<random>

using namespace std;

#define N 10000000
#define M 100000000
#define GRAIN_SIZE 32

struct Edge
{
    int u, v;
    float w;
} E[M];


// Takes the graph as input in an edge-list data structure

void input_graph(int &n, int &m)
{
    scanf("%d %d", &n, &m);
    for(int i = 1; i <= m; ++i)
        scanf("%d %d %f", &E[i].u, &E[i].v, &E[i].w);
}

void print_array(int *A, int n, const char *message)
{
    puts(message);

    for(int i = 1; i <= n; ++i)
        printf("%d\t", A[i]);

    printf("\n\n");
}


// Insertion sort for the array segment A[q : r]

void insertion_sort(int *A, int q, int r)
{
    for(int j = q + 1; j <= r; ++j)
    {
        int key = A[j], i = j - 1;
        while(i >= q && A[i] > key)
        {
            A[i + 1] = A[i];
            i--;
        }

        A[i + 1] = key;
    }
}


// Computes the prefix sum of the n-length array X into the array S

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


// Rearrange the elements of A[q : r] and return an index k in [q, r] such that all elements in A[q : k - 1] are smaller than x,
// and all elements in A[k + 1 : r] are larger than x, and A[k] = x

int parallel_partition(int *A, int q, int r, int x)
{
    int n = r - q + 1;

    if(n == 1)
        return q;

    int *B = new int[n + 1], *LT = new int[n + 1], *GT = new int[n + 1];

    for(int i = 1; i <= n; ++i)
    {
        B[i] = A[q + i - 1];
        LT[i] = (B[i] < x);
        GT[i] = (B[i] > x);
    }

    parallel_prefix_sum(LT, n, LT),
    parallel_prefix_sum(GT, n, GT);

    int k = q + LT[n];
    A[k] = x;

    for(int i = 1; i <= n; ++i)
        if(B[i] < x)
            A[q + LT[i] - 1] = B[i];
        else if(B[i] > x)
            A[k + GT[i]] = B[i];


    delete B, delete LT, delete GT;

    return k;
}

// Generate a thread-safe random number integer in the bound [low, high]
// Ref: https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers8
int random_int(const int &low, const int &high)
{
    static thread_local mt19937 generator;
    uniform_int_distribution<int> distribution(low, high);

    return distribution(generator);
}


// Randomized quicksort for the array segment A[q : r]

void parallel_randomized_quicksort(int *A, int q, int r)
{
    int n = r - q + 1;

    if(n <= GRAIN_SIZE)
        insertion_sort(A, q, r);
    else
    {
        int randomIdx = random_int(q, r);
        int x = A[randomIdx];

        //print_array(A + q - 1, n, "Un-partitioned array");

        int k = parallel_partition(A, q, r, x);

        //printf("Position of element %d is %d\n", x, k);
        //print_array(A + q - 1, n, "Partitioned array");

        //print_array(A, 16, "Full array:");

        parallel_randomized_quicksort(A, q, k - 1);
        parallel_randomized_quicksort(A, k + 1, r);
    }
}


bool is_non_decreasing(int *A, int n)
{
    for(int i = 2; i <= n; ++i)
        if(A[i] < A[i - 1])
            return false;

    return true;
}

void generate_random_array(int *A, int n, int low, int high)
{
    for(int i = 1; i <= n; ++i)
    {
        A[i] = random_int(low, high);
        for(int j = 1; j < i; ++j)
            if(A[j] == A[i])
            {
                i--;
                break;
            }
    }
}

void test()
{
    int X[100001], n = 9800;
    int low = 0, high = 1000000;

    generate_random_array(X, n, low, high);
    print_array(X, n, "Initial array:");

    //parallel_partition(X, 1, n, X[1]);
    //print_array(X, n, "Same array but now partitioned around the 1st element");


    parallel_randomized_quicksort(X, 1, n);
    print_array(X, n, "Sorted array:");

    if(!is_non_decreasing(X, n))
        puts("Incorrect sorting\n\n\n");
    else
        puts("Correct sorting\n\n\n");


    parallel_prefix_sum(X, n, X);

    print_array(X, n, "Prefix sums:");

}

int main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    int n, m;

    input_graph(n, m);

    test();

    return 0;
}
