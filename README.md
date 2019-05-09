# Parallel-Minimum-Spanning-Forest

## Compile
----------

### For small graphs:

To work with static array sizes < 2GB:

`icc Parallel-MSF.cpp -o parallel-msf -std=c++11`

This will work for the current version of *Parallel-MSF.cpp*. Specifically, use this for graphs with edge-counts up-to ≈ 67M. For larger graphs, use the following compilation.

### For large graphs:

To work with static array sizes > 2GB:

`icc Parallel-MSF.cpp -o parallel-msf -std=c++11 -mcmodel medium -shared-intel`

Note that at the *Parallel-MSF.cpp*, you need to modify the macro `M` to suit your edge-list length. Specifically, modify your edge-count at line 13, `#define M 34690007`.

## Run
------

`./parallel-msf -p <process-count> -cw <concurrent-write-method> -f <input-file>`

- The `cw` parameter value should be in {0: Radix sort, 1: Radix sort with ranking by counting sort, 2: Binary search}.

- The input file should be inside a directory named **turn-in**, which should be inside the same directory as the executable *parallel-msf*.

- There should also exist a directory named **outputs** at the same directory as the executable *parallel-msf*, inside which the output files will be generated.


For example, to run the executable *parallel-msf* for a file named *com-lj-in.txt* using 48 processes and using the binary search based parallel concurrent-write, execute the following.

`./parallel-msf -p 48 -cw 2 -f turn-in/com-lj-in.txt`

The output file will be *outputs/com-lj-MST-search-out.txt*.
