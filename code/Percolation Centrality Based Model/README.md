## Code : percolation_centrality.cu
### Author : Sayantan Jana 
compile instruction : 
```sh
nvcc percolation_centrality.cu
```
run instruction : 
```sh
./a.out < <input_file> > <output_file>
```
Note that the input would be of the following form :
- First line contains 2 space separated integers N and M, denoting the count of nodes and edges
- M lines follow describing the edges containing 2 space separated integers u and v, denoting there is an edge present between u and v.
- Sample structure of input :
```sh
N M
u1 v1
u2 v2
.
.
.
uM vM
```
Note that the percolation has been assumed as a function 1/i for the simplicity, it can be changed as desired in line 153.