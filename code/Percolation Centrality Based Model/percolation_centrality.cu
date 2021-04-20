#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <vector>
#include <set>
#include <time.h>
#include <iostream>
#include <algorithm>
#include <assert.h>

using namespace std;
/*
compile instruction : nvcc percolation_centrality.cu
run instruction : ./a.out < <input_file> > <output_file>
Note that the input would be of the following form :
First line contains 2 space separated integers N and M, denoting the count of nodes and edges
M lines follow describing the edges containing 2 space separated integers u and v, denoting
there is an edge present between u and v.
Sample structure of input :
N M
u1 v1
u2 v2
.
.
.
uM vM
Note that the percolation has been assumed as a function 1/i for the simplicity, it can be changed
as desired in line 153.
*/

#define NUM_THREADS 32
#define NUM_BLOCKS 1024

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__global__ void bfs(int V, int E, int *dColumn, int *dRow, int *Distance, int *Queue,
float *Paths, int *Dist, int *Sigma)
{
	__shared__ int arr[1];
	int *QLen = arr;
	
	int rootIndex = blockIdx.x + 1;
	int *Q = Queue + (blockIdx.x)*(V+1);
	float *dPaths = Paths + (blockIdx.x)*(V+1);
	int *dDistance = Distance + (blockIdx.x)*(V+1);

	while(rootIndex <= V)
	{
		for(int i=threadIdx.x; i<=V; i+=NUM_THREADS) dPaths[i] = 0;
		for(int i=threadIdx.x; i<=V; i+=NUM_THREADS) dDistance[i] = -1;

		if(threadIdx.x==0)
		{
			*QLen = 1;
			int root = rootIndex;
			Q[0] = root;
			dPaths[root] = 1.0f;
			dDistance[root] = 0;
		}
		__syncthreads();

		int oldQLen = 0;
		while(oldQLen < *QLen)
		{
			int id = threadIdx.x;
			int	source = Q[oldQLen++];
			int degree = dRow[source+1] - dRow[source];

			while(id < degree)
			{
				int neighbour = dColumn[dRow[source]+id];
				if(dDistance[neighbour] == -1)
				{
					dDistance[neighbour] = dDistance[source]+1;
					Q[atomicAdd(QLen, 1)] = neighbour;
				}
				if(dDistance[neighbour] == dDistance[source]+1)
					dPaths[neighbour] += dPaths[source];
					
				id += NUM_THREADS;
			}
			__syncthreads();
		}
		__syncthreads();
		
		for(int i=threadIdx.x; i<=V; i+=NUM_THREADS)
		{
			if(i == 0) continue;
			Dist[(rootIndex-1)*V+i] = dDistance[i];
			Sigma[(rootIndex-1)*V+i] = dPaths[i];
		}
		rootIndex += NUM_BLOCKS;
	}
}

__global__ void percolation_estimate(int V, int E, float *dCentrality,int *sDist,int *sSigma,float *pc)
{
	int k = blockIdx.x + 1;
	while(k<=V)
	{
		float v[NUM_THREADS]; 
		for(int i=threadIdx.x; i<=V; i+=NUM_THREADS)
		{
			if(i == 0) continue;
			if(i == k) continue;
			for(int j=1; j<=V; j++)
			{
				if(i == j) continue;
				if(k == j) continue;
				v[threadIdx.x] = 0;
				if(sDist[(k-1)*V+i]+sDist[(i-1)*V+j] == sDist[(k-1)*V+j])
				{
					v[threadIdx.x] = sSigma[(k-1)*V+i]*sSigma[(i-1)*V+j];
					v[threadIdx.x] = v[threadIdx.x]/(float)(sSigma[(k-1)*V+j]);
					v[threadIdx.x] = v[threadIdx.x]*max(0.0f,pc[k]-pc[j]);
					atomicAdd(&dCentrality[i], v[threadIdx.x]);
				}
				// printf("%d %d %d %f\n",k,i,j,v[threadIdx.x]);
			}
		}
		k += NUM_BLOCKS;
	}
}

int main()
{
	int V, E;
	cin >> V >> E;

	vector <vector <int> > graph(V+1);
	for(int i=0; i<E; ++i)
	{
		int u, v;
		cin >> u >> v;
		if(u == v) continue;
		graph[u].push_back(v);
		graph[v].push_back(u);
	}

	int *hColumn = new int[2*E];
	int *hRow	 = new int[V+2];
	float *perc  = new float[V+2];

	for(int i=1;i<=V;++i)
		perc[i] = 1.0/(float)(i);
	perc[0] = perc[V+1] = 1.0;

	for(int index=0, i=1; i<=V; ++i) 
	{
		for(int j=0;j<(int)graph[i].size();++j)
		{
			int n = graph[i][j]; 
			hColumn[index++] = n;
		}
	}
	
	// Filling row array
	long count = 0;
	for(int i=0; i<=V;)
	{
		for(int j=0;j<(int)graph.size();++j)
		{
			vector<int> v = graph[i];
			hRow[i++] = count;
			count += v.size();
		}
	}
	hRow[V+1] = count;

	float *Paths;
	int *Dist, *Sigma;
	int *dColumn, *dRow, *Distance, *Queue;

	cudaMalloc((void**)&dRow,    		sizeof(int)*(V+2));
	cudaMalloc((void**)&Dist,			sizeof(int)*(V*V+2));
	cudaMalloc((void**)&Sigma,			sizeof(int)*(V*V+2));
	cudaMalloc((void**)&dColumn, 		sizeof(int)*(2*E));
	cudaMalloc((void**)&Queue,    		sizeof(int)*(V+1)*NUM_BLOCKS);
	cudaMalloc((void**)&Distance,		sizeof(int)*(V+1)*NUM_BLOCKS);
	cudaMalloc((void**)&Paths,			sizeof(float)*(V+1)*NUM_BLOCKS);

	cudaMemcpy(dRow, hRow, sizeof(int)*(V+2),cudaMemcpyHostToDevice);
	cudaMemcpy(dColumn, hColumn, sizeof(int)*(2*E), cudaMemcpyHostToDevice);
	gpuErrchk( cudaPeekAtLastError() );

	bfs <<<NUM_BLOCKS, NUM_THREADS, 32>>> (V, E, dColumn, dRow, Distance, Queue, Paths, Dist, Sigma);
	cudaDeviceSynchronize();
	gpuErrchk( cudaPeekAtLastError() );
	
	cudaDeviceSynchronize();

	int *GetDist = new int[V*V+2];
	int *GetSigma = new int[V*V+2];
	cudaMemcpy(GetDist, Dist, sizeof(int)*(V*V+2), cudaMemcpyDeviceToHost);
	cudaMemcpy(GetSigma, Sigma, sizeof(int)*(V*V+2), cudaMemcpyDeviceToHost);
	/*
	for(int i=1; i<=V; ++i)
	{
		for(int j=1;j<=V; ++j)
		{
			printf("%d ",GetSigma[(i-1)*V+j]);
		}
		printf("\n");
	}
	*/
	float *dCentrality, *pc;
	int *sSigma,*sDist;

	cudaMalloc((void**)&dCentrality,	sizeof(float)*(V+2));
	cudaMalloc((void**)&sDist,			sizeof(int)*(V*V+2));
	cudaMalloc((void**)&sSigma,			sizeof(int)*(V*V+2));
	cudaMalloc((void**)&pc,				sizeof(float)*(V+2));

	cudaMemcpy(sDist, GetDist, sizeof(int)*(V*V+2),cudaMemcpyHostToDevice);
	cudaMemcpy(sSigma, GetSigma, sizeof(int)*(V*V+2),cudaMemcpyHostToDevice);
	cudaMemcpy(pc, perc, sizeof(float)*(V+2),cudaMemcpyHostToDevice);
	gpuErrchk( cudaPeekAtLastError() );

	percolation_estimate <<<NUM_BLOCKS, NUM_THREADS, 32>>> (V, E, dCentrality, sDist, sSigma, pc);
	cudaDeviceSynchronize();
	gpuErrchk( cudaPeekAtLastError() );
	
	cudaDeviceSynchronize();

	vector<pair<float,int> > perc_pair(V+1);
	vector<float> contrib(V+1);
	perc_pair[0].first = 0;
	perc_pair[0].second = 0;
    for(int i=1;i<=V;++i)
    {
		perc_pair[i].first = perc[i];
		perc_pair[i].second = i;
    }
	sort(perc_pair.begin(),perc_pair.end());
	float carry = 0,sum_x = 0;
	for(int i=1;i<=V;++i)
	{
		contrib[perc_pair[i].second] = (float)(i-1)*perc_pair[i].first-carry;
		carry += perc_pair[i].first;
		sum_x += contrib[perc_pair[i].second];
	}
	carry = 0;
	for(int i=V;i>=1;i--)
	{
		contrib[perc_pair[i].second] += carry-(float)(V-i)*perc_pair[i].first;
		carry += perc_pair[i].first;
	}

	float *Centrality = new float[V+2];
	cudaMemcpy(Centrality, dCentrality, sizeof(float)*(V+2), cudaMemcpyDeviceToHost);
	
	for(int i=1; i<=V; ++i)
	{
		printf("%f\n", Centrality[i]/(sum_x-contrib[i]));
	}
	
	delete[] hRow;
	delete[] hColumn;

	cudaFree(Queue);
	cudaFree(dRow);
	cudaFree(dColumn);
	cudaFree(Distance);
	return 0;
}


void printfGraph(vector <vector<int> > &graph)
{
	printf("Graph is\n");
	for(int i=1; i<graph.size(); ++i)
	{
		printf("\n%d\t", i);
		for(int j=0;j<(int)graph[i].size();++j)
		{
			int n = graph[i][j];
			printf("%d ", n);
		}
	}
}

__global__ void dPrintGraph(int V, int E, int *dRow, int *dColumn)
{
	printf("printing from device:\nRow is\n");
	for(int i=0; i<=V+1; ++i) printf("%d ", dRow[i]);
	printf("\nCol is\n");
	for(int i=0; i<2*E; ++i) printf("%d ", dColumn[i]);
	printf("\n");
}

__global__ void dPrintDist(int V, int E, int *dDistance)
{
	printf("Distances \n");
	for(int i=1; i<(V+1)*V; ++i)
	{
		if(i%(V+1)==0) { ++i; printf("\n"); }
		printf("%d ", dDistance[i]);
	}
}