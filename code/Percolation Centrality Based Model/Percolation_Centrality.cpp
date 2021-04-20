#include<iostream> 
#include<vector>
#include<queue>

double R(double u){
	if(u>0) return u;
	return 0;
}

void bfs(int src, std::vector<int> &dist, std::vector<int> &sigma, std::vector<std::vector<int> > &g){
	dist[src] = 0;
	sigma[src] = 1;
	std::queue<int> q;
	q.push(src);
	while(!q.empty()){
		int u = q.front();
		q.pop();
		for(int v:g[u]){
			if(dist[v] < 0){
				q.push(v);
				dist[v] = dist[u] + 1;
			}
			if(dist[v] == dist[u] + 1){
				sigma[v] += sigma[u];
			}
		}
	}	
}

int main(){
	std::vector<std::vector<int> > g;
	std::vector<std::vector<int> > dist;
	std::vector<std::vector<int> > sigma;
	std::vector<double> perc;
	int n,m;
	std::cin>>n>>m;
	perc.resize(n);
	g.resize(n);
	sigma.resize(n,std::vector<int>(n));
	dist.resize(n,std::vector<int>(n,-1));
	for(int i=0;i<n;i++){
		perc[i] = 1.0/(double)(i+1);
	}
	for(int i=0;i<m;i++){
		int u,v;
		std::cin>>u>>v;
		u--;
		v--;
		g[u].push_back(v);
		g[v].push_back(u);
	}
	// std::cerr<<"ok"<<std::endl;

	const int INF = n+1;

	for(int i=0;i<n;i++){
		bfs(i,dist[i], sigma[i], g);
	}

	// for(int u=0;u<n;u++){
	// 	for(int v=u+1;v<n;v++){
	// 		std::cout<<u<<" "<<v<<" "<<sigma[u][v]<<std::endl;
	// 	}
	// }

	double denominator = 0;
	std::vector<double> PC(n);

	for(int s=0;s<n;s++){
		for(int r=0;r<n;r++){
			denominator += R(perc[s]-perc[r]);
		}
	}

	for(int v=0;v<n;v++){
		double denom_v = denominator;
		for(int s=0;s<n;s++){
			denom_v -= R(perc[v]-perc[s]);
			denom_v -= R(perc[s]-perc[v]);
			for(int r=0;r<n;r++){
				if(s==r || s==v || r==v) continue;
				// std::cout<<s<<" "<<r<<" "<<denominator<<std::endl;
				if(dist[s][r] == dist[s][v] + dist[v][r]){
					double ratio = sigma[s][v] * sigma[v][r];
					ratio /= sigma[s][r];
					PC[v] += ratio * R(perc[s]-perc[r]);
				}
			}
		}
		PC[v] /= denom_v;
	}

	for(int i=0;i<n;i++){
		std::cout<<PC[i]<<std::endl;
	}


	return 0;
}