#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <cstring>

using namespace std;

class Graph
{
	int V;
	list<int> *adj;
	void find_path(int v,bool visited[],int parent[]);
public:
	Graph(int V);
    void addEdge(int v, int w);
    void print_path(int s,int d);
};

Graph::Graph(int V)
{
	this->V = V;
    adj = new list<int>[V];
}

void Graph::addEdge(int v,int w)
{
	adj[v].push_back(w);  // Add w to v’s list.
}

void Graph::find_path(int v,bool visited[],int parent[])
{
	visited[v] = true;

	list<int>::iterator i;
    for(i = adj[v].begin(); i != adj[v].end(); ++i)
        if(!visited[*i]){
        	parent[*i] = v;
            find_path(*i, visited, parent);
        }
}

void Graph::print_path(int s,int d)
{
	bool *visited = new bool[V];
	int *parent = new int[V];
	for(int i=0;i<V;i++){
		visited[i] = false;
	}

	find_path(s,visited,parent);

	if(visited[d] == true){
		cout << "Path between the given vertices is" << endl;
		vector<int> path_s_d;
		//path_s_d.reserve(V);
		path_s_d.push_back(d);

		int temp = d;

		while(1){
			int n = parent[temp];
			path_s_d.push_back(n);
			if(n == s)
				break;
			temp = n;
		}

		reverse(path_s_d.begin(),path_s_d.end());

		for(vector<int>::iterator i=path_s_d.begin();i != path_s_d.end()-1;i++){
			cout << *i << "-->";
		}
		cout << *(path_s_d.end()-1) << endl;
	}else{
		cout << "There is no path between the given vertices" << endl;
	}
}

int main()
{
	int **edges,V,start,end;
    // Use this file for the firt maze testcase
    ifstream inp ("input_graph_1.txt");
    
    // Use this file for the second maze testcase and uncomment this line
    // ifstream inp ("input_graph_2.txt");
    if(inp.is_open())
    {
        inp >> V;
        V = sqrt(V);
        edges = new int*[V];
        cout << V << endl;
        for(int i=0;i<V;i++){
            edges[i] = new int[V];
            for(int j=0;j<V;j++)
                inp >> edges[i][j];
        }
        inp >> start >> end;
        inp.close();
    }
    
    Graph g(V*V);
    for(int i=0;i<V;i++)
    {
    	for(int j=0;j<V;j++)
    	{
    		if(edges[i][j] == 1){
    			if(i-1<V && i-1>-1){
    				if(edges[i-1][j] == 1){
    					g.addEdge(V*i+j,V*(i-1)+j);
    				}
    			}
    			if(i+1<V && i+1>-1){
    				if(edges[i+1][j] == 1){
    					g.addEdge(V*i+j,V*(i+1)+j);
    				}
    			}
    			if(j-1<V && j-1>-1){
    				if(edges[i][j-1] == 1){
    					g.addEdge(V*i+j,V*i+(j-1));
    				}
    			}
    			if(j+1<V && j+1>-1){
    				if(edges[i][j+1] == 1){
    					g.addEdge(V*i+j,V*i+(j+1));
    				}
    			}
    		}
    	}
    }

    g.print_path(start,end);
	return 0;
}