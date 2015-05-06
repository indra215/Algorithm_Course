#include <iostream>
#include <list>
#include <vector>
#include <algorithm>

using namespace std;

bool iscycle = false;

class Graph
{
	int V;
	list<int> *adj;
	void get_cycle(int v,bool visited[],int parent[]);
	void print_cycle(int v,int u,int parent[]);

public:
	Graph(int V);
    void addEdge(int v, int w);
    void find_cycle();
};

Graph::Graph(int V)
{
	this->V = V;
    adj = new list<int>[V];
}

void Graph::addEdge(int v,int w)
{
	adj[v].push_back(w);  // Add w to v’s list.
	adj[w].push_back(v);  // Add v to w's list
}

void Graph::print_cycle(int u,int v,int parent[])
{
	cout << v << "-->";
	while(1){
		int n = parent[v];
		if(n == u)
			break;
		cout << n << "-->";
		v = n;
	}
	cout << u << endl;
}

void Graph::get_cycle(int v,bool visited[],int parent[])
{
	visited[v] = true;

	list<int>::iterator i;
    for(i = adj[v].begin(); i != adj[v].end(); ++i){
    	if(!visited[*i]){
        	parent[*i] = v;
            get_cycle(*i, visited, parent);
        }else if(visited[*i]==true && parent[v] != *i){
        	if(!iscycle){
        		cout << "Cycle Detected" << endl;
        		print_cycle(*i,v,parent);
        		iscycle = true;
        	}
        }
    }
}

void Graph::find_cycle()
{
	bool *visited = new bool[V];
	int *parent = new int[V];

	for(int i=0;i<V;i++)
	{
		visited[i] = false;
	}
	
	for(int i=0;i<V;i++)
	{
		if(visited[i] == false){
			get_cycle(i,visited,parent);
		}
	}
}

int main()
{
	// Testcase 1
	int V = 8;
	Graph g(V);
    g.addEdge(0,3);
    g.addEdge(0,4);
    g.addEdge(0,2);
    g.addEdge(1,0);
    g.addEdge(2,1);
    g.addEdge(2,3);
    g.addEdge(3,6);
    g.addEdge(4,5);
    g.addEdge(5,4);
    g.addEdge(6,3);
    g.addEdge(6,7);
    g.addEdge(7,7);

    // Testcase 2
    // Uncomment this to run for this testcase
    /*int V = 8;
	Graph g(V);
    g.addEdge(0,3);
    g.addEdge(0,4);
    g.addEdge(1,0);
    g.addEdge(1,3);
    g.addEdge(2,1);
    g.addEdge(3,2);
    g.addEdge(3,6);
    g.addEdge(4,5);
    g.addEdge(4,6);
    g.addEdge(6,4);
    g.addEdge(7,6);*/


    g.find_cycle();

	return 0;
}