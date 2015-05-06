#include <iostream>
#include <list>
#include <vector>
#include <algorithm>

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
	adj[w].push_back(v);  // Add v to w's list
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
		path_s_d.reserve(V);
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
	// Testcse 1
	int V = 8;
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
    g.addEdge(7,6);

    // Testcse 2
    // Uncomment to run for this testcase
    /*int V = 5
	Graph g(V);
    g.addEdge(1, 0);
    g.addEdge(0, 2);
    g.addEdge(2, 1);
    g.addEdge(0, 3);
    g.addEdge(3, 4);*/


    int s,d;

    cout << "Enter the source and destination vertices " << endl;
    cin >> s >> d;
    g.print_path(s,d);

	return 0;
}