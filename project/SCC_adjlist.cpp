// Implementation of Kosaraju's algorithm to print all SCCs
#include <iostream>
#include <list>
#include <stack>

using namespace std;
 
class Graph
{
    int V;    // No. of vertices
    list<int> *adj;    // An array of adjacency lists
 
    // Fills Stack with vertices (in increasing order of finishing times)
    // The top element of stack has the maximum finishing time
    void fill_vertex_finishtime(int v, bool visited[], stack<int> &Stack);
 
    // A recursive function to print DFS starting from v
    void DFS_invert_graph(int v, bool visited[]);
public:
    Graph(int V);
    void addEdge(int v, int w);
 
    // The main function that finds and prints strongly connected components
    void print_SCCs();
 
    // Function that returns reverse (or transpose) of this graph
    Graph invert_graph();
};
 
Graph::Graph(int V)
{
    this->V = V;
    adj = new list<int>[V];
}
 
// A recursive function to print DFS starting from v
void Graph::DFS_invert_graph(int v, bool visited[])
{
    // Mark the current node as visited and print it
    visited[v] = true;
    cout << v << " ";
 
    // Recur for all the vertices adjacent to this vertex
    list<int>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
        if (!visited[*i])
            DFS_invert_graph(*i, visited);
}
 
Graph Graph::invert_graph()
{
    Graph g(V);
    for (int v = 0; v < V; v++)
    {
        // Recur for all the vertices adjacent to this vertex
        list<int>::iterator i;
        for(i = adj[v].begin(); i != adj[v].end(); ++i)
        {
            g.adj[*i].push_back(v);
        }
    }
    return g;
}
 
void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w); // Add w to v’s list.
}
 
void Graph::fill_vertex_finishtime(int v, bool visited[], stack<int> &Stack)
{
    // Mark the current node as visited and print it
    visited[v] = true;
 
    // Recur for all the vertices adjacent to this vertex
    list<int>::iterator i;
    for(i = adj[v].begin(); i != adj[v].end(); ++i)
        if(!visited[*i])
            fill_vertex_finishtime(*i, visited, Stack);
 
    // All vertices reachable from v are processed by now, push v to Stack
    Stack.push(v);
}
 
// The main function that finds and prints all strongly connected components
void Graph::print_SCCs()
{
    stack<int> Stack;
 
    // Mark all the vertices as not visited (For first DFS)
    bool *visited = new bool[V];
    for(int i = 0; i < V; i++)
        visited[i] = false;
 
    // Fill vertices in stack according to their finishing times
    for(int i = 0; i < V; i++)
        if(visited[i] == false)
            fill_vertex_finishtime(i, visited, Stack);
 
    // Create a reversed graph
    Graph gr = invert_graph();
 
    // Mark all the vertices as not visited (For second DFS)
    for(int i = 0; i < V; i++)
        visited[i] = false;
 
    // Now process all vertices in order defined by Stack
    while (Stack.empty() == false)
    {
        // Pop a vertex from stack
        int v = Stack.top();
        Stack.pop();
 
        // Print Strongly connected component of the popped vertex
        if (visited[v] == false)
        {
            gr.DFS_invert_graph(v, visited);
            cout << endl;
        }
    }
}
 
// Driver program to test above functions
int main()
{
    // Create a graph given in the above diagram
    // Testcase 1
    int V = 5
    Graph g(V);
    g.addEdge(1, 0);
    g.addEdge(0, 2);
    g.addEdge(2, 1);
    g.addEdge(0, 3);
    g.addEdge(3, 4);


    // Testcase 2
    // Uncomment this to run for this case
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

    // Testcase 3
    // Uncomment this to run for this case
    /*int V = 8;
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
    g.addEdge(7,7);*/

 
    cout << "Following are strongly connected components in given graph \n";
    g.print_SCCs();
 
    return 0;
}