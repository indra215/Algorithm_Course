#include <iostream>
#include <vector>

using namespace std;

int main()
{
	vector <float> f,f1;
	vector <float> g[2];
	f.push_back(0.3);
	f.push_back(0.2);
	f.push_back(0.1);
	f1.push_back(0.4);
	f1.push_back(0.2);
	f1.push_back(0.1);
	g[0].push_back(f);
	//g[1].push_back(f1);
	for(int i=0;i<=2;i++){
		cout << g[0][i] << endl;
	}
	return 0;
}
