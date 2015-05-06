x	//MTH524A Assignment 1 part (a)
// This algorithm takes input R from the file "inputR.txt" 
// and then outputs all the points Rl i.e. R intersection Cl
// and writes the output to the file "output.txt".


#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <stdlib.h>
#include <vector>

using namespace std;

double a=1;				
int K=3;

//structure of coordinate where we store x,y,z coordinates of a particular point
struct coord{
	int x,y,z;
};

//PARTITION function where given a r,L pair, gives all possible l's for which r belongs to Cl
vector<struct coord> PARTITION(double r[],int L)
{
	vector<struct coord> l;
	for(int i=0;i<L;i++){
		for(int j=0;j<L;j++){
			for(int k=0;k<L;k++){
				if(i*a/L <= r[0] && (i+1)*a/L >= r[0]){
					if(j*a/L <= r[1] && (j+1)*a/L >= r[1]){
						if(k*a/L <= r[2] && (k+1)*a/L >= r[2]){
							struct coord new_coord;
							new_coord.x = i;
							new_coord.y = j;
							new_coord.z = k;
							l.push_back(new_coord);
						}
					}
				}
			}
		}
	}
	return l;
}

//PARTITION-ALL function where given a set R of points, gives a set Rl which is intersection of R and Cl 
bool partition_all(int L, char* inputFileName, char* outputFileName)
{
	int Rl[L][L][L][K];
	memset(Rl,0,sizeof(Rl));
	double R[K][3];
	string line;
	ifstream inpfile (inputFileName);
	int i=0;
	if(inpfile.is_open())
	{
		while (getline(inpfile, line))
    	{
	        double value;
	        int k = 0;
	        stringstream ss(line);
	
	        while (ss >> value)
	        {
	            R[i][k] = value;
	            ++k;
	        }
	        ++i;
    	}
	}
	K = i;
	for(int i=0;i<K;i++){
		cout << R[i][0] << R[i][1] << R[i][2] << endl;
	}
	inpfile.close();
	ofstream ofile (outputFileName);
	for(int i=0;i<K;i++){
		vector<struct coord> ord;
		ord = PARTITION(R[i],L);
		for(int j=0;j<ord.size();j++){
			Rl[ord[j].x][ord[j].y][ord[j].z][i] = 1;
		}
	}
	for(int i=0;i<L;i++){
		for(int j=0;j<L;j++){
			for(int k=0;k<L;k++){
				int ans = 0;
				for(int p=0;p<K;p++){
					if(Rl[i][j][k][p] == 1)
						ans++;
				}
				ofile << i << " " << j << " " << k << endl;
				ofile << ans << endl;
				for(int p=0;p<K;p++){
					if(Rl[i][j][k][p] == 1) {
							ofile << R[p][0] << " " << R[p][1] << " " << R[p][2] << endl;	
					}
				}
			}
		}
	}
	ofile.close();
	return true;
}

int main()
{
	int L=2;
	bool succ = partition_all(L,"inputR.txt","output.txt");
	return 0;
}
