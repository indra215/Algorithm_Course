// MTH524A Assignment 1 part (b)
// This algorithm takes input P from the file "inputP.txt", input Q from the file "inputQ.txt" 
// and then computes the matrix ATA and writes it to the file "outputATA.txt"


#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;

int N=1,M=2;
float a=1;

//function to check whether a given point lies on the boundary of Cl
bool checkValidCube(float point[],int L,int l[])
{
	if((l[0]*a/L <= point[0] && (l[0]+1)*a/L >= point[0])) 
		return true;
	return false;
}

//function to check whether a given point lies on the boundary of Sl
bool checkValidSphere(float point[],int l[],int L)
{
	float dist = pow((point[0] - (2*l[0]+1)*a/(2*L)),2) + pow((point[1] - (2*l[1]+1)*a/(2*L)),2) + pow((point[2] - (2*l[2]+1)*a/(2*L)),2);
	if(dist == pow((2*a/L),2))
		return true;
	return false;
}

//function to calculate ATA given set of points P and Q
float **createATA(float *P,float *Q)
{
	float **A;
	A = new float*[M];
	float num;
	for(int i=0;i<M;i++)
	{
		A[i] = new float[N];
		for(int j=0;j<N;j++)
		{
			num = sqrt(pow((*(Q+i*3)-*(P+j*3)),2) + pow((*(Q+i*3+1)-*(P+j*3+1)),2) + pow((*(Q+i*3+2)-*(P+j*3+2)),2));
			A[i][j] = sin(num)/num;
		}
	}
	float A_trans[N][M];
	float **ATA;
	ATA = new float*[N];
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<N;j++)
		{
			A_trans[j][i] = A[i][j];
		}
	}
	for(int i=0;i<N;i++)
    {
    	for(int j=0;j<N;j++)
        {
        	ATA[i] = new float[N];
        	ATA[i][j]=0;
            for(int k=0;k<M;k++)
            {
            	ATA[i][j]=ATA[i][j]+A_trans[i][k]*A[k][j];
        	}
        }
    }
    return ATA;
}

//function which generates ATA and writes it to the output file
bool matrix(int L, int l[3], char* fileForP, char* fileForQ,char* fileForMatrix)
{
	float P[N][3],Q[M][3];
	float temp[3],tempS[3];
	int K;
	string line;
	ifstream Pfile (fileForP);
	ifstream Qfile (fileForQ);
	int i=0,k=0;
	if(Pfile.is_open())
	{
		while (getline(Pfile, line))
    	{
	        float value;
	        int k = 0;
	        int yes=0;
	        stringstream ss(line);
	
	        while (ss >> value)
	        {
	        	temp[k] = value;
	        	if(checkValidCube(temp,L,l)){
	        		yes += 1;	
				}
	            ++k;
	        }
	        if(yes == 3){
	        	P[i][0] = temp[0];
	        	P[i][1] = temp[1];
	        	P[i][2] = temp[2];
	        	++i;
			}
    	}
    	N = i;
    	cout << N << endl;
		for(int i=0;i<N;i++)
		{
			cout << P[i][0] << " " << P[i][1] << " " << P[i][2] << endl;
		}
	}
	Pfile.close();
	i=0;k=0;
	if(Qfile.is_open())
	{
		while (getline(Qfile, line))
    	{
	        float value;
	        int k = 0;
	        stringstream ss(line);
	
	        while (ss >> value)
	        {
	            tempS[k] = value;
	            ++k;
	        }
	        if(checkValidSphere(tempS,l,L)){
	        	Q[i][0] = tempS[0];
	        	Q[i][1] = tempS[1];
	        	Q[i][2] = tempS[2];
				++i;
			}
    	}
		
	}
	M = i;
		for(int i=0;i<M;i++)
		{
			cout << Q[i][0] << " " << Q[i][1] << " " << Q[i][2] << endl;
		}
	float **ATA = createATA(P[0],Q[0]);	
	ofstream ofile (fileForMatrix);
	if(ofile.is_open())
	{
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
			{
				ofile << i << "\t" << j << "\t" << ATA[i][j];
				ofile << "\n";
			}
		}
	}
	return true;
}

int main()
{
	int L=2;
	int l[3]={0,0,0};
	bool succ = matrix(L,l,"inputP.txt","inputQ.txt","outputATA.txt");
	return 0;
}

