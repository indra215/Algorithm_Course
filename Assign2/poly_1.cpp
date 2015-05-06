#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>

using namespace std;

double *gradient(double *coeff,int m)
{
	int ind=0;
	double new_coeff[m+1][m+1][m+1];
	double coeff_p[m+1][m+1][m+1];
	memset(new_coeff,0,sizeof(new_coeff));
	for(int i=0;i<=m;i++){
		for(int j=0;j<=m;j++){
			for(int k=0;k<=m;k++){
				coeff_p[i][j][k] = coeff[ind++];
			}
		}
	}
	for(int i=0;i<=m;i++){
		for(int j=0;j<=m;j++){
			for(int k=0;k<=m;k++){
				if(i > 1)
					new_coeff[i-2][j][k] += coeff_p[i][j][k] * i * i-1;
				if(j > 1)
					new_coeff[i][j-2][k] += coeff_p[i][j][k] * j * j-1;
				if(k > 1)
					new_coeff[i][j][k-2] += coeff_p[i][j][k] * k * k-1;	
			}
		}
	}
	ind = 0;
	for(int i=0;i<=m;i++){
		for(int j=0;j<=m;j++){
			for(int k=0;k<=m;k++){
				coeff[ind++] = new_coeff[i][j][k];
			}
		}
	}
	return coeff;
}

int main()
{
	
	return 0;
}
