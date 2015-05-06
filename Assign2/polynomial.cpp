#include <iostream>
#include <cmath>
#include <cstring>

using namespace std;

double gradient_value(double *coeff,int m,int iter,double x[])
{
	double new_coeff[m+1][m+1][m+1];
	if(iter == 0){
		int ind=0;
		double coeff_p[m+1][m+1][m+1];
		for(int i=0;i<=m;i++){
			for(int j=0;j<=m;j++){
				for(int k=0;k<=m;k++){
					coeff_p[i][j][k] = coeff[ind++];
				}
			}
		}
		double grad_value = 0;
		for(int i=0;i<=m;i++){
			for(int j=0;j<=m;j++){
				for(int k=0;k<=m;k++){
						grad_value += coeff_p[i][j][k]*pow(x[0],i)*pow(x[1],j)*pow(x[2],k);
				}
			}
		}
		return grad_value;
	}
	for(int i=0;i<iter;i++){
		int ind=0;
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
	}
	double grad_value = 0;
	for(int i=0;i<=m;i++){
		for(int j=0;j<=m;j++){
			for(int k=0;k<=m;k++){
					grad_value += new_coeff[i][j][k]*pow(x[0],i)*pow(x[1],j)*pow(x[2],k);
			}
		}
	}
	return grad_value;
}

int abk(int a,int b,int k)
{
	int abk = 1;
	if(k==0){
		return 1;
	}
	for(int i=0;i<k;i++)
	{
		abk *= (a+i*b);
	}
	return abk;
}

double mod_x_power(double x[],int k)
{
	double sum=0;
	for(int i=0;i<3;i++){
		sum += pow(abs(x[i]),2*k);	
	}
	return sum;
}

double compute_q(int m,double *cp,double *x)
{
	double mod_x2 = 0;
	int n = 3;
	mod_x2 = mod_x_power(x,1);
	double q_x = 0;
	for (int k = 0; k<=m/2; k++)
	{
		q_x += (pow(-1,k)*pow(mod_x2,k)*gradient_value(cp,m,k,x))/((abk(2,2,k+1)*abk(n+2*m-2*k,2,k+1)));
	}
	q_x *= mod_x2;
	return q_x;
}

int main()
{
	//degree of the polynomial
	int m=4;
	int size_arr = pow(m+1,3);
	//coefficients of the polynomial p(x)
	double co_effp[size_arr];
	memset(co_effp,0,sizeof(co_effp));
	//point of evaluation
	double x[3]={1,1,2};
	cout << "Value of Q(x) at given x is " << endl;
	cout << compute_q(m,co_effp,x) << endl;
	return 0;
}
