int q_x(double  x[],int m, double cp[])
{
	mod_x2 = 0;
	int n = 3;
	for (int i = 0; i<n; i++)
	{
		mod_x2 += x[i]*x[i]; 
	}
	q_x = 0;
	int n = 3;
	for (int k = 0; k<=m/2; k++ )
	{
		q_x += mod_x2*(pow(-1,k)*pow(mod_x2, k)*delta_p(k)/(abk[2,2,k]*abk[n+2m-2k,2,k]))  // Indra find value of delta
	}
}
