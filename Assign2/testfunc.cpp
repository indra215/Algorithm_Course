double* differentiate(double* coeff)
{ int i,j,k,x=0; /* Easier to deal with 3-D array than pointers therefore conversion */ 
for(i=x=0;i<=m;i++)
{ for(j=0;j<=m;j++)
{ for(k=0;k<=m;k++) 
cont1[i][j][k]=0.0,
cont[i][j][k]=coeff[x++]; 
} 
} 
for(i=0;i<=m;i++)
{ for(j=0;j<=m;j++)
{ for(k=0;k<=m;k++)
{ if(i>1) cont1[i-2][j][k]+=(cont[i][j][k]*i*(i-1)); 
if(j>1) cont1[i][j-2][k]+=(cont[i][j][k]*j*(j-1)); 
if(k>1) cont1[i][j][k-2]+=(cont[i][j][k]*k*(k-1)); } } } 
for(i=x=0;i<=m;i++)
{ for(j=0;j<=m;j++)
{ for(k=0;k<=m;k++) coeff[x++]=cont1[i][j][k]; } } 
return coeff; 
}

