#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

// Dimensions of each of the set of points
int K=1,N=1,M=2,L=2;
float a=1;

// All the matrices
float **R,**P,**Q,*Gamma,*F_q,**A,**A_trans,**ATA,*ATF_l_q;
// ATF_l_q = new float[N];
// ATA = new float*[N];
// A_trans = new float*[N];
// A = new float*[M];
// R = new float*[K];
// P = new float*[N];
// Q = new float*[M];
// gamma = new float*[K];
// F_q = new float[M];

float calc_norm(float *p,float *q)
{
    float norm;
    for(int i=0;i<3;i++){
        norm += pow((p[i]-q[i]),2);
    }
    norm = sqrt(norm);
    return norm;
}

vector<int> check_points_l(float **R,int i,int j,int k)
{
    vector<int> indices;
    for(int p=0;p<K;i++){
        if(i*a/L <= R[p][0] && (i+1)*a/L >= R[p][0] && j*a/L <= R[p][1] && (j+1)*a/L >= R[p][1] && k*a/L <= R[p][2] && (k+1)*a/L >= R[p][2]){
            indices.push_back(p);
        }
    }
    return indices;
}

void compute_ATF_l_q(vector<int> indices,int i,int j,int k)
{
    F_q = new float[M];
    float num,f=0,q[3];
    for(int p=0;p<M;p++){
        q[0] = Q[p][0]+(i/L);
        q[1] = Q[p][1]+(j/L);
        q[2] = Q[p][2]+(k/L);
        for(int s=0;s<indices.size();s++){
            num = calc_norm(q,R[indices[s]]);
            float g = Gamma[indices[s]];
            f += (sin(num)/num)*g;
        }
        F_q[p] = f;
    }
    ATF_l_q = new float[N];
    for(int i=0;i<N;i++){
        for(int j=0;j<M;j++){
            ATF_l_q[i] += A_trans[i][j]*F_q[j]; 
        }
    }
}

double *solve_LS(float **B,float *b) 
{
    int n = N;
    float **A;
    A = new float*[n];
    for(int i=0;i<n;i++){
        A[i] = new float[n+1];
        for(int j=0;j<n;j++){
            A[i][j] = B[i][j];
        }
    }
    for(int j=0;j<n;j++){
            A[j][n] = b[j];
        }
    for (int i=0; i<n; i++) {
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }
        for (int k=i; k<n+1;k++) {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }
        for (int k=i+1; k<n; k++) {
            double c = -A[k][i]/A[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }
    double *x;
    x = new double[n];
    for (int i=n-1; i>=0; i--) {
        x[i] = A[i][n]/A[i][i];
        for (int k=i-1;k>=0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    return x;
}

void PARTITION_ALL(float **R,int L,char* outputFileName)
{
    vector<int> indices;
    ofstream out (outputFileName);
    double *solution;
    for(int i=0;i<L;i++){
        for(int j=0;j<L;j++){
            for(int k=0;k<L;k++){
                indices = check_points_l(R,i,j,k);
                if(indices.size()){
                    compute_ATF_l_q(indices,i,j,k);
                    solution = solve_LS(ATA,ATF_l_q);
                    out << i << " " << j << " " << k << endl;
                    for(int i=0;i<N;i++){
                        out << solution[i] << " ";
                    }
                    out << endl;
                }
            }
        }
    }
}

void **buildmatrix(float **P,float **Q)
{
    ATA = new float*[N];
    A_trans = new float*[N];
    A = new float*[M];
    float num;
    for(int i=0;i<M;i++)
    {
        A[i] = new float[N];
        for(int j=0;j<N;j++)
        {
            num = calc_norm(P[i],Q[j]);
            A[i][j] = sin(num)/num;
        }
    }
    for(int i=0;i<M;i++)
    {
        A_trans[i] = new float[M];
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
}

bool solve(int L, char* inputFileR, char* inputFileP, char* inputFileQ, char* inputFileGamma, char* outputFile)
{
    ifstream Rfile (inputFileR);
    ifstream Pfile (inputFileP);
    ifstream Qfile (inputFileQ);
    ifstream Gammafile (inputFileGamma);

    //float **P,**Q,**R,**gamma;
    R = new float*[K];
    P = new float*[N];
    Q = new float*[M];
    Gamma = new float[K];

    if(Rfile.is_open() || Pfile.is_open() || Qfile.is_open() || Gammafile.is_open()){
        for(int i=0;i<K;i++){
            R[i] = new float[3];
            Rfile >> R[i][0] >> R[i][1] >> R[i][2];
        }
        for(int i=0;i<K;i++){
            Gammafile >> Gamma[i];
        }
        for(int i=0;i<N;i++){
            P[i] = new float[3];
            Pfile >> P[i][0] >> P[i][1] >> P[i][2];
        }
        for(int i=0;i<M;i++){
            Q[i] = new float[3];
            Qfile >> Q[i][0] >> Q[i][1] >> Q[i][2];
        }
        PARTITION_ALL(R,L,outputFile);
        buildmatrix(P,Q);
    }else{
        return false;
    }


    return true;
}

int main()
{
    int L=2;
    bool succ = solve(L,"inputR.txt","inputP.txt","inputQ.txt","inputGamma.txt","output.txt");
    if(succ){

    }else{
        cout << "Error in the input files::File Empty" << endl;
    }
    return 0;
}
