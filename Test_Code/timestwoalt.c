#include "mex.h"

void timestwo_alt(double *y, double x)
{
    *y = 2.0*x;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *M;
    int m,n;
    M = mxGetPr(prhs[0]);
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    plhs[0] = mxCreateDoubleMatrix(m, n, mxINT32_CLASS, mxREAL);
    double *A;
    A = mxGetPr(plhs[0]);
    for(int i=0; i<m*n; i++)
    {
        timestwo_alt(A+i, *(M+i));
    }
}