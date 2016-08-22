/*
 *  Usage: [W H objKL timeKL] = ccd_KL(V, k, max_iter, Winit, Hinit, trace);
 *
 * Given the nonnegative input matrix V, this code solves the following KL-NMF problem to find the low-rank approximation WH for V. 
 *
 *  min_{W>=0,H>=0} sum_{i,j} V_{ij}*log(V_{ij}/(WH)_{ij})
 *
 *  Input arguments
 *  	V: n by m nonnegative input matrix.
 *  	k: rank of output matrices W and H. 
 *  	max_iter: maximum iteration. 
 *  	Winit: k by n initial matrix for W. 
 *  	Hinit: k by m initial matrix for H. 
 *  	trace: 1: compute objective value per iteration. 
 *  		   0: do not compute objective value per iteration. (default)
 *
 *  Output arguments
 *  	W: k by n dense matrix.
 *  	H: k by m dense matrix.
 *  	objKL: objective values.
 *  	timeKL: time taken by this algorithm. 
 *
 */

#include "math.h"
#include "mex.h" 
#include <time.h>

int gradDescent(double *neighVecs, int numNeighs, double *currentCoord, int k, double *sum, int maxIter, double eps)
{
	double normCurrentVec = 0;
    double *tmp_i = (double *)malloc(sizeof(double)*k);
    double normTmp_i = 0;
    
    for(int i = 0; i < k ; i++)
    {
        tmp_i[i]=0;
    }
    for(int iter = 0; iter < maxIter ; iter++)
    {
//         normCurrentVec = 0;
//         normTmp_i = 0;
//         for(int i = 0; i < k ; i++)
//         {
//             normCurrentVec += pow(currentCoord[i]-tmp_i[i],2);
//             normTmp_i += pow(tmp_i[i],2);
//         }
//         normCurrentVec = sqrt(normCurrentVec);
//         normTmp_i = sqrt(normTmp_i);
//         if(normCurrentVec/normTmp_i < eps)
//         {
//             break;
//         }
        
        int maxinner = 2;
        double g, h, s;
        for ( int q=0 ; q<k ; q++ )
        {
            g=0; h=0; s=0;
            for (int inneriter =0 ; inneriter<maxinner ; inneriter++)
            {
                double oldCurQ, newCurQ, diff;
                g = sum[q];
                for(int i = 0; i < numNeighs ; i++)
                {
                    double dotProd = 0;
                    for(int j = 0; j < k ; j++)
                    {
                        dotProd += neighVecs[j*numNeighs+i]*currentCoord[j];
                    }
                    dotProd += s*neighVecs[q*numNeighs+i];
                    g -= neighVecs[q*numNeighs+i]*(1/dotProd);
                    h += (neighVecs[q*numNeighs+i]*neighVecs[q*numNeighs+i])*(1/(dotProd*dotProd));
                }
                s = -g/h;
                oldCurQ=currentCoord[q];
                newCurQ=currentCoord[q] + s;
                if(newCurQ < 1e-15)
                    newCurQ = 1e-15;
                diff = newCurQ-oldCurQ;
                currentCoord[q] = newCurQ;
                if ( fabs(diff) < fabs(oldCurQ)*0.5 )
                    break;
            }
        }
        for(int j = 0; j < k ; j++)
        {
            tmp_i[j] = currentCoord[j];
        }
    }
//     for(int j = 0; j < k ; j++)
//     {
//         resNorm[j] = 0;
//     }
    for(int i = 0; i < numNeighs ; i++)
    {
        double dotProd = 0;
        for(int j = 0; j < k ; j++)
        {
            dotProd += neighVecs[j*numNeighs+i]*currentCoord[j];
        }
//         for(int j = 0; j < k ; j++)
//         {
//             resNorm[j] -= neighVecs[j*numNeighs+i]*(1/dotProd);
//         }
    }


    free(tmp_i);
	return 0;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *neighVecs;
    int numNeighs;
    double *currentCoord;
    double *sum;
    int maxIter;
    double eps;
//     double *resNorm;
    

	// Check input/output number of arguments
	if ( nlhs > 3 || nrhs != 5)
	{
		printf("--Number of input or output arguments are not correct.\n");
		return;
	}

	neighVecs = mxGetPr(prhs[0]);
	numNeighs = mxGetM(prhs[0]);
    
    currentCoord = mxGetPr(prhs[1]);
    int k  =mxGetN(prhs[1]);
    sum = mxGetPr(prhs[2]);
    
	maxIter = mxGetScalar(prhs[3]);
	eps = mxGetScalar(prhs[4]);

    double *outCurrentCoord;
    plhs[0] = mxCreateDoubleMatrix(1,k,mxREAL);
    outCurrentCoord = mxGetPr(plhs[0]);
    
    double *outResNorm;
    plhs[1] = mxCreateDoubleMatrix(1,k,mxREAL);
    outResNorm = mxGetPr(plhs[1]);
    int outTisValCount;
    
	gradDescent(neighVecs, numNeighs, currentCoord, k, sum, maxIter, eps);
    
    for(int j = 0; j < k ; j++)
    {
        outCurrentCoord[j] = currentCoord[j];
    }
    
    plhs[2] = mxCreateDoubleScalar((double)outTisValCount);

	return;
}
