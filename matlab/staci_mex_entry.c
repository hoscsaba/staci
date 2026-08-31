#include "mex.h"

void staciMexFunction(int nlhs, mxArray *plhs[], int nrhs,
                      const mxArray *prhs[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    staciMexFunction(nlhs, plhs, nrhs, prhs);
}
