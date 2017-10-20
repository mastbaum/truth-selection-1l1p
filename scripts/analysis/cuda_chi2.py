'''CUDA-accelerated chi2 calculation.

A work in progress, performing a Monte Carlo chi2 calculation for a toy model
using GPU acceleration via PyCUDA.

A. Mastbaum <mastbaum@uchicago.edu>, 2017/10
'''

import numpy as np
import pycuda.autoinit
import pycuda.driver as drv
import pycuda.curandom
from pycuda import gpuarray
from pycuda.compiler import SourceModule
import scipy.stats
from matplotlib import pyplot as plt

# Define CV and correlations
cv = np.array([50, 100, 75], dtype=np.float32)
cor = np.array([[1.00, 0.00, 0.00],
                [0.00, 1.00, 0.00],
                [0.00, 0.00, 1.00]], dtype=np.float32)

# Compute actual covariance matrix
cov = cor * np.sqrt(np.outer(cv, cv))
print(cov)

# Sample with np.random.multivariate_normal and reconstitute a sample
# covariance matrix
nu = 10000
samples = np.random.multivariate_normal(cv, cov, size=nu)
d = samples - cv
cov2 = np.sum(np.einsum('...i,...j->...ij', d, d), axis=0) / nu
print(cov2)

# Do our own Cholesky decomposition and sample
c = np.linalg.cholesky(cov)
print(c)

xs = np.empty_like(samples)
for i in range(nu):
    r = np.random.normal(size=c.shape[1])
    xs[i] = cv + np.sum(c * r, axis=1)

d = xs - cv
cov3 = np.sum(np.einsum('...i,...j->...ij', d, d), axis=0) / nu
print(cov3)

# CUDA kernels
code = '''
    #include <curand_kernel.h>

    const int nstates = %(NGENERATORS)s;
    __device__ curandState_t* states[nstates];

    extern "C" {

      __global__ void init_device_rngs(int seed) {
        int idx = threadIdx.x + blockIdx.x * blockDim.x;
        if (idx < nstates) {
          curandState_t* s = new curandState_t;
          if (s != 0) {
            curand_init(seed, idx, 0, s);
          }
          states[idx] = s;
        }
      }

      __device__ void sample(curandState* state,
                             const int n,
                             float* v,
                             float* cv,
                             float* dcov) {
        for (unsigned j=0; j<n; j++) {
          float r = curand_normal(state);
          float w = 0;
          for (unsigned int i=0; i<n; i++) {
            w += dcov[j * n + i] * r;
          }
          //v[j] = cv[j] - curand_poisson(state, cv[j] + w);
          //v[j] = cv[j] - curand_normal(state, cv[j] + w);
          //v[j] = cv[j] - cv[j] + w;
          v[j] = 1.0 * w;// + cv[j];
          //v[j] = curand_uniform(state);
        }
      }

      __global__ void chi2(const int ns,
                           const int n,
                           float* v,
                           float* cv,
                           float* icov,
                           float* dcov) {
        int idx = threadIdx.x + blockIdx.x * blockDim.x;
        if (idx < nstates) {
          curandState_t s = *states[idx];
          for(int i=idx; i<ns; i+=blockDim.x*gridDim.x) {
            float* pv = new float[n];
            sample(&s, n, pv, cv, dcov);
            //for (int j=0; j<n; j++) {
            //  v[n*i+j] = pv[j];
            //}
            //v[n*i+0] = pv[0];
            //v[n*i+1] = blockDim.x * gridDim.x;
            //v[n*i+2] = i;
            v[i] = 0;
            for (int ii=0; ii<n; ii++) {
              for (int jj=0; jj<n; jj++) {
                v[i] += pv[ii] * icov[jj * n + ii] * pv[jj];
              }
            }
          }
          *states[idx] = s;
        }
      }

    }
'''

print('\n')

N = 384
nblock = 32

bx = int(N / nblock)
gx = int((N + bx - 1) / bx)
print(

mod = SourceModule(code % {'NGENERATORS' : N}, no_extern_c=True, arch='sm_30')
init_func = mod.get_function('init_device_rngs')
fill_func = mod.get_function('chi2')

nvalues = int(nblock * N)
g_data = gpuarray.zeros(nvalues, dtype=np.float32)
g_cv = gpuarray.to_gpu(cv)
g_dcov = gpuarray.to_gpu(c)
g_icov = gpuarray.to_gpu(np.linalg.inv(cov))

init_func(np.int32(42),
          block=(192,1,1), grid=(1,1,1))


h, e = None, None

fill_func = mod.get_function('chi2')
fill_func(np.int32(nvalues), np.int32(len(cv)),
          g_data, g_cv, g_icov, g_dcov,
          block=(192,1,1),
          grid=(nblock,1,1))

data = g_data.get()
print(data)
print(data.shape)

if h is None:
    h, e = np.histogram(data, bins=100)
else:
    hh, _e = np.histogram(data, bins=e)
    h += hh

h = h.astype(np.float32) / np.sum(h) / np.diff(e)[0]

plt.bar(e[:-1]+np.diff(e)/2, h, width=np.diff(e)[0], label='MC')
plt.plot(e[:-1], scipy.stats.chi2.pdf(e[:-1], df=len(cv)),
                 color='black', label='$\chi^2(dof=N)$')
plt.yscale('log')
plt.show()

