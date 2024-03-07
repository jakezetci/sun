import numpy as np
import time
b = np.array([2,1,2])
print(b**2)
tic = time.perf_counter()
for i in range(1000000):
    a = np.inner(b,b)
toc = time.perf_counter()
print(f'{toc-tic:.2}')
tic = time.perf_counter()
for i in range(1000000):
    a = np.linalg.norm(b)**2
toc = time.perf_counter()
print(f'{toc-tic:.2}')

tic = time.perf_counter()
for i in range(1000000):
    a = np.sum(b**2)
toc = time.perf_counter()
print(f'{toc-tic:.2}')
