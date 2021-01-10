# etfe

The class in `ETFE.hpp` emulates output from MATLAB's `tfestimate`, `pwelch`, and `cpsd` functions. It calculates the experimental transfer function estimate `txy`, the power spectral density of the input x and output y `psdx` and `psdy`, and the cross spectral density of x and y `csdxy`. By default, it behaves exactly as MATLAB's functions, but similarly can be provided with specified FFT windows, overlaps, and FFT sizes.

