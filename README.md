# etfe

`ETFE.hpp` emulates MATLAB's `tfestimate`, `pwelch`, and `cpsd` functions. It calculates the experimental transfer function estimate between input x and output y `txy`, the power spectral densities `pxx` and `pyy`, and the cross spectral density `pxy`. By default, it behaves exactly as MATLAB's functions, and similarly can be provided with specified windows, overlap, and FFT sizes. The output has been tested to precisely match that of MATLAB's (see `matlab\test.m`). 

**Note:** *Currently, only real-valued 1D sample inputs are supported*.

## Usage

MATLAB
```matlab
[pxx,f] = pwelch(x,[],[],[],fs);
[pxy,f] = cpsd(x,y,[],[],[],fs);
[txy,f] = tfestimate(x,y,[],[],[],fs);
mag     = 20*log10(abs(txy));
phase   = 180/pi*angle(txy);
```

ETFE.hpp
```cpp
ETFE etfe(nsamples,fs);
ETFE::Result& result = etfe.estimate(x,y);

result.f;
result.txy;
result.pxx;
result.pxy;
result.mag;
result.phase;
```

## Example Apps

- `etfe_cl.cpp`    - command line interface to compute ETFE from CSV files
- `filter_toy.cpp` - interactive spectrum and Bode plots for filters

![filter_toy](https://user-images.githubusercontent.com/29577475/104275075-c9446700-5467-11eb-8d86-19c2688b1951.gif)

## Dependencies

- [KISS FFT](https://github.com/mborgerding/kissfft) - other FFT packages should be easily substitutable 
- [mahi-gui](https://github.com/mahilab/mahi-gui) - for filter demo only
- [IIR](https://github.com/berndporr/iir1) - for filter demo only

## Resources

https://www.mathworks.com/matlabcentral/answers/29641-quetsions-about-matlab-pwelch-implementation
https://www.gaussianwaves.com/2015/11/interpreting-fft-results-complex-dft-frequency-bins-and-fftshift/
https://stackoverflow.com/questions/14536950/applying-kiss-fft-on-audio-samples-and-getting-nan-output
http://www.iowahills.com/FFTCode.html
https://community.sw.siemens.com/s/article/window-correction-factors
https://www.mathworks.com/matlabcentral/answers/372516-calculate-windowing-correction-factor
