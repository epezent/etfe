// MIT License

// Copyright (c) 2021 Evan Pezent

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once
#include <kiss_fftr.h>
#include <complex>
#include <vector>

namespace etfe {

/// PI
static constexpr double pi = 3.14159265358979323846264338327950288;
/// Complex type
typedef std::complex<double> complex;

//-----------------------------------------------------------------------------
// FFT
//-----------------------------------------------------------------------------

/// Fast Fourier Transfrom Wrapper (KISSFFT)
class FFT {
public:

    // Constructor
    FFT(std::size_t nfft) :
        m_nfft(nfft), 
        m_in(nfft),
        m_out(nfft/2+1)  // may need to be size nfft for other FFT packages
    { 
        m_cfg = kiss_fftr_alloc((int)nfft, 0, NULL, NULL);
    }

    // Destructor
    ~FFT() { free(m_cfg); }

    // Perform FFT Transform
    const std::vector<complex>& operator()(const double* in, std::size_t n) {
        kiss_fft_cpx* out = reinterpret_cast<kiss_fft_cpx*>(&m_out[0]);
        if (n < m_nfft) {
            std::fill(m_in.begin(), m_in.end(), 0);
            for (std::size_t i = 0; i < n; ++i)
                m_in[i] = in[i];
            kiss_fftr(m_cfg, m_in.data(), out);            
        }
        else {
            kiss_fftr(m_cfg, in, out);
        }
        return m_out;
    }

private:
    std::size_t          m_nfft; ///< FFT size
    std::vector<double>  m_in;   ///< local padded input buffer for when n < m_nfft
    std::vector<complex> m_out;  ///< FFT output
    kiss_fftr_cfg        m_cfg;  ///< KISSFFT plan

};

//-----------------------------------------------------------------------------
// WINDOWS
//-----------------------------------------------------------------------------

/// Generates hamming window of size W
inline std::vector<double> hamming(std::size_t W) {
    std::vector<double> window(W);
    for (std::size_t w = 0; w < W; ++w)
        window[w] = (0.54 - 0.46 * std::cos(2*pi*w/(W-1)));
    return window;
}

// Generates hann window of size W
inline std::vector<double> hann(std::size_t W) {
    std::vector<double> window(W);
    for (std::size_t w = 0; w < W; ++w)
        window[w] = (0.5 - 0.5 * std::cos(2*pi*w/(W-1)));
    return window;
}

// Generates winrect window of size W
inline std::vector<double> winrect(std::size_t W) {
    return std::vector<double>(W,1.0);
}

//-----------------------------------------------------------------------------
// ETFE
//-----------------------------------------------------------------------------

// Empirical Transfer-Function Estimate
class ETFE {
public:   

    /// Results from estimating the experimental transfer function
    class Result {
    public:

        Result(std::size_t n) :  f(n), psdx(n), psdy(n), csdxy(n), txy(n),  mag(n), phase(n), ampx(n), ampy(n), m_size(n) { }
        Result()                   { }
        std::size_t size() const   {return m_size; }

        std::vector<double>  f;     ///< freqeuncies [Hz]
        std::vector<double>  psdx;  ///< power spectral density of input x
        std::vector<double>  psdy;  ///< power spectral density of output y
        std::vector<complex> csdxy; ///< cross power spectral density of input x and output y
        std::vector<complex> txy;   ///< transfer function estimate for input x and output y (i.e. conj(csdxy) ./ psdx)
        std::vector<double>  mag;   ///< transfer function magnitude (i.e 20*log10(abs(txy)))
        std::vector<double>  phase; ///< transfer function phase (i.e. 180/pi*angle(txy))
        std::vector<double>  ampx;  ///< amplitudes of x
        std::vector<double>  ampy;  ///< amplitudes of y

    private:
        friend class ETFE;
        void reset() {
            std::fill(ampx.begin(), ampx.end(), 0);
            std::fill(ampy.begin(), ampy.end(), 0);
            std::fill(psdx.begin(), psdx.end(), 0);
            std::fill(psdy.begin(), psdy.end(), 0);
            std::fill(csdxy.begin(), csdxy.end(), complex(0,0));
        }
        std::size_t m_size;
    };

    /// Constructor 1
    ETFE(std::size_t nsamples, double fs) :
        ETFE(nsamples, fs, hamming(2*nsamples/(8+1)), nsamples/(8+1), std::max(256, (int)std::pow(2,std::ceil(std::log2(2*nsamples/(8+1))))))
    { }

    /// Constructor 2
    ETFE(std::size_t nsamples, double fs, const std::vector<double>& window, std::size_t noverlap, std::size_t nfft) :
        N(nsamples), 
        W(window.size()),
        O(noverlap), 
        K(O == 0 ? N / W : (O - W + N) / O),
        Fs(fs), 
        dF(fs/nfft),
        Nfft(nfft),
        m_fftx(nfft), 
        m_ffty(nfft),
        m_in(W),
        m_r(nfft/2+1)
    {         
        // generate frequencies
        for (std::size_t i = 0; i < Nfft/2+1; ++i)
            m_r.f[i] = i * dF;   
        // set window
        setWindow(window);
    }

    /// Estimate the experimental transfer function given input x and output y, each of length n provided to constructor
    const Result& operator()(const double* x, const double* y) {
        // reset results
        m_r.reset();
        // iterate window segments
        for (std::size_t k = 0; k < K; ++k) {

            // window kth segment of x
            for (std::size_t w = 0; w < W; ++w)
                m_in[w] = x[O*k+w] * m_win[w];
            // fft kth segment of x
            auto& woutx = m_fftx(m_in.data(), W);

            // window kth segment of y
            for (std::size_t w = 0; w < W; ++w)
                m_in[w] = y[O*k+w] * m_win[w];
            // fft kth segment of y
            auto& wouty = m_ffty(m_in.data(), W);

            // calc window contributions to average estimates
            for (std::size_t i = 0; i < Nfft/2+1; ++i)  {

                // power spectral densities
                const double psdScale = (1 + (i > 0 && i < Nfft/2+1)) * m_pcf / K;
                m_r.psdx[i]  += std::pow(std::abs(woutx[i]),2) * psdScale;          
                m_r.psdy[i]  += std::pow(std::abs(wouty[i]),2) * psdScale;  
                m_r.csdxy[i] += woutx[i] * std::conj(wouty[i]) * psdScale;      

                // amplitudes
                const double ampScale = 1.0 / (Nfft/2) * (1.0/K) * m_acf; 
                m_r.ampx[i] += std::abs(woutx[i]) * ampScale;
                m_r.ampy[i] += std::abs(wouty[i]) * ampScale;
            }
        }

        // compute ETFE
        for (std::size_t i = 0; i < Nfft/2+1; ++i) {
            m_r.txy[i]   = std::conj(m_r.csdxy[i]) / m_r.psdx[i];
            m_r.mag[i]   = 20*log10(std::abs(m_r.txy[i]));
            m_r.phase[i] = 180/pi * std::arg(m_r.txy[i]);
        }
 
        // return result
        return m_r;
    }

    /// Get the most recent result
    const Result& result() const {
        return m_r;
    }

    /// Set windowing type
    void setWindow(const std::vector<double>& window) {
        W     = window.size();
        m_win = window;
        double sum   = 0;
        double sumsq = 0;
        for (auto& w : m_win) {
            sum   += w;
            sumsq += w*w;
        }
        m_acf = m_win.size() / sum;
        m_pcf = 1.0 / (Fs * sumsq); // may need to be Nfft
    }

private:

public:
    std::size_t    N;    ///< number of samples expected in u and y
    std::size_t    W;    ///< number of samples per window segment
    std::size_t    O;    ///< number of overlapping samples between segments
    std::size_t    K;    ///< number of segments
    double         Fs;   ///< sampling frequency
    double         dF;   ///< frequency resolution
    std::size_t    Nfft; ///< FFT size

private:
    FFT                  m_fftx;  ///< FFT object for x
    FFT                  m_ffty;  ///< FFT object for y
    std::vector<double>  m_win;   ///< window coefficients
    double               m_acf;   ///< amplitude correction factor
    double               m_pcf;   ///< power correction factor
    std::vector<double>  m_in;    ///< windowed input buffer
    Result               m_r;     ///< result buffers
};

} // namespace etfe

//-----------------------------------------------------------------------------
// RESOURES
//-----------------------------------------------------------------------------

// Reproducing MATLAB's pwelch
// https://www.mathworks.com/matlabcentral/answers/29641-quetsions-about-matlab-pwelch-implementation
// https://www.gaussianwaves.com/2015/11/interpreting-fft-results-complex-dft-frequency-bins-and-fftshift/
// https://stackoverflow.com/questions/14536950/applying-kiss-fft-on-audio-samples-and-getting-nan-output
// http://www.iowahills.com/FFTCode.html
// https://community.sw.siemens.com/s/article/window-correction-factors
// https://www.mathworks.com/matlabcentral/answers/372516-calculate-windowing-correction-factor