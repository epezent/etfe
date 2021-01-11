#include <random>
#include <Mahi/Gui.hpp>
#include <Mahi/Util.hpp>
#include <kiss_fftr.h>
#include <Iir.h>
#include "ETFE.hpp"

using namespace mahi::gui;
using namespace mahi::util;
using namespace etfe;

class FilterToy : public Application {
public:

    FilterToy() : 
        Application()
    { 
        ImPlot::StyleColorsDark();
    }

    void update() override {
        static bool init = true;
        static bool open = true;
        constexpr int N = 10000;
        constexpr double Fs = 1000;
        constexpr double dt = 1.0 / Fs;
        static double f1 = 3;
        static double f2 = 7;
        static double a1 = 1;
        static double a2 = 1;
        static double ng = 1;
        static double Fc = 250;
        static std::vector<double> noise(N);
        static Iir::Butterworth::LowPass<2> butt;
        static std::default_random_engine generator;
        static std::normal_distribution<double> distribution(0,1);

        if (init) {
            // generate static noise
            for (auto& n : noise)
                n = distribution(generator);
            // initialize filter
            butt.setup(Fs,Fc);
            init = false;
        }

        // gui inputs
        ImGui::Begin("Filter Toy",&open);
        ImGui::Text("x(t) = a1*sin(2*pi*f1*t) + a2*sin(2*pi*f2*t) + noise");
        ImGui::Separator();
        ImGui::DragDouble("f1 [Hz]",&f1,0.5f,0,500);
        ImGui::DragDouble("a1",&a1,0.1f,0,10);
        ImGui::DragDouble("f2 [Hz]",&f2,0.5f,0,500);
        ImGui::DragDouble("a2",&a2,0.1f,0,10);
        ImGui::DragDouble("noise",&ng,0.01f,0,10);
        if (ImGui::DragDouble("Fc [Hz]",&Fc,0.5f,0,500))
            butt.setup(Fs,Fc);
        ImGui::Separator();
        // reset filter
        butt.reset();
        // generate waveform
        static double x[N];
        static double y[N];
        static double t[N];
        for (int i = 0; i < N; ++i) {            
            t[i] = i * dt ;
            x[i] = a1*sin(2*PI*f1*t[i]) + a2*sin(2*PI*f2*t[i]) + ng * noise[i];
            y[i] = butt.filter(x[i]);
        }

        // plots
        ImPlot::SetNextPlotLimits(0,10,-10,10);
        if (ImPlot::BeginPlot("##Filter","Time [s]","Signal")) {
            ImPlot::PlotLine("x(t)", t, x, N);
            ImPlot::PlotLine("y(t)", t, y, N);
            ImPlot::EndPlot();
        }

        // perform ETFE
        static ETFE etfe(N, Fs, hamming(2000), 1000, 2048);
        auto& result = etfe.estimate(x,y);    

        ImGui::Separator();
        static int mode = 0;
        ImGui::ModeSelector(&mode, {"Amplitude","Power","Magnitude","Phase"});
        if (mode == 0) {
            ImPlot::SetNextPlotLimits(0,500,0,1);
            if (ImPlot::BeginPlot("##Amp","Frequency [Hz]","Amplitude")) {
                ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.75f);
                ImPlot::PlotShaded("u(f)",result.f.data(),result.ampx.data(),(int)result.f.size());
                ImPlot::PlotLine("u(f)",result.f.data(),result.ampx.data(),(int)result.f.size());
                ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.75f);
                ImPlot::PlotShaded("y(f)",result.f.data(),result.ampy.data(),(int)result.f.size());
                ImPlot::PlotLine("y(f)",result.f.data(),result.ampy.data(),(int)result.f.size());
                ImPlot::EndPlot();
            }
        }
        else if (mode == 1) {
            static std::vector<double> psdx10, psdy10;
            psdx10.resize(result.psdx.size());
            psdy10.resize(result.psdy.size());
            for (int i = 0; i < psdx10.size(); ++i) {
                psdx10[i] = 10*std::log10(result.psdx[i]);
                psdy10[i] = 10*std::log10(result.psdy[i]);
            }
            ImPlot::SetNextPlotLimits(0,500,-100,0);
            if (ImPlot::BeginPlot("##Power","Frequency [Hz]","PSD (dB/Hz)")) {
                ImPlot::PlotLine("x(f)",result.f.data(),psdx10.data(),(int)result.f.size());
                ImPlot::PlotLine("y(f)",result.f.data(),psdy10.data(),(int)result.f.size());
                ImPlot::EndPlot();
            }
        }
        else if (mode == 2) {
            ImPlot::SetNextPlotLimits(1,500,-100,10);
            if (ImPlot::BeginPlot("##Bode1","Frequency [Hz]","Magnitude [dB]",ImVec2(-1,0),0,ImPlotAxisFlags_LogScale)) {
                ImPlot::NextColormapColor();
                ImPlot::PlotLine("##Mag",result.f.data(),result.mag.data(),(int)result.f.size());
                ImPlot::EndPlot();
            }
        }            
        else if (mode == 3) {
            ImPlot::SetNextPlotLimits(1,500,-180,10);
            if (ImPlot::BeginPlot("##Bode2","Frequency [Hz]","Phase Angle [deg]",ImVec2(-1,0),0,ImPlotAxisFlags_LogScale)) {
                ImPlot::NextColormapColor();
                ImPlot::PlotLine("##Phase",result.f.data(),result.phase.data(),(int)result.f.size());
                ImPlot::EndPlot();
            }
        }        

        ImGui::End();
        if (!open)
            quit();
    }
};

int main(int, char**) { 
    FilterToy toy;
    toy.run();
}