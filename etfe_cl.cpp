#include <Mahi/Util.hpp>
#include "ETFE.hpp"

using namespace mahi::util;
using namespace etfe;

/// Column adapter for csv_read_rows
struct CsvData {
    struct Adapter {
        CsvData* p;
        std::size_t r;
        double& operator[](std::size_t c) { return p->data[c][r];  }
        std::size_t size() const  { return p->data.size(); }
    };
    CsvData(std::size_t rows, std::size_t cols) : data(cols, std::vector<double>(rows, 0)) { }
    Adapter operator[](std::size_t r) { return {this, r}; }
    std::size_t size() const { return data[0].size(); }
    std::vector<std::vector<double>> data;
};

int main(int argc, char *argv[])
{
    Options options("etfe.exe","Emperical Transfer Function Estimate");
    options.add_options()
        ("h,help",      "Produces help message")
        ("d,data",      "Filepath to {n x 3} data in CSV format with column headers [Time][Input][Output]", value<std::string>())
        ("n,nsammples", "Number of samples in the input and output (i.e. number of rows excluding header)", value<int>())
        ("f,fs",        "Sampling frequency in Hz",                                                         value<double>());
    auto cl = options.parse(argc, argv);
    if (cl.count("h")) {
        print("{}",options.help());
        return 0;
    }
    else if (cl.count("d")  && cl.count("n") && cl.count("f")) {
        std::string data = cl["d"].as<std::string>();
        int n            = cl["n"].as<int>();
        double fs        = cl["f"].as<double>();
        CsvData csv(n,3);
        if (csv_read_rows(data,csv,1)) {
            ETFE etfe(n,fs);      
            auto& result = etfe.estimate(csv.data[1].data(), csv.data[2].data());
            std::cout << result.size();
            for (int i = 0; i < result.size(); ++i) 
                print("{: .4f} {} {:.4f}i",result.txy[i].real(),result.txy[i].imag() < 0 ? "-" : "+",std::abs(result.txy[i].imag()));  
        }
        else {
            print("Failed to open data file: {}", data);
        }     
    }
    else {
        print("Not enough command line arguments provided");
        print("{}",options.help());
    }
    return 0;
}
