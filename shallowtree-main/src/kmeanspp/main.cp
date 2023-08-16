#include "fastKMeans.h"
#include "../utils/data_utils.h"

#include <chrono>

using namespace fastKMeans;
int main()
{


    
    std::shared_ptr<blaze::DynamicMatrix<double>> data;
    data = parse("../input/input_insertAll.txt",60000,2);



    auto start = std::chrono::high_resolution_clock::now();
    fastKMeans::fastKMeansAlgo f (*data);

    auto ret = f.fastKMeansPP(500);
    auto stop = std::chrono::high_resolution_clock::now();

    std::cout << "Runtime: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count() << "\n";
    for (const auto& x : ret)
        std::cout << x << "\n";

}
