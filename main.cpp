#include <iostream>
#include "min_matching.h"
#include "rural_postman.h"

auto StartTime() {
    return std::chrono::high_resolution_clock::now();
}

void EndTime(std::chrono::high_resolution_clock::time_point start) {
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
}

void RunMockTestOne() {
    int n = 4;
    int m = 6;
    std::vector<int> u = {1, 1, 1, 2, 2, 3}; // starting node
    std::vector<int> v = {2, 3, 4, 3, 4, 4}; // ending node
    std::vector<int> w = {3, 1, 1, 1, 1, 1}; // weight
    std::vector<int> chosen_edges = {0, 5}; // indexes of chosen edges
    auto ans = RuralPostman().Run(n, m, u,
                                  v, w, chosen_edges);
    std::cout << "Cost: " << ans.first << "\n";
    std::cout << "Cycle traversal: ";
    for (int i = 0; i < ans.second.size(); ++i) {
        std::cout << ans.second[i] << " ";
    }
    std::cout << "\n";
}
