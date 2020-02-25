#include "../include/ACO.h"

#include <iostream>

int main(int argc, char const *argv[]) {
    AntSystem ACO;
    ACO.runACOFromFile("data/djibouti");
    double bestLength = ACO.getBestLength();
    std::vector<int> bestPath = ACO.getBestPath();

    std::cout << "Best Length is: " << bestLength << std::endl;
    std::cout << "Best Path is: " << bestPath.size() << std::endl;
    for (int i = 0; i < bestPath.size(); ++i) {
        std::cout << " " << bestPath[i];
    }
    std::cout << std::endl;
    return 0;
}
