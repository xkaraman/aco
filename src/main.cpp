#include "../include/ACO.h"

#include <iostream>

int main(int argc, char const *argv[]) {
    AntSystem ACO;

    ACO.runACOFromFile("../../data/qatar");
    auto bestLength = ACO.getBestLength();
    auto bestPath = ACO.getBestPath();

    std::cout << "Best Length is: " << bestLength << '\n';
    std::cout << "Best Path is: " << '\n';
    for (size_t i = 0; i < bestPath.size(); i++) {
        std::cout << bestPath[i] << ' ';
    }

    return 0;
}
