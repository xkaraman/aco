#include "AntSystemSimple.h"
#include "utils.h"

#include <iostream>

int main(int argc, char const *argv[]) {


	AntSystemSimple antSystem;

	std::vector<std::vector<double>> coordinates;
	std::vector<std::vector<double>> distances;

	coordinates = readCoords("data/qatar");
	distances = calcDistances(coordinates);

	antSystem.init();
	antSystem.setInputDataMatrix(distances);
	antSystem.setParameters();

	antSystem.run();



//    AntSystem ACO;
//    ACO.readDataFromFile("data/qatar");
//    ACO.setParameters(50,50,0.99,0.35,0.9,0.9,0.1);
//    ACO.runACO();
//
//    double bestLength = ACO.getBestLength();
//    std::vector<int> bestPath = ACO.getBestPath();
//
//    std::cout << "Best Length is: " << bestLength << std::endl;
//    std::cout << "Best Path is: " << bestPath.size() << std::endl;
//    for (int i = 0; i < bestPath.size(); ++i) {
//        std::cout << " " << bestPath[i];
//    }
//    std::cout << std::endl;
   return 0;
}
