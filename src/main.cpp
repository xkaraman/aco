#include "AntSystemSimple.h"
#include <iostream>
#include "../include/utils.hpp"

clock_t START_TIMER;

clock_t tic()
{
    return START_TIMER = clock();
}

void toc(clock_t start = START_TIMER)
{
    std::cout
        << "Elapsed time: "
        << (clock() - start) / (double)CLOCKS_PER_SEC << "s"
        << std::endl;
}

int main(int argc, char const *argv[]) {

	std::vector<std::vector<double>> coordinates;
	std::vector<std::vector<double>> distances;

	coordinates = readCoords("data/djibouti");
	distances = calcDistances(coordinates);

//	printVector("MAIN::coordinates", coordinates);
//	printVector("MAIN::distances", distances);

	AntSystemSimple antSystem;
	double sum = 0;
	int repeat = 20;

	for (int i = 0; i < repeat; ++i) {

		antSystem.init();
		antSystem.setInputDataMatrix(distances);
		antSystem.setParameters();

		tic();
		antSystem.run();
		toc();

		std::cout << "Best Length is: " << antSystem.getBestLength() << std::endl;
		printVectorInt("Best Path is: ", antSystem.getBestTour());
		sum += antSystem.getBestLength();
	}
	std::cout << "Average Best Length is: " << sum/repeat<< std::endl;



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
