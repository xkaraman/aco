#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "../include/ACO.h"

AntSystem::AntSystem():
        mWidth()
        ,mHeight()
        ,mOptimalLength()
        ,mGen(mRandomDevice())
{

}

/*********************************************************************
* Comment
*********************************************************************/
// AntSystem::AntSystem(std::vector<std::vector<double> > distances):
// mWidth()
// ,mHeight()
// ,mOptimalLength()
// ,mGen(mRandomDevice())
// {
//     mDistances = distances;
// }
/*********************************************************************
* Comment
*********************************************************************/

std::vector<int> AntSystem::getBestPath()
{
    return mOptimalPath;
}

/*********************************************************************
* Comment
*********************************************************************/
double AntSystem::getBestLength()
{
    return mOptimalLength;
}
/*********************************************************************
* Comment
*********************************************************************/
void AntSystem::runACOFromFile(const std::string& distancesFilename)
{
    mDestinations = readCoords(distancesFilename);
    mDistances  = calcDistances(mDestinations);

    runACO();
}
/*********************************************************************
* Comment
*********************************************************************/
void AntSystem::runACO()
{
	double sumP;
	int nextMove;
	double length;

	if( mDistances.empty() ){
		std::cout << "Exiting... Distance Matrix is Empty" << "\n";
		return;
	}

	int numberDestinations = mDistances.size();
	std::vector<std::vector<double>> pheromoneDeposistedT(numberDestinations,std::vector<double>(numberDestinations,0));
	std::vector<std::vector<double>> desirabilityTransitionH(numberDestinations,std::vector<double>(numberDestinations,0));

	double t0 = 0;
	std::vector<int> bestTour(numberDestinations + 1);

//	printVector("Coordinates:",mDestinations);
//	printVector("Matrix Table:",mDistances);
	for (int i = 0; i < pheromoneDeposistedT.size(); ++i) {
		for (int j = 0; j < pheromoneDeposistedT[i].size(); ++j) {
			if ( i == j) {
				desirabilityTransitionH[i][j] = 1 / mDistances[i][j];
			}
		}
	}

	int nextNode = 0;
	double bestLength = 0.0;

	std::vector<double> results;
	results.resize(mMaxIterations);

	for (int i = 0; i < numberDestinations - 1; ++i) {
		bestLength = bestLength + mDistances[i][i+1];
	}

	std::uniform_int_distribution<> dis(0, numberDestinations - 1);
	double min = std::numeric_limits<double>::max();
	int startingNode = dis(mGen);
	std::vector<int> neighborUnvisited;
	double nearNeighbor;

	bestTour[0] = startingNode;

	for (int i = 0; i < numberDestinations; ++i) {
		neighborUnvisited.push_back(i);
	}

    // REMOVE ELEMENT WITH VALUE 'Startingnode' FROM VECTOR NBUnvisited
	neighborUnvisited.erase(std::remove(neighborUnvisited.begin(),neighborUnvisited.end(),startingNode),neighborUnvisited.end());

	for (int i = 0; i < neighborUnvisited.size(); ++i) {
		if ( min > mDistances[startingNode][neighborUnvisited[i] ]) {
			min = mDistances[startingNode][neighborUnvisited[i] ];
			nextNode = neighborUnvisited[i];
		}
	}

	nearNeighbor = nearNeighbor + mDistances[startingNode][nextNode];
	neighborUnvisited.erase(std::remove(neighborUnvisited.begin(),neighborUnvisited.end(), nextNode), neighborUnvisited.end());
	bestTour[1] = nextNode;
	startingNode = nextNode;

	int count = 1;
//	bool listEmpty = false;

	while (!neighborUnvisited.empty()) {
		count++;
		min = std::numeric_limits<double>::max();
		for (int i = 0; i < neighborUnvisited.size(); ++i) {
			if (min > mDistances[startingNode][neighborUnvisited[i]]) {
				min = mDistances[startingNode][neighborUnvisited[i]];
				nextNode = neighborUnvisited[i];
			}
		}
		nearNeighbor = nearNeighbor + mDistances[startingNode][nextNode];
	    neighborUnvisited.erase(std::remove(neighborUnvisited.begin(),neighborUnvisited.end(),nextNode),neighborUnvisited.end());
		bestTour[count] = nextNode;
		startingNode = nextNode;
//		std::cout << neighborUnvisited.size() << std::endl;
	}
//	std::cout << "Out of While" << std::endl;

	bestTour[count + 1] = bestTour[0];
	nearNeighbor = nearNeighbor + mDistances[bestTour[count]][bestTour[count + 1]];

	for (int i = 0; i < numberDestinations; i++) {
		for (int j = 0; j < numberDestinations; j++) {
			t0 = (1 / ((nearNeighbor * numberDestinations)));
			pheromoneDeposistedT[i][j] = t0;
		}
	}

	int iteration;
	double tMax = 0;
	tMax = (1 / ((1 - mExhaustRatePheromones))) * (1 / nearNeighbor);
	double tMin	= 0;
	tMin = tMax * (1 - std::pow(0.05, 1 / numberDestinations)) / ((numberDestinations / 2 - 1) * std::pow(0.05, 1 / numberDestinations));
	double totalRandomLength[500];
	double randomLength = 0;

	// TODO: Add 500 as parameter
	// Represents how many random solutions to be used to initialise TO
	int randomSolutions = 500;
	for (int g = 0; g < randomSolutions; ++g) {
		randomLength = 0;
		std::vector<int> randomUnvisited;
		int start = dis(mGen);
		std::vector<int> randomTour(numberDestinations+1);
		randomTour[0] =  start;

		std::iota(randomUnvisited.begin(), randomUnvisited.end(), 0);
	    randomUnvisited.erase(std::remove(randomUnvisited.begin(),randomUnvisited.end(),nextNode),randomUnvisited.end());

	    int countRandom = 1;

	    while (!randomUnvisited.empty()) {
			std::uniform_int_distribution<int> dist(0, randomUnvisited.size() - 1);
			int next = dist(mGen);
			randomTour[countRandom] = randomUnvisited[next];
		    randomUnvisited.erase(std::remove(randomUnvisited.begin(),randomUnvisited.end(),randomUnvisited[next]),randomUnvisited.end());
		    randomLength = randomLength + mDistances[randomTour[countRandom-1]][randomTour[countRandom]];
		    countRandom++;
	    }

	    randomTour[countRandom] = randomTour[0];
	    randomLength = randomLength + mDistances[randomTour[countRandom-1]][randomTour[countRandom]];

	    totalRandomLength[g] = randomLength;
	}

	double DC = 0;
	double SDC = 0;
	for (int i = 0; i < randomSolutions; ++i) {
		DC = DC + std::abs(totalRandomLength[i] - totalRandomLength[i + 1]);
		SDC = SDC + std::pow((totalRandomLength[i] - totalRandomLength[i + 1]) - DC, 2);
	}

	DC = DC / double(randomSolutions);
	SDC = std::sqrt(float(SDC / randomSolutions));

	double temperature = (DC + 3 * SDC) / (std::log(1.0f / 0.1f));
	double activeLength = nearNeighbor;
	std::vector<int> activeSolution(numberDestinations + 1);
	activeSolution = bestTour;

	// ACO ITERATIONS
	while (iteration < mMaxIterations) {
		if(mStopped){
			return;
		}

		std::vector<int> tourIteration(numberDestinations + 1);
		double lengthIteration = std::pow(nearNeighbor,10);

	}






}

/*********************************************************************
* Comment
*********************************************************************/
std::vector< std::vector<double> > AntSystem::readCoords(const std::string& distancesFilename){
	int line_no = 0;
	int dim = 0;

	std::vector<std::vector<double> > coords;

	std::string line;
	std::ifstream myfile(distancesFilename);

	if (!myfile.is_open()) {
	    std::perror(("Error while opening file " + distancesFilename).c_str());
	    std::exit(-1); // @suppress("Function cannot be resolved")
	}

	while(std::getline(myfile,line)){
		line_no++;
//		std::cout << line << "\n";
		if(line_no == 5) {
			std::stringstream ss(line);
			std::string temp;
			ss >> temp >> temp >> dim;
//			std::cout << temp << " " << dim << "\n";
			coords.resize(dim);
			for (int i = 0; i < dim; ++i) {
				coords[i].resize(2);
			}
//			printVector("Vector: ",coords);
		}
		if (line_no > 7 && !myfile.eof()) {
//			std::cout << line << "\n";
			processLine(line,coords);
		}
	}

    if (myfile.bad())
        perror(("error while reading file " + distancesFilename).c_str());

    myfile.close();
//	printVector("Vector: ",coords);

  return coords;
}

void AntSystem::setParameters(double maxIterations, double numAnts,
		double explorationProbability, double exhaustRatePheromones,
		double transitionProbabilityA, double transitionProbabilityB,
		double transitionProbabilityX) {
}

/*********************************************************************
* Comment
*********************************************************************/
std::vector<std::vector<double> > AntSystem::calcDistances(const std::vector<std::vector<double> > coords){
  std::vector<std::vector<double> > result(coords.size());
  for (size_t i = 0; i < result.size(); i++) {
    result[i].resize(result.size());
  }

  for (size_t i = 0; i < coords.size(); i++) {
    for (size_t j = i+1; j < coords.size(); j++) {
      result[i][j] = std::sqrt( std::pow( coords[i][0] - coords[j][0], 2.0 ) + std::pow(coords[i][1] - coords[j][1], 2.0 ) ); // @suppress("Ambiguous problem")
      result[j][i] = result[i][j];
    }
  }


  return result;
}

void AntSystem::printVector(const std::string title,
		const std::vector<std::vector<double> > vector) const {
	std::cout << title << std::endl;
	for (int i = 0; i < vector.size(); i++) {
	std::cout<< i <<":\t";
	for (int j = 0; j < vector[i].size(); j++) {
	std::cout<< vector[i][j] << "\t";
	}
	std::cout << std::endl;
   }
}

void AntSystem::processLine(const std::string line,
		std::vector<std::vector<double> > &vector) const {
	std::stringstream ss(line);
//	std::cout << line.c_str() << " \n";
	int count;
	ss >> count;
	ss >> vector[count-1][0] >> vector[count-1][1];
//	std::cout << count << vector[count-1][0] << vector[count-1][1];
//	std::cout << "Exit\n";
}
