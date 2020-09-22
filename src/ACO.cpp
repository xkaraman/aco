#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

#include "../include/twooptswap.hpp"
#include "../include/ACO.h"

AntSystem::AntSystem():
        mDistances()
		,mDestinations()
		,mWidth()
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
void AntSystem::readDataFromFile(const std::string& distancesFilename)
{
//    printVector("Destinations", mDestinations);
	mDestinations = readCoords(distancesFilename);
//    printVector("Destinations", mDestinations);

//    printVector("Distances", mDistances);
    mDistances  = calcDistances(mDestinations);
//    printVector("Distances", mDistances);

    setParameters(150,50,0.98,0.3);
//    std::cout << mMaxIterations << mNumAnts << std::endl;
//    runACO();
}
/*********************************************************************
* Comment
*********************************************************************/
void AntSystem::runACO()
{
	double sumP = 0.0;
	int nextMove = 0;
	double length = 0.0;

	if( mDistances.empty() ){
		std::cout << "Exiting... Distance Matrix is Empty" << std::endl;
		return;
	}

	int numberDestinations = mDestinations.size();

	// Initialize with zeros
	std::vector<std::vector<double>> pheromoneDeposistedT(numberDestinations,std::vector<double>(numberDestinations,0));
	std::vector<std::vector<double>> desirabilityTransitionH(numberDestinations,std::vector<double>(numberDestinations,0));

	double t0 = 0.0;
	std::vector<int> bestTour(numberDestinations + 1);

	for (int i = 0; i < desirabilityTransitionH.size(); ++i) {
		for (int j = 0; j < desirabilityTransitionH[i].size(); ++j) {
			if (i == j) {
				desirabilityTransitionH[i][j] = 1 / std::numeric_limits<double>::max();
			}
			desirabilityTransitionH[i][j] = 1 / mDistances[i][j];
		}
	}

	int nextNode = 0;
	double bestLength = 0.0;

	std::vector<double> results(mMaxIterations);

	// Find best Length by neighbors
	for (int i = 0; i < numberDestinations - 1; ++i) {
		bestLength = bestLength + mDistances[i][i+1];
	}

	std::uniform_int_distribution<> dis(0, numberDestinations - 1);
	double min = std::numeric_limits<double>::max();
	int startingNode = dis(mGen);
	double nearNeighbor = 0.0;
	int count = 0;
	bestTour[count] = startingNode;

	std::vector<int> neighborUnvisited(numberDestinations);
	std::iota(neighborUnvisited.begin(), neighborUnvisited.end(), 0);

    // REMOVE ELEMENT WITH VALUE 'startingNode' FROM VECTOR neighborUnvisited
	neighborUnvisited.erase(std::remove(neighborUnvisited.begin(),neighborUnvisited.end(),startingNode),neighborUnvisited.end());


	// Find min between startingNode and all other neighbours that are unvisited
//	for (int i = 0; i < neighborUnvisited.size(); ++i) {
//		if ( min > mDistances[startingNode][neighborUnvisited[i] ]) {
//			min = mDistances[startingNode][neighborUnvisited[i] ];
//			nextNode = neighborUnvisited[i];
//		}
//	}
//
//	nearNeighbor = nearNeighbor + mDistances[startingNode][nextNode];
//	neighborUnvisited.erase(std::remove(neighborUnvisited.begin(),neighborUnvisited.end(), nextNode), neighborUnvisited.end());
//	bestTour[1] = nextNode;
//	startingNode = nextNode;

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
//	std::cout << "Out of While Count= " << count << std::endl;
//	std::cout << "Best Tour size= " << bestTour.size() << std::endl;

	// End where we started
	bestTour[count + 1] = bestTour[0];
	nearNeighbor = nearNeighbor + mDistances[bestTour[count]][bestTour[count + 1]];

	t0 = (1 / ((nearNeighbor * mNumAnts)));
	for (int i = 0; i < pheromoneDeposistedT.size(); i++) {
		for (int j = 0; j < pheromoneDeposistedT.size(); j++) {
			pheromoneDeposistedT[i][j] = t0;
		}
	}

	int iteration = 1;
	double tMax = 0;
	tMax = (1 / ((1 - mExhaustRatePheromones))) * (1 / nearNeighbor);
	double tMin	= 0;
	tMin = tMax * (1 - std::pow(0.05, 1 / numberDestinations)) / ((numberDestinations / 2 - 1) * std::pow(0.05, 1 / numberDestinations));
	double randomLength = 0;

	// TODO: Add 500 as parameter
	// Represents how many random solutions to be used to initialise temperature for simulated annealing
	int randomSolutions = 1000;
	std::vector<double> totalRandomLength(randomSolutions);
	for (int g = 0; g < randomSolutions; ++g) {
		randomLength = 0;
		std::vector<int> randomTour(numberDestinations+1);

		std::uniform_int_distribution<> dis(0, numberDestinations - 1);
		int start = dis(mGen);

		std::vector<int> randomUnvisited(numberDestinations);
		std::iota(randomUnvisited.begin(), randomUnvisited.end(), 0);
	    randomUnvisited.erase(std::remove(randomUnvisited.begin(),randomUnvisited.end(),start),randomUnvisited.end());

		randomTour[0] =  start;

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

	double DC = 0.0;
	double SDC = 0.0;
	for (int i = 0; i < randomSolutions - 1; ++i) {
//		std::cout << "TotalRandomLength= "<< totalRandomLength[i] << std::endl;
//		std::cout << "TotalRandomLength= "<< totalRandomLength[i+1] << std::endl;

		DC = DC + std::abs(totalRandomLength[i] - totalRandomLength[i + 1]);
		SDC = SDC + std::pow((totalRandomLength[i] - totalRandomLength[i + 1]) - DC, 2.0);
	}

	DC = DC / randomSolutions;
//	std::cout << "DC is= " << DC << std::endl;
	SDC = std::sqrt(SDC / randomSolutions); // @suppress("Ambiguous problem")
//	std::cout << "SDC is=" << SDC << std::endl;

	double temperature = (DC + 3 * SDC) / (std::log(1.0f / 0.1f));
	double activeLength = nearNeighbor;
	std::vector<int> activeSolution(numberDestinations + 1);
	activeSolution = bestTour;

	std::cout << "Starting ACO Iterations..." << std::endl;

//	printVector("Coordinates:",mDestinations);
//	printVector("Matrix Table:",mDistances);
//	printVector("Pheromone", pheromoneDeposistedT);
//	printVector("Desirabilty", desirabilityTransitionH);

	// ACO ITERATIONS
	while (iteration < mMaxIterations) {
//		std::cout << "Iteration " << iteration << std::endl;
		if(mStopped){
			return;
		}

		std::vector<int> tourIteration(numberDestinations + 1);
		double lengthIteration = std::pow(nearNeighbor,10.0);

		int moves = 0;
		for (int ant = 0; ant < mNumAnts; ++ant) {
			moves = 0;
			std::vector<int> tour(numberDestinations + 1);
			std::uniform_int_distribution<int> dist(0,numberDestinations - 1);

			tour[moves] = dist(mGen);

			std::vector<int> unvisited(numberDestinations);
			std::iota(unvisited.begin(),unvisited.end(),0);
	        unvisited.erase(std::remove(unvisited.begin(),unvisited.end(),tour[0]),unvisited.end());

			std::vector<double> choice;
	        for (int trip = 0; trip < numberDestinations - 1; ++trip) {
				int c = tour[trip];
				choice.clear();

//				std::cout << "-----Entering for unvisited---- Size: " << unvisited.size() << std::endl;
				for (int i = 0; i < unvisited.size(); ++i) {
					int j = unvisited[i];
//					std::cout << "J:= " << j << " C:= " << c << std::endl;
//					std::cout << "PheromeT: " << pheromoneDeposistedT[c][j];
//					std::cout <<" DesirabilityH: " << desirabilityTransitionH[c].size() <<std::endl;
					double value = std::pow(pheromoneDeposistedT[c][j], mTransitionProbabilityA ) * std::pow(desirabilityTransitionH[c][j], mTransitionProbabilityB);
					choice.push_back(value);
				}
//				std::cout << "-----Exiting for unvisited---- Size: " << unvisited.size() << std::endl;

				std::uniform_real_distribution<double> unif(0,1);
				double random1 = unif(mGen);
//				std::cout << "Random double= " << random1 << std::endl;

				if (random1 < mExplorationProbability) {
					int maxIndex = std::max_element(choice.begin(),choice.end()) - choice.begin();

//					std::cout << "-----Choice List ----" << std::endl;
//					for (int index = 0; index < choice.size(); ++index) {
//						std::cout <<  choice[index]  << " ";
//					}
//					std::cout << "Max Index= " << maxIndex << std::endl;
//					double maxValue = *std::max_element(choice.begin(), choice.end());
					nextMove = unvisited[maxIndex];
				} else {
					std::vector<double> p;
					sumP = 0;
					for (int i = 0; i < unvisited.size(); ++i) {
						int j = unvisited[i];
						double value = (std::pow( pheromoneDeposistedT[c][j], mTransitionProbabilityA) * std::pow(desirabilityTransitionH[c][j], mTransitionProbabilityB));
						sumP = sumP + value;
						p.push_back(value);
					}

					double cumulativeSum = 0;
					double randomNum = unif(mGen);

					for (int i = 0; i < p.size(); ++i) {
						p[i] = p[i] / sumP;
						p[i] = cumulativeSum + p[i];
						cumulativeSum = p[i];
					}

					for (int j = 0; j < p.size() - 1; ++j) {
						if (randomNum >= p[j] && randomNum < p[j+1]) {
							nextMove = unvisited[j];
							break;
						}

					}

				}

				if (nextMove == c) {
					nextMove = unvisited[0];
				}

				tour[trip + 1] = nextMove;
		        unvisited.erase(std::remove(unvisited.begin(),unvisited.end(),tour[trip + 1]),unvisited.end());

		        pheromoneDeposistedT[c][tour[trip+1]] = std::max(pheromoneDeposistedT[c][tour[trip + 1]] * (1 - mTransitionProbabilityX) + mTransitionProbabilityX * t0, tMin);
			}

	        tour[tour.size()-1] = tour[0];
	        length = 0;
	        for (int i = 0; i < tour.size() -1 ; ++i) {
	        	length =+ mDistances[tour[i]][tour[i+1]];
			}

	        if (length < lengthIteration) {
				tourIteration = tour;
				lengthIteration = length;
			}
		}

		int improve = 0;
		while (improve <= randomSolutions) {
//			std::cout << "Improve = " << improve << std::endl;
			double newDistance = 0.0;
			for (int i = 0; i < tourIteration.size() - 1; ++i) {
				newDistance = newDistance + mDistances[tourIteration[i]][tourIteration[i+1]];
			}

			for (int i = 1; i < tourIteration.size() - 2; ++i) {
				for (int j = i + 1; j < tourIteration.size() - 1; ++j) {
					std::vector<int> newRoute;
					newRoute = TwoOptSwap::optSwap(tourIteration,i,j);
					double newLength = 0.0;
					for (int l = 0; l < tourIteration.size() - 1; ++l) {
						newLength += mDistances[newRoute[l]][newRoute[l+1]];
					}

					if ( newLength < newDistance) {
						tourIteration = newRoute;
						newDistance = newLength;
						improve = 0;
						break;
					} else {
						improve++;
					}
				}

			}
		}

		lengthIteration = 0.0;
		for (int i = 0; i < tourIteration.size() - 1; ++i) {
			lengthIteration = lengthIteration + mDistances[tourIteration[i]][tourIteration[i+1]];
		}

		if (lengthIteration < bestLength) {
			bestLength = lengthIteration;
			bestTour = tourIteration;
		}

		if (activeLength > lengthIteration) {
			activeSolution = tourIteration;
			activeLength = lengthIteration;
		} else {
			double C = lengthIteration - activeLength;
			std::uniform_real_distribution<double> unif(0,1);
			if (unif(mGen) < std::exp(-(float)C / (float)temperature)) {
				activeSolution = tourIteration;
				activeLength = lengthIteration;
			}
		}
//
		temperature = temperature * 0.9999;

		for (size_t i = 0; i < pheromoneDeposistedT.size(); i++) {
			for (size_t j = 0; j < pheromoneDeposistedT[i].size(); j++) {
				pheromoneDeposistedT[i][j] =  std::max(pheromoneDeposistedT[i][j] * (1 - pheromoneDeposistedT[i][j] / (tMin + tMax)), tMin);
			}
		}

		tMax = (1 / ((1 - mExhaustRatePheromones))) * (1 / bestLength);
	    tMin = tMax * (1 - std::pow(0.05, 1 / numberDestinations)) / ((numberDestinations / 2 - 1) * std::pow(0.05, 1 / numberDestinations));

		for (int i = 0; i < activeSolution.size()-1; ++i) {
			pheromoneDeposistedT[activeSolution[i]][activeSolution[i+1]] = std::min(pheromoneDeposistedT[activeSolution[i]][activeSolution[i + 1]] + (pheromoneDeposistedT[activeSolution[i]][activeSolution[i + 1]] / (tMax + tMin)) * (1 / activeLength), tMax);
		}

		results[iteration] = bestLength;
		iteration++;
	}

//	std::cout << "Best Length= " << bestLength << std::endl;
	mOptimalPath = bestTour;
	mOptimalLength = bestLength;

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
//	printVector("Vector: ",coords);

    if (myfile.bad())
        perror(("error while reading file " + distancesFilename).c_str());

    myfile.close();

  return coords;
}

void AntSystem::setParameters(double maxIterations, double numAnts,
		double explorationProbability, double exhaustRatePheromones,
		double transitionProbabilityA, double transitionProbabilityB,
		double transitionProbabilityX) {

	mMaxIterations = maxIterations;
	mNumAnts = numAnts;
	mExplorationProbability = explorationProbability;
	mExhaustRatePheromones = exhaustRatePheromones;
	mTransitionProbabilityA = transitionProbabilityA;
	mTransitionProbabilityB = transitionProbabilityB;
	mTransitionProbabilityX = transitionProbabilityX;
}

/*********************************************************************
* Comment
*********************************************************************/
std::vector<std::vector<double> > AntSystem::calcDistances(const std::vector<std::vector<double> >& coords){
  std::vector<std::vector<double> > result(coords.size(),std::vector<double>(coords.size(),0));
//  printVector("Distances",result);

//  for (size_t i = 0; i < result.size(); i++) {
//    result[i].resize(result.size());
//  }

  for (int i = 0; i < coords.size(); i++) {
    for (int j = i; j < coords.size(); j++) {
      result[i][j] = std::sqrt( std::pow( coords[i][0] - coords[j][0], 2.0 ) + std::pow(coords[i][1] - coords[j][1], 2.0 ) ); // @suppress("Ambiguous problem")
      result[j][i] = result[i][j];
    }
  }

//  printVector("Distances", result);

  return result;
}

void AntSystem::printVector(const std::string& title,
		const std::vector<std::vector<double> >& vector) const {
	std::cout << title << std::endl;
	for (int i = 0; i < vector.size(); i++) {
		std::cout<< i <<":\t";
		for (int j = 0; j < vector[i].size(); j++) {
			std::cout<< vector[i][j] << '\t';
		}
		std::cout << std::endl;
	}
}

void AntSystem::processLine(const std::string& line,
		std::vector<std::vector<double> > &vector) const {
	std::stringstream ss(line);
//	std::cout << line.c_str() << " \n";
	int count;
	ss >> count;
	ss >> vector[count-1][0] >> vector[count-1][1];
//	std::cout << count << vector[count-1][0] << vector[count-1][1];
//	std::cout << "Exit\n";
}
