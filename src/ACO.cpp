#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

#include "../include/twooptswap.hpp"
#include "../include/ACO.h"
#include "../include/utils.hpp"

AntColonyOptimizationSystem::AntColonyOptimizationSystem():
        mDistances()
		,mDestinations()
		,mWidth()
        ,mHeight()
//        ,mOptimalLength()
        ,mGen(mRandomDevice())
		,mMaxCycles()
		,mNumOfDestinations()
{
	this->setParameters();
}

AntColonyOptimizationSystem::~AntColonyOptimizationSystem() {
	// TODO Auto-generated destructor stub
}

void AntColonyOptimizationSystem::setInputDataMatrix(
		const std::vector<std::vector<double> > &inputMatrix) {

	if (inputMatrix.empty()) {
		std::cerr << "Input Matrix is Empty";
		exit(EXIT_FAILURE);
	}

	mDistances = inputMatrix;
	mNumOfDestinations = inputMatrix.size();

	mIntensityOfTrail.resize(mNumOfDestinations);
	for (auto& row : mIntensityOfTrail ) {
		row.resize(mNumOfDestinations);
	}

	mDifIntensityOfTrail.resize(mNumOfDestinations);
	for (auto& row : mDifIntensityOfTrail ) {
		row.resize(mNumOfDestinations);
	}
	mVisibility.resize(mNumOfDestinations);
	for (auto& row : mVisibility ) {
		row.resize(mNumOfDestinations);
	}
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

//std::vector<int> AntColonyOptimizationSystem::getBestPath()
//{
//    return mOptimalPath;
//}

/*********************************************************************
* Comment
*********************************************************************/
//double AntColonyOptimizationSystem::getBestLength()
//{
//    return mOptimalLength;
//}

/*********************************************************************
* Comment
*********************************************************************/
void AntColonyOptimizationSystem::runACO()
{
	double sumP = 0.0;
	int nextMove = 0;
	double length = 0.0;

	if( mDistances.empty() ){
		std::cout << "Exiting... Distance Matrix is Empty" << std::endl;
		return;
	}

	int numberDestinations = mNumOfDestinations;

	// Initialize with zeros
	std::vector<std::vector<double>> pheromoneDeposistedT(numberDestinations,std::vector<double>(numberDestinations,0));
//	std::vector<std::vector<double>> desirabilityTransitionH(numberDestinations,std::vector<double>(numberDestinations,0));

	double t0 = 0.0;
	std::vector<int> bestTour(numberDestinations + 1);

//	for (int i = 0; i < desirabilityTransitionH.size(); ++i) {
//		for (int j = 0; j < desirabilityTransitionH[i].size(); ++j) {
//			if (i == j) {
//				desirabilityTransitionH[i][j] = 1 / std::numeric_limits<double>::max();
//			}
//			desirabilityTransitionH[i][j] = 1 / mDistances[i][j];
//		}
//	}

//	printVector("desirability: ",desirabilityTransitionH);

	int nextNode = 0;
	double bestLength = 0.0;

	std::vector<double> results(mMaxCycles);

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
	int randomSolutions = 10;
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

//	printVector("Matrix Table:",mDistances);
//	printVector("Pheromone", pheromoneDeposistedT);
//	printVector("Desirabilty", desirabilityTransitionH);

	// ACO ITERATIONS

//	int cycleCounter = 0;
//	while ( cycleCounter < mMaxCycles) {
//		for (Ant& ant : mAnts) {
//			ant.reset(mNumOfDestinations);
//			ant.computeTour(mIntensityOfTrail,mVisibility,mTransitionProbabilityA,mTransitionProbabilityB,mExplorationProbability);
//			ant.computeTourLength(mDistances);
//		}
//	}

	while (iteration < mMaxCycles) {
//		std::cout << "Iteration " << iteration << std::endl;
		if(mStopped){
			return;
		}

		std::vector<int> tourIteration;
		double lengthIteration = std::pow(nearNeighbor,10.0);

		for (Ant& ant : mAnts) {
			ant.reset(mNumOfDestinations);
			ant.computeTour(mIntensityOfTrail,mVisibility,mTransitionProbabilityA,mTransitionProbabilityB,mExplorationProbability);
			ant.computeTourLength(mDistances);
			if (length < lengthIteration) {
				tourIteration = ant.getTour();
				lengthIteration = ant.getTourLength();
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
		tourIteration.push_back(tourIteration[0]);
//		printVectorInt("Tour Iterations: ",tourIteration);
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
	mOptimalTour = bestTour;
	mOptimalLength = bestLength;

}

void AntColonyOptimizationSystem::init() {
	mAnts.resize(mNumAnts);
	mOptimalLength = std::numeric_limits<double>::max();

	this->initialize();

}

void AntColonyOptimizationSystem::initialize() {
	// Init phase
	for (auto& row : mIntensityOfTrail) {
			for (auto& col: row) {
				col = 0.00001;
		}
	}
//	printVector("Initialize::Intensity",mIntensityOfTrail);

	for (auto& row : mDifIntensityOfTrail) {
			for (auto& col: row) {
				col = 0.0;
		}
	}
//	printVector("Initialize::DifIntensity",mDifIntensityOfTrail);

//	for (auto i = 0 ; i < mVisibility.size() ; ++i) {
//			for (auto j = 0 ; j < mVisibility.size() ; ++j) {
//				mVisibility[i][j] = 1.0 / mDistances[i][j];
//		}
//	}

	calculateVisibility();
//	printVector("Initialize::Visibility",mVisibility);

	for(Ant& ant : mAnts){
		ant.init(mNumOfDestinations);
//		ant.startNode();
	}
}

void AntColonyOptimizationSystem::calculateVisibility() {
	for (auto i = 0 ; i < mVisibility.size() ; ++i) {
			for (auto j = 0 ; j < mVisibility.size() ; ++j) {
				mVisibility[i][j] = 1.0 / mDistances[i][j];
		}
	}
}



void AntColonyOptimizationSystem::run() {
	runACO();
}

void AntColonyOptimizationSystem::setParameters(double maxIterations, double numAnts,
		double explorationProbability, double exhaustRatePheromones,
		double transitionProbabilityA, double transitionProbabilityB,
		double transitionProbabilityX) {

	mNumAnts = numAnts;
	mMaxCycles = maxIterations;
	mExplorationProbability = explorationProbability;
	mExhaustRatePheromones = exhaustRatePheromones;
	mTransitionProbabilityA = transitionProbabilityA;
	mTransitionProbabilityB = transitionProbabilityB;
	mTransitionProbabilityX = transitionProbabilityX;
}
