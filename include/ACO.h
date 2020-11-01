#ifndef ACO_H
#define ACO_H

#include "../src/AbstractAntSystem.h"

#include <vector>
#include <string>
#include <random>

class AntColonyOptimizationSystem : public AbstractAntSystem {
public:
	AntColonyOptimizationSystem();
	virtual ~AntColonyOptimizationSystem();

	void setInputDataMatrix(const std::vector<std::vector<double>> &inputMatrix);

	// This can be called as
	// .setParameters() for default Parameter settings
	void setParameters(double maxIterations = 100,
					   double numAnts = 40,
					   double explorationProbability = 0.9,
					   double exhaustRatePheromones = 0.2,
					   double transitionProbabilityA = 2.5,
					   double transitionProbabilityB = 3.5,
					   double transitionProbabilityX = 0.1);

	void init();
	void run();

//	std::vector<int> getBestPath();
//	double getBestLength();



//	void readDataFromFile(const std::string& distancesFilename);
	void runACO();
	// void runACO(const std::string& distancesFilename);


protected:
	void initialize();
	void calculateVisibility();

private:
	std::vector< std::vector<double> > mDistances;
	std::vector< std::vector<double> > mDestinations;
	std::vector< std::vector<double> > mVisibility;
	std::vector< std::vector<double> > mIntensityOfTrail;
	std::vector< std::vector<double> > mDifIntensityOfTrail;

	double	mNumAnts = 150;
	double	mMaxCycles = 50;
	int		mNumOfDestinations;

	// ACO Parameters
	double mExplorationProbability = 0.7;
	double mExhaustRatePheromones = 0.5;
	double mTransitionProbabilityA = 2.5;
	double mTransitionProbabilityB = 3.5;
	double mTransitionProbabilityX = 0.1;


	int mWidth = 0;
	int mHeight = 0;
	bool mStopped = false;

//	double mOptimalLength = 0.0;
	bool mCapacitated = false;

	//std::vector<int> mDemand;
	//std::vector<std::vector<int> > mBestList;

	std::random_device mRandomDevice;
	std::mt19937 mGen;

	// Helper functions
	//std::vector< std::vector<double> > readCoords(const std::string& distancesFilename);
	//std::vector< std::vector<double> > calcDistances(const std::vector<std::vector<double> >& coords);
	//
	//void printVector(const std::string& title,const std::vector<std::vector<double>>& vector) const;
	//void processLine(const std::string& line, std::vector<std::vector<double>>& vector) const;
};

#endif
