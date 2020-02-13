#ifndef ACO_H
#define ACO_H

#include <vector>
#include <string>
#include <random>

class AntSystem {
public:
AntSystem();
// AntSystem(std::vector<std::vector<double> > distances);


std::vector<int> getBestPath();
double getBestLength();
void setParameters(double maxIterations = 3000,
                   double numAnts = 20,
                   double explorationProbability = 0.9,
                   double exhaustRatePheromones = 0.1,
                   double transitionProbabilityA = 1,
                   double transitionProbabilityB = 2,
                   double transitionProbabilityX = 0.1);

void runACOFromFile(const std::string& distancesFilename);
void runACO();
// void runACO(const std::string& distancesFilename);

private:
std::vector<std::vector<double> > mDistances;
std::vector<std::vector<double> > mDestinations;

// ACO Parameters
double mMaxIterations = 3000;
double mNumAnts = 20;
double mExplorationProbability = 0.9;
double mExhaustRatePheromones = 0.1;
double mTransitionProbabilityA = 1;
double mTransitionProbabilityB = 2;
double mTransitionProbabilityX = 0.1;

int mWidth;
int mHeight;
bool mStopped = false;

std::vector<int> mOptimalPath;
double mOptimalLength;
bool mCapacitated = false;

std::vector<int> mDemand;
std::vector<std::vector<int> > mBestList;

std::random_device mRandomDevice;
std::mt19937 mGen;

// Helper functions
std::vector< std::vector<double> > readCoords(const std::string& distancesFilename);
std::vector< std::vector<double> > calcDistances(const std::vector<std::vector<double> > coords);


};

#endif
