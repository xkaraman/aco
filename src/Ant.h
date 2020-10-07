/*
 * Ant.h
 *
 *  Created on: Sep 19, 2020
 *      Author: xenofon
 */

#ifndef SRC_ANT_H_
#define SRC_ANT_H_

#include <vector>
#include <random>

class Ant {
public:
	Ant();
	Ant(const Ant&) = default;
	virtual ~Ant();

	void init();
	void init(int numOfDestinations);

	void startNode();
	void computeTour(const std::vector<std::vector<double>> &intensity,const std::vector<std::vector<double>> &visibility,const double &a,const double &b);
	double computeTourLength(const std::vector<std::vector<double>> &distances);
	bool edgeExists(const int i,const int j) const;
	void reset(const int numOfDestinations);

	std::vector<int>	getTour() const;
	double 				getTourLength() const;

	float calculateTransisionProbability();

	static std::random_device randomDevice;

protected:

private:
	std::vector<int> 	unvisited;
	std::vector<int> 	visited;
	std::vector<double> probability;

	int currentNode;
	double tourLength;

	std::mt19937_64 gen;
};
#endif /* SRC_ANT_H_ */
