/*
 * AbstractAntSystem.h
 *
 *  Created on: Oct 8, 2020
 *      Author: xenofon
 */

#ifndef SRC_ABSTRACTANTSYSTEM_H_
#define SRC_ABSTRACTANTSYSTEM_H_

#include "Ant.h"

#include <vector>

class AbstractAntSystem {
public:
	AbstractAntSystem();
	virtual ~AbstractAntSystem();

	virtual void init() = 0;
	virtual void run() = 0 ;
	virtual void setInputDataMatrix(const std::vector<std::vector<double>> &inputMatrix);

	virtual const double getBestLength() const;
	virtual const std::vector<int> &getBestPath() const;

protected:
	std::vector<Ant> mAnts;

	double mOptimalLength;
	std::vector<int> mOptimalTour;

private:

};

#endif /* SRC_ABSTRACTANTSYSTEM_H_ */
