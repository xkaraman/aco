/*
 * AbstractAntSystem.h
 *
 *  Created on: Oct 8, 2020
 *      Author: xenofon
 */

#ifndef SRC_ABSTRACTANTSYSTEM_H_
#define SRC_ABSTRACTANTSYSTEM_H_

#include <vector>

class AbstractAntSystem {
public:
	AbstractAntSystem();
	virtual ~AbstractAntSystem();

	virtual void init() = 0;
	virtual void run() = 0 ;

protected:
	double mOptimalLength;
	std::vector<int> mOptimalTour;

private:

};

#endif /* SRC_ABSTRACTANTSYSTEM_H_ */
