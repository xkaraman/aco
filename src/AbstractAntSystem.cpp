/*
 * AbstractAntSystem.cpp
 *
 *  Created on: Oct 8, 2020
 *      Author: xenofon
 */

#include "AbstractAntSystem.h"

#include <limits>

AbstractAntSystem::AbstractAntSystem():
mOptimalLength(std::numeric_limits<double>::max()){
	// TODO Auto-generated constructor stub
	mOptimalTour.resize(1);
}

AbstractAntSystem::~AbstractAntSystem() {
	// TODO Auto-generated destructor stub
}

