/*
 * AbstractAntSystem.cpp
 *
 *  Created on: Oct 8, 2020
 *      Author: xenofon
 */

#include "AbstractAntSystem.h"

#include <limits>
#include <iostream>
#include <cstdlib>

AbstractAntSystem::AbstractAntSystem():
mOptimalLength(std::numeric_limits<double>::max()){
	// TODO Auto-generated constructor stub
	mOptimalTour.resize(1);
}

AbstractAntSystem::~AbstractAntSystem() {
	// TODO Auto-generated destructor stub
}

const double AbstractAntSystem::getBestLength() const {
	return mOptimalLength;
}

const std::vector<int> &AbstractAntSystem::getBestPath() const {
	return mOptimalTour;
}

void AbstractAntSystem::setInputDataMatrix(const std::vector<std::vector<double>> &inputMatrix){
	std::cerr << "Implement this functions in subclass";
	exit(EXIT_FAILURE);
};
