/*
 * Ant.cpp
 *
 *  Created on: Sep 19, 2020
 *      Author: xenofon
 */

#include "Ant.h"

#include <numeric>
#include <algorithm>

#include "../include/utils.hpp"
#include <cmath>

Ant::Ant():
currentNode(){
	// TODO Auto-generated constructor stub

}

Ant::~Ant() {
	// TODO Auto-generated destructor stub
}

void Ant::init() {
}

void Ant::init(int numOfDestinations) {
	unvisited.resize(numOfDestinations);
	visited.resize(numOfDestinations);
	probability.resize(numOfDestinations);
	std::iota (std::begin(unvisited), std::end(unvisited), 0);
	std::random_shuffle(unvisited.begin(),unvisited.end());

	printVectorInt("ANT::init", unvisited);

}
void Ant::nextNode(const std::vector<std::vector<double>> &intensity,const std::vector<std::vector<double>> &visibility,const double &a,const double &b) {
	// First node to visit
	if(visited.empty()){
		visited.push_back(unvisited.back());
		unvisited.pop_back();
		currentNode = visited.front();
	}
	else{
		double sum = 0.0;

		// Calculate sum of probability from currentNode to allowedNodes
		for (auto& allowed : unvisited ) {
			sum += std::pow(intensity[currentNode][allowed],a) * std::pow(visibility[currentNode][allowed],b);
		}
		// Calcualte proability of each allowed Node
		for (int i = 0 ; i < unvisited.size(); ++i) {
			probability[i] = std::pow(intensity[currentNode][unvisited[i]],a) * std::pow(visibility[currentNode][unvisited[i]],b);;
			probability[i] = probability[i] / sum ;
		}

		// Choose next node according to its probability

	}
}
