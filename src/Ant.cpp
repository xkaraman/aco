/*
 * Ant.cpp
 *
 *  Created on: Sep 19, 2020
 *      Author: xenofon
 */

#include "Ant.h"

#include <numeric>
#include <algorithm>

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
	std::iota (std::begin(unvisited), std::end(unvisited), 0);
	std::random_shuffle(unvisited.begin(),unvisited.end());

}
void Ant::nextNode() {
	// First node to visit
	if(visited.empty()){
		visited.push_back(unvisited.front());
	}
}
