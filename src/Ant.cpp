/*
 * Ant.cpp
 *
 *  Created on: Sep 19, 2020
 *      Author: xenofon
 */

#include "Ant.h"

#include <numeric>
#include <algorithm>
#include <random>

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
//	visited.resize(numOfDestinations);
	probability.resize(numOfDestinations);
	std::iota (std::begin(unvisited), std::end(unvisited), 0);
	std::random_shuffle(unvisited.begin(),unvisited.end());

	printVectorInt("\nANT::init", unvisited);

}
void Ant::nextNode(const std::vector<std::vector<double>> &intensity,const std::vector<std::vector<double>> &visibility,const double &a,const double &b) {
	// First node to visit
	if(visited.empty()){
		visited.push_back(unvisited.back());
		unvisited.pop_back();
		currentNode = visited.front();
		printVectorInt("\nANT::nextNode::empty()", unvisited);

	}
	else{
		double sum = 0.0;
		std::vector<double> accumulatedSum;

		// Calculate sum of probability from currentNode to allowedNodes
		for (auto& allowed : unvisited ) {
			sum += std::pow(intensity[currentNode][allowed],a) * std::pow(visibility[currentNode][allowed],b);
			accumulatedSum.push_back(sum);
		}

		// Calcualte proability of each allowed Node
		for (int i = 0 ; i < unvisited.size(); ++i) {
			probability[i] = std::pow(intensity[currentNode][unvisited[i]],a) * std::pow(visibility[currentNode][unvisited[i]],b);;
			probability[i] = probability[i] / sum ;
			accumulatedSum[i] = accumulatedSum[i] / sum;
		}

		printVector("\nANT:nextNode::notempty()",probability);
		printVector("\nANT:nextNode::accu()",accumulatedSum);
		// Choose next node according to its probability
	    std::random_device rd;  //Will be used to obtain a seed for the random number engine
	    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	    std::uniform_real_distribution<double> dis(0.0, 1.0);

	    double exceed = dis(gen);
	    int ind = 0;
	    std::cout << "\nSum to exceed: "<< exceed;
	    std::cout << "\nAccumulated[index]: ";

	    while(accumulatedSum[ind] < exceed){
		    std::cout << "\nAccumulated[" << ind << "]: " << accumulatedSum[ind];
	    	++ind;

//	    	if(ind > accumulatedSum.size()){
//	    		ind = 0;
//	    	}
	    }
	    std::cout << "\nSelecting Node at index: " << ind << " with acccumulated sum " << accumulatedSum[ind];
	}
}
