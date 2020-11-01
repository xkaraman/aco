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
#include <chrono>

Ant::Ant():
currentNode(),
tourLength(),
gen(randomDevice())
{
	// TODO Auto-generated constructor stub
}

Ant::~Ant() {
	// TODO Auto-generated destructor stub
}

void Ant::init() {
}

void Ant::init(int numOfDestinations) {
	unvisited.resize(numOfDestinations);
	visited.reserve(numOfDestinations);
	probability.reserve(numOfDestinations);
	std::iota (std::begin(unvisited), std::end(unvisited), 0);
	std::random_shuffle(unvisited.begin(),unvisited.end());
//	printVectorInt("\nANT::init", unvisited);
}

void Ant::reset(const int numOfDestinations){
	unvisited.clear();
	unvisited.resize(numOfDestinations);
	visited.clear();
	visited.reserve(numOfDestinations);
	probability.clear();
	probability.reserve(numOfDestinations);
	std::iota (std::begin(unvisited), std::end(unvisited), 0);
	std::random_shuffle(unvisited.begin(),unvisited.end());
	startNode();
}

void Ant::startNode(){
	if(visited.empty()){
		visited.push_back(unvisited.back());
		unvisited.pop_back();
		currentNode = visited.front();
//		printVectorInt("\nANT::startNode::empty()", unvisited);

	}
}

// ACO version
void Ant::computeTour(const std::vector<std::vector<double>> &intensity,const std::vector<std::vector<double>> &visibility,const double &a,const double &b,const double &exploration) {
//	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(randomDevice()); //Standard mersenne_twister_engine seeded with rd()
	//		// Choose next node according to its probability
	std::uniform_real_distribution<double> dis(0.0, 1.0);
	while (!unvisited.empty()) {
		double sum = 0.0;
		std::vector<double> accumulatedSum;
		probability.clear();

		// Calculate sum of probability from currentNode to allowedNodes
		for (auto& allowed : unvisited ) {
			sum += std::pow(intensity[currentNode][allowed],a) * std::pow(visibility[currentNode][allowed],b);
			accumulatedSum.push_back(sum);
		}

		for (auto& allowed : unvisited ) {
			probability.push_back((std::pow(intensity[currentNode][allowed],a) * std::pow(visibility[currentNode][allowed],b)) / sum);
		}

		for (auto& acc : accumulatedSum ) {
			acc = acc / sum;
		}
		// Calcualte proability of each allowed Node
		//		for (int i = 0 ; i < unvisited.size(); ++i) {
		//			probability[i] = std::pow(intensity[currentNode][unvisited[i]],a) * std::pow(visibility[currentNode][unvisited[i]],b);;
		//			probability[i] = probability[i] / sum ;
		//			accumulatedSum[i] = accumulatedSum[i] / sum;
		//		}

		//			printVector("\nANT:nextNode::notempty()",probability);
		//			printVector("\nANT:nextNode::accu()",accumulatedSum);

		// For ACO version
		double random1 = dis(gen);
		int nextMoveIndex;
		if (random1 < exploration) {
			nextMoveIndex = std::max_element(probability.begin(),probability.end()) - probability.begin();
		} else {
		double exceed = dis(gen);
		int ind = 0;
//			    std::cout << "\nSum to exceed: "<< exceed;
		//	    std::cout << "\nAccumulated[index]: ";

		while(accumulatedSum[ind] < exceed){
			//		    std::cout << "\nAccumulated[" << ind << "]: " << accumulatedSum[ind];
			++ind;

		}
			//	    	if(ind > accumulatedSum.size()){
			//	    		ind = 0;
			//	    	}
			nextMoveIndex = ind;
		}
//		std::cout << "\nSelecting Node at index: " << ind << " with acccumulated sum " << accumulatedSum[ind];
//
//		printVectorInt("\nANT::computeTour::notempty()::unvisitedBeforeChose", unvisited);
//		printVectorInt("\nANT::computeTour::notempty()::visitedBeforeChose", visited);

		currentNode = unvisited[nextMoveIndex];
		visited.push_back(unvisited[nextMoveIndex]);
		// REMOVE ELEMENT WITH VALUE 'ind' FROM VECTOR unvisited
		unvisited.erase(std::remove(unvisited.begin(),unvisited.end(),unvisited[nextMoveIndex]),unvisited.end());

//		printVectorInt("\nANT::computeTour::notempty()::unvisitedAfterChose", unvisited);
//		printVectorInt("\nANT::computeTour::notempty()::visitedAfterchose", visited);
	}
}

void Ant::computeTour(const std::vector<std::vector<double>> &intensity,const std::vector<std::vector<double>> &visibility,const double &a,const double &b) {
//	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(randomDevice()); //Standard mersenne_twister_engine seeded with rd()
	//		// Choose next node according to its probability
	std::uniform_real_distribution<double> dis(0.0, 1.0);
	while (!unvisited.empty()) {
		double sum = 0.0;
		std::vector<double> accumulatedSum;

		// Calculate sum of probability from currentNode to allowedNodes
		for (auto& allowed : unvisited ) {
			sum += std::pow(intensity[currentNode][allowed],a) * std::pow(visibility[currentNode][allowed],b);
			accumulatedSum.push_back(sum);
		}

		for (auto& allowed : unvisited ) {
			probability.push_back((std::pow(intensity[currentNode][allowed],a) * std::pow(visibility[currentNode][allowed],b)) / sum);
		}

		for (auto& acc : accumulatedSum ) {
			acc = acc / sum;
		}
		// Calcualte proability of each allowed Node
		//		for (int i = 0 ; i < unvisited.size(); ++i) {
		//			probability[i] = std::pow(intensity[currentNode][unvisited[i]],a) * std::pow(visibility[currentNode][unvisited[i]],b);;
		//			probability[i] = probability[i] / sum ;
		//			accumulatedSum[i] = accumulatedSum[i] / sum;
		//		}

		//			printVector("\nANT:nextNode::notempty()",probability);
		//			printVector("\nANT:nextNode::accu()",accumulatedSum);

		double exceed = dis(gen);
		int ind = 0;
//			    std::cout << "\nSum to exceed: "<< exceed;
		//	    std::cout << "\nAccumulated[index]: ";

		while(accumulatedSum[ind] < exceed){
			//		    std::cout << "\nAccumulated[" << ind << "]: " << accumulatedSum[ind];
			++ind;

			//	    	if(ind > accumulatedSum.size()){
			//	    		ind = 0;
			//	    	}
		}
//		std::cout << "\nSelecting Node at index: " << ind << " with acccumulated sum " << accumulatedSum[ind];
//
//		printVectorInt("\nANT::computeTour::notempty()::unvisitedBeforeChose", unvisited);
//		printVectorInt("\nANT::computeTour::notempty()::visitedBeforeChose", visited);

		currentNode = unvisited[ind];
		visited.push_back(unvisited[ind]);
		// REMOVE ELEMENT WITH VALUE 'ind' FROM VECTOR unvisited
		unvisited.erase(std::remove(unvisited.begin(),unvisited.end(),unvisited[ind]),unvisited.end());

//		printVectorInt("\nANT::computeTour::notempty()::unvisitedAfterChose", unvisited);
//		printVectorInt("\nANT::computeTour::notempty()::visitedAfterchose", visited);
	}
}

double Ant::computeTourLength(
		const std::vector<std::vector<double> > &distances) {
	double length = 0;

//	printVectorInt("\nANT::computeTourLength", visited);
	for (auto it = visited.begin(); it != visited.end(); ){
	    int current = *it;
	    int next = *it;

	    if( it == visited.end()-1){
		    next = *visited.begin();
		    ++it;
	    }
	    else {
	    	next = (*++it);
	    }
//		std::cout << "\nDistance at Nodes " << current <<" and " <<  next << " is: " << distances[current][next];
		length += distances[current][next];
	}

	tourLength = length;
	return length;
}

std::random_device Ant::randomDevice;

bool Ant::edgeExists(const int i,const int j) const{
	auto it = std::find(visited.begin(),visited.end(),i);
	return (*++it) == j;
}
std::vector<int> Ant::getTour() const {
	return visited;
}

double Ant::getTourLength() const {
	return tourLength;
}
