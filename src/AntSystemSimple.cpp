/*
 * AntSystemSimple.cpp
 *
 *  Created on: Sep 19, 2020
 *      Author: xenofon
 */

#include "AntSystemSimple.h"

#include <iostream>
#include <utils.hpp>
#include <limits>

AntSystemSimple::AntSystemSimple()
:mNumOfAnts(),
 mNumOfDestinations(),
 mMaxCycles(),
 mImportanceOfTrailA(),
 mImportanceOfVisibilityB(),
 mEvaporationOfTrailR(),
 mQuantityOfTrailRelatedQ(),
 mBestLength(std::numeric_limits<double>::max())
{
	// TODO Auto-generated constructor stub
	mAnts.resize(1);
	mDistances.resize(1,std::vector<double>(1));
	mIntensityOfTrail.resize(1,std::vector<double>(1));
	mDifIntensityOfTrail.resize(1,std::vector<double>(1));
	mVisibility.resize(1,std::vector<double>(1));
	mBestTour.resize(1);
}

AntSystemSimple::~AntSystemSimple() {
	// TODO Auto-generated destructor stub
}

void AntSystemSimple::init() {
	this->setParameters();
	mAnts.resize(mNumOfAnts);
	mBestLength = std::numeric_limits<double>::max();

}

void AntSystemSimple::run() {

	this->initialize();
//	std::cout<< "======Init done=======\n";

	int cycleCounter = 0;
	while (cycleCounter < mMaxCycles) {
		for(Ant& ant : mAnts){
			ant.reset(mNumOfDestinations);
			ant.computeTour(mIntensityOfTrail,mVisibility,mImportanceOfTrailA,mImportanceOfVisibilityB);
			ant.computeTourLength(mDistances);
		}
//		std::cout << "======Paths found for each ant=====\n";

		double length;
		for(const Ant& ant:mAnts){
			length = (ant.getTourLength());
			if(length < mBestLength){
				mBestLength = length;
				mBestTour = ant.getTour();
			}
		}
//		std::cout << "======Best tour found =====\n";
		//	printVector("Run::DifIntensity",mDifIntensityOfTrail);

		double antDifIntensityOfTrail;
		for (int i = 0; i < mDifIntensityOfTrail.size(); ++i) {
			for (int j = 0; j < mDifIntensityOfTrail[i].size(); ++j) {
				for(const Ant& ant : mAnts){
		//				printVectorInt("Run::Ant",ant.getTour());
		//				std::cout << "Edge between nodes " << i << " and " << j << " exists: " << ant.edgeExists(i,j) << "\n";
					if (i == j){
						antDifIntensityOfTrail = 0.0;

					}
					else if (ant.edgeExists(i,j)){
						antDifIntensityOfTrail = mQuantityOfTrailRelatedQ / ant.getTourLength();
					}
					else {
						antDifIntensityOfTrail = 0.0;
					}
					mDifIntensityOfTrail[i][j] += antDifIntensityOfTrail;
				}
			}
		}
		//	printVector("Run::DifIntensity",mDifIntensityOfTrail);

		for (int i = 0; i < mIntensityOfTrail.size(); ++i) {
				for (int j = 0; j < mIntensityOfTrail[i].size(); ++j){
					mIntensityOfTrail[i][j] = mEvaporationOfTrailR * mIntensityOfTrail[i][j] + mDifIntensityOfTrail[i][j];
					mDifIntensityOfTrail[i][j] = 0.0;
				}
		}
		cycleCounter++;
	}
}

void AntSystemSimple::setParameters() {
	mNumOfAnts = 60;
	mMaxCycles = 40;
	mImportanceOfTrailA = 2.5;
	mImportanceOfVisibilityB = 1.0;
	mQuantityOfTrailRelatedQ = 50.0;
	mEvaporationOfTrailR = 0.5;
}

void AntSystemSimple::setInputDataMatrix(const std::vector<std::vector<double>> &inputMatrix){
	mDistances = inputMatrix;
	mNumOfDestinations = inputMatrix.size();

	mIntensityOfTrail.resize(mNumOfDestinations);
	for (auto& row : mIntensityOfTrail ) {
		row.resize(mNumOfDestinations);
	}

	mDifIntensityOfTrail.resize(mNumOfDestinations);
	for (auto& row : mDifIntensityOfTrail ) {
		row.resize(mNumOfDestinations);
	}
	mVisibility.resize(mNumOfDestinations);
	for (auto& row : mVisibility ) {
		row.resize(mNumOfDestinations);
	}

}

void AntSystemSimple::calculateVisibility() {
}

void AntSystemSimple::initialize() {
	// Init phase
	for (auto& row : mIntensityOfTrail) {
			for (auto& col: row) {
				col = 0.00001;
		}
	}
//	printVector("Initialize::Intensity",mIntensityOfTrail);

	for (auto& row : mDifIntensityOfTrail) {
			for (auto& col: row) {
				col = 0.0;
		}
	}
//	printVector("Initialize::DifIntensity",mDifIntensityOfTrail);

	for (auto i = 0 ; i < mVisibility.size() ; ++i) {
			for (auto j = 0 ; j < mVisibility.size() ; ++j) {
				mVisibility[i][j] = 1.0 / mDistances[i][j];
		}
	}

	for(Ant& ant : mAnts){
		ant.init(mNumOfDestinations);
//		ant.startNode();
	}
}

double AntSystemSimple::getBestLength() const {
	return mBestLength;
}

const std::vector<int>& AntSystemSimple::getBestTour() const {
	return mBestTour;
}
