/*
 * AntSystemSimple.cpp
 *
 *  Created on: Sep 19, 2020
 *      Author: xenofon
 */

#include "AntSystemSimple.h"

AntSystemSimple::AntSystemSimple()
:mNumOfAnts(),
 mNumOfDestinations(),
 mMaxCycles(),
 mImportanceOfTrailA(),
 mImportanceOfVisibilityB(),
 mEvaporationOfTrailR(),
 mQuantityOfTrailRelated()
{
	// TODO Auto-generated constructor stub
	mAnts.resize(1);
	mDistances.resize(1,std::vector<double>(1));
	mIntensityOfTrail.resize(1,std::vector<double>(1));

}

AntSystemSimple::~AntSystemSimple() {
	// TODO Auto-generated destructor stub
}

void AntSystemSimple::init() {
	this->setParameters();
	mAnts.resize(mNumOfAnts);
}

void AntSystemSimple::run() {

	this->initialize();

	for(Ant& ant : mAnts){
		ant.nextNode();
	}
}

void AntSystemSimple::setParameters() {
	mNumOfAnts = 20;
	mMaxCycles = 100;
	mImportanceOfTrailA = 1.0;
	mImportanceOfVisibilityB = 5.0;
	mQuantityOfTrailRelated = 100.0;
	mEvaporationOfTrailR = 0.5;
}

void AntSystemSimple::setInputDataMatrix(const std::vector<std::vector<double>> &inputMatrix){
	mDistances = inputMatrix;
	mNumOfDestinations = inputMatrix.size();
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

	for (auto& row : mDifIntensityOfTrial) {
			for (auto& col: row) {
				col = 0.0;
		}
	}

	for (auto i = 0 ; i < mVisibility.size() ; ++i) {
			for (auto j = 0 ; j < mVisibility.size() ; ++j) {
				mVisibility[i][j] = 1.0 / mDistances[i][j];
		}
	}

	for(Ant& ant : mAnts){
		ant.init(mNumOfDestinations);
		ant.nextNode();
	}

}
