/*
 * AntSystemSimple.h
 *
 *  Created on: Sep 19, 2020
 *      Author: xenofon
 */

#ifndef SRC_ANTSYSTEMSIMPLE_H_
#define SRC_ANTSYSTEMSIMPLE_H_

/*
 * Based on http://www.cs.unibo.it/babaoglu/courses/cas05-06/tutorials/Ant_Colony_Optimization.pdf
 */
#include "Ant.h"
#include "AbstractAntSystem.h"

#include <string>
#include <vector>

class AntSystemSimple : public AbstractAntSystem {
public:
	AntSystemSimple();
	virtual ~AntSystemSimple();

	void init();
	void run();

	void setParameters();

	void setInputDataMatrix(const std::vector<std::vector<double>> &inputMatrix);

	double getBestLength() const;
	const std::vector<int>& getBestTour() const;

protected:
	void calculateVisibility();
	void initialize();

private:
	std::vector<std::vector<double>> mDistances;
	std::vector<std::vector<double>> mVisibility;
	std::vector<std::vector<double>> mIntensityOfTrail;
	std::vector<std::vector<double>> mDifIntensityOfTrail;
	std::vector<Ant> mAnts;

	// Parameters for AS
	double mEvaporationOfTrailR;
	double mImportanceOfTrailA;
	double mImportanceOfVisibilityB;
	double mQuantityOfTrailRelatedQ;

	int mNumOfDestinations;
	int mNumOfAnts;
	int mMaxCycles;

//	double mOptimalLength;
//	std::vector<int> mOptimalTour;
};

#endif /* SRC_ANTSYSTEMSIMPLE_H_ */
