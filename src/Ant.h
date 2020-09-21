/*
 * Ant.h
 *
 *  Created on: Sep 19, 2020
 *      Author: xenofon
 */

#ifndef SRC_ANT_H_
#define SRC_ANT_H_

#include <vector>


class Ant {
public:
	Ant();
	virtual ~Ant();

	void init();
	void init(int numOfDestinations);

	void pickNode();

	void nextNode();
	float calculateTransisionProbability();

private:
	std::vector<int> unvisited;
	std::vector<int> visited;
	int currentNode;
};

#endif /* SRC_ANT_H_ */
