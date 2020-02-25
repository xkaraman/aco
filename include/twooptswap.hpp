/*
 * twooptswap.hpp
 *
 *  Created on: 25 Feb 2020
 *      Author: xenofon
 */

#ifndef INCLUDE_TWOOPTSWAP_HPP_
#define INCLUDE_TWOOPTSWAP_HPP_

#include <vector>

namespace TwoOptSwap
{
	std::vector<int> optSwap(const std::vector<int>& tempTour,const int& i,const int& v){
		int size = tempTour.size();
		std::vector<int> newTour(size);

		for (int j = 0; j <= i - 1; ++j) {
			newTour[j] = tempTour[j];
		}

		int eff = 0;

		for (int j = i; j <= v; ++j) {
			newTour[j] = tempTour[v - eff];
			eff++;
		}

		for (int j = v + 1; j < size; ++j) {
			newTour[j] = tempTour[j];
		}

		return newTour;
	};

}



#endif /* INCLUDE_TWOOPTSWAP_HPP_ */
