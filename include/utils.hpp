/*
 * utils.h
 *
 *  Created on: Sep 20, 2020
 *      Author: xenofon
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#ifndef INCLUDE_UTILS_H_
#define INCLUDE_UTILS_H_

// Helper functions
inline void
printVector(const std::string& title,
		const std::vector<std::vector<double> >& vector) {
	std::cout << title << "\n";
	for (auto row : vector){
		for (auto col : row) {
			std::cout << col << " ";
		}
		std::cout << "\n";
	}

}

inline void
printVector(const std::string& title,
		const std::vector<double> & vector) {
	std::cout << title << "\n";
	for (auto row : vector){
			std::cout << row << " ";
		}
}

inline void
printVectorInt(const std::string& title,
		const std::vector<int>& vector) {
	std::cout << title << "\n";
	for (auto item : vector){
		std::cout << item << " ";
	}
	std::cout << "\n";
}

inline void
processLine(const std::string& line,
		std::vector<std::vector<double> > &vector) {
	std::stringstream ss(line);
//	std::cout << line.c_str() << " \n";
	int count;
	ss >> count;
	ss >> vector[count-1][0] >> vector[count-1][1];
//	std::cout << count << vector[count-1][0] << vector[count-1][1];
//	std::cout << "Exit\n";
}

inline std::vector< std::vector<double> >
readCoords(const std::string& distancesFilename){
		int line_no = 0;
		int dim = 0;

		std::vector<std::vector<double> > coords;

		std::string line;
		std::ifstream myfile(distancesFilename);

		if (!myfile.is_open()) {
		    std::perror(("Error while opening file " + distancesFilename).c_str());
		    std::exit(-1); // @suppress("Function cannot be resolved")
		}

		while(std::getline(myfile,line)){
			line_no++;
	//		std::cout << line << "\n";
			if(line_no == 5) {
				std::stringstream ss(line);
				std::string temp;
				ss >> temp >> temp >> dim;
	//			std::cout << temp << " " << dim << "\n";
				coords.resize(dim);
				for (int i = 0; i < dim; ++i) {
					coords[i].resize(2);
				}
	//			printVector("Vector: ",coords);
			}
			if (line_no > 7 && !myfile.eof()) {
	//			std::cout << line << "\n";
				processLine(line,coords);
			}
		}
	//	printVector("Vector: ",coords);

	    if (myfile.bad())
	        perror(("error while reading file " + distancesFilename).c_str());

	    myfile.close();

	  return coords;
}
inline std::vector<std::vector<double> >
calcDistances(const std::vector<std::vector<double> >& coords){
	std::vector<std::vector<double> > result(coords.size(),std::vector<double>(coords.size(),0));
	//  printVector("Distances",result);

	//  for (size_t i = 0; i < result.size(); i++) {
	//    result[i].resize(result.size());
	//  }

	for (int i = 0; i < coords.size(); i++) {
	for (int j = i; j < coords.size(); j++) {
	  result[i][j] = std::sqrt( std::pow( coords[i][0] - coords[j][0], 2.0 ) + std::pow(coords[i][1] - coords[j][1], 2.0 ) ); // @suppress("Ambiguous problem")
	  result[j][i] = result[i][j];
	}
	}

	//  printVector("Distances", result);

	return result;
}




#endif /* INCLUDE_UTILS_H_ */
