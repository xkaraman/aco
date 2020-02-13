#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "../include/ACO.h"

AntSystem::AntSystem():
        mWidth()
        ,mHeight()
        ,mOptimalLength()
        ,mGen(mRandomDevice())
{

}

/*********************************************************************
* Comment
*********************************************************************/
// AntSystem::AntSystem(std::vector<std::vector<double> > distances):
// mWidth()
// ,mHeight()
// ,mOptimalLength()
// ,mGen(mRandomDevice())
// {
//     mDistances = distances;
// }
/*********************************************************************
* Comment
*********************************************************************/

std::vector<int> AntSystem::getBestPath()
{
    return mOptimalPath;
}

/*********************************************************************
* Comment
*********************************************************************/
double AntSystem::getBestLength()
{
    return mOptimalLength;
}
/*********************************************************************
* Comment
*********************************************************************/
void AntSystem::runACOFromFile(const std::string& distancesFilename)
{
    mDestinations = readCoords(distancesFilename);
    mDistances  = calcDistances(mDestinations);

    runACO();
}
/*********************************************************************
* Comment
*********************************************************************/
void AntSystem::runACO()
{
	double bestLength;
	int iteration;
	double sumP;
	int nextMove;
	double length;

	if( mDistances.empty() ){
		std::cout << "Exiting" << ' \n';
		return;
	}

	int numberDestinations = mDistances.size();
	std::vector<std::vector<double>> pheromoneDeposistedT(numberDestinations,std::vector<double>(numberDestinations,0));
	std::vector<std::vector<double>> desirabilityTransitionH(numberDestinations,std::vector<double>(numberDestinations,0));

	double t0 = 0;
	double nearNb;

	std::vector<int> bestTour(numberDestinations + 1);

	for (int i = 0; i < pheromoneDeposistedT.size(); ++i) {
		for (int j = 0; j < pheromoneDeposistedT[i].size(); ++j) {

		}
	}



}

/*********************************************************************
* Comment
*********************************************************************/
std::vector< std::vector<double> > AntSystem::readCoords(const std::string& distancesFilename){
  std::string line;
  std::ifstream myfile(distancesFilename);

  int line_no = 0;
  int dim;

  // std::vector<int> customers;
  std::vector<std::vector<double> > coords;

  if (myfile.is_open()){
    while ( getline (myfile,line) ){
      line_no++;
      //diavazoume tin grammi 5 apo to txt arxeio pou exei ta dedomena  https://stackoverflow.com/questions/26288145/how-to-read-a-specific-line-from-file-using-fstream-c
      if (line_no == 5) {
        // sLine contains the fifth line in the file.
        // std::cout << line << '\n';
        std::stringstream ss;
        std::string temp;
        ss << line;
        while ( !ss.eof() ) {
          ss >> temp;
          if (std::isdigit(temp[0]) ) {
            dim = std::stoi(temp);
            // std::cout << dim;
            // customers.resize(dim);
            coords.resize(dim);
            for (size_t i = 0; i < coords.size(); i++) {
              coords[i].resize(2);
            }
          }
        }
      }

      if ( (line_no >= 7) && (line_no < 7 + dim) ) {
        // std::cout << line << '\n';
        int a;
        myfile >> a;
        // std::cout << a << std::endl;
        myfile >> coords[a-1][0] >> coords[a-1][1];
        // std::cout << coords[a-1][0] << "  " << coords[a-1][1] << std::endl;
      }
    }
    myfile.close();
  }
  else {
    std::cout << "Unable to open file";
  }

  return coords;
}
/*********************************************************************
* Comment
*********************************************************************/
std::vector<std::vector<double> > AntSystem::calcDistances(const std::vector<std::vector<double> > coords){
  std::vector<std::vector<double> > result(coords.size());
  for (size_t i = 0; i < result.size(); i++) {
    result[i].resize(result.size());
  }

  for (size_t i = 0; i < coords.size(); i++) {
    for (size_t j = i+1; j < coords.size(); j++) {
      result[i][j] = std::sqrt( std::pow( coords[i][0] - coords[j][0], 2.0 ) + std::pow(coords[i][1] - coords[j][1], 2.0 ) ); // @suppress("Ambiguous problem")
      result[j][i] = result[i][j];
    }
  }

  // for (size_t i = 0; i < result.size(); i++) {
  //   for (size_t j = 0; j < result.size(); j++) {
  //     if (i==j) {
  //       result[i][j] = 0;
  //     }
  //     result[j][i] = result[i][j];
  //     std::cout<<result[i][j] << " ";
  //   }
  //   std::cout << std::endl;
  // }
  return result;
}
