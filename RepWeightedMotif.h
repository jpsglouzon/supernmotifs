/*
 * RepWeightedMotif.h
 *
 *  Created on: Jul 28, 2014
 *      Author: Séhi
 */

#ifndef REPWEIGHTEDMOTIF_H_
#define REPWEIGHTEDMOTIF_H_

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <math.h>       /* fabs */


#include "Eigen/Dense"

using namespace std;

class RepWeightedMotif {

    Eigen::ArrayXXf matALLMotifWeigthed;
    float minOcc;
    map<string,float> allStructureFeatureForWeigthOfMotifs;
    vector< map<string,float> > nmotifsForEachStructure;


public:

	RepWeightedMotif();
    RepWeightedMotif(map<string,float>, vector<map<string, float> > , float);
	virtual ~RepWeightedMotif();

    Eigen::ArrayXXf filterMotifs(float);

	//Getters
    const map<string, float>& getAllStructureFeatureForWeigthOfMotifs() const {
		return allStructureFeatureForWeigthOfMotifs;
	}

    const Eigen::ArrayXXf& getMatAllMotifWeigthed() const {
		return matALLMotifWeigthed;
	}
};


#endif /* REPWEIGHTEDMOTIF_H_ */
