/*
 * RepWeightedMotif.cpp
 *
 *  Created on: Jul 28, 2014
 *      Author: Séhi
 */

#include "RepWeightedMotif.h"

/*********************************************************************************/

RepWeightedMotif::RepWeightedMotif() {
	// TODO Auto-generated constructor stub
}

/*********************************************************************************/

RepWeightedMotif::RepWeightedMotif(map<string,float> allStructureFeature, vector< map<string,float> > nmotifsForEachStructure, float minOcc)
{
	// TODO Auto-generated constructor stub
this->allStructureFeatureForWeigthOfMotifs=allStructureFeature;
this->nmotifsForEachStructure=nmotifsForEachStructure;
this->minOcc=minOcc;
matALLMotifWeigthed=filterMotifs(this->minOcc);
}


/*********************************************************************************/

RepWeightedMotif::~RepWeightedMotif() {
	// TODO Auto-generated destructor stub
}

/*******************************************************************************************************/
Eigen::ArrayXXf RepWeightedMotif::filterMotifs(float minOcc){


    map<string,float> allStructureFeatureForWeigthOfMotifs_filtered;
    std::map<string,float>::iterator it = allStructureFeatureForWeigthOfMotifs.begin();
    int sumAllNmotifs=0;
    int colMatALLMotifFiltered=0;

    for (std::map<string, float>::iterator it2 = allStructureFeatureForWeigthOfMotifs.begin(); it2 != allStructureFeatureForWeigthOfMotifs.end(); ++it2)
    {sumAllNmotifs=sumAllNmotifs+it2->second;}

    if (minOcc==0) //Automatic minimum occurrence
    {minOcc=float(sumAllNmotifs)/float(allStructureFeatureForWeigthOfMotifs.size());}

    for (std::map<string, float>::iterator it2 = allStructureFeatureForWeigthOfMotifs.begin(); it2 != allStructureFeatureForWeigthOfMotifs.end(); ++it2)
    {
        if(it2->second>=minOcc)
        {
          allStructureFeatureForWeigthOfMotifs_filtered.insert ( std::pair<string,int>(it2->first,it2->second) );
        }
    }
    allStructureFeatureForWeigthOfMotifs=allStructureFeatureForWeigthOfMotifs_filtered;

    unsigned int nrows=nmotifsForEachStructure.size();
    Eigen::ArrayXXf matALLMotifWeigthed(nrows,allStructureFeatureForWeigthOfMotifs_filtered.size());
    for (unsigned int i = 0; i < nrows; i++) {
        colMatALLMotifFiltered=0;
        for (std::map<string, float>::iterator it2 = allStructureFeatureForWeigthOfMotifs_filtered.begin();it2 != allStructureFeatureForWeigthOfMotifs_filtered.end(); ++it2)
        {
            it = nmotifsForEachStructure[i].find(it2->first);
            if(it!=nmotifsForEachStructure[i].end())
            {matALLMotifWeigthed(i, colMatALLMotifFiltered)=log2(it->second+1);}
            else{matALLMotifWeigthed(i, colMatALLMotifFiltered)=0;}
            colMatALLMotifFiltered++;
        }
    }
   return matALLMotifWeigthed;
}
