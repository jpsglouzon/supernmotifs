/*
 * RepSupermotif.h
 *
 *  Created on: Jul 28, 2014
 *      Author: Séhi
 */

#ifndef REPSUPERMOTIF_H_
#define REPSUPERMOTIF_H_

#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include "Eigen/Dense"

using namespace std;

class RepSupernmotif {
private:

    Eigen::ArrayXXf matS_supernmotifs;
    Eigen::ArrayXXf matS_supernmotifsFull;
    Eigen::ArrayXXf matNM_supernmotifs;
    Eigen::ArrayXXf matNM_supernmotifsFull;
    Eigen::ArrayXXf matDissimS2S;
    Eigen::ArrayXXf matDissimS2NM;
    Eigen::RowVectorXf singularVal;
    Eigen::RowVectorXf singularValFull;


public:

	RepSupernmotif();

    RepSupernmotif(Eigen::ArrayXXf, int, int);

	virtual ~RepSupernmotif();

    Eigen::ArrayXXf repSupermotifComputeMatrixDissimCosS2S(Eigen::ArrayXXf);
    Eigen::ArrayXXf repSupermotifComputeMatrixDissimCosS2NM(Eigen::ArrayXXf,Eigen::ArrayXXf);

	//Getters
    const Eigen::ArrayXXf& getMatDissimS2Nm() const {
		return matDissimS2NM;
	}

    const Eigen::ArrayXXf& getMatDissimS2S() const {
		return matDissimS2S;
	}

    const Eigen::ArrayXXf& getMatNmSupernmotifs() const {
        return matNM_supernmotifs;
    }

    const Eigen::ArrayXXf& getMatSSupernmotifs() const {
		return matS_supernmotifs;
	}

    const Eigen::RowVectorXf& getSingularVal() const {
		return singularVal;
	}

    const Eigen::ArrayXXf& getMatSSupernmotifsFull() const{
        return matS_supernmotifsFull;
    }

    const Eigen::RowVectorXf& getSingularValFull() const{
        return singularValFull;
    }

    const Eigen::ArrayXXf& getMatNM_supernmotifsFull() const
    {
        return matNM_supernmotifsFull;
    }
};

#endif /* REPSUPERMOTIF_H_ */
