/*
 * RepSupermotif.cpp
 *
 *  Created on: Jul 28, 2014
 *      Author: Séhi
 */

#include "RepSupermotif.h"


RepSupernmotif::RepSupernmotif() {
    // TODO Auto-generated constructor stub
}

RepSupernmotif::RepSupernmotif(Eigen::ArrayXXf matALLMotifWeigthed,int nbSupNmotifs,int p_outputOption) {
	// TODO Auto-generated constructor stub

        if((p_outputOption==0)||(p_outputOption==1))//Compute the super n-motifs representation
        {
			int maxDim;
			if(matALLMotifWeigthed.rows()>matALLMotifWeigthed.cols())
			{maxDim=matALLMotifWeigthed.cols();}
			else{maxDim=matALLMotifWeigthed.rows();}

            if (nbSupNmotifs>maxDim)
			{
                cerr<<"warning: the number of super n-motifs selected (-k) is greater than the maximum number of computed super n-motifs. -k is set to the maximum number of computed super n-motifs."<<endl;
                nbSupNmotifs=maxDim;
			}


			//compute SVD//
            Eigen::JacobiSVD<Eigen::MatrixXf> svd(matALLMotifWeigthed.matrix(), Eigen::ComputeThinU);
            singularVal=svd.singularValues();

            //Brocken stick model

            if(nbSupNmotifs==0)
            {
                Eigen::RowVectorXf singularVal_precent=singularVal/singularVal.sum();
                Eigen::RowVectorXf brModel=singularVal;
                float tempbrModel;

                for(int i=0; i<brModel.size(); i++)
                {
                    tempbrModel=0;
                    for(int j=i; j<brModel.size(); j++)
                    {
                        tempbrModel+=1/float(j+1);
                    }
                    brModel(i)=tempbrModel/brModel.size();
                }

                int nbSupNmotifs2=0;
                while(singularVal_precent(nbSupNmotifs2)>brModel(nbSupNmotifs2))
                {nbSupNmotifs2++;}

                nbSupNmotifs=nbSupNmotifs2;

                if (nbSupNmotifs<2)
                {
                    nbSupNmotifs=2;
                }

            }

            matS_supernmotifs=svd.matrixU().leftCols(nbSupNmotifs);
            for(int i=0;i<nbSupNmotifs;i++)
			{matS_supernmotifs.col(i)=matS_supernmotifs.col(i)*(singularVal(i));}

            //Compute the SS dissimilarities using the cosine dissimilarity
            //if (p_outputOption==0)
            //{matDissimS2S=repSupermotifComputeMatrixDissimCosS2S(matS_supernmotifs);}

		}
        else if(p_outputOption==2)
        {

            int maxDim;
            if(matALLMotifWeigthed.rows()>matALLMotifWeigthed.cols())
            {maxDim=matALLMotifWeigthed.cols();}
            else{maxDim=matALLMotifWeigthed.rows();}

            if (nbSupNmotifs>maxDim)
            {
                cerr<<"warning: the number of super n-motifs selected (-k) is greater than the maximum number of computed super n-motifs. -k is set to the maximum number of computed super n-motifs."<<endl;
                nbSupNmotifs=maxDim;
            }

            //*******compute SVD*************************//
            Eigen::JacobiSVD<Eigen::MatrixXf> svd(matALLMotifWeigthed.matrix(), Eigen::ComputeThinU | Eigen::ComputeThinV);
            singularValFull=svd.singularValues();

            //Brocken stick model

            if(nbSupNmotifs==0)
            {
                Eigen::RowVectorXf singularVal_precent=singularValFull/singularValFull.sum();
                Eigen::RowVectorXf brModel=singularValFull;
                float tempbrModel;

                for(int i=0; i<brModel.size(); i++)
                {
                    tempbrModel=0;
                    for(int j=i; j<brModel.size(); j++)
                    {
                        tempbrModel+=1/float(j+1);
                    }
                    brModel(i)=tempbrModel/brModel.size();
                }

                int nbSupNmotifs2=0;
                while(singularVal_precent(nbSupNmotifs2)>brModel(nbSupNmotifs2))
                {nbSupNmotifs2++;}

                nbSupNmotifs=nbSupNmotifs2;

                if (nbSupNmotifs<2)
                {
                    nbSupNmotifs=2;
                }

            }

            singularVal=singularValFull.head(nbSupNmotifs);
            matS_supernmotifs=svd.matrixU().leftCols(nbSupNmotifs);
            matNM_supernmotifs=svd.matrixV().leftCols(nbSupNmotifs);

            for(int i=0;i<nbSupNmotifs;i++)
            {
             matS_supernmotifs.col(i)=matS_supernmotifs.col(i)*(singularVal(i));
             matNM_supernmotifs.col(i)=matNM_supernmotifs.col(i)*(singularVal(i));
            }

            //Compute the SS dissimilarities using the cosine dissimilarity
            //matDissimS2S=repSupermotifComputeMatrixDissimCosS2S(matS_supernmotifs);
            matDissimS2NM=repSupermotifComputeMatrixDissimCosS2NM(matS_supernmotifs,matNM_supernmotifs);

        }
        else if(p_outputOption==3)//Compute the secondary structure dissimilarities based on the n-motifs representation
        {
         cout<<"No super n-motifs computed (-p 3) ..."<<endl;
         //matDissimS2S=repSupermotifComputeMatrixDissimCosS2S(matALLMotifWeigthed);
        }
        else if ((p_outputOption==4))
        {
        cout<<"No super n-motifs computed (-p 4) ..."<<endl;
        }
        else if (p_outputOption==5)
		{
			int maxDim;
			if(matALLMotifWeigthed.rows()>matALLMotifWeigthed.cols())
			{maxDim=matALLMotifWeigthed.cols();}
			else{maxDim=matALLMotifWeigthed.rows();}

			//compute SVD
            if (nbSupNmotifs>maxDim)
			{
                cerr<<"warning: the number of super n-motifs selected (-k) is greater than the maximum number of computed super n-motifs. -k is set to the maximum number of computed super n-motifs."<<endl;
                nbSupNmotifs=maxDim;
			}

			//compute SVD//
            Eigen::JacobiSVD<Eigen::MatrixXf> svd(matALLMotifWeigthed.matrix(), Eigen::ComputeThinU | Eigen::ComputeThinV);
            singularValFull=svd.singularValues();
            //Brocken stick model
            if(nbSupNmotifs==0)
            {
                Eigen::RowVectorXf singularVal_precent=singularValFull/singularValFull.sum();
                Eigen::RowVectorXf brModel=singularValFull;
                float tempbrModel;

                for(int i=0; i<brModel.size(); i++)
                {
                    tempbrModel=0;
                    for(int j=i; j<brModel.size(); j++)
                    {
                        tempbrModel+=1/float(j+1);
                    }
                    brModel(i)=tempbrModel/brModel.size();
                }

                int nbSupNmotifs2=0;
                while(singularVal_precent(nbSupNmotifs2)>brModel(nbSupNmotifs2))
                {nbSupNmotifs2++;}

                nbSupNmotifs=nbSupNmotifs2;

                if (nbSupNmotifs<2)
                {
                    nbSupNmotifs=2;
                }

            }


            matS_supernmotifsFull=svd.matrixU();
            matNM_supernmotifsFull=svd.matrixV();

            singularVal=singularValFull.head(nbSupNmotifs);
            matS_supernmotifs=matS_supernmotifsFull.leftCols(nbSupNmotifs);
            matNM_supernmotifs=matNM_supernmotifsFull.leftCols(nbSupNmotifs);;


            //Compute super n-motifs representation of n-motifs

            for(int i=0;i<nbSupNmotifs;i++)
			{
				matS_supernmotifs.col(i)=matS_supernmotifs.col(i)*(singularVal(i));
                matNM_supernmotifs.col(i)=matNM_supernmotifs.col(i)*(singularVal(i));
			}

            //Compute the dissimilarities between SS and n-motifs in the super n-motif space using the cosine dissimilarity
            matDissimS2NM=repSupermotifComputeMatrixDissimCosS2NM(matS_supernmotifs,matNM_supernmotifs);
            //Compute the SS dissimilarities using the cosine dissimilarity
            //matDissimS2S=repSupermotifComputeMatrixDissimCosS2S(matS_supernmotifs);
		}


}

RepSupernmotif::~RepSupernmotif() {
	// TODO Auto-generated destructor stub
}
/*
Eigen::ArrayXXf RepSupernmotif::repSupermotifComputeMatrixDissimCosS2S(Eigen::ArrayXXf matSuperMotif){

    Eigen::MatrixXf matDistS2S = Eigen::MatrixXf::Zero(matSuperMotif.rows(),matSuperMotif.rows());
	int n=matSuperMotif.rows();
    Eigen::VectorXf allSquaredNorm,allNorm,vectTempi,vectTempj;
	allNorm= matSuperMotif.matrix().rowwise().norm();
	allSquaredNorm= matSuperMotif.matrix().rowwise().squaredNorm(); //à retravaller

	//Cosine dissimilarity
	for (int i=0; i<n;i++)
	{
		for (int j=0; j<i;j++)
		{
            matDistS2S(i,j)=1-(matSuperMotif.row(i).matrix().dot(matSuperMotif.row(j).matrix())/(allNorm[i]*allNorm[j]));
			if(matDistS2S(i,j)<0){matDistS2S(i,j)=0;}
		}
	}
return matDistS2S;}
*/

Eigen::ArrayXXf RepSupernmotif::repSupermotifComputeMatrixDissimCosS2NM(Eigen::ArrayXXf matS_supernmotifs,Eigen::ArrayXXf matNM_supernmotifs){

	//compute Cosine dissimilarity between matS_supernmotifs and matNM_supernmotifs//
		int n=matS_supernmotifs.rows();
		int m=matNM_supernmotifs.rows();
        Eigen::MatrixXf matDistS2NM = Eigen::MatrixXf::Zero(matS_supernmotifs.rows(),matNM_supernmotifs.rows());
        Eigen::VectorXf allNormMatS_supernmotifs,allNormMatNM_supernmotifs,vectTempi,vectTempj;
		allNormMatS_supernmotifs= matS_supernmotifs.matrix().rowwise().norm();
		allNormMatNM_supernmotifs= matNM_supernmotifs.matrix().rowwise().norm();
		for (int i=0; i<n;i++)
		{
			for (int j=0; j<m;j++)
			{
				vectTempi=matS_supernmotifs.row(i);
				vectTempj=matNM_supernmotifs.row(j);
				matDistS2NM(i,j)=1-(vectTempi.dot(vectTempj)/(allNormMatS_supernmotifs[i]*allNormMatNM_supernmotifs[j]));
				if(matDistS2NM(i,j)<0){matDistS2NM(i,j)=0;}
			}
		}

return matDistS2NM;
}


