//============================================================================
// Name        : Supernmotifs program
// Author      : Jean-Pierre Sehi Glouzon
// Copyright   : GNU/GPL
// Description : Supernmotifs algorithm in C++, Ansi-style
//============================================================================

#include<string>
#include<stdio.h>
#include <iostream>
#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <vector>
#include <ctype.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <queue>
#include <stdexcept>
#include "RepMotif.h"
#include "RepWeightedMotif.h"
#include "RepSupermotif.h"
#include "Eigen/Dense"

using namespace std;

string i_PathFileInput;//input file
string o_PathDirOutput;//Output directory
int p_outputOption;// Output options
int nbSupNmotifs;//number of super nmotifs
int maxLevelNmotifs;//level of nmotifs
float minOcc;//minimum occurrence
int ngram;//ngram
string filenameTemp;

void initParameters(string&, string&, int&,int&,int&,float&,int&);
void setParameters(int, char*[] ,string&, string&,int&,int&,int&,float&,int&);

void writeMatALLMotif(const vector<string> ,map<string,float>, vector<map<string, float> > ,const string, const string );
void writeMatALLMotifWeighted(const vector<string>,const Eigen::ArrayXXf, map<string,float>,const string ,const string );
void writeMatS_Supernmotifs(const vector<string>,const Eigen::ArrayXXf,const string ,string );
void writeMatNmotifs_Supernmotifs( map<string,float> , const Eigen::ArrayXXf,const string ,const string );
void computeNwriteMatDissimS2S( const Eigen::ArrayXXf, const string ,const string);
void writeMatDissimS2SFull(const vector<string>, const Eigen::ArrayXXf, const string ,const string);
void writeMatDissimS2NM(const vector<string>, const Eigen::ArrayXXf, map<string,float>, const string ,const string);
void writeSingularVal(const Eigen::RowVectorXf,const string ,const string );
void writeMatNmotifsPositionsInSS(const vector<string> , vector< map<string,vector<int>> > , map<string,float>  ,const string , const string );
void writeStat(const vector<string>,map<string,float>,map<string,float>,const string, const string, const int,const Eigen::ArrayXXf );
string help();

int main(int argc, char *argv[])
{

    initParameters(i_PathFileInput, o_PathDirOutput,p_outputOption,nbSupNmotifs,maxLevelNmotifs,minOcc,ngram);
    setParameters(argc, argv, i_PathFileInput, o_PathDirOutput, p_outputOption,nbSupNmotifs, maxLevelNmotifs,minOcc,ngram);

    cout<<endl<<"Running the super n-motifs program..."<<endl;

    clock_t t1,t2,t3,t4;
    t1=clock();t2=clock();t3=clock();t4=clock();

    //Extract n-motifs.
    cout<<endl<<"Extract n-motifs..."<<endl;
    RepNmotif repnMotif(i_PathFileInput,maxLevelNmotifs,ngram);

    t1=clock()-t1;
    cout << "Took:" << ((float)t1)/CLOCKS_PER_SEC << " sec."<<endl;

    //Filter and weight n-motifs : n-motifs representation.
    cout<<"Filter and weight n-motifs to build the n-motifs representation ..."<<endl;
    RepWeightedMotif repFilterWeightNmotif(repnMotif.getNmotifsAllStructure(),repnMotif.getNmotifsForEachStructure(),minOcc);

    t2=clock()-t2;
    cout << "Took:" << ((float)t2)/CLOCKS_PER_SEC << " sec."<<endl;

    //Compute the super n-motifs representation and the dissimilarity matrix using the cosine dissimilarity.
    //Alternatively compute the dissimilarity matrix based on the n-motifs representation.
    if ((p_outputOption==0))
    {cout<<"Compute the super n-motif representation..."<<endl;}
    else if(p_outputOption==1) {cout<<"Compute the SS*super n-motifs matrix that is the super n-motif representation of SS..."<<endl;}
    else if(p_outputOption==2) {cout<<"Compute the SS*SS dissimilarity matrix, the SS*n-motifs dissimilarity matrix, "
                                      "the SS*super n-motifs matrix that is the super n-motif representation of SS and "
                                      "the n-motifs*super n-motifs matrix that is the super n-motifs representation of n-motifs  "<<endl;}
    else if(p_outputOption==3) {cout<<"Compute SS*SS dissimilarity matrix based on the n-motifs representation..."<<endl;}
    else if(p_outputOption==4) {cout<<"Compute SS*nm matrix that is the n-motif representation of SS..."<<endl;}
    else if(p_outputOption==5) {cout<<"Compute All... "<<endl;}
    RepSupernmotif repSupernmotif(repFilterWeightNmotif.getMatAllMotifWeigthed(),nbSupNmotifs,p_outputOption);

    t3=clock()-t3;
    cout << "Took:" << ((float)t3)/CLOCKS_PER_SEC << " sec."<<endl;

    cout<<endl<<"Compute and Write output files :"<<endl;

    if (p_outputOption==0)
    {
        //write the SS*SS dissimilarity matrix (SS:Secondary Structures)
        cout<<"SS*SS dissimilarity matrix, the lower part (matDissim_SSbySS.csv) ..."<<endl;
        filenameTemp="matDissim_SSbySS.csv";

        computeNwriteMatDissimS2S(repSupernmotif.getMatSSupernmotifs(),o_PathDirOutput,filenameTemp);

        //write statistics on data.
        cout<<"Statistics: nb. n-motifs, nb. relevantor filter n-motifs etc. (stat.csv) ..."<<endl;
        filenameTemp="stat.csv";
        writeStat(repnMotif.getHeaders(),repnMotif.getNmotifsAllStructure(),repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs(),o_PathDirOutput,filenameTemp,nbSupNmotifs,repSupernmotif.getMatSSupernmotifs());
        cout<<"Ending of writing."<<endl<<endl;

    }
    else if (p_outputOption==1)
    {
        //write the SS*supernmotifs matrix that is the super n-motif representation of SS
        cout<<"SS*supernmotifs matrix (matSnmRep_SSbySnm.csv) ..."<<endl;
        filenameTemp="matSnmRep_SSbySnm.csv";
        writeMatS_Supernmotifs(repnMotif.getHeaders(),repSupernmotif.getMatSSupernmotifs(),o_PathDirOutput,filenameTemp);
    }
    else if (p_outputOption==2)
    {
        //write the SS*SS dissimilarity matrix
        cout<<"SS*SS dissimilarity matrix, the lower part (matDissim_SSbySS.csv) ..."<<endl;
        filenameTemp="matDissim_SSbySS.csv";
        computeNwriteMatDissimS2S(repSupernmotif.getMatSSupernmotifs(),o_PathDirOutput,filenameTemp);

        //write the SS*n-motifs dissimilarity matrix
        cout<<"SS*n-motifs dissimilarity matrix (matDissim_SSbynm.csv) ..."<<endl;
        filenameTemp="matDissim_SSbynm.csv";
        writeMatDissimS2NM(repnMotif.getHeaders(), repSupernmotif.getMatDissimS2Nm(), repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs() , o_PathDirOutput,filenameTemp);

        //write the SS*super n-motifs matrix that is the super n-motif representation of SS
        cout<<"SS*supernmotifs matrix (matSnmRep_SSbySnm.csv) ..."<<endl;
        filenameTemp="matSnmRep_SSbySnm.csv";
        writeMatS_Supernmotifs(repnMotif.getHeaders(),repSupernmotif.getMatSSupernmotifs(),o_PathDirOutput,filenameTemp);

        //write the n-motifs*super n-motifs matrix that is the super n-motif representation of n-motifs
        cout<<"n-motifs*supernmotifs matrix (matSnmRep_nmbySnm.csv) ..."<<endl;
        filenameTemp="matSnmRep_nmbySnm.csv";
        writeMatNmotifs_Supernmotifs(repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs(),repSupernmotif.getMatNmSupernmotifs(),o_PathDirOutput,filenameTemp);

        //write the full singular values associated with the super n-motifs.
        cout<<"Full singular values (singularValuesFull_supernmotifs.csv) ..."<<endl;
        filenameTemp="singularValuesFull_supernmotifs.csv";
        writeSingularVal(repSupernmotif.getSingularValFull(),o_PathDirOutput,filenameTemp);

        //write the nucleotide position associates with n-motifs in the n-motifs representation of SS
        cout<<"n-motifsPosition matrix (matnmPos.csv) ..."<<endl;
        filenameTemp="matnmPos.csv";
        writeMatNmotifsPositionsInSS(repnMotif.getHeaders(),repnMotif.getNmotifsForEachStructureWithPosNucOfnmotifs(), repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs() ,o_PathDirOutput, filenameTemp);

    }
    else if (p_outputOption==3)
    {
        //write the SS*SS dissimilarity matrix based on the n-motifs representation
        cout<<"SS*SS dissimilarity matrix, the lower part, based on the n-motif representation (matDissimNmotifsBased_SSbySS.csv) ..."<<endl;
        filenameTemp="matDissimNmotifsBased_SSbySS.csv";
        computeNwriteMatDissimS2S(repFilterWeightNmotif.getMatAllMotifWeigthed(),o_PathDirOutput,filenameTemp);
    }
    else if (p_outputOption==4)
    {
        //write the SS*n-motifs matrix
        cout<<"SS*n-motifs matrix that is the n-motif representation of SS (matNmRep_SSbyNm.csv) ..."<<endl;
        filenameTemp="matNmRep_SSbyNm.csv";
        writeMatALLMotifWeighted(repnMotif.getHeaders(), repFilterWeightNmotif.getMatAllMotifWeigthed(), repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs(), o_PathDirOutput,filenameTemp);
    }
    else if (p_outputOption==5)
    {
        //write the SS*SS dissimilarity matrix (SS:Secondary Structures)
        cout<<"SS*SS dissimilarity matrix, the lower part (matDissim_SSbySS.csv) ..."<<endl;
        filenameTemp="matDissim_SSbySS.csv";
        computeNwriteMatDissimS2S(repSupernmotif.getMatSSupernmotifs(),o_PathDirOutput,filenameTemp);

        //write the full SS*SS dissimilarity matrix.
        cout<<"Full SS*SS dissimilarity matrix (matDissimFull_SSbySS.csv) ..."<<endl;
        filenameTemp="matDissimFull_SSbySS.csv";
        writeMatDissimS2SFull(repnMotif.getHeaders(),repSupernmotif.getMatDissimS2S(), o_PathDirOutput,filenameTemp);

        //write the SS*n-motifs dissimilarity matrix
        cout<<"SS*n-motifs dissimilarity matrix (matDissim_SSbynm.csv) ..."<<endl;
        filenameTemp="matSnmRep_SSbySnm.csv";
        writeMatDissimS2NM(repnMotif.getHeaders(), repSupernmotif.getMatDissimS2Nm(), repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs() , o_PathDirOutput,filenameTemp);

        //write the raw SS*n-motifs matrix.
        cout<<"Raw SS*n-motifs matrix (matNmRepRaw_SSbyNm.csv) ..."<<endl;
        filenameTemp="matNmRepRaw_SSbyNm.csv";
        writeMatALLMotif(repnMotif.getHeaders(), repnMotif.getNmotifsAllStructure(),repnMotif.getNmotifsForEachStructure(),o_PathDirOutput,filenameTemp);

        //write the SS*n-motifs matrix that is the n-motif representation of SS.
        cout<<"SS*n-motifs matrix that is the n-motif representation of SS (matNmRep_SSbyNm.csv) ..."<<endl;
        filenameTemp="matNmRep_SSbyNm.csv";
        writeMatALLMotifWeighted(repnMotif.getHeaders(), repFilterWeightNmotif.getMatAllMotifWeigthed(), repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs(), o_PathDirOutput,filenameTemp);

        //write the SS*supernmotifs matrix that is the super n-motif representation of SS
        cout<<"SS*super n-motifs matrix that is the super n-motif representation of SS (matSnmRep_SSbySnm.csv) ..."<<endl;
        filenameTemp="matSnmRep_SSbySnm.csv";
        writeMatS_Supernmotifs(repnMotif.getHeaders(),repSupernmotif.getMatSSupernmotifs(),o_PathDirOutput,filenameTemp);

        //write the full SS*supernmotifs matrix that is the super n-motif representation of SS
        cout<<"Full SS*super n-motifs matrix that is the super n-motif representation of SS (matSnmRepFull_SSbySnm.csv) ..."<<endl;
        filenameTemp="matSnmRepFull_SSbySnm.csv";
        writeMatS_Supernmotifs(repnMotif.getHeaders(),repSupernmotif.getMatSSupernmotifsFull(),o_PathDirOutput,filenameTemp);

        //write the n-motifs*super n-motifs matrix that is the super n-motif representation of n-motifs
        cout<<"n-motifs*supernmotifs matrix (matSnmRep_nmbySnm.csv) ..."<<endl;
        filenameTemp="matSnmRep_nmbySnm.csv";
        writeMatNmotifs_Supernmotifs(repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs(),repSupernmotif.getMatNmSupernmotifs(),o_PathDirOutput,filenameTemp);

        //write the full n-motifs*super n-motifs matrix that is the super n-motif representation of n-motifs
        cout<<"n-motifs*supernmotifs matrix (matSnmRepFull_nmbySnm.csv) ..."<<endl;
        filenameTemp="matSnmRepFull_nmbySnm.csv";
        writeMatNmotifs_Supernmotifs(repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs(),repSupernmotif.getMatNM_supernmotifsFull(),o_PathDirOutput,filenameTemp);

        //write the singular values associated with the super n-motifs.
        cout<<"Singular values (singularValues_supernmotifs.csv) ..."<<endl;
        filenameTemp="singularValues_supernmotifs.csv";
        writeSingularVal(repSupernmotif.getSingularVal(),o_PathDirOutput,filenameTemp);

        //write the full singular values associated with the super n-motifs.
        cout<<"Full singular values (singularValuesFull_supernmotifs.csv) ..."<<endl;
        filenameTemp="singularValuesFull_supernmotifs.csv";
        writeSingularVal(repSupernmotif.getSingularValFull(),o_PathDirOutput,filenameTemp);

        //write the nucleotide position associates with n-motifs in the n-motifs representation of SS
        cout<<"n-motifsPosition matrix (matnmPos.csv) ..."<<endl;
        filenameTemp="matnmPos.csv";
        writeMatNmotifsPositionsInSS(repnMotif.getHeaders(),repnMotif.getNmotifsForEachStructureWithPosNucOfnmotifs(), repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs() ,o_PathDirOutput, filenameTemp);

        //write statistics on data.
        cout<<"Statistics: nb. n-motifs, nb. relevantor filter n-motifs etc. (stat.csv) ..."<<endl;
        filenameTemp="stat.csv";
        writeStat(repnMotif.getHeaders(),repnMotif.getNmotifsAllStructure(),repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs(),o_PathDirOutput,filenameTemp,nbSupNmotifs,repSupernmotif.getMatSSupernmotifs());
        cout<<"Ending of writing."<<endl<<endl;
    }
    cout<<"Execution of super n-motifs program completed"<<endl;
     t4=clock()-t4;
    cout << "Took:" << ((float)t4)/CLOCKS_PER_SEC << " sec."<<endl;
    return 0;
}

/******************************************************************************/
void initParameters(string& i_PathFileInput,string& o_PathDirOutput, int& p_outputOption,int& nbSupNmotifs,int& maxLevelNmotifs, float& minOcc, int& ngram){
    i_PathFileInput=""; //input file
    o_PathDirOutput=""; //Output directory
    p_outputOption=0; //Output the dissimilarity matrix
    nbSupNmotifs=0; //number of super n-motifs -->0: automatic determination of dimension: brocken stick model
    maxLevelNmotifs=1; //n-motifs parameters
    minOcc=0;//automatic minimum support
    ngram=0;
}

/******************************************************************************/
void setParameters(int argc,char* argv[],string& i_PathFileInput,string& o_PathDirOutput,int& p_outputOption, int& nbSupNmotifs,int& maxLevelNmotifs, float& minOcc, int& ngram){

    ifstream viennaFile;
    //string filename="helpTerminal.txt";
    //string filename="README.md";

    //ifstream helpFile (filename.c_str());

    string lineHelpFile;
    bool paramRequiredinput=false;
    bool paramRequiredoutput=false;

    struct stat sb;
    for (int i=0;i<argc;i++){
        if (argv[i][0]=='-'){
            switch ( argv[i][1] ) {
            case 'i':
                if (argv[i+1]!=nullptr)
                {
                    i_PathFileInput=string(argv[i+1]);
                    viennaFile.open(i_PathFileInput.c_str());
                      if (!(viennaFile.is_open()))
                      {throw invalid_argument("Unable to open vienna file. Please check the input path file.");}
                    viennaFile.close();
                    paramRequiredinput=true;
                }
                else
                {   throw invalid_argument("Empty value for parameter -i. Please enter the input path file.");}
            break;
            case 'o':
                if (argv[i+1]!=nullptr)
                {
                 o_PathDirOutput=string(argv[i+1]);
                 if (!(stat(o_PathDirOutput.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)))
                  {
                    cerr<<"The ouput path directory doesn't exits. It will be created."<<endl;
                    #if defined(_WIN32)
                        mkdir(o_PathDirOutput.c_str());
                         #else
                        mkdir(o_PathDirOutput.c_str(), 0700);
                         #endif
                  }
                 paramRequiredoutput=true;
                }
                else
                {   throw invalid_argument("Empty value for parameter -o. Please enter the output path directory.");}
              break;
            case 'p':
                if (argv[i+1]!=nullptr)
                {
                    p_outputOption=stoi(argv[i+1]);
                    if (!(p_outputOption==0) && !(p_outputOption==1) && !(p_outputOption==2)&& !(p_outputOption==3)&& !(p_outputOption==4)&& !(p_outputOption==5) )
                    {throw invalid_argument("Output option (-p) must be "
                                            "0 to output the SS*SS dissimilarity matrix (matDissim_SSbySS.csv), "
                                            "1 to output the SS*super n-motifs matrix that is the super n-motif representation of SS (matSnmRep_SSbySnm.csv), "
                                            "2 to output the SS*SS dissimilarity matrix (matDissim_SSbySS.csv), the SS*n-motifs dissimilarity matrix (matDissim_SSbynm.csv), "
                                            "the SS*super n-motifs matrix that is the super n-motif representation of SS (matSnmRep_SSbySnm.csv), "
                                            "the n-motifs*super n-motifs matrix that is the super n-motifs representation of n-motifs (matSnmRep_nmbySnm.csv),  "
                                            "3 to output the SS*SS dissimilarity matrix based on the n-motif representation (matDissimNmotifsBased_SSbySS.csv), "
                                            "4 to output the SS*n-motifs matrix that is the n-motif representation of SS (matNmRep_SSbyNm.csv), "
                                            "or 5 to get all outputs with their full version including the singular values associated with the super n-motifs (singularValues_supernmotifs.csv) "
                                            "and  some statistics (raw nb. n-motifs, nb. n-motifs etc.) (stat.csv).");}
                }
                else
                {   throw invalid_argument("Empty value for parameter -p. "
                                           "Output option (-p) must be "
                                           "0 to output the SS*SS dissimilarity matrix (matDissim_SSbySS.csv), "
                                           "1 to output the SS*super n-motifs matrix that is the super n-motif representation of SS (matSnmRep_SSbySnm.csv), "
                                           "2 to output the SS*SS dissimilarity matrix (matDissim_SSbySS.csv), the SS*n-motifs dissimilarity matrix (matDissim_SSbynm.csv), "
                                           "the SS*super n-motifs matrix that is the super n-motif representation of SS (matSnmRep_SSbySnm.csv), "
                                           "the n-motifs*super n-motifs matrix that is the super n-motifs representation of n-motifs (matSnmRep_nmbySnm.csv),  "
                                           "3 to output the SS*SS dissimilarity matrix based on the n-motif representation (matDissimNmotifsBased_SSbySS.csv), "
                                           "4 to output the SS*n-motifs matrix that is the n-motif representation of SS (matNmRep_SSbyNm.csv), "
                                           "or 5 to get all outputs with their full version including the singular values associated with the super n-motifs (singularValues_supernmotifs.csv) "
                                           "and  some statistics (raw nb. n-motifs, nb. n-motifs etc.) (stat.csv).");}
              break;

            case 'n':
                if (argv[i+1]!=nullptr)
                {
                    maxLevelNmotifs=stoi(argv[i+1]);
                    if (!(maxLevelNmotifs==0)&&!(maxLevelNmotifs==1)&& !(maxLevelNmotifs==2))
                    {throw invalid_argument("Maximum level of n-motifs parameter (-n) is used to extract the n-motifs. It must be 0, 1 or 2. "
                                            "For instance: when it is set to 0, 0-motifs will be extracted. "
                                            "when it is set to 1, 0-motifs and 1-motifs will be extracted. etc.");}
                }
                else
                {   throw invalid_argument("Empty value for parameter -n. "
                                           "The maximum level of n-motifs parameter (-n) is used to extract the n-motifs. It must be 0, 1 or 2."
                                           "For instance: when it is set to 0, 0-motifs will be extracted. "
                                           "when it is set to 1, 0-motifs and 1-motifs will be extracted. etc.");};
              break;
            case 's':
                if (argv[i+1]!=nullptr)
                {

                    nbSupNmotifs=stoi(argv[i+1]);
                    if (!(nbSupNmotifs>=2)&&!(nbSupNmotifs==0))
                    { throw invalid_argument("Super n-motifs parameter (-s) must be s=0 or s>=2 respectively for the automatic determination of the number of super-n-motifs or retaining the first s super-n-motifs.");}
                }
                else
                {   throw invalid_argument("Empty value for parameter -s. "
                                           "Super n-motifs parameter (-s) must be s=0 or s>=2 respectively for the automatic determination of the number of super-n-motifs or retaining the first s super-n-motifs.");};
              break;

            case 'm':
                if (argv[i+1]!=nullptr)
                {
                    minOcc=stod(argv[i+1]);
                    if (!(minOcc>=0))
                    {throw invalid_argument("Minimum occurrence of motifs parameter (-m) must be positive (>=0)."
                                            "N-motifs with occurrence below the minimum occurrence are removed."
                                            "When it is set to 0, it computes the automatic n-motif minimum occurrence.");}
                }
                else
                {   throw invalid_argument("Empty value for parameter -m."
                                           "Minimum occurrence of motifs parameter (-m) must be positive (>=0). "
                                           "N-motifs with occurrence below the minimum occurrence are removed."
                                           "When it is set to 0, it computes the automatic n-motif minimum occurrence.");};

              break;
            case 'g':
                if (argv[i+1]!=nullptr)
                {
                    ngram=stoi(argv[i+1]);
                    if (!(ngram==0) && !(ngram==1) && !(ngram==2)&& !(ngram==3)&& !(ngram==4)&& !(ngram==5))
                    {throw invalid_argument("Ngram parameter (-g) must be : 0, 1, 2 , 3, 4 or 5.");}
                }
                else
                {   throw invalid_argument("Empty value for parameter -g. "
                                           "Ngram parameter (-g) must be : 0, 1, 2 , 3, 4 or 5.");};

              break;
            case 'h': //print help and exit the program
                  //if (helpFile.is_open())
                  //{
                  //  while ( getline (helpFile,lineHelpFile) )
                  //  { cout << lineHelpFile << '\n';}
                  //  helpFile.close();
                  //}
                  //else cout << "Unable to open file"<<endl;
                  cout<< help();

                exit(EXIT_SUCCESS);
              break;
            default:
                cerr<<"No parameters found. The program requires minimally the input vienna file (-i) and the ouput folder (-o). For help, type -h."<<endl;
                exit(EXIT_FAILURE);
            break;
            }
            if ((argv[i][1]!='h')&&(argv[i][1]!='i')&&(argv[i][1]!='o')&&(argv[i][1]!='p')&&(argv[i][1]!='n')&&(argv[i][1]!='s')&&(argv[i][1]!='m')&&(argv[i][1]!='g'))
            {cerr<<"Wrong parameters. Please check parameter list in the help (-h).";exit(EXIT_FAILURE);}
        }
    }
    if((paramRequiredinput==false || paramRequiredoutput==false))
        {cerr<<"The program requires the input vienna file (-i) and the ouput folder (-o)."<<endl;exit(EXIT_FAILURE);}
}

/******************************************************************************/
void writeMatALLMotif(const vector<string> headers, map<string,float> allStructureFeature, vector<map<string, float> > nmotifsForEachStructure  ,const string o_PathDirOutput, const string filename)
{
    ofstream outputMatMotifs;
    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    outputMatMotifs.open(tempfilename.c_str());

    Eigen::ArrayXXi matS_nmotifs = Eigen::ArrayXXi::Zero(nmotifsForEachStructure.size(), allStructureFeature.size());
    std::map<string, float>::iterator it;
    int ind = 0;
    for (unsigned int i = 0; i < nmotifsForEachStructure.size(); i++) {
        for (std::map<string, float>::iterator it2 =nmotifsForEachStructure[i].begin();it2 != nmotifsForEachStructure[i].end(); ++it2)
        {
            it = allStructureFeature.find(it2->first);
            ind = distance(allStructureFeature.begin(), it);
            matS_nmotifs(i, ind) = it2->second;
        }
    }

    //write labels
    map<string,float>::iterator it1;
    outputMatMotifs<<" ";
    for(it1=allStructureFeature.begin();it1!=allStructureFeature.end();++it1)
    {outputMatMotifs<<","<< it1->first;}
    outputMatMotifs<<"\n";

    //write data
    Eigen::IOFormat Comma(Eigen::StreamPrecision, 0, ",","","","");
    for (unsigned int i=0;i<matS_nmotifs.rows();i++)
    {
     outputMatMotifs<<headers[i]<<",";
     outputMatMotifs<<matS_nmotifs.row(i).format(Comma)<<"\n";
    }

    outputMatMotifs.close();
}

/******************************************************************************/
void writeMatALLMotifWeighted(const vector<string> headers, const Eigen::ArrayXXf matALLMotifWeighted, map<string,float> allStructureFeatureForWeigthOfMotifs ,const string o_PathDirOutput, const string filename){

    ofstream outputMatMotifs;
    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    outputMatMotifs.open(tempfilename.c_str());

    //write n-motifs labels
    map<string,float>::iterator it;
    outputMatMotifs<<" ";
    for(it=allStructureFeatureForWeigthOfMotifs.begin();it!=allStructureFeatureForWeigthOfMotifs.end();++it)
    {outputMatMotifs<<"," << it->first;}
    outputMatMotifs<<"\n";

    //write data
    Eigen::IOFormat Comma(Eigen::StreamPrecision, 0, ",","\n","","");
    for (unsigned int i=0;i<matALLMotifWeighted.rows();i++)
    {outputMatMotifs<<headers[i]<<","<<matALLMotifWeighted.row(i).format(Comma)<<endl;}

    outputMatMotifs.close();
}

/*********************************************************************************************/
void writeMatS_Supernmotifs(const vector<string> headers, const Eigen::ArrayXXf matSuperMotif,const string o_PathDirOutput,const string filename){

    ofstream outputMatMotifs;
    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    outputMatMotifs.open(tempfilename.c_str());

    outputMatMotifs<<" ";
    for(int i=0;i<matSuperMotif.cols();i++)
    {outputMatMotifs<<","<< "snm_"<<i+1;}
    outputMatMotifs<<"\n";

    Eigen::IOFormat Comma(Eigen::StreamPrecision, 0, ",","\n","","");
    for (unsigned int i=0;i<matSuperMotif.rows();i++)
    {outputMatMotifs<<headers[i]<<","<<matSuperMotif.row(i).format(Comma)<<endl;}

    outputMatMotifs.close();
}

/*********************************************************************************************/
void writeMatNmotifs_Supernmotifs( map<string,float> allStructureFeatureForWeigthOfMotifs , const Eigen::ArrayXXf matnm_SuperMotif,const string o_PathDirOutput,const string filename){

    ofstream outputMatMotifs;
    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    outputMatMotifs.open(tempfilename.c_str());

    outputMatMotifs<<" ";
    for(int i=0;i<matnm_SuperMotif.cols();i++)
    {outputMatMotifs<<","<< "snm_"<<i+1;}
    outputMatMotifs<<"\n";

    map<string,float>::iterator it;
    Eigen::IOFormat Comma(Eigen::StreamPrecision, 0, ",","\n","","");
    int i=0;
    for(it=allStructureFeatureForWeigthOfMotifs.begin();it!=allStructureFeatureForWeigthOfMotifs.end();++it)
    {
      outputMatMotifs<<it->first<<","<<matnm_SuperMotif.row(i).format(Comma)<<endl;
      i++;
    }
    outputMatMotifs.close();
}


/*********************************************************************************************/
void computeNwriteMatDissimS2S(const Eigen::ArrayXXf matSuperMotif,const string o_PathDirOutput,const string filename){

    ofstream outputMatMotifs;
    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    outputMatMotifs.open(tempfilename.c_str());

    //for (int i=0; i<matDissimS2S.rows();i++)
    //{
    //   for (int j=0; j<i;j++)
    //       {outputMatMotifs<<matDissimS2S(i,j)<<",";}
    //  outputMatMotifs<<"\n";
    //}

    float tempDist=0;
	int n=matSuperMotif.rows();
    Eigen::VectorXf allSquaredNorm,allNorm,vectTempi,vectTempj;
	allNorm= matSuperMotif.matrix().rowwise().norm();
	allSquaredNorm= matSuperMotif.matrix().rowwise().squaredNorm(); //à retravaller

	/*Cosine dissimilarity*/
	for (int i=0; i<n;i++)
	{
		for (int j=0; j<i;j++)
		{
            tempDist=1-(matSuperMotif.row(i).matrix().dot(matSuperMotif.row(j).matrix())/(allNorm[i]*allNorm[j]));
			if(tempDist<0){tempDist=0; }
			outputMatMotifs<<tempDist<<",";
		}
		outputMatMotifs<<"\n";
	}

    outputMatMotifs.close();
}

/*********************************************************************************************/
void writeMatDissimS2SFull(const vector<string> headers, const Eigen::ArrayXXf matDissimS2S,const string o_PathDirOutput,const string filename){

    ofstream outputMatMotifs;
    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    outputMatMotifs.open(tempfilename.c_str());

    //write n-motifs labels
    outputMatMotifs<<" ";
    for(unsigned int i=0;i<headers.size();i++)
    {outputMatMotifs<<","<< headers[i];}
    outputMatMotifs<<"\n";

    //write data
    Eigen::ArrayXXf matDissimS2STempURight=matDissimS2S.matrix().transpose();
    Eigen::ArrayXXf matDissimS2STempULeft=matDissimS2S;
    matDissimS2STempULeft=matDissimS2STempULeft+matDissimS2STempURight;
    Eigen::IOFormat Comma(Eigen::StreamPrecision, 0, ",","\n","","");
    for (unsigned int i=0;i<matDissimS2STempULeft.rows();i++)
    {outputMatMotifs<<headers[i]<<","<<matDissimS2STempULeft.row(i).format(Comma)<<endl;}

    outputMatMotifs.close();
}

/*********************************************************************************************/
void writeMatDissimS2NM(const vector<string> headers, const Eigen::ArrayXXf matDissimS2NM, map<string, float> allStructureFeatureForWeigthOfMotifs , const string o_PathDirOutput, const string filename){

    ofstream outputMatMotifs;
    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    outputMatMotifs.open(tempfilename.c_str());

    //write n-motifs labels
    map<string,float>::iterator it;
    outputMatMotifs<<" ";
    for(it=allStructureFeatureForWeigthOfMotifs.begin();it!=allStructureFeatureForWeigthOfMotifs.end();++it)
    {outputMatMotifs<<","<< it->first;}
    outputMatMotifs<<"\n";

    //write data
    Eigen::IOFormat Comma(Eigen::StreamPrecision, 0, ",","\n","","");
    for (unsigned int i=0;i<matDissimS2NM.rows();i++)
    {outputMatMotifs<<headers[i]<<","<<matDissimS2NM.row(i).format(Comma)<<endl;}

    outputMatMotifs.close();
}

/******************************************************************************/
void writeSingularVal( const Eigen::RowVectorXf singularVal,const string o_PathDirOutput,const string filename){

    ofstream outputMatMotifs;
    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    outputMatMotifs.open(tempfilename.c_str());

    //write labels
    for(unsigned i=0;i<singularVal.cols();i++)
    {outputMatMotifs <<" ,snm_"<<i+1<<"";}
    outputMatMotifs<<"\n";

    outputMatMotifs<<" ,";
    Eigen::IOFormat Comma(Eigen::StreamPrecision, 0, ",","\n","","");
    outputMatMotifs<<singularVal.format(Comma)<<endl;

    outputMatMotifs.close();
}

/******************************************************************************/
void writeMatNmotifsPositionsInSS(const vector<string> headers, vector< map<string,vector<int>> > FeatureForEachStructureWithPosNuc, map<string,float> allStructureFeatureForWeigthOfMotifs ,const string o_PathDirOutput, const string filename){


    ofstream outputMatMotifs;
    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    outputMatMotifs.open(tempfilename.c_str());

    //write n-motifs labels
    map<string,float>::iterator it;
    outputMatMotifs<<" ";
    for(it=allStructureFeatureForWeigthOfMotifs.begin();it!=allStructureFeatureForWeigthOfMotifs.end();++it)
    {outputMatMotifs<<"," << it->first;}
    outputMatMotifs<<"\n";

    //write data
    map<string,vector<int>>::iterator it2;
    for(unsigned int i=0;i<headers.size();i++)
    {
        outputMatMotifs<<headers[i]<<",";
        for(it=allStructureFeatureForWeigthOfMotifs.begin();it!=allStructureFeatureForWeigthOfMotifs.end();++it)
        {

            it2=FeatureForEachStructureWithPosNuc[i].find(it->first);
            if (it2!=FeatureForEachStructureWithPosNuc[i].end())
            {
                for(unsigned int j=0;j<it2->second.size();j++)
                {
                    outputMatMotifs<<it2->second[j]+1<<"|";
                }
            }
            else
            {
                outputMatMotifs<<"x";
            }
            outputMatMotifs<<",";
        }
        outputMatMotifs<<"\n";
    }

    outputMatMotifs.close();
}

/************************************************************************************/
void writeStat(const vector<string> headers, map<string, float> allStructureFeature, map<string, float> allStructureFeatureForWeigthOfMotifs, const string o_PathDirOutput, const string filename, const int nbSupNmotifs, const Eigen::ArrayXXf repSupernmotif){

    ofstream outputMatMotifs;
    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    outputMatMotifs.open(tempfilename.c_str());

    if (nbSupNmotifs==0)
   {
       outputMatMotifs<<"Nb. of Structures, Raw nb. of nmotifs, Nb. of relevant n-motifs, Nb. of automatically determined super-n-motifs "<<endl;


       outputMatMotifs<<headers.size()<<", "<<allStructureFeature.size()<<", "<<allStructureFeatureForWeigthOfMotifs.size()<<", "<<repSupernmotif.cols()<<endl;
   }
   else{
       outputMatMotifs<<"Nb. of Structures, Raw nb. of nmotifs, Nb. of relevant n-motifs, Nb. of selected super-n-motifs "<<endl;
       outputMatMotifs<<headers.size()<<", "<<allStructureFeature.size()<<", "<<allStructureFeatureForWeigthOfMotifs.size()<<", "<<repSupernmotif.cols()<<endl;
   }

    outputMatMotifs.close();
}

string help()
{
    //From README.md
    string help;

    help=
"\n### Super n-motifs model ###\n"
"Usage:\n"
"\n"
"* supernmotifs [Parameters]...\n"
"\n"
"Examples :\n"
"\n"
"* Execute super n-motifs in command line\n"
"\n"
"  ./pathOfSupernmotifsProgram/supernmotifs -i /PathtoDbFile\n"
"  -o /OutputDirectoryPath/\n"
"\n"
"Important notes:\n"
"\n"
"* Circular RNA\n"
"\n"
"  For the processing of secondary structures of circular RNA, adding 'c_'\n"
"   at the beginning of the header of each circular RNA is required.\n"
"\n"
"* Pseudoknots and Gquadruplexes\n"
"\n"
"  Special characters '{}','<>','[]', and alphabets such as 'Aa','Bb','Zz'\n"
"  are used to represent base pairs involved in pseudoknots. '+' is used to\n"
"  represent each guanine involved in the Gquadraplexes formation.\n"
"\n"
"### Parameters ###\n"
"\n"
"  -h\n"
"\n"
"	Print help and exit the program.\n"
"\n"
"  -i [input vienna file]\n"
"\n"
"    Input file of RSS in vienna/db format (required).\n"
"    >strucID\n"
"    AAAAAUU\n"
"    ((...))\n"
"\n"
"  -o [output folder]\n"
"\n"
"    Results folder (required).\n"
"\n"
"  -n [0|1|2]\n"
"\n"
"    Specify the maximum level of n-motifs used to\n"
"    extract the n-motifs. It must be 0, 1 or 2. For instance: when\n"
"    it is set to 0, 0-motifs will be extracted. When it is set\n"
"    to 1, 0-motifs and 1-motifs will be extracted. etc.\n"
"    (default : -n 1)\n"
"\n"
"  -s [s=0|s>=2]\n"
"\n"
"    Specify  the number of super n-motifs.\n"
"    (default : -s 0 for automatic determination of the number of\n"
"    super n-motifs using the broken stick model)\n"
"\n"
"  -m [m>=0]\n"
"\n"
"    Specify the minimum occurrence of n-motifs. N-motifs with occurrence\n"
"    below the minimum occurrence are removed. When it is set to 0,\n"
"    it computes the automatic n-motif minimum occurrence.\n"
"    (default : -m 0)\n"
"\n"
"  -g [0|1|2|3|4|5]\n"
"\n"
"    Specify the length of sequence patterns (ngrams).\n"
"    (default : -g 0)\n"
"\n"
"  -p [0|1|2|3|4]\n"
"\n"
"     Specify the output options.\n"
"     (default: -p 0)\n"
"     '-p 0' to output the SS*SS dissimilarity matrix (matDissim_SSbySS.csv)\n"
"     and some statistics (stat.csv)\n"
"     '-p 1' to output the SS*super n-motifs matrix that is the super n-motif\n"
"     representation of SS (matSnmRep_SSbySnm.csv)\n"
"     '-p 2' to output the SS*SS dissimilarity matrix (matDissim_SSbySS.csv),\n"
"     the SS*n-motifs dissimilarity matrix (matDissim_SSbynm.csv), the\n"
"     SS*supern-motifs matrix that is the super n-motif representation of SS\n"
"     (matSnmRep_SSbySnm.csv), the n-motifs*super n-motifs matrix that is the\n"
"     super n-motifs representation of n-motifs (matSnmRep_nmbySnm.csv),\n"
"     the singular values (singularValuesFull_supernmotifs.csv) and the matrix\n"
"     associating with each n-motifs its position in each SS (matnmPos.csv)\n"
"     '-p 3' to output the SS*SS dissimilarity matrix based on the n-motif\n"
"     representation (matDissimNmotifsBased_SSbySS.csv)\n"
"     '-p 4' to output the SS*n-motifs matrix that is the n-motif representation\n"
"     of SS (matNmRep_SSbyNm.csv)\n"
"     '-p 5' to get all outputs including the full version of all the matrices\n\n";

return help;
}
