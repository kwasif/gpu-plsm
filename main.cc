/*
 * Copyright (c) 2013 - Idiap Research Institute - http://www.idiap.ch/
 * 
 * See file COPYING for the licence associated to this software
 * 
 * This file is a part of the TopicModels4Video software version: academic-1.0
*/
#include "CmdLine.h"
#include <iostream>
#include <string> //string, getline()
#include <fstream>
#include <sstream>
#include "TAM.h"

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <limits>
#include <sys/stat.h>
#include <cmath>

using namespace std;

bool userTouchedTheStopFileForSince(const char *datafile, const time_t &since) {
    char *file = new char[strlen(datafile) + 10];
    sprintf(file, "%s.stop", datafile);
    bool res = false;
    struct stat b;
    if (!stat(file, &b)) {
        time_t &modtime = b.st_mtime;
        res = modtime - since > 0;
    }
    delete[] file;
    return res;
}

int main(int argc, char * argv[]) {
    const time_t startupDate = time(NULL);

    int seed = (int) time(NULL);

    int docLength;
    char *datafile;
    int n_Aspects, nTz;
    float end_accuracy;
    int n_iter;
    char *modelfile = NULL;
    char *initfile = NULL;
    char *initfilePtrwz = NULL;
    float initscalePtrwz = 1.f;
    int initype = 1;
    float scaleOnRead = 1.f;
    int multiRunCount = 1;

    int trainFlag;
    struct Options opt;
    opt.z_prior = false;
    opt.tr_prior= false;
    opt.prune_topics = false;
    opt.z_thresh = 0.95;
    bool d_overlap, read_several;
    bool align_topics = false;
    time_t start,end;

    char *writePerInstantLL; // will be set to "" by default by the cmd
    
    Torch::CmdLine cmd;
    cmd.info("Temporal Aspect Modeling.\n");

    cmd.addText("\nArguments:");

    cmd.addSCmdArg("Data file", &datafile, "data file");
    cmd.addICmdArg("DocumentLength", &docLength, "no of frames in a Document (use -1 if unknown)");
    cmd.addICmdArg("Topics", &n_Aspects, "no of topics");
    cmd.addICmdArg("nTz", &nTz, "Topic duration");
    cmd.addICmdArg("train flag", &trainFlag, "train(1)/infer(0)");

    cmd.addText("\n options");
    cmd.addICmdOption("-seed", &seed, seed, "Force random seed");
    cmd.addICmdOption("-multi", &multiRunCount, multiRunCount, "Do multiple runs and keep the best (number of runs)");
    cmd.addBCmdOption("-inList", &read_several, false, "Flag to indicate whether the input file contains a list of filenames/not[false]");
    cmd.addBCmdOption("-do", &d_overlap, false, "Flag to indicate document Overlap required/not[false]");
    cmd.addRCmdOption("-sor", &scaleOnRead, scaleOnRead, "scale factor when reading document (before converting the count to an integer)");
    cmd.addSCmdOption("-mo", &modelfile, "", "model files to obtain pwz and ptrwz");
    cmd.addSCmdOption("-init", &initfile, "", "initialization to obtain pwz and ptszd");
    cmd.addSCmdOption("-initPtrwz", &initfilePtrwz, "", "initialization to from existing pwz");
    cmd.addRCmdOption("-initPtrwzScaling", &initscalePtrwz, initscalePtrwz, "scaling factor for the -initPtrwz file (weight of the dirichlet prior)");
    cmd.addRCmdOption("-acc", &end_accuracy, 0.0001, "convergence criteria[0.0001]");
    cmd.addICmdOption("-iter", &n_iter, 50, "No of iterations[50]");
    cmd.addRCmdOption("-lts", &opt.lambdaTsSparsity, 0.0, "Lambda factor for the p(ts) thresholding[0.0]");
    cmd.addRCmdOption("-trwt", &opt.tr_wt, 1.0, "Weight for tr prior while learning[1.0]");
    cmd.addRCmdOption("-zth", &opt.z_thresh, 0.95, "Threshold for prunning topics[0.95]");
    cmd.addSCmdOption("-writeTaLL", &writePerInstantLL, "", "write the per instant loglikelihood to this file");

    //bools
    cmd.addBCmdOption("-prn", &opt.prune_topics, false, "Prune topics using a threshold[0.95] for inference[false]");
    cmd.addBCmdOption("-zp", &opt.z_prior, false, "Flag to use or not z prior in inference[false]");
    cmd.addBCmdOption("-trp", &opt.tr_prior, false, "Flag to use or not tr prior in learning[false]");
    cmd.addBCmdOption("-alignz", &align_topics, false, "Flag to use or not tr prior in learning[false]");

    cmd.addText("\n Fit modes:");
    cmd.addText("       0: table");
    cmd.addText("       1: gaussian (fit1: stddevToAdd, fit2: addUniformProportion)");
    cmd.addText("       2: gaussian/uniform (fit1: stddevToAdd, fit2: initialPi_u)");
    cmd.addICmdOption("-fitMode", (int*)&opt.fitMode, (int)opt.fitMode, "Fit Mode");
    cmd.addRCmdOption("-fit1", &opt.fitParam1, opt.fitParam1, "Fit parameter #1");
    cmd.addRCmdOption("-fit2", &opt.fitParam2, opt.fitParam2, "Fit parameter #2");



    cmd.read(argc, argv);
    //printf("opt.prune_topics: %s\n",(opt.prune_topics)?"true":"false");
    //printf("opt.z_prior: %s\n",(opt.z_prior)?"true":"false");
    //printf("opt.tr_prior: %s\n",(opt.tr_prior)?"true":"false");


    printf("data file: %s \n",datafile);
    printf("Topics: %d \n",n_Aspects);
    printf("nTz: %d \n",nTz);
    printf("provided doc length: %d \n",docLength);
    //printf("b_train: %d \n",b_train);*/

    if (!trainFlag) {
        if (!strcmp(modelfile, "")) {
            printf("please provide the model files \n");
            exit(0);
        }
    }
    if (strcmp(initfile, "") && trainFlag) {
        cout << strcmp(initfile, "") << "\n";
        initype = 2;
        printf("initype: %d %s \n", initype, initfile);
    }
    if (strcmp(initfilePtrwz, "") && trainFlag) {
        cout << strcmp(initfile, "") << "\n";
        initype = 3;
        printf("initype: %d %s \n", initype, initfilePtrwz);
    }


    Corpus *data = new Corpus(0);
    Document src_doc, doc;
    Frame frame;

    int length, count, word, actualDocLength = 0, currentDocLength = 0;

    if (read_several) {// Data loading
        cout << "Houston" << endl;
        string linelist;
        ifstream infilelist(datafile, ios_base::in);
	while (getline(infilelist, linelist, '\n')) {
            if (linelist[0] == '#')
                continue;
            // Trim
            linelist.erase (linelist.find_last_not_of (" \t\r\n") + 1);
            linelist.erase (0, linelist.find_first_not_of (" \t\r\n"));
            
	    string line;
	    ifstream infile(linelist.c_str (), ios_base::in);
	    while (getline(infile, line, '\n')) {
		if (line[0] == '#')
		    continue;
		size_t firstSpace = line.find(" ");
		size_t firstColon = line.find(":");
		istringstream inLine(line);
		if (firstSpace == string::npos || firstSpace > firstColon) {
		    length = std::count(line.begin(), line.end(), ':');
		} else {
		    // old style
		    inLine >> length;
		}
		for (int n = 0; n < length; n++) {
		    int w;
		    char c;
		    inLine >> word;
		    inLine.get(c);
		    if (c != ':') {
			printf
			    ("Error in file interpretation (expected ':' but got '%c')\n",
			     c);
			throw "File interpretation error";
		    }
		    float fcount;
		    inLine >> fcount;
		    count = (int) (fcount * scaleOnRead);
		    if (!inLine.eof() && inLine.get(c) && c != ' ') {
			printf
			    ("Error in file interpretation (expected ' ' but got '%c')\n",
			     c);
			throw "File interpretation error";
		    }
		    if (count != 0 && !isnan(count)) {
			frame.push_back(WordCount(word, count));
			data->m_totWords += count;
			if (word >= data->m_corpusLength) {
			    data->m_corpusLength = word + 1;
			}
		    }
		}
		actualDocLength++;
                currentDocLength++;
		src_doc.push_back(frame);
		frame.clear();
	    }
	    data->add (src_doc);
            data->m_docLength.push_back (currentDocLength);
            //cout << currentDocLength << " ";
            currentDocLength = 0;
	    src_doc.clear ();
	    data->m_numDocs++;
	}
        //cout << endl;
        //for (vector<int>::iterator it = data->m_docLength.begin() ; it != data->m_docLength.end() ; it++)
        //    cout << (*it) << " ";
        //cout << endl;
    } else {
        cout << "OK" << endl;
	{			// Data loading (from a single file)
	    string line;
	    ifstream infile(datafile, ios_base::in);
	    while (getline(infile, line, '\n')) {
		if (line[0] == '#')
		    continue;
		size_t firstSpace = line.find(" ");
		size_t firstColon = line.find(":");
		istringstream inLine(line);
		if (firstSpace == string::npos || firstSpace > firstColon) {
		    length = std::count(line.begin(), line.end(), ':');
		} else {
		    // old style
		    inLine >> length;
		}
		for (int n = 0; n < length; n++) {
		    int w;
		    char c;
		    inLine >> word;
		    inLine.get(c);
		    if (c != ':') {
			printf
			    ("Error in file interpretation (expected ':' but got '%c')\n",
			     c);
			throw "File interpretation error";
		    }
		    float fcount;
		    inLine >> fcount;
		    count = (int) (fcount * scaleOnRead);
		    if (!inLine.eof() && inLine.get(c) && c != ' ') {
			printf
			    ("Error in file interpretation (expected ' ' but got '%c')\n",
			     c);
			throw "File interpretation error";
		    }
		    if (count != 0 && !isnan(count)) {
			frame.push_back(WordCount(word, count));
			data->m_totWords += count;
			if (word >= data->m_corpusLength) {
			    data->m_corpusLength = word + 1;
			}
		    }
		}
		actualDocLength++;
		src_doc.push_back(frame);
		frame.clear();
	    }
	}

        if (docLength < 0)
            docLength = actualDocLength;
	// split doc
	int f_count = 0;
	int n_line = 0;
	data->m_numDocs = 0;
	while (f_count < src_doc.size()) {
	    //cout<<data->m_numDocs<<"\n";
	    if ((n_line + 1) % docLength == 0) {
		doc.push_back(src_doc.at(f_count));
		f_count++;
		n_line++;
		data->add(doc);
                data->m_docLength.push_back (docLength);
		doc.clear();
		data->m_numDocs++;
		n_line = 0;
		if (d_overlap) {
		    f_count = f_count - nTz + 1;
		    //cout<<"Document Overlap :"<<d_overlap<<endl;
		}

	    }			// push the frame into the document
	    else {

		doc.push_back(src_doc.at(f_count));
		n_line++;
		f_count++;
		frame.clear();
		//cout<<f_count<<"\t"<<n_line<<"\t";
	    }
	    //getchar();
	}
	if (doc.size() > 0) {
	    cout << doc.size() << endl;
	    data->add(doc);
            data->m_docLength.push_back (docLength);
	    doc.clear();
	    data->m_numDocs++;
	}
    }
    
    //data->m_numDocs = data->m_docs.size(); //update the numDocs
    //data->m_docLength = docLength; //update docLength
    printf("sum of all doc lengths: %d \n",actualDocLength);
    printf("> loaded %s (%d docs and %d terms)\n", datafile, data->m_numDocs, data->m_corpusLength);

    double bestll = -numeric_limits<double>::max();
    for (int run = 0; run < multiRunCount; run++) {
        printf("######### Learning ########\n");
        printf("> Random seed is: %d\n", seed);
        srand(seed);

        double diff = 0.0;
        double ll;
        
        TAModel *train_model;
        
        if (trainFlag) {
            if (initype == 1){
                // random initialization
                train_model = new TAModel(data, n_Aspects, nTz, trainFlag, &opt);
                train_model->setTrPriors();
                if(align_topics) {
                    train_model->align_topics(); }
                //cout<<"out of initialization :"<<train_model->m_trainFlag<<endl;
            }
            if (initype == 2){
                // in case of initializing from plsa
                train_model = new TAModel(data, initfile, n_Aspects, nTz, trainFlag, &opt);
                train_model->setTrPriors();
                if(align_topics) {
                    train_model->align_topics(); }
            }
            if (initype == 3){
                // in case of initializing from plsa
                train_model = TAModel::createWithPtrwz(data, initscalePtrwz, initfilePtrwz, n_Aspects, nTz, trainFlag, &opt);
                train_model->setTrPriors();
                if(align_topics) {
                    train_model->align_topics(); }
            }
        } else {
            //it is inference parameter order is changed just to make the difference
            train_model = new TAModel(data, trainFlag, n_Aspects, nTz, modelfile, &opt);
            //train_model->setTrPriors();
            if(opt.prune_topics) train_model->prune_topics(n_Aspects);
            cout<<"m_trainFlag " << train_model->m_trainFlag << "\n";
        }

       
        train_model->compute_pwtad();
        
        double old_ll = train_model->LogLikelihood(data);

        if (multiRunCount>1) printf("You can stop plsm with \"touch %s.stop\"\n", datafile);
        for (int iter = 0; iter < n_iter; iter++) {
            
            train_model->EM(data);
            train_model->normalize();
            
            
            if(align_topics && trainFlag) {
                train_model->align_topics();
            }
            
            time (&start);		
            train_model->compute_pwtad();
            ll = train_model->LogLikelihood(data);
            time (&end);
            diff = difftime (end,start);
            //printf ("It took %.2lf seconds to complete calc likelihoods.\n", diff );
            
            diff = (ll - old_ll) / (abs(ll) + EPS);
            printf("Run %d/%d: %f %f ", run+1, multiRunCount, ll, diff);
            if ((iter>5 && (abs(diff) < end_accuracy)) || (iter > n_iter)) { // 30 iter minimum
                printf("Early stopping...\n");
                break;
            }
            old_ll = ll;
            printf("\n");
        }
        if (ll > bestll) {
            train_model->LogLikelihood_doc(data);
            train_model->saveModel(datafile, 1, 1, 1, 1);
            printf(" > model saved to file\n");
            bestll = ll;
            if (writePerInstantLL && strlen(writePerInstantLL)>0) {
                train_model->writePerInstantLogLikelihood(writePerInstantLL, "%f\n", data);
            }
        } else {
            printf(" > dropped model as it is not better than previous ones\n");
        }

        delete train_model;
        seed = seed +100;
        if (userTouchedTheStopFileForSince(datafile, startupDate)) {
            printf("Stopping plsm because user asked for it by touching %s.stop\n", datafile);
            break;
        }
    }
    delete data;

}


