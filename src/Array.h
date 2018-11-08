/*
 * Copyright (c) 2013 - Idiap Research Institute - http://www.idiap.ch/
 * 
 * See file COPYING for the licence associated to this software
 * 
 * This file is a part of the TopicModels4Video software version: academic-1.0
*/
#include <limits>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>

#ifndef ARRAY
#define ARRAY


// review by remi: this might cause problem
#define EPS std::numeric_limits<double>::epsilon()

class Array {

private:


public:
    //int m_width;
    //int m_height;
   // int m_depth;
    //double ***mat;
    //std::vector<std::vector<std::vector<double> > > mat;


    //Array();

    //Array(int r, int c, int d);
    
    //Array(const Array&);
    //Array& operator=(const Array&);

    //Array(int r,int c);
    //void set(int r,int c, double val);
    //double get(int r,int c);

    //void set(int r, int c, int d, double val);
    //double get(int r, int c, int d);
    // ^ not implemented, actually existing code uses this->mat[i][j][k] directly
    void file_initialize(const char *filename, double *input, int ncols, int nrows, int ndepth);
    void file_initialize_3d1(const char *filename, double *input, int ncols, int nrows, int ndepth);
    void file_initialize_3d(const char *filename, double *input, int ncols, int nrows, int ndepth);
    void file_initialize_2d(const char *filename, double *input, int ncols, int nrows);
    void multiply(double *input, double val, int width, int height, int depth);
    void add(double *input, double val, int width, int height, int depth);
    void random_initialize(double* input,int width,int height,int depth, int order);
    double sum(double* input,int width,int height,int depth);

    void val_initialize(double* input,int width,int height,int depth,double val);
    void normalize(double* input,int width,int height,int depth);
    void normalize_y(double* input,int width,int height,int depth, int dim);
    void normalize2(double* input,int width,int height,int depth);
    void saveArray(const char* datafile, const char *ext, int n_Aspects, int n_Tz, double* output, int width,int height,int depth, int order );
    void readArray(const char* datafile, const char *ext, int n_Aspects, int n_Tz, double* output, int width,int height,int depth ) ;
    void tr_prior(double* input,int height); 
    void print(double* input,int width,int height,int depth);
    void compute_pd(double *h_Pd, double *h_Doc, int nDocs, std::vector<int>& m_nTd , int nVoc);
    void compute_pz(double *Pz, double *train_pzd, double *train_pd,int TrainDocs , int nZ);
    void compute_initial_ptrwz_pwz(double* h_Ptrwz,double* h_Pwz,double*p_w_tr_given_z,double* p_w_given_tr_z,double *trp, int nZ,int nVoc,int nTr);
    void setTrPriors(double *log_trp, int nTz);
    //void align_topics(double *ptrwz,double* ptszd,double *pwz,double *log_trp,double *trp, int nTz,int nAspects, int nVoc);
    void align_topics(double *ptrwz,double* ptszd,double *pwz,double *log_trp,double *trp, std::vector<int>& m_nTd, std::vector<int>& m_nTs, int nTz,int nAspects, int nVoc, int nDocs,int nTs);
    //void copy(Array *);
    void copy(double *output,double *input,int x, int y, int z );
 	
    
   
	
   

    ~Array();

};


#endif

