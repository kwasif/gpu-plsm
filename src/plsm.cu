 /* Copyright (c) 2016 - Advanced Digital Sciences Centre - http://www.adsc.com.sg
 * 
 * See file COPYING for the licence associated to this software
 * 
 * 
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <stdio.h>
using namespace std; 

#include <float.h>




// plsm library
#include "CmdLine.h"
#include "Array.h"
#include <sys/stat.h>
#include <limits>


//#include "cuPrintf.cu"

// includes CUDA
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// includes, project
#include <helper_cuda.h>
#include <helper_functions.h> // helper functions for SDK examples

typedef std::numeric_limits <double> dbl;

typedef std::pair<int, int> WordCount;
typedef std::vector<WordCount> Frame;
typedef std::vector<Frame> Document;


enum FitMode {DEFAULT =0,NORMAL,NORMAL_UNIFORM};

struct Options{
    bool debug_flag;
    bool tr_prior;
    bool z_prior;
    bool prune_topics;
    float z_thresh;
    float lambdaTsSparsity;
    float tr_wt;
    int nvoc;
    FitMode fitMode;
    float fitParam1;
    float fitParam2;
    bool usePriorPtrwz;
    bool align_topics;

    Options() : fitMode(DEFAULT), fitParam1(0), fitParam2(0) {
    }
};





// Global Arrays
double *h_Priorwz;
double *h_Priortrwz;
std::vector<int> m_nTs;
std::vector<int> m_nTd;

Array p_obj;

const int Max_iterations =50;
const double epsilon = 0.0001;
const int tsSparsityFlag=1;
const int totlhood=1;
const double lambda=0.25;
////////////////////////////////////////////////////////////////////////////////
// declaration, forward
double plsm_func( double *h_Doc,double *h_Pd,double* h_Ptszd, double* h_Ptrwz, double *h_Pzd, double *h_Pwz, double *h_Pwtad, double *h_Pz,double *h_trp,double *log_trp, int nDocs, int max_nTd, int nTz, int nVoc,int nZ, int nTs,int totWords, float end_accuracy,bool trainflag, struct Options &opt);

//void random_initialize(double* input,int width,int height,int depth);



__global__ void calc_pwtad (double *g_Pwtad,double *g_Ptszd,double *g_Pzd,double *g_Pwz,double *g_Ptrwz,double *g_Pd, int g_nDocs,int g_nZ, int g_nTs,int g_nVoc, int g_nTr)
{

// nTr,Nvoc
//__shared__ double tmp_pwtad[100];
//y is nz  x is nts
const int block_width_x=128;
__shared__ double tmp_pd;
__shared__ double tmp_pwtad[block_width_x];
__shared__ double tmp_pzd;
__shared__ double tmp_pwz;
__shared__ double tmp_ptrwz;
int pwtad_index;
int Td_index;
__shared__ double tmp_Ptszd[block_width_x];
__shared__ int w;
tmp_pd=g_Pd[blockIdx.x];
w=blockIdx.y;
__syncthreads();

int nTs_counter=ceil((float)g_nTs/(float)block_width_x);

int index_Ptszd;
for (int tr=0;tr<g_nTr;tr++)
{   

 // for (int w=0;w<g_nVoc;w++)
	{
           for (int ts=0;ts<nTs_counter;ts++)
              {

		tmp_pwtad[threadIdx.x]=0;
		__syncthreads();
               for (int z=0;z<g_nZ;z++)
                   {
                    
                   tmp_pwz=g_Pwz[z*g_nVoc+w];
		   tmp_ptrwz=g_Ptrwz[(z*g_nVoc+w)*g_nTr+tr];
		   tmp_Ptszd[threadIdx.x]=0;
                    
		    tmp_pzd=g_Pzd[blockIdx.x*g_nZ+z];
                   
                    
                    __syncthreads();
		index_Ptszd=(blockIdx.x*g_nZ+z)*g_nTs+block_width_x*ts+threadIdx.x;	 
			
                        if ((block_width_x*ts+threadIdx.x)<g_nTs )
                           {
			
                        	tmp_Ptszd[threadIdx.x]= g_Ptszd[index_Ptszd];
                            }
			__syncthreads();

                     

                          
                        
			
			if 	 ((block_width_x*ts+threadIdx.x)< g_nTs)
			

			{
			
			tmp_pwtad[threadIdx.x]+=tmp_pd*tmp_pzd*tmp_Ptszd[threadIdx.x]*tmp_pwz*tmp_ptrwz;
			
			}
                        __syncthreads();
                        
                       

                    } // nz loop  

			Td_index=tr+block_width_x*ts+threadIdx.x;
			
                          if(Td_index< (g_nTs+g_nTr-1))
                	     {
				
			
			        pwtad_index=(blockIdx.x*(g_nTs+g_nTr-1)+Td_index)*g_nVoc+w;
                             
                               g_Pwtad[pwtad_index]+=tmp_pwtad[threadIdx.x];
				
			    }	
                             __syncthreads();
						
			
				
                 		
                	     

                        
                   
               } // nts counter loop
 

	}// tr loop
	}// w loop



}


//1024(nvoc) block size and number of blocks is nDocs 
__global__ void calc_sum_doc(double *g_Docsum, double *g_Doc,int g_nDocs, int g_nTs,int g_nVoc, int g_nTr,int totWords)
{
const int block_width=1024;
__shared__ double shared_block[block_width];

int iterations_nvoc = ceil((float)g_nVoc/(float)block_width);
int doc_index;
double sum=0;

int nTd=g_nTs+g_nTr-1;
// (g_nTs+g_nTr)
for (unsigned int ntd_counter=0;ntd_counter<nTd;ntd_counter++)
    {
	 for (unsigned int nvoc_counter=0;nvoc_counter<iterations_nvoc;nvoc_counter++)
	     {
                doc_index=((blockIdx.x*nTd+ntd_counter)*g_nVoc+nvoc_counter*block_width+threadIdx.x);
                 // initialize shared memory to zero  
	        shared_block[threadIdx.x]=0;
               
                
                __syncthreads();
              
		if ( nvoc_counter*block_width+threadIdx.x< g_nVoc)
                      { 
                         shared_block[threadIdx.x]=g_Doc[doc_index];
                         //shared_block[threadIdx.x]=logf(2);
                        
	         	
                      
		       
                      
                 
                     __syncthreads();
                       }
               
	        for (unsigned int s=block_width/2;s>=1;s>>=1)
                    { 
	            if (threadIdx.x<s)
			{
	            	  shared_block[threadIdx.x]+=shared_block[threadIdx.x+s];
               	    	 
			}
                    __syncthreads(); 
                    }

                
                if(threadIdx.x==0)
                {
                sum=sum+shared_block[0];
                
               
                }
                __syncthreads();
	     }

      }

//last iteration do it here to avoid padding 

if(threadIdx.x==0)
                {
g_Docsum[blockIdx.x]=sum/totWords;




}
__syncthreads();
}




//1024(nvoc) block size and number of blocks is nDocs 
__global__ void calc_logmatsum(double *g_Logmat_sum,double *g_Docsum, double *g_Doc,double *g_Pwtad,int g_nDocs, int g_nTs,int g_nVoc, int g_nTr)
{
const int block_width=1024;
//__shared__ double shared_block[block_width];
__shared__ double shared_logmat_block[block_width];
int iterations_nvoc = ceil((float)g_nVoc/(float)block_width);
int doc_index;
//double sum=0;
double logmat_sum=0;
int nTd=g_nTs+g_nTr-1;
// (g_nTs+g_nTr)
for (unsigned int ntd_counter=0;ntd_counter<nTd;ntd_counter++)
    {
	 for (unsigned int nvoc_counter=0;nvoc_counter<iterations_nvoc;nvoc_counter++)
	     {
                doc_index=((blockIdx.x*nTd+ntd_counter)*g_nVoc+nvoc_counter*block_width+threadIdx.x);
                 // initialize shared memory to zero  
	       // shared_block[threadIdx.x]=0;
                shared_logmat_block[threadIdx.x]=0;
                
                __syncthreads();
                
		if ( nvoc_counter*block_width+threadIdx.x< g_nVoc)
                      { 
                        
	         	
                      
		       
                       shared_logmat_block[threadIdx.x]=g_Doc[doc_index]*log(g_Pwtad[doc_index]+DBL_EPSILON);
                   
                     __syncthreads();
                       }
            
	        for (unsigned int s=block_width/2;s>=1;s>>=1)
                    { 
	            if (threadIdx.x<s)
			{
	            	
               	    	  shared_logmat_block[threadIdx.x]+=shared_logmat_block[threadIdx.x+s];
			}
                    __syncthreads(); 
                    }

                
                if(threadIdx.x==0)
                {
                
                logmat_sum=logmat_sum+shared_logmat_block[0];
              
                }
                __syncthreads();
	     }

      }

//last iteration do it here to avoid padding 

if(threadIdx.x==0)
                {
//g_Docsum[blockIdx.x]=sum;
g_Logmat_sum[blockIdx.x]=logmat_sum;//g_Docsum[blockIdx.x];


}
__syncthreads();
}




__global__ void calc_observ_lhood(double L,double *g_Logmat_sum,int g_nDocs)
{
//change this based on number of documents
for (int i=0;i<g_nDocs;i++)
L=L+g_Logmat_sum[i];

printf("Lilkelihood is :: %0.8f \n",L);
}

__global__ void calc_lhood(double E,double *g_Pzd,double *g_Ptszd,double *g_Pwz, double *g_Ptrwz,double *g_Pwtad,double *g_Doc,int g_nDocs, int g_nTs,int g_nVoc, int g_nTr,int g_nZ)
{
int Td=g_nTs+g_nTr-1;
double temp=0;
double temp1=0;
	for (int d=0;d<g_nDocs;d++)
	     {
		for (int w=0;w<g_nVoc;w++)
	    	    {   

  	     	      for (int tr=0;tr<g_nTr;tr++)
		          {
           	   
		           for (int z=0;z<g_nZ;z++)
              		       {
               	          
			        for (int ts=0;ts<g_nTs;ts++) 
				    {
                                     temp1=g_Pzd[d*g_nZ+z]*g_Ptszd[(d*g_nZ+z)*g_nTs+ts]*g_Pwz[z*g_nVoc+w]*g_Ptrwz[(z*g_nVoc+w)*g_nTr+tr];
                                     temp=temp1*log(temp1+DBL_EPSILON/g_Pwtad[(d*Td+ts+tr-((ts+tr)/(Td-1)))*g_nVoc+w]+DBL_EPSILON);
				      E=E+g_Doc[(d*Td+ts+tr-((ts+tr)/(Td-1)))*g_nVoc+w]	*temp;



				     }
			       }
		          }
		    }
		}
	

printf("E is ::%0.8f \n",E);
}







 //dim3  grid_tszd(nZ, nDocs,1);

 //   dim3  threads_tszd(32, 16,1);
__global__ void calc_tszd(double *g_zd, double *g_tszd,double *g_Ptszd,double *g_Pzd,double *g_Pwz,double *g_Ptrwz,double *g_Pwtad,double *g_Doc,double *pd,double *pz, int g_nDocs,int g_nZ, int g_nTs,int g_nVoc, int g_nTr,float lambdaTsSparsity, bool z_prior , bool trainflag)
{

const int block_width_x=16;
const int block_height_y=16;
__shared__ double partial_tszd[block_height_y*block_width_x];
double tmp_tszd=0;
double tmp_zd=0;
int gdoc_index;
int iterations_nvoc = ceil((float)g_nVoc/(float)block_width_x);
int iterations_ntr = ceil((float)g_nTr/(float)block_height_y);
__shared__ double decrease;


        
	for (int ts=0;ts<g_nTs;ts++) 
            {
                //tmp_sum=0;
                  tmp_tszd=0;
		    for (unsigned int nvoc_counter=0;nvoc_counter<iterations_nvoc;nvoc_counter++)
		     {
			 for (unsigned int ntr_counter=0;ntr_counter<iterations_ntr;ntr_counter++)
		     		{
			partial_tszd[threadIdx.y*block_width_x+threadIdx.x]=0;
			__syncthreads();	
			if ((ntr_counter*blockDim.y+threadIdx.y)<g_nTr && (nvoc_counter*blockDim.x+threadIdx.x)<g_nVoc )
			{	gdoc_index=(blockIdx.y*(g_nTs+g_nTr-1)+ts+ntr_counter*blockDim.y+threadIdx.y-(ts+ntr_counter*blockDim.y+threadIdx.y)/(g_nTs+g_nTr-1))*g_nVoc+nvoc_counter*blockDim.x+threadIdx.x;
                           
			partial_tszd[threadIdx.y*block_width_x+threadIdx.x]=(pd[blockIdx.y]*g_Doc[gdoc_index]*g_Pzd[(blockIdx.y)*g_nZ+blockIdx.x]*g_Ptszd[(blockIdx.y*g_nZ+blockIdx.x)*g_nTs+ts]*g_Pwz[(blockIdx.x)*g_nVoc+nvoc_counter*blockDim.x+threadIdx.x]*g_Ptrwz[(blockIdx.x*g_nVoc+nvoc_counter*blockDim.x+threadIdx.x)*g_nTr+ntr_counter*blockDim.y+threadIdx.y])/(g_Pwtad[gdoc_index]+DBL_EPSILON);
	                 
		
                         } 
			__syncthreads();
                        for (unsigned int s=block_width_x*block_height_y/2;s>=1;s>>=1)
       				{//g_tszd[s]=s;
       				 
				  if (threadIdx.y*block_width_x+threadIdx.x<s)
	    			  {
	     			  
             			    partial_tszd[threadIdx.y*block_width_x+threadIdx.x]+=partial_tszd[threadIdx.y*block_width_x+threadIdx.x+s];
            			    
	    			   }
				__syncthreads();
       
        			}
                       if (threadIdx.x==0 && threadIdx.y==0)
                       {

                         
			tmp_tszd+=partial_tszd[0];
			
			}      
			__syncthreads();
		     } //ntr counter
		      }	//nvoc counter
		   
		if (threadIdx.x==0 && threadIdx.y==0)	
		    {	
			
			if (lambdaTsSparsity!=0){
			decrease = lambdaTsSparsity*pd[blockIdx.y]/(g_nZ*g_nTs);
			tmp_tszd-=decrease;
			tmp_tszd=max(DBL_EPSILON,tmp_tszd);
			}
			
			g_tszd[(blockIdx.y*g_nZ+blockIdx.x)*g_nTs+ts]=tmp_tszd;
			tmp_zd+=tmp_tszd; 
			
			
		     }

	        __syncthreads();
	    } // nts loop
if (threadIdx.x==0 && threadIdx.y==0)	
{
if (z_prior && !trainflag)
tmp_zd+= pd[blockIdx.y]/(g_nZ*pz[blockIdx.x]);
	
g_zd[blockIdx.y*g_nZ+blockIdx.x]=tmp_zd;
}

}



//thread 32 normalize 2d matrix
//grid (1,g_height,1)
__global__ void
normalize2d_gpu(double *g_input, int g_width,int g_height)
{

// 2d matrix
const int block_width_x=32;
int x_iter = ceil((float(g_width)/(float)block_width_x));
__shared__ double partial_sum[block_width_x];

__shared__ double tmp_sum;
tmp_sum=0;
__syncthreads();
		for (unsigned int x_counter=0;x_counter<x_iter;x_counter++)
	     	    {   partial_sum[threadIdx.x]=0;
			__syncthreads();
			if (x_counter*block_width_x+threadIdx.x < g_width)
              		   {
				partial_sum[threadIdx.x]=g_input[blockIdx.y*g_width+x_counter*block_width_x+threadIdx.x];
			        __syncthreads();
			   }
			 for (unsigned int s=block_width_x/2;s>=1;s>>=1)
       			     {
       				 
				 if (threadIdx.x<s)
	    			    {
	     			   
             			    partial_sum[threadIdx.x]+=partial_sum[threadIdx.x+s];
            			      __syncthreads();
	    			    }
				
       
        		     }
			 if(threadIdx.x==0)
			   {
				
                                 tmp_sum+=partial_sum[0];
				
		           	__syncthreads();
			   }
	            }
		 if (threadIdx.x==0 )
                     {
			
			tmp_sum=tmp_sum+(tmp_sum==0);
			
			
		   
		      }
		__syncthreads();
         	for (unsigned int x_counter1=0;x_counter1<x_iter;x_counter1++)
	    	    {  
		      if (x_counter1*block_width_x+threadIdx.x<g_width)
			 {		      
			
			   g_input[blockIdx.y*g_width+x_counter1*block_width_x+threadIdx.x]=g_input[blockIdx.y*g_width+x_counter1*block_width_x+threadIdx.x]/tmp_sum;

		          
			  }
			__syncthreads();
                     }
		      
		
}

//thread 64 normalize 3d matrix

//grid (g_Width,g_height)
// tszd  along dim ts
//dim3  grid_norm3d(nZ,nDocs,1);
//dim3  threads_norm3d(1,1, 64);
__global__ void
normalize3d_gpu(double *g_input, int g_width,int g_height, int g_depth)
{

// 2d matrix
const int block_depth_x=64;
int x_iter = ceil((float(g_depth)/(float)block_depth_x));
__shared__ double partial_sum[block_depth_x];
__shared__ double row_sum;
__shared__ double tmp_sum;
row_sum=0;
//double depth_sum=0;
tmp_sum=0;
__syncthreads();
	
		for (unsigned int x_counter=0;x_counter<x_iter;x_counter++)
	     	    {   partial_sum[threadIdx.x]=0;
			__syncthreads();
			if (x_counter*block_depth_x+threadIdx.x < g_depth)
              		   {
				partial_sum[threadIdx.x]=g_input[( blockIdx.y*g_width+blockIdx.x)*g_depth+x_counter*block_depth_x+threadIdx.x];
			       
			   }
                         __syncthreads();
			 for (unsigned int s=block_depth_x/2;s>=1;s>>=1)
       			     {
       				 
				 if (threadIdx.x<s)
	    			    {
	     			   
             			    partial_sum[threadIdx.x]+=partial_sum[threadIdx.x+s];
            			    
	    			    }
				  __syncthreads();
       
        		     }
			 if(threadIdx.x==0)
			   {
				//g_sum[blockIdx.y]+=partial_sum[0];
                                 tmp_sum+=partial_sum[0];
				
		           	
			   }
			  __syncthreads();
	            }
		 if (threadIdx.x==0 )
                     {
			
			row_sum=tmp_sum+(tmp_sum==0);
		
		      }
                 //depth_sum=tmp_sum+(tmp_sum==0);
		__syncthreads();
		//row_sum=depth_sum
         	for (unsigned int x_counter1=0;x_counter1<x_iter;x_counter1++)
	    	    {  
		      if (x_counter1*block_depth_x+threadIdx.x<g_depth)
			 {		      
			 
 			g_input[( blockIdx.y*g_width+blockIdx.x)*g_depth+x_counter1*block_depth_x+threadIdx.x]=g_input[( blockIdx.y*g_width+blockIdx.x)*g_depth+x_counter1*block_depth_x+threadIdx.x]/row_sum;

		          __syncthreads();
			  }
                     }
		      
		
}

// grid nz,ndocs,1
//threads 32  // loop over nts
//tszd_sparsityflag<<<grid_tszd,threads_tszd>>>(d_tszd,d_Docsum,lambda,nDocs,nZ,nTs,nVoc);
__global__ void tszd_sparsityflag(double *g_tszd,double *g_Docsum,const double lambda,int nDocs, int nZ,int nTs,int nVoc)
{
__shared__ double thresh;
const int block_depth_x=32;
__shared__ double shared_sparsity[block_depth_x];
thresh=(lambda*g_Docsum[blockIdx.y])/(nZ*nTs);
//__syncthreads();


int x_iter = ceil((float(nTs)/(float)block_depth_x));
//__shared__ double partial_ts[block_depth_x];

	
		for (unsigned int x_counter=0;x_counter<x_iter;x_counter++)
		    {  // partial_ts[threadIdx.x]=0;
			shared_sparsity[threadIdx.x]=0;
			__syncthreads();
			
			if ((x_counter*block_depth_x+threadIdx.x) < nTs)                      
			   {			   
				
                        shared_sparsity[threadIdx.x]=max(DBL_EPSILON,g_tszd[(blockIdx.y*nZ+blockIdx.x)*nTs+x_counter*block_depth_x+threadIdx.x]-thresh);
                      
			
			}
			__syncthreads();
			if ((x_counter*block_depth_x+threadIdx.x) < nTs)      
			g_tszd[(blockIdx.y*nZ+blockIdx.x)*nTs+x_counter*block_depth_x+threadIdx.x]=shared_sparsity[threadIdx.x]	;
			__syncthreads();
			
		     }// nTs loop
}



//dim3  grid(nVoc,nZ,1);
    //dim3  threads(32,16)(nts, nTr,1);
//calc_tszd1<<<grid_tszd,threads_tszd>>>(d_zd, d_tszd,d_Ptszd,d_Pzd,d_Pwz,d_Ptrwz,d_Pwtad,d_Doc,nDocs,nZ, nTs,nVoc, nTr);
__global__ void calc_wz(double *g_wz,double *g_Ptszd,double *g_Pzd,double *g_Pwz,double *g_Ptrwz,double *g_Pwtad,double *g_Doc,int g_nDocs,int g_nZ, int g_nTs,int g_nVoc, int g_nTr)
{

const int block_width_x=32;
const int block_height_y=16;

__shared__ double partial_wz[block_width_x*block_height_y];

__shared__ double partial_ptszd[block_width_x];
__shared__ double partial_ptrwz[block_height_y];
__shared__ double partial_zd;
__shared__ double partial_pwz;

int Td;
 double tmp_sum;

 double tmp_sum_outerloop;
tmp_sum_outerloop=0;
int gdoc_index;
int iterations_nTs = ceil((float)g_nTs/(float)block_width_x);
int iterations_nTr = ceil((float)g_nTr/(float)block_height_y);
partial_pwz=g_Pwz[blockIdx.y*g_nVoc+blockIdx.x];



__syncthreads();



 
        for (int d=0;d<g_nDocs;d++)
	{    	         tmp_sum=0;
			partial_zd= g_Pzd[d*g_nZ+blockIdx.y];
	               __syncthreads();
                	

		    for (unsigned int nts_counter=0;nts_counter<iterations_nTs;nts_counter++)
		     {
			


			 for (unsigned int ntr_counter=0;ntr_counter<iterations_nTr;ntr_counter++)
		         {
                                
				partial_ptszd[threadIdx.x]=0;
			        partial_ptrwz[threadIdx.y]=0;
			        partial_wz[threadIdx.y*block_width_x+threadIdx.x]=0;
			        __syncthreads();


   
				if (ntr_counter*blockDim.y+threadIdx.y< g_nTr)
				partial_ptrwz[threadIdx.y]=g_Ptrwz[(blockIdx.y*g_nVoc+blockIdx.x)*g_nTr+ntr_counter*blockDim.y+threadIdx.y];
				__syncthreads();	
			if ((ntr_counter*blockDim.y+threadIdx.y)<g_nTr && (nts_counter*blockDim.x+threadIdx.x)<g_nTs)
				
			{	
				Td=(nts_counter*blockDim.x+threadIdx.x)+ntr_counter*blockDim.y+threadIdx.y;
				gdoc_index=(d*(g_nTs+g_nTr-1)+Td-(Td/(g_nTs+g_nTr-1)))*g_nVoc+ blockIdx.x;
                                 
                                partial_ptszd[threadIdx.x]=g_Ptszd[(d*g_nZ+blockIdx.y)*g_nTs+nts_counter*block_width_x+threadIdx.x];
				//__syncthreads();
                         
                          
			partial_wz[threadIdx.y*block_width_x+threadIdx.x]=(g_Doc[gdoc_index]*partial_zd*partial_ptszd[threadIdx.x]*partial_pwz*partial_ptrwz[threadIdx.y])/(g_Pwtad[gdoc_index]+DBL_MIN);
			
                         }
			__syncthreads();
			
			 
                    
                        for (unsigned int s=block_width_x*block_height_y/2;s>=1;s>>=1)
       				{//g_tszd[s]=s;
       				 
				  if (threadIdx.y*block_width_x+threadIdx.x<s)
	    			  {
	     			   
             			    partial_wz[threadIdx.y*block_width_x+threadIdx.x]+=partial_wz[threadIdx.y*block_width_x+threadIdx.x+s];
            			    
	    			   }
				__syncthreads();
       
        			}
                      
		      	//__syncthreads();
			 if (threadIdx.x==0 && threadIdx.y==0)
                       {
			tmp_sum += partial_wz[0];
			
                       }
		       __syncthreads();
			
		} //ntr counter
		} // nts counter loop
       if (threadIdx.x==0 && threadIdx.y==0)
	{	tmp_sum_outerloop+=tmp_sum;
		
	}
	__syncthreads();

	

	} //nDocs for loop
if (threadIdx.x==0 && threadIdx.y==0)	
g_wz[blockIdx.y*g_nVoc+blockIdx.x]=tmp_sum_outerloop;

__syncthreads();
}




// loop over ndocs
//dim3  grid_trwz(nVoc,nZ,1);
//dim3  threads_trwz(64, 1,1); loop over nts
__global__ void calc_trwz(double *g_trwz,double *g_wz,double *g_Ptszd,double *g_Pzd,double *g_Pwz,double *g_Ptrwz,double *g_Pwtad,double *g_Doc,double *pd,double *g_trp,double *priortrwz,double *priorwz,int g_nDocs,int g_nZ, int g_nTs,int g_nVoc, int g_nTr,int totWords,float tr_wt,bool usePriorPtrwz,bool tr_prior)
{

const int block_width_x=64;
int gdoc_index;
int iterations_nTs = ceil((float)g_nTs/(float)block_width_x);

 __shared__ double partial_trwz[block_width_x];
int Td;
 __shared__ double tmp_sum;

 __shared__ double tmp_sum_outerloop;



__shared__ double partial_zd;
__shared__ double partial_pwz;
__shared__ double partial_ptrwz;

__shared__ double more;
__shared__ double tmp_pd;
 __shared__ double tmp_sum_wz;
__shared__ double  word_sum_count[block_width_x];
double sum_count;

tmp_sum_wz=0;
tmp_sum=0;
tmp_sum_outerloop=0;
sum_count=0;
partial_pwz=g_Pwz[blockIdx.y*g_nVoc+blockIdx.x];




__syncthreads();


       for(int tr=0;tr<g_nTr;tr++)
	{
         partial_ptrwz=g_Ptrwz[(blockIdx.y*g_nVoc+blockIdx.x)*g_nTr+tr];
	 sum_count=0;
 	 __syncthreads();
        for (int d=0;d<g_nDocs;d++)
	{    	         partial_zd=g_Pzd[d*g_nZ+blockIdx.y];
			tmp_sum=0;
                        tmp_pd=pd[d];
			__syncthreads();
                	

		    for (unsigned int nts_counter=0;nts_counter<iterations_nTs;nts_counter++)
		     {
			partial_trwz[threadIdx.x]=0;
			word_sum_count[threadIdx.x]=0;
	               __syncthreads();
				
			if (nts_counter*blockDim.x+threadIdx.x <g_nTs)
				
			{	//gdoc_index=(blockIdx.y*(g_nTs+g_nTr-1)+ts+threadIdx.y-(ts+threadIdx.y)/(g_nTs+g_nTr-2))*g_nVoc+nvoc_counter*blockDim.x+threadIdx.x+ *g blockIdx.x;
				Td=(nts_counter*block_width_x+threadIdx.x)+tr;

				gdoc_index=(d*(g_nTs+g_nTr-1)+Td-(Td/(g_nTs+g_nTr-1)))*g_nVoc+ blockIdx.x;
                                 
                                word_sum_count[threadIdx.x]=g_Doc[gdoc_index];
				partial_trwz[threadIdx.x]=(g_Doc[gdoc_index]*tmp_pd*partial_zd*g_Ptszd[(d*g_nZ+blockIdx.y)*g_nTs+nts_counter*block_width_x+threadIdx.x]*partial_pwz*partial_ptrwz)/(g_Pwtad[gdoc_index]+DBL_MIN);
                             
			
                         }
			__syncthreads();
			
                        for (unsigned int s=block_width_x/2;s>=1;s>>=1)
       				{//g_tszd[s]=s;
       				 
				  if (threadIdx.x<s)
	    			  {
	     			   partial_trwz[threadIdx.x]+=partial_trwz[threadIdx.x+s];
            			   word_sum_count[threadIdx.x]+=word_sum_count[threadIdx.x+s];
	    			   }
				__syncthreads();
       
        			}
                      
		      	//__syncthreads();
			 if (threadIdx.x==0 )
                       {
			
			tmp_sum += partial_trwz[0];
			sum_count+=word_sum_count[0];
		 
                        	  }
		       __syncthreads();
			//g_tszd[(blockIdx.y*g_nZ+blockIdx.x)*g_nTs+ts]+=partial_tszd[0];
                       	//printf("g_tszd value for index %d is %0.8f\n",(blockIdx.y*g_nZ+blockIdx.x)*g_nTs+ts,partial_tszd[0]);
			    
			//__syncthreads();
		     } // nts counter loop
       		if (threadIdx.x==0 )
		{	tmp_sum_outerloop+=tmp_sum;
			// if ( blockIdx.x==1 && blockIdx.y==0 && blockIdx.z==0 )
			//printf(" d is %d and tmp_sum is  %0.8f and  tmp_sum_outerloop  is %0.8f \n",d,tmp_sum,tmp_sum_outerloop);
		}
		__syncthreads();

	

	} //nDocs for loop
if (threadIdx.x==0 )	
{

 	if (usePriorPtrwz)
	   {   
		tmp_sum_outerloop+=priortrwz[(blockIdx.y*g_nVoc+blockIdx.x)*g_nTr+tr]*priorwz[blockIdx.y*g_nVoc+blockIdx.x]*totWords;
	    }
	else if (tr_prior){
		
		tmp_sum_outerloop+=tr_wt*sum_count*g_trp[tr]/(g_nVoc*g_nZ);// fill it
	    }
	g_trwz[(blockIdx.y*g_nVoc+blockIdx.x)*g_nTr+tr]=tmp_sum_outerloop;
	tmp_sum_wz+=tmp_sum_outerloop;
	tmp_sum_outerloop=0;
	tmp_sum=0;
 
}
__syncthreads();
  }// ntr loop
if (threadIdx.x==0 )
{
	
g_wz[blockIdx.y*g_nVoc+blockIdx.x]=tmp_sum_wz;
}
__syncthreads();
}



//dim3  grid(nTr,nZ,1);
    //dim3  threads(16)(Nvoc, 1,1);
// loop over ndocs

__global__ void MakeJoint(double *g_Topics,double *g_Pwz,double *g_Ptrwz,int g_nZ,int g_nVoc, int g_nTr)
{

const int block_width_x=16;

int iterations_nvoc = ceil((float)g_nVoc/(float)block_width_x);




 		for (unsigned int nvoc_counter=0;nvoc_counter<iterations_nvoc;nvoc_counter++)
		     {
			
				
			if (nvoc_counter*blockDim.x+threadIdx.x < g_nVoc)
                       g_Topics[(blockIdx.y*g_nTr+blockIdx.x)*g_nVoc+nvoc_counter*blockDim.x+threadIdx.x]=g_Pwz[blockIdx.y*g_nVoc+nvoc_counter*blockDim.x+threadIdx.x]*g_Ptrwz[(blockIdx.y*g_nVoc+nvoc_counter*blockDim.x+threadIdx.x)*g_nTr+blockIdx.x];
			__syncthreads();
                     }
}


double* normalize_this(double *orig, int orig_x, int orig_y,int orig_z, int dim) 
{
 double *orig_norm = (double *) malloc(sizeof(double)*orig_x*orig_y*orig_z);
// one dimensional input array
if (dim==0 && (orig_x==1 || orig_y==1))
	{   
	    int max_dim=max(orig_x,orig_y);
	    
	    
            double orig_sum=0;
		for (int i=0; i<max_dim;i++)
		orig_sum+=orig[i];
                if (orig_sum==0)
                   orig_sum=1;
                for (int j=0; j<max_dim;j++)
		orig_norm[j]=orig[j]/orig_sum;  
	}
// 2/3 dimensional input array and normalizing along columns
else if (dim==1)
	{
	    unsigned int mem_size = sizeof(double) * orig_x*orig_z;

    	   // allocate host memory
           
            double *orig_sum = (double *) malloc(mem_size);
	
	    	for (int i=0;i<orig_x;i++)
			{
			for (int k=0;k<orig_z;k++)
	    		{				
				for (int j=0;j<orig_y;j++)
 				orig_sum[i*orig_z+k]+=orig[(i+orig_x*j)*orig_z+k];
		         
			if (orig_sum[i*orig_z+k]==0)
			   orig_sum[i*orig_z+k]+=1;
			}
             		}

            
	    	for (int i=0;i<orig_x;i++)
			{
				 for (int k=0;k<orig_z;k++)
				{
				for (int j=0;j<orig_y;j++)
 				{
				orig_norm[(i+orig_x*j)*orig_z+k]=orig[(i+orig_x*j)*orig_z+k]/orig_sum[i*orig_z+k];
		                
				}
		
				}
	    		}
         }  
return (orig_norm);
}

// function to copy vector to array

void vecToarray2d(vector<vector<double> > vals, double ** temp, int N,int M)
{
  //double** temp;
 // temp = new double[N][M];

 for(unsigned i=0; (i < N); i++)
 {
    for(unsigned j=0; (j < M); j++)
    {
        temp[i][j] = vals[i][j];
    }
 }
 //return temp;
}


void vecToarray3d(vector<vector<vector<double> >  >vals, double *** temp, int N, int M, int P)
{
  

 for(unsigned i=0; (i < N); i++)
 {
    for(unsigned j=0; (j < M); j++)
    {
	for(unsigned k=0; (k < P); k++)        
	temp[i][j][k] = vals[i][j][k];
    }
 }
 //return temp;
}


////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
	   //Options 
    struct Options opt;
    opt.debug_flag=false;  // useful to print stuff
    opt.tr_prior=false;
    opt.z_prior=false;
    opt.prune_topics=false;
    opt.z_thresh=0.95;
    opt.align_topics=false;

   const time_t startupDate = time(NULL);

    int seed = (int) time(NULL);
  

    int docLength;
    char *datafile;
    int n_Aspects, nTz;
    float end_accuracy=0.0001;
    int n_iter;
    char *modelfile = NULL;
    char *initfile = NULL;
    char *initfilePtrwz = NULL;
    float initscalePtrwz = 1.f;
    int initype = 1;
    float scaleOnRead = 1.f;
    int multiRunCount = 1;

    int trainFlag;
    
    bool d_overlap, read_several;
    //bool align_topics = false;
    time_t start,end;

    // 

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
    cmd.addBCmdOption("-alignz", &opt.align_topics, false, "Flag to use or not tr prior in learning[false]");

    cmd.addText("\n Fit modes:");
    cmd.addText("       0: table");
    cmd.addText("       1: gaussian (fit1: stddevToAdd, fit2: addUniformProportion)");
    cmd.addText("       2: gaussian/uniform (fit1: stddevToAdd, fit2: initialPi_u)");
    cmd.addICmdOption("-fitMode", (int*)&opt.fitMode, (int)opt.fitMode, "Fit Mode");
    cmd.addRCmdOption("-fit1", &opt.fitParam1, opt.fitParam1, "Fit parameter #1");
    cmd.addRCmdOption("-fit2", &opt.fitParam2, opt.fitParam2, "Fit parameter #2");



    cmd.read(argc, argv);
   
    printf("data file: %s \n",datafile);
    printf("Topics: %d \n",n_Aspects);
    printf("nTz: %d \n",nTz);
    printf("provided doc length: %d \n",docLength);
    printf("Paramter d_overlap: %d \n",d_overlap);
    printf("Paramter opt.lambdaTsSparsity: %0.8f \n",opt.lambdaTsSparsity);
    printf("Paramter opt.tr_prior: %d \n",opt.tr_prior);

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


    //Corpus *data = new Corpus(0);
    Document src_doc, doc;
    Frame frame;
    long m_totWords=0;
    int m_corpusLength=0;
    int m_numDocs=0;
    double *h_Doc;

   
    int max_nTd=0; 
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

			//data->m_totWords += count;
			m_totWords += count;

			//if (word >= data->m_corpusLength) {
			if (word >= m_corpusLength) {
			    //data->m_corpusLength = word + 1;
			    m_corpusLength = word + 1;	
			}
		    }
		}
		actualDocLength++;
                currentDocLength++;
                src_doc.push_back(frame);
		frame.clear();
	    }
	    
	   
           
            currentDocLength = 0;
	    src_doc.clear ();
	    //data->m_numDocs++;
	    m_numDocs++;
	}
     
    } else {
        cout << "OK" << endl;
	{			// Data loading (from a single file)
	    string line;
	    ifstream infile(datafile, ios_base::in);
            int curr_line=0;
	    while (getline(infile, line, '\n')) {curr_line++;
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

			//data->m_totWords += count;
			m_totWords += count;
			//if (word >= data->m_corpusLength) {
			
			if (word >= m_corpusLength) {
			    //data->m_corpusLength = word + 1;
			    m_corpusLength = word + 1;
			}
			 
		    }
		}
		actualDocLength++;
		src_doc.push_back(frame);
		frame.clear();
	    }
	}

cout<<"actualDocLength::"<<actualDocLength<<endl;
cout<<"m_totWords::"<<m_totWords<<endl;
cout<<"m_corpusLength::"<<m_corpusLength<<endl;
cout<<"m_numDocs::"<<m_numDocs<<endl;
cout<<"docLength::"<<docLength<<endl;
// ensure there is only one document.

//docLength=actualDocLength;
int total_docs=ceil((float)actualDocLength/docLength);
cout<<"total_docs is "<<total_docs<<endl;

if (d_overlap) 
//total_docs=ceil((float)(total_docs*(nTz-1)+actualDocLength)/docLength)+1;
total_docs=ceil((float)(actualDocLength)/(float)(docLength-nTz+1));
cout<<"total_docs now is "<<total_docs<<endl;
	dim3 dims_Doc(docLength,total_docs,m_corpusLength);
        unsigned int mem_size_Doc = dims_Doc.x * dims_Doc.y *dims_Doc.z* sizeof(double);

        h_Doc = (double *) malloc(mem_size_Doc);
        memset(h_Doc,0,mem_size_Doc);

          if (docLength < 0)
            docLength = actualDocLength;
	// split doc
	int f_count = 0;
	int n_line = 0;
	int overlap_ctr=0;
	//data->m_numDocs = 0;
	m_numDocs = 0;
	while (f_count < src_doc.size()) {
	    //cout<<data->m_numDocs<<"\n";
	    if ((n_line + 1) % docLength == 0) {
		doc.push_back(src_doc.at(f_count));
		for(std::vector<WordCount>::const_iterator it=src_doc.at(f_count).begin();it!=src_doc.at(f_count).end();++it )
		//h_Doc[(it->first*docLength+n_line)*total_docs+m_numDocs]=it->second;
		{h_Doc[(docLength*m_numDocs+n_line)*m_corpusLength+it->first]+=it->second;
		//cout<<"n_line for if loop is:: "<<n_line<<endl;
		}
		f_count++;
		n_line++;
		
		//data->add(doc);
		   
                 m_nTd.push_back (docLength);
		 if (max_nTd<docLength)
		     max_nTd=docLength;
                 m_nTs.push_back (docLength - nTz + 1);
    
                //data->m_docLength.push_back (docLength);
		doc.clear();
		//data->m_numDocs++;
		m_numDocs++;
		
		n_line = 0;
		if (d_overlap) {
		    f_count = f_count - nTz + 1;
		   overlap_ctr=0;
		    //cout<<"Document Overlap :"<<d_overlap<<endl;
		}

	    }			// push the frame into the document
	    else {
               // cout<<"source doc size at line:: "<<n_line<<" is:: "<<src_doc.at(f_count).size()<<endl;
		doc.push_back(src_doc.at(f_count));
		
		for(std::vector<WordCount>::const_iterator it=src_doc.at(f_count).begin();it!=src_doc.at(f_count).end();++it )
		{
		h_Doc[(docLength*m_numDocs+n_line)*m_corpusLength+it->first]+=it->second;
		
		}
		
		n_line++;
		f_count++;
		frame.clear();
		//cout<<f_count<<"\t"<<n_line<<"\t";
	    }// else loop
	    //getchar();
	}  //while loop
	if (doc.size() > 0) {
	    cout <<"doc size is:: "<<doc.size() <<"and src_docsize is "<<src_doc.size()<< endl;
		int i=0;
		

	   


	    m_nTd.push_back (docLength);
            m_nTs.push_back (docLength - nTz + 1);	
	    //m_docLength.push_back (docLength);
	    doc.clear();
	    //data->m_numDocs++;
	      m_numDocs++;
 
	}

cout<<"actualDocLength::"<<actualDocLength<<endl;
cout<<"m_totWords::"<<m_totWords<<endl;
cout<<"m_corpusLength::"<<m_corpusLength<<endl;
cout<<"m_numDocs::"<<m_numDocs<<endl;
cout<<"docLength::"<<docLength<<endl;
  //dim3 dims_Doc(nTd,nDocs,nVoc);

} //end of else 

double bestll = -numeric_limits<double>::max();
double L_plsm;
int nDocs=m_numDocs,nZ=n_Aspects,nVoc=m_corpusLength;
int nTr=nTz;//docLength
int nTs=docLength-nTz+1;    
int nTd=nTs+nTr-1;
	
//initializations

dim3 dims_d(nDocs,1,1);
unsigned int mem_size_d = dims_d.x * dims_d.y *dims_d.z* sizeof(double);
double *h_Pd = (double *) malloc(mem_size_d);
double *h_TrainPd;
	
dim3 dims_z(nZ,1,1);
unsigned int mem_size_z = dims_z.x * dims_z.y *dims_z.z* sizeof(double);
double *h_Pz = (double *) malloc(mem_size_z);

dim3 dims_trp(1,nTz,1);
unsigned int mem_size_trp = dims_trp.x * dims_trp.y *dims_trp.z* sizeof(double);
double *trp = (double *) malloc(mem_size_trp);

dim3 dims_logtrp(1,nTz+2,1);
unsigned int mem_size_logtrp = dims_logtrp.x * dims_logtrp.y *dims_logtrp.z* sizeof(double);
double *log_trp = (double *) malloc(mem_size_logtrp);

dim3 dims_pwtr_z(nZ,nVoc,nTr);
unsigned int mem_size_pwtr_z = dims_pwtr_z.x * dims_pwtr_z.y *dims_pwtr_z.z* sizeof(double);
double *p_w_tr_given_z = (double *) malloc(mem_size_pwtr_z);

double *p_w_given_tr_z = (double *) malloc(mem_size_pwtr_z);

dim3 dims_Ptszd(nZ,nDocs,nTs);
unsigned int mem_size_Ptszd = dims_Ptszd.x * dims_Ptszd.y *dims_Ptszd.z* sizeof(double);
double *h_Ptszd = (double *) malloc(mem_size_Ptszd);

dim3 dims_Pzd(nZ,nDocs,1);
unsigned int mem_size_Pzd = dims_Pzd.x * dims_Pzd.y *dims_Pzd.z* sizeof(double);
double *h_Pzd = (double *) malloc(mem_size_Pzd);
double *h_TrainPzd;


dim3 dims_Pwz(nVoc,nZ,1);
unsigned int mem_size_Pwz = dims_Pwz.x * dims_Pwz.y *dims_Pwz.z* sizeof(double);
double *h_Pwz = (double *) malloc(mem_size_Pwz);


//double *
h_Priorwz = (double *) malloc(mem_size_Pwz);

dim3 dims_Ptrwz(nVoc,nZ,nTr);
unsigned int mem_size_Ptrwz = dims_Ptrwz.x * dims_Ptrwz.y *dims_Ptrwz.z* sizeof(double);
double *h_Ptrwz = (double *) malloc(mem_size_Ptrwz);

//double *
h_Priortrwz = (double *) malloc(mem_size_Ptrwz);

dim3 dims_Pwtad(nTd,nDocs,nVoc);
unsigned int mem_size_Pwtad = dims_Pwtad.x * dims_Pwtad.y *dims_Pwtad.z* sizeof(double);
double *h_Pwtad = (double *) malloc(mem_size_Pwtad);


seed=1487059346; // temporary to test cpu result remove later

    for (int run = 0; run < multiRunCount; run++) 
	{
        printf("######### Learning ########\n");
	
	 
        printf("> Random seed is: %d\n", seed);
	 
        srand(seed);
	
        double diff = 0.0;
        double ll;

	// initializations 
        
    	// doing mem set
    	//memset(h_Ptszd,1,mem_size_Ptszd);   
    	// loading from txt file
         

    	  if (trainFlag) 
	      {
            		  if (initype == 1)
				{
				// Tammodel functionality initialize all the arrays  
                                 p_obj.compute_pd(h_Pd,h_Doc,nDocs, m_nTd , nVoc);
				 				
				 p_obj.normalize_y(h_Pd,nDocs,1,1,1);	
						
                                
				
			        p_obj.random_initialize(h_Pzd,nZ,nDocs,1,3);
				
				
				
				p_obj.normalize_y(h_Pzd,nZ,nDocs,1,1);
				
				p_obj.random_initialize(h_Ptszd,nZ,nDocs,nTs,2);
				
				
				p_obj.normalize_y(h_Ptszd,nDocs,nZ,nTs,3);
				
                                
 				
				memset(h_Pwtad,0,mem_size_Pwtad);
                                 if (opt.tr_prior)
				    {
					p_obj.tr_prior(trp,nTz) ;
					p_obj.normalize_y(trp,1,nTz,1,2) ;
					
					
				    }

				 else
				     {
  					
                                        p_obj.val_initialize(trp,1, nTz, 1,1);
					
					p_obj.normalize_y(trp,1,nTz,1,2) ;
					
				     }	
                                 
				p_obj.random_initialize(p_w_given_tr_z,nZ,nVoc,nTr,4);
				
				p_obj.normalize_y(p_w_given_tr_z,nZ,nVoc,nTr,2);
				
				
				
				const char *ext_w_trwz="w_trz.txt";
				//p_obj.file_initialize_3d1( ext_w_trwz, p_w_given_tr_z,nZ,nVoc,nTz );
				

				p_obj.compute_initial_ptrwz_pwz(h_Ptrwz,h_Pwz,p_w_tr_given_z,p_w_given_tr_z,trp,nZ,nVoc,nTr);
                                
				
				p_obj.setTrPriors(log_trp, nTz);
			
                                if(opt.align_topics) // check for answers when disp is -1 or +1 fix error
				p_obj.align_topics(h_Ptrwz,h_Ptszd,h_Pwz,log_trp,trp, m_nTd, m_nTs,nTz,nZ,nVoc, nDocs, nTs);
				
				
				
				
							
				
  			      	}
			  if (initype == 2) // 
				{  // read pwz and ptszd for initializations
				 char *filename;
					
    				 sprintf(filename, "%s_%d_%d.pwz", initfile, nZ, nTz);
					
			     
                                  p_obj.file_initialize_2d(filename, h_Pwz,nVoc,nZ);
				
			    
			        // initializing ptszd from modelfile
			        sprintf(filename, "%s_%d_%d.ptszd", initfile, nZ, nTz);
				
				p_obj.file_initialize_3d(filename, h_Ptszd, nZ, nDocs, nTs);
			       
				
				p_obj.random_initialize(h_Pzd,nZ,nDocs,1,0);
				p_obj.normalize_y(h_Pzd,nZ,nDocs,1,1);
				//p_obj.normalize(h_Pzd,nZ,nDocs,1);				
				p_obj.random_initialize(h_Ptrwz,nVoc,nZ,nTr,0);
                                p_obj.normalize_y(h_Ptrwz,nVoc,nZ,nTr,3);
				//p_obj.normalize(h_Ptrwz,nVoc,nZ,nTr);
								
				memset(h_Pwtad,0,mem_size_Pwtad);
				
				p_obj.tr_prior(trp,nTz) ;
				p_obj.normalize_y(trp,1,nTz,1,2) ;
				//p_obj.normalize(trp,1,nTz,1) ;


				p_obj.normalize_y(h_Ptszd,nZ,nDocs,nTs,3);
				p_obj.normalize_y(h_Pwz,nZ,nVoc,1,2) ;
				

				p_obj.setTrPriors(log_trp, nTz);
                                 if(opt.align_topics)
				p_obj.align_topics(h_Ptrwz,h_Ptszd,h_Pwz,log_trp,trp, m_nTd, m_nTs,nTz,nZ,nVoc, nDocs, nTs);
				
  			      	}

			  if (initype == 3)

				{
				     p_obj.compute_pd(h_Pd,h_Doc,nDocs, m_nTd , nVoc);				
				     p_obj.normalize_y(h_Pd,nDocs,1,1,1);

					
					// initializing pwz from model files with a scale
				    char filename[8095];
				    sprintf(filename, "%s_%02d_%02d.pwz", initscalePtrwz, nZ, nTz);
				    p_obj.file_initialize_2d(filename,h_Pwz, nVoc,nZ);
				    p_obj.add(h_Pwz,(double)(0.1/nVoc),nZ,nVoc,1);
				    p_obj.normalize_y(h_Pwz,nZ,nVoc,1,2);

				    p_obj.file_initialize_2d(filename,h_Priorwz, nVoc,nZ);
				    p_obj.multiply(h_Priorwz,(double)(initscalePtrwz/nTz),nZ,nVoc,1);
				    p_obj.normalize_y(h_Priorwz,nZ,nVoc,1,2);
				   
                                     	
				    // initializing ptrwz from model file
				    sprintf(filename, "%s_%02d_%02d.ptrwz", initscalePtrwz, nZ, nTz);

			            p_obj.file_initialize_3d(filename,h_Ptrwz, nVoc,nZ, nTz);
				    p_obj.add(h_Ptrwz,(double)(0.1/nTz),nVoc,nZ,nTz);
			            p_obj.normalize_y(h_Ptrwz,nVoc,nZ,nTz,3);
				    //that->ptrwz->initialize(filename);
				    //that->ptrwz->add(.1 / that->m_nTz);
				    p_obj.file_initialize_3d(filename,h_Priortrwz, nVoc,nZ, nTz);
				    p_obj.normalize_y(h_Priortrwz,nVoc,nZ,nTz,3);
				    //that->priorPtrwz->initialize(filename);


				    p_obj.random_initialize(h_Pzd,nZ,nDocs,1,0);
				    //p_obj.normalize(h_Pzd,nZ,nDocs,1);
				    p_obj.normalize_y(h_Pzd,nZ,nDocs,1,1);
				
				    p_obj.random_initialize(h_Ptszd,nZ,nDocs,nTs,0);
			            p_obj.normalize_y(h_Ptszd,nDocs,nZ,nTs,3);



				    memset(h_Pwtad,0,mem_size_Pwtad);
                                 
				    p_obj.tr_prior(trp,nTz) ;
			            p_obj.normalize_y(trp,1,nTz,1,2) ;
			            //p_obj.normalize(trp,1,nTz,1) ;
				  

				}


	      }
	   else {
            //it is inference parameter order is changed just to make the difference
		//compute pz matrix by using trained pzd file
				    
								//read parameters file
				    FILE *fid;
				    char filename[8095];
				    char str[2];
				    int t;
				    int file_nAspects,file_nTz,file_nVoc,file_maxnTd;
				    sprintf(filename, "%s_%02d_%02d.%s", modelfile, n_Aspects, nTz, "parameters");
				    fid = fopen(filename, "r");
				    if (!fid) {
					printf("cannot open file: %s\n", filename);
					exit(1);
				    }
				    int nTrainDocs;
				
				    t = fscanf(fid, "%d %d %d %d %d", &file_nVoc, &nTrainDocs, &file_maxnTd, &file_nAspects, &file_nTz);
				    if (file_nAspects != n_Aspects || file_nTz != nTz) {
					cout << "given parameters dont match with the modelfile parameters \n";
					exit(1);
				    }
				    if (m_corpusLength != file_nVoc){
				    	// reset m_corpusLength using m_nVoc
				    	m_corpusLength = file_nVoc;
				    }
				  
				
				   // initializing pwz from model file
				    sprintf(filename, "%s_%02d_%02d.pwz", modelfile, file_nAspects, file_nTz);
				  
                                    p_obj.file_initialize_2d(filename,h_Pwz, nVoc,nZ);
				
				    
				    // initializing ptrwz from modelfile
				    sprintf(filename, "%s_%02d_%02d.ptrwz", modelfile, file_nAspects,file_nTz);
				    //this->ptrwz->initialize(filename);
                                    p_obj.file_initialize_3d(filename,h_Ptrwz, nVoc,nZ, nTz);

				    // read pzd and pd from model file////////////////////////////////
				    //Array train_pzd(m_nAspects, nTrainDocs, 1);
                                    dim3 dims_trainzd(file_nAspects,nTrainDocs,1);
				   unsigned int mem_size_trainzd = dims_trainzd.x * dims_trainzd.y *dims_trainzd.z* sizeof(double);
				    h_TrainPzd = (double *) malloc(mem_size_trainzd);
				    sprintf(filename, "%s_%02d_%02d.pzd", modelfile, file_nAspects, file_nTz);
				    //train_pzd.initialize(filename);
				    p_obj.file_initialize_2d( filename, h_TrainPzd,nZ,nDocs );

				    //Array train_pd(nTrainDocs, 1, 1);
				    dim3 dims_traind(nTrainDocs,1,1);
				   unsigned int mem_size_traind = dims_traind.x * dims_traind.y *dims_traind.z* sizeof(double);
				   h_TrainPd = (double *) malloc(mem_size_traind);
                                   
				    sprintf(filename, "%s_%02d_%02d.pd", modelfile, file_nAspects, file_nTz);
				    //train_pd.initialize(filename);
				     p_obj.file_initialize_2d( filename, h_TrainPd,nTrainDocs,1 );

				    memset(h_Pz,0,mem_size_z);

				    // compute 
				    p_obj.compute_pz(h_Pz, h_TrainPzd, h_TrainPd,nTrainDocs , nZ);


				  





		
        	}
	// call plsm function here
	cout<<"Paramters are nDocs: "<<nDocs<<"and max_nTd is: "<<max_nTd<<"nTz is: "<<nTz<<"nVoc is :" <<nVoc<<"nZ is "<<nZ<<"nTs is:"<<nTs<<endl;
          
				    const char *ext_ptszd="ptszd.txt";
				    const char *ext_trwz="trwz.txt";
				    const char *ext_zd="zd.txt";
				    const char *ext_pwz="pwz.txt";
                              



			
	L_plsm=plsm_func(h_Doc,h_Pd,h_Ptszd, h_Ptrwz, h_Pzd, h_Pwz, h_Pwtad,h_Pz,trp,log_trp, nDocs, max_nTd, nTz, nVoc,nZ, nTs,m_totWords,end_accuracy,trainFlag, opt);	
 
    if (L_plsm > bestll)
	{
	// save the model
	printf("saving model now \n");
	FILE *fpt;
    	char filename_param[8095];
    	//sprintf(filename, "%s.parameters",datafile);
    	sprintf(filename_param, "%s_%02d_%02d.%s", datafile, n_Aspects, nTz, "parameters");
    	fpt = fopen(filename_param, "w");
	int tmp=75;
    	fprintf(fpt, "%d %d %d %d %d\n", tmp, nDocs, max_nTd, n_Aspects, nTz);
        printf("train flag is %d\n",trainFlag);
	if (trainFlag)fprintf(fpt, "Mode \t\t\t :%s\n", "Training");
	    if (!trainFlag)fprintf(fpt, "Mode \t\t\t :%s\n", "Inference");
	    fprintf(fpt, "Number of documents \t :%d\n", nDocs);
	    fprintf(fpt, "Number of Words \t :%d\n", nVoc);
	    fprintf(fpt, "Number of Topics \t :%d\n", n_Aspects);
	    fprintf(fpt, "Topic Length \t\t :%d\n", nTz);
	    fprintf(fpt, "Document Length \t :%d\n", max_nTd);
	    fprintf(fpt, "sparsity# \t\t : %f\n", opt.lambdaTsSparsity);
	    fprintf(fpt, "tr_prior# \t\t : %d\n", opt.tr_prior);
	    fprintf(fpt, "fit mode \t\t : %d\n", opt.fitMode);
	    fprintf(fpt, "fit params \t\t : %f %f\n", opt.fitParam1, opt.fitParam2);
	    if (!trainFlag) {
		fprintf(fpt, "pruned topics# \t\t : %d\n", opt.prune_topics);
		fprintf(fpt, "pruned topics# \t\t : %d\n", n_Aspects);
		fprintf(fpt, "z_thresh# \t\t : %f\n", opt.z_thresh);
		fprintf(fpt, "z_prior# \t\t : %d\n", opt.z_prior);
	    }
	    fclose(fpt);

	     p_obj.saveArray(datafile, (const char * )"pzd", n_Aspects, nTz, h_Pzd,nZ,nDocs,1,2 );
	     //pzd->saveArray(datafile, (char *) "pzd", m_nAspects, m_nTz);
    	     
	     p_obj.saveArray(datafile, (const char *)"ptszd", n_Aspects, nTz, h_Ptszd,nZ,nDocs,nTs,3 ); 
	     //ptszd->saveArray(datafile, (char *) "ptszd", m_nAspects, m_nTz);
    	     
	     //saveArray(datafile, const char *"pzd", n_Aspects, nTz, h_Pzd,nZ,nDocs,1 );
	     //pd->saveArray(datafile, (char *) "pd", m_nAspects, m_nTz);
             //doc_lhood->saveArray(datafile, (char *) "doc_lhood", m_nAspects, m_nTz);
	    if (trainFlag) {
		p_obj.saveArray(datafile, (const char *)"ptrwz", n_Aspects, nTz, h_Ptrwz,nVoc,nZ,nTz,3);
		//ptrwz->saveArray(datafile, (char *) "ptrwz", m_nAspects, m_nTz);
		p_obj.saveArray(datafile, (const char *)"pwz", n_Aspects, nTz, h_Pwz,nVoc,nZ,1,2 );
		//pwz->saveArray(datafile, (char *) "pwz", m_nAspects, m_nTz);
	    }
	} 
seed=seed+100;
	} // for loop	

  
}

////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
double 
plsm_func( double *h_Doc,double *h_Pd,double* h_Ptszd, double* h_Ptrwz, double *h_Pzd, double *h_Pwz, double *h_Pwtad,double *h_Pz, double *h_trp,double *log_trp, int nDocs, int max_nTd, int nTz, int nVoc,int nZ, int nTs,int totWords, float end_accuracy,bool trainflag, struct Options &opt)

{



     bool bTestResult = true;
     // update parameters here after reading command line arguements
    //int nDocs=5
    int nTr=nTz;
   // int nTd=max_nTd;
     int nTd=nTs+nTr-1;
    //int nTz=10;
    int iterations=0;
  int converged=0;
    //double L=0,
    double L_prev=0,L_Diff;

   double L_plsm=0;


   //5*5
    vector<vector<double> > Pzd(nZ,vector<double>(nDocs,rand()));
    
   
   

  

    unsigned int num_threads = 32;
    unsigned int mem_size = sizeof(double) * num_threads;
 dim3 dims_Doc(nTd,nDocs,nVoc);
unsigned int mem_size_Doc = dims_Doc.x * dims_Doc.y *dims_Doc.z* sizeof(double);


     // initialize plsm attributes
   
    dim3 dims_Ptszd(nZ,nDocs,nTs);
    unsigned int mem_size_Ptszd = dims_Ptszd.x * dims_Ptszd.y *dims_Ptszd.z* sizeof(double);
    
   

    dim3 dims_Pzd(nZ,nDocs,1);
    unsigned int mem_size_Pzd = dims_Pzd.x * dims_Pzd.y *dims_Pzd.z* sizeof(double);
  
  


    dim3 dims_Pwz(nVoc,nZ,1);
    unsigned int mem_size_Pwz = dims_Pwz.x * dims_Pwz.y *dims_Pwz.z* sizeof(double);
   


    dim3 dims_Ptrwz(nVoc,nZ,nTr);
    unsigned int mem_size_Ptrwz = dims_Ptrwz.x * dims_Ptrwz.y *dims_Ptrwz.z* sizeof(double);
   
    
   



    dim3 dims_Pwtad(nTd,nDocs,nVoc);
    unsigned int mem_size_Pwtad = dims_Pwtad.x * dims_Pwtad.y *dims_Pwtad.z* sizeof(double);

   
    
    dim3 dims_Pd(nDocs,1,1);
    unsigned int mem_size_pd = dims_Pd.x * dims_Pd.y *dims_Pd.z* sizeof(double);

    dim3 dims_Pz(nZ,1,1);
    unsigned int mem_size_pz = dims_Pz.x * dims_Pz.y *dims_Pz.z* sizeof(double);

	dim3 dims_trp(1,nTz,1);
	unsigned int mem_size_trp = dims_trp.x * dims_trp.y *dims_trp.z* sizeof(double);


   
    double *h_tszd = (double *) malloc(mem_size_Ptszd);
   
    double *h_zd = (double *) malloc(mem_size_Pzd);
   
    // Device variables and copying of data
     double *d_Ptszd;
     checkCudaErrors(cudaMalloc((void **) &d_Ptszd, mem_size_Ptszd));
     

     double *d_Pzd;
     checkCudaErrors(cudaMalloc((void **) &d_Pzd, mem_size_Pzd));

     
     double *d_Pwz;
     checkCudaErrors(cudaMalloc((void **) &d_Pwz, mem_size_Pwz));

     double *d_Priorwz;
     checkCudaErrors(cudaMalloc((void **) &d_Priorwz, mem_size_Pwz));
     
 
     double *d_Ptrwz;
     checkCudaErrors(cudaMalloc((void **) &d_Ptrwz, mem_size_Ptrwz));

     double *d_Priortrwz;
     checkCudaErrors(cudaMalloc((void **) &d_Priortrwz, mem_size_Ptrwz));
     

     double *d_Pwtad;
     checkCudaErrors(cudaMalloc((void **) &d_Pwtad, mem_size_Pwtad));
     

     double *d_Doc;
     checkCudaErrors(cudaMalloc((void **) &d_Doc, mem_size_Doc));

    double *d_Pd;
     checkCudaErrors(cudaMalloc((void **) &d_Pd, mem_size_pd));

   double *d_Pz;
     checkCudaErrors(cudaMalloc((void **) &d_Pz, mem_size_pz));
     
	double *d_trp;
     checkCudaErrors(cudaMalloc((void **) &d_trp, mem_size_trp));
     //Output     
     double *d_tszd;
     checkCudaErrors(cudaMalloc((void **) &d_tszd, mem_size_Ptszd));

     double *d_zd;
     checkCudaErrors(cudaMalloc((void **) &d_zd, mem_size_Pzd));

   

     double *d_o;
     
    double *d_idata;
    checkCudaErrors(cudaMalloc((void **) &d_idata, mem_size));
   

    // allocate device memory for result
    double *d_odata;
    checkCudaErrors(cudaMalloc((void **) &d_odata, mem_size));


   double *h_wz = (double *) malloc(mem_size_Pwz);
double *d_wz;
checkCudaErrors(cudaMalloc((void **) &d_wz, mem_size_Pwz));
double *h_wz_norm = (double *) malloc(mem_size_Pwz);


// calc trwz
double *h_trwz = (double *) malloc(mem_size_Ptrwz);
double *d_trwz;
checkCudaErrors(cudaMalloc((void **) &d_trwz, mem_size_Ptrwz));
double *h_trwz_norm = (double *) malloc(mem_size_Ptrwz);


double *h_zd_norm = (double *) malloc(mem_size_Pzd);



double *h_tszd_norm = (double *) malloc(mem_size_Ptszd);

    
    double *h_Pwtad_new = (double *) malloc(mem_size_Pwtad);
    
    //double *d_Pwtad_new;
    //checkCudaErrors(cudaMalloc((void **) &d_Pwtad_new, mem_size_Pwtad));
   	

     // Host to Device memory transfer
     checkCudaErrors(cudaMalloc((void **) &d_o, mem_size_Ptszd));

     checkCudaErrors(cudaMemcpy(d_Ptszd, h_Ptszd, mem_size_Ptszd,
                               cudaMemcpyHostToDevice));
     checkCudaErrors(cudaMemcpy(d_Pzd, h_Pzd, mem_size_Pzd,
                               cudaMemcpyHostToDevice));
     checkCudaErrors(cudaMemcpy(d_Pwz, h_Pwz, mem_size_Pwz,
                               cudaMemcpyHostToDevice));
     checkCudaErrors(cudaMemcpy(d_Ptrwz, h_Ptrwz, mem_size_Ptrwz,
                               cudaMemcpyHostToDevice));

     checkCudaErrors(cudaMemcpy(d_Pwtad, h_Pwtad, mem_size_Pwtad,
                               cudaMemcpyHostToDevice));
     checkCudaErrors(cudaMemcpy(d_Doc, h_Doc, mem_size_Doc,
                               cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMemcpy(d_tszd, h_tszd, mem_size_Ptszd,
                               cudaMemcpyHostToDevice));

        checkCudaErrors(cudaMemcpy(d_zd, h_zd, mem_size_Pzd,
                               cudaMemcpyHostToDevice));
	 checkCudaErrors(cudaMemcpy(d_Pd, h_Pd, mem_size_pd,
                               cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Pz, h_Pz, mem_size_pz,
                               cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_trp, h_trp, mem_size_trp,
                               cudaMemcpyHostToDevice));


// kernel to compute sum of doc matrix along ndoc dimension
unsigned int mem_size_logmatsum = nDocs* sizeof(double);
dim3  grid_sumdoc(nDocs, 1,1);
// loop over nTs along blockidx
dim3  threads_sumdoc(1024,1,1);



double *h_Logmat_sum = (double *) malloc(mem_size_logmatsum);
double *d_Logmat_sum;
checkCudaErrors(cudaMalloc((void **) &d_Logmat_sum, mem_size_logmatsum));



StopWatchInterface *timer3 = 0;
    sdkCreateTimer(&timer3);
	sdkStartTimer(&timer3);



dim3 grid_pwtad(nDocs,nVoc,1);
dim3 threads_pwtad(128,1,1);	
calc_pwtad<<<grid_pwtad,threads_pwtad>>>(d_Pwtad,d_Ptszd,d_Pzd,d_Pwz,d_Ptrwz,d_Pd,nDocs,nZ, nTs,nVoc,nTr);
cudaDeviceSynchronize();
checkCudaErrors(cudaMemcpy(h_Pwtad, d_Pwtad, mem_size_Pwtad,cudaMemcpyDeviceToHost));
//(nTd,nVoc,nDocs)
double sum_pwtad=0;









// kernel to compute logmat sum 
//calc_sum_doc<<<grid_sumdoc,threads_sumdoc>>>(d_Logmat_sum,d_Docsum,d_Doc,d_Pwtad,nDocs,nTs,nVoc, nTr);


calc_logmatsum<<<grid_sumdoc,threads_sumdoc>>>(d_Logmat_sum,d_Pd,d_Doc,d_Pwtad,nDocs,nTs,nVoc, nTr);
cudaDeviceSynchronize();
//checkCudaErrors(cudaMemcpy(h_Docsum, d_Docsum, mem_size_logmatsum,cudaMemcpyDeviceToHost));
checkCudaErrors(cudaMemcpy(h_Logmat_sum, d_Logmat_sum, mem_size_logmatsum,cudaMemcpyDeviceToHost));



dim3  grid_1(1, 1,1);
dim3  block_1(1,1,1);

for (int i=0;i<nDocs;i++)
L_prev=L_prev+h_Logmat_sum[i];




//cout<<"L_prev is "<<L_prev<<endl;
 printf("L_prev is ::  %0.8f \n",L_prev);





// kernel calculate likelihood


if (totlhood==1)
{
// not being used as of now
dim3  grid_lhood(1, 1,1);
dim3  threads_lhood(1,1,1);
//calc_lhood<<<grid_sumdoc,threads_sumdoc>>>(E,d_Docsum,d_Pzd,d_Ptszd,d_Pwz,d_Ptrwz,d_Pwtad,d_Doc,nDocs,nTs,nVoc, nTr,nZ);
//calc_lhood<<<grid_lhood,threads_lhood>>>(E,d_Pzd,d_Ptszd,d_Pwz, d_Ptrwz,d_Pwtad,d_Doc,nDocs, nTs,nVoc, nTr,nZ);
}



while (converged!=1)

   {   
     memset(h_Pwtad_new,0,mem_size_Pwtad);
	StopWatchInterface *timer2 = 0;
    sdkCreateTimer(&timer2);
	sdkStartTimer(&timer2);

StopWatchInterface *timerzd = 0;
    sdkCreateTimer(&timerzd);
	sdkStartTimer(&timerzd);
	
    // setup execution parameters
 
    dim3  grid_tszd(nZ, nDocs,1);
//dim3  grid_tszd(1, 1,1);
    //dim3  grid1(1, 1,1);
    dim3  threads_tszd(16, 16,1);
  
calc_tszd<<<grid_tszd,threads_tszd>>>(d_zd, d_tszd,d_Ptszd,d_Pzd,d_Pwz,d_Ptrwz,d_Pwtad,d_Doc,d_Pd,d_Pz,nDocs,nZ, nTs,nVoc, nTr,opt.lambdaTsSparsity, opt.z_prior, trainflag);
cudaDeviceSynchronize();

  sdkStopTimer(&timerzd);
        printf("Processing time for tszd n zd kernel : %f (ms)\n", sdkGetTimerValue(&timerzd));
 checkCudaErrors(cudaMemcpy(h_tszd, d_tszd, mem_size_Ptszd,cudaMemcpyDeviceToHost));



// normalize tszd matrix





  if ( trainflag )  // only in this case we compute pwz and ptrwz
    {	



dim3  grid_trwz(nVoc,nZ,1);
dim3  threads_trwz(64, 1,1);
//  m_Options.fitMode == NORMAL || m_Options.fitMode == NORMAL_UNIFORM) currently working for default mode


calc_trwz<<<grid_trwz,threads_trwz>>>(d_trwz,d_wz,d_Ptszd,d_Pzd,d_Pwz,d_Ptrwz,d_Pwtad,d_Doc,d_Pd,d_trp,d_Priortrwz,d_Priorwz,nDocs,nZ,nTs,nVoc,nTr,totWords,opt.tr_wt,opt.usePriorPtrwz,opt.tr_prior);

cudaDeviceSynchronize();
//checkCudaErrors(cudaMemcpy(h_trwz,d_trwz, mem_size_Ptrwz,cudaMemcpyDeviceToHost));
checkCudaErrors(cudaMemcpy(h_wz, d_wz, mem_size_Pwz,cudaMemcpyDeviceToHost));








dim3  grid_normwz(1, nZ,1);
dim3  threads_normwz(32, 1,1);
normalize2d_gpu<<<grid_normwz,threads_normwz>>>(d_wz,nVoc, nZ);
cudaDeviceSynchronize();
  checkCudaErrors(cudaMemcpy(h_wz_norm,d_wz, mem_size_Pwz,cudaMemcpyDeviceToHost));










dim3  grid_norm3dtrwz(nVoc,nZ,1);
dim3  threads_norm3dtrwz(64,1,1 );
//normalize3d_gpu(double *g_input,double*g_output,double *g_sum, int g_width,int g_height, int g_depth)
normalize3d_gpu<<<grid_norm3dtrwz,threads_norm3dtrwz>>>(d_trwz,nVoc,nZ,nTr); 
cudaDeviceSynchronize(); 
// ts z d    tr w z
   // 
 checkCudaErrors(cudaMemcpy(h_trwz_norm, d_trwz, mem_size_Ptrwz,cudaMemcpyDeviceToHost));
//for (int i=0;i<nTr*nVoc*nZ;i++)
//printf("trwz norm:: h_trwz_norm[%d]:: %0.8f \n",i,h_trwz_norm[i]);


dim3  grid_norm(1, nDocs,1);
dim3  threads_norm(32, 1,1);
normalize2d_gpu<<<grid_norm,threads_norm>>>(d_zd, nZ,nDocs);
   checkCudaErrors(cudaMemcpy(h_zd_norm, d_zd, mem_size_Pzd,cudaMemcpyDeviceToHost));
cudaDeviceSynchronize();

dim3  grid_norm3d(nZ,nDocs,1);
dim3  threads_norm3d(64,1,1 );
//normalize3d_gpu(double *g_input,double*g_output,double *g_sum, int g_width,int g_height, int g_depth)
normalize3d_gpu<<<grid_norm3d,threads_norm3d>>>(d_tszd,nZ,nDocs,nTs);
  cudaDeviceSynchronize();
 checkCudaErrors(cudaMemcpy(h_tszd_norm, d_tszd, mem_size_Ptszd,cudaMemcpyDeviceToHost));



// calling align topics
if (opt.align_topics){

p_obj.align_topics(h_trwz_norm,h_tszd_norm,h_wz_norm,log_trp,h_trp, m_nTd, m_nTs,nTz,nZ,nVoc, nDocs, nTs);



}
checkCudaErrors(cudaMemcpy(d_tszd, h_tszd_norm, mem_size_Ptszd,cudaMemcpyHostToDevice));
checkCudaErrors(cudaMemcpy(d_trwz, h_trwz_norm, mem_size_Ptrwz,cudaMemcpyHostToDevice));
checkCudaErrors(cudaMemcpy(d_wz, h_wz_norm, mem_size_Pwz,cudaMemcpyHostToDevice));

//p_obj.align_topics(h_Ptrwz,h_Ptszd,h_Pwz,log_trp,h_trp, m_nTd, m_nTs,nTz,nZ,nVoc, nDocs, nTs);

//nTr*nVoc*nZ
//for (int i=0;i<nTr*nVoc*nZ;i++)
//printf("trwz:: h_trwz[%d]:: %0.8f \n",i,h_trwz[i]);


checkCudaErrors(cudaMemcpy(d_Pwtad, h_Pwtad_new, mem_size_Pwtad,cudaMemcpyHostToDevice));
StopWatchInterface *timer_pwtad = 0;
    sdkCreateTimer(&timer_pwtad);
                sdkStartTimer(&timer_pwtad);

dim3 grid_pwtad(nDocs,nVoc,1);
dim3 threads_pwtad(128,1,1);	
calc_pwtad<<<grid_pwtad,threads_pwtad>>>(d_Pwtad,d_tszd,d_zd,d_wz,d_trwz,d_Pd,nDocs,nZ, nTs,nVoc,nTr);

cudaDeviceSynchronize();
  sdkStopTimer(&timer_pwtad);
            printf("Processing time for pwtad kernel : %f (ms)\n", sdkGetTimerValue(&timer_pwtad));

checkCudaErrors(cudaMemcpy(h_Pwtad_new, d_Pwtad, mem_size_Pwtad,cudaMemcpyDeviceToHost));
checkCudaErrors(cudaMemcpy(h_Ptszd, d_tszd, mem_size_Ptszd,cudaMemcpyDeviceToHost));

}  // if loop computing ptrwz and pwz  end of training flag if
else
{
dim3  grid_norm(1, nDocs,1);
dim3  threads_norm(32, 1,1);
normalize2d_gpu<<<grid_norm,threads_norm>>>(d_zd, nZ,nDocs);
   checkCudaErrors(cudaMemcpy(h_zd_norm, d_zd, mem_size_Pzd,cudaMemcpyDeviceToHost));
cudaDeviceSynchronize();

dim3  grid_norm3d(nZ,nDocs,1);
dim3  threads_norm3d(64,1,1 );
//normalize3d_gpu(double *g_input,double*g_output,double *g_sum, int g_width,int g_height, int g_depth)
normalize3d_gpu<<<grid_norm3d,threads_norm3d>>>(d_tszd,nZ,nDocs,nTs);
  cudaDeviceSynchronize();
 checkCudaErrors(cudaMemcpy(h_tszd_norm, d_tszd, mem_size_Ptszd,cudaMemcpyDeviceToHost));


calc_pwtad<<<grid_pwtad,threads_pwtad>>>(d_Pwtad,d_tszd,d_zd,d_Pwz,d_Ptrwz,d_Pd,nDocs,nZ, nTs,nVoc,nTr);
cudaDeviceSynchronize();
checkCudaErrors(cudaMemcpy(h_Pwtad_new, d_Pwtad, mem_size_Pwtad,cudaMemcpyDeviceToHost));
}


double *h_Logmat_sum_new = (double *) malloc(mem_size_logmatsum);
double *d_Logmat_sum_new;
checkCudaErrors(cudaMalloc((void **) &d_Logmat_sum_new, mem_size_logmatsum));

calc_logmatsum<<<grid_sumdoc,threads_sumdoc>>>(d_Logmat_sum_new,d_Pd,d_Doc,d_Pwtad,nDocs,nTs,nVoc, nTr);
cudaDeviceSynchronize();
checkCudaErrors(cudaMemcpy(h_Logmat_sum_new, d_Logmat_sum_new, mem_size_logmatsum,cudaMemcpyDeviceToHost));

for (int i=0;i<nDocs;i++)
L_plsm=L_plsm+h_Logmat_sum_new[i];


sdkStopTimer(&timer2);
        printf("Processing time array : %f (ms)\n", sdkGetTimerValue(&timer2));

printf(" L is %0.8f and iteration number %d\n",L_plsm,iterations);


L_Diff=(L_plsm-L_prev)/(abs(L_plsm)+EPS);

//printf("EPS is %lg and ldiff is %lg and end_Accuracy is %0.8f\n",EPS,L_Diff,end_accuracy);
//converged=1;
//if (L_Diff < epsilon || iterations > Max_iterations)
if ( ( iterations>5 && (abs(L_Diff) < end_accuracy) ) || iterations >= Max_iterations-1)
	{
	converged=1;

	printf("Iterations taken %d\n",iterations);

	checkCudaErrors(cudaMemcpy(h_Ptszd, d_tszd, mem_size_Ptszd,cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_Pzd, d_zd, mem_size_Pzd,cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_Ptrwz, d_trwz, mem_size_Ptrwz,cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_Pwz, d_wz, mem_size_Pwz,cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaFree(d_tszd));
	checkCudaErrors(cudaFree(d_zd));
        checkCudaErrors(cudaFree(d_trwz));
	checkCudaErrors(cudaFree(d_wz));
	checkCudaErrors(cudaFree(d_Pwtad));
	checkCudaErrors(cudaFree(d_Ptrwz));
	checkCudaErrors(cudaFree(d_Pwz));
	checkCudaErrors(cudaFree(d_Ptszd));

	return(L_plsm);
        
	}
	// copying the new variables
else
	{
	L_prev=L_plsm;
	L_plsm=0;
        iterations++;
	checkCudaErrors(cudaMemcpy(d_Ptszd, h_tszd_norm, mem_size_Ptszd,cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Pzd, h_zd_norm, mem_size_Pzd,cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Pwz, h_wz_norm, mem_size_Pwz,cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Ptrwz, h_trwz_norm, mem_size_Ptrwz,cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Pwtad, h_Pwtad_new, mem_size_Pwtad,cudaMemcpyHostToDevice));
	//printf("Iteration # %d\n",iterations);
	
	}
}
  

   
    cudaDeviceReset();
    exit(bTestResult ? EXIT_SUCCESS : EXIT_FAILURE);
}
