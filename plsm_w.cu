// includes, system


#include <stdio.h>
using namespace std; 

#include <float.h>





#include <sys/stat.h>
#include <limits>


//#include "cuPrintf.cu"

// includes CUDA
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// includes, project
//#include <helper_cuda.h>
//#include <helper_functions.h> // helper functions for SDK examples



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

tmp_pd=g_Pd[blockIdx.x];
/*
__shared__ double tmp_pwz[block_height_y];
__shared__ double tmp_pzd[block_height_y];
__shared__ double tmp_ptrwz[block_height_y];
*/
int nTs_counter=ceil((float)g_nTs/(float)block_width_x);

int index_Ptszd;
for (int tr=0;tr<g_nTr;tr++)
{   

  for (int w=0;w<g_nVoc;w++)
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
                    //tmp_pwz[threadIdx.y]=0;
                   // tmp_pzd[threadIdx.y]=0;
		  //  tmp_ptrwz[threadIdx.y]=0;
                    
                    __syncthreads();
		index_Ptszd=(blockIdx.x*g_nZ+z)*g_nTs+block_width_x*ts+threadIdx.x;	 
			/*
                                           
			tmp_pwz[threadIdx.y]=g_Pwz[w*g_nz+z*blockDim.y+threadIdx.y];
                        tmp_pzd[threadIdx.y]=g_Pzd[blockIdx.x*g_nz+z*blockDim.y+threadIdx.y];
                        tmp_ptrwz[threadIdx.y]=g_Ptrwz[((z*blockDim.y+threadIdx.y)*g_nVoc+w)*g_nTr+tr];
			
			*/
                        if ((block_width_x*ts+threadIdx.x)<g_nTs )
                           {
			
                        	tmp_Ptszd[threadIdx.x]= g_Ptszd[index_Ptszd];
                            }
			__syncthreads();

                     

                          
                        
			//tmp_pwtad[threadIdx.y][threadIdx.x]=tmp_pzd[threadIdx.y]*tmp_Ptszd[threadIdx.y][threadIdx.x]*tmp_pwz[threadIdx.y]*tmp_ptrwz[threadIdx.y];
                        //__syncthreads();
			if 	 ((block_width_x*ts+threadIdx.x)< g_nTs)
			

			{
			
			tmp_pwtad[threadIdx.x]+=tmp_pd*tmp_pzd*tmp_Ptszd[threadIdx.x]*tmp_pwz*tmp_ptrwz;
			
			}
                        __syncthreads();
                          //if (threadIdx.x==0 && threadIdx.y==0)
 			//printf("tmp_pwtad[0][0]::%0.8f \n",tmp_pwtad[0][0]);  
                       

                    } // nz loop  

			Td_index=tr+block_width_x*ts+threadIdx.x;
			//-((tr+block_width_x*ts+threadIdx.x)/(g_nTr+g_nTs-1));
    
                          if(Td_index< (g_nTs+g_nTr-1))
                	     {
				
			
			        pwtad_index=(blockIdx.x*(g_nTs+g_nTr-1)+Td_index)*g_nVoc+w;
                              //g_Pwtad[(w*(g_nTr+g_nTs-1)+(tr+blockDim.x*ts+threadIdx.x))*g_nDocs+blockIdx.x]+=tmp_pwtad[0][threadIdx.x];
				//g_Pwtad[(blockIdx.x*(g_nTr+g_nTs-1)+tr+blockDim.x*ts+threadIdx.x -((tr+blockDim.x*ts+threadIdx.x)/(g_nTr+g_nTs-2)))*g_nVoc+w]+=tmp_pwtad[0][threadIdx.x];
                               g_Pwtad[pwtad_index]+=tmp_pwtad[threadIdx.x];
				
			    }	
                             __syncthreads();
						
			
				
                 		
                	     

                        
                   
               } // nts counter loop
 

	}// tr loop
	}// w loop



}




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
__shared__ double shared_ptszd;
__shared__ double shared_pzd;
__shared__ double shared_pd;
        
	for (int ts=0;ts<g_nTs;ts++) 
            {
                //tmp_sum=0;
                  tmp_tszd=0;
		   shared_ptszd=g_Ptszd[(blockIdx.y*g_nZ+blockIdx.x)*g_nTs+ts];
		   shared_pzd=g_Pzd[(blockIdx.y)*g_nZ+blockIdx.x];
		   shared_pd=pd[blockIdx.y];
		  __syncthreads();
		    for (unsigned int nvoc_counter=0;nvoc_counter<iterations_nvoc;nvoc_counter++)
		     {
			 for (unsigned int ntr_counter=0;ntr_counter<iterations_ntr;ntr_counter++)
		     		{
			partial_tszd[threadIdx.y*block_width_x+threadIdx.x]=0;
			__syncthreads();	
			if ((ntr_counter*blockDim.y+threadIdx.y)<g_nTr && (nvoc_counter*blockDim.x+threadIdx.x)<g_nVoc )
			{	gdoc_index=(blockIdx.y*(g_nTs+g_nTr-1)+ts+ntr_counter*blockDim.y+threadIdx.y-(ts+ntr_counter*blockDim.y+threadIdx.y)/(g_nTs+g_nTr-1))*g_nVoc+nvoc_counter*blockDim.x+threadIdx.x;
                           
			partial_tszd[threadIdx.y*block_width_x+threadIdx.x]=(shared_pd*g_Doc[gdoc_index]*shared_pzd* shared_ptszd*g_Pwz[(blockIdx.x)*g_nVoc+nvoc_counter*blockDim.x+threadIdx.x]*g_Ptrwz[(blockIdx.x*g_nVoc+nvoc_counter*blockDim.x+threadIdx.x)*g_nTr+ntr_counter*blockDim.y+threadIdx.y])/(g_Pwtad[gdoc_index]+DBL_EPSILON);
	                 
		
                         } 
			__syncthreads();
                        for (unsigned int s=block_width_x*block_height_y/2;s>=1;s>>=1)
       				{//g_tszd[s]=s;
       				 
				  if (threadIdx.y*block_width_x+threadIdx.x<s)
	    			  {
	     			  
             			    partial_tszd[threadIdx.y*block_width_x+threadIdx.x]+=partial_tszd[threadIdx.y*block_width_x+threadIdx.x+s];
            			    //shared_partial_tszd[0][][threadIdx.y]+=shared_partial_tszd[idx][threadIdx.z][threadIdx.y];
             			    //sdata[threadIdx.x]+=sdata[threadIdx.x+s];
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
                             
			// if ( blockIdx.x==0 && blockIdx.y==0 && blockIdx.z==0 && d==0 && nts_counter==0 )
			//	printf("tid x is  %d and g_Doc[gdoc_index] is  %0.8f and partial_trwz[threadIdx.x] is %0.8f \n",threadIdx.x,g_Doc[gdoc_index],partial_trwz[threadIdx.x]);
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
	   {    //double more = priorPwz->mat[w][z][0] * priorPtrwz->mat[tr][w][z] * m_totWords;
		tmp_sum_outerloop+=priortrwz[(blockIdx.y*g_nVoc+blockIdx.x)*g_nTr+tr]*priorwz[blockIdx.y*g_nVoc+blockIdx.x]*totWords;
	    }
	else if (tr_prior){
		//temp1.mat[tr][w][z] = sum_trwz + m_Options.tr_wt * (sum_count * trp->mat[tr][0][0]) / (m_nVoc * m_nAspects);
 		//sum_wz += temp1.mat[tr][w][z];
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



