


#include "Array.h"
#include <vector>


void Array::random_initialize(double* input,int width,int height,int depth, int order)
{
	


	 if (order==1) 
	    {
		for (int d=0;d<depth;d++)
		    {
		    
		      for (int r=0;r<height;r++)
			{
		
			    
				for (int c=0;c<width;c++)
			    	input[(r*width+c)*depth+d]=1+rand()%9;
			}
		    }
	    }

	else if (order==2)
		{
		for (int r=0;r<height;r++)
		    {
		    
		      for (int d=0;d<depth;d++)
			{
		
			  for (int c=0;c<width;c++)
			    	input[(r*width+c)*depth+d]=1+rand()%9;
			}
		    }
                 } 
	else if (order==3)
		{
		for (int d=0;d<depth;d++)
		    {
		    
		      for (int c=0;c<width;c++)
			{
		
			   for (int r=0;r<height;r++) 
			    	input[(r*width+c)*depth+d]=1+rand()%9;
			}
		    }
                 }
	else if (order==4)
		{
		for (int c=0;c<width;c++)
		    {
		    
		      for (int r=0;r<height;r++) 
			{
		
			  for (int d=0;d<depth;d++) 
			    	input[(r*width+c)*depth+d]=1+rand()%9;
			}
		    }
                 }

		
}





void Array::val_initialize(double* input,int width,int height,int depth,double val)
{
	for (int r=0;r<height;r++)
	    {
	    
	      for (int c=0;c<width;c++)
		{
		
		    for (int d=0;d<depth;d++)
		    	input[(r*width+c)*depth+d]=val;
		}
	    }
		
}

void Array::add(double *input, double val, int width, int height, int depth) {
  for (int r=0;r<height;r++)
	    {
	    
	      for (int c=0;c<width;c++)
		{
		
		    for (int d=0;d<depth;d++)
		    	input[(r*width+c)*depth+d]+=val;
		}
	    }
}


void Array::multiply(double *input, double val, int width, int height, int depth) {
  for (int r=0;r<height;r++)
	    {
	    
	      for (int c=0;c<width;c++)
		{
		
		    for (int d=0;d<depth;d++)
		    	input[(r*width+c)*depth+d]*=val;
		}
	    }
}


void Array::file_initialize(const char *filename, double *input, int ncols, int nrows, int ndepth) {

    int t;
    double temp;
    FILE *fid;
    //int nrows = this->m_height;
    //int ncols = this->m_width;
    //int ndepth = this->m_depth;

    printf("reading file: %s\n", filename);
    fid = fopen(filename, "r");
    if (!fid) {
        printf("cannot open file: %s\n", filename);
        exit(1);
    }
    for (int d = 0; d < ndepth; d++) {
        for (int r = 0; r < nrows; r++) {
            for (int c = 0; c < ncols; c++) {
                t = fscanf(fid, "%lg ", &temp);
                //this->mat[r][c][d] = temp;
		input[(r*ncols+c)*ndepth+d]=temp;
            }
        }
    }

    fclose(fid);


}


void Array::file_initialize_3d1(const char *filename, double *input, int ncols, int nrows, int ndepth) {

    int t;
    double temp;
    FILE *fid;
    //int nrows = this->m_height;
    //int ncols = this->m_width;
    //int ndepth = this->m_depth;

    printf("reading file: %s\n", filename);
    fid = fopen(filename, "r");
    if (!fid) {
        printf("cannot open file: %s\n", filename);
        exit(1);
    } 
	for (int c = 0; c < ncols; c++) {

    
        for (int r = 0; r < nrows; r++) {
		
            
			for (int d = 0; d < ndepth; d++) {
                t = fscanf(fid, "%lg ", &temp);
                //this->mat[r][c][d] = temp;
		input[(r*ncols+c)*ndepth+d]=temp;
            }
        }
    }

    fclose(fid);


}




void Array::file_initialize_3d(const char *filename, double *input, int ncols, int nrows, int ndepth) {

    int t;
    double temp;
    FILE *fid;
    //int nrows = this->m_height;
    //int ncols = this->m_width;
    //int ndepth = this->m_depth;

    printf("reading file: %s\n", filename);
    fid = fopen(filename, "r");
    if (!fid) {
        printf("cannot open file: %s\n", filename);
        exit(1);
    }
    
        for (int r = 0; r < nrows; r++) {
		for (int d = 0; d < ndepth; d++) {
            for (int c = 0; c < ncols; c++) {
                t = fscanf(fid, "%lg ", &temp);
                //this->mat[r][c][d] = temp;
		input[(r*ncols+c)*ndepth+d]=temp;
            }
        }
    }

    fclose(fid);


}

void Array::file_initialize_2d(const char *filename, double *input, int ncols, int nrows) {

    int t;
    double temp;
    FILE *fid;
    //int nrows = this->m_height;
    //int ncols = this->m_width;
    //int ndepth = this->m_depth;

    printf("reading file: %s\n", filename);
    fid = fopen(filename, "r");
    if (!fid) {
        printf("cannot open file: %s\n", filename);
        exit(1);
    }
    
        
	
            for (int c = 0; c < ncols; c++) {
		for (int r = 0; r < nrows; r++) {
                t = fscanf(fid, "%lg ", &temp);
                //this->mat[r][c][d] = temp;
		input[(r*ncols+c)]=temp;
           
        }
    }

    fclose(fid);


}


Array::~Array() {

}


void Array::normalize(double* input,int width,int height,int depth) {


    for (int d = 0; d < depth; d++) {
        for (int c = 0; c < width; c++) {
            double sum = 0.0;
            for (int r = 0; r < height; r++) {
                sum += input[(r*width+c)*depth+d];
            }
            for (int r = 0; r < height; r++) {
                input[(r*width+c)*depth+d]/=(sum+EPS);
		//this->mat[r][c][d] /= (sum + EPS);
            }
        }
    }
    
}


void Array::normalize_y(double* input,int width,int height,int depth, int dim) {

if (dim==2){

	    for (int d = 0; d < depth; d++) {
		for (int c = 0; c < width; c++) {
		    double sum = 0.0;
		    for (int r = 0; r < height; r++) {
		        sum += input[(r*width+c)*depth+d];
			
		    }
		    for (int r = 0; r < height; r++) {
		        input[(r*width+c)*depth+d]/=(sum+EPS);
			//this->mat[r][c][d] /= (sum + EPS);
			
		    }
		}
	    }

 	}

if (dim==1){

	    for (int d = 0; d < depth; d++) {
		for (int r = 0; r < height; r++) {
		    double sum = 0.0;
		    for (int c = 0; c < width; c++) {
		        sum += input[(r*width+c)*depth+d];
		    }
		   for (int c = 0; c < width; c++) {
		        input[(r*width+c)*depth+d]/=(sum+EPS);
			//this->mat[r][c][d] /= (sum + EPS);
		    }
		}
	    }

 	}

if (dim==3){

	    for (int c = 0; c < width; c++) {
		for (int r = 0; r < height; r++) {
		    double sum = 0.0;
		   for (int d = 0; d < depth; d++)  {
		        sum += input[(r*width+c)*depth+d];
		    }
		   for (int d = 0; d < depth; d++) {
		        input[(r*width+c)*depth+d]/=(sum+EPS);
			//this->mat[r][c][d] /= (sum + EPS);
		    }
		}
	    }

 	}


    
}






void Array::normalize2(double* input,int width,int height,int depth) {
    for (int d = 0; d < depth; d++) {
        double sum = 0.0;
        for (int c = 0; c < width; c++) {
            for (int r = 0; r < height; r++) {
                sum += input[(r*width+c)*depth+d];
            }
        }
        for (int c = 0; c < width; c++) {
            for (int r = 0; r < height; r++) {
                input[(r*width+c)*depth+d] /= (sum + EPS);
            }
        }
    }
    
}


void Array::saveArray(const char* datafile, const char *ext, int n_Aspects, int n_Tz, double* output, int width,int height,int depth, int order ) {

    FILE* fid;
    char filename[8095];

    int nrows = height;
    int ncols = width;
    int ndepth = depth;

    sprintf(filename, "%s_%02d_%02d.%s", datafile, n_Aspects, n_Tz, ext);
    fid = fopen(filename, "w");
    printf("%s \n", filename);

    if (order == 2)
 	for (int c = 0; c < ncols; c++) {
    	    for (int d = 0; d < ndepth; d++) {
       
            for (int r = 0; r < nrows; r++) {
                fprintf(fid, "%g ", output[(r*ncols+c)*ndepth+d]);
		//fprintf(fid, "%g ", this->mat[r][c][d]);  column elements are written first
            }
            fprintf(fid, "\n");
        }
        //fprintf(fid,"\n");
    }

   else if (order ==3)
   for (int r = 0; r < nrows; r++) {
    	   for (int d = 0; d < ndepth; d++) {
       
             for (int c = 0; c < ncols; c++) {
                fprintf(fid, "%g ", output[(r*ncols+c)*ndepth+d]);
		//fprintf(fid, "%g ", this->mat[r][c][d]);  column elements are written first
            }
            fprintf(fid, "\n");
        }
        //fprintf(fid,"\n");
    }

    fclose(fid);

}



void Array::readArray(const char* datafile, const char *ext, int n_Aspects, int n_Tz, double* output, int width,int height,int depth ) {

    FILE* fid;
    char filename[8095];

    int nrows = height;
    int ncols = width;
    int ndepth = depth;

    sprintf(filename, "%s_%02d_%02d.%s", datafile, n_Aspects, n_Tz, ext);
    fid = fopen(filename, "r");
    printf("%s \n", filename);

    for (int d = 0; d < ndepth; d++) {
        for (int r = 0; r < nrows; r++) {
            for (int c = 0; c < ncols; c++) {
                fscanf(fid, "%lf ", &output[(r*ncols+c)*ndepth+d]);
		//fprintf(fid, "%g ", this->mat[r][c][d]);  column elements are written first
            }
            fprintf(fid, "\n");
        }
        //fprintf(fid,"\n");
    }

    fclose(fid);

}


void Array::print(double* input,int width,int height,int depth) {

    

    for (int d = 0; d < depth; d++) {
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                printf("%lf ", input[(r*width+c)*depth+d]);
            }
            printf("\n");
        }
        printf("\n");
    }
}



double Array::sum(double* input,int width,int height,int depth) {

    double sum=0;

    for (int d = 0; d < depth; d++) {
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                sum=sum+ input[(r*width+c)*depth+d];
            }
            //printf("\n");
        }
        //printf("\n");
    }
   return(sum);
}

void Array::compute_pz(double *Pz, double *train_pzd, double *train_pd,int TrainDocs , int nZ)
	{
	 double temp;
				   // pz->initialize(0.0);
				    for (int k = 0; k < nZ; k++) {
					temp = 0.0;
					for (int d = 0; d < TrainDocs; d++) {
					   // temp = temp + train_pzd.mat[k][d][0] * train_pd.mat[d][0][0];
					    temp = temp + train_pzd[d*nZ+k] * train_pd[d];
					}
					//pz->mat[k][0][0] = temp;
					Pz[k]=temp;
				    }
				    //pz->normalize();
				    normalize_y(Pz,nZ,1,1,1);

        
	}




void Array::compute_pd(double *h_Pd, double *h_Doc, int nDocs, std::vector<int>& m_nTd , int nVoc)
	{
     
	for (int d = 0; d < nDocs; d++) {
        h_Pd[d]=0;
	
        for (int ta = 0; ta < m_nTd[d]; ta++) {
		
            for (int w = 0; w < nVoc; w++) {
		//printf("d is %d and ta is %d and w is %d and m_nTd[d] is %d\n",d,ta,w,m_nTd[d]);
		//if (d==nDocs-1)		
		//printf("h_Doc[%d] is %f and h_Pd[d] is %f\n",(d*m_nTd[d]+ta)*nVoc+w,h_Doc[(d*m_nTd[d]+ta)*nVoc+w],);
		h_Pd[d]+=h_Doc[(d*m_nTd[d]+ta)*nVoc+w];                
		//pd->mat[d][0][0] = pd->mat[d][0][0] + (double) corpus->count(w, ta, d);
            }
        }
        //sum = sum + pd->mat[d][0][0];
	
    	}

        
	}


void Array::tr_prior(double* input,int height) 
	{


    double uni = 1;
    double gamma = 0.3;
    double base = 3;
    int sizeRT = height;
    input[0] = base;
    for (int i = 1; i < sizeRT; i++) {
//        this->mat[i][0][0] = (base + uni + gamma * (sizeRT - i));
          input[i] = (base + uni + gamma * (sizeRT - i));
    }

}


void Array::compute_initial_ptrwz_pwz(double* h_Ptrwz,double* h_Pwz,double*p_w_tr_given_z,double* p_w_given_tr_z,double *trp, int nZ,int nVoc,int nTr)
{
     // h_Ptrwz x,y,z  dims_pwtr_z(nVoc,nZ,nTr)
     // h_Pwz dims_Pwz(nZ,nVoc,1);
     // p_w_tr_given_z dims_pwtr_z(nZ,nVoc,nTr)
     //	p_w_given_tr_z	dims_pwtr_z(nZ,nVoc,nTr)		

printf("nZ is %d, nTr is %d, nVoc is %d\n",nZ,nTr,nVoc);
		for (int z = 0; z < nZ; z++) {
        for (int w = 0; w < nVoc; w++) {
            double sum = 0.0;
            for (int tr = 0; tr < nTr; tr++) {
		//p_w_tr_given_z[ (z*nVoc+w)*nTr+tr]= p_w_given_tr_z[(z*nVoc+w)*nTr+tr]* trp[tr];
		p_w_tr_given_z[ (z+w*nZ)*nTr+tr]= p_w_given_tr_z[(z+w*nZ)*nTr+tr]* trp[tr];
		sum+= p_w_tr_given_z[ (z+w*nZ)*nTr+tr];
		//sum+= p_w_tr_given_z[ (z*nVoc+w)*nTr+tr];
		//printf("p_w_tr_given_Z[%d][%d][%d]  is %lg and trp[tr] is %lg\n",z,w,tr, p_w_tr_given_z[(z+w*nZ)*nTr+tr],trp[tr] );
		//p_w_tr_given_z.mat[w][tr][z] = p_w_given_tr_z.mat[w][tr][z] * trp->mat[tr][0][0];
                //sum += p_w_tr_given_z.mat[w][tr][z];
            }
             //h_Pwz[w*nZ+z]=sum; 
	    h_Pwz[w+z*nVoc]=sum; 
	    // printf("h_Pwz[%d]:: %0.8f\n",w+z*nVoc,h_Pwz[w+z*nVoc]);	
	     //this->pwz->mat[w][z][0] = sum;
        }
    }
       // this->pwz->normalize();
	
	normalize_y(h_Pwz,nVoc,nZ,1,1) ;
	//normalize(h_Pwz,nVoc,nZ,1) ;
    for (int z = 0; z < nZ; z++) {
        for (int w = 0; w < nVoc; w++) {
            for (int tr = 0; tr < nTr; tr++) {
                //h_Ptrwz[(z*nVoc+w)*nTr+tr]=p_w_tr_given_z[(w*nZ+z)*nTr+tr]/(h_Pwz[w*nZ+z]+EPS);
		h_Ptrwz[(z*nVoc+w)*nTr+tr]=p_w_tr_given_z[(w*nZ+z)*nTr+tr]/(h_Pwz[z*nVoc+w]+EPS);
		//this->ptrwz->mat[tr][w][z] = p_w_tr_given_z.mat[w][tr][z]/(this->pwz->mat[w][z][0]+EPS);
		
		}	        
	}
    }
	
normalize_y(h_Ptrwz,nVoc,nZ,nTr,3);

//normalize(h_Ptrwz,nVoc,nZ,nTr);
 }




void Array::setTrPriors(double *log_trp, int nTz){

    double epsilon = 1;
    double base = 3;
    double uni = 1;
    double gamma = 0.3;
    log_trp[0]=epsilon;	
    double sum=base+2*epsilon;
    //this->log_trp->mat[-1 + 1][0][0] = epsilon;

    log_trp[1]=base;
    //this->log_trp->mat[0 + 1][0][0] = base;
    for (int i = 1; i < nTz; i++) {
	log_trp[i +1] = (base + uni + gamma * (nTz + 2 - i));
	sum+=log_trp[i +1];
        //this->log_trp->mat[i + 1][0][0] = (base + uni + gamma * (m_nTz + 2 - i));
    }
    log_trp[nTz + 1] = epsilon;
    //this->log_trp->mat[m_nTz + 1][0][0] = epsilon;
    //this->log_trp->print();
    //normalize(log_trp,1,nTz+2,1);
    //normalize_y(log_trp,1,nTz+2,1,2);
    //this->log_trp->normalize();
    //this->log_trp->print();

	// not calling normalization function 
    for (int i = 0; i < nTz + 2; i++) {
	
	log_trp[i] = log((log_trp[i]/sum)+EPS);
	//printf("log_trp val is %0.8f\n",log_trp[i]);
        //this->log_trp->mat[i][0][0] = log(this->log_trp->mat[i][0][0]+EPS);
    }
    
}


void Array::align_topics(double *ptrwz,double* ptszd,double *pwz,double *log_trp,double *trp, std::vector<int>& m_nTd, std::vector<int>& m_nTs, int nTz,int nAspects, int nVoc, int nDocs, int nTs) {

    

   // Array ptrz(m_nTz, m_nAspects, 1);
   // Array p_w_tr_given_z(m_nVoc,m_nTz,m_nAspects);

   unsigned int mem_size_ptrz = nAspects * nTz *1* sizeof(double);
   double *ptrz = (double *) malloc(mem_size_ptrz);
   //dim3 dims_Pwz(nZ,nVoc,1); 
   //dim3 dims_Ptrwz(nVoc,nZ,nTr);	
   //dim3 dims_pwtr_z(nZ,nVoc,nTz/ntr);
   unsigned int mem_size_pwtrz = nAspects* nVoc *nTz * sizeof(double);
   double *p_w_tr_given_z = (double *) malloc(mem_size_pwtrz);
   

    //marginalize to get ptrz
    for (int tr = 0; tr < nTz; tr++) {
        for (int z = 0; z < nAspects; z++) {
            //ptrz.mat[tr][z][0] = 0;
	      ptrz[tr*nAspects+z]=0;	
            for (int w = 0; w < nVoc; w++){

		ptrz[tr*nAspects+z] +=  ptrwz[(z*nVoc+w)*nTz+tr] * pwz[w+z*nVoc];
                p_w_tr_given_z[(w*nAspects+z)*nTz+tr] = ptrwz[(z*nVoc+w)*nTz+tr] * pwz[w+z*nVoc];
		
                //ptrz.mat[tr][z][0] +=  ptrwz->mat[tr][w][z] * pwz->mat[w][z][0];
                //p_w_tr_given_z.mat[w][tr][z] = ptrwz->mat[tr][w][z] * pwz->mat[w][z][0];
            }
          
        }
    }
    //ptrz.normalize();
    //normalize(ptrz,nTz,nAspects,1);
 
    normalize_y(ptrz,nAspects,nTz,1,2);
 //printf("ptrz norm is %0.8f\n",ptrz[0]);
//printf("ptrz norm is %0.8f\n",ptrz[1]); 
	
  
    //dim3 dims_pwtr_z(nZ,nVoc,nTz/ntr);
    //Array temp_pwtrz = p_w_tr_given_z;
    //Array *temp_ptszd = ptszd;
//    double *temp_pwtrz=p_w_tr_given_z;
double *temp_pwtrz = (double *) malloc(mem_size_pwtrz);
copy(temp_pwtrz,p_w_tr_given_z,nAspects, nVoc, nTz );
//temp_pwtrz=p_w_tr_given_z;
//double *temp_ptszd = (double *) malloc( nAspects* nDocs *nTz * sizeof(double));
    double *temp_ptszd=ptszd;
    //dim3 dims_Ptszd(nZ,nDocs,nTs);



    
    double corr[3];
    for (int z = 0; z < nAspects; z++) {

        // find the displacements for each disp
        for (int disp = -1; disp <= +1; disp++) {
            corr[disp + 1] = 0;
            int ofs = disp+1;
            for (int tr = 0; tr < nTz; tr++) {// -1
                //corr[disp + 1] +=  ptrz.mat[tr][z][0] * log_trp->mat[tr + ofs][0][0];
		corr[disp + 1] +=  ptrz[tr*nAspects+z] * log_trp[tr + ofs];
		//printf("ptrz is %0.8f and log_trp is  %0.8f\n",ptrz[tr*nAspects+z], log_trp[tr + ofs]);
            }
            
        }
        //find max in corr
	int final_disp=-1;
        float m = corr[0]; 
        if (corr[1]>=corr[0]){
                m = corr[1];
                final_disp=0;
        }
        if (corr[2]>corr[1]){
                m = corr[2];
                final_disp=1;
        }
        //cout<<corr[0]<<"\t"<<corr[1]<<"\t"<<corr[2]<<endl;
        //cout<<"max :"<<m<<"\t"<<"disp :"<<final_disp;
	//printf("z is %d and corr 0 is %0.8f corr 1 is %0.8f and corr 2 is %0.8f\n",z,corr[0],corr[1],corr[2]);
        
        if (final_disp == -1) {
           //cout<<"final_disp :"<<final_disp<<endl;
		printf("final_disp %d\n",final_disp);
           for(int w = 0; w < nVoc; w ++) {
                //temp_pwtrz.mat[w][m_nTz - 1][z] = (1/m_nVoc)*trp->mat[m_nTz-1][0][0];
		temp_pwtrz[(w*nAspects+z)*nTz+nTz-1] = (1/nVoc)*trp[nTz-1];
                for (int tr = 0; tr < nTz - 1; tr++) {
                    //temp_pwtrz.mat[w][tr][z] = p_w_tr_given_z.mat[w][tr + 1][z];
		    temp_pwtrz[(w*nAspects+z)*nTz+tr] = p_w_tr_given_z[(w*nAspects+z)*nTz+tr+1]; 	
                }
           }
            for(int d =0; d < nDocs; d++){
                for (int ts = m_nTs[d]-1; ts > 0; ts--){
		     //dim3 dims_Ptszd(nZ,nDocs,nTs);	
                     //temp_ptszd->mat[ts][z][d] =  temp_ptszd->mat[ts-1][z][d];
		     temp_ptszd[(d*nAspects+z)*nTs+ts] =  temp_ptszd[(d*nAspects+z)*nTs+ts-1];
                }
            }

        } else if (final_disp == +1) {
            //cout<<"final_disp :"<<final_disp<<endl;
	    printf("final_disp %d\n",final_disp);
	  // printf("pw_tr_given_z is %0.8f",p_w_tr_given_z[0]);
            for(int w = 0; w < nVoc; w ++) {
                 //temp_pwtrz.mat[w][0][z] = (1/m_nVoc)*trp->mat[0][0][0];
		 temp_pwtrz[(w*nAspects+z)*nTz] = (1/nVoc)*trp[0];
		// temp_pwtrz[(w*nAspects+z)*nTz+1]=p_w_tr_given_z[(w*nAspects+z)*nTz];
                for (int tr = nTz-1; tr > 0 ; tr--) {
		// for (int tr = 1; tr <nTz ; tr++) {
                    //temp_pwtrz.mat[w][tr][z] = p_w_tr_given_z.mat[w][tr-1][z];
		    temp_pwtrz[(w*nAspects+z)*nTz+tr] = p_w_tr_given_z[(w*nAspects+z)*nTz+tr-1];
			
                }
            }
		
		//printf("pw_tr_given_z is %0.8f",p_w_tr_given_z[0]);
		//printf("temp_pwtrz is %0.8f",temp_pwtrz[1]);

            for(int d =0; d < nDocs; d++){
		//printf("m_nTs[%d] is %d\n",d,m_nTs[d]);
                for (int ts = 0; ts < m_nTs[d]-1; ts++){
                    //temp_ptszd->mat[ts][z][d] =  temp_ptszd->mat[ts+1][z][d];
			
		    temp_ptszd[(d*nAspects+z)*nTs+ts] =  temp_ptszd[(d*nAspects+z)*nTs+ts+1];
                }
            }
	//for(int d =0; d < 5; d++)
	//printf("temp_pwtrz at ts %d is %0.8f\n",d,temp_pwtrz[d]);
		
        }
   }
   ptszd = temp_ptszd;
//printf("ptszd at 0  is %0.8f\n",ptszd[(0*nAspects+0)*nTs]);
//printf("ptszd at 1 is %0.8f\n",ptszd[(0*nAspects+1)*nTs]);
//printf("ptszd at 2 is %0.8f\n",ptszd[(0*nAspects+2)*nTs]);
   
   normalize_y(ptszd,nAspects,nDocs,nTs,3);
   //normalize(ptszd,nAspects,nDocs,nTs);
   //ptszd->normalize();

   for (int z = 0; z < nAspects; z++) {
        for (int w = 0; w < nVoc; w++) {
            double sum = 0.0;
            for (int tr = 0; tr < nTz; tr++) {
   
               // sum += temp_pwtrz.mat[w][tr][z];
		sum += temp_pwtrz[(w*nAspects+z)*nTz+tr];
            }
            //this->pwz->mat[w][z][0] = sum;
	    //printf("sum at z %d and w %d is %0.8f\n",z,w,sum);
	    pwz[w+z*nVoc] = sum;
           
        }
    }
    //this->pwz->normalize(); 
    //normalize(pwz,nAspects,nVoc,1);
    normalize_y(pwz,nVoc,nAspects,1,1);
    
    	

    //dims_Ptrwz(nVoc,nAspects,nTz);
    //mem_size_pwtrz = nApsects* nVoc *nTz	
	
    for (int z = 0; z < nAspects; z++) {
        for (int w = 0; w < nVoc; w++) {
            for (int tr = 0; tr < nTz; tr++)
                //this->ptrwz->mat[tr][w][z] = temp_pwtrz.mat[w][tr][z]/(this->pwz->mat[w][z][0]+EPS);
		ptrwz[(z*nVoc+w)*nTz+tr]=temp_pwtrz[(w*nAspects+z)*nTz+tr]/(pwz[w+z*nVoc]+EPS);
        }
    }
    //ptrwz->normalize();
    //normalize(ptrwz,nAspects,nVoc,nTz);
    normalize_y(ptrwz,nVoc,nAspects,nTz,3);


//    printf(" ptszd 0 is %lf \n",ptszd[0]); 
//printf(" ptszd 1 is %lf \n",ptszd[1]);
}


void Array:: copy(double *output,double *input,int x, int y, int z )
{

for (int i=0;i<x*y*z;i++)
output[i]=input[i];
}



