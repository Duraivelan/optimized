# include "structure_definitions.h"
#include <sstream>

const int n_vars = 1 ; 									// number of different properties you want to autocorrelate
int pcor = 64 ;
int p2 = pcor/2 ; 								//  p2 = pcor/2
const int Max_level = 22;
double acor[n_vars][Max_level][64]={-2.0};
double fcor[n_vars][Max_level][32]={0};
double ncor[n_vars][Max_level][32]={0};
int pointcor[n_vars][Max_level]={0};

void zerocor() {

	for (int  i=0; i < n_vars; i++ ) 
		{       
			for (int  j=0; j < Max_level ; j++ )
				{
					for (int  k=0; k < p2; k++ )
						{
							acor[i][j][k] = -2.0e9 ; 
							fcor[i][j][k] =  0.0 ; 
							ncor[i][j][k] =  0.0 ; 
						}
						
					for (int  k=p2; k < pcor; k++ )
						{
							acor[i][j][k] = -2.0e9 ; 
						}	
										
					pointcor[i][j] = 0 ; 	
				}
		}
}		


void addcor(double new_data, int nf, int k)		// nf is data_type identifer , k is the correlator_level
{
	int i, point, j;
	
// shift pointer and put in new_data
	point 				= ( pointcor[nf][k] ) % pcor ;				// modulo 
	pointcor[nf][k]		= point + 1 ;
	acor[nf][k][point]	= new_data ;

// do the correlation
	for(i=0; i<p2; i++)
		{
			j = (pcor+point-i-p2)%pcor  ;
			if(acor[nf][k][j]>-2.0e9) 				// array must be filled to do correlation
				{
					fcor[nf][k][i]=fcor[nf][k][i]+new_data*acor[nf][k][j] ;
					ncor[nf][k][i]=ncor[nf][k][i]+1.0 ;
				}
		}
      
      
// in 0 we put the first 7 values
	if (k==1)
		{
			for(i=0; i<p2-1; i++)
				{
					j = (pcor+point-i-1)%pcor  ;
					if(acor[nf][1][j]>-2.0e9)
						{
							fcor[nf][0][i]=fcor[nf][0][i]+acor[nf][1][point]*acor[nf][1][j] ;
							ncor[nf][0][i]=ncor[nf][0][i]+1.0 ;
						}
				}		
       }

// if needed, push down to the next correlator
	if ( ((point+1)%2) == 0 )
		{
			if ( k+1 < Max_level)						// max number of correlator levels 
				{
					addcor( (acor[nf][k][point] + acor[nf][k][point-1] ) / 2., nf, k+1 ) ;
				}
		}
} 



void writecor() {

	std::ofstream outFile1("cor.dat");
	
	int i,k ;

	for ( i=0; i < p2-1; i++ )  			// first 7 dt=dt*1
		{
			if(ncor[0][0][i]>0) 
				{
					outFile1<<fcor[0][0][i]/ncor[0][0][i] << '\t';
				}
		}	

	for ( k=1; k < Max_level; k++ ) 
		{       
			for ( i=0; i < p2; i++ )
				{
					if( ncor[0][k][i] > 0)
						{
							outFile1<<fcor[0][k][i]/ncor[0][k][i]  << '\t'; 
						}
				}
		}
outFile1.close();	
}
 
int main() {

vctr3D new_data;
int data_lngth = 10000000 ;	      
		
zerocor();										// intialize the array with zeros
/*
	for (int  i=0; i < n_vars; i++ ) 
		{       
			for (int  j=0; j < Max_level ; j++ )
				{
					for (int  k=0; k < p2; k++ )
						{
							cout<< acor[i][j][k] <<'\t' ; 
							cout<< fcor[i][j][k] <<'\t' ; 
							cout<< ncor[i][j][k] <<'\t' ; 
						}
						
					for (int  k=p2; k < pcor; k++ )
						{
							cout<< acor[i][j][k] <<'\t' ; 
						}	
										
					cout<< pointcor[i][j] <<'\t' ; 
				}
		}

double test[3][3] = {-2.0};

	for (int  k=0; k < 3; k++ ) 
		{       
			for (int  i=0; i < 3; i++ )
				{
					cout << test[k][i]<<'\n'; 
					test[k][i] = {-2.0} ;
					cout << test[k][i]<<'\n'; 
				}
		}

*/           
std::string fileName="com.dat";
std::ifstream dataFile;
dataFile.open(fileName);
if(!dataFile.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
		std::string line;
  	
		for (int i=0;i<data_lngth;i++) 			// read data file line by line
			{
				std::getline(dataFile,line);
				std::istringstream currentLine(line);    
				currentLine >> new_data.comp[0];
				currentLine >> new_data.comp[1];
				currentLine >> new_data.comp[2];   	
				
				// call autocorr function
				addcor(new_data.comp[0],0,1);	// format (new_data,data_type,correlator_level) 		
		
			}
			
   	}
		
 writecor(); 
 	
}
