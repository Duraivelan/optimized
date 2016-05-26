/* 
 * File:   force.h
 * Author: duraivelan
 *
 * Created on 12 December, 2014, 11:12 AM
 */

# include "rigid_force.h"
# include "defs.h"
#include <algorithm>
#include <functional>
#include <array>

static inline double sign(double a,double b) { return a= fabs(a),(b<0)?-a:a; }

using namespace std;

 struct RowSort {
        bool operator()(vector<int> a, vector<int>  b)
        {   
            return a[0] < b[0];
        }   
    } ;

void forceUpdate( vector<SubData>& particle,  double *p_energy, int* combine_now , int combine[][4], int* step) {

  int    NrCells[3],MaxNrCells;
  double scale[3];

	int    i,j,i_seg,j_seg;
	int    ii,jj;
	int    mi[3],m,mj[3];
	double F,r2, r;
	vctr3D dr,dR, Fij;
	double r2inv, r6inv, r12inv;

  for ( i = 0 ; i < 3 ; i++ )
  {
    NrCells[i] = floor ( box.comp[i]*Nsegm / (r_cut*apct_rt) ); // cellnr runs from 0 to NrCells-1
    scale  [i] = NrCells[i] * rbox.comp[i];
    if ( NrCells[i] < 3 ) { cout << "*** NrCells[" << i << "] = " << NrCells[i] << endl ; abort(); }
  }

// periodic boundary conditions

  MaxNrCells = max( max( NrCells[x], NrCells[y] ), NrCells[z]);
  int    periodN[ MaxNrCells + 2 ][3];
  double periodR[ MaxNrCells + 2 ][3];

  for ( j = 0 ; j < 3 ; j++ )
  {
    periodN[0][j] = NrCells[j] - 1;          // left neighbour of leftmost cell
    periodR[0][j] = -box.comp[j];       // correction to add to particle j in rij = ri - rj
    for ( i = 1 ;  i < NrCells[j] + 1 ; i++ )
    {
      periodN[i][j] = i - 1; // same cell
      periodR[i][j] = 0.;
    } // i
    periodN[NrCells[j] + 1][j] = 0;          // right neigbour of rightmost cell
    periodR[NrCells[j] + 1][j] = +box.comp[j];
  } // j
  

// generate grid list
 	int grid[NrCells[x]][NrCells[y]][NrCells[z]][MaxPerCell+1];
 	int seg_grid[NrCells[x]][NrCells[y]][NrCells[z]][MaxPerCell+1];

  for ( mi[x] = 0 ; mi[x] < NrCells[x] ; mi[x]++ )
  {
    for ( mi[y] = 0 ; mi[y] < NrCells[y] ; mi[y]++ )
    {
      for ( mi[z] = 0 ; mi[z] < NrCells[z] ; mi[z]++ )
      {

        grid[mi[x]][mi[y]][mi[z]][0] = 0;
        seg_grid[mi[x]][mi[y]][mi[z]][0] = 0;
      } // miz
    } // miy
  } // mix

for ( int i = 0 ; i < NrParticles ; i ++ )
  {

vctr3D r_segm ; 
	for (int j =0 ; j < Nsegm ; j++ )
		{

            r_segm = particle[i].pos + particle[i].dir* lh * ( (double(j) + 0.5) * NsegmINV - 0.5 ) ;
			
			r_segm.PBC(box,rbox);
			
			mi[x] = int ( (r_segm.comp[x]+havbox.comp[x]) * scale[x] );
			mi[y] = int ( (r_segm.comp[y]+havbox.comp[y]) * scale[y] );
			mi[z] = int ( (r_segm.comp[z]+havbox.comp[z]) * scale[z] );       
		
			if ( int (grid[mi[x]][mi[y]][mi[z]][0]) >= MaxPerCell-1 )
				{
					cout << "*** cell overfull" << endl;
					cout << mi[x] << "  " << mi[y] << "  " << mi[z] << endl;
					abort();
				}

			grid[mi[x]][mi[y]][mi[z]][0] ++ ;
			//  cout << i << "  " << mix << "  " << miy << "  " << miz << "  " << grid[mix][miy][miz][0] << endl;
			grid[mi[x]][mi[y]][mi[z]][ int (grid[mi[x]][mi[y]][mi[z]][0])] = i;
			seg_grid[mi[x]][mi[y]][mi[z]][ int (grid[mi[x]][mi[y]][mi[z]][0])] = j;
		}
} // i
	//		 cout<<"exit grid "<<endl;

// zero forces

for ( int i = 0 ; i < NrParticles ; i ++ )
  {
	particle[i].frc=null3D;
}

// calculate energy and forces

  for ( mi[x] = 0 ; mi[x] < NrCells[x] ; mi[x]++ )
  {
    for ( mi[y] = 0 ; mi[y] < NrCells[y] ; mi[y]++ )
    {
      for ( mi[z] = 0 ; mi[z] < NrCells[z] ; mi[z]++ )
      {
        for ( ii = 1 ; ii <= grid[mi[x]][mi[y]][mi[z]][0] ; ii++ )
        {
          i = grid[mi[x]][mi[y]][mi[z]][ii];
          i_seg = seg_grid[mi[x]][mi[y]][mi[z]][ii];

          // particle j in same cell as i
          dR = null3D;
          for ( jj = ii + 1 ; jj <= grid[mi[x]][mi[y]][mi[z]][0] ; jj++ )
          {
			j = grid[mi[x]][mi[y]][mi[z]][jj];
			j_seg = seg_grid[mi[x]][mi[y]][mi[z]][jj];
		//	if (particle[i].cluster!=particle[j].cluster)
			//	{
					#include "pairforce_rod.h"
			//	}
          } // jj

          // particle j in neighbour cell to i
          for ( m = 0 ; m < 13 ; m++ )
          {
            mj[x]      = periodN[ mi[x] + dm[m][x] + 1 ][x];
            mj[y]      = periodN[ mi[y] + dm[m][y] + 1 ][y];
            mj[z]      = periodN[ mi[z] + dm[m][z] + 1 ][z];
            dR.comp[x] = periodR[ mi[x] + dm[m][x] + 1 ][x];
            dR.comp[y] = periodR[ mi[y] + dm[m][y] + 1 ][y];
            dR.comp[z] = periodR[ mi[z] + dm[m][z] + 1 ][z];
            for ( jj = 1 ; jj <= grid[mj[x]][mj[y]][mj[z]][0] ; jj++ )
            {
				j = grid[mj[x]][mj[y]][mj[z]][jj];
				j_seg = seg_grid[mj[x]][mj[y]][mj[z]][jj];
		//	if (particle[i].cluster!=particle[j].cluster)
			//	{
					#include "pairforce_rod.h"
			//	}
            } // jj
          } // m
        } // ii
      } // miz
    } // miy
  } // mix

}



