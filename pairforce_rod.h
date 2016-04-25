//  Calculate the force between particles i and j
if (particle[i].cluster!=particle[j].cluster)
	{
		
/* Revision of Carlos Vega & Santiago Lago Computers Chem. 18, 55-59, 1994

 Subrutine to evaluate the shortest distance between two rods of
 different length

 The original code did not give the symmetry property of the distance for almost parallel rods.
 The coordinates of the centers of the rods should be given in a periodic system

 r1,r2: centers of rods
 w1,w2: unit orientation vectors of rods
 lh1,lh2: halves of the length of rods
 Lv.x,Lv.y,Lv.z the edges of the periodic simulation cell
*/

// Minimum distance in the periodic system:
	
//dr=particle[i].pos-(particle[j].pos+dR);
    dr = particle[i].pos - particle[j].pos ;
    dr.PBC(box,rbox);
//---------------- Distance of two rods: -------------------------------------
  double lamda1, lamda2, lamdai, lamdaj, 
  r2= dr.norm2(),
  rw1= dr*particle[i].dir,
  rw2= dr*particle[j].dir,
  uu = particle[i].dir*particle[j].dir,
  uu2 = 1.0/(1.0-(uu)*(uu));
  vctr3D dt; 

            lamda1 = ( -rw1 + rw2*uu ) * uu2 ;
            lamda2 = (  rw2 - rw1*uu ) * uu2 ;
            dt = dr + particle[i].dir*lamda1 - particle[j].dir*lamda2 ;
            r2 =  dt.norm2();
            if ((fabs(lamda1) < lh && fabs(lamda2) < lh ) ||  r2 > DIAMpairlistSQ ) 
			{
				r2+=0.0;
				}
			else {

            lamda1 = sign(1.0,lamda1) * min( lh,fabs(lamda1) ) ;
            lamda2 = sign(1.0,lamda2) * min( lh,fabs(lamda2) ) ;

            lamdaj = particle[j].dir *( dr + particle[i].dir*lamda1 ); 

            if (fabs(lamdaj) < lh ) {
               dt = dr + particle[i].dir *lamda1 - particle[j].dir*lamdaj ; 
               r2 =  dt.norm2();
		   } 
            else
            {
               lamdai = particle[i].dir *(particle[j].dir*lamda2 - dr ); 
               if (fabs(lamdai) < lh ) {
				  dt = dr + particle[i].dir *lamdai - particle[j].dir*lamda2 ;                   
                  r2 =  dt.norm2();
			  }
               else {
                  lamdai = sign(1.0,lamdai) * min( lh,fabs(lamdai) ) ;
                  lamdaj = sign(1.0,lamdaj) * min( lh,fabs(lamdaj) ) ;
                  dt = dr + particle[i].dir *lamdai - particle[j].dir*lamdaj ; 
                  r2 =  dt.norm2();
			  }
			}
		}
/*
 double  xla,xmu,
  r2= dr.norm2(),
  rw1= dr*particle[i].dir,
  rw2= dr*particle[j].dir,
  w1w2= particle[i].dir*particle[j].dir,
  cc= 1.0-(w1w2)*(w1w2);

// Checking whether the rods are or not parallel:
// The original code is modified to have symmetry:

 if(cc<1.0e-6) {
  if(rw1 && rw2) {
   xla= rw1/2.0;
   xmu= -rw2/2.0;
  }
 }

 else {

// Step 1

  xla= (rw1-w1w2*rw2)/cc;
  xmu= (-rw2+w1w2*rw1)/cc;
 }

// Step 2

if( fabs(xla)>lh || fabs(xmu)>lh ) {

// Step 3 - 7

  if((fabs(xla)-lh) > (fabs(xmu)-lh)) {
   xla= sign(lh,xla);
   xmu= xla*w1w2-rw2;
   if( fabs(xmu)>lh ) xmu= sign(lh,xmu);
  }
  else {
   xmu= sign(lh,xmu);
   xla= xmu*w1w2+rw1;
   if( fabs(xla)>lh ) xla= sign(lh,xla);
  }
 }

// Step 8

 r2+= (xla)*(xla)+(xmu)*(xmu) + 2.0 *(xmu*rw2 -xla*(rw1+xmu*w1w2));	
 */
 	 
		r2-=r_min2*xxshift;

		if (r2<(r_cut2)) 
		{
		//	cout<<r2<<endl;
			if(xxclustering && r2 < (r_min2)) 
				{
					*combine_now += 1;
					combine[	*combine_now	]	[ 0 ] = min(particle[i].cluster,particle[j].cluster);
					combine[	*combine_now	]	[ 1 ] = max(particle[i].cluster,particle[j].cluster);
					combine[	*combine_now	]	[ 2 ] = min(i,j);
					combine[	*combine_now	]	[ 3 ] = max(i,j);
				//	cout<<combine[*combine_now][0]<<'\t'<<combine[*combine_now][1]<<"insdide pairforce"<<'\t'<<*step<<endl;
				}

			 if (r2<(rs2)) {
						    	r=std::sqrt(r2);								
												
	// 							exponential potential from PHYSICAL REVIEW E VOLUME 50, NUMBER 3 SEPTEMBER 1994 Browniian dynamics simulations of self-difFusion and shear viscosity of near-hard-sphere colloids
	//							F = n*epsilon*pow(sigma*sigma/rs,n)/rs; //  n is the exponent of the potential function
            					Fij=dr*Fs;
								particle[i].frc+=Fij;
								particle[j].frc-=Fij;

							//	if(cell[neighborList[i][j][k][m][0]][neighborList[i][j][k][m][1]][neighborList[i][j][k][m][2]][lC2] > cell[i][j][k][lC1]) {
								*p_energy+= phis - phicut + Fs*(rs-r)*(rs-r);
							} 
			else {				
									
				r2inv=1.0/r2;
        		r6inv 	= 	r2inv*r2inv*r2inv;
        		r12inv 	= 	r6inv*r6inv;
	//			simple potential
	//			F=2*epsilon*(sigma-r)*rx/r;
	//			*p_energy+=2*epsilon*(sigma-r)*(sigma-r);
	
	//			exponential potential from PHYSICAL REVIEW E VOLUME 50, NUMBER 3 SEPTEMBER 1994 Browniian dynamics simulations of self-difFusion and shear viscosity of near-hard-sphere colloids
	//			F = n*epsilon*pow(sigma*sigma/r2,n/2)/r2; //  n is the exponent of the potential function
	//			F = 4*epsilon*(12*pow(sigma*sigma/r2,6)-6*pow(sigma*sigma/r2,3))/r2; 	WCA potential
				F = 4.0*epsilon*(12.0*sigma12*r12inv-6.0*sigma6*r6inv)*r2inv;
        		Fij=dr*F;
				particle[i].frc+=Fij;
				particle[j].frc-=Fij;
				*p_energy+=4.0*epsilon*(sigma12*r12inv-sigma6*r6inv) - phicut;
			} 
			
		}

}  
