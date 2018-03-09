//  Calculate the force between particles i and j
if (particle[i].cluster!=particle[j].cluster)
	{


// dr=particle[i].pos-(particle[j].pos+dR+shift*xxshift);

dr=particle[i].pos-particle[j].pos;
dr.comp[0] -= ( round( dr.comp[1] * rbox.comp[1] ) * (*DEL_BOX)  ) ;						
dr.PBC(box,rbox);
r2=dr.norm2();

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
		/*	else {				
									
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
			} */
			
		}

} /* else {

if (r2<(r_cut2)) 
		{

			 if (r2<(rs2)) {
						    	r=std::sqrt(r2);								
								*p_energy+= phis - phicut + Fs*(rs-r);
							} 
			else {									
				r2inv=1.0/r2;
        		r6inv 	= 	r2inv*r2inv*r2inv;
        		r12inv 	= 	r6inv*r6inv;
				*p_energy+=4.0*epsilon*(sigma12*r12inv-sigma6*r6inv) - phicut;

			} 
			
		}
} */
