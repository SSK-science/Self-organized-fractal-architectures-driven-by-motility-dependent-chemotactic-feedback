#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>

#define Nx 100  //Chemical grid in x 
#define Ny 100	//Chemical grid in y
#define Np 1500		//Total # of particles
#define pi 3.14159265359
#define ux 0.5 
#define uy 0.86602540378

void Chem_neighbours(int P[Nx*Ny][6]);
void Sublattice(int Cx[2*Nx][Ny][2],int Cy[2*Nx][Ny][2]);
void Chemical_directions(double C_dirx[6],double C_diry[6]);
int Hex_grid(double x, double y,int Cx[2*Nx][Ny][2],int Cy[2*Nx][Ny][2]);
void init ( double * x, double * y, double * xw, double * yw,double * xold, double * yold, int L ,double g[Ny*Nx],int n[Np][6],int P[Nx*Ny][6],int Cx[2*Nx][Ny][2],int Cy[2*Nx][Ny][2]);
void init_conf( double * x, double * y, double * xw, double * yw,double g[Ny*Nx],int n[Np][6],int P[Nx*Ny][6],int Cx[2*Nx][Ny][2],int Cy[2*Nx][Ny][2], double decay_rate,double deposit_rate);
double Dynamics(double * x, double * y, double * xw, double * yw,double * xold,double * yold,int n[Np][6],double g[Ny*Nx], int P[Nx*Ny][6],int Cx[2*Nx][Ny][2],int Cy[2*Nx][Ny][2],double C_dirx[6],double C_diry[6],double persistence,double deposit_rate,double decay_rate,double stddev,double dr,double rc2);
double total_e ( double * rx, double * ry, int i, int L, double rc2,int Ns[Np],int Nid[Np][Np]); 
void verlet_list(double * rx, double * ry, int L, double rs,int Ns[Np],int Nid[Np][Np]);
int clustering(double *x,double *y, int time, double deposit_rate, double decay_rate);
double gaussian(double mean, double stddev);
double Pr(double x);
void decay(double decay_rate,double g[Ny*Nx],int G);
double PBC(double x,double L);

int Ns[Np], Nid[Np][Np]={}; double rs=8;
int T,i,j,t,ts;//T=Total time
int P[Nx*Ny][6], Cx[2*Nx][Ny][2],Cy[2*Nx][Ny][2],n[Np][6]; double C_dirx[6], C_diry[6];
double * x, * y,* xw, * yw,* xold, * yold;
double deposit_rate,decay_rate,persistence,stddev,rc2,dr,incr;
//Initialize with 0 chemicals 
double g[Ny*Nx] = {},G; 

int main(int argc, char *argv[]){
	
	srand48(time(0));
	//defining chemical neighbouring sites
	Chem_neighbours(P); 
	Sublattice(Cx,Cy);
		
	// Allocate the position arrays 
	x = (double*)malloc(Np*sizeof(double));	y = (double*)malloc(Np* sizeof(double));
	xw = (double*)malloc(Np*sizeof(double));	yw = (double*)malloc(Np* sizeof(double));
	xold = (double*)malloc(Np*sizeof(double));  yold = (double*)malloc(Np*sizeof(double));

	//Parameters
	deposit_rate=5.00; decay_rate=atof(argv[1]);
	incr=0.05; persistence=0.00; stddev= pi/4.0; rc2=pow(2.0,1.0/6.0); dr=1.0; ts =(int)(1.0/decay_rate); T= 100000;

	char filename[50]; 
	FILE *fw, *fwc; 
	sprintf(filename,"Pos_hex_%d_%2.4f_%2.4f.d",Np,deposit_rate,decay_rate);
	fw = fopen(filename,"a");
	sprintf(filename,"Chem_hex_%d_%2.4f_%2.4f.d",Np,deposit_rate,decay_rate);
	fwc = fopen(filename,"a");
		
	init(x,y,xw,yw,xold,yold,Nx,g,n,P,Cx,Cy);	//init_conf(x,y,xw,yw,g,n,P,Cx,Cy,decay_rate,deposit_rate); //To start from a given configuration
	verlet_list(x,y,Nx,rs,Ns,Nid);

	//Time evolution
	for(t=0;t<T+1;t++){//t is time
		printf("time = %d\tdeposit= %2.4f\tdecay= %2.4f\t Np= %d\n",t,deposit_rate,decay_rate,Np );
		//Montecarlo step at time t
		Dynamics(x,y,xw,yw,xold,yold,n,g,P,Cx,Cy,C_dirx,C_diry,persistence,deposit_rate,decay_rate,stddev,dr,rc2);
		
		//Stochastic decay of chemical
		decay(decay_rate,g,G);

		//Measurement
		if(t%(1000)==0){	
			//Chemical environment
			for(i=0;i<Nx;i++){
				for (j=0;j<Ny;j++) {
					fprintf(fwc,"%d %lf %lf %f \n",t,i+ux*(j%2),j*uy,deposit_rate*g[j*Nx+i]);
				}
			}
			// Particle positions
			for(i=0;i<Np;i++){
				fprintf(fw,"%d\t%d\t%lf\t%lf\t%lf\t%lf\n",t,i,x[i],y[i],xw[i],yw[i]);
			}
			
			//Clustering 
			//j=clustering(x,y,t,deposit_rate,decay_rate); 
			//fprintf(fw,"%d\t%d\n",t,j);
			
		}
		if(t%(10)==0) verlet_list(x,y,Nx,rs,Ns,Nid);
		
	}//time loop over
	//Chemical environment
	/*for(i=0;i<Nx;i++){
		for (j=0;j<Ny;j++) {
			fprintf(fwc,"%d %lf %lf %f \n",t,i+ux*(j%2),j*uy,deposit_rate*g[j*Nx+i]);
		}
	}*/	
	fclose(fw);	fclose(fwc);
	
return 0;
}

void Chem_neighbours(int P[Nx*Ny][6]){
	int i,j;
	for(i=0; i<Nx; i++){
		for(j=0; j<Ny; j++){
			P[i+j*Nx][0] =  i+j%2 + (j+1)*Nx ;
			P[i+j*Nx][1] =  i-(j+1)%2 + (j+1)*Nx ;
			P[i+j*Nx][2] =  i-1 + j*Nx ;
			P[i+j*Nx][3] =  i-(j+1)%2 + (j-1)*Nx;
			P[i+j*Nx][4] =  i+j%2 + (j-1)*Nx ;
			P[i+j*Nx][5] =  i+1 + j*Nx ;
		}
	} 
		//Periodic Boundary neighbours
	for(i=0; i<=Nx-1; i++){
		P[i][3] =  i-(Ny+1)%2 + (Ny-1)*Nx;				P[i][4] = i+Ny%2 + (Ny-1)*Nx ;	 	
		P[i+(Ny-1)*Nx][0] =  i+1 ;	P[i+(Ny-1)*Nx][1] =  i ;
	}
	for(j=0; j<=Ny-1; j++){
		P[j*Nx][1]= (j%2)*(0+(j+1)*Nx) + ((j+1)%2)*(Nx-1 + (j+1)*Nx) ; 
		P[j*Nx][2]=  Nx-1 + j*Nx ; 
		P[j*Nx][3]= (j%2)*(0 + (j-1)*Nx) + (Nx-1 + (j-1)*Nx)*((j+1)%2)  ;
		P[j*Nx + Nx-1][0]= (j%2)*(j+1)*Nx + ((j+1)%2)*(Nx-1 + (j+1)*Nx);
		P[j*Nx+ Nx-1][5]=  j*Nx ;  
		P[j*Nx+ Nx-1][4]= (j%2)*(j-1)*Nx + ((j+1)%2)*(Nx-1 + (j-1)*Nx);
	}
	P[0][3] = Nx-1+ (Ny-1)*Nx; 	P[Nx-1][4] = (Ny-1)*Nx+ Nx-1; 	P[Nx-1+(Ny-1)*Nx][0] = 0; P[(Ny-1)*Nx][1] = 0;

}
void Chemical_directions(double C_dirx[6],double C_diry[6]){
	
	C_dirx[0]= ux; C_diry[0]= uy; C_dirx[1]= -ux; C_diry[1]= uy;	 
	C_dirx[2]= -1; C_diry[2]= 0; C_dirx[3]= -ux; C_diry[3]= -uy;
	C_dirx[4]= ux; C_diry[4]= -uy; C_dirx[5]= 1;  C_diry[5]= 0; 
}

// Defining the sub-lattice
void Sublattice(int Cx[2*Nx][Ny][2],int Cy[2*Nx][Ny][2]){
	int i,j; 

	for(j=0; j<Ny; j++){
		for(i=0; i<2*Nx; i++){
			Cx[i][j][0]= i; 	Cy[i][j][0]= ((j+1)%2)*(j+ i%2) + (j%2)*(j+(i+1)%2)  ;
			Cx[i][j][1]= i+1; 	Cy[i][j][1]= ((j+1)%2)*(j+ (i+1)%2) + (j%2)*(j+i%2);
		}
	}
	for(j=0;j<Ny;j++){Cx[2*Nx-1][j][1]= 0; }
	for(i=0;i<2*Nx;i++){Cy[i][Ny-1][0]= (Ny-1)*(i%2); Cy[i][Ny-1][1]= (Ny-1)*((i+1)%2);}
	
}
// Initialize particle positions by assigning them on a cubic grid, then scaling positions to achieve a given box size and thereby, volume and density 
void init ( double * x, double * y, double * xw, double * yw,double * xold, double * yold, int L ,double g[Ny*Nx],int n[Np][6],int P[Nx*Ny][6],int Cx[2*Nx][Ny][2],int Cy[2*Nx][Ny][2]) {
	int i,k,ix,iy,ij;
	int n3=2; L = L/1;
	// Find the lowest perfect square, n3, greater than or equal to the number of particles  
	while ((n3*n3)<Np) n3++;

	ix=iy=0;
	//Assign particle positions 
	for (i=0;i<Np;i++) {
		
		x[i] =  ((double)ix+0.5)*L/n3; 
		y[i] = 	((double)iy+0.5)*L*uy/n3;
		xw[i]=x[i];	yw[i]=y[i];

		ix++;
		if (ix==n3) {
			ix=0;
			iy++;
    		}
		if (iy==n3) {
			iy=0;
      		}
		
		xold[i] = x[i] - 1.*(drand48()-0.5);	yold[i] = y[i] - 1.*(drand48()-0.5);
		ij = Hex_grid(x[i],y[i],Cx,Cy); 
		for(k=0;k<6; k++){
			n[i][k] = P[ij][k]; g[n[i][k]] +=1;
		}
	}
	
}
//Predefined initial configuration
void init_conf( double * x, double * y, double * xw, double * yw,double g[Ny*Nx],int n[Np][6],int P[Nx*Ny][6],int Cx[2*Nx][Ny][2],int Cy[2*Nx][Ny][2], double decay_rate,double deposit_rate) {
	int i,j,ij,k,m,c,t,nc,np; float rx,ry,cc; 
	char filename[50]; 
	FILE *fr, *frc;
	printf("hi");
	sprintf(filename,"Configuration_%d_%2.4f_%2.4f.d",Np,deposit_rate,decay_rate);
	fr = fopen(filename,"r");
	sprintf(filename,"Chem_hex_%d_%3.4f_%2.4f.d",Np,deposit_rate,decay_rate);
	frc = fopen(filename,"r");
	printf("hello");
	for(i=0;i<Np; i++){
		fscanf(fr,"%d %d %f %f %d %d",&t,&m,&rx,&ry,&nc,&np);
		x[m]= rx; y[m]=ry;	xw[m]=x[m];	yw[m]=y[m];		
		ij = Hex_grid(x[m],y[m],Cx,Cy); 
		for(k=0;k<6; k++){
			n[m][k] = P[ij][k]; //g[n[m][k]] +=1;
		}
	}
	for(i=0;i<Nx*Ny; i++){
		fscanf(frc,"%d %f %f %f",&t,&rx,&ry,&cc);
		j=(int)(ry/uy); k=(int)(rx-ux*(j%2)); c=(int)(cc/deposit_rate);
		g[k+j*Nx] =c;		
	}
	fclose(fr); fclose(frc);
}


//Particles movement
double Dynamics(double * x, double * y, double * xw, double * yw,double * xold,double * yold,int n[Np][6],double g[Ny*Nx],int P[Nx*Ny][6],int Cx[2*Nx][Ny][2],int Cy[2*Nx][Ny][2],double C_dirx[6],double C_diry[6],double persistence,double deposit_rate,double decay_rate,double stddev,double dr,double rc2){
	int i,j,k,ij;
	double dx,dy,E_old,E_new,xo,yo,xow,yow,xdrn=0,ydrn=0,theta,theta0,temp=1.0,r=1.0; char boundary[]="PBC"; int res = strcmp(boundary,"PBC");
	double avgdx=0,avgdy=0;
	
	Chemical_directions(C_dirx,C_diry); 
	for(j=0; j<Np; j++){//1 time unit corresponds to Np particle movement on an average
		i=Np*drand48();
		//E_old = total_e(x,y,i,Nx,rc2,Ns,Nid);
		xo = x[i];	yo = y[i];	xow = xw[i];	yow = yw[i];

		//printf("x %lf %lf \n",x[i],y[i]);
		xdrn = 0.001*(drand48()-0.5) +persistence*(x[i] - xold[i]); 
		ydrn = 0.001*(drand48()-0.5) +persistence*(y[i] - yold[i]);	 	
		for(k=0; k<6; k++){
			xdrn += Pr(deposit_rate*g[n[i][k]])*C_dirx[k]; ydrn += Pr(deposit_rate*g[n[i][k]])*C_diry[k];
			//printf("k %d n %d g %d\n",k,n[i][k],g[n[i][k]]);
		}
		//printf("%f %f\n",xdrn,ydrn);
		theta0 = atan2(ydrn,xdrn) ;
		//printf("theta0 %lf\n",theta0*180/pi);
		//gaussian distribution(mean=theta0,std. dev = stddev)
		theta = theta0 + gaussian(0.0,stddev);	
		
		r = 1.0*drand48();
		dx = dr*r*cos(theta);	dy = dr*r*sin(theta);
		x[i] += dx;	y[i] += dy; 	xw[i] += dx;	yw[i] += dy;	//Particle displaced in the direction of highest concentration most probably
		x[i] =PBC(x[i],1.0*Nx);	y[i] =PBC(y[i],uy*Ny);
		
		E_new = total_e(x,y,i,Nx,rc2,Ns,Nid);
		if( drand48() < exp(-(E_new-E_old)*temp ) ){
			xold[i] = xo;	yold[i] = yo;
			ij = Hex_grid(x[i],y[i],Cx,Cy); //Add chemical after moving to a new aite
			//ij = Hex_grid(xo,yo,Cx,Cy); 	// Add chemical before leaving a site
			for(k=0;k<6; k++){
				n[i][k] = P[ij][k]; g[n[i][k]] +=1; //Adding chemical to the accepted hexagonal grid
			}
			
		}
		else{
			x[i] = xo;	y[i] = yo;	xw[i] = xow;	yw[i] = yow;			
		}
		
		/*ij = Hex_grid(x[i],y[i],Cx,Cy); 
		for(k=0;k<6; k++){
			n[i][k] = P[ij][k]; g[n[i][k]] +=1; //Adding chemical to its hexagonal grid
		}*/

	}//Np particle movement at a given time done
}

//verlet list
void verlet_list(double * rx, double * ry, int L, double rs,int Ns[Np],int Nid[Np][Np]){
	int i,j;
	double dx, dy, r2, r6i,Lx=L,Ly=L*uy;
	double hLx=Lx/2.0,hLy=Ly/2.0;
	rs = rs*rs;
   	for(i=0;i<Np;i++){
		Ns[i]=0;
		for (j=0;j<Np;j++) {
			dx  = (rx[i]-rx[j]);
			dy  = (ry[i]-ry[j]);
			
			if (dx>hLx)       dx-=Lx;
			else if (dx<-hLx) dx+=Lx;
			if (dy>hLy)       dy-=Ly;
			else if (dy<-hLy) dy+=Ly;
			
			r2 = dx*dx + dy*dy ;
			if (r2<rs && j!=i ) {
				Nid[i][Ns[i]] = j;	Ns[i] += 1;
				
			}
		}
	}
}
// An algorithm for computing the total energy. 
double total_e ( double * rx, double * ry, int i, int L, double rc2,int Ns[Np],int Nid[Np][Np]) {
	int j;
	double dx, dy, r2, r6i,Lx=L,Ly=L*uy;
	double e = 0.0, hLx=Lx/2.0,hLy=Ly/2.0;
	rc2 = rc2*rc2;
   
	for (j=0;j<Ns[i];j++) {
		dx  = (rx[i]-rx[ Nid[i][j] ] );
		dy  = (ry[i]-ry[ Nid[i][j] ]);
		if (dx>hLx)       dx-=Lx;
		else if (dx<-hLx) dx+=Lx;
		if (dy>hLy)       dy-=Ly;
		else if (dy<-hLy) dy+=Ly;
		
		r2 = dx*dx + dy*dy ;
		if (r2<rc2 && j!=i ) {
			r6i = 1.0/(r2*r2*r2);
			e += 4*(r6i*r6i - r6i + 0.25) ;
		}
	}
   	return e;
}
//Clustering algorithm 
int clustering(double *x,double *y, int time, double deposit_rate, double decay_rate){
	
	FILE *fw, *fwc; char filename[50];
	sprintf(filename,"cluster_hex_%d_%2.4f_%2.4f.d",Np,deposit_rate,decay_rate);
	fw = fopen(filename,"a");
	sprintf(filename,"Ns_hex_%d_%2.4f_%2.4f.d",Np,deposit_rate,decay_rate);
	fwc = fopen(filename,"a");

	int N=Np-1,j,i,k,kc,l,r,s,t,lmxo,lmxn;double rij,dx,dy,rc=4,Lb=Nx,hL = Lb/2.0;
	int NC[Np]={}, C[2][Np]={};

	for(i=0; i<=N; i++){	C[0][i] = i; }
	C[1][0] = C[0][N]; N=N-1;k=1;kc=1; lmxo=0;lmxn=1;s=1; j=C[1][0]; 
	fprintf(fw,"%d\t%d\t%2.5f\t%2.5f\t%d\t%d\n",time,j,x[j],y[j],kc,lmxn);  // N=Number of particles -1 
	
	while(N>=0){  //while no particle remained for clustering
		while(s>0){ // while no particle in a given cluster k remained to check
			s=0; 
			for(l=lmxo;l<lmxn; l++){ //for newly cluseterd particles to check with unclustered ones
				r=0;
				while(r<N+1 ){ // while checking is done for a given particle 
					i=C[k][l]; j= C[0][r];
					//printf("\nk= %d\tl= %d\tr= %d\ti= %d\tj= %d, x= %2.4f,%2.4f\tN= %d",k,l,r,i,j,x[i],x[j],N);
					dx  = x[i]-x[j]; dy  = y[i]-y[j];
					if (dx>hL)       dx-=Lb;
					else if (dx<-hL) dx+=Lb;
					if (dy>hL*uy)       dy-=Lb*uy;
					else if (dy<-hL*uy) dy+=Lb*uy;
					rij = dx*dx + dy*dy ; 
					if(rij < rc){
						s+=1; C[k][lmxn+s-1] = j; C[0][r]=C[0][N]; N=N-1; //if clstrd add j to kth clstr nd bring the last unclstred to the first   
						fprintf(fw,"%d\t%d\t%2.5f\t%2.5f\t%d\t%d\n",time,j,x[j],y[j],kc,lmxn+s); 
					}	
					else r = r+1;// if not go to the next unclusterd one  
				}
			}	
			lmxo = lmxn; lmxn = lmxn +s;//update for the newly added particles  
		}
		NC[kc] = lmxn;//printf("\n%d\t%d\n",NC[k],k);
		fprintf(fwc,"%d\t%d\t%d\n",time,kc,NC[kc]);		
 
		if(N>-1){
			kc = kc+1; C[k][0] = C[0][N]; N=N-1;lmxo=0;lmxn=1;s=1;NC[kc] = lmxn;
			j=C[k][0]; fprintf(fw,"%d\t%d\t%2.5f\t%2.5f\t%d\t%d\n",time,j,x[j],y[j],kc,lmxn); 
			
		} //Form the next clstr and repeat the process 
	}
	//printf("\n%d\n",kc);
	fprintf(fw,"\n"); fclose(fw);
	fclose(fwc);
	return kc;
}
//Finding out the hexagonal lattice sites
int Hex_grid(double x, double y,int Cx[2*Nx][Ny][2],int Cy[2*Nx][Ny][2]){
	int ix,iy;
	double x1,y1,x2,y2,d1,d2;

	ix = (int)(x/ux); iy=(int)(y/uy);
	x1=ux*Cx[ix][iy][0];x2=ux*Cx[ix][iy][1]; y1=uy*Cy[ix][iy][0];y2=uy*Cy[ix][iy][1]; 
	d1 = (x-x1)*(x-x1) + (y-y1)*(y-y1); d2 = (x-x2)*(x-x2) + (y-y2)*(y-y2);
	
	if(d1<d2){  	iy=(int)(y1/uy); ix = (int)(x1-ux*(iy%2)); 	}
	else {  	iy=(int)(y2/uy); ix = (int)(x2-ux*(iy%2)); 	}
	
	return iy*Nx+ix;	
}
//Pr(x)=Probability of following a chemical having concentration x
double Pr(double x){  
double y;
y = (1.0)*x/(1.0+x);
return y;
}

double gaussian(double mean, double stddev){
double u1,u2,x; 
u1 = sqrt( -2.0*log( drand48() ) ); u2 = 2*pi*drand48();
x = mean + u1*cos(u2)*stddev;
return x;
}
//Chemical decay
void decay(double decay_rate,double g[Ny*Nx],int G){
	int i,j,k; double q;
	G=0;
	for(i=0; i<Nx; i++){
		for(j=0; j<Ny; j++){
			//Decay 			
			k = g[j*Nx+i] ; 
			while(k>0){
				q = (int)(1.0 - decay_rate + drand48());
				g[j*Nx+i] = g[j*Nx+i] - 1*(1-q); 
				k = k-1;
			}
			G += g[j*Nx+i];
		}
	}
}

//Periodic boundary condition
double PBC(double x,double L){
if(x>L) x = x - L;
else if(x<0) x = x + L;
return x;
}
