#include<stdio.h>
#include<math.h>

//compute short range energy of the system (pseudo 2d electrostatic and short range)
//this version computes energies as an NxN array where element i,j is the interaction energy between particles i and j
//only the upper triangle of the matrix has energy values
void calc_pair_energies(int N, double x[N], double y[N], double ch[N], double rad[N], double Ly, double cutoff, double epsilon, double invTstar, double xsi, double energies[N][N]) {
	double sigma; //sigma value of WCA interaction
	double lj_frac; //sigma divided by dist
	double dx,dy,dist;
	
	for (int i=0;i<N;i++) {
		//calculate pair interaction energies
		for (int j=i+1;j<N;j++) {
			energies[i][j]=0.;
			//calculate particle distances from one another (assuming 2d system with periodicity in y only)
			dx=x[i]-x[j];
			dy=y[i]-y[j];
			dy=dy-round(dy/Ly)*Ly;
			dist=sqrt(dx*dx+dy*dy);
				
			//calculate WCA energies between particles
			sigma=rad[i]+rad[j];
			if (dist<cutoff*sigma) {
				lj_frac = sigma/dist;
				energies[i][j]+=4.*epsilon*(pow(lj_frac,12)-pow(lj_frac,6));
			}
				
			//calculate short ranged electrostatic energies between particles
			energies[i][j]-=ch[i]*ch[j]*invTstar*log(dist/(dist+xsi));
		}
	}
}

//update short range energy of the system (pseudo 2d electrostatic and short range)
//this version computes energies as an NxN array where element i,j is the interaction energy between particles i and j
//only the upper triangle of the matrix has energy values
void update_pair_energies(int N, int num, double testx, double testy, double x[N], double y[N], double ch[N], double rad[N], double Ly, double cutoff, double epsilon, double invTstar, double xsi, double energies[N][N]) {
	double sigma; //sigma value of WCA interaction
	double lj_frac; //sigma divided by dist
	double dx,dy,dist;
	
	for (int i=0;i<N;i++) {
		//calculate pair interaction energies
		for (int j=i+1;j<N;j++) {
			if (i==num) {
				energies[i][j]=0.;
				//calculate particle distances from one another (assuming 2d system with periodicity in y only)
				dx=testx-x[j];
				dy=testy-y[j];
				dy=dy-round(dy/Ly)*Ly;
				dist=sqrt(dx*dx+dy*dy);
					
				//calculate WCA energies between particles
				sigma=rad[i]+rad[j];
				if (dist<cutoff*sigma) {
					lj_frac = sigma/dist;
					energies[i][j]+=4.*epsilon*(pow(lj_frac,12)-pow(lj_frac,6));
				}
					
				//calculate short ranged electrostatic energies between particles
				energies[i][j]-=ch[i]*ch[j]*invTstar*log(dist/(dist+xsi));
			}
			if (j==num) {
				energies[i][j]=0.;
				//calculate particle distances from one another (assuming 2d system with periodicity in y only)
				dx=x[i]-testx;
				dy=y[i]-testy;
				dy=dy-round(dy/Ly)*Ly;
				dist=sqrt(dx*dx+dy*dy);
					
				//calculate WCA energies between particles
				sigma=rad[i]+rad[j];
				if (dist<cutoff*sigma) {
					lj_frac = sigma/dist;
					energies[i][j]+=4.*epsilon*(pow(lj_frac,12)-pow(lj_frac,6));
				}
					
				//calculate short ranged electrostatic energies between particles
				energies[i][j]-=ch[i]*ch[j]*invTstar*log(dist/(dist+xsi));
			}
		}
	}
}

//this function copies data from arr1 to arr2
void copy_data(int N, double arr1[N][N], double arr2[N][N]) {
	for (int i=0;i<N;i++) {
		for (int j=i+1;j<N;j++) {
			arr2[i][j]=arr1[i][j];
		}
	}
}

//calculates energy from the applied field
void calc_field_energies(int N, double x[N], double ch[N], double efield, double energies[N]) {
	for (int i=0;i<N;i++) {
		energies[i]=-efield*ch[i]*x[i];
	}
}

//updates energy from the applied field
void update_field_energy(int N, int num, double testx, double ch[N], double efield, double energies[N]) {
	energies[num]=-efield*ch[num]*testx;
}

//sum of elements in double array of N elements
double sum(int N, double array[N]) {
	//returns the sum
	double sum=0.;
	for (int i=0;i<N;i++) {
		sum+=array[i];
	}
	return sum;
}

//sum of the upper triangle of elements in an NxN array of doubles
double upper_triangular_sum(int N, double array[N][N]) {
	//returns the sum
	double sum=0.;
	for (int i=0;i<N;i++) {
		for (int j=i+1;j<N;j++) {
			sum+=array[i][j];
		}
	}
	return sum;
}

//returns the mathematical (always positive) modulus
double math_mod(double a, double b) {
    double r = fmod(a,b);
    return (r < 0) ? r + b : r;
}

//prints NxN array to the terminal in a readable way
void printarray_2d(int N, double array[N][N]) {
	for (int i=0;i<N;i++) {
		for (int j=0;j<N;j++) {
			if (j==N-1) printf("%lf\n",array[i][j]);
			else printf("%lf\t",array[i][j]);
		}
	}
}

//prints 1D array of length N to the terminal in a readable way
void printarray_1d(int N, double array[N]) {
	for (int i=0;i<N;i++) {
		if (i==N-1) printf("%lf\n",array[i]);
		else printf("%lf\t",array[i]);
	}
}