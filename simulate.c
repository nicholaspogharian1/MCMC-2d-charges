#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include "energies.c"
#include "gaussian.c"

/*
Command line arguments:
1: value of e field in x
2: length of box in x
3: length of box in y
4: run number
*/

//compile: gcc -o MCSIM simulate.c -lm -std=c99

void main(int argc, char* argv[]) {
	
	srand( (unsigned)time(NULL) + atoi(argv[4]) ) ; //set random seed to cpu time plus the argv1 value (so parallel jobs have different random seeds)
	
	FILE *file0;
	char adr[100];

	///define important simulation variables
	double pi=acos(-1);
	//int n_steps=20000000; //number of MC moves to make
	int n_steps=400000000; //number of MC moves to make
	int between_steps=1000000; //number of steps between recorded steps
	int N=40; //number of particles
	double vp=2; //valence of cations
	double vn=-2; //valence of anions
	int Np=20; //number of cations
	int Nn=20; //number of anions
	double sp=1; //diameter of cations
	double sn=1.2; //diameter of anions
	double mean_step_size=0.1; //mean MC step size to take
	double step_size_sd=mean_step_size/10; //standard deviation of MC step size
	double cutoff=pow(2.,1./6.); //cutoff multiplier for WCA interactions
	double epsilon=1.;
	double invTstar=2.14; //strength of electrostatic interactions in pseudo 2d (defined as 1 / T* where T*=0.12) (this parameter is actually equal to 1/T*/Z^2 since we multiply the potential by ch[i]*ch[j])
	double xsi=(80.*1.9)/(2.*2.); //parameter defining range of pseudo 2d electrostatic interactions
	double theta,step_size,testx,testy,old_field_energy,old_energy,new_energy,move_prob; //parameters to update during MC steps
	int moving_particle; //particle we are currently moving, also to be updated during MC steps

	double efield=atof(argv[1]); //magnitude of applied e field (kT/e/d)
	double Lx=atof(argv[2]); //length of box in x (aperiodic)
	double Ly=atof(argv[3]); //length of simulation cell in y (periodic)
	
	int successful_moves=0; //keep track of the number of MC moves which actually work

	///set up the particles and their relevant properties
	double x[N],y[N],rad[N],ch[N],pair_energies[N][N],new_pair_energies[N][N],field_energies[N],dx,dy,dist;
	memset(pair_energies,0.,sizeof(pair_energies));
	memset(new_pair_energies,0.,sizeof(new_pair_energies));
	
	//particle charges and radii
	for (int i=0;i<Np;i++) {
		ch[i]=vp;
		rad[i]=sp/2.;
	}
	for (int i=Np;i<N;i++) {
		ch[i]=vn;
		rad[i]=sn/2.;
	}
	
	//reset energies file
	sprintf(adr,"./efield_%s_lx_%s_ly_%s/energies_%s.dat",argv[1],argv[2],argv[3],argv[4]);
	file0=fopen(adr,"w");
	fclose(file0);
	
	//particle positions
	//cations can go anywhere in the system, anions are limited in x: 0<x<240
	for (int i=0;i<N;i++) {
		checkpoint: ;
		if (i<Np) {
			x[i]=(double)rand()/RAND_MAX*(Lx-2*rad[i])+rad[i];
		}
		else {
			x[i]=(double)rand()/RAND_MAX*(240.-2*rad[i])+rad[i]; //negative charges can go between x=0 and x=240 (due to size constraints in the real system)
		}
		y[i]=(double)rand()/RAND_MAX*Ly;
		//check to make sure there are no overlapping particles
		for (int j=0;j<i;j++) {
			dx=x[i]-x[j];
			dy=y[i]-y[j];
			dy=dy-round(dy/Ly)*Ly;
			dist=sqrt(dx*dx+dy*dy);
			if (dist<(rad[i]+rad[j])) {
				goto checkpoint; //reposition the particle if it is too close to any of the other particles
			}
		}
	}

	//calculate initial system energies and total energy
	calc_pair_energies(N,x,y,ch,rad,Ly,cutoff,epsilon,invTstar,xsi,pair_energies);
	calc_field_energies(N,x,ch,efield,field_energies);
	copy_data(N,pair_energies,new_pair_energies);

	///make monte carlo moves
	for (int s=0;s<n_steps;s++) {
		bool selected_movepos=false;
		while (!selected_movepos) {
			step_size=gaussran()*step_size_sd+mean_step_size;
			theta=(double)rand()/RAND_MAX*2.*pi;
			moving_particle=rand()%N;
			testx=x[moving_particle]+step_size*cos(theta);
			testy=y[moving_particle]+step_size*sin(theta);
			testy=math_mod(testy,Ly);
			//see if the moved particle remains in the bounds of the system
			if (moving_particle<Np) {
				if (testx>rad[moving_particle]&testx<Lx-rad[moving_particle]) {
					selected_movepos=true;
				}
			}
			else {
				if (testx>rad[moving_particle]&testx<240.-rad[moving_particle]) {
					selected_movepos=true;
				}
			}
		}
		//printf("theta: %lf\n",theta);
		/*
		printf("s: %d\ndata before update\n",s);
		printf("moving particle: %d\n",moving_particle);
		printarray_1d(N,field_energies);
		printf("\n");
		printarray_2d(N,pair_energies);
		printf("\n");
		printarray_2d(N,new_pair_energies);
		printf("\n\n");
		*/
		old_field_energy=field_energies[moving_particle];
		old_energy=upper_triangular_sum(N,pair_energies)+old_field_energy;
		update_pair_energies(N,moving_particle,testx,testy,x,y,ch,rad,Ly,cutoff,epsilon,invTstar,xsi,new_pair_energies);
		update_field_energy(N,moving_particle,testx,ch,efield,field_energies);
		new_energy=upper_triangular_sum(N,new_pair_energies)+field_energies[moving_particle];
		/*
		printf("data mid update\n");
		printarray_1d(N,field_energies);
		printf("\n");
		printarray_2d(N,pair_energies);
		printf("\n");
		printarray_2d(N,new_pair_energies);
		printf("\n\n");
		*/
		if (new_energy<old_energy) { //move the particle if the free energy goes down
			//printf("%d\t%lf\t%lf\n",s,upper_triangular_sum(N,pair_energies)+old_field_energy+sum(N,field_energies)-field_energies[moving_particle],upper_triangular_sum(N,new_pair_energies)+sum(N,field_energies));
			x[moving_particle]=testx;
			y[moving_particle]=testy;
			copy_data(N,new_pair_energies,pair_energies);
			successful_moves+=1;
		}
		else {
			move_prob=exp(old_energy-new_energy);
			//printf("%d\t%lf\t%lf\t%lf\n",s,upper_triangular_sum(N,pair_energies)+old_field_energy+sum(N,field_energies)-field_energies[moving_particle],upper_triangular_sum(N,new_pair_energies)+sum(N,field_energies),move_prob);
			if ((double)rand()/RAND_MAX<=move_prob) { //move the particle if the random number is within the acceptance window
				//printf("move accepted\n");
				x[moving_particle]=testx;
				y[moving_particle]=testy;
				copy_data(N,new_pair_energies,pair_energies);
				successful_moves+=1;
			}
			else { //otherwise, do not move the particle
				//printf("move not accepted\n");
				field_energies[moving_particle]=old_field_energy;
			}
		}
		/*
		printf("data post update\n");
		printarray_1d(N,field_energies);
		printf("\n");
		printarray_2d(N,pair_energies);
		printf("\n");
		printarray_2d(N,new_pair_energies);
		printf("\n\n");
		*/
		//print out some results every so often (energy,positions)
		//sample every 40000 steps
		if (s%between_steps==0) {
			//write positions to file
			
			sprintf(adr,"./efield_%s_lx_%s_ly_%s/run_%s/positions_%d.dat",argv[1],argv[2],argv[3],argv[4],(int)(s/between_steps));
			file0=fopen(adr,"w");
			for (int i=0;i<N;i++) {
				fprintf(file0,"%lf\t%lf\t%d\n",x[i],y[i],ch[i]);
			}
			fclose(file0);
			
			//write energies to a file
			sprintf(adr,"./efield_%s_lx_%s_ly_%s/energies_%s.dat",argv[1],argv[2],argv[3],argv[4]);
			file0=fopen(adr,"a");
			fprintf(file0,"%lf\n",upper_triangular_sum(N,pair_energies)+sum(N,field_energies));
			//printf("%lf\n",upper_triangular_sum(N,pair_energies)+sum(N,field_energies));
			fclose(file0);
			//write info to a file
			/*
			sprintf(adr,"./efield_%s_lx_%s_ly_%s/info_%s.dat",argv[1],argv[2],argv[3],argv[4]);
			file0=fopen(adr,"a");
			fprintf(file0,"%lf\t%lf\t%lf\n",new_energy-old_energy,exp(old_energy-new_energy),upper_triangular_sum(N,pair_energies)+field_energies[moving_particle]-old_energy);
			fclose(file0);
			*/
		}
	}
	
	//print the move success rate to a file at the end
	sprintf(adr,"./efield_%s_lx_%s_ly_%s/success_rate_%s.dat",argv[1],argv[2],argv[3],argv[4]);
	file0=fopen(adr,"w");
	fprintf(file0,"%lf\n",(double)successful_moves/n_steps);
	//printf("%lf\n",(double)successful_moves/n_steps);
	fclose(file0);
	
	//print run info to a file
	sprintf(adr,"./efield_%s_lx_%s_ly_%s/run_%s/info",argv[1],argv[2],argv[3],argv[4]);
	file0=fopen(adr,"w");
	fprintf(file0,"%lf %lf %d %d %d\n",Lx,Ly,Np,Nn,n_steps/between_steps);
	fclose(file0);
}