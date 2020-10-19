#include <stdlib.h>
#include <stdio.h>
#include <time.h>
//#include <mex.h>
#include <math.h>
//#include "matrix.h"

#define delta 1e-10


/* 
Revision History
  
  First Version available on October 10, 2009 

  A runnable version on October 15, 2009

  Major revision on October 29, 2009
  (Some functions appearing in a previous version have deleted, please refer to the previous version for the old functions.
   Some new functions have been added as well)

*/

/*

Files contained in this header file sfa.h:

1. Algorithms for solving the linear system A A^T z0 = Av (see the description of A from the following context)

  int PTRANS1(double *zMax, double *z0, double *Av, int nn)

  void Thomas(double *zMax, double *z0, 
              double * Av, int nn)

  void Rose(double *zMax, double *z0, 
            double * Av, int nn)

  int supportSet(double *x, double *v, double *z, 
                 double *g, int * S, double lambda, int nn, int d)

  void dualityGap(double *gap, double *z, 
                  double *g, double *s, double *Av, 
				  double lambda, int nn)

  void dualityGap2(double *gap, double *z, 
                  double *g, double *s, double *Av, 
				  double lambda, int nn)


2. The Subgraident Finding Algorithm (SFA) for solving problem (4) (refer to the description of the problem for detail) 
  
  int sfa(double *x,     double *gap,
		 double *z,     double *z0,   double * v,   double * Av, 
		 double lambda, int nn,       int maxStep,
		 double *s,     double *g,
		 double tol,    int tau,       int flag)

  int sfa_special(double *x,     double *gap,
		 double *z,     double * v,   double * Av, 
		 double lambda, int nn,       int maxStep,
		 double *s,     double *g,
		 double tol,    int tau)

  int sfa_one(double *x,     double *gap,
		 double *z,     double * v,   double * Av, 
		 double lambda, int nn,       int maxStep,
		 double *s,     double *g,
		 double tol,    int tau)


*/


/*

  Some mathematical background.

  In this file, we discuss how to solve the following subproblem,

        min_x  1/2 \|x-v\|^2  + lambda \|A x\|_1,                 (1)

  which is a key problem used in the Fused Lasso Signal Approximator (FLSA).

  Also, note that, FLSA is a building block for solving the optimation problmes with fused Lasso penalty.
  
  In (1), x and v are n-dimensional vectors, 
        and A is a matrix with size (n-1) x n, and is defined as follows (e.g., n=4):
		                                         A= [ -1  1  0  0;
												       0  -1 1  0;
													   0  0  -1 1]

  The above problem can be reformulated as the following equivalent min-max optimization problem

        min_x  max_z  1/2 \|x-v\|^2  + <A x, z>
		subject to   \|z\|_{infty} \leq lambda                     (2)


  It is easy to get that, at the optimal point

                         x = v - AT z,                             (3)

  where z is the optimal solution to the following optimization problem

        min_z  1/2  z^T A AT z - < z, A v>,
		subject to  \|z\|_{infty} \leq lambda                      (4)


  
  Let B=A A^T. It is easy to get that B is a (n-1) x (n-1) tridiagonal matrix.
  When n=5, B is defined as:
                                                B= [ 2  -1   0    0;
												     -1  2   -1   0;
													 0  -1   2    -1;
													 0   0   -1   2]

  Let z0 be the solution to the linear system:

                               A A^T * z0 = A * v                  (5)

  The problem (5) can be solve by the Thomas Algorithm, in about 5n multiplications and 4n additions.

  It can also be solved by the Rose's Algorithm, in about 2n multiplications and 2n additions.

  Moreover, considering the special structure of the matrix A (and B), 
       it can be solved in about n multiplications and 3n additions

  If lambda \geq \|z0\|_{infty}, x_i= mean(v), for all i, 
                the problem (1) admits near analytical solution


  We have also added the restart technique, please refer to our paper for detail!

*/


/*
///////////////    Solving the Linear System via PTRANS-1 Algorithm (Askar and Karawia 2015) \\\\\\\\\\\\\\\\\\
*/

int PTRANS1(double *zMax, double *z0, double *Av, int nn)
{

/*

	We apply the PTRANS-1 algorithm for solving the following linear system
	                   B * z0 = Av


    B is a banded pentadiagonal matrix

    B= [ 6  -4   1   0   0
        -4   6  -4   1   0
         1  -4   6  -4   1
         0   1  -4   6  -4
         0   0   1  -4   6 ]

    z0 is the result,  Av is unchanged after the computation

    z0 is an nn dimensional vector
    
	*/

    int i;
    int j;
    double tt, z_max;

    /* INPUT */
    double *a = (double *)malloc(sizeof(double)*nn);
    double *b = (double *)malloc(sizeof(double)*nn);
    double *c = (double *)malloc(sizeof(double)*nn);
    double *d = (double *)malloc(sizeof(double)*nn);
    double *e = (double *)malloc(sizeof(double)*nn);

    for(i = 0; i < nn; i++) {
        a[i] = -4.0;
        b[i] = 1.0;
        c[i] = -4.0;    
        d[i] = 6.0;
        e[i] = 1.0;
    }
    a[nn-1] = 0.0;
    b[nn-2] = 0.0;    
    b[nn-1] = 0.0;
    c[0] = 0.0;
    e[0] = 0.0;
    e[1] = 0.0;

    double *alpha = (double *)malloc(sizeof(double)*(nn-1));
    double *beta = (double *)malloc(sizeof(double)*(nn-2));
    double *z = (double *)malloc(sizeof(double)*nn);
    double *gamma = (double *)malloc(sizeof(double)*nn);    
    double *mu = (double *)malloc(sizeof(double)*nn);
    
    /* Step 3 */
    mu[0] = d[0];
    if(mu[0] == 0.0) {
        free(a);
        free(b);
        free(c);
        free(d);
        free(e);
        free(alpha);
        free(beta);
        free(z);
        free(gamma);
        free(mu);
        return -1;
    }
    alpha[0] = a[0] / mu[0];
    beta[0] = b[0] / mu[0];
    z[0] = Av[0] / mu[0];
    
    /* Step 4 */
    gamma[1] = c[1];
    mu[1] = d[1] - alpha[0] * gamma[1];
    if(mu[1] == 0.0) {
        free(a);
        free(b);
        free(c);
        free(d);
        free(e);
        free(alpha);
        free(beta);
        free(z);
        free(gamma);
        free(mu);
        return -1;
    }
    alpha[1] = (a[1] - beta[0] * gamma[1]) / mu[1];
    beta[1] = b[1] / mu[1];
    z[1] = (Av[1] - z[0] * gamma[1]) / mu[1]; 

    /* Step 5 */
    for(j = 3; j <= nn-2; j++) {
        i = j-1;

        gamma[i] = c[i] - alpha[i-2] * e[i];
        mu[i] = d[i] - beta[i-2] * e[i] - alpha[i-1] * gamma[i];
        if(mu[i] == 0.0) {
            free(a);
            free(b);
            free(c);
            free(d);
            free(e);
            free(alpha);
            free(beta);
            free(z);
            free(gamma);
            free(mu);
            return -1;
        }
        alpha[i] = ( a[i] - beta[i-1] * gamma[i] ) / mu[i];
        beta[i] = b[i] / mu[i];
        z[i] = ( Av[i] - z[i-2] * e[i] - z[i-1] * gamma[i] ) / mu[i];  
    }   
    gamma[nn-2] = c[nn-2] - alpha[nn-4] * e[nn-2];
    mu[nn-2] = d[nn-2] - beta[nn-4] * e[nn-2] - alpha[nn-3] * gamma[nn-2];
    if(mu[nn-2] == 0.0) {
        free(a);
        free(b);
        free(c);
        free(d);
        free(e);
        free(alpha);
        free(beta);
        free(z);
        free(gamma);
        free(mu);
        return -1;
    }
    alpha[nn-2] = ( a[nn-2] - beta[nn-3] * gamma[nn-2] ) / mu[nn-2];
    gamma[nn-1] = c[nn-1] - alpha[nn-3] * e[nn-1];
    mu[nn-1] = d[nn-1] - beta[nn-3] * e[nn-1] - alpha[nn-2] * gamma[nn-1];
    if(mu[nn-1] == 0.0) {
        free(a);
        free(b);
        free(c);
        free(d);
        free(e);
        free(alpha);
        free(beta);
        free(z);
        free(gamma);
        free(mu);
        return -1;
    }
//    z[nn-2] = ( Av[nn-2] - z[nn-3] * e[nn-2] - z[nn-3] * gamma[nn-2] ) / mu[nn-2];        // Algorithm was incorrect in paper
//    z[nn-1] = ( Av[nn-1] - z[nn-2] * e[nn-1] - z[nn-2] * gamma[nn-1] ) / mu[nn-1];        // Algorithm was incorrect in paper
    z[nn-2] = ( Av[nn-2] - z[nn-4] * e[nn-2] - z[nn-3] * gamma[nn-2] ) / mu[nn-2];
    z[nn-1] = ( Av[nn-1] - z[nn-3] * e[nn-1] - z[nn-2] * gamma[nn-1] ) / mu[nn-1];

    /* Step 6 */
    z0[nn-1] = z[nn-1];
    z0[nn-2] = z[nn-2] - alpha[nn-2] * z0[nn-1];
    for(j = nn - 2; j >= 1; j--) {
        i = j-1;

        z0[i] = z[i] - alpha[i] * z0[i+1] - beta[i] * z0[i+2];
    }

    z_max= fabs(z0[0]);
    for(i = 1; i < nn; i++) {
        tt = fabs(z0[i]);        
        if(z_max < tt) {
            z_max = tt;
        }
    }
    *zMax=z_max;

    free(a);
    free(b);
    free(c);
    free(d);
    free(e);
    free(alpha);
    free(beta);
    free(z);
    free(gamma);
    free(mu);
    return 0;
}


/*
///////////////    Solving the Linear System via Thomas's Algorithm \\\\\\\\\\\\\\\\\\
*/

void Thomas(double *zMax, double *z0, double * Av, int nn){

	/*

	We apply the Tomas algorithm for solving the following linear system
	                   B * z0 = Av


    Thomas algorithm is also called the tridiagonal matrix algorithm

  B=[ 2  -1   0    0;
	  -1  2   -1   0;
	  0  -1   2    -1;
	  0   0   -1   2]

    z0 is the result,  Av is unchanged after the computation


    c is a precomputed nn dimensional vector
	c=[-1/2, -2/3, -3/4, -4/5, ..., -nn/(nn+1)]

    c[i]=- (i+1) / (i+2)
	c[i-1]=- i / (i+1)

    z0 is an nn dimensional vector
    
	*/

	int i;
	double tt, z_max;

	/*
	Modify the coefficients in Av (copy to z0)
	*/
	z0[0]=Av[0]/2;
	for (i=1;i < nn; i++){
		tt=Av[i] + z0[i-1];
		z0[i]=tt - tt / (i+2);
	}

	/*z0[i]=(Av[i] + z0[i-1]) * (i+1) / (i+2);*/
		
	/*z0[i]=(Av[i] + z0[i-1])/ ( 2 - i / (i+1));*/

	
	/*
	Back substitute (obtain the result in z0)
	*/
	z_max= fabs(z0[nn-1]);

	for (i=nn-2; i>=0; i--){

		z0[i]+=  z0[i+1] -  z0[i+1]/ (i+2);

		/*z0[i]+=  z0[i+1] * (i+1) / (i+2);*/

		tt=fabs(z0[i]);

		if (tt > z_max)
			z_max=tt;

	}
	*zMax=z_max;
	
}



			
/*
///////////////    Solving the Linear System via Rose's Algorithm \\\\\\\\\\\\\\\\\\
*/

void Rose(double *zMax, double *z0,	double * Av, int nn){

	/*
	We use the Rose algorithm for solving the following linear system
	                   B * z0 = Av


  B=[ 2  -1   0    0;
	  -1  2   -1   0;
	  0  -1   2    -1;
	  0   0   -1   2]

    z0 is the result,  Av is unchanged after the computation

    z0 is an nn dimensional vector
    
	*/

	int i, m;
	double s=0, z_max;


	/*
	We follow the style in CLAPACK
	*/
	m= nn % 5;
	if (m!=0){
		for (i=0;i<m; i++)
			s+=Av[i] * (i+1);
	}
	for(i=m;i<nn;i+=5)
		s+=   Av[i]   * (i+1) 
		    + Av[i+1] * (i+2) 
		    + Av[i+2] * (i+3) 
			+ Av[i+3] * (i+4) 
			+ Av[i+4] * (i+5);
	s/=(nn+1);


	/*
    from nn-1 to 0
	*/
	z0[nn-1]=Av[nn-1]- s;
	for (i=nn-2;i >=0; i--){
		z0[i]=Av[i] + z0[i+1];
	}

	/*
    from 0 to nn-1
	*/
	z_max= fabs(z0[0]);
	for (i=0; i<nn; i++){

		z0[i]+=  z0[i-1];

		s=fabs(z0[i]);

		if (s > z_max)
			z_max=s;

	}
	*zMax=z_max;
	
}



/*
////////////////    compute x for restarting \\\\\\\\\\\\\\\\\\\\\\\\\

x=omega(z)

v: the vector to be projected
z: the approximate solution
g: the gradient at z (g should be computed before calling this function

nn: the length of z, g, and S (maximal length for S)

n:  the length of x and v

S: records the indices of the elements in the support set
*/

int supportSet(double *x, double *v, double *z, double *g, int * S, double lambda, int nn, int d){

	int i, j, n=nn+d, numS=0;
	double temp;


	/*
	we first scan z and g to obtain the support set S
	*/

	/*numS: number of the elements in the support set S*/
	for(i=0;i<nn; i++){
		if ( ( (z[i]==lambda) && (g[i] < delta) ) || ( (z[i]==-lambda) && (g[i] >delta) )){
			S[numS]=i;
			numS++;
		}
	}
	
	/*
	printf("\n %d",numS);
	*/

	if (numS==0){ /*this shows that S is empty*/
		temp=0;
		for (i=0;i<n;i++)
			temp+=v[i];

		temp=temp/n;
		for(i=0;i<n;i++)
			x[i]=temp;

		return numS;
	}


    /*
	 Next, we deal with numS >=1
     */

	/*process the first block
	 
	   j=0
	*/
	temp=0;
	for (i=0;i<=S[0]; i++)
		temp+=v[i];
	/*temp =sum (v [0: s[0] ]*/
	temp=( temp + z[ S[0] ] ) / (S[0] +1);
	for (i=0;i<=S[0]; i++)
		x[i]=temp;


	/*process the middle blocks
	
	  If numS=1, it belongs the last block
	*/
	for (j=1; j < numS; j++){
		temp=0;
		for (i= S[j-1] +1; i<= S[j]; i++){
			temp+=v[i];
		}

		/*temp =sum (v [ S[j-1] +1: s[j] ]*/

		temp=(temp - z[ S[j-1] ] + z[ S[j] ])/ (S[j]- S[j-1]);              /*  #### NOTE HERE THAT THERE SEEMS TO BE THE FUSION PENALTY LCOATED HERE   (I DON'T THINK IT IS RELEVANT THOUGH)*/

		for (i= S[j-1] +1; i<= S[j]; i++){
			x[i]=temp;
		}
	}

	/*process the last block
	j=numS-1;
	*/
	temp=0;
	for (i=S[numS-1] +1 ;i< n; i++)
		temp+=v[i];
	/*temp =sum (v [  (S[numS-1] +1): (n-1) ]*/

	temp=( temp - z[ S[numS-1] ] ) / (nn - S[numS-1]); /*S[numS-1] <= nn-1*/

	for (i=S[numS-1] +1 ;i< n; i++)
		x[i]=temp;
	
	return numS;

}



/*

////////////  Computing the duality gap \\\\\\\\\\\\\\\\\\\\\\\\\\

we compute the duality corresponding the solution z

z: the approximate solution
g: the gradient at z (we recompute the gradient)
s: an auxiliary variable
Av: A*v

nn: the lenght for z, g, s, and Av

The variables g and s shall be revised.

The variables z and Av remain unchanged.
*/

//void dualityGap(double *gap, double *z, double *g, double *s, double *Av, double lambda, int nn){
void dualityGap(double *gap, double *z, double *g, double *s, double *Av, double lambda, int nn, int d){

	int i, m;
	double temp;

	if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
        g[0]=z[0] + z[0] - z[1] - Av[0];                              
        for (i=1;i<nn-1;i++){
            g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
        }
        g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];        
    }
    else if(d == 2) {   /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
        g[0] = 6.0 * z[0] - 4.0 * z[1] + z[2] - Av[0];                  
        g[1] = -4.0 * z[0] + 6.0 * z[1] - 4.0 * z[2] + z[3] - Av[1];		
        for (i=2;i<nn-2;i++){
            g[i] = z[i-2] - 4.0 * z[i-1] + 6.0 * z[i] - 4.0 * z[i+1] + z[i+2] - Av[i];              
        }
        g[nn-2] = z[nn-4] - 4.0 * z[nn-3] + 6.0 * z[nn-2] - 4.0 * z[nn-1] - Av[nn-2]; 
        g[nn-1] = z[nn-3] - 4.0 * z[nn-2] + 6.0 * z[nn-1] - Av[nn-1]; 
    }
    else {
        exit(-1);
    }
	
	for (i=0;i<nn;i++)
		if (g[i]>0)
			s[i]=lambda + z[i];
		else
			s[i]=-lambda + z[i];

		
	temp=0;					
	m=nn%5;
	
	if (m!=0){
		for(i=0;i<m;i++)
			temp+=s[i]*g[i];
	}
	
	for(i=m;i<nn;i+=5)
		temp=temp + s[i]  *g[i]
		          + s[i+1]*g[i+1]
		          + s[i+2]*g[i+2]
		          + s[i+3]*g[i+3]
		          + s[i+4]*g[i+4];
	*gap=temp;
}


/*
Similar to dualityGap,

  The difference is that, we assume that g has been computed.
*/

void dualityGap2(double *gap, double *z, double *g, double *s, double *Av, double lambda, int nn){

	int i, m;
	double temp;


	/*
	g[0]=z[0] + z[0] - z[1] - Av[0];
	for (i=1;i<nn-1;i++){
		g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
	}	
	g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];



     #### IN THEORY, IT SHOULD BE BASED ON SOMETHING LIKE THIS
    if(d == 1) { 
        g[0]=z[0] + z[0] - z[1] - Av[0];                              
        for (i=1;i<nn-1;i++){
            g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
        }
        g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];        
    }
    else if(d == 2) {
        g[0] = 6.0 * z[0] - 4.0 * z[1] + z[2] - Av[0];                  
        g[1] = -4.0 * z[0] + 6.0 * z[1] - 4.0 * z[2] + z[3] - Av[1];		
        for (i=2;i<nn-2;i++){
            g[i] = z[i-2] - 4.0 * z[i-1] + 6.0 * z[i] - 4.0 * z[i+1] + z[i+2] - Av[i];              
        }
        g[nn-2] = z[nn-4] - 4.0 * z[nn-3] + 6.0 * z[nn-2] - 4.0 * z[nn-1] - Av[nn-2]; 
        g[nn-1] = z[nn-3] - 4.0 * z[nn-2] + 6.0 * z[nn-1] - Av[nn-1]; 
    }
    else {
        exit(-1);
    }

    */
	
	for (i=0;i<nn;i++)
		if (g[i]>0)
			s[i]=lambda + z[i];
		else
			s[i]=-lambda + z[i];

		
	temp=0;					
	m=nn%5;
	
	if (m!=0){
		for(i=0;i<m;i++)
			temp+=s[i]*g[i];
	}
	
	for(i=m;i<nn;i+=5)
		temp=temp + s[i]  *g[i]
		          + s[i+1]*g[i+1]
		          + s[i+2]*g[i+2]
		          + s[i+3]*g[i+3]
		          + s[i+4]*g[i+4];
	*gap=temp;
}


/*
generateSolution:

  generate the solution x based on the information of z and g 
  (!!!!we assume that g has been computed as the gradient of z!!!!)

*/

/*
int generateSolution(double *x, double *z, double *gap,
					  double *v, double *Av,
					  double *g, double *s, int *S,
					  double lambda, int nn){
*/
int generateSolution(double *x, double *z, double *gap,
					  double *v, double *Av,
					  double *g, double *s, int *S,
					  double lambda, int nn, int d, int gamma1, int gamma2){

	int i, m, numS, n=nn+d;
	double temp, funVal1, funVal2;
	               
	/*
	z is the appropriate solution,
	and g contains its gradient
	*/


   /*
		We assume that n>=3, and thus nn>=2
		
		  We have two ways for recovering x. 
		  The first way is x = v - A^T z
		  The second way is x =omega(z)
  */

	temp=0;
    m=nn%5;
    if (m!=0){
        for (i=0;i<m;i++)
            temp+=z[i]*(g[i] + Av[i]);
    }
    for (i=m;i<nn;i+=5)
        temp=temp + z[i]  *(g[i]   + Av[i])
                  + z[i+1]*(g[i+1] + Av[i+1])
                  + z[i+2]*(g[i+2] + Av[i+2])
                  + z[i+3]*(g[i+3] + Av[i+3])
                  + z[i+4]*(g[i+4] + Av[i+4]);
    funVal1=temp /2;
    /* funVal1 = SUM_{i=0, 1, ..., nn-1} z[i](g[i] + Av[i]) */

    
    temp=0;
    m=nn%5;
    if(gamma1 == 1) {
        if (m!=0){
            for (i=0;i<m;i++)
                temp+=fabs(g[i]);
        }
        for (i=m;i<nn;i+=5)
            temp=temp + fabs(g[i])
            + fabs(g[i+1])
            + fabs(g[i+2])
            + fabs(g[i+3])
            + fabs(g[i+4]);
        funVal1=funVal1+ temp*lambda;
        /* funVal1 = (1/2) * SUM_{i=0,1,...,nn-1} z[i](g[i] + Av[i]) + LAMBDA * SUM_{i=0,1,...,nn-1} |g[i]| */
    }
    else if(gamma1 == 2) {
        if (m!=0){
            for (i=0;i<m;i++)
                temp+=g[i]*g[i];
        }
        for (i=m;i<nn;i+=5)
            temp=temp + g[i]*g[i]
            + g[i+1]*g[i+1]
            + g[i+2]*g[i+2]
            + g[i+3]*g[i+3]
            + g[i+4]*g[i+4];
        funVal1=funVal1+ temp*lambda;
        /* funVal1 = (1/2) * SUM_{i=0,1,...,nn-1} z[i](g[i] + Av[i]) + LAMBDA * SUM_{i=0,1,...,nn-1} g[i]^2 */
    }
    else {
        printf("ERROR: gamma1 not equal to 1 or 2\n");
        exit(-1);
    }
    
    
    /*
           we compute the solution by the second way
    */

    numS= supportSet(x, v, z, g, S, lambda, nn, d);
    
	/*
        we compute the objective function of x computed in the second way
    */
    
    temp=0;
    m=n%5;
    if (m!=0){
        for (i=0;i<m;i++)
            temp+=(x[i]-v[i]) * (x[i]-v[i]);
    }
    for (i=m;i<n;i+=5)
        temp=temp + (x[i]-  v[i]) * (  x[i]-  v[i])
                  + (x[i+1]-v[i+1]) * (x[i+1]-v[i+1])
                  + (x[i+2]-v[i+2]) * (x[i+2]-v[i+2])
                  + (x[i+3]-v[i+3]) * (x[i+3]-v[i+3])
                  + (x[i+4]-v[i+4]) * (x[i+4]-v[i+4]);
    funVal2=temp/2;
    /* funVal2 = (1/2) * SUM_{i=0,1,...,n-1} (x[i]-v[i])^2 */    

    if(d == 1) {                       /*  THIS MAY NOT BE CORRECT. BE CAREFUL AS IT IS POSSIBLE THAT THIS WAS ACTUALLY A DIFFERENT PENALTY PART. NEED TO CHECK CAREFULLY.*/
        if(gamma2 == 1) {        
            temp=0;
            m=nn%5;
            if (m!=0){
                for (i=0;i<m;i++)
                    temp+=fabs( x[i+1]-x[i] );                  /*  #### NOTE HERE THAT THERE SEEMS TO BE THE FUSION PENALTY LCOATED HERE   */
            }
            for (i=m;i<nn;i+=5)
                temp=temp + fabs( x[i+1]-x[i] )                 /*  #### NOTE HERE THAT THERE SEEMS TO BE THE FUSION PENALTY LCOATED HERE   */
                          + fabs( x[i+2]-x[i+1] )               /*  #### NOTE HERE THAT THERE SEEMS TO BE THE FUSION PENALTY LCOATED HERE   */
                          + fabs( x[i+3]-x[i+2] )               /*  #### NOTE HERE THAT THERE SEEMS TO BE THE FUSION PENALTY LCOATED HERE   */
                          + fabs( x[i+4]-x[i+3] )               /*  #### NOTE HERE THAT THERE SEEMS TO BE THE FUSION PENALTY LCOATED HERE   */
                          + fabs( x[i+5]-x[i+4] );              /*  #### NOTE HERE THAT THERE SEEMS TO BE THE FUSION PENALTY LCOATED HERE   */
            funVal2=funVal2 + lambda * temp;
            /* funVal2 = (1/2) * SUM_{i=0,1,...,n-1} (x[i]-v[i])^2 + LAMBDA * SUM_{i=0,1,...,nn-1} |x[i+1]-x[i]| */
        }
        else if(gamma2 == 2) {        
            temp=0;
            m=nn%5;
            if (m!=0){
                for (i=0;i<m;i++)
                    temp+=( x[i+1]-x[i] )*( x[i+1]-x[i] );                  /*  #### NOTE HERE THAT THERE SEEMS TO BE THE FUSION PENALTY LCOATED HERE   */
            }
            for (i=m;i<nn;i+=5)
                temp=temp + ( x[i+1]-x[i] )*( x[i+1]-x[i] )                 /*  #### NOTE HERE THAT THERE SEEMS TO BE THE FUSION PENALTY LCOATED HERE   */
                          + ( x[i+2]-x[i+1] )*( x[i+2]-x[i+1] )               /*  #### NOTE HERE THAT THERE SEEMS TO BE THE FUSION PENALTY LCOATED HERE   */
                          + ( x[i+3]-x[i+2] )*( x[i+3]-x[i+2] )               /*  #### NOTE HERE THAT THERE SEEMS TO BE THE FUSION PENALTY LCOATED HERE   */
                          + ( x[i+4]-x[i+3] )*( x[i+4]-x[i+3] )               /*  #### NOTE HERE THAT THERE SEEMS TO BE THE FUSION PENALTY LCOATED HERE   */
                          + ( x[i+5]-x[i+4] )*( x[i+5]-x[i+4] );              /*  #### NOTE HERE THAT THERE SEEMS TO BE THE FUSION PENALTY LCOATED HERE   */
            funVal2=funVal2 + lambda * temp;
            /* funVal2 = (1/2) * SUM_{i=0,1,...,n-1} (x[i]-v[i])^2 + LAMBDA * SUM_{i=0,1,...,nn-1} (x[i+1]-x[i])^2 */
        } 
        else {
            printf("ERROR: gamma2 not equal to 1 or 2\n");
            exit(-1);
        }
    }   
    else if(d == 2) {                       /*  THIS MAY NOT BE CORRECT. BE CAREFUL AS IT IS POSSIBLE THAT THIS WAS ACTUALLY A DIFFERENT PENALTY PART. NEED TO CHECK CAREFULLY.*/
        if(gamma2 == 1) {        
            temp=0;
            for (i = 0; i < nn-2; i++) {
                temp = temp + fabs( x[i] - 2.0*x[i+1] + x[i+2] );
            }
            funVal2=funVal2 + lambda * temp;
            /* funVal2 = (1/2) * SUM_{i=0,1,...,n-1} (x[i]-v[i])^2 + LAMBDA * SUM_{i=0,1,...,nn-2} |x[i] - 2*x[i+1] + x[i+2]| */ 
        }
        else if(gamma2 == 2) {        
            temp=0;
            for (i = 0; i < nn-2; i++) {
                temp = temp + ( x[i] - 2.0*x[i+1] + x[i+2] )*( x[i] - 2.0*x[i+1] + x[i+2] );
            }
            funVal2=funVal2 + lambda * temp;
            /* funVal2 = (1/2) * SUM_{i=0,1,...,n-1} (x[i]-v[i])^2 + LAMBDA * SUM_{i=0,1,...,nn-2} (x[i] - 2*x[i+1] + x[i+2])^2 */ 
        }
        else {
            printf("ERROR: gamma2 not equal to 1 or 2\n");
            exit(-1);
        }
    }      
    else {
        exit(-1);
    } 


    /*
	printf("\n    funVal1=%e, funVal2=%e, diff=%e\n", funVal1, funVal2, funVal1-funVal2);
	*/


  
    
    if (funVal2 > funVal1){  /*
                                  we compute the solution by the first way
                              */
        /* NOTE TO SELF. THIS IS x = v - A^T * z (z is p-d dimensional vector), x and v are p dimensional vectors. */
        if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
            x[0]=v[0] + z[0];
	        for(i=1;i<n-1;i++)
		        x[i]= v[i] - z[i-1] + z[i];
	        x[n-1]=v[n-1] - z[n-2];       
        }
        else if(d == 2) {   /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
	        x[0] = v[0] - z[0];
            x[1] = v[1] + 2.0 * z[0] - z[1];		    
            for(i = 2; i < n-2; i++) {
		        x[i] = v[i] - z[i-2] + 2.0 * z[i-1] - z[i];
            }
	        x[n-2] = v[n-2] - z[n-4] + 2.0 * z[n-3];
            x[n-1] = v[n-1] - z[n-2];     
        }
        else {
            exit(-1);
        }


    }
    else{
        
        /*
        the solution x is computed in the second way
        the gap can be further reduced
        (note that, there might be numerical error)
         */
        
        *gap=*gap - (funVal1- funVal2);
        if (*gap <0)
            *gap=0;
	}

	return (numS);
}


/*
void restartMapping(double *g, double *z,  double * v, 
		 double lambda, int nn)
*/
void restartMapping(double *g, double *z,  double * v, 
		 double lambda, int nn, int d)
{

	int i, n=nn+d;
	double temp;
	int* S=(int *) malloc(sizeof(int)*nn);
	double *x=(double *)malloc(sizeof(double)*n);
	double *s=(double *)malloc(sizeof(double)*nn);
	double *Av=(double *)malloc(sizeof(double)*nn);
	int numS=-1;    
    int status;
     
	/*
	for a given input z, 
	we compute the z0 after restarting

    The elements in z lie in [-lambda, lambda]

    The returned value is g
	*/


        if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
            for (i=0;i<nn; i++) {
        		Av[i]=v[i+1]-v[i];      
            }
        }
        else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
            for (i = 0; i < nn; i++) {
        		Av[i] = v[i] - 2.0 * v[i+1] + v[i+2];      
            }
        } 
        else {
            exit(-1);
        }


	
		if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
            g[0]=z[0] + z[0] - z[1] - Av[0];                              
            for (i=1;i<nn-1;i++){
	            g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
            }
            g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];        
        }
        else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
            g[0] = 6.0 * z[0] - 4.0 * z[1] + z[2] - Av[0];                  
            g[1] = -4.0 * z[0] + 6.0 * z[1] - 4.0 * z[2] + z[3] - Av[1];		
            for (i=2;i<nn-2;i++){
	            g[i] = z[i-2] - 4.0 * z[i-1] + 6.0 * z[i] - 4.0 * z[i+1] + z[i+2] - Av[i];              
            }
            g[nn-2] = z[nn-4] - 4.0 * z[nn-3] + 6.0 * z[nn-2] - 4.0 * z[nn-1] - Av[nn-2]; 
            g[nn-1] = z[nn-3] - 4.0 * z[nn-2] + 6.0 * z[nn-1] - Av[nn-1];        
        }
        else {
            exit(-1);
        }


		numS = supportSet(x, v, z, g, S, lambda, nn, d);
		
		
		/*With x, we compute z via
		AA^T z = Av - Ax
		 */
	
		/*
		compute s= Av -Ax
		*/
		
        if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
            for (i=0;i<nn; i++) {
        		s[i]=Av[i] - x[i+1] + x[i];      
            }
        }
        else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
            for (i=0;i < nn; i++) {
                s[i]=Av[i] - (x[i] - 2.0 * x[i+1] + x[i+2]);
            }
        }
        else {
            exit(-1);
        }				
	

        if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
            /*
		    Apply Rose Algorithm for solving z
		    */
					
		    Thomas(&temp, g, s, nn); 					
		    /*
		    Rose(&temp, g, s, nn);
		    */      
        }
        else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
            /*
		    Apply PTRANS-1 Algorithm from Askar and Karawia (20150 for solving z
		    */
					
		    status = PTRANS1(&temp, g, s, nn);
            if(status == -1) {
                printf("ERROR: unable to solve system because mu[i]=0 with PTRANS-1 algorithm\n");
            }      
        }
        else {
            exit(-1);
        }
	

	   /*
	   project g to [-lambda, lambda]
	   */
	
		for(i=0;i<nn;i++){		
			if (g[i]>lambda)
				g[i]=lambda;
			else
				if (g[i]<-lambda)
					g[i]=-lambda;
		}

	
		free (S);
		free (x);
		free (s);
		free (Av);

}



	/*

     /////////////////////////////////////// Explanation for the function sfa \\\\\\\\\\\\\\\\\\\\\\\\\\\\

     Our objective is to solve the fused Lasso signal approximator (flsa) problem:

             min_x  g(x) 1/2 \|x-v\|^2  + lambda \|A x\|_1,                      (1)

     Let x* be the solution (which is unique), it satisfies
	    
		              0 in  x* - v +  A^T * lambda *SGN(Ax*)                     (2)

     To solve x*, it suffices to find
	                
					  y*  in A^T * lambda *SGN(Ax*)                              (3)
	 that satisfies
	    
		              x* - v + y* =0                                             (4)
	 which leads to
	                  x*= v - y*                                                 (5)

     Due to the uniqueness of x*, we conclude that y* is unique. 
	 
	 As y* is a subgradient of lambda \|A x*\|_1, 
	         we name our method as Subgradient Finding Algorithm (sfa).

     y* in (3) can be further written as
	                  
					  y*= A^T * z*                                               (6)
	 where

	                  z* in lambda* SGN (Ax*)                                    (7)

     From (6), we have
	                  z* = (A A^T)^{-1} A * y*                                   (8)

     Therefore, from the uqniueness of y*, we conclude that z* is also unique.
	 Next, we discuss how to solve this unique z*.

     The problem (1) can reformulated as the following equivalent problem:	 
	   
		 min_x  max_z  f(x, z)= 1/2 \|x-v\|^2  + <A x, z>
		 subject to   \|z\|_{infty} \leq lambda                                  (9)

     At the saddle point, we have
              
				        x = v - AT z,                                            (10)

	 which somehow concides with (5) and (6)

     Plugging (10) into (9), we obtain the problem
	        
			  min_z  1/2  z^T A AT z - < z, A v>,
		      subject to  \|z\|_{infty} \leq lambda,                             (11)

     In this program, we apply the Nesterov's method for solving (11).


    Duality gap:
	
	At a given point z0, we compute x0= v - A^T z0.
	It is easy to show that
	                  min_x f(x, z0) = f(x0, z0) <= max_z f(x0, z)               (12)

    Moreover, we have
	                  max_z f(x0, z) - min_x f(x, z0) 
					     <= lambda * \|A x0\|_1 - < z0, Av - A A^T z0>           (13)

    It is also to get that
	     
		              f(x0, z0) <= f(x*, z*) <= max_z f(x0, z)                   (14)

					  g(x*)=f(x*, z*)                                            (15)

                      g(x0)=max_z f(x0, z)                                       (17)

    Therefore, we have

	                  g(x0)-g(x*) <= lambda * \|A x0\|_1 - < z0, Av - A A^T z0>  (18)


    We have applied a restarting technique, which is quite involved; and thus, we do not explain here.

     /////////////////////////////////////// Explanation for the function sfa \\\\\\\\\\\\\\\\\\\\\\\\\\\\
	*/


/*
////////////               sfa              \\\\\\\\\\\\\\\\\\\\\

For sfa, the stepsize of the Nesterov's method is fixed to 1/4, so that no line search is needed.


  
    Explanation of the parameters:
    
	Output parameters
	x:    the solution to the primal problem
	gap:  the duality gap (pointer)

    Input parameters
	z:    the solution to the dual problem (before calling this function, z contains a starting point)
               !!!!we assume that the starting point has been successfully initialized in z !!!!
	z0:   a variable used for multiple purposes:
			  1) the previous solution z0
			  2) the difference between z and z0, i.e., z0=z- z0

	lambda:   the regularization parameter (and the radius of the infity ball, see (11)).
	nn:       the length of z, z0, Av, g, and s
	maxStep:  the maximal number of iterations

	v:    the point to be projected (not changed after the program)
	Av:   A*v (not changed after the program)

	s:        the search point (used for multiple purposes)
	g:        the gradient at g (and it is also used for multiple purposes)

	tol:      the tolerance of the gap
	tau:  the duality gap or the restarting technique is done every tau steps
	flag: if flag=1,  we apply the resart technique
	         flag=2,  just run the SFA algorithm, terminate it when the absolution change is less than tol
		     flag=3,  just run the SFA algorithm, terminate it when the duality gap is less than tol
			 flag=4,  just run the SFA algorithm, terminate it when the relative duality gap is less than tol


     We would like to emphasis that the following assumptions 
	       have been checked in the functions that call this function:
		   1) 0< lambda < z_max
		   2) nn >=2
		   3) z has been initialized with a starting point
		   4) z0 has been initialized with all zeros
		   
The termination condition is checked every tau iterations.

  For the duality gap, please refer to (12-18)
*/

/*
int sfa(double *x,     double *gap, int * activeS,
		 double *z,     double *z0,   double * v,   double * Av, 
		 double lambda, int nn,       int maxStep,
		 double *s,     double *g,
		 double tol,    int tau,       int flag){
*/
int sfa(double *x,     double *gap, int * activeS,
		 double *z,     double *z0,   double * v,   double * Av, 
		 double lambda, int nn,       int maxStep,
		 double *s,     double *g,
		 double tol,    int tau,       int flag, int d, int gamma1, int gamma2){

	int i, iterStep, m, tFlag=0, n=nn+d;
	double alphap=0, alpha=1, beta=0, temp;
	int* S=(int *) malloc(sizeof(int)*nn);
	double gapp=-1, gappp=-1;	/*gapp denotes the previous gap*/
	int numS=-1, numSp=-2, numSpp=-3;;
    int status;
                     
	/*
	numS denotes the number of elements in the Support Set S
	numSp denotes the number of elements in the previous Support Set S
	*/

	*gap=-1; /*initial a value -1*/

	/*
	The main algorithm by Nesterov's method

     B is an nn x nn tridiagonal matrix.

     The nn eigenvalues of B are 2- 2 cos (i * PI/ n), i=1, 2, ..., nn
	*/

	for (iterStep=1; iterStep<=maxStep; iterStep++){


		/*-------------   Step 1 ---------------------*/

		beta=(alphap -1 ) / alpha;
		/*
		compute search point
		    
			  s= z + beta * z0

        We follow the style of CLAPACK
		*/
		m=nn % 5;
		if (m!=0){
			for (i=0;i<m; i++)
				s[i]=z[i]+ beta* z0[i];
		}
		for (i=m;i<nn;i+=5){
			s[i]   =z[i]   + beta* z0[i];
			s[i+1] =z[i+1] + beta* z0[i+1];
			s[i+2] =z[i+2] + beta* z0[i+2];
			s[i+3] =z[i+3] + beta* z0[i+3];			
			s[i+4] =z[i+4] + beta* z0[i+4];
		}

		/*
		s and g are of size nn x 1

		compute the gradient at s

        g= B * s - Av,
		
		  where B is an nn x nn tridiagonal matrix. and is defined as

                                                B= [ 2  -1   0    0;
												     -1  2   -1   0;
													 0  -1   2    -1;
													 0   0   -1   2]

        We assume n>=3, which leads to nn>=2
		*/
		if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
            g[0]=z[0] + z[0] - z[1] - Av[0];                              
            for (i=1;i<nn-1;i++){
	            g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
            }
            g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];        
        }
        else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
            g[0] = 6.0 * z[0] - 4.0 * z[1] + z[2] - Av[0];                  
            g[1] = -4.0 * z[0] + 6.0 * z[1] - 4.0 * z[2] + z[3] - Av[1];		
            for (i=2;i<nn-2;i++){
	            g[i] = z[i-2] - 4.0 * z[i-1] + 6.0 * z[i] - 4.0 * z[i+1] + z[i+2] - Av[i];              
            }
            g[nn-2] = z[nn-4] - 4.0 * z[nn-3] + 6.0 * z[nn-2] - 4.0 * z[nn-1] - Av[nn-2]; 
            g[nn-1] = z[nn-3] - 4.0 * z[nn-2] + 6.0 * z[nn-1] - Av[nn-1];        
        }
        else {
            exit(-1);
        }

		/* 
		z0 stores the previous -z 
		*/
		m=nn%7;
		if (m!=0){
			for (i=0;i<m;i++)
				z0[i]=-z[i];
		}
		for (i=m; i <nn; i+=7){
			z0[i]   = - z[i];
			z0[i+1] = - z[i+1];
			z0[i+2] = - z[i+2];
			z0[i+3] = - z[i+3];
			z0[i+4] = - z[i+4];
			z0[i+5] = - z[i+5];
			z0[i+6] = - z[i+6];
		}
		

		/* 
		do a gradient step based on s to get z
		*/
		m=nn%5;
		if (m!=0){
			for(i=0;i<m; i++)
				z[i]=s[i] - g[i]/4;
		}
		for (i=m;i<nn; i+=5){			
			z[i]   = s[i]   -  g[i]  /4;
			z[i+1] = s[i+1] -  g[i+1]/4;
			z[i+2] = s[i+2] -  g[i+2]/4;
			z[i+3] = s[i+3] -  g[i+3]/4;
			z[i+4] = s[i+4] -  g[i+4]/4;
		}

		/*
		project z onto the L_{infty} ball with radius lambda

        z is the new approximate solution
		*/			
		for (i=0;i<nn; i++){
			if (z[i]>lambda)
				z[i]=lambda;
			else
				if (z[i]<-lambda)
					z[i]=-lambda;
		}
				
		/*
		compute the difference between the new solution 
		   and the previous solution (stored in z0=-z_p)

        the difference is written to z0
		*/
		
		m=nn%5;
		if (m!=0){
			for (i=0;i<m;i++)
				z0[i]+=z[i];
		}
		for(i=m;i<nn; i+=5){
			z0[i]  +=z[i];
			z0[i+1]+=z[i+1];
			z0[i+2]+=z[i+2];
			z0[i+3]+=z[i+3];
			z0[i+4]+=z[i+4];
		}

			
		alphap=alpha;
		alpha=(1+sqrt(4*alpha*alpha+1))/2;		

		/*
		check the termination condition
		*/
		if (iterStep%tau==0){


			/*
			The variables g and s can be modified

            The variables x, z0 and z can be revised for case 0, but not for the rest
			*/
			switch (flag){
			case 1:

			   /*

                terminate the program once the "duality gap" is smaller than tol

				compute the duality gap:
				
				       x= v - A^T z
				       Ax = Av - A A^T z = -g, 
				where
				       g = A A^T z - A v 

  
	            the duality gap= lambda * \|Ax\|-1 - <z, Ax>
		                       = lambda * \|g\|_1 + <z, g>

                In fact, gap=0 indicates that,
		             if g_i >0, then z_i=-lambda
					 if g_i <0, then z_i=lambda
				*/
				
				gappp=gapp;
				gapp=*gap;  /*record the previous gap*/
				numSpp=numSp;
				numSp=numS; /*record the previous numS*/

				dualityGap(gap, z, g, s, Av, lambda, nn, d);
				/*g is computed as the gradient of z in this function*/

				
				/*
				printf("\n Iteration: %d, gap=%e, numS=%d", iterStep, *gap, numS);
				*/
				
				/*
				If *gap <=tol, we terminate the iteration
				Otherwise, we restart the algorithm
				*/

				if (*gap <=tol){
					tFlag=1;
					break;

				} /* end of *gap <=tol */
				else{

					/* we apply the restarting technique*/

					/*
					we compute the solution by the second way
					*/
					numS = supportSet(x, v, z, g, S, lambda, nn, d);	
					/*g, the gradient of z should be computed before calling this function*/

					/*With x, we compute z via
					     AA^T z = Av - Ax
				    */

					/*
					printf("\n iterStep=%d, numS=%d, gap=%e",iterStep, numS, *gap);
					*/


					m=1;
					if (nn > 1000000)
						m=10;
					else
						if (nn > 100000)
							m=5;

					if ( abs(numS-numSp) < m){

						numS=generateSolution(x, z, gap, v, Av,
							                  g, s, S, lambda, nn, d, gamma1, gamma2);
						/*g, the gradient of z should be computed before calling this function*/

					
						if (*gap <tol){
							tFlag=2;	 /*tFlag =2 shows that the result is already optimal
										   There is no need to call generateSolution for recomputing the best solution
							              */					
							break;
						}

						if ( (*gap ==gappp) && (numS==numSpp) ){
					
							tFlag=2;
							break;

						}

						/*we terminate the program is *gap <1
						                          numS==numSP
												  and gapp==*gap
						*/
					}

					/*
                    compute s= Av -Ax
					*/
					if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
                        for (i=0;i<nn; i++) {
                    		s[i]=Av[i] - x[i+1] + x[i];      
                        }
                    }
                    else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
                        for (i=0;i < nn; i++) {
                            s[i]=Av[i] - (x[i] - 2.0 * x[i+1] + x[i+2]);
                        }
                    }
                    else {
                        exit(-1);
                    }

					if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
                        /*
		                Apply Rose Algorithm for solving z
		                */
					
		                Thomas(&temp, z, s, nn); 					
		                /*
		                Rose(&temp, z, s, nn);
		                */      
                    }
                    else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
                        /*
		                Apply PTRANS-1 Algorithm from Askar and Karawia (20150 for solving z
		                */
					
		                status = PTRANS1(&temp, z, s, nn);
                        if(status == -1) {
                            printf("ERROR: unable to solve system because mu[i]=0 with PTRANS-1 algorithm\n");
                        }      
                    }
                    else {
                        exit(-1);
                    }

					/*
					printf("\n Iteration: %d, %e", iterStep, temp);
					*/

					/*
					project z to [-lambda2, lambda2]
					*/
					for(i=0;i<nn;i++){
						if (z[i]>lambda)
							z[i]=lambda;
						else
							if (z[i]<-lambda)
								z[i]=-lambda;
					}

			
					
					m=nn%7;
					if (m!=0){
						for (i=0;i<m;i++)
							z0[i]=0;
					}
					for (i=m; i<nn; i+=7){
						z0[i]   = z0[i+1] 
							    = z0[i+2]
								= z0[i+3]
								= z0[i+4]
								= z0[i+5]
								= z0[i+6]
								=0;
					}

					
					alphap=0; alpha=1;

					/*
					we restart the algorithm
					*/

				}

				break; /*break case 1*/ 

			case 2: 

				/*
				The program is terminated either the summation of the absolution change (denoted by z0)
				    of z (from the previous zp) is less than tol * nn,
					     or the maximal number of iteration (maxStep) is achieved
				Note: tol indeed measures the averaged per element change.
				*/
				temp=0;
				m=nn%5;
				if (m!=0){
					for(i=0;i<m;i++)
						temp+=fabs(z0[i]);
				}
				for(i=m;i<nn;i+=5)
					temp=temp + fabs(z0[i])
					          + fabs(z0[i+1])
							  + fabs(z0[i+2])
							  + fabs(z0[i+3])
							  + fabs(z0[i+4]);
				*gap=temp / nn;

				if (*gap <=tol){

					tFlag=1;
				}

				break;

			case 3:

				/*

                terminate the program once the "duality gap" is smaller than tol

				compute the duality gap:
				
				       x= v - A^T z
				       Ax = Av - A A^T z = -g, 
				where
				       g = A A^T z - A v 

  
	            the duality gap= lambda * \|Ax\|-1 - <z, Ax>
		                       = lambda * \|g\|_1 + <z, g>

                In fact, gap=0 indicates that,
		             if g_i >0, then z_i=-lambda
					 if g_i <0, then z_i=lambda
				*/
				
								
				if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
		            g[0]=z[0] + z[0] - z[1] - Av[0];                              
		            for (i=1;i<nn-1;i++){
			            g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
		            }
		            g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];        
                }
                else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
                    g[0] = 6.0 * z[0] - 4.0 * z[1] + z[2] - Av[0];                  
                    g[1] = -4.0 * z[0] + 6.0 * z[1] - 4.0 * z[2] + z[3] - Av[1];		
                    for (i=2;i<nn-2;i++){
			            g[i] = z[i-2] - 4.0 * z[i-1] + 6.0 * z[i] - 4.0 * z[i+1] + z[i+2] - Av[i];              
		            }
		            g[nn-2] = z[nn-4] - 4.0 * z[nn-3] + 6.0 * z[nn-2] - 4.0 * z[nn-1] - Av[nn-2]; 
                    g[nn-1] = z[nn-3] - 4.0 * z[nn-2] + 6.0 * z[nn-1] - Av[nn-1];        
                }
                else {
                    exit(-1);
                }

				for (i=0;i<nn;i++)
					if (g[i]>0)
						s[i]=lambda + z[i];
					else
						s[i]=-lambda + z[i];

				temp=0;					
				m=nn%5;
				if (m!=0){
					for(i=0;i<m;i++)
						temp+=s[i]*g[i];
				}					
				for(i=m;i<nn;i+=5)
					temp=temp + s[i]  *g[i]
					          + s[i+1]*g[i+1]
				              + s[i+2]*g[i+2]
				              + s[i+3]*g[i+3]
				              + s[i+4]*g[i+4];
				*gap=temp;

				/*
				printf("\n %e", *gap);
				*/

					
				if (*gap <=tol)
					tFlag=1;

				break;

			case 4:

				/*

                terminate the program once the "relative duality gap" is smaller than tol
			       

				compute the duality gap:
				
				       x= v - A^T z
				       Ax = Av - A A^T z = -g, 
				where
				       g = A A^T z - A v 

  
	            the duality gap= lambda * \|Ax\|-1 - <z, Ax>
		                       = lambda * \|g\|_1 + <z, g>

                In fact, gap=0 indicates that,
		             if g_i >0, then z_i=-lambda
					 if g_i <0, then z_i=lambda

                
                Here, the "relative duality gap" is defined as:
				      duality gap / - 1/2 \|A^T z\|^2 + < z, Av>

                We efficiently compute - 1/2 \|A^T z\|^2 + < z, Av> using the following relationship

                      - 1/2 \|A^T z\|^2 + < z, Av>
					  = -1/2 <z, A A^T z - Av -Av>
					  = -1/2 <z, g - Av>
				*/
				
				
				if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
		            g[0]=z[0] + z[0] - z[1] - Av[0];                              
		            for (i=1;i<nn-1;i++){
			            g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
		            }
		            g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];        
                }
                else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
                    g[0] = 6.0 * z[0] - 4.0 * z[1] + z[2] - Av[0];                  
                    g[1] = -4.0 * z[0] + 6.0 * z[1] - 4.0 * z[2] + z[3] - Av[1];		
                    for (i=2;i<nn-2;i++){
			            g[i] = z[i-2] - 4.0 * z[i-1] + 6.0 * z[i] - 4.0 * z[i+1] + z[i+2] - Av[i];              
		            }
		            g[nn-2] = z[nn-4] - 4.0 * z[nn-3] + 6.0 * z[nn-2] - 4.0 * z[nn-1] - Av[nn-2]; 
                    g[nn-1] = z[nn-3] - 4.0 * z[nn-2] + 6.0 * z[nn-1] - Av[nn-1];        
                }
                else {
                    exit(-1);
                }

				for (i=0;i<nn;i++)
					if (g[i]>0)
						s[i]=lambda + z[i];
					else
						s[i]=-lambda + z[i];

				temp=0;					
				m=nn%5;
				if (m!=0){
					for(i=0;i<m;i++)
						temp+=s[i]*g[i];
				}					
				for(i=m;i<nn;i+=5)
					temp=temp + s[i]  *g[i]
					          + s[i+1]*g[i+1]
				              + s[i+2]*g[i+2]
				              + s[i+3]*g[i+3]
				              + s[i+4]*g[i+4];
				*gap=temp;
				/*
				Now, *gap contains the duality gap
				Next, we compute
				   - 1/2 \|A^T z\|^2 + < z, Av>
				   =-1/2 <z, g - Av>
				*/

				temp=0;
				m=nn%5;
				if (m!=0){
					for(i=0;i<m;i++)
						temp+=z[i] * (g[i] - Av[i]);
				}					
				for(i=m;i<nn;i+=5)
					temp=temp + z[i]  * (g[i] -  Av[i])
					          + z[i+1]* (g[i+1]- Av[i+1])
				              + z[i+2]* (g[i+2]- Av[i+2])
				              + z[i+3]* (g[i+3]- Av[i+3])
				              + z[i+4]* (g[i+4]- Av[i+4]);
				temp=fabs(temp) /2; 

				if (temp <1)
					temp=1;

				*gap/=temp;
				/*
				*gap now contains the relative gap
				*/

					
				if (*gap <=tol){
					tFlag=1;
				}

				break;

			default:

				/*
				The program is terminated either the summation of the absolution change (denoted by z0)
				    of z (from the previous zp) is less than tol * nn,
					     or the maximal number of iteration (maxStep) is achieved
				Note: tol indeed measures the averaged per element change.
				*/
				temp=0;
				m=nn%5;
				if (m!=0){
					for(i=0;i<m;i++)
						temp+=fabs(z0[i]);
				}
				for(i=m;i<nn;i+=5)
					temp=temp + fabs(z0[i])
					          + fabs(z0[i+1])
							  + fabs(z0[i+2])
							  + fabs(z0[i+3])
							  + fabs(z0[i+4]);
				*gap=temp / nn;

				if (*gap <=tol){

					tFlag=1;
				}

				break;

			}/*end of switch*/
			
			if (tFlag)
				break;
			
		}/* end of the if for checking the termination condition */

		/*-------------- Step 3 --------------------*/
		
	}

	/*
	for the other cases, except flag=1, compute the solution x according the first way (the primal-dual way)
	*/
	
    /* NOTE TO SELF. THIS IS x = v - A^T * z (z is p-d dimensional vector), x and v are p dimensional vectors. */
    if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
	    if ( (flag !=1) || (tFlag==0) ){
		    x[0]=v[0] + z[0];
		    for(i=1;i<n-1;i++)
			    x[i]= v[i] - z[i-1] + z[i];
		    x[n-1]=v[n-1] - z[n-2];
	    }        
    }
    else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
        if ( (flag !=1) || (tFlag==0) ){
		    x[0] = v[0] - z[0];
            x[1] = v[1] + 2.0 * z[0] - z[1];		    
            for(i = 2; i < n-2; i++) {
			    x[i] = v[i] - z[i-2] + 2.0 * z[i-1] - z[i];
            }
		    x[n-2] = v[n-2] - z[n-4] + 2.0 * z[n-3];
            x[n-1] = v[n-1] - z[n-2];
	    }        
    }
    else {
        exit(-1);
    }

    

	if ( (flag==1) && (tFlag==1)){
	
	    /*
		We assume that n>=3, and thus nn>=2
		
		  We have two ways for recovering x. 
		  The first way is x = v - A^T z
		  The second way is x =omega(z)
		*/
		
		/*
         We first compute the objective function value of the first choice in terms f(x), see our paper
		*/
				
		/*
		for numerical reason, we do a gradient descent step
		*/

		/*
		---------------------------------------------------
		  A gradient step  begins
		*/
		if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
		    g[0]=z[0] + z[0] - z[1] - Av[0];                              
		    for (i=1;i<nn-1;i++){
			    g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
		    }
		    g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];        
        }
        else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
            g[0] = 6.0 * z[0] - 4.0 * z[1] + z[2] - Av[0];                  
            g[1] = -4.0 * z[0] + 6.0 * z[1] - 4.0 * z[2] + z[3] - Av[1];		
            for (i=2;i<nn-2;i++){
			    g[i] = z[i-2] - 4.0 * z[i-1] + 6.0 * z[i] - 4.0 * z[i+1] + z[i+2] - Av[i];              
		    }
		    g[nn-2] = z[nn-4] - 4.0 * z[nn-3] + 6.0 * z[nn-2] - 4.0 * z[nn-1] - Av[nn-2]; 
            g[nn-1] = z[nn-3] - 4.0 * z[nn-2] + 6.0 * z[nn-1] - Av[nn-1];        
        }
        else {
            exit(-1);
        }

		
		/* 
		do a gradient step based on z to get the new z
		*/
		m=nn%5;
		if (m!=0){
			for(i=0;i<m; i++)
				z[i]=z[i] - g[i]/4;
		}
		for (i=m;i<nn; i+=5){			
			z[i]   = z[i]   -  g[i]  /4;
			z[i+1] = z[i+1] -  g[i+1]/4;
			z[i+2] = z[i+2] -  g[i+2]/4;
			z[i+3] = z[i+3] -  g[i+3]/4;
			z[i+4] = z[i+4] -  g[i+4]/4;
		}

		/*
		project z onto the L_{infty} ball with radius lambda

        z is the new approximate solution
		*/			
		for (i=0;i<nn; i++){
			if (z[i]>lambda)
				z[i]=lambda;
			else
				if (z[i]<-lambda)
					z[i]=-lambda;
		}

		/*
		---------------------------------------------------
		  A gradient descent step ends
		*/

		/*compute the gradient at z*/
		
		if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
		    g[0]=z[0] + z[0] - z[1] - Av[0];                              
		    for (i=1;i<nn-1;i++){
			    g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
		    }
		    g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];        
        }
        else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
            g[0] = 6.0 * z[0] - 4.0 * z[1] + z[2] - Av[0];                  
            g[1] = -4.0 * z[0] + 6.0 * z[1] - 4.0 * z[2] + z[3] - Av[1];		
            for (i=2;i<nn-2;i++){
			    g[i] = z[i-2] - 4.0 * z[i-1] + 6.0 * z[i] - 4.0 * z[i+1] + z[i+2] - Av[i];              
		    }
		    g[nn-2] = z[nn-4] - 4.0 * z[nn-3] + 6.0 * z[nn-2] - 4.0 * z[nn-1] - Av[nn-2]; 
            g[nn-1] = z[nn-3] - 4.0 * z[nn-2] + 6.0 * z[nn-1] - Av[nn-1];        
        }
        else {
            exit(-1);
        }

		numS=generateSolution(x, z, gap, v, Av,
							                  g, s, S, lambda, nn, d, gamma1, gamma2);
		/*g, the gradient of z should be computed before calling this function*/

	}

	free (S);
	/*
	free the variables S
	*/

	*activeS=numS;
	return (iterStep);

}


/*

Refer to sfa for the defintions of the variables  

In this file, we restart the program every step, and neglect the gradient step.

  It seems that, this program does not converge.

  This function shows that the gradient step is necessary.
*/

/*
int sfa_special(double *x,     double *gap,  int * activeS,
		 double *z,     double * v,   double * Av, 
		 double lambda, int nn,       int maxStep,
		 double *s,     double *g,
		 double tol,    int tau){
*/
int sfa_special(double *x,     double *gap,  int * activeS,
		 double *z,     double * v,   double * Av, 
		 double lambda, int nn,       int maxStep,
		 double *s,     double *g,
		 double tol,    int tau, int d){

	int i, iterStep, tFlag=0, n=nn+d;
	double temp;
	int* S=(int *) malloc(sizeof(int)*nn);
	double gapp=-1;	/*gapp denotes the previous gap*/
	int numS=-1, numSp=-1;  
    int status;  
                                                                       
	/*
	numS denotes the number of elements in the Support Set S
	numSp denotes the number of elements in the previous Support Set S
	*/

	*gap=-1; /*initialize *gap a value*/

	for (iterStep=1; iterStep<=maxStep; iterStep++){

		if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
		    g[0]=z[0] + z[0] - z[1] - Av[0];                              
		    for (i=1;i<nn-1;i++){
			    g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
		    }
		    g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];        
        }
        else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
            g[0] = 6.0 * z[0] - 4.0 * z[1] + z[2] - Av[0];                  
            g[1] = -4.0 * z[0] + 6.0 * z[1] - 4.0 * z[2] + z[3] - Av[1];		
            for (i=2;i<nn-2;i++){
			    g[i] = z[i-2] - 4.0 * z[i-1] + 6.0 * z[i] - 4.0 * z[i+1] + z[i+2] - Av[i];              
		    }
		    g[nn-2] = z[nn-4] - 4.0 * z[nn-3] + 6.0 * z[nn-2] - 4.0 * z[nn-1] - Av[nn-2]; 
            g[nn-1] = z[nn-3] - 4.0 * z[nn-2] + 6.0 * z[nn-1] - Av[nn-1];        
        }
        else {
            exit(-1);
        }

		numSp=numS; /*record the previous numS*/
		numS = supportSet(x, v, z, g, S, lambda, nn, d);
		
		
		/*With x, we compute z via
		AA^T z = Av - Ax
		 */
	
		/*
		compute s= Av -Ax
		*/
		
		if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
            for (i=0;i<nn; i++) {
        		s[i]=Av[i] - x[i+1] + x[i];      
            }
        }
        else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
            for (i=0; i < nn; i++) {
                s[i]=Av[i] - (x[i] - 2.0 * x[i+1] + x[i+2]);
            }
        }
        else {
            exit(-1);
        }			
	

	    if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
            /*
		    Apply Rose Algorithm for solving z
		    */
					
		    Thomas(&temp, z, s, nn); 					
		    /*
		    Rose(&temp, z, s, nn);
		    */      
        }
        else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
            /*
		    Apply PTRANS-1 Algorithm from Askar and Karawia (20150 for solving z
		    */
					
		    status = PTRANS1(&temp, z, s, nn); 		
            if(status == -1) {
                printf("ERROR: unable to solve system because mu[i]=0 with PTRANS-1 algorithm\n");
            }      
        }
        else {
            exit(-1);
        }

	   /*
	   project z to [-lambda, lambda]
	   */
	
		for(i=0;i<nn;i++){		
			if (z[i]>lambda)
				z[i]=lambda;
			else
				if (z[i]<-lambda)
					z[i]=-lambda;
		}


		if (iterStep%tau==0){
			gapp=*gap;  /*record the previous gap*/

			dualityGap(gap, z, g, s, Av, lambda, nn, d);

			/*
			printf("\n iterStep=%d, numS=%d, gap=%e, diff=%e",iterStep, numS, *gap, *gap -gapp);

            */
			
			if (*gap <=tol){
				tFlag=1;
				break;
			}

			if ( (*gap <1) && (numS==numSp) && fabs(gapp == *gap) ){
				tFlag=1;			
				break;
			/*we terminate the program is *gap <1
			                         numS==numSP
								   and gapp==*gap
		    */
			}

		}/*end of if tau*/
				
	}/*end for */		
	
	free (S);

	* activeS=numS;
	return(iterStep);

}


/*

  We do one gradient descent, and then restart the program
*/


/*
int sfa_one(double *x,     double *gap, int * activeS,
		 double *z,     double * v,   double * Av, 
		 double lambda, int nn,       int maxStep,
		 double *s,     double *g,
		 double tol,    int tau){
*/
int sfa_one(double *x,     double *gap, int * activeS,
		 double *z,     double * v,   double * Av, 
		 double lambda, int nn,       int maxStep,
		 double *s,     double *g,
		 double tol,    int tau, int d, int gamma1, int gamma2){

	int i, iterStep, m, tFlag=0, n=nn+d;
	double temp;
	int* S=(int *) malloc(sizeof(int)*nn);
	double gapp=-1, gappp=-2;	/*gapp denotes the previous gap*/
	int numS=-100, numSp=-200, numSpp=-300;    
    int status;
              
	/*
	numS denotes the number of elements in the Support Set S
	numSp denotes the number of elements in the previous Support Set S
	*/

	*gap=-1; /*initialize *gap a value*/

	/*
	The main algorithm by Nesterov's method

     B is an nn x nn tridiagonal matrix.

     The nn eigenvalues of B are 2- 2 cos (i * PI/ n), i=1, 2, ..., nn
	*/


	/*
	we first do a gradient step based on z
	*/


	/*
		---------------------------------------------------
		  A gradient step  begins
	*/

    if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
		g[0]=z[0] + z[0] - z[1] - Av[0];                              
		for (i=1;i<nn-1;i++){
			g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
		}
		g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];        
    }
    else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
        g[0] = 6.0 * z[0] - 4.0 * z[1] + z[2] - Av[0];                  
        g[1] = -4.0 * z[0] + 6.0 * z[1] - 4.0 * z[2] + z[3] - Av[1];		
        for (i=2;i<nn-2;i++){
			g[i] = z[i-2] - 4.0 * z[i-1] + 6.0 * z[i] - 4.0 * z[i+1] + z[i+2] - Av[i];              
		}
		g[nn-2] = z[nn-4] - 4.0 * z[nn-3] + 6.0 * z[nn-2] - 4.0 * z[nn-1] - Av[nn-2]; 
        g[nn-1] = z[nn-3] - 4.0 * z[nn-2] + 6.0 * z[nn-1] - Av[nn-1];        
    }
    else {
        exit(-1);
    }

		
		/* 
		do a gradient step based on z to get the new z
		*/
		m=nn%5;
		if (m!=0){
			for(i=0;i<m; i++)
				z[i]=z[i] - g[i]/4;
		}
		for (i=m;i<nn; i+=5){			
			z[i]   = z[i]   -  g[i]  /4;
			z[i+1] = z[i+1] -  g[i+1]/4;
			z[i+2] = z[i+2] -  g[i+2]/4;
			z[i+3] = z[i+3] -  g[i+3]/4;
			z[i+4] = z[i+4] -  g[i+4]/4;
		}

		/*
		project z onto the L_{infty} ball with radius lambda

        z is the new approximate solution
		*/			
		for (i=0;i<nn; i++){
			if (z[i]>lambda)
				z[i]=lambda;
			else
				if (z[i]<-lambda)
					z[i]=-lambda;
		}

		/*
		---------------------------------------------------
		  A gradient descent step ends
		*/


	/*compute the gradient at z*/

    if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
		g[0]=z[0] + z[0] - z[1] - Av[0];                              
		for (i=1;i<nn-1;i++){
			g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
		}
		g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];        
    }
    else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
        g[0] = 6.0 * z[0] - 4.0 * z[1] + z[2] - Av[0];                  
        g[1] = -4.0 * z[0] + 6.0 * z[1] - 4.0 * z[2] + z[3] - Av[1];		
        for (i=2;i<nn-2;i++){
			g[i] = z[i-2] - 4.0 * z[i-1] + 6.0 * z[i] - 4.0 * z[i+1] + z[i+2] - Av[i];              
		}
		g[nn-2] = z[nn-4] - 4.0 * z[nn-3] + 6.0 * z[nn-2] - 4.0 * z[nn-1] - Av[nn-2]; 
        g[nn-1] = z[nn-3] - 4.0 * z[nn-2] + 6.0 * z[nn-1] - Av[nn-1];        
    }
    else {
        exit(-1);
    }
    

	for (iterStep=1; iterStep<=maxStep; iterStep++){


		/*
		---------------------------------------------------
		restart the algorithm with x=omega(z)
		*/
				
		numSpp=numSp;
		numSp=numS; /*record the previous numS*/
		numS = supportSet(x, v, z, g, S, lambda, nn, d);
		
		
		/*With x, we compute z via
		AA^T z = Av - Ax
		 */
	
		/*
		compute s= Av -Ax
		*/
		
		if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
            for (i=0;i<nn; i++) {
        		s[i]=Av[i] - x[i+1] + x[i];      
            }
        }
        else if(d == 2) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
            for (i=0;i<nn; i++) {
                s[i]=Av[i] - (x[i] - 2.0 * x[i+1] + x[i+2]);
            }
        }
        else {
            exit(-1);
        }			
	

		if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
            /*
		    Apply Rose Algorithm for solving z
		    */
					
		    Thomas(&temp, z, s, nn); 					
		    /*
		    Rose(&temp, z, s, nn);
		    */      
        }
        else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
            /*
		    Apply PTRANS-1 Algorithm from Askar and Karawia (20150 for solving z
		    */
					
		    status = PTRANS1(&temp, z, s, nn); 			
            if(status == -1) {
                printf("ERROR: unable to solve system because mu[i]=0 with PTRANS-1 algorithm\n");
            }      
        }
        else {
            exit(-1);
        }

	    /*
	    project z to [-lambda, lambda]
	    */
	
		for(i=0;i<nn;i++){		
			if (z[i]>lambda)
				z[i]=lambda;
			else
				if (z[i]<-lambda)
					z[i]=-lambda;
		}

		/*
		---------------------------------------------------
		restart the algorithm with x=omega(z)

        we have computed a new z, based on the above relationship
		*/


		/*
		---------------------------------------------------
		  A gradient step  begins
		*/       
        if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
		    g[0]=z[0] + z[0] - z[1] - Av[0];                              
		    for (i=1;i<nn-1;i++){
			    g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
		    }
		    g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];        
        }
        else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
            g[0] = 6.0 * z[0] - 4.0 * z[1] + z[2] - Av[0];                  
            g[1] = -4.0 * z[0] + 6.0 * z[1] - 4.0 * z[2] + z[3] - Av[1];		
            for (i=2;i<nn-2;i++){
			    g[i] = z[i-2] - 4.0 * z[i-1] + 6.0 * z[i] - 4.0 * z[i+1] + z[i+2] - Av[i];              
		    }
		    g[nn-2] = z[nn-4] - 4.0 * z[nn-3] + 6.0 * z[nn-2] - 4.0 * z[nn-1] - Av[nn-2]; 
            g[nn-1] = z[nn-3] - 4.0 * z[nn-2] + 6.0 * z[nn-1] - Av[nn-1];        
        }
        else {
            exit(-1);
        }
		
		/* 
		do a gradient step based on z to get the new z
		*/
		m=nn%5;
		if (m!=0){
			for(i=0;i<m; i++)
				z[i]=z[i] - g[i]/4;
		}
		for (i=m;i<nn; i+=5){			
			z[i]   = z[i]   -  g[i]  /4;
			z[i+1] = z[i+1] -  g[i+1]/4;
			z[i+2] = z[i+2] -  g[i+2]/4;
			z[i+3] = z[i+3] -  g[i+3]/4;
			z[i+4] = z[i+4] -  g[i+4]/4;
		}

		/*
		project z onto the L_{infty} ball with radius lambda

        z is the new approximate solution
		*/			
		for (i=0;i<nn; i++){
			if (z[i]>lambda)
				z[i]=lambda;
			else
				if (z[i]<-lambda)
					z[i]=-lambda;
		}

		/*
		---------------------------------------------------
		  A gradient descent step ends
		*/

		/*compute the gradient at z*/
        if(d == 1) {    /* #### HERE IS THE CONSTANT TREND PENALTY #### */
		    g[0]=z[0] + z[0] - z[1] - Av[0];                              
		    for (i=1;i<nn-1;i++){
			    g[i]= - z[i-1] + z[i] + z[i] - z[i+1] - Av[i];
		    }
		    g[nn-1]= -z[nn-2] + z[nn-1] + z[nn-1] - Av[nn-1];        
        }
        else if(d == 2) {    /* #### HERE WE MODIFIED TO INCLUDE THE LINEAR TREND PENALTY #### */
            g[0] = 6.0 * z[0] - 4.0 * z[1] + z[2] - Av[0];                  
            g[1] = -4.0 * z[0] + 6.0 * z[1] - 4.0 * z[2] + z[3] - Av[1];		
            for (i=2;i<nn-2;i++){
			    g[i] = z[i-2] - 4.0 * z[i-1] + 6.0 * z[i] - 4.0 * z[i+1] + z[i+2] - Av[i];              
		    }
		    g[nn-2] = z[nn-4] - 4.0 * z[nn-3] + 6.0 * z[nn-2] - 4.0 * z[nn-1] - Av[nn-2]; 
            g[nn-1] = z[nn-3] - 4.0 * z[nn-2] + 6.0 * z[nn-1] - Av[nn-1];        
        }
        else {
            exit(-1);
        }

		if (iterStep % tau==0){
			gappp=gapp;
			gapp=*gap;  /*record the previous gap*/

			dualityGap2(gap, z, g, s, Av, lambda, nn);
			/*g, the gradient of z should be computed before calling this function*/


			/*
			printf("\n iterStep=%d, numS=%d, gap=%e",iterStep, numS, *gap);
			*/

		
			/*
			printf("\n  %d  & %d   &  %2.0e \\\\ \n \\hline ",iterStep, numS, *gap);
			*/
			

			/*
			printf("\n %e",*gap);
			*/

			/*		

			printf("\n %d",numS);

			*/
			
			if (*gap <=tol){
				tFlag=1;
				break;
			}

			m=1;
			if (nn > 1000000)
				m=5;
			else
				if (nn > 100000)
					m=3;

			if ( abs( numS-numSp) <m ){

				/*
				printf("\n numS=%d, numSp=%d",numS,numSp);
				*/

				m=generateSolution(x, z, gap, v, Av,
					                  g, s, S, lambda, nn, d, gamma1, gamma2);
				/*g, the gradient of z should be computed before calling this function*/

				if (*gap < tol){

					numS=m;
					tFlag=2;					
					break;
				}


				if ( (*gap ==gappp) && (numS==numSpp) ){
					
					tFlag=2;
					break;
							
				}
				
            } /*end of if*/

		}/*end of if tau*/


	} /*end of for*/



	if (tFlag!=2){
		numS=generateSolution(x, z, gap, v, Av, g, s, S, lambda, nn, d, gamma1, gamma2);
       /*g, the gradient of z should be computed before calling this function*/
    }

	free(S);

	*activeS=numS;
	return(iterStep);
}

