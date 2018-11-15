//THIS CODE HAS BEEN STRONGLY INSPIRED BY SZABOS BOOK ON HF THEORY
//COMPILE WITH gcc quantum.c -o quantum -llapack -lblas -lm
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#define DEBUG    0
#define PI       3.1415926535898
#define SIZE     2
#define max(a,b)     ((a)>(b)?(a):(b))
#define SQR(a)       ((a)*(a))

typedef struct elem {
  int size;
  double H[SIZE][SIZE];                 //HAMILTONIAN
  double S[SIZE][SIZE];                 //OVERLAP MATRIX
  double TT[SIZE][SIZE][SIZE][SIZE];    //TWOE MATRIX
  double X[SIZE][SIZE];
  double XT[SIZE][SIZE];
  double G[SIZE][SIZE];
  double P[SIZE][SIZE];
  int choice;
} elem;

typedef struct CONTAINER {
  elem dat;
} ELEMENT;

typedef ELEMENT *SYSM;

void tqli(double*, double*, int, double*);
void diag2(double*, int, double*, double*);
void HFcalc(int, int, float, float*, int*, SYSM);
void integrals(int, int, float, float*, int*, SYSM);
void scf(int, int, float, float*, int*, SYSM);
void scf_2(int, int, float, float*, int*, SYSM);
void matout(double[SIZE][SIZE],char*,int,int);
void vectout(double*,char*,int);
void ttout(double[SIZE][SIZE][SIZE][SIZE],char*, int);
void formG(SYSM);
void diag(double[SIZE][SIZE],double[SIZE][SIZE],double[SIZE][SIZE]);
void mult(double[SIZE][SIZE],double[SIZE][SIZE],double[SIZE][SIZE],int);
void test(void);
void copy(double[SIZE][SIZE],double[SIZE][SIZE], int, int);

float pythag(float, float);
double S(double,double,double);
double T(double,double,double);
double V(double,double,double,double,float);
double E2I(double,double,double,double,double,double,double);
double F0(double);

int main()
{
  int IOP=0,N=2,Z[2]={2,1};
  int val=0;
  float R=1.4632,ZETA[2]={2.0925,1.24};
  SYSM HS;  
  FILE *fptr;

  fptr=fopen("inp","r");
  if(!(fptr==NULL)){
    fscanf(fptr,"%d %d %d %d %d\n",&IOP,&N,&Z[0],&Z[1],&val);
    fscanf(fptr,"%f %f %f",&R,&ZETA[0],&ZETA[1]);
  }
  fclose(fptr);

  HS = malloc(sizeof(ELEMENT));
  HS->dat.choice=val;
  HFcalc(IOP,N,R,ZETA,Z,HS);

  return 0;
}

void HFcalc(int IOP, int N, float R, float ZETA[],int Z[], SYSM HS)
{
  //Do the one and two electron integral calculations
  integrals(IOP,N,R,ZETA,Z,HS);
  scf(IOP,N,R,ZETA,Z,HS);
}

void integrals(int IOP, int N, float R, float ZETA[],int Z[],SYSM HS)
{
  int M=3;
  double RAP[2],RBP[2],RAQ[2],RBQ[2],RPQ[2],R2;
  double RES[SIZE*SIZE];
  double coef[3*3] ={1.0,0.0,0.0,0.678914,0.430129,0.0,0.444635,0.535328,0.154329};
  double expon[3*3]={0.270950,0.0,0.0,0.151623,0.851819,0.0,0.109818,0.4057710,2.22766};
  double D1[3],A1[3],D2[3],A2[3];
  int i,j,k,l,s,m=M-1,n=N-1,s2,s3;
  double** Aptr;
  R2=R*R;

  //D are unnormalized contraction coefficients 
  //A are unnormalized contraction exponents
  if(IOP<2){
    printf("\n\nCOEFFICIENTS AND EXPONENTS\n");
    printf("---------------------------\n");   
  }
  for(i=0;i<N;i++){
    k=(N-1)*M+i;
    A1[i]=expon[k]*pow(ZETA[0],2);
    D1[i]=coef[k]*pow(2.0*A1[i]/PI,0.75);
    A2[i]=expon[k]*pow(ZETA[1],2);
    D2[i]=coef[k]*pow(2.0*A2[i]/PI,0.75);
    if(IOP<2){
      printf("   %f13 %f13\n",coef[k],expon[k]);
    }
  }
  if(IOP<2)
    printf("---------------------------\n");

  HS->dat.size=SIZE;
  s=HS->dat.size;
  s2=s*s;
  s3=s2*s;
  //INIT
  for(i=0;i<s;i++){
    for(j=0;j<s;j++){
      HS->dat.H[i][j]=0.0;
      for(k=0;k<s;k++)
	for(l=0;l<s;l++)
	  HS->dat.TT[i][j][k][l]=0.0;
    }
  }

  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      RAP[0]=A2[j]*R/(A1[i]+A2[j]);
      RAP[1]=RAP[0]*RAP[0];
      RBP[1]=(R-RAP[0])*(R-RAP[0]);
      HS->dat.S[0][1]+=S(A1[i],A2[j],R2)*D1[i]*D2[j];
      HS->dat.H[0][0]+=T(A1[i],A1[j],0.0)*D1[i]*D1[j];
      HS->dat.H[0][1]+=T(A1[i],A2[j],R2)*D1[i]*D2[j];
      HS->dat.H[1][1]+=T(A2[i],A2[j],0.0)*D2[i]*D2[j];
      HS->dat.H[0][0]+=V(A1[i],A1[j],0.0,0.0,Z[0])*D1[i]*D1[j];
      HS->dat.H[0][1]+=V(A1[i],A2[j],R2,RAP[1],Z[0])*D1[i]*D2[j];
      HS->dat.H[1][1]+=V(A2[i],A2[j],0.0,R2,Z[0])*D2[i]*D2[j];
      HS->dat.H[0][0]+=V(A1[i],A1[j],0.0,R2,Z[1])*D1[i]*D1[j];
      HS->dat.H[0][1]+=V(A1[i],A2[j],R2,RBP[1],Z[1])*D1[i]*D2[j];
      HS->dat.H[1][1]+=V(A2[i],A2[j],0.0,0.0,Z[1])*D2[i]*D2[j];
    }
  }

  //2^3 2^2 2^1 
  //2e- integrals
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
	for(l=0;l<N;l++){
	  RAP[0]=A2[i]*R/(A2[i]+A1[j]);
	  RAP[1]=RAP[0]*RAP[0];
	  RBP[0]=R-RAP[0];RBP[1]=RBP[0]*RBP[0];
	  RAQ[0]=A2[k]*R/(A2[k]+A1[l]);RAQ[1]=RAQ[0]*RAQ[0];
	  RBQ[0]=R-RAQ[0]; RBQ[1]=RBQ[0]*RBQ[0];
	  RPQ[0]=RAP[0]-RAQ[0];RPQ[1]=RPQ[0]*RPQ[0];
	  HS->dat.TT[0][0][0][0]+=E2I(A1[i],A1[j],A1[k],A1[l],0.0,0.0,0.0)*D1[i]*D1[j]*D1[k]*D1[l];
	  HS->dat.TT[1][0][0][0]+=E2I(A2[i],A1[j],A1[k],A1[l],R2,0.0,RAP[1])*D2[i]*D1[j]*D1[k]*D1[l];
	  HS->dat.TT[1][0][1][0]+=E2I(A2[i],A1[j],A2[k],A1[l],R2,R2,RPQ[1])*D2[i]*D1[j]*D2[k]*D1[l];
	  HS->dat.TT[1][1][0][0]+=E2I(A2[i],A2[j],A1[k],A1[l],0.0,0.0,R2)*D2[i]*D2[j]*D1[k]*D1[l];
	  HS->dat.TT[1][1][1][0]+=E2I(A2[i],A2[j],A2[k],A1[l],0.0,R2,RBQ[1])*D2[i]*D2[j]*D2[k]*D1[l];
	  HS->dat.TT[1][1][1][1]+=E2I(A2[i],A2[j],A2[k],A2[l],0.0,0.0,0.0)*D2[i]*D2[j]*D2[k]*D2[l];
	}
      }
    }
  }
  HS->dat.TT[0][1][0][0]=HS->dat.TT[1][0][0][0];
  HS->dat.TT[0][0][1][0]=HS->dat.TT[1][0][0][0];
  HS->dat.TT[0][0][0][1]=HS->dat.TT[1][0][0][0];
  HS->dat.TT[0][1][1][0]=HS->dat.TT[1][0][1][0];
  HS->dat.TT[1][0][0][1]=HS->dat.TT[1][0][1][0];
  HS->dat.TT[0][1][0][1]=HS->dat.TT[1][0][1][0];
  HS->dat.TT[0][0][1][1]=HS->dat.TT[1][1][0][0];
  HS->dat.TT[1][1][0][1]=HS->dat.TT[1][1][1][0];
  HS->dat.TT[1][0][1][1]=HS->dat.TT[1][1][1][0];
  HS->dat.TT[0][1][1][1]=HS->dat.TT[1][1][1][0];

  HS->dat.H[1][0]=HS->dat.H[0][1];

  HS->dat.S[0][0]=1.0;
  HS->dat.S[1][0]=HS->dat.S[0][1];
  HS->dat.S[1][1]=1.0;

  HS->dat.X[0][0]=1.0/sqrt(2.0*(1.0+HS->dat.S[0][1]));
  HS->dat.X[1][0]=HS->dat.X[0][0];
  HS->dat.X[0][1]=1.0/sqrt(2.0*(1.0-HS->dat.S[0][1]));
  HS->dat.X[1][1]=-HS->dat.X[0][1];

  HS->dat.XT[0][0]=HS->dat.X[0][0];
  HS->dat.XT[0][1]=HS->dat.X[1][0];
  HS->dat.XT[1][0]=HS->dat.X[0][1];
  HS->dat.XT[1][1]=HS->dat.X[1][1];

  if(IOP<2){
    matout(HS->dat.S,"S",SIZE,SIZE);
    matout(HS->dat.H,"H",SIZE,SIZE);
    matout(HS->dat.X,"X",SIZE,SIZE);
    ttout(HS->dat.TT,"TWO ELECTRON",SIZE);
  }
}

void scf_2(int IOP, int N, float R, float ZETA[],int Z[],SYSM HS)
{
  double C[SIZE][SIZE],CPR[SIZE][SIZE],OLDP[SIZE][SIZE],E[SIZE][SIZE];
  double F[SIZE][SIZE],FPR[SIZE][SIZE],H[SIZE][SIZE],EN,ENT;
  double FP[SIZE][SIZE],ONE[4]={1,0,0,1};
  double WORK[26*SIZE];
  double tol[3]={1e-6,1e-6,1e-6},delta=1.0;
  int ISUPPZ[SIZE][SIZE],IWORK[10*SIZE];
  int MAXIT[3]={100,1000,1},choice,iter=0;  
  int i,j,k,l,s,M;
  choice=HS->dat.choice;
  s=SIZE;

  for(i=0;i<s;i++)
    for(j=0;j<s;j++)
      HS->dat.P[i][j]=0.0;
  while(iter<MAXIT[choice] && delta > tol[choice]){
    iter++;
    if(IOP<2)
      fprintf(stderr,"::::::::::::::::::::::::>  ITERATION:: %d\n",iter);
    formG(HS);
    if(IOP<2){
      matout(HS->dat.P,"P",s,s);
      matout(HS->dat.G,"G",s,s);
    }
    //create fock matrix
    for(i=0;i<s;i++)
      for(j=0;j<s;j++)
	F[i][j]=HS->dat.G[i][j]+HS->dat.H[i][j];
    EN=0.0;
    for(i=0;i<s;i++)
      for(j=0;j<s;j++)
	EN+=0.5*(HS->dat.P[i][j])*(F[i][j]+HS->dat.H[i][j]);
    if(IOP<2)
      matout(F,"F",SIZE,SIZE);
    fprintf(stderr,"\n:::Energy::: %e\n",EN);
    mult(F,HS->dat.X,HS->dat.G,s);
    mult(HS->dat.XT,HS->dat.G,FPR,s);
    // DIAGONALIZE FPR 
    // EIGENVECT IN CPR 
    // EIGENV. IN E
    diag(FPR,CPR,E);
    mult(HS->dat.X,CPR,C,s);
    for(i=0;i<s;i++)
      for(j=0;j<s;j++){
	OLDP[i][j]=HS->dat.P[i][j];
	HS->dat.P[i][j]=0.0;
	for(k=0;k<1;k++)
	  HS->dat.P[i][j]+=2.0*C[i][k]*C[j][k];
      }
    if(IOP<2){
      matout(FPR,"F'",s,s);
      matout(CPR,"C'",s,s);
      matout(E,"E",s,s);
      matout(C,"C",s,s);
      matout(HS->dat.P,"P",s,s);
    }
    delta=0.0;
    for(i=0;i<s;i++)
      for(j=0;j<s;j++)
	delta+=pow(HS->dat.P[i][j]-OLDP[i][j],2.0);
    delta=sqrt(delta/4);
    fprintf(stderr,"CONVERGENCE OF DENSITY: %e\n",delta);
  }
  //CONVERGENCE
  ENT=EN+Z[0]*Z[1]/R;
  fprintf(stderr,"::::::::::::::::::::::::>  FINISHED ITERATIONS\n");
  if(delta > tol[choice])
    fprintf(stderr,"\nWARNING::METHOD DID NOT CONVERGE\n");
  fprintf(stderr,"\nTotal energy:: %e\nElectronic energy:: %e\n\n",ENT,EN);
  if(IOP>1){
     matout(HS->dat.G,"G",SIZE,SIZE);
     matout(F,"F",SIZE,SIZE);
     matout(E,"E",SIZE,SIZE);
     matout(C,"C",SIZE,SIZE);
     matout(HS->dat.P,"P",SIZE,SIZE);
  }
  //MULLIKEN CHARGES
  mult(HS->dat.P,HS->dat.S,OLDP,SIZE);
  if(IOP>0)
    matout(OLDP,"MULLIKEN",SIZE,SIZE);
}

void scf(int IOP, int N, float R, float ZETA[],int Z[],SYSM HS)
{
  double C[SIZE][SIZE],CPR[SIZE][SIZE],OLDP[SIZE][SIZE],E[SIZE][SIZE];
  double F[SIZE][SIZE],FPR[SIZE][SIZE],H[SIZE][SIZE],EN,ENT;
  double FP[SIZE][SIZE],ONE[4]={1,0,0,1};
  double WORK[26*SIZE];
  double tol[3]={1e-8,1e-8,1e-8},delta=1.0;
  int ISUPPZ[SIZE][SIZE],IWORK[10*SIZE];
  int MAXIT[3]={1000,1000,1000},choice,iter=0;  
  int i,j,k,l,s,M;
  choice=HS->dat.choice;
  s=SIZE;

  for(i=0;i<s;i++)
    for(j=0;j<s;j++)
      HS->dat.P[i][j]=0.0;
  while(iter<MAXIT[choice] && delta > tol[choice]){
    iter++;
    if(IOP<2)
      fprintf(stderr,"::::::::::::::::::::::::>  ITERATION:: %d\n",iter);
    formG(HS);
    if(IOP<2){
      matout(HS->dat.P,"P",s,s);
      matout(HS->dat.G,"G",s,s);
    }
    //create fock matrix
    for(i=0;i<s;i++)
      for(j=0;j<s;j++)
	F[i][j]=HS->dat.G[i][j]+HS->dat.H[i][j];
    EN=0.0;
    for(i=0;i<s;i++)
      for(j=0;j<s;j++)
	EN+=0.5*(HS->dat.P[i][j])*(F[i][j]+HS->dat.H[i][j]);
    if(IOP<2)
      matout(F,"F",SIZE,SIZE);
    fprintf(stderr,"\n:::Energy::: %e\n",EN);
    mult(F,HS->dat.X,HS->dat.G,s);
    mult(HS->dat.XT,HS->dat.G,FPR,s);
    // DIAGONALIZE FPR 
    // EIGENVECT IN CPR 
    // EIGENV. IN E
    diag(FPR,CPR,E);
    /* get the eigenvalues and eigenvectors
    dsyevr('V', 'A', 'U', SIZE, FPR, SIZE, 0, 0, 0, 0, dlamch('S'), &M,
           E, CPR, SIZE, ISUPPZ, WORK, 26*SIZE, IWORK, 10*SIZE);
    */
    mult(HS->dat.X,CPR,C,s);
    for(i=0;i<s;i++)
      for(j=0;j<s;j++){
	OLDP[i][j]=HS->dat.P[i][j];
	HS->dat.P[i][j]=0.0;
	for(k=0;k<1;k++)
	  HS->dat.P[i][j]+=2.0*C[i][k]*C[j][k];
      }
    if(IOP<2){
      matout(FPR,"F'",s,s);
      matout(CPR,"C'",s,s);
      matout(E,"E",s,s);
      matout(C,"C",s,s);
      matout(HS->dat.P,"P",s,s);
    }
    delta=0.0;
    for(i=0;i<s;i++)
      for(j=0;j<s;j++)
	delta+=pow(HS->dat.P[i][j]-OLDP[i][j],2.0);
    delta=sqrt(delta/4.0);
    fprintf(stderr,"CONVERGENCE OF DENSITY: %f\n",delta);
  }
  //CONVERGENCE
  ENT=EN+Z[0]*Z[1]/R;
  fprintf(stderr,"::::::::::::::::::::::::>  FINISHED ITERATIONS\n");
  if(delta > tol[choice])
    fprintf(stderr,"\nWARNING::METHOD DID NOT CONVERGE\n");
  fprintf(stderr,"\nTotal energy:: %e\nElectronic energy:: %e\n\n",ENT,EN);
  if(IOP>1){
     matout(HS->dat.G,"G",SIZE,SIZE);
     matout(F,"F",SIZE,SIZE);
     matout(E,"E",SIZE,SIZE);
     matout(C,"C",SIZE,SIZE);
     matout(HS->dat.P,"P",SIZE,SIZE);
  }
  //MULLIKEN CHARGES
  mult(HS->dat.P,HS->dat.S,OLDP,SIZE);
  if(IOP>0)
    matout(OLDP,"MULLIKEN",SIZE,SIZE);
}

void diag(double F[SIZE][SIZE], double C[SIZE][SIZE], double E[SIZE][SIZE])
{
  double theta,temp,cost,sint;
  int s;
  s=SIZE;

//  if(abs(F[0][0]-F[1][1])>1e-200)
    theta=0.5*atan(2.0*F[0][1]/(F[0][0]-F[1][1]));
//  else
//    theta=PI/4.0;
  cost=cos(theta);
  sint=sin(theta);
  C[0][0]=cost; 
  C[1][0]=sint;
  C[0][1]=sint; 
  C[1][1]=-cost;
  E[0][0]=F[0][0]*SQR(cost)+F[1][1]*SQR(sint)+F[0][1]*sin(2.0*theta);
  E[1][1]=F[1][1]*SQR(cost)+F[0][0]*SQR(sint)-F[0][1]*sin(2.0*theta);
  E[1][0]=0.0;
  E[0][1]=0.0;

  if(E[1][1]<E[0][0]){
    temp=E[1][1];
    E[1][1]=E[0][0];
    E[0][0]=temp;
    temp=C[0][1];
    C[0][1]=C[0][0];
    C[0][0]=temp;
    temp=C[1][1];
    C[1][1]=C[1][0];
    C[1][0]=temp;
  }
  
}

void formG(SYSM HS)
{
  int i,j,k,l,s,s2,s3;
  
  s=SIZE;s2=s*s;s3=s2*s;
  for(i=0;i<s;i++)
    for(j=0;j<s;j++){
      HS->dat.G[i][j]=0.0;
      for(k=0;k<s;k++)
	for(l=0;l<s;l++)
	  HS->dat.G[i][j]+=HS->dat.P[k][l]*(HS->dat.TT[i][j][k][l]-HS->dat.TT[i][l][k][j]*0.5);
    }	  
}

void matout(double A[SIZE][SIZE], char* name, int N, int M)
{
  int i,j;

  //do stuff with matrix here.

  fprintf(stderr,"\nMATRIX::%s\n",name);
  for(i=0;i<N;i++){
    for(j=0;j<M;j++){
      fprintf(stderr,"% 12.5E ",A[i][j]);
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");
}

void vectout(double A[],char* name,int N) // also need vectout
{
  int i,j;

  fprintf(stderr,"\nVECTOR::%s\n",name);
  for(i=0;i<N;i++){
      fprintf(stderr,"% 12.5E ",A[i]);
  }
  fprintf(stderr,"\n");
}

void ttout(double TT[SIZE][SIZE][SIZE][SIZE],char* name, int N)
{
  int i,j,k,l,n1,n2,n3;
  n1=N;n2=n1*n1;n3=n2*n1;

  fprintf(stderr,"\n%s ELEMENTS:\n",name);
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      for(k=0;k<N;k++)
	for(l=0;l<N;l++)
	  fprintf(stderr,"TT(%d,%d,%d,%d) =%12.5E\n",i+1,j+1,k+1,l+1,TT[i][j][k][l]);
  fprintf(stderr,"\n");
}

double S(double A, double B, double RAB2)
{
  double SR=0.0;
  SR=pow(PI/(A+B),1.5)*exp(-A*B*RAB2/(A+B));
  return SR;
}

double T(double A, double B, double RAB2)
{
  double TR=0.0;
  TR=A*B/(A+B)*(3.0-2.0*A*B*RAB2/(A+B))*(pow((PI/(A+B)),1.5))*exp(-A*B*RAB2/(A+B));
  return TR;
}

double V(double A, double B, double RAB2, double RCP2, float Z)
{
  double VR=0.0;
  VR=-2.0*PI/(A+B)*F0((A+B)*RCP2)*exp(-A*B*RAB2/(A+B))*Z;
  if(DEBUG){
    printf("\n========V INTEGRAL:: %f\n%f %f %f %f %f\n---------------------------\n",VR,A,B,RAB2,RCP2,Z);
    printf("F0 PART::%f\n========\n",F0((A+B)*RCP2));
  }
  return VR;
}

double E2I(double A, double B, double C, double D,
	 double RAB2, double RCD2, double RPQ2)
{
  double TWOE;
  TWOE=2.0*pow(PI,2.5)/((A+B)*(C+D)*sqrt(A+B+C+D));
  TWOE=TWOE*F0((A+B)*(C+D)*RPQ2/(A+B+C+D));
  TWOE=TWOE*exp(-A*B*RAB2/(A+B)-C*D*RCD2/(C+D));
  return TWOE;
}

double F0(double ARG)
{
  double FO;
  if(ARG>1e-6){
    FO=sqrt(PI/ARG)*erf(sqrt(ARG))/2.0;
    if(DEBUG)
      printf("F0-PARTS::%f %f\n",sqrt(PI/ARG),erf(sqrt(ARG))/2.0);
  }
  else
    FO=1.0-ARG/3.0;
  return FO;
}

void test(void)
{
  double A[SIZE][SIZE]={0,0,1,1};
  double B[SIZE][SIZE]={1,2,3,4};
  double C[SIZE][SIZE]={0,0,0,0};

  mult(A,B,C,SIZE);
  matout(C,"MULT [[0 0];[1 1]] with [[1 2];[3 4]]",SIZE,SIZE);
}

void mult(double A[SIZE][SIZE], double B[SIZE][SIZE], double C[SIZE][SIZE], int N)
{
  C[0][0]=A[0][0]*B[0][0]+A[0][1]*B[1][0];
  C[0][1]=A[0][0]*B[0][1]+A[0][1]*B[1][1];
  C[1][0]=A[1][0]*B[0][0]+A[1][1]*B[1][0];
  C[1][1]=A[1][0]*B[0][1]+A[1][1]*B[1][1];
}

void copy(double A[SIZE][SIZE], double B[SIZE][SIZE], int N, int M)
{
  int i,j;

  for(i=0;i<N;i++)
    for(j=0;j<M;j++)
      A[i][j]=B[i][j];
}
