/*Compile with -Ofast flag*/
#include <iostream>
#include <cmath>
#include <cstdlib>

const double DX = 1.0/10000.; //Grid with denominator cells
const double LX = 1.0;
const int NX =  LX/DX + 1;
const double C = 0.3; //Courant Number
const double v = 374.17;
const double DT = DX*C/v;
const double g = 1.4;
const double A = 0.5; //Point of the discontinuity
const double T = 0.2*0.65 ; // Real time of the simulation
/*------Initial data-----*/
const double rhoL = 1.0; 
const double uL = 0.75;
const double pL = 1.0;

const double rhoD = 0.125;
const double uD = 0.0; 
const double pD = 0.1;
/*---------------------*/
	
class roeNewtonian
{
public:

	double u[NX];
	double rho[NX];
	double p[NX];
	double E[NX];
	double H[NX];
	double e[NX];
	double c[NX]; //Local speed of sound
	double mach[NX]; //Mach number
	
	/*-------------------*/
	double uAve;
	double HAve;
	double aAve;
	double uAveMinus;
	double HAveMinus;
	double aAveMinus;
	/*Average eigenvalue*/
	double lambda1;
	double lambda2;
	double lambda3;
	double lambdaMinus1;
	double lambdaMinus2;
	double lambdaMinus3;
	/*-------------------*/
	
	double deltaU1;
	double deltaU2;
	double deltaU3;
	double deltaU1Minus;
	double deltaU2Minus;
	double deltaU3Minus;

	double sigmaRho[NX];
	double sigmaRhoU[NX];
	double sigmaE[NX];
	double sigma[NX];

	void begin (void);
	void boundaries (void);
	void averageValues (void);
	void averageEigenvalues (void);
	void finiteVolume (void);
	void print(void);
	void gnuplot(void);	
	void minmod (void);

};

void roeNewtonian::begin (void){ 
	for (int i = 1; i<NX*A; i++ ){
		u[i]=  uL;		rho[i]=  rhoL;
		p[i]=  pL;		E[i]= p[i]/(g-1.) + 0.5*rho[i]*pow(u[i],2);
		H[i]=(E[i]+p[i])/rho[i]; 
		e[i]=p[i]/((g-1.)*rho[i]);
		c[i]=sqrt(g*p[i]/rho[i]);
		mach[i]=u[i]/c[i];
	}

    for (int i = NX*A; i<NX-1; i++ ){
		u[i]=  uD;		rho[i]=  rhoD;
		p[i]= pD;		E[i] =  p[i]/(g-1.) + 0.5*rho[i]*pow(u[i],2);
		H[i]=(E[i]+p[i])/rho[i];
		e[i]=p[i]/((g-1.)*rho[i]);
		c[i]=sqrt(g*p[i]/rho[i]);
		mach[i]=u[i]/c[i];
	}
	for (int i = 1; i < NX-1; ++i){
		sigmaRho[i] = rho[i];
		sigmaRhoU[i] = u[i]*rho[i];
		sigmaE[i] = E[i];
	}
}

void roeNewtonian::boundaries (void){
	u[0]= uL;	u[NX-1]= uD;  
	rho[0]=  rhoL;	rho[NX-1]= rhoD;
	p[0]=pL;	p[NX-1]=pD;
	e[0]=p[0]/((g-1.)*rho[0]);	e[NX-1]=p[NX-1]/((g-1.)*rho[NX-1]);	
	E[0]= p[0]/(g-1.) + 0.5*rho[0]*pow(u[0],2);	E[NX-1]= p[NX-1]/(g-1.) + 0.5*rho[NX-1]*pow(u[NX-1],2);
	H[0]=(E[0]+p[0])/rho[0];	H[NX-1]=(E[NX-1]+p[NX-1])/rho[NX-1];
	c[0]=sqrt(g*p[0]/rho[0]);	c[NX-1]=sqrt(g*p[NX-1]/rho[NX-1]);
	mach[0]=u[0]/c[0];	mach[NX-1]=u[NX-1]/c[NX-1];
}


void roeNewtonian::averageValues (void){
	/*u Average*/
	double UL, ULPlus, ULMinus;
	double rhoSquareL, rhoSquareLPlus,rhoSquareLMinus;
	for(int ii=1; ii<NX*A; ii++){
		UL+=u[ii];
		ULPlus+=u[ii+1];
		ULMinus+=u[ii-1];
		rhoSquareL+=sqrt(rho[ii]);
		rhoSquareLPlus+=sqrt(rho[ii+1]);
		rhoSquareLMinus+=sqrt(rho[ii-1]);
	}
		
	 UL=UL/(NX*A);	ULPlus=ULPlus/(NX*A);	ULMinus=ULMinus/(NX*A);
	 rhoSquareL=rhoSquareL/(NX*A);	rhoSquareLPlus=rhoSquareLPlus/(NX*A);	rhoSquareLMinus=rhoSquareLMinus/(NX*A);

	double UR, URMinus, URPlus;
	double rhoSquareR, rhoSquareRPlus,	rhoSquareRMinus;
	for(int ii=NX*A; ii<NX-1; ii++){
		UR+=u[ii];
		URPlus+=u[ii+1];
		URMinus+=u[ii-1];
		rhoSquareR+=sqrt(rho[ii]);
		rhoSquareRPlus+=sqrt(rho[ii+1]);
		rhoSquareRMinus+=sqrt(rho[ii-1]);
	}
	URMinus=URMinus/(NX*(1.-A));
	URPlus=URPlus/(NX*(1.-A));
	UR=UR/(NX*(1.-A));
	rhoSquareR=rhoSquareR/(NX*(1.-A));
	rhoSquareRPlus=rhoSquareRPlus/(NX*(1.-A));
	rhoSquareRMinus=rhoSquareRMinus/(NX*(1.-A));

	uAve = (rhoSquareL*UL+rhoSquareR*UR)/(rhoSquareR+rhoSquareL);
	
	uAveMinus = (rhoSquareLMinus*ULMinus+rhoSquareRMinus*URMinus)/(rhoSquareRMinus+rhoSquareLMinus);
	/*H Average*/
	double HL,	HLMinus;
	double HR,	HRMinus;
	for(int ii=1; ii<NX*A; ii++){
		HL+=H[ii+1];
		HLMinus+=H[ii-1];
	}
	for(int ii=NX*A; ii<NX-1; ii++){
		HR+=H[ii+1];
		HRMinus+=H[ii-1];
	}
	HL=HL/(NX*A);
	HLMinus=HLMinus/(NX*A);
	HR=HR/(NX*(1.-A));
	HRMinus=HRMinus/(NX*(1.-A));
	HAve=(rhoSquareL*HL+rhoSquareR*HR)/(rhoSquareR+rhoSquareL);
	HAveMinus=(rhoSquareLMinus*HLMinus+rhoSquareRMinus*HRMinus)/(rhoSquareRMinus+rhoSquareLMinus);
	/*a Average*/
	aAve=sqrt((g-1.)*(HAve-0.5*pow(uAve,2))); 
	uAveMinus=sqrt((g-1.)*(HAveMinus-0.5*pow(uAveMinus,2)));
	/*Deltas in order to find alphas*/
	deltaU1=pow(rhoSquareR,2)-pow(rhoSquareL,2);
	deltaU1Minus=pow(rhoSquareRMinus,2)-pow(rhoSquareLMinus,2);
	deltaU2=pow(rhoSquareR,1)*UR-pow(rhoSquareL,1)*UL;
	deltaU2Minus=pow(rhoSquareRMinus,1)*URMinus-pow(rhoSquareLMinus,1)*ULMinus;
	deltaU3=pow(rhoSquareR,2)*HR - pow(rhoSquareL,2)*HL;
	deltaU3Minus=pow(rhoSquareRMinus,2)*HRMinus - pow(rhoSquareLMinus,2)*HLMinus;
	
}
void roeNewtonian::averageEigenvalues(void){
	lambda1=abs(uAve - aAve);
	lambda2=abs(uAve);
	lambda3=abs(uAve + aAve);
	lambdaMinus1=abs(uAveMinus - aAveMinus);
	lambdaMinus2=abs(uAveMinus);
	lambdaMinus3=abs(uAveMinus + aAveMinus);
}

void roeNewtonian::finiteVolume (void){
	double rhoNew[NX] = {};	double uNew[NX] = {};
	double pNew[NX] = {};	double ENew[NX] = {}; 
	double lambda=0.5*DT/DX;
	double lambdaMayor = lambda1+lambda2+lambda3;
	double lambdaMenor = lambdaMinus1+lambdaMinus2+lambdaMinus3;

	for (int i = 0; i < NX; ++i)		rhoNew[i] = rho[i];
	for (int i = 0; i < NX; ++i)		uNew[i] = u[i];
	for (int i = 0; i < NX; ++i)		pNew[i] = p[i];
	for (int i = 0; i < NX; ++i)		ENew[i] = E[i];
	 
	for (int ii = 1; ii < NX-1; ++ii)
		rho[ii] = rhoNew[ii] + lambda*(rhoNew[ii-1]*uNew[ii-1] + (rhoNew[ii+1]-rhoNew[ii])*lambdaMayor +(rhoNew[ii-1]-rhoNew[ii])*lambdaMenor
			 - rhoNew[ii+1]*uNew[ii+1])-(0.5*lambda)*(DX-DT)*(sigmaRho[ii]-sigmaRho[ii-1]);

	for (int ii = 1; ii < NX-1; ++ii)
		u[ii] = rhoNew[ii]*uNew[ii] + lambda*(rhoNew[ii-1]*pow(uNew[ii-1],2)+pNew[ii-1] 
			+(rhoNew[ii+1]*uNew[ii+1]-rhoNew[ii]*uNew[ii])*lambdaMayor +(-rhoNew[ii]*uNew[ii]+rhoNew[ii-1]*uNew[ii-1])*lambdaMenor
			-(rhoNew[ii+1]*pow(uNew[ii+1],2)+pNew[ii+1]))- (0.5*DT/DX)*(DX-DT)*(sigmaRhoU[ii]-sigmaRhoU[ii-1]);
	for (int ii = 0; ii < NX; ++ii)
		u[ii] = u[ii]/rho[ii];		
	for (int ii = 1; ii < NX-1; ++ii)
		E[ii] = ENew[ii] + lambda*(uNew[ii-1]*(ENew[ii-1]+pNew[ii-1]) +(ENew[ii+1]-ENew[ii])*lambdaMayor+(-ENew[ii]+ENew[ii-1])*lambdaMenor 
			-(uNew[ii+1]*(ENew[ii+1]+pNew[ii+1])))- (0.5*lambda)*(DX-DT)*(sigmaE[ii]-sigmaE[ii-1]);		
	for (int ii = 0; ii < NX; ++ii)
		p[ii]=(E[ii]-0.5*rho[ii]*pow(u[ii],2))*(g-1.);
	for (int ii = 0; ii < NX; ++ii)	
		e[ii]=p[ii]/((g-1.)*rho[ii]);	
	for (int i = 0; i < NX; ++i)
		c[i]=sqrt(g*p[i]/rho[i]);
	for (int i = 0; i < NX; ++i)
		mach[i]=u[i]/c[i];
	for (int i = 0; i < NX; ++i)	H[i] = (E[i] + p[i])/rho[i];

}	


void roeNewtonian::print(void)
{
 //std::cout << "plot '-' with lines linecolor rgb 'black' lw 5 " << std::endl;
 //std::cout << "plot '-' with points pt 7 lc rgb 1 using xticlabels using yticlabels" << std::endl;
    std::cout << "plot '-' w l lw 2 " << std::endl;
    double x;
    for (int ii = 0; ii < NX-1; ++ii)
    {
        x = ii * DX;
        std::cout << x << "  " << p[ii]<<  std::endl;
        //std::cout << x << "  " << rho[ii]/10.<< "  " << p[ii]/1000.<< "  " << u[ii]<< std::endl;
    }
    std::cout << "e" << std::endl;
   //std::cout << "set yrange [-0.5:1.5] " << std::endl;
}

void roeNewtonian::gnuplot(void)
{
    std::cout << "set terminal pdf" << std::endl;
    std::cout << "unset key" << std::endl;
    //std::cout << "set xtics ( 0.0, 0.2, 0.4, 0.6, 0.8, 1.0)" << std::endl;  
    //std::cout << "set ytics ( 0.0, 0.4, 0.6, 0.8, 1.0)" << std::endl;
    //std::cout << "set ylabel 'Density'" << std::endl;
    //std::cout << "set xlabel 'Position'" << std::endl;
   // std::cout << "set size  1,1" << std::endl;
    //std::cout << "set yrange[0.0:1.0]" << std::endl;
    //std::cout << " set style line 7 lc rgb '#bf000a' lt 1 lw 5.5 " << std::endl;

    //std::cout << "set out 'anyTest.pdf'" << std::endl;
 	std::cout << "set contour base" << std::endl;
  	//std::cout << "set pm3d" << std::endl;
   // std::cout << "set xlabel 'x' " << std::endl;

	//std::cout << "set terminal gif animate" << std::endl;
    std::cout << "set out 'example1.pdf'" << std::endl;
    //std::cout << "set contour base" << std::endl;
    
}	

void roeNewtonian::minmod (void){
	
	double alpha[NX] = {};
	double beta[NX] = {};
	for (int i = 0; i < NX; ++i)
		alpha[i] = (u[i]-u[i-1])/DX;
	for (int i = 0; i < NX-1; ++i)
		beta[i] = (u[i+1]-u[i])/DX;	
	for (int i = 0; i < NX-1; ++i)
	{	
		if (alpha[i]*beta[i]<=0.)
			sigma[i] = 0.;
			else{

				if (abs(alpha[i]) < abs(beta[i]))
					sigma[i] = alpha[i];
				if (abs(alpha[i]) > abs(beta[i]))
					sigma[i] = beta[i];
		}
	}
}

int main(void)
{
	roeNewtonian roeNewtonian;
    int nit;
    nit=T/DT;
    roeNewtonian.begin();
    roeNewtonian.boundaries();
    roeNewtonian.gnuplot();
    for (int t = 0; t <= nit; t++)
    {	
    	roeNewtonian.averageValues();
    	roeNewtonian.averageEigenvalues();
       	roeNewtonian.minmod();
    	roeNewtonian.finiteVolume();	
        roeNewtonian.boundaries();
        //roeNewtonian.print();
    }
    roeNewtonian.print();
    return 0;
}
