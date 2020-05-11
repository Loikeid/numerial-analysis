#include <stdio.h>
#include <math.h>

#define EPSILON 1e-50
#define SECANTRATE 9e-1 //suppose between 0 and 1, though other value do not cause wrong.
#define TIMES 100
#define X_0 3

double Newtoon(double (*Nfunc)(double ),double (*Nderfunc)(double ),double x_0,double epsilon){
    double x=x_0;int i=0;
    printf(" x_0=%.15f\n",x_0);
    do{
        x_0=x;x=x_0-Nfunc(x_0)/Nderfunc(x_0);i+=1;
        printf(" x_%d=%.15f\n",i,x);
    }while(i<TIMES&&fabs(x-x_0)>epsilon);
    return x;
}
double NewtoonSecant(double (*NSfunc)(double),double x_0,double epsilon){
    double secant_rate=SECANTRATE;
    double NSderfunc(double x){
        return (NSfunc(x)-NSfunc(x*secant_rate))/(x*(1-secant_rate));
    }
    return Newtoon(NSfunc,NSderfunc,x_0,epsilon);
}
double func(double );
double derfunc(double);

int main(){
    double x_0;x_0=X_0;
    printf("\n %.15f\n",NewtoonSecant(func,x_0,EPSILON));
    return 0;
}

double func(double x){
    return pow(x,5.0)-225;
}
double derfunc(double x){
    return 5.0*pow(x,4.0);
}
