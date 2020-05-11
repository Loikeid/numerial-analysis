#include <stdio.h>
#define INF 1e20
#define EPSILON 1e-50
#define TIMES 100

double fabs(double x)
{
    return (x>0)?x:-x;
}
double Max(double *array,int n){
    int i=0;double max=-INF;
    for(;i<n;i++){
        max=(max<*(array+i))?*(array+i):max;
    }return max;
}
double Min(double *array,int n){
    int i=0;double min=INF;
    for(;i<n;i++){
        min=(min>*(array+i))?*(array+i):min;
    }return min;
}
double sum(double list[],double weight[],int n){
    double result=0.0;int i=0;
    for(;i<n;i++){
        result+=list[i]*weight[i];
    }return result;
}
double powerEigenValue(double *A, double *x,int n, double epsilon){
    int i=0,j=0;double max=0,max2=0,y[n];
    do{
        for(i=0;i<n;i++){
            double nA[n];int k;
            for(k=0;k<n;k++)nA[k]=A[i*n+k];
            y[i]=sum(x,nA,n);
        }j++;
        printf(" &y_{%d}=(%.4f,%.4f,%.4f)&\\quad",j,y[0],y[1],y[2]);
        max2=max;max=Max(y,n);
        for(i=0;i<n;i++){
            x[i]=y[i]/max;
        }
        printf("\\lambda_{%d}=%.4f\\quad&x_{%d}=(%.4f,%.4f,%.4f)\\\\\n",j,max,j,x[0],x[1],x[2]);
    }while(j<TIMES&&fabs(max-max2)>epsilon);
    return max;
}
int main(){
    double A[9]={1,2,0,2,-1,1,0,1,3};
    double x[3]={1,0,0}; double lambda;
    lambda=powerEigenValue(A,x,3,EPSILON);
    printf(" %g\n",lambda);
    return 0;
}
