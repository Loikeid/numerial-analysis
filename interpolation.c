#include <stdio.h>

double L(double t, double x[], double y[], int n)
{
    int i,j;double basic,result=0;
    for(j=0;j<n;j++){
        basic=1;
        for(i=0;i<n;i++){
            basic*=(j==i)?1:(t-x[i])/(x[j]-x[i]);
        }
        result+=y[j]*basic;
    }return result;
}

double difference(double a[], double b[], int n)
{
    if(n<=0)return 0xFFFFFFFF;
    else if (n==1)return b[0];
    else if(n==2)return (b[0]-b[1])/(a[0]-a[1]);
    else return (difference(a,b,n-1)-difference(a+1,b+1,n-1))/(a[0]-a[n-1]);
}

double N(double t, double x[], double y[], int n)
{
    int i,j;double temp,result=0;
    for(i=1;i<n+1;i++){
        temp=difference(x,y,i);
        for(j=0;j<i-1;j++){
            temp*=t-x[j];
        }result+=temp;
    }return result;
}

int main()
{
    double x[3]={0,4,6},y[3]={1,3,2},t=2;
    return 0;
}
