#define N 5
#include <stdio.h>

void input(float * a,int n){
    int i;
    for(i=0;i<n;i++){
        scanf("%f",a+i);
    }
}

void getPolynomial(float * a, int n){
    int i;
    printf("\n\tp(x)=");
    for(i=0;i<n-1;i++){
        if (*(a+i)!=1)printf("%gx^%d+",*(a+i),n-1-i);
        else printf("x^%d+",n-i-1);
    }
    printf("%g\n\n",*(a+n-1));
}

int main(){
	int k=1;
	float a[N+1],x;
	printf("Please input %d values as 'a' array:\n",N+1);
	input(a,N+1);
	getPolynomial(a,N+1);
	printf("Please input compute value: ");
	input(&x,1);
	float s=a[0];
	for(;k<=N;k++)s=s*x+a[k];
	printf("\n\n\tp(%g)=%g\n",x,s);
	return 0;
}
