#include<stdio.h>

#define TIMES 100
#define EPSILON 1e-10

//assist functions
double abs(double x)
{
    return (x>0)?x:-x;
}
void output(double *x, int n, int l)
{
    int i=0,col;col=(n%l==0)?(n/l):(n/l+1);
    for(;i<n;i++)
        if((i+1)%col!=0)printf("%f\t",*(x+i));
        else printf("%f\n",*(x+i));
    printf("\n");
}

// algorithms
void Jacobi(double A[], double b[], double x[],int n)
{
    int i,j,k=0;double temp,norm=0,norm1,x1[n];
    do{
        norm1=norm;
        for(i=0;i<n;i++){temp=0;
            for(j=0;j<n;j++){
                temp+=A[i*n+j]*x[j]*(i!=j);
            }
            x1[i]=1/A[i*n+i]*(-temp+b[i]);
        }norm=0;
        for(j=0;j<n;j++){
            norm+=abs(x1[j]);
            x[j]=x1[j];
        }
        k++;output(x,n,n);
    }while(k<TIMES&&abs(norm-norm1)>EPSILON);
    printf("\n");
}
void Gauss_Seidel(double A[], double b[], double x[],int n)
{
    int i,j,k=0;double temp,norm=0,norm1;
    do{
        norm1=norm;
        for(i=0;i<n;i++){temp=0;
            for(j=0;j<n;j++){
                temp+=A[i*n+j]*x[j]*(i!=j);
            }
            x[i]=1/A[i*n+i]*(-temp+b[i]);
        }norm=0;
        for(j=0;j<n;j++)norm+=abs(x[j]);
        k++;output(x,n,n);
    }while(k<TIMES&&abs(norm-norm1)>EPSILON);
}
void Gauss(double A_bar[], int n)
{
    double l;int i=0,j,k,col=n+1;
    printf("原矩阵为：\n");
    output(A_bar,n*col,n);
    printf("消元所得结果：\n");
    for(;i<n;i++){
        for(j=i+1;j<n;j++){
            l=-A_bar[j*col+i]/A_bar[i*col+i];
            for(k=i;k<col;k++)
                A_bar[j*col+k]+=l*A_bar[i*col+k];
        }
    }
    output(A_bar,n*col,n);
    printf("回代所的结果：\n");
    for(i=n-1;i>0;i--){
        for(j=i-1;j>=0;j--){
            l=-A_bar[j*col+i]/A_bar[i*col+i];
            A_bar[j*col+n]+=l*A_bar[i*col+n];
            A_bar[j*col+i]=0;
        }
    }output(A_bar,n*col,n);
    printf("最终结果：\n");
    for(i=0;i<n;i++){
        A_bar[i*col+n]/=A_bar[i*col+i];
        A_bar[i*col+i]=1;
    }output(A_bar,n*col,n);
}
void Gauss_Jordan(double A_bar[], int n)
{
    double l;int i=0,j,k,col=n+1;
    printf("原矩阵为：\n");
    output(A_bar,n*col,n);
    printf("消元所得结果：\n");
    for(;i<n;i++){
        for(j=0;j<n;j++){
            if(j==i)continue;
            l=-A_bar[j*col+i]/A_bar[i*col+i];
            for(k=i;k<col;k++)
                A_bar[j*col+k]+=l*A_bar[i*col+k];
        }
    }
    output(A_bar,n*col,n);
    printf("最终结果：\n");
    for(i=0;i<n;i++){
        A_bar[i*col+n]/=A_bar[i*col+i];
        A_bar[i*col+i]=1;
    }output(A_bar,n*col,n);
}
void Gauss_principal(double A_bar[], int n)
{
    double l;int i=0,j,k,p_index,col=n+1;
    printf("原矩阵为：\n");
    output(A_bar,n*col,n);
    for(;i<n;i++){
        p_index=i;printf("交换主元：\n");
        for(j=i;j<n;j++){
            p_index=(abs(A_bar[j*col+i])>abs(A_bar[p_index*col+i]))?j:p_index;
        }for(k=i;k<col;k++){
            l=A_bar[i*col+k];
            A_bar[i*col+k]=A_bar[p_index*col+k];
            A_bar[p_index*col+k]=l;
        }output(A_bar,n*col,n);
        for(j=i+1;j<n;j++){
            l=-A_bar[j*col+i]/A_bar[i*col+i];
            for(k=i;k<col;k++)
                A_bar[j*col+k]+=l*A_bar[i*col+k];
        }
    }
    output(A_bar,n*col,n);
    printf("回代所的结果：\n");
    for(i=n-1;i>0;i--){
        for(j=i-1;j>=0;j--){
            l=-A_bar[j*col+i]/A_bar[i*col+i];
            A_bar[j*col+n]+=l*A_bar[i*col+n];
            A_bar[j*col+i]=0;
        }
    }output(A_bar,n*col,n);
    printf("最终结果：\n");
    for(i=0;i<n;i++){
        A_bar[i*col+n]/=A_bar[i*col+i];
        A_bar[i*col+i]=1;
    }output(A_bar,n*col,n);
}
double Doolittle(double A[],int n)
{
    int i,j,k,m,s=n*n;double d=1;
    for(k=0;k<n;k++){
        for(j=k;j<n;j++){
            for(m=0;m<k;m++){
                A[k*n+j]-=A[k*n+m]*A[m*n+j];
            }
        }
        if(A[k*n+k]==0)return 0;
        for(i=k+1;i<n;i++){
            for(m=0;m<k;m++){
                A[i*n+k]-=A[i*n+m]*A[m*n+k];
            }
            A[i*n+k]/=A[k*n+k];
        }
        d*=A[k*n+k];
    }output(A,s,n);
    return d;
}

// main function
int main()
{
    double A[9]={1,-2,2,-1,1,-1,-2,-2,1};
    return Doolittle(A,3);
}
