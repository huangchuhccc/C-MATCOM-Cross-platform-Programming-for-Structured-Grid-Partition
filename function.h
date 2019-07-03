#include <iostream>
#include "matlib.h"
#include <cmath>
using namespace std;
double wing(double x); //计算翼形在对应点的函数值
//Mm bound(Mm& sita); //计算中心射线与边界的交点
double getborder(double sita);
void coefficence(double& a,double& b,double& c,double i,double j,Mm& x, Mm& y,double dx,double dy);
void Cacl_p(Mm& x,Mm& y ,Mm& x_temp,Mm& y_temp,double dx,double dy,int& NUM,double& omega);
void Cacl_l(Mm& x,Mm& y ,Mm& x_temp,Mm& y_temp,double dx,double dy,int NUM,double omega);


double wing(double x){    //计算翼型的边界方程
return 0.2969*sqrt(x)-0.1260*x-0.3516*pow(x,2)+0.2843*pow(x,3)-0.1015*pow(x,4);
}

double getborder(double sita){   
double x[3], y[3]; 
x[0]=0.5;  //利用二分法，计算射线与翼型边界的交点
if (cos(sita)>0)
x[2]=1;
else
x[2]=0;
while (abs(x[0]-x[2])>1e-6)
{
x[1]=(x[0]+x[2])/2;
if (sin(sita)>=0)
y[2]=wing(x[2])-tan(sita) * (x[2]-0.5),
y[1]=wing(x[1])-tan(sita) * (x[1]-0.5),
y[0]=wing(x[0])-tan(sita) * (x[0]-0.5);
else
y[2]=-wing(x[2])-tan(sita) * (x[2]-0.5),
y[1]=-wing(x[1])-tan(sita) * (x[1]-0.5),
y[0]=-wing(x[0])-tan(sita) * (x[0]-0.5);
if (y[1]*y[0]<0)x[2]=x[1];
else
x[0]=x[1];
}
return x[1];
}

void coefficence(double& a,double& b,double& c,double i,double j,Mm& x, Mm& y,double dx,double dy){
double x_yita,y_yita,x_yip,y_yip;
int i0=i,i1=i;
if (i0==1)
{
i0=x.size(1)-1;
}
else
i0=i-1;
if (i1==x.size(1))
{
i1=2;
}
else
i1=i+1;
x_yita=(x.r(i,j+1)-x.r(i,j-1))/(2*dy);
x_yip =(x.r(i1,j)-x.r(i0,j))/(2*dx);
y_yita=(y.r(i,j+1)-y.r(i,j-1))/(2*dy);
y_yip =(y.r(i1,j)-y.r(i0,j))/(2*dx);
a=pow(x_yita,2)+pow(y_yita,2);
b=x_yip*x_yita+y_yip*y_yita;
c=pow(x_yip,2)+pow(y_yip,2);
}
void Cacl_p(Mm& x,Mm& y ,Mm& x_temp,Mm& y_temp,double dx,double dy,int& NUM,double& omega){ //Gauss-Siedel点迭代。充分利用新值
double a,b,c;
for (int i=2;i<=x.size(1);i++)
{
for (int j=2;j<=x.size(2)-1;j++)
{
int i0=i-1,i1=i;
if (i1==x.size(1))
{
i1=2;
}
else
i1=i+1;
coefficence(a,b,c,i,j,x,y,dx,dy);
x.r(i,j)=(a/(dx*dx)*(x.r(i1,j)+x.r(i0,j)) - b/(2*dx*dy)*(x.r(i1,j+1)-x.r(i0,j+1)-x.r(i1,j-1)+x.r(i0,j-1))+c/(dy*dy)*(x.r(i,j+1)+x.r(i,j-1)))/(2*a/(dx*dx)+2*c/(dy*dy));
y.r(i,j)=(a/(dx*dx)*(y.r(i1,j)+y.r(i0,j)) - b/(2*dx*dy)*(y.r(i1,j+1)-y.r(i0,j+1)-y.r(i1,j-1)+y.r(i0,j-1))+c/(dy*dy)*(y.r(i,j+1)+y.r(i,j-1)))/(2*a/(dx*dx)+2*c/(dy*dy));
if (NUM==3) //超松弛法点迭代
{
x.r(i,j)=x_temp.r(i,j)+omega*(x.r(i,j)-x_temp.r(i,j));
y.r(i,j)=y_temp.r(i,j)+omega*(y.r(i,j)-y_temp.r(i,j));
}
}
}
x(1,c_p)=x(x.size(1),c_p);
y(1,c_p)=y(x.size(1),c_p);
}
void Cacl_l(Mm& x,Mm& y ,Mm& x_temp,Mm& y_temp,double dx,double dy,int NUM,double omega){
int JM=x.size(2),IM=x.size(1);
double a,b,c;
Mm A ,bx, by;
for (int i=2;i<=IM;i++)
{
	Mm BX,BY;
A=diag(2*ones(1,JM-2))-diag(ones(1,JM-3),-1)-diag(ones(1,JM-3),1);
BX=zeros(JM-2,1),BY=zeros(JM-2,1);
for (int j=2;j<=JM-1;j++)
{
int i1=i;
if (i1==IM)
{
i1=2;
}
else
i1=i+1;
coefficence(a,b,c,i,j,x,y,dx,dy);
BX.r(j-1,1)=(4*a*dy*dy*(x.r(i1,j)+x.r(i-1,j))-2*b*dx*dy*(x.r(i1,j+1)-x.r(i1,j+1)-x.r(i1,j-1)+x.r(i-1,j-1)))/(4*dx*dx*c);
BY.r(j-1,1)=(4*a*dy*dy*(y.r(i1,j)+y.r(i-1,j))-2*b*dx*dy*(y.r(i1,j+1)-y.r(i1,j+1)-y.r(i1,j-1)+y.r(i-1,j-1)))/(4*dx*dx*c);
A.r(j-1,j-1)=A.r(j-1,j-1)+pow(dy,2)/c*a/pow(dx,2)*2;
}
BX.r(1,1)=BX.r(1,1)+x.r(i,1);
BX.r(JM-2,1)=BX.r(JM-2,1)+x.r(i,JM);
BY.r(1,1)=BY.r(1,1)+y.r(i,1);
BY.r(JM-2,1)=BY.r(JM-2,1)+y.r(i,JM);
x(i,colon(2,JM-1))=mldivide(A,BX);
y(i,colon(2,JM-1))=mldivide(A,BY);
}
x(1,c_p)=x(IM,c_p);
y(1,c_p)=y(IM,c_p);
}

