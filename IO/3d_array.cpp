#include<iostream>
#include<cmath>
#include<time.h>
#define m 500

using namespace std;

int main(){
  
double*** arr=new double** [m];
for(int i=0;i<m;i++) arr[i]=new double* [m];
for(int i=0;i<m;i++){
	for(int j=0;j<m;j++){
		arr[i][j]=new double [m];
	}
}


  clock_t t0,t1,t2,t3,t4,t5;
   
  t0=clock();
  for(int i=0;i<m;i++){
    for(int j=0;j<m;j++){
      for(int k=0;k<m;k++){
	
	arr[i][j][k]=i+j+k;
	//cout<<arr[i][j][k]<<endl;
      }
    }
 }
 
  t0=clock()-t0;
  cout<<"i,j,k :"<<(float(t0))/CLOCKS_PER_SEC<<endl; 
 
 t1=clock(); 
 for(int i=0;i<m;i++){
    for(int k=0;k<m;k++){
      for(int j=0;j<m;j++){
	
	arr[i][j][k]=i+j+k;
      }
    }
 }
  t1=clock()-t1;
  cout<<"i,k,j :"<<(float(t1))/CLOCKS_PER_SEC<<endl;
  
  t2=clock();
  for(int k=0;k<m;k++){
    for(int j=0;j<m;j++){
      for(int i=0;i<m;i++){
	
	arr[i][j][k]=i+j+k;
      }
    }
 }
 
 t2=clock()-t2;
 cout<<"k,j,i :"<<(float(t2))/CLOCKS_PER_SEC<<endl;
 
  t3=clock();
  for(int k=0;k<m;k++){
    for(int i=0;i<m;i++){
      for(int j=0;j<m;j++){
	
	arr[i][j][k]=i+j+k;
      }
    }
 }
 
 t3=clock()-t3;
 cout<<"k,i,j :"<<(float(t3))/CLOCKS_PER_SEC<<endl; 
 
  t4=clock();
  for(int j=0;j<m;j++){
    for(int i=0;i<m;i++){
      for(int k=0;k<m;k++){
	
	arr[i][j][k]=i+j+k;
      }
    }
 }
 
 t4=clock()-t4;
 cout<<"j,i,k :"<<(float(t4))/CLOCKS_PER_SEC<<endl;   
 
  t5=clock();
  for(int j=0;j<m;j++){
    for(int k=0;k<m;k++){
      for(int i=0;i<m;i++){
	
	arr[i][j][k]=i+j+k;
      }
    }
 }
 
 t5=clock()-t5;
 cout<<"j,k,i :"<<(float(t5))/CLOCKS_PER_SEC<<endl;    
  
 return 0; 
}

/* m = 500
   i,j,k :1.23416
   i,k,j :2.24749
   k,j,i :7.57492
   k,i,j :5.92052
   j,i,k :1.25228
   j,k,i :5.16701
*/