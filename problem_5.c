//Newton-Raphson method for non linear system of equations

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void function_value(double x1,double x2);
void tangent_matrix(double x1,double x2);
void print_output(int n,double x1,double x2,double R);
void solving_U(void);
double norm_R(double x01,double x02,double x1,double x2);

int main()
{
	int i,n,n_max;
	double x1,x2,x01,x02,fx1,fx2,TOL,R,U[2];
	
	x1=0.5;x2=0.5;	//setting initial values
	n=0;n_max=25;TOL=0.0000001;R=1;	//setting tolerance and maximum iterations
	
	function_value(x1,x2);
	print_output(n,x1,x2,R);
	
	while((R>TOL)&&(n<n_max))
	{
		tangent_matrix(x1,x2);	//to find tangent matrix
		function_value(x1,x2);	//to find fuction values
		solving_U();	//solving U
		{
			FILE *fp;
			fp=fopen("U_matrix.txt","r");
			for(i=0;i<2;i++)
			fscanf(fp,"%lf",&U[i]);
			fclose(fp);	
		}
		//updating xi values
		x01=x1;x02=x2;
		x1=x1+U[0];
		x2=x2+U[1];
		R=norm_R(x01,x02,x1,x2); //error
		n++;
		print_output(n,x1,x2,R);
	}
	
	if(n<n_max)
	printf("The point of intersection of circle and parabola\n (x1,y1)=(%4.7lf,%4.7lf) and (x2,y2)=(%4.7lf,%4.7lf)\n",x1,x2,-x1,x2);
	else
	printf("Newton Raphson method did not converge\n");
	return 0;
}

void function_value(double x1,double x2)	//function to find function values
{
	double fx1,fx2,F[2];
	int i;
	fx1=x1*x1+x2*x2-1;
	fx2=x1*x1-x2;
	F[0]=-fx1;F[1]=-fx2;
	
	FILE *fp;
	fp=fopen("func_value.txt","w");
	for(i=0;i<2;i++)
	fprintf(fp,"%4.7lf\t",F[i]);
	fclose(fp);
}

void tangent_matrix(double x1,double x2)	//for finding tangent matrix
{
	double dfx1,dfx2,K[2][2];
	int i,j;
	dfx1=2*x1;
	dfx2=2*x2;
	for(i=0;i<2;i++)
	{
		K[i][0]=dfx1;
	}
	K[0][1]=dfx2;K[1][1]=-1;
	
	FILE *fp;
	fp=fopen("jacobian.txt","w");
	for(i=0;i<2;i++)
	{
		for(j=0;j<2;j++)
		{
			fprintf(fp,"%4.7lf\t",K[i][j]);
		}
		fprintf(fp,"\n");
	}
	
	fclose(fp);
}

double norm_R(double x01,double x02,double x1,double x2)	//function for finding norm
{
	double a,b,d,norm;
	a=x1-x01;
	b=x2-x02;
	d=a*a+b*b;
	norm= sqrt(d);
	return(norm);
}
	


void print_output(int n,double x1,double x2,double R)	//output printing function
{
	FILE *fp;
	fp=fopen("output_problem_5.txt","a");
	fprintf(fp,"\n k=%d \t%4.7lf \t%4.7lf \t%4.7lf",n,x1,x2,R);
	fclose(fp);
}

void solving_U(void)	//function for solving U using gauss elimination
{
	int i,j,k,l,q,n=2;
	double c,sum,temp;
	double mat[n][n+1],U[n];
	
	FILE *f1,*f2;
	f1 = fopen("jacobian.txt","r");
	for(i=0;i<2;i++)
	{
		for(j=0;j<2;j++)
		{
			fscanf(f1,"%lf",&mat[i][j]);
		}
	
	}
	fclose(f1);
	
	f2= fopen("func_value.txt","r");
	for(i=0;i<2;i++)
	{
		fscanf(f1,"%lf",&mat[i][n]);
	}
	fclose(f1);
	
	
	for(j=0;j<n;j++)
	{
	
			if(mat[j][j]==0)	//swaping rows having pivot element equals zero
			{
				
				for(l=j+1;l<n;l++)
				{
					if(mat[l][j]!=0)
					{
						for(q=j;q<n+1;q++)
							{
								temp=mat[j][q];
								mat[j][q]=mat[l][q];
								mat[l][q]=temp;
							}
						break;
					}
					else
					{
						if(j<n)
						{
							j=j+1;
						}
			
						else
						{
							printf("matrix is singular\n");
							exit(0);
						}
					}
				}
			}
	
		//elimination step
		for(i=0;i<n;i++)
		{
			if(i>j)
			{
				c=mat[i][j]/mat[j][j];
				for(k=0;k<n+1;k++)
				{
					mat[i][k]=mat[i][k]-c*mat[j][k];
				}
			}
		}
	
	}

	//back substitution phase
	U[n-1]=mat[n-1][n]/mat[n-1][n-1];
	for(i=n-2;i>=0;i--)
	{
		sum=0;
		for(j=i+1;j<n;j++)
		{
			sum=sum+mat[i][j]*U[j];
		}
		U[i]=(mat[i][n]-sum)/mat[i][i];
	}
	
	FILE *fp;
	fp=fopen("U_matrix.txt","w");
	for(j=0;j<2;j++)
	{
		fprintf(fp,"%4.7lf\t",U[j]);
	}
	fclose(fp);
}
