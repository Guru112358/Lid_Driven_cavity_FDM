void init_variables(dmatrix &u,dmatrix&v,dmatrix &p,dmatrix &ustar,dmatrix&vstar,dmatrix &pstar);
void solve_x_mom(dmatrix&u,dmatrix&unp1,dmatrix	&v,dmatrix &p);
void solve_y_mom(dmatrix&v,dmatrix &vnp1,dmatrix &u,dmatrix &p);
void apply_bc(dmatrix&u,dmatrix &v,dmatrix &p);
double compute_residual(dmatrix&u,dmatrix &v);
void compute_collocated_values(dmatrix &u,dmatrix &v ,dmatrix &p,dmatrix &uc,dmatrix &vc ,dmatrix &pc);
void swap_variables(dmatrix &u,dmatrix &v ,dmatrix &p,dmatrix &unp1,dmatrix &vnp1,dmatrix &pnp1,dmatrix &ustar,dmatrix &vstar,dmatrix &pstar);






void init_variables(dmatrix &u,dmatrix&v,dmatrix &p,dmatrix &ustar,dmatrix&vstar,dmatrix &pstar)
{


for(int i=0;i<=(nx-1);i++)
{
	for(int j=0;j<=ny;j++)
	{
		ustar(i,j)=0;
		u(i,j)=0.0;
		
	}
}


	
for(int i=0;i<=nx;i++)
{
	for(int j=0;j<=(ny-1);j++)
	{	
		v(i,j)=0.0;
		vstar(i,j)=0;
		
	}
}


for(int i=0;i<=nx;i++)
{
	for(int j=0;j<=ny;j++)
	{	
		p(i,j)=0;
		pstar(i,j)=1;
		
	}
}
}

void solve_x_mom(dmatrix&u,dmatrix&unp1,dmatrix	&v,dmatrix &pstar)
{

double Ax,Ay,Dx,Dy;
double gradp_x;

#pragma omp parallel for schedule(static) private(Ax,Ay,Dx,Dy,gradp_x)
for(int i=1;i<=(nx-2);i++)
{
	for(int j=1;j<=(ny-1);j++)
	{

	    gradp_x=pstar(i+1,j)-pstar(i,j);

		Ax=(u(i+1,j)*u(i+1,j) - u(i-1,j)*u(i-1,j)) / (2.0*dx);
		
		Ay=0.25 *(((u(i,j) + u(i,j+1)) * (v(i,j) + v(i+1,j)) - (u(i,j) + u(i,j-1)) * (v(i+1,j-1) + v(i,j-1)) ) / dy );
		
		Dx=(u(i+1,j)+u(i-1,j)-2*u(i,j))/(dx*dx);
		
		Dy=(u(i,j+1)+u(i,j-1)-2*u(i,j))/(dy*dy);
				
		unp1(i,j)=u(i,j)-(Ax+Ay)*dt +(dt/Re)*(Dx+Dy)-(dt/dx)*gradp_x;

	
	
	}

}

}


void solve_y_mom(dmatrix&v,dmatrix &vnp1,dmatrix &u,dmatrix &pstar)
{

double Ax,Ay,Dx,Dy;
double gradp_y;

#pragma omp parallel for schedule(static) private(Ax,Ay,Dx,Dy,gradp_y)

for (int i = 1; i <= nx - 1; i++)
 {
    for (int j = 1; j <= ny - 2; j++)
	 {

		gradp_y=pstar(i,j+1)-pstar(i,j);

		Ax=0.25*(((u(i,j) + u(i,j+1)) * (v(i,j) + v(i+1,j)) - (u(i-1,j) + u(i-1,j+1)) * (v(i,j) + v(i-1,j)))/ dx);

		Ay=(v(i,j+1) * v(i,j+1) - v(i,j-1) * v(i,j-1)) / (2.0 * dy);

		Dx=(v(i+1,j)+v(i-1,j)-2*v(i,j))/(dx*dx);
		
		Dy=(v(i,j+1)+v(i,j-1)-2*v(i,j))/(dy*dy);	
		
		vnp1(i,j)=v(i,j)-(Ax+Ay)*dt +(dt/Re)*(Dx+Dy) -(dt/dy)*gradp_y;

		
    }
}

}




void solve_pressure_correction_equation(dmatrix&u,dmatrix &v,dmatrix&p_prime,int niter)
{


double a,b,c,d,rhs;

a=2.0*( ((dt)/(dx*dx))  + ((dt)/(dy*dy)) );
b=dt/(dx*dx);
c=dt/(dy*dy);
double norm;

p_prime.setZero(nx+1,ny+1);

for(int iter_count=0;iter_count<niter;iter_count++)
{
	//norm=0;
	d=0;
	rhs=0;
	
	#pragma omp parallel for schedule(static) private(d,rhs) 
	
	
	for(int i = 1; i <= (nx - 1); i++)
	{
 	  for(int j = 1; j <= (ny - 1); j++)
		{

			d=( (u(i,j) - u(i-1,j)) / dx + (v(i,j) - v(i,j-1)) / dy );  //mass balance source term;

			rhs=(b*(p_prime(i+1,j)+p_prime(i-1,j))+c*(p_prime(i,j+1)+p_prime(i,j-1))-d)/a;

			p_prime(i,j)=rhs;
			
		}

}

}
}

void correct_u(dmatrix &unp1,dmatrix &pprime,dmatrix &ustar)
{

double gradp_x;

for(int i=1;i<=(nx-2);i++)
{
	for(int j=1;j<=(ny-1);j++)
	{	
       gradp_x=pprime(i+1,j)-pprime(i,j);
	   unp1(i,j)=ustar(i,j)-(dt/dx)*gradp_x;

	}
}

}

void correct_v(dmatrix &vnp1,dmatrix &pprime,dmatrix &vstar)
{

double gradp_y;
for (int i = 1; i <= nx - 1; i++)
 {
    for (int j = 1; j <= ny - 2; j++)
	 {
	 
       gradp_y=pprime(i,j+1)-pprime(i,j);
	   vnp1(i,j)=vstar(i,j)-(dt/dy)*gradp_y;

	}
}

}

void correct_p(dmatrix &p,dmatrix &pstar,dmatrix&p_prime)
{

for(int i = 1; i <= (nx - 1); i++)
{
	for(int j = 1; j <= (ny - 1); j++)
	{
		p(i,j)=pstar(i,j)+urf_p*p_prime(i,j);
	}
}


}



void apply_bc(dmatrix&u,dmatrix &v,dmatrix &p)
{


//u velocity 

for(int j = 0; j <= ny - 1; j++)
{
	u(0,j) = 0.0;
	u(nx - 1,j) = 0.0;
}

for(int i = 0; i <= nx - 1; i++)
{
	u(i,0) = -u(i,1);
	u(i,ny) = 2 - u(i,ny - 1);
	
}
//v velocity

for(int j = 1; j <= ny - 2; j++)
{
	v(0,j) = -v(1,j);
	v(nx,j) = -v(nx - 1,j);	
}

for(int i = 0; i <= nx; i++)
{
	v(i,0) = 0.0;
	v(i,ny - 1) = 0.0;
}


//pressure
for (int i = 1; i <= (nx - 1); i++)
{
	p(i,0) = p(i,1);
	p(i,ny) = p(i,ny - 1);
}

for (int j = 0; j <= ny; j++)
{
	p(0,j) = p(1,j);
	p(nx,j) = p(nx - 1,j);
}



}


double compute_residual(dmatrix&u,dmatrix &v)
{
double res=0;

double mb;

#pragma omp parallel for schedule(static) private(mb) reduction(+:res)
for (int i = 1; i <= (nx - 1); i++)
{
	for (int j = 1; j <= (ny - 1); j++)
	{
		mb = ( (u(i,j) - u(i-1,j)) / dx + (v(i,j) - v(i,j-1)) / dy );
		res+= fabs(mb);
		
	}
}

return res;

}



void compute_collocated_values(dmatrix &u,dmatrix &v ,dmatrix &p,dmatrix &uc,dmatrix &vc ,dmatrix &pc)
{


for (int i = 0; i <= (nx - 1); i++)
{
	for (int j = 0; j <= (ny - 1); j++)
	{	
		uc(i,j) = 0.5 * (u(i,j) + u(i,j+1));
		vc(i,j) = 0.5 * (v(i,j) + v(i+1,j));
		pc(i,j) = 0.25 * (p(i,j) + p(i+1,j) + p(i,j+1) + p(i+1,j+1));
	}
}

}




void swap_variables(dmatrix &u,dmatrix &v ,dmatrix &p,dmatrix &unp1,dmatrix &vnp1,dmatrix &pstar)
{

#pragma omp parallel for schedule(static)
for (int i = 0; i <= (nx - 1); i++)
{
	for (int j = 0; j <= ny; j++)
	{
		u(i,j)=unp1(i,j);
	}
}

#pragma omp parallel for schedule(static)
for (int i = 0; i <= nx; i++)
{
	for (int j = 0; j <= (ny - 1); j++)
	{	
		v(i,j)=vnp1(i,j);
	}
}

#pragma omp parallel for schedule(static)

for (int i = 0; i <= nx; i++)
{
	for (int j = 0; j <= ny; j++)
	{
		pstar(i,j)=p(i,j);
	}
}

}




