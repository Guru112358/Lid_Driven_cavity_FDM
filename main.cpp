#include "params.cpp"
#include "discequations.cpp"
#include "postproc.cpp"
#include "allocate.cpp"


int main()
{


double residual=1;

int count=0;
bool loop_switch=true;

init_variables(u,v,p,ustar,vstar,pstar);
p_prime.setZero(nx+1,ny+1);

while(loop_switch)
{

 	solve_x_mom(u,ustar,v,pstar);
 	solve_y_mom(v,vstar,u,pstar);

	solve_pressure_correction_equation(ustar,vstar, p_prime,pressure_iters);

	correct_p(p,pstar,p_prime);
	correct_u(unp1,p_prime,ustar);
	correct_v(vnp1,p_prime,vstar);
	
 	apply_bc(unp1,vnp1,p);
	

	swap_variables(u,v,p,unp1,vnp1,pstar);

	residual=compute_residual(ustar,vstar);
	
	if(residual<tol&&(count!=0))	
	{
	compute_collocated_values(u,v,p,uc,vc,pc);
	write_file(uc,vc,pc,nx,ny);
	std::cout<<"converged to a tolerance of: "<<tol<<" in "	<<count<<" Iterations"<<std::endl;
	loop_switch=false;
	
	}
	
	if(count%print_interval==0)
	{
	compute_collocated_values(u,v,p,uc,vc,pc);
	write_file(uc,vc,pc,nx,ny);
	std::cout<<"|| iteration is: "<<count<<", continuity residual : " <<residual<<" ||"<<std::endl;
	
	}
	
	
	count+=1;	

}

return 0;

}
