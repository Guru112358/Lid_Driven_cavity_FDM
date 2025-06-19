//staggered variables
dmatrix u(nx,ny+1),unp1(nx,ny+1),ustar(nx,ny+1);
dmatrix v(nx+1,ny),vnp1(nx+1,ny),vstar(nx+1,ny);
dmatrix p(nx+1,ny+1),pnp1(nx+1,ny+1),p_prime(nx+1,ny+1),pstar(nx+1,ny+1);
dmatrix b(nx+1,ny+1);
//collocated final variables
dmatrix uc(nx,ny);
dmatrix vc(nx,ny);
dmatrix pc(nx,ny);

