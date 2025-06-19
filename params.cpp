#include<iostream>
#include<cmath>
#include<Eigen/Dense>
#include<fstream>
#include<omp.h>


using dmatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;   //the loops are accessed in row major so this is a dirty way to improve the performance without upending the code

//grid:
const double lx=1.0;
const double ly=1.0;
const int nx=128;
const int ny=128;


const double dx=lx/nx;
const double dy=ly/ny;

//simulation parameters:
const double dt=0.001;
const double Re=1000;
const double tol=0.0000001;
const double pressure_tol=0.000001;

const int pressure_iters=1;

//data parameter:
const int print_interval=10000;



double urf_p=0.8;
