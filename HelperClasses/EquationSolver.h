#ifndef _Equation_Solver_h_
#define _Equation_Solver_h_


#include<cmath>
#include<string>
#include<iostream>
#include<vector>
#include<complex>


/**
 *\Function EquationSolver:
 *
 *Solves 3rd degree equations
 *
 * \Author A. Orso M. Iorio
 * 
 *
 *\version  $Id: EquationSolver.h,v 1.1 2013/02/27 12:18:42 degrutto Exp $
 *
 *
*/



template <class T>
std::vector< T > const EquationSolve(const T & a, const T & b,const T & c,const T & d){


  std::vector<T> result;


  std::complex<T> x1;
  std::complex<T> x2;
  std::complex<T> x3;

  if (a != 0) {
    
    T q = (3*a*c-b*b)/(9*a*a);
    T r = (9*a*b*c - 27*a*a*d - 2*b*b*b)/(54*a*a*a);
    T Delta = q*q*q + r*r;

    std::complex<T> s;
    std::complex<T> t;

    T rho=0;
    T theta=0;
    
    if( Delta<=0){
      rho = sqrt(-(q*q*q));

      theta = acos(r/rho);

      s = std::polar<T>(sqrt(-q),theta/3.0); 
      t = std::polar<T>(sqrt(-q),-theta/3.0); 
    }
    
    if(Delta>0){ 
      s = std::complex<T>(cbrt(r+sqrt(Delta)),0);
      t = std::complex<T>(cbrt(r-sqrt(Delta)),0);
    }
  
    std::complex<T> i(0,1.0); 
    
    
     x1 = s+t+std::complex<T>(-b/(3.0*a),0);
     x2 = (s+t)*std::complex<T>(-0.5,0)-std::complex<T>(b/(3.0*a),0)+(s-t)*i*std::complex<T>(sqrt(3)/2.0,0);
     x3 = (s+t)*std::complex<T>(-0.5,0)-std::complex<T>(b/(3.0*a),0)-(s-t)*i*std::complex<T>(sqrt(3)/2.0,0);

    if(fabs(x1.imag())<0.0001)result.push_back(x1.real());
    if(fabs(x2.imag())<0.0001)result.push_back(x2.real());
    if(fabs(x3.imag())<0.0001)result.push_back(x3.real());

    return result;
  }
  else{return result;}


  return result;
}

#endif 
