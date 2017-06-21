// g++ density_function.cc -o density_function

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

double h;
double pi=M_PI;

double F1_3d(double phi, double r0, double R_0, double B1) {

  double integral;
  double mu, a, logs, invtan, u;
  double I0, I1, I_1, I_2, I_3, I_4, I_5;
  double a2, cosp2, r03, r0h2, r0h3, r0h_2, r0h_3;

  a = R_0/r0;
  mu = cos(phi)/a / sqrt(1. + pow(cos(phi)/a,2));

  a2 = pow(a,2);
  cosp2 = pow(cos(phi),2);
  r03 = pow(r0,3);
  r0h2 = pow(r0/h,2);
  r0h3 = pow(r0/h,3);
  r0h_2 = pow(r0/h,-2);
  r0h_3 = pow(r0/h,-3);

  I0  = phi;
  I1  = - asin( sqrt( (1.+cosp2/a2)/(1.+1./a2) ) );

  I_2 = phi +   a2 * tan(phi);
  I_4 = phi + 2*a2 * tan(phi) + 1./3.*pow(a,4) * tan(phi)*(2. + 1./cosp2);

  //u = sqrt(1.-(1.+a2)*pow(mu,2));
  u = sin(phi)*sqrt(1-mu*mu);
  logs = log(1+u) - log(1-u);
  invtan = atan(u/a);
  I1  = invtan;

  I_1 = a/2.*logs + invtan;
  I_3 = I_1 + a*(1.+a2)/4. *(2*u/(1-pow(u,2)) + logs);
  I_5 = I_3 + a*pow((1.+a2),2)/16. *( (10*u - 6*pow(u,3))/pow(1-pow(u,2),2) + 3.*logs);

  integral =  r0h3/pi  * (1./6. *I_2 - 3./40.*r0h2 *I_4 + 1./40.*r0h3 *I_5 + B1/r03 *I0);

  /*
  printf("B1: %g\n", B1);
  printf("phi: %g\n", phi);
  printf("I_2: %g\n", I_2);
  printf("I_4: %g\n", I_4);
  printf("I_5: %g\n", I_5);
  printf("I0: %g\n\n", I0);
  */

  return integral;
}


double F2_3d(double phi, double r0, double R_0, double B2) {

  double integral;
  double mu, a, logs, invtan, u;
  double I0, I1, I_1, I_2, I_3, I_4, I_5;
  double a2, cosp2, r03, r0h2, r0h3, r0h_2, r0h_3;

  a = R_0/r0;
  mu = cos(phi)/a / sqrt(1. + pow(cos(phi)/a,2));

  a2 = pow(a,2);
  cosp2 = pow(cos(phi),2);
  r03 = pow(r0,3);
  r0h2 = pow(r0/h,2);
  r0h3 = pow(r0/h,3);
  r0h_2 = pow(r0/h,-2);
  r0h_3 = pow(r0/h,-3);

  I0  = phi;
  I1  = - asin( sqrt( (1.+cosp2/a2)/(1.+1./a2) ) );

  I_2 = phi +   a2 * tan(phi);
  I_4 = phi + 2*a2 * tan(phi) + 1./3.*pow(a,4) * tan(phi)*(2. + 1./cosp2);

  //u = sqrt(1.-(1.+a2)*pow(mu,2));  The two expressions for u are equivallent. The bottom one is chosen for numerical reasons.
  u = sin(phi)*sqrt(1-mu*mu);
  logs = log(1+u) - log(1-u);
  invtan = atan(u/a);
  I1  = invtan;

  I_1 = a/2.*logs + invtan;
  I_3 = I_1 + a*(1.+a2)/4. *(2*u/(1-pow(u,2)) + logs);
  I_5 = I_3 + a*pow((1.+a2),2)/16. *( (10*u - 6*pow(u,3))/pow(1-pow(u,2),2) + 3.*logs);

  integral =  r0h3/pi  * (0.25 * (4./3. *I_2 - (r0/h) *I_3 + 0.3*r0h2 *I_4 - 
					  1./30.*r0h3 *I_5 + 1./15. *r0h_3 *I1) + B2/r03 *I0);

  /*
  printf("B2: %g\n", B2);
  printf("phi: %g\n", phi);
  printf("I_2: %g\n", I_2);
  printf("I_3: %g\n", I_3);
  printf("I_4: %g\n", I_4);
  printf("I_5: %g\n", I_5);
  printf("I0: %g\n", I0);
  printf("I1: %g\n\n", I1);
  */

  return integral;
}

double F3_3d(double phi, double r0, double R_0, double B3) {

  double integral;
  double I0, I1;
  double a, a2, cosp2, r03, r0h3, r0h_3, u, invtan, mu;

  a = R_0/r0;
  mu = cos(phi)/a / sqrt(1. + pow(cos(phi)/a,2));

  a2 = pow(a,2);
  cosp2 = pow(cos(phi),2);
  r03 = pow(r0,3);
  r0h3 = pow(r0/h,3);
  r0h_3 = pow(r0/h,-3);

  I0  = phi;
  I1  = - asin( sqrt( (1.+cosp2/a2)/(1.+1./a2) ) );

  u = sin(phi)*sqrt(1-mu*mu);
  invtan = atan(u/a);
  I1  = invtan;

  integral = r0h3/pi  * (-0.25*r0h_3 *I1 + B3/r03 *I0);

  /*
  printf("B3: %g\n", B3);
  printf("phi: %g\n", phi);
  printf("I0: %g\n", I0);
  printf("I1: %g\n\n", I1);
  */

  return integral;
}


double full_integral(double phi, double r0, double R_0) {

  double B1, B2, B3, mu, a, logs, invtan, u;
  double full_int;
  double a2, cosp2, r03, r0h2, r0h3, r0h_2, r0h_3;
  double r, R, linedist, phi1, phi2;

  if(r0==0.0) return 0.0;
  if(R_0==0.0) return 0.0;
  if(phi==0.0) return 0.0;

  a = R_0/r0;
  mu = cos(phi)/a / sqrt(1. + pow(cos(phi)/a,2));

  a2 = pow(a,2);
  cosp2 = pow(cos(phi),2);
  r03 = pow(r0,3);
  r0h2 = pow(r0/h,2);
  r0h3 = pow(r0/h,3);
  r0h_2 = pow(r0/h,-2);
  r0h_3 = pow(r0/h,-3);


  if(r0 >= 2.0*h) {
    B3 = pow(h,3) /4.;
  }
  else if(r0 > h) {
    B3 = r03/4. *(-4./3. + (r0/h) - 0.3*r0h2 + 1./30.*r0h3 - 1./15. *r0h_3+ 8./5.*r0h_2);
    B2 = r03/4. *(-4./3. + (r0/h) - 0.3*r0h2 + 1./30.*r0h3 - 1./15. *r0h_3);
  }
  else {
    B3 = r03/4. *(-2./3. + 0.3*r0h2 - 0.1*r0h3 + 7./5.*r0h_2);
    B2 = r03/4. *(-2./3. + 0.3*r0h2 - 0.1*r0h3 - 1./5.*r0h_2);
    B1 = r03/4. *(-2./3. + 0.3*r0h2 - 0.1*r0h3);
  }


  /*
  if(r0/mu < 1.0*h) {
    full_int = F1_3d(phi,r0,R_0,B1);
  }
  else if(r0/mu < 2.0*h) {
    full_int = F2_3d(phi,r0,R_0,B2);
  }
  else {
    full_int = F3_3d(phi,r0,R_0,B3);
  }

  */

  linedist = sqrt(r0*r0 + R_0*R_0);
  R = R_0/cos(phi);
  r = sqrt(r0*r0 + R*R);


    full_int = 0.0;

    if(linedist <= 1.0*h) {
      phi1 = acos(R_0/sqrt(h*h-r0*r0));
      phi2 = acos(R_0/sqrt(4.0*h*h-r0*r0));

      if(r <= 1.0*h) {
        full_int = full_int + F1_3d(phi, r0, R_0, B1) - F1_3d(0.0, r0, R_0, B1);
      }
      else if(r <= 2.0*h) {
        full_int = full_int + F1_3d(phi1, r0, R_0, B1) - F1_3d(0.0, r0, R_0, B1);
        full_int = full_int + F2_3d(phi, r0, R_0, B2) - F2_3d(phi1, r0, R_0, B2);
      }
      else { 
        full_int = full_int + F1_3d(phi1, r0, R_0, B1) - F1_3d(0.0, r0, R_0, B1);
        full_int = full_int + F2_3d(phi2, r0, R_0, B2) - F2_3d(phi1, r0, R_0, B2);
        full_int = full_int + F3_3d(phi, r0, R_0, B3) - F3_3d(phi2, r0, R_0, B3);
      }
    }


    if((linedist>1.0*h) && (linedist<=2.0*h)) {
      phi2 = acos(R_0/sqrt(4.0*h*h-r0*r0));

      if(r <= 2.0*h) {
        full_int = full_int + F2_3d(phi, r0, R_0, B2) - F2_3d(0.0, r0, R_0, B2);
      }
      else {
        full_int = full_int + F2_3d(phi2, r0, R_0, B2) - F2_3d(0.0, r0, R_0, B2);
        full_int = full_int + F3_3d(phi, r0, R_0, B3) - F3_3d(phi2, r0, R_0, B3);
      }
    }



    if(linedist > 2.0*h) {
      full_int = full_int + F3_3d(phi, r0, R_0, B3) - F3_3d(0.0, r0, R_0, B3);
    }

  return full_int;
}






double full_integral_new(double phi, double r0, double R_0) {

  double B1, B2, B3, mu, a, logs, u;
  double full_int;
  double a2, cosp, cosp2, r03, r0h2, r0h3, r0h_2, r0h_3, tanp;
  double r, R, linedist, phi1, phi2;
  double I0, I1, I_1, I_2, I_3, I_4, I_5;
  double D1, D2, D3;

  if(r0==0.0) return 0.0;
  if(R_0==0.0) return 0.0;
  if(phi==0.0) return 0.0;

  r03 = r0*r0*r0;
  r0h2 = r0/h*r0/h;
  r0h3 = r0h2*r0/h;
  r0h_2 = h/r0*h/r0;
  r0h_3 = r0h_2*h/r0;


  if(r0 >= 2.0*h) {
    B3 = h*h*h /4.;
  }
  else if(r0 > h) {
    B3 = r03/4. *(-4./3. + (r0/h) - 0.3*r0h2 + 1./30.*r0h3 - 1./15. *r0h_3+ 8./5.*r0h_2);
    B2 = r03/4. *(-4./3. + (r0/h) - 0.3*r0h2 + 1./30.*r0h3 - 1./15. *r0h_3);
  }
  else {
    B3 = r03/4. *(-2./3. + 0.3*r0h2 - 0.1*r0h3 + 7./5.*r0h_2);
    B2 = r03/4. *(-2./3. + 0.3*r0h2 - 0.1*r0h3 - 1./5.*r0h_2);
    B1 = r03/4. *(-2./3. + 0.3*r0h2 - 0.1*r0h3);
  }


  a = R_0/r0;
  a2 = a*a;

  linedist = sqrt(r0*r0 + R_0*R_0);
  R = R_0/cos(phi);
  r = sqrt(r0*r0 + R*R);


  full_int = 0.0;
  D2 = 0.0;
  D3 = 0.0;

  if(linedist <= 1.0*h) {
    ////// phi1 business /////
    phi1 = acos(R_0/sqrt(h*h-r0*r0));
    
    cosp = cos(phi1);
    cosp2 = cosp*cosp;
    mu = cosp/a / sqrt(1. + cosp2/a2);
    
    tanp = tan(phi1);
    
    I0  = phi1;
    I_2 = phi1 +   a2 * tanp;
    I_4 = phi1 + 2*a2 * tanp + 1./3.*a2*a2 * tanp*(2. + 1./cosp2);
    
    u = sin(phi1)*sqrt(1-mu*mu);
    logs = log(1+u) - log(1-u);
    I1 = atan(u/a);
    
    I_1 = a/2.*logs + I1;
    I_3 = I_1 + a*(1.+a2)/4. *(2*u/(1-u*u) + logs);
    I_5 = I_3 + a*(1.+a2)*(1.+a2)/16. *( (10*u - 6*u*u*u)/(1-u*u)/(1-u*u) + 3.*logs);
    
    D2 = -1./6.*I_2 + 0.25*(r0/h) *I_3 - 0.15*r0h2 *I_4 + 1./30.*r0h3 *I_5 - 1./60. *r0h_3 *I1 + (B1-B2)/r03 *I0;
    
    
    ////// phi2 business /////
    phi2 = acos(R_0/sqrt(4.0*h*h-r0*r0));
    
    cosp = cos(phi2);
    cosp2 = cosp*cosp;
    mu = cosp/a / sqrt(1. + cosp2/a2);
    
    tanp = tan(phi2);
    
    I0  = phi2;
    I_2 = phi2 +   a2 * tanp;
    I_4 = phi2 + 2*a2 * tanp + 1./3.*a2*a2 * tanp*(2. + 1./cosp2);
    
    u = sin(phi2)*sqrt(1-mu*mu);
    logs = log(1+u) - log(1-u);
    I1 = atan(u/a);
    
    I_1 = a/2.*logs + I1;
    I_3 = I_1 + a*(1.+a2)/4. *(2*u/(1-u*u) + logs);
    I_5 = I_3 + a*(1.+a2)*(1.+a2)/16. *( (10*u - 6*u*u*u)/(1-u*u)/(1-u*u) + 3.*logs);
    
    D3 = 1./3.*I_2 - 0.25*(r0/h) *I_3 + 3./40.*r0h2 *I_4 - 1./120.*r0h3 *I_5 + 4./15. *r0h_3 *I1 + (B2-B3)/r03 *I0 + D2;
  }
  else if(linedist <= 2.0*h) {
    ////// phi2 business /////
    phi2 = acos(R_0/sqrt(4.0*h*h-r0*r0));
    
    cosp = cos(phi2);
    cosp2 = cosp*cosp;
    mu = cosp/a / sqrt(1. + cosp2/a2);
    
    tanp = tan(phi2);
    
    I0  = phi2;
    I_2 = phi2 +   a2 * tanp;
    I_4 = phi2 + 2*a2 * tanp + 1./3.*a2*a2 * tanp*(2. + 1./cosp2);
    
    u = sin(phi2)*sqrt(1-mu*mu);
    logs = log(1+u) - log(1-u);
    I1 = atan(u/a);
    
    I_1 = a/2.*logs + I1;
    I_3 = I_1 + a*(1.+a2)/4. *(2*u/(1-u*u) + logs);
    I_5 = I_3 + a*(1.+a2)*(1.+a2)/16. *( (10*u - 6*u*u*u)/(1-u*u)/(1-u*u) + 3.*logs);
    
    D3 = 1./3.*I_2 - 0.25*(r0/h) *I_3 + 3./40.*r0h2 *I_4 - 1./120.*r0h3 *I_5 + 4./15. *r0h_3 *I1 + (B2-B3)/r03 *I0 + D2;
  }
  
  
  printf("phi1, phi2: %g %g\n", phi1, phi2);
  printf("D2, D3: %g %g\n", D2, D3);
  printf("linedist/h, r/h: %g %g\n", linedist/h, r/h);
  printf("%g %g %g %g %g %g %g\n", I0, I1, I_1, I_2, I_3, I_4, I_5);


  //////////////////////////////


  cosp = cos(phi);
  cosp2 = cosp*cosp;
  mu = cosp/a / sqrt(1. + cosp2/a2);

  tanp = tan(phi);

  I0  = phi;
  I_2 = phi +   a2 * tanp;
  I_4 = phi + 2*a2 * tanp + 1./3.*a2*a2 * tanp*(2. + 1./cosp2);

  u = sin(phi)*sqrt(1-mu*mu);
  logs = log(1+u) - log(1-u);
  I1 = atan(u/a);

  I_1 = a/2.*logs + I1;
  I_3 = I_1 + a*(1.+a2)/4. *(2*u/(1-u*u) + logs);
  I_5 = I_3 + a*(1.+a2)*(1.+a2)/16. *( (10*u - 6*u*u*u)/(1-u*u)/(1-u*u) + 3.*logs);

  
  //if(r0/mu < 1.0*h) {
  if(r < 1.0*h) {
    //full_int = F1_3d(phi,r0,R_0,B1);
    full_int = r0h3/pi  * (1./6. *I_2 - 3./40.*r0h2 *I_4 + 1./40.*r0h3 *I_5 + B1/r03 *I0);
  }
  //else if(r0/mu < 2.0*h) {
  else {
    if(r < 2.0*h) {
      //full_int = F2_3d(phi,r0,R_0,B2);
      full_int=  r0h3/pi  * (0.25 * (4./3. *I_2 - (r0/h) *I_3 + 0.3*r0h2 *I_4 - 
				     1./30.*r0h3 *I_5 + 1./15. *r0h_3 *I1) + B2/r03 *I0 + D2);
    }
    else {
      //full_int = F3_3d(phi,r0,R_0,B3);
      full_int = r0h3/pi  * (-0.25*r0h_3 *I1 + B3/r03 *I0 + D3);
    }
  }  

  //printf("mu: %g, r0/mu: %g\n",mu,r0/mu);
  //printf("I1: %g %g\n", I1, pi/2.- asin( sqrt( (1.+cosp2/a2)/(1.+1./a2) ) ));

  return full_int;
}


///////////////////////////////////
/// Numerical solution ///////////
/////////////////////////////////

double kernel(double y, double x, double z) {
  double r, W;

  r = sqrt(x*x+y*y+z*z);
  W = 1.0/h/h/h/pi;

  if(r < h) {
    W = W * (1 - 1.5*(r/h)*(r/h) + 0.75*(r/h)*(r/h)*(r/h) );
  }
  else {
    if(r < 2*h) {
      W = W * 0.25*(2.0-r/h)*(2.0-r/h)*(2.0-r/h);
    }
    else {
      W = 0.0;
    }
  }

  return W;
}

double line(double x, double z, double phi, int N) {
  int i;
  double dy, y1, y2, ym, f1, f2, fm, side;

  dy = x*tan(phi)/double(N);
  y1 = 0.0;
  f1 = kernel(y1,x,z);
  side = 0.0;

  for(i=1;i<=N;i++){
    y2 = double(i)*dy;
    ym = (y1+y2)/2.0;
    f2 = kernel(y2,x,z);
    fm = kernel(ym,x,z);
    side = side + (f1+4.0*fm+f2)/6.0*(y2-y1);
    y1 = y2;
    f1 = f2;
  }

  return side;
}

double area(double z, double r0, double R_0, double phi, int N) {
  int i;
  double dx, x1, x2, xm, f1, f2, fm, R, S;

  R = z/r0 * R_0;
  dx = R/double(N);
  x1 = 0.0;
  f1 = line(x1,z,phi,N);
  S = 0.0;

  for(i=1;i<=N;i++){
    x2 = double(i)*dx;
    xm = (x1+x2)/2.0;
    f2 = line(x2,z,phi,N);
    fm = line(xm,z,phi,N);
    S = S + (f1+4.0*fm+f2)/6.0*(x2-x1);
    x1 = x2;
    f1 = f2;
  }

  return S;
}

double volume(double r0, double R_0, double phi, int N) {
  int i;
  double dz, z1, z2, zm, f1, f2, fm, vol;

  dz = r0/double(N);
  z1 = 0.0;
  f1 = area(r0-z1,r0,R_0,phi,N);
  vol = 0.0;

  for(i=1;i<=N;i++){
    z2 = double(i)*dz;
    zm = (z1+z2)/2.0;
    f2 = area(r0-z2,r0,R_0,phi,N);
    fm = area(r0-zm,r0,R_0,phi,N);
    vol = vol + (f1+4.0*fm+f2)/6.0*(z2-z1);
    z1 = z2;
    f1 = f2;
  }

  return vol;
}



int main() {
  double phi, r0, R_0;
  double mu, a, invtan, u, I1;
  double a2, cosp2;
  int i, j, k;

  h=0.25;
  r0=0.1*h;
  R_0=1.0*h;
  phi=0.0;

  printf("(%g, %g)\n", full_integral(phi, r0, R_0), full_integral_new(phi, r0, R_0));

  for(i=0;i<10;i++) {
    for(j=0;j<10;j++) {
      for(k=0;k<10;k++) {
	r0 = (3.0/10.0)*i*h;
	R_0 = (3.0/10.0)*j*h;
	phi = (pi/2.0/10.0)*k;

	a = R_0/r0;
	mu = cos(phi)/a / sqrt(1. + pow(cos(phi)/a,2));

	a2 = pow(a,2);
	cosp2 = pow(cos(phi),2);
	
	I1  = - asin( sqrt( (1.+cosp2/a2)/(1.+1./a2) ) );
	
	u = sin(phi)*sqrt(1-mu*mu);
	invtan = atan(u/a);

	//printf("(%g %g %g) -> %g %g\n",r0,R_0,phi,I1+pi/2.,invtan);
	//printf("(%g %g %g) -> %g %g\n",r0,R_0,phi,r0*R_0*R_0*tan(phi)/6.0, volume(r0,R_0,phi,10));
	//printf("(%g %g %g) -> %g \n",r0,R_0,phi, full_integral_new(phi, r0, R_0) - volume(r0,R_0,phi,30));
      }
    }
  }

  //phi = 1.41861;
  h=1.5;
  r0=0.3*h;
  R_0=0.3*h;

  for(i=1;i<=20;i++) {
    phi = pi/2./20.*double(i-1);
    printf("(%g %g %g) -> %g %g %g\n\n",r0/h,R_0/h,phi, full_integral(phi, r0, R_0), full_integral_new(phi, r0, R_0), volume(r0,R_0,phi,20));
  }
  return 0;
}

