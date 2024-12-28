#include <iostream>
#include <iomanip>
#include<cmath>
#include <fstream>
#include<cstdlib>

using namespace std;
const double ti = 0.0;
const double tf = 1.0;
const double accuracy = 0.00000001;

double analytic_fn(double t){
    double x = -1 + exp(pow(t,3)/3 - 2*pow(t,2) + 4*t); 
    return x;
}


// all doubles below are the functions that I will need in the for loop that i will iterate using t0 and y0
double f(double t,double x){ // this is the function for the first derivative or slope of x (dx/dt)
    double a = pow((t-2),2)*(x+1);
    return a;
}

// all ki terms are the different orders of the runge-kutta method. doing these seperately so that the 
// code is much easier to read and implement
// all these ki are for the 4th order runge-kutta method
double k1(double t, double x) {
    double k1 = f(t,x);
    return k1;
}

double k2(double t,double x,double h){
    double k2 = f(t + (h/2), x + h*(k1(t,x)/2));
    return k2;
}

double k3(double t,double x,double h){
    double k3 = f(t + (h/2), x + h*(k2(t,x,h)/2));
    return k3;
}

double k4(double t,double x,double h){
    double k4 = f(t + h, x + h*k3(t,x,h));
    return k4;
}

double x4_next(double t,double x, double h){
    double x4_next = x + (h/6)*(k1(t,x) + 2*k2(t,x,h) + 2*k3(t,x,h) + k4(t,x,h));
    return x4_next;
}


int main(){
    const int len = 10;
    double discrepency[len] = {};
    double H[len] = {1.0,0.5,0.25, 0.1,0.05,0.025, 0.01,0.005,0.0025, 0.001};
    // for the error of the numerical value at t=1
    for(int j=0; j<len; j++){
        double h = H[j];
        double x_true = 0;
        double x4_0 = 0.0; // initial y term for the 4th order rk method
        double t0 = ti; // time does not need to be seperate for the different methods
        int steps = (tf-ti)/h;
        for(int i=0; i<steps; i++){
            // t_array[i] = t0;
            // x4_array[i] = x4_0;
            // xanalytic_array[i] = analytic_fn(t0);
            // discrepency[i] = xanalytic_array[i] - x4_array[i];
            x4_0 = x4_next(t0,x4_0,h);
            t0 = t0+h;
            x_true = analytic_fn(t0);
        }
    discrepency[j] = x_true - x4_0;
    cout << setprecision(9) << "h = " << h << ", t = " << t0 << ", x4(t=1) = " << x4_0 << ", x_true(t=1) = " << x_true << ", steps = " << steps << "\n";
    // RK4 is accurate to 8 significant figures for values of h below but not equal to 0.005

    }    
    

    // for the plot of rk4 with h = 0.1
    const int N = 1+(tf-ti)/0.1;
    double x4_array[N] = {};
    double t_array[N] = {};
    double xanalytic_array[N] = {};
    const double h1 = H[3];
    double x4_01 = 0.0;
    double t01 = 0.0;
    for(int i=0;i<N;i++){
        t_array[i] = t01;
        x4_array[i] = x4_01;
        xanalytic_array[i] = analytic_fn(t01);
        x4_01 = x4_next(t01,x4_01,h1);
        t01 = t01 + h1;

    }

    ofstream dataFile("PNS_RK_method.txt");
    ofstream dataFile2("please.txt");
    ofstream functionFile("PNS_RK_func_file.txt");
    for(int i=0;i<len;i++){
        dataFile << H[i] << " " << discrepency[i] << " " << accuracy << "\n";
    }
     dataFile.close();

    for(int i=0;i<N;i++){
        dataFile2 << t_array[i] << " " << x4_array[i] << " " << xanalytic_array[i] << " " << xanalytic_array[i] - x4_array[i] << "\n";
    }
    dataFile2.close();


    // plotting error of RK4 at t=1 for multiple step sizes
    system("gnuplot -p -e \"set xlabel 'h size'; "
    "set ylabel 'Error at t=1'; "
    "set logscale x;"
    "set logscale y;"
    "set title 'Error between RK4 and Analytic Solution at t=1 for Different h Size'; "
    "plot 'PNS_RK_method.txt' using 1:2 with linespoints title 'Error', ""'PNS_RK_method.txt' using 1:3 with lines title '8 s.f.'\""); 

    // plotting RK4 for step size h = 0.1
    system("gnuplot -p -e \"set xlabel 't'; "
    "set ylabel 'x'; "
    "set title 'Plot of RK4 for Step Size h = 0.1'; "
    "plot 'please.txt' using 1:3 with linespoints title 'Analytic Solution',""'please.txt' using 1:2 with linespoints title 'RK4'\""); 

    // plotting the error of the RK4 for step size h=0.1
    system("gnuplot -p -e \"set xlabel 't'; "
    "set ylabel 'Error'; "
    "set title 'Error of RK4 for Step Size h = 0.1'; "
    "plot 'please.txt' using 1:4 with linespoints title 'Error'\""); 

    return 0;
}

// something to note it appear that the second order rk method struggles when it is approximating 
// an area of extreme slope and the 4th order does much better there. this is assuming that the 4th 
// order rk is more accurate though. have no proof of this considering there is no analytical solution for 
// the approximated equation