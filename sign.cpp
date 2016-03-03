#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <tuple>
static inline double sign(double a,double b) { return a= fabs(a),(b<0)?-a:a; }
int main() {
double a =2, b= -3;
std::cout<<sign(a,b)<<'\t'<<sign(b,a)<<std::endl;

}
