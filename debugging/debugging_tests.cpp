
#include <iostream>
#include "../custom_datatypes.h"
#include <assert.h>
#include <random>
//returnss random integer in interval with uniform distribution
double rand01() {
    // Create a random number engine based on Mersenne Twister algorithm
    std::mt19937 rng(std::random_device{}());

    // Create a uniform distribution for real numbers between 0 and 1
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Generate a random number between 0 and 1
    return dist(rng);
}

int main(){
for(int i=0;i<100;i++)
{
double f=rand01();
std::cout<<f<<",";
}
}
