
#include <iostream>
#include <limits>
#include <random>
#include "custom_datatypes.h"
#include "slimemould_func_prototypes.h"
#include <cassert>
#include <algorithm>
//#include "custom_math_operations.cpp"

using xycirc = std::vector<std::vector<double>>;
void SMA(const model_struct *model, double bestPositions[6], double X[180],xycirc xCircle,xycirc yCircle)
{
    //This code is a C++ translation of a modified version of the algorithm presented
    //in the following work: DOI: https://doi.org/10.1016/j.future.2020.03.055
    SMAsol_struct sol1;
    double weight[180]={}; //weight of each entry. 180 because we have 30 agents for 6 values
    double AllFitness[30]; //Fitness of each agent
    double y[30];
    double Destination_fitness;//best fitness found thus far
    double lb = model->lb;;//lower bound
    double ub = model->ub;;//upper bound
    double model_tmp= ub - lb;;//space between upper and lower bound
    unsigned max_iter=200; //set to 200 before release
    double z=0.03;
    //at start of algorithm, all entries of bestPositions should be 0 
    for (int i = 0; i < 6; i++) {bestPositions[i] == 0;}

    //set fitness to infinity as a start
    double inf = std::numeric_limits<double>::infinity();
    Destination_fitness=inf;
    std::fill(AllFitness, AllFitness + 6, inf);
    
    for (int j{0}; j < 180; j++){
      weight[j] = 1.0;
      X[j]=rand01()*model_tmp+lb; //initialization function for quadratic map; we handle non-quadratic shape with obstacle
      } //initialize weights and X
    
//    for(int i=0;i<8;i++){X[i]=0;} //not sure what this is for
//    std::cout<<"SMA: Variables initialized\n";

  //  Main loop
  for (int it=0; it < max_iter; it++) {
    std::cout<<"\nSMA: entered main loop with it="<<it<<"\n";
    double guide_points[6];
    double S;
    double bestFitness;
    int k;
    // sort the fitness
    for (int i{0}; i < 30; i++) 
    {
      //  Check if solutions go outside the search space and bring them back
      for (int j{0}; j < 6; j++) {
        k = i + 30 * j;
        if(X[k]>ub){
        //  std::cout<<"Bringing back X["<<k<<"] from "<<X[k]<<" to "<<ub<<std::endl;
          X[k]=ub;
        }
        else if(X[k]<lb){
        //  std::cout<<"Bringing back X["<<k<<"] from "<<X[k]<<" to "<<lb<<std::endl;
          X[k]=lb;
        }
        guide_points[j] = X[k]; 
      }
      //std::cout<<"guide points: ";
      //for(int w=0; w<6;w++){std::cout<<guide_points[w]<<",";} //in the Matlab Code, X is 30x6, which makes sense: 30 candidates,6 points
      //std::cout<<std::endl;
      //std::cout<<"\nF00 started on i="<<i<<"; ";
      bool nancond = std::isnan(guide_points[0])|| std::isnan(guide_points[1])|| std::isnan(guide_points[2])|| std::isnan(guide_points[3])|| std::isnan(guide_points[4])|| std::isnan(guide_points[5]);//check if the guide points are nan and skip F00 if yes
      if(!nancond)
      {
        //AllFitness[i] = F00(guide_points, model, &sol1, &xCircle,&yCircle);
        //use another fitness function to make sure that SMA works
        AllFitness[i] = (guide_points[0]-46)*(guide_points[0]-46)+(guide_points[1]-58)*(guide_points[1]-58)+(guide_points[2]-48)*(guide_points[2]-48)+(guide_points[3]-58)*(guide_points[3]-58)+(guide_points[4]-21)*(guide_points[4]-21)+(guide_points[5]-5)*(guide_points[5]-5);
        //assert(AllFitness[i]>=0 && !std::isnan(AllFitness[i]));
      }
     //std::cout<<"F00 finished on i="<<i<<"\n";
    }

    //sort(AllFitness, smell_index); // Eq.(2.6)
    //for sorting according to Eq. 2.6, we'll call std::sort
    std::vector<int> smell_index(30);
    std::iota(smell_index.begin(), smell_index.end(), 0); //vector is filled with int values from 0 to 29
    std::sort(smell_index.begin(), smell_index.end(), [&AllFitness](int i, int j) {
    return AllFitness[i] < AllFitness[j];
    });
    //assert(AllFitness[smell_index[0]]<=AllFitness[smell_index[1]] && AllFitness[smell_index[0]]<=AllFitness[smell_index[30]]);

    for (int i = 0; i < 30; i++) {
      y[i] = AllFitness[smell_index[i]];
      //std::cout<<"y["<<i<<"]: "<<y[i]<<"\n";
      }
    assert(y[10]<=y[11] && y[11]<=y[12]); //make sure values are sorted

    bestFitness = y[0];
    S = (y[0] - y[29]) + 2.2204460492503131E-16; //  plus eps to avoid denominator zero
    //std::cout<<"SMA: Values sorted by Fitness\n";

    // calculate the fitness weight of each slime mold
    //I didn't test this part, but it looks correct to me
    for (int i{0}; i < 30; i++) {
      double log_arg=std::log10((y[0] - y[i]) / S + 1);
      //std::cout<<"log_arg: "<<std::log10((y[0] - y[i]) / S + 1)<<"\n";
      assert(log_arg < 1 && log_arg>-1);
      for (int j{0}; j < 6; j++) {
        if (i < 15) {
          // Eq.(2.5)
          weight[(smell_index[i] + 30 * j) - 1] = 1 + rand01() * log_arg;
          //std::cout<<"weight "<<smell_index[i] + 30 * j<<": "<<weight[(smell_index[i] + 30 * j) - 1]<<std::endl;
        } else {
          weight[(smell_index[i] + 30 * j) - 1] = 1 - rand01() * log_arg;
          //std::cout<<"weight "<<smell_index[i] + 30 * j<<": "<<weight[(smell_index[i] + 30 * j) - 1]<<std::endl;
        }
      }
    }

for(int k=0;k<180;k++){assert(weight[k]>=0);} //seems to ve true for first iteration
    //std::cout<<"SMA: fitness weights calculated\n";
    //std::cout<<"bestFitness:"<<bestFitness<<" Destination_fitness"<<Destination_fitness<<"\n";
    // update the best fitness value and best position
    if (bestFitness < Destination_fitness) {
      for (int j{0}; j < 6; j++) {
        int best_index=(smell_index[0] + 30 * j) - 1;
        //std::cout<<"X["<<best_index<<"]: "<<X[best_index]<<"\n";
        bestPositions[j] = X[best_index];
        } //best position
      //std::cout<<"updating destination fitness. bestFitness: "<<bestFitness<<", Destination Fitness: "<<Destination_fitness<<std::endl;
      Destination_fitness = bestFitness; //best fitness
    }
    //printing out best Positions for current iteration
    //std::cout<<"best Positions:\n";
    //for(int w=0; w<6;w++){std::cout<<bestPositions[w]<<",";}
    std::cout<<std::endl;
    //Eq.(2.4):
    double b=1-static_cast<double>(it)/max_iter;
    double a=std::atanh(b-2.2204460492503131E-16); //-eps to avoid atanh(1)=inf
    //std::cout<<"it: "<<it<<", max_iter: "<<max_iter<<", b: "<<b<<", a: "<<a<<"\n";
    assert(b>=0);

//    std::cout<<"SMA: updated best fitness value and best position\n";
//save X of old iteration so we wouldn't make a mess when writing onto XA and XB
    double Xold[180];
    for (int i = 0; i < 180; ++i) {
        Xold[i] = X[i];
    }
    //  Update the Position of search agents
    for (int i{0}; i < 30; i++) {
      if (rand01() < z) {
        // Eq.(2.7)
        //std::cout<<"agent "<<i<<": random<z\n";
        for (int j{0}; j < 6; j++) { X[i + 30 * j] = model_tmp * rand01() + lb;//std::cout<<"random X computed to "<<X[i+30*j]<<"\n";
        }
      } 
      else {
        double vc[6]={}, vb[6]={};
        double p = std::tanh(std::abs(AllFitness[i] - Destination_fitness)); // Eq.(2.2)
        assert(p>=0);
        for (int j{0}; j < 6; j++) {
          vb[j]=rand01()*a;// used to be (rand01()-0.5)*2*a;, but neg values don't make sense in this application
          vc[j]=rand01()*b;// used to be vc[j]=(rand01()-0.5)*2*b;, but negative values don't make sense in this application
          assert(abs(vb[j])<=a);
          assert(abs(vc[j])<=a);

          double r = rand01();
          int A=randi(1,30); //index of XA => value of random agent
          int B=randi(1,30);
          if(r<p)
          {
            int index_A=A+30*j - 1; //this one might be wrong
            int index_B=B+30*j - 1;
            assert(index_A>=0 && index_A<180);
            assert(index_B>=0 && index_B<180);
            assert(i+30*j<180);
            X[i + 30 * j] = bestPositions[j] + vb[j] * (weight[i + 30 * j] * Xold[index_A] - Xold[index_B]);
            //std::cout<<"r<p: X["<<i+30*j<<"]: "<<X[i+30*j]<<std::endl;
            //std::cout<<"vb: "<<vb[j]<<", weight: "<<weight[i+30*j]<<", XA: "<<Xold[index_A]<<", XB: "<<Xold[index_B]<<std::endl;
            //if(std::isnan(X[i+30*j])){std::cout<<"X=nan detected for k="<<i+30*j<<".\n";
            //std::cout<<"bestPositions:"<<bestPositions[0]<<","<<bestPositions[1]<<","<<bestPositions[2]<<","<<bestPositions[3]<<","<<bestPositions[4]<<","<<bestPositions[5]<<"\n";
            //std::cout<<"vb["<<j<<"]: "<<vb[j]<<"\n";
            //std::cout<<"weight[k]: "<<weight[i+30*j]<<"\n";
            //std::cout<<"XA: "<<X[index_A]<<",XB: "<<X[index_B]<<"\n";
            //}
          }
          else
          {
            X[i+30*j] *=vc[j];
            //std::cout<<"r>p: X["<<i+30*j<<"]: "<<X[i+30*j]<<" with vc[j]="<<vc[j]<<std::endl;
          }
        //          if ~isvalid(X(i,:)) %prüfen ob Spline zwischen den Punkte in
        //          den Hindernissen liegt
        //              X(i,:) = unifrnd(lb,ub,dim); %Wenn Punkte nicht korrekt
        //              liegen werden sie zufällig auf dem Gebiet gleichverteilt
        //          end
        }
      }
    }
  //std::cout<<"Updated the position of search agents. Iteration "<<it<<" done\n";
  //std::cout<<"printing X:\n";
  //for(k=0;k<180;k++){std::cout<<k<<","<<X[k]<<"\n";//assert(X[k]>=0 && !std::isnan(X[k]));
  //}
  }
std::cout<<"SMA: SMA done\n";
}