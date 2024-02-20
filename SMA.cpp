
#include <iostream>
#include <limits>
#include <random>
#include "custom_datatypes.h"
#include "slimemould_func_prototypes.h"
#include <cassert>
#include <algorithm>
//#include "custom_math_operations.cpp"

void showSpinner(int iteration) {
    const char spinner[] = {'/', '-', '\\', '|'};
    std::cout <<"\r"<< spinner[iteration % 4] << std::flush;
}

output_SMA SMA(const model_struct *model, SMAsol_struct* sol1, kinematics_time* kin, const unsigned max_iter)
{
    std::cout<<"SMA started\n";
    //This code is a C++ translation of a modified version of the algorithm presented
    //in the following work: DOI: https://doi.org/10.1016/j.future.2020.03.055
    double weight[30][6]={}; //weight of each entry. 180 because we have 30 agents for 6 values
    double bestPositions[6]; 
    double X[30][6];
    double AllFitness[30], Violation[30]; //Fitness of each agent
    double y[30];
    double Destination_fitness, Destination_Vio;//best fitness found thus far
    const double lb = model->lb;;//lower bound
    const double ub = model->ub;;//upper bound
    const double model_tmp= ub - lb;;//space between upper and lower bound
    const double z=0.03;
    double bestVio;
    struct_fitness fitness;
    //at start of algorithm, all entries of bestPositions should be 0 
    for (int i = 0; i < 6; i++) {bestPositions[i] == 0;}

    //set fitness to infinity as a start
    const double inf = std::numeric_limits<double>::infinity();
    Destination_fitness=inf;
    Destination_Vio=inf;
    std::fill(AllFitness, AllFitness + 30, inf);
    
    for (int i = 0; i < 30; i++) {
        for (int j = 0; j < 6; j++) {
            weight[i][j] = 1.0;
            X[i][j] = rand01()*model_tmp+lb;
        }
    }

  //  Main loop
  for (int it=0; it < max_iter; it++) {
    //std::cout<<"SMA: iteration "<<it;
    if(it%21==0){showSpinner(it);} //showing a spinner in the terminal in case the user thinks that the algorithm isn't running
    double guide_points[6];
    double FitnessRange;
    double bestFitness;
    int k;
    // sort the fitness
    for (int i{0}; i < 30; i++) 
    {
      //  Check if solutions go outside the search space and bring them back
      // for first and last waypoint, we need to check for trespassing ub by calculating the vector to the end point
      for (int j{0}; j < 6; j++) {
        //check for upper bound
        if(j>=1 && j<=4){ //waypoints 1 to 4
          k = i + 30 * j;
          if(X[i][j]>(ub-model->tolerance)){ //ub - tolerance bcs path represents middle line of AGV, it also shouldn't bump against the wall on the sides
            X[i][j]=ub-model->tolerance;
          }
        }
        else if(j==0)// waypoint 0
        {
          point w;
          w[0] = model->xs + X[i][j] * b_cosd(model->SA); //vector to waypoint
          w[1] = model->ys + X[i][j] * b_sind(model->SA); //todo: find a way not to do this 30 times every iteration
          if(w[0]>(ub-model->tolerance)){X[i][j]-=(model->tolerance + w[0] - ub) / b_cosd(model->SA);} //set X back inside the search area
          if(w[1]>(ub-model->tolerance)){X[i][j]-=(model->tolerance + w[1] - ub) / b_sind(model->SA);}
        }
        else// waypoint 5
        {
          point w;
          w[0] = model->xt + X[i][j] * b_cosd(model->EA);
          w[1] = model->yt + X[i][j] * b_sind(model->EA);
          if(w[0]>(ub-model->tolerance)){X[i][j]-=(model->tolerance + w[0] - ub) / b_cosd(model->EA);} //set X back inside the search area
          if(w[1]>(ub-model->tolerance)){X[i][j]-=(model->tolerance + w[1] - ub) / b_sind(model->EA);}
        }

        if(X[i][j]<(lb+model->tolerance)){ //check for lower bound => waypoint shouldn't be closer to wall than tolerance
          X[i][j]=lb+model->tolerance; //F00 function cannot deal with 0
        }
        guide_points[j] = X[i][j]; 
      }
      //std::cout<<"iteration "<<it<<": F00 started on i="<<i<<"; ";
      bool nancond = std::isnan(guide_points[0])|| std::isnan(guide_points[1])|| std::isnan(guide_points[2])|| std::isnan(guide_points[3])|| std::isnan(guide_points[4])|| std::isnan(guide_points[5]);//check if the guide points are nan and skip F00 if yes
      if(!nancond)
      {
        fitness = F00(guide_points, model, sol1, kin, (it>max_iter/2));
        AllFitness[i]=fitness.res;
        Violation[i]=fitness.Vio;
      }
    }

    //for sorting according to Eq. 2.6, we'll call std::sort
    std::vector<int> smell_index(30);
    std::iota(smell_index.begin(), smell_index.end(), 0); //vector is filled with int values from 0 to 29
    std::sort(smell_index.begin(), smell_index.end(), [&AllFitness](int i, int j) {
    return AllFitness[i] < AllFitness[j];
    });

    for (int i = 0; i < 30; i++) {
      y[i] = AllFitness[smell_index[i]];
      }
    assert(y[10]<=y[11] && y[11]<=y[12]); //make sure values are sorted

    bestFitness = y[0];
    FitnessRange = (y[0] - y[29]) + 2.2204460492503131E-16; //  plus eps to avoid denominator zero

    // calculate the fitness weight of each slime mold
    for (int i{0}; i < 30; i++) {
      double log_arg=std::log10((y[0] - y[i]) / FitnessRange + 1);
      assert(log_arg < 1 && log_arg>-1);
      for (int j{0}; j < 6; j++) {
        if (i < 15) {
          // Eq.(2.5)
          weight[smell_index[i]][j] = 1 + rand01() * log_arg;
        } else {
          weight[smell_index[i]][j] = 1 - rand01() * log_arg;
        }
      }
    }

    //logic for determining whether new solution is better than old solution
    bool bettersol;
    bestVio=Violation[smell_index[0]];
    if(it<max_iter){bettersol=bestFitness<Destination_fitness;}
    //todo: feel free to use a different logic here. For instance, only allow improvements with lower Violation.

    // update the best fitness value and best position
    if (bettersol) { //in second half of iterations: Violation needs to decrease
      for (int j{0}; j < 6; j++) {
        int best_index[2]={smell_index[0],j};
        bestPositions[j] = X[best_index[0]][best_index[1]];
        } //best position
      Destination_fitness = bestFitness; //best fitness
      Destination_Vio=bestVio;
      //std::cout<<"iteration "<<it<<": updating best solution. fitness: "<<bestFitness<<", violation: "<<bestVio<<"\n";
    }
    //printing out best Positions for current iteration
//    std::cout<<"\nbest Positions: ";
//    for(int w=0; w<6;w++){std::cout<<bestPositions[w]<<",";}
    //std::cout<<std::endl;

    //Eq.(2.4):
    double b=1-static_cast<double>(it)/max_iter;
    double a=std::atanh(b-2.2204460492503131E-16); //-eps to avoid atanh(1)=inf
    assert(b>=0);

//let me consider an experiment:
//z*=b;
//todo: Analyze whether changing z over the course of the iterations changes the algorithm's performance.

//save X of old iteration so we wouldn't make a mess when writing onto XA and XB
    double Xold[30][6];
    for (int i = 0; i < 30; ++i) {
      for(int j=0;j<6;j++){
        Xold[i][j] = X[i][j];
      }
    }
    //  Update the Position of search agents
    for (int i{0}; i < 30; i++) {
      if (rand01() < z) {
        // Eq.(2.7)
        for (int j{0}; j < 6; j++) { X[i][j] = model_tmp * rand01() + lb;//std::cout<<"random X computed to "<<X[i+30*j]<<"\n";
        }
      } 
      else {
        double vc[6]={}, vb[6]={};
        double p = std::tanh(std::abs(AllFitness[i] - Destination_fitness)); // Eq.(2.2)
        assert(p>=0);
        for (int j{0}; j < 6; j++) {
          vb[j]=(rand01()-0.5)*2*a;
          vc[j]=(rand01()-0.5)*2*b;
          assert(abs(vb[j])<=a);
          assert(abs(vc[j])<=a);

          double r = rand01();
          int A=randi(0,29); //index of XA => value of random agent
          int B=randi(0,29);
          if(r<p)
          {
            assert(A>=0 && A<30);
            assert(B>=0 && B<30);
            assert(i+30*j<180);
            X[i][j] = bestPositions[j] + vb[j] * (weight[i][j] * Xold[A][j] - Xold[B][j]);
          }
          else
          {
            X[i][j] *=vc[j];
          }
        }
      }
    }
  }
  sol1->Violation=bestVio;
  output_SMA out;
  for(int i=0;i<6;i++)
  {out.waypoints[i]=bestPositions[i];}
  out.Vio=Destination_Vio;
  out.res=Destination_fitness;
  //std::cout<<"SMA: SMA done. best Violation: "<<bestVio<<"\n";
  return out;
}