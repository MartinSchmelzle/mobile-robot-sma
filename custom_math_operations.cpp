#include <cmath>
#include <iostream>
#include <random>
#include <cstring>
#include<algorithm>

const double PI = 3.14159265358979323846;

double b_cosd(double x)
{
    x=x*PI/180;
    x=std::cos(x);
    return x;
}

double b_sind(double x)
{
    x=x*PI/180;
    x=std::sin(x);
    return x;
}

// Function to generate random numbers between 0 and 1
double rand01() {
    // Create a random number engine based on Mersenne Twister algorithm
    std::mt19937 rng(std::random_device{}());

    // Create a uniform distribution for real numbers between 0 and 1
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Generate a random number between 0 and 1
    return dist(rng);
}

//Random numbers generator for 6 numbers

void c_rand(double r[6])
{
  for (int k{0}; k < 6; k++) {
    r[k] =rand01(); 
  }
}
//Sorting algorithm copied from MATLAB Coder output
//depends on function merge

void merge(int idx[30], double x[30], int offset, int np, int nq, int iwork[30], double xwork[30])
{
  if (nq != 0) {
    int iout;
    int n_tmp;
    int p;
    int q;
    n_tmp = np + nq;
    for (int j{0}; j < n_tmp; j++) {
      iout = offset + j;
      iwork[j] = idx[iout];
      xwork[j] = x[iout];
    }
    p = 0;
    q = np;
    iout = offset - 1;
    int exitg1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork[p] <= xwork[q]) {
        idx[iout] = iwork[p];
        x[iout] = xwork[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx[iout] = iwork[q];
        x[iout] = xwork[q];
        if (q + 1 < n_tmp) {
          q++;
        } else {
          q = iout - p;
          for (int j{p + 1}; j <= np; j++) {
            iout = q + j;
            idx[iout] = iwork[j - 1];
            x[iout] = xwork[j - 1];
          }
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}
//sorting algorithm
void sort(double x[30], int idx[30])
{
  double xwork[30];
  double x4[4];
  int iwork[30];
  int b_i1;
  int i;
  int i1;
  int i3;
  int ib;
  int nNaNs;
  int quartetOffset;
  signed char idx4[4];
  x4[0] = 0.0;
  idx4[0] = 0;
  x4[1] = 0.0;
  idx4[1] = 0;
  x4[2] = 0.0;
  idx4[2] = 0;
  x4[3] = 0.0;
  idx4[3] = 0;
  std::memset(&idx[0], 0, 30U * sizeof(int));
  std::memset(&xwork[0], 0, 30U * sizeof(double));
  nNaNs = 0;
  ib = 0;
  for (int k{0}; k < 30; k++) {
    if (std::isnan(x[k])) {
      idx[29 - nNaNs] = k + 1;
      xwork[29 - nNaNs] = x[k];
      nNaNs++;
    } else {
      ib++;
      idx4[ib - 1] = static_cast<signed char>(k + 1);
      x4[ib - 1] = x[k];
      if (ib == 4) {
        double d;
        double d1;
        int i4;
        quartetOffset = k - nNaNs;
        if (x4[0] <= x4[1]) {
          i1 = 1;
          ib = 2;
        } else {
          i1 = 2;
          ib = 1;
        }
        if (x4[2] <= x4[3]) {
          i3 = 3;
          i4 = 4;
        } else {
          i3 = 4;
          i4 = 3;
        }
        d = x4[i3 - 1];
        d1 = x4[i1 - 1];
        if (d1 <= d) {
          d1 = x4[ib - 1];
          if (d1 <= d) {
            i = i1;
            b_i1 = ib;
            i1 = i3;
            ib = i4;
          } else if (d1 <= x4[i4 - 1]) {
            i = i1;
            b_i1 = i3;
            i1 = ib;
            ib = i4;
          } else {
            i = i1;
            b_i1 = i3;
            i1 = i4;
          }
        } else {
          d = x4[i4 - 1];
          if (d1 <= d) {
            if (x4[ib - 1] <= d) {
              i = i3;
              b_i1 = i1;
              i1 = ib;
              ib = i4;
            } else {
              i = i3;
              b_i1 = i1;
              i1 = i4;
            }
          } else {
            i = i3;
            b_i1 = i4;
          }
        }
        idx[quartetOffset - 3] = idx4[i - 1];
        idx[quartetOffset - 2] = idx4[b_i1 - 1];
        idx[quartetOffset - 1] = idx4[i1 - 1];
        idx[quartetOffset] = idx4[ib - 1];
        x[quartetOffset - 3] = x4[i - 1];
        x[quartetOffset - 2] = x4[b_i1 - 1];
        x[quartetOffset - 1] = x4[i1 - 1];
        x[quartetOffset] = x4[ib - 1];
        ib = 0;
      }
    }
  }
  if (ib > 0) {
    signed char perm[4];
    perm[1] = 0;
    perm[2] = 0;
    perm[3] = 0;
    if (ib == 1) {
      perm[0] = 1;
    } else if (ib == 2) {
      if (x4[0] <= x4[1]) {
        perm[0] = 1;
        perm[1] = 2;
      } else {
        perm[0] = 2;
        perm[1] = 1;
      }
    } else if (x4[0] <= x4[1]) {
      if (x4[1] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 2;
        perm[2] = 3;
      } else if (x4[0] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 3;
        perm[2] = 2;
      } else {
        perm[0] = 3;
        perm[1] = 1;
        perm[2] = 2;
      }
    } else if (x4[0] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 3;
    } else if (x4[1] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 3;
      perm[2] = 1;
    } else {
      perm[0] = 3;
      perm[1] = 2;
      perm[2] = 1;
    }
    i = static_cast<unsigned char>(ib);
    for (int k{0}; k < i; k++) {
      quartetOffset = ((k - nNaNs) - ib) + 30;
      b_i1 = perm[k];
      idx[quartetOffset] = idx4[b_i1 - 1];
      x[quartetOffset] = x4[b_i1 - 1];
    }
  }
  ib = (nNaNs >> 1) + 30;
  for (int k{0}; k <= ib - 31; k++) {
    i1 = (k - nNaNs) + 30;
    quartetOffset = idx[i1];
    idx[i1] = idx[29 - k];
    idx[29 - k] = quartetOffset;
    x[i1] = xwork[29 - k];
    x[29 - k] = xwork[i1];
  }
  if ((nNaNs & 1) != 0) {
    i = ib - nNaNs;
    x[i] = xwork[i];
  }
  if (30 - nNaNs > 1) {
    std::memset(&iwork[0], 0, 30U * sizeof(int));
    i3 = (30 - nNaNs) >> 2;
    quartetOffset = 4;
    while (i3 > 1) {
      if ((i3 & 1) != 0) {
        i3--;
        ib = quartetOffset * i3;
        i1 = 30 - (nNaNs + ib);
        if (i1 > quartetOffset) {
          merge(idx, x, ib, quartetOffset, i1 - quartetOffset, iwork, xwork);
        }
      }
      ib = quartetOffset << 1;
      i3 >>= 1;
      for (int k{0}; k < i3; k++) {
        merge(idx, x, k * ib, quartetOffset, quartetOffset, iwork, xwork);
      }
      quartetOffset = ib;
    }
    if (30 - nNaNs > quartetOffset) {
      merge(idx, x, 0, quartetOffset, 30 - (nNaNs + quartetOffset), iwork,
            xwork);
    }
  }
}

//linspace from Matlab
std::vector<double> linspace(double start, double end, int numPoints) {
    std::vector<double> result;
    if (numPoints <= 1) {
        result.push_back(start);
        return result;
    }

    double step = (end - start) / (numPoints - 1);
    for (int i = 0; i < numPoints; ++i) {
        double value = start + i * step;
        result.push_back(value);
    }
    return result;
}

//diff from Matlab
std::vector<double> diff(const std::vector<double>& input) {
    std::vector<double> differences;
    if (input.size() < 2) {
        // If input size is less than 2, cannot compute differences
        std::cerr << "Input vector should have at least two elements." << std::endl;
        return differences;
    }

    for (size_t i = 1; i < input.size(); ++i) {
        differences.push_back(input[i] - input[i - 1]);
    }

    return differences;
}

//equivalent to sign function in Matlab
int sign(double n)
{
  if(n>1e-5){return 1;}
  else if(n<=1e-5 && n>=-1e-5){return 0;}
  else {return -1;}
}

//returnss random integer in interval with uniform distribution
int randi(int lower_bound, int upper_bound) {
    // Seed for the random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Create a uniform distribution for integers
    std::uniform_int_distribution<int> distribution(lower_bound, upper_bound);

    // Generate a random integer
    return distribution(gen);
}