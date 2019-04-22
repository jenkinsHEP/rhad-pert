#ifndef _VARIABLE_H_
#define _VARIABLE_H_

#include <iostream>
#include <random>

/******************************************************************************/

/* Single variable sampling method */

class Variable1d
  {
  public:
    // Constructors
    Variable1d();
    Variable1d(double mean, double error);

    // Returns mean or standard deviation
    double getMean() const;
    double getError() const;

    // Sets mean or standard deviation
    void setMean(double mean);
    void setError(double error);

    // Samples the distribution
    double sample();

  private:
    // Mean and standard deviation of the sample
    double mean;
    double error;

    // Unit normal distribution and random number generator
    std::normal_distribution<> dist;
    std::mt19937 gen{ std::random_device{}() };
  };

/******************************************************************************/

/* Compute the mean or standard deviation of a list */

double mean(double *list, int size);
double stdev(double *list, int size);

#endif
