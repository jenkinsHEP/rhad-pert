#include <iostream>

#include "Variable.h"

/******************************************************************************/

/* Single variable sampling method */

// Default constructor
Variable1d::Variable1d()
  : Variable1d(0.0, 0.0)
  {}

// Loads the mean and standard deviation at construction time
Variable1d::Variable1d(double mean_in, double error_in)
  {
    mean = mean_in;
    error = error_in;
  }

// Returns the mean or error
double Variable1d::getMean() const
  {
    return mean;
  }

double Variable1d::getError() const
  {
    return error;
  }

// Sets the mean or error
void Variable1d::setMean(double mean_in)
  {
    mean = mean_in;
  }

void Variable1d::setError(double error_in)
  {
    error = error_in;
  }

// Samples the variable
double Variable1d::sample()
  {
    return mean + error * dist(gen);
  }

/******************************************************************************/

/* Compute mean or standard deviation of a list */

double mean(double *list, int size)
  {
    double mean;
    for (int i=0; i<size; i++)
      {
        mean += list[i];
      }
    mean /= size;

    return mean;
  }

double stdev(double *list, int size)
  {
    double avg = mean(list, size);
    double stdev;
    for (int i=0; i<size; i++)
      {
        stdev += pow(avg - list[i], 2.0);
      }
    stdev = sqrt( stdev / (size - 1.0) );

    return stdev;
  }

/******************************************************************************/
