#ifndef STATS_H
#define STATS_H

#include <numeric>
#include <cmath>
#include <algorithm>
#include <functional>

/*
#include <valarray>
#include <iostream>
*/

using namespace std;

struct OutOfRange {
  OutOfRange(int min, int max) : min_(min), max_(max) {}
  bool operator()(int x) {
    return (x < min_) || (x > max_);
  }
  int min_;
  int max_;
};

// Median
//
//  - It would be preferable to just send the iterators, but how
//    to retrieve the median from only iterators?
//
template<class T>
T findMedian(valarray<T> x, bool sorted = false) {
  valarray<T> newx(x);
  if (! sorted) sort(&newx[0],&newx[0] + newx.size());
  if (newx.size() % 2)
    return newx[(newx.size()-1)/2];
  else
    return 0.5*(newx[newx.size()/2-1] + newx[newx.size()/2]);
}

// Clipped statistics
//
template<int N, int I, class T>
void clippedStats(valarray<T> x, valarray<T> m,
		  T& median, T& mean, T& std_dev) {

  double tmedian, tmean, tstd_dev;

  // Mask out values and calculate initial statistics
  //
  valarray<T> xcopy(x.size());
  double sum = 0.0;
  double sum2 = 0.0;
  unsigned int nelem = 0;
  for (unsigned int i=0; i<x.size(); i++) {
    if (m[i] == 0) {
      xcopy[nelem++] = x[i];
      sum += x[i];
      sum2 += x[i] * x[i];
    }
  }
  tmean = sum/nelem;
  tstd_dev = sqrt( sum2/nelem - tmean*tmean );
  //  tstd_dev = computeStdDev(&xcopy[0], &xcopy[0]+xcopy.size(), tmean);

  // There has to be a better way to just shrink xcopy without going
  // through xtmp.
  //
  valarray<T> xtmp(nelem);
  for (unsigned int i=0; i<nelem; i++)
    xtmp[i] = xcopy[i];
  xcopy.resize(nelem);
  xcopy = xtmp;

  // Sort and find median
  sort(&xcopy[0],&xcopy[0]+xcopy.size());
  tmedian = findMedian(xcopy,true);
  
  
  // Do iterative clipping
  //
  for (unsigned int i=0; i<I; i++) {
    
    // Determine range of sorted xcopy to keep
    //
    int i1 = 0;
    int i2 = xcopy.size()-1;
    while (i1 < int(xcopy.size()) && xcopy[i1] <= tmedian - N*tstd_dev)
      i1++;
    while (i2 >= 0 && xcopy[i2] >= tmedian + N*tstd_dev)
      i2--;
  
    if (i1 < i2) {
      // Create a temporary array (I'm sure this can be solved in
      // a better way).
      //
      double sum = 0.0;
      double sum2 = 0.0;
      unsigned int nelem = i2-i1;
      valarray<T> xtmp(nelem);
      for (unsigned int n=0; n<nelem; n++) {
	xtmp[n] = xcopy[i1+n];
	sum += xtmp[n];
	sum2 += xtmp[n] * xtmp[n];
      }
      xcopy.resize(nelem);
      xcopy = xtmp;
      
      tmedian = findMedian(xcopy,true);
      tmean = sum/nelem;
      tstd_dev = sqrt( sum2/nelem - tmean*tmean );
      //      std_dev = computeStdDev(&xcopy[0], &xcopy[0]+xcopy.size(), mean);
    } else {
      tmedian  = 0.0;
      tmean    = 0.0;
      tstd_dev = 0.0;
      break;
    }
  }

  median = T(tmedian);
  mean = T(tmean);
  std_dev = T(tstd_dev);
}

template<int N, int I, class T>
void clippedStats(valarray<T> x, T& median, T& mean, T& std_dev) {
  valarray<T> m(x.size(),0);
  clippedStats<N,I>(x, m, median, mean, std_dev);
}

template<class T>
void runningStats(valarray<T> x, valarray<T> m, int N,
		  T& median, T& mean, T& std_dev, T& nuse) {
  median = mean = std_dev = NULL;
  nuse = NULL;
  
  double tmedian = NULL, tmean = NULL, tstd_dev = NULL;

  // Set the masked values to the min value
  //
  valarray<T> xcopy(x);
  int nmask = 0;
  T xmin = xcopy.min();
  for (unsigned int i=0; i<xcopy.size(); i++)
    if (m[i]) {
      xcopy[i] = xmin;
      nmask++;      // Determine the number of masked elements
    }

  sort(&xcopy[0],&xcopy[0]+xcopy.size());

  int len = xcopy.size()-nmask-2*N;
  if (len > 0) {
    valarray<T> xcrop(len);
    double sum = 0.0;
    double sum2 = 0.0;
    for (int i=0; i<len; i++) {
      xcrop[i] = xcopy[nmask+N+i];
      sum += xcrop[i];
      sum2 += xcrop[i] * xcrop[i];
    }
    nuse = len;
    tmedian = findMedian(xcrop,true);
    tmean = sum/len;
    tstd_dev = sqrt( sum2/len - tmean*tmean );
    //    std_dev = computeStdDev(&xcrop[0],&xcrop[0]+xcrop.size(),tmean);
  } else {
    nuse = 0;
  }

  median = T(tmedian);
  mean = T(tmean);
  std_dev = T(tstd_dev);
}


//  I am not sure why I wrote this putting N in the template definition.
//  The new version ov runningStats defined above does not do this.
//
template<int N, class T>
void runningStats(valarray<T> x, valarray<T> m,
		  T& median, T& mean, T& std_dev, T& nuse) {
  median = mean = std_dev = NULL;
  nuse = NULL;
  
  double tmedian = NULL, tmean = NULL, tstd_dev = NULL;

  // Set the masked values to the min value
  //
  valarray<T> xcopy(x);
  int nmask = 0;
  T xmin = xcopy.min();
  for (unsigned int i=0; i<xcopy.size(); i++)
    if (m[i]) {
      xcopy[i] = xmin;
      nmask++;      // Determine the number of masked elements
    }

  sort(&xcopy[0],&xcopy[0]+xcopy.size());

  int len = xcopy.size()-nmask-2*N;
  if (len > 0) {
    valarray<T> xcrop(len);
    double sum = 0.0;
    double sum2 = 0.0;
    for (int i=0; i<len; i++) {
      xcrop[i] = xcopy[nmask+N+i];
      sum += xcrop[i];
      sum2 += xcrop[i] * xcrop[i];
    }
    nuse = len;
    tmedian = findMedian(xcrop,true);
    tmean = sum/len;
    tstd_dev = sqrt( sum2/len - tmean*tmean );
    //    std_dev = computeStdDev(&xcrop[0],&xcrop[0]+xcrop.size(),tmean);
  } else {
    nuse = 0;
  }

  median = T(tmedian);
  mean = T(tmean);
  std_dev = T(tstd_dev);
}

template<int N, class T>
void runningStats(valarray<T> x, T& median, T& mean, T& std_dev, T& nuse) {
  valarray<T> m(x.size(),0.0);
  runningStats<N>(x,m,median,mean,std_dev,nuse);
}


template<int N, class T>
T nthPower(T x) {
  T ret = x;
  for (int i=1; i < N; i++) {
    ret *= x;
  }
  return ret;
}

template<class T, int N>
struct SumDiffNthPower {
  SumDiffNthPower(T x) : mean_(x) {};
  T operator()(T sum, T current) {
    return sum + nthPower<N>(current - mean_);
  }
  T mean_;
};

template<class T, int N, class Iter_T>
T nthMoment(Iter_T first, Iter_T last, T mean) {
  size_t cnt = distance(first, last);
  return accumulate(first, last, T(), SumDiffNthPower<T, N>(mean)) / cnt;
}

template<class T, class Iter_T>
T computeVariance(Iter_T first, Iter_T last, T mean) {
  return nthMoment<T, 2>(first, last, mean);
}

template<class T, class Iter_T>
T computeStdDev(Iter_T first, Iter_T last, T mean) {
  return sqrt(computeVariance(first, last, mean));
}

template<class T, class Iter_T>
T computeSkew(Iter_T begin, Iter_T end, T mean) {
  T m3 = nthMoment<T, 3>(begin, end, mean);
  T m2 = nthMoment<T, 2>(begin, end, mean);
  return m3 / (m2 * sqrt(m2));
}

template<class T, class Iter_T>
T computeKurtosisExcess(Iter_T begin, Iter_T end, T mean) {
  T m4 = nthMoment<T, 4>(begin, end, mean);
  T m2 = nthMoment<T, 2>(begin, end, mean);
  return m4 / (m2 * m2) - 3;
}

template<class T, class Iter_T>
void computeStats(Iter_T first, Iter_T last, T& sum, T& mean,
		  T& var, T& std_dev, T& skew, T& kurt) {
  size_t cnt = distance(first, last);
  sum = accumulate(first, last, T());
  mean = sum / cnt;
  var = computeVariance(first, last, mean);
  std_dev = sqrt(var);
  skew = computeSkew(first, last, mean);
  kurt = computeKurtosisExcess(first, last, mean);
}

/*
int main() {
  valarray<float> v(6);
  v[0] = 2;
  v[1] = 4;
  v[2] = 8;
  v[3] = 10;
  v[4] = 99;
  v[5] = 1;

  valarray<float> m(6);
  m[0] = 0;
  m[1] = 0;
  m[2] = 0;
  m[3] = 0;
  m[4] = 0;
  m[5] = 0;

  float sum, mean, median, clippedmedian, var, dev, skew, kurt;
  median = findMedian(v);
  //  clippedmedian = clippedMedian<1,2>(v);
  float nuse;
  runningStats<0>(v,m,median,mean,dev,nuse);

  cout << "median = " << median << "\n";
  cout << "sum = " << nuse << "\n";
  cout << "mean = " << mean << "\n";
  cout << "standard deviation = " << dev << "\n";
  cout << endl;


//  computeStats(&v[0], &v[0]+v.size(), sum, mean, var, dev, skew, kurt);
//  cout << "count = " << v.size() << "\n";
//  cout << "sum = " << sum << "\n";
//  cout << "median = " << median << "\n";
//  cout << "clippedmedian = " << clippedmedian << "\n";
//  cout << "mean = " << mean << "\n";
//  cout << "variance = " << var << "\n";
//  cout << "standard deviation = " << dev << "\n";
//  cout << "skew = " << skew << "\n";
//  cout << "kurtosis excess = " << kurt << "\n";
//  cout << endl;
}
*/

#endif
