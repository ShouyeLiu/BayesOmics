// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BayesOmics, a statistical genetics software package
// developed by Shouye Liu.
//
// Portions of this file are adapted from GCTB,
// originally licensed under the MIT License. See below for the original license.
//
// Original source: https://cnsgenomics.com/software/gctb/#Download (accessed 20 June 2024)
//
// Modifications and additional code:
// Copyright (C) 2025 Shouye Liu <syliu.xue@foxmail.com>
//
// This file is licensed under the GNU General Public License v3.0 or (at your option)
// any later version. You may redistribute and/or modify it under the terms of the GPL.
//
// BayesOmics is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the LICENSE file for details.


#ifndef toolbox_hpp
#define toolbox_hpp

// #include <sys/time.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <memory>
#include <stdexcept>
#include <array>
#ifdef _WIN32
#include <windows.h>
#endif
#include <sys/stat.h>
#include <unistd.h>
// #include <map>
#include <Eigen/Eigen>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

#include "Logger.hpp"
#include "Macro.hpp"
#include "Stat.hpp"


using namespace std;
using namespace Eigen;
// using namespace arma;

namespace Gadget {

class Timer {
    time_t prev, curr;
public:
    Timer(){
        setTime();
    };
    void setTime(void);
    time_t getTime(void);
    time_t getElapse(void);
    string format(const time_t time);
    string getDate(void);
    void printElapse(void);
};

class Tokenizer : public vector<string> {
    // adopted from matvec
public:
    void getTokens(const string &str, const string &sep);
    int  getIndex(const string &str);
};

template <class T> class Recoder : public map<T,unsigned> {
    // adopted from matvec
    unsigned count;
public:
    Recoder(void){count=0;}
    unsigned code(T s){
        typename map<T,unsigned>::iterator mapit = this->find(s);
        if(mapit == this->end()){
            (*this)[s] = ++count;
            return count;
        }
        else {
            return (*mapit).second;
        }
    }
    void display_codes(ostream & os = cout){
        typename Recoder::iterator it;
        for (it=this->begin(); it!=this->end();it++){
            os << (*it).first << " " << (*it).second << endl;
        }
    }
};

// file process functions
string getFileName(const string &file);
string getFileSuffix(const string &file);
void fileExist(const string &filename);
bool directoryExist(const string &dirname);
bool createDirectory(const string& dirname);

// statistics functions
double calcMean(const VectorXd &vec);
double calcMean(const vector<double> &vec);
double calcMean(const vector<VectorXd> &vec);
VectorXd calColMeans(const MatrixXd &mat);
VectorXd calColMeans(const vector<VectorXd> &vecVec);
VectorXd calRowMeans(const MatrixXd &mat);
double calcVariance(const VectorXd &vec);
double calcVariance(const vector<double> &vec);
double calcVariance(const vector<VectorXd> &vec);
double calcCovariance(const VectorXd &vec1, const VectorXd &vec2);
double calcCorrelation(const VectorXd &vec1, const VectorXd &vec2);
double calcRegression(const VectorXd &y, const VectorXd &x);
double calcPvalFromZ(double z);
double findMedian(const VectorXd &vec);

/// check 
bool checkScaleMatrix(const arma::dmat& Psi, const double nu,const bool messageBool);

bool isSymmetric(const Matrix2d& A, double tol = 1e-6);
Matrix2d makeSymmetricPositiveDefinite(Matrix2d A, double epsilon = 1e-6);
bool isSymmetricPositiveDefinite(const arma::mat& X, double tol = 1e-8);

// shuffle
vector<int> shuffle_index(const int start, const int end);

void shuffle_vector(vector<int> &vec);

std::string formatTime(int seconds);
// Function to display the progress bar with estimated time remaining and color
void showProgressBar(int current, int total, std::chrono::steady_clock::time_point startTime, const std::string& taskName, const std::string& color = BLUE_COLOR);
size_t getTotalLineNumFromFile(const std::string &filename);
size_t getTotalLineNumFromGzFile(const std::string &filename);
// calculate
bool isNumber(const std::string& s);
// is positive number 
bool isPositiveInteger(const std::string& str);
////// other parameter
std::string getSSEvar();
std::string getOSName();
std::string getHostName();

// show some usefull words
std::string getCurrentQuotes();
// git current git hash
std::string getGitCommitHash();

// Function to calculate the mean of the map values
template<typename K, typename V>
double calculateMean(const std::map<K, V>& dataMap) {
    static_assert(std::is_arithmetic<V>::value, "Value type must be numeric");
    double sum = 0.0;
    for (const auto& pair : dataMap) {
        sum += pair.second;
    }
    return sum / dataMap.size();
}

// Function to calculate the standard deviation of the map values
template<typename K, typename V>
double calculateStandardDeviation(const std::map<K, V>& dataMap, double mean) {
    static_assert(std::is_arithmetic<V>::value, "Value type must be numeric");
    double sumOfSquaredDifferences = 0.0;
    for (const auto& pair : dataMap) {
        double difference = pair.second - mean;
        sumOfSquaredDifferences += difference * difference;
    }
    return std::sqrt(sumOfSquaredDifferences / dataMap.size());
}


// jl_array_t* eigen_matrix_to_julia(const Eigen::MatrixXd& matrix);
// Eigen::MatrixXd julia_array_to_eigen_matrix(jl_array_t* result);
}

#endif /* gadgets.hpp */
