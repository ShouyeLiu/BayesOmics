//
//  gadgets.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "Gadgets.hpp"
#include <random>


void Gadget::Tokenizer::getTokens(const string &str, const string &sep){
    clear();
    string::size_type begidx,endidx;
    begidx = str.find_first_not_of(sep);
    while (begidx != string::npos) {
        endidx = str.find_first_of(sep,begidx);
        if (endidx == string::npos) endidx = str.length();
        push_back(str.substr(begidx,endidx - begidx));
        begidx = str.find_first_not_of(sep,endidx);
    }
}

int Gadget::Tokenizer::getIndex(const string &str){
    for (unsigned i=0; i<size(); i++){
        if((*this)[i]==str){
            return i;
        }
    }
    return -1;
}

void Gadget::Timer::setTime(){
    prev = curr = time(0);
}

time_t Gadget::Timer::getTime(){
    return curr = time(0);
}

time_t Gadget::Timer::getElapse(){
    return curr - prev;
}

string Gadget::Timer::format(const time_t time){
    return to_string((long long)(time/3600)) + ":" + to_string((long long)((time % 3600)/60)) + ":" + to_string((long long)(time % 60));
}

string Gadget::Timer::getDate(){
    // return ctime(&curr);
    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);
    auto tm = std::put_time(std::localtime(&now_c), "%T %Z on %a %b %d %Y");
    std::ostringstream oss;
    oss << tm;
    return oss.str();
}

void Gadget::Timer::printElapse(){
    getTime();
    LOGGER << "Time elapse: " << format(getElapse()) << endl;
}

string Gadget::getFileName(const string &file){
    size_t start = file.rfind('/');
    size_t end   = file.rfind('.');
    start = start==string::npos ? 0 : start+1;
    return file.substr(start, end-start);
}

string Gadget::getFileSuffix(const string &file){
    size_t start = file.rfind('.');
    return file.substr(start);
}

void Gadget::fileExist(const string &filename){
    ifstream file(filename.c_str());
    if(!file) LOGGER.e(0," can not open the file ["+filename+"] to read.");
}

bool Gadget::directoryExist(const string& dirname)
{
    struct stat info;

    if (stat(dirname.c_str(), &info) != 0)
        return false;
    else if (info.st_mode & S_IFDIR)
        return true;
    else
        return false;
}

bool Gadget::createDirectory(const string& dirname)
{
    int status = mkdir(dirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (status == -1) {
        std::string errorMsg = strerror(errno);
        if (errno == EEXIST) {
            // The directory already exists
            return true;
        } else if (errno == EACCES) {
            // Permission denied
            LOGGER.e(0,"Permission denied: Unable to create directory '" + dirname + "': " + errorMsg);
            return false;
        } else {
            // Another error occurred
            LOGGER.e(0, "Error: Unable to create directory '" + dirname + "': " + errorMsg);
            return false;
        }
    } else {
        // The directory was created successfully
        return true;
    }
}

vector<int> Gadget::shuffle_index(const int start, const int end){
    vector<int> vec;
    for (unsigned i = start; i <= end; i++) {
        vec.push_back(i);
    }
    
    Gadget::shuffle_vector(vec);

    return vec;
}

void Gadget::shuffle_vector(vector<int> &vec){
    // Get a random integer using Boost random generator (which has been seeded)
    constexpr int max_integer = std::numeric_limits<int>::max();
    int random_integer = static_cast<int>(Stat::ranf() * max_integer);
    
    // Create a thread-local random number generator
    thread_local std::mt19937 rng(random_integer);

    // Shuffle using the thread-local RNG
    std::shuffle(vec.begin(), vec.end(), rng);
}

double Gadget::calcMean(const VectorXd &vec){
    int totalSize = vec.size();
    return vec.mean();
}

double Gadget::calcMean(const vector<double> &vec){
    Eigen::VectorXd vecXd = Eigen::Map<const Eigen::VectorXd>(vec.data(), vec.size());
    return vecXd.mean();
}

double Gadget::calcMean(const vector<VectorXd> &vec){
    int totalSize = 0;
    for (const auto& svec : vec) {
        totalSize += svec.size();
    }
    VectorXd tvec(totalSize); // total vector;
    int currentIndex = 0;
    for (const auto& svec : vec) {
        tvec.segment(currentIndex, svec.size()) = svec;
        currentIndex += svec.size();
    }
    return calcMean(tvec);
}

VectorXd Gadget::calColMeans(const MatrixXd &mat){
    return mat.colwise().mean();
}

VectorXd Gadget::calRowMeans(const MatrixXd &mat){
    return mat.rowwise().mean();
}
VectorXd Gadget::calColMeans(const vector<VectorXd> &vecVec){
    VectorXd means;
    means.setZero(vecVec.size());
    for(unsigned i = 0; i < vecVec.size(); i++ ){
        means(i) = vecVec[i].mean();
    }
    return means;
}

double Gadget::calcVariance(const VectorXd &vec){
    return (vec.array() - vec.mean()).square().sum()/ (vec.size() -1 );
}

double Gadget::calcVariance(const vector<double> &vec){
    Eigen::VectorXd vecXd = Eigen::Map<const Eigen::VectorXd>(vec.data(), vec.size());
    return (vecXd.array() - vecXd.mean()).square().sum()/ (vecXd.size() -1 );
}

double Gadget::calcVariance(const vector<VectorXd> &vec){
    int totalSize = 0;
    for (const auto& svec : vec) {
        totalSize += svec.size();
    }
    VectorXd tvec(totalSize); // total vector;
    int currentIndex = 0;
    for (const auto& svec : vec) {
        tvec.segment(currentIndex, svec.size()) = svec;
        currentIndex += svec.size();
    }
    return calcVariance(tvec);
}

double Gadget::calcCovariance(const VectorXd &vec1, const VectorXd &vec2){
    if (vec1.size() != vec2.size()) {
        LOGGER.e(0," Gadget::calcCovariance: the two vectors have different sizes.");
    }
    VectorXd vec1_double = vec1.cast<double>();
    VectorXd vec2_double = vec2.cast<double>();
    return (vec1_double.array()-vec1_double.mean()).cwiseProduct(vec2_double.array()-vec2_double.mean()).sum()/vec1_double.size();
}

double Gadget::calcCorrelation(const VectorXd &vec1, const VectorXd &vec2){
    double cov = calcCovariance(vec1, vec2);
    double var1 = calcVariance(vec1);
    double var2 = calcVariance(vec2);
    return cov/sqrt(var1*var2);
}

double Gadget::calcRegression(const VectorXd &y, const VectorXd &x){
    double cov = calcCovariance(y, x);
    double varx = calcVariance(x);
    return cov/varx;
}

double Gadget::calcPvalFromZ(double z){
    // double p = 0.5 * erfc(-z / sqrt(2));
    double p = 1.0 - 0.5 * (1.0 + erf(std::abs(z) / sqrt(2.0)));
    p = 2.0 * p;
    return p;
}

double Gadget::findMedian(const VectorXd &vec){
    VectorXd tmp = vec;
    std::sort(tmp.data(), tmp.data() + tmp.size());
    // return tmp[tmp.size()/2];
    return tmp.size() % 2 == 0 ? tmp.segment( (tmp.size()-2)/2, 2 ).mean() : tmp( tmp.size()/2 );
}

bool Gadget::checkScaleMatrix(const arma::dmat& Psi,const double nu, const bool messageBool) {
    // // Check if Psi is square
    // if (Psi.n_rows != Psi.n_cols) {
    //     if(messageBool) LOGGER.e(0,"Scale matrix must be square.");
    //     return false;
    // }
    // // Check if Psi is symmetric
    // if (!Psi.is_symmetric()) {
    //     if(messageBool) LOGGER.e(0,"Scale matrix must be symmetric.");
    //     return false;
    // }

    // // Check if Psi is positive definite
    // if (!Psi.is_sympd()) {
    //     if(messageBool) LOGGER.e(0,"Scale matrix must be positive definite.");
    //     return false;
    // }

    // // Check degrees of freedom
    // if (nu <= Psi.n_rows - 1) {
    //     if(messageBool) LOGGER.e(0,"Degrees of freedom must be greater than the number of rows/columns of the scale matrix minus one.");
    //     return false;
    // }

    // If all checks pass
    return true;
}

// Function to check if a matrix is symmetric
bool Gadget::isSymmetric(const Matrix2d& A, double tol) {
    return A.isApprox(A.transpose(), tol);
}


using namespace arma;
using namespace std;

// Function to check if a matrix is symmetric positive definite
bool Gadget::isSymmetricPositiveDefinite(const arma::mat& X, double tol) {
    // Step 1: Check if the matrix is symmetric
    if (!arma::approx_equal(X, X.t(), "absdiff", tol)) {
        return false;
    }

    // Step 2: Attempt Cholesky decomposition to check positive definiteness
    arma::mat L;
    bool pd = arma::chol(L, X);
    return pd;
}



// Function to make a matrix symmetric positive definite (SPD)
Matrix2d Gadget::makeSymmetricPositiveDefinite(Matrix2d A, double epsilon) {
   // Step 1: Make the matrix symmetric by averaging with its transpose
    if (!A.isApprox(A.transpose(), epsilon)) {
        A = (A + A.transpose()) / 2;
    }

    // Step 2: Perform SVD decomposition
    JacobiSVD<Matrix2d> svd(A, ComputeFullU | ComputeFullV);
    Vector2d singularValues = svd.singularValues();
    Matrix2d U = svd.matrixU();
    Matrix2d V = svd.matrixV();

    // Step 3: Adjust singular values to ensure they are all positive
    for (int i = 0; i < singularValues.size(); ++i) {
        if (singularValues(i) <= 0) {
            singularValues(i) = epsilon;
        }
    }

    // Step 4: Reconstruct the SPD matrix
    Matrix2d A_spd = U * singularValues.asDiagonal() * V.transpose();
    return A_spd;
}



// value
bool Gadget::isNumber(const std::string& s)
{
    return !s.empty() && std::find_if(s.begin(), 
        s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
}
// positive interger
bool Gadget::isPositiveInteger(const std::string& str) {
    // A positive integer must not be empty and must consist of digits only.
    // It also must not start with '0' unless it is exactly "0".
    return !str.empty() && 
           std::all_of(str.begin(), str.end(), ::isdigit) &&
           (str.size() == 1 || str[0] != '0');
}

// Function to convert seconds into a human-readable format
std::string Gadget::formatTime(int seconds) {
    if (seconds < 60) {
        return std::to_string(seconds) + "s";
    } else if (seconds < 3600) { // Less than an hour
        int minutes = seconds / 60;
        int remainingSeconds = seconds % 60;
        return std::to_string(minutes) + "m " + std::to_string(remainingSeconds) + "s";
    } else { // 1 hour or more
        int hours = seconds / 3600;
        int remainingMinutes = (seconds % 3600) / 60;
        int remainingSeconds = seconds % 60;

        std::string formattedTime = std::to_string(hours) + "h";
        if (remainingMinutes > 0 || remainingSeconds > 0) {
            formattedTime += " " + std::to_string(remainingMinutes) + "m";
        }
        if (remainingSeconds > 0) {
            formattedTime += " " + std::to_string(remainingSeconds) + "s";
        }
        return formattedTime;
    }
}

// Function to display the progress bar with estimated time remaining and color
void Gadget::showProgressBar(int current, int total, std::chrono::steady_clock::time_point startTime, const std::string& taskName, const std::string& color) {
    int barWidth = 50;
    float progress = static_cast<float>(current) / total;

    // Calculate elapsed time
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - startTime).count();

    // Estimate total time and remaining time
    int estimatedTotalTime = static_cast<int>(elapsed / progress);
    int remainingTime = estimatedTotalTime - elapsed;

    // Print task name with color
    std::cout << color << taskName << ": ";

    // Display the progress bar with color
    std::cout << "[";
    int pos = static_cast<int>(barWidth * progress);
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "% ";
    // Display elapsed and remaining time
    std::cout << "Elapsed: " << formatTime(elapsed) << "s"
              << " | Remaining: " << formatTime(remainingTime) << "s" << RESET_COLOR << "\r";
    std::cout.flush();
}

size_t Gadget::getTotalLineNumFromFile(const std::string &filename) {
    FILE *file = fopen(filename.c_str(), "rb");
    if (!file) {
        LOGGER.e(0,"Error: Cannot open file [ "+ filename + "].");
        return 0;
    }
    const size_t bufferSize = 65536;  // 64 KB buffer size for faster reading
    char *buffer = new char[bufferSize]; // Dynamically allocate the buffer
    size_t totalLines = 0;

    while (size_t bytesRead = fread(buffer, 1, bufferSize, file)) {
        // Process each character in the buffer
        for (size_t i = 0; i < bytesRead; ++i) {
            if (buffer[i] == '\n') {
                ++totalLines;
            }
        }
    }
    fclose(file);  // Close the file after reading
    delete[] buffer; // Free the dynamically allocated buffer
    return totalLines;
}

size_t Gadget::getTotalLineNumFromGzFile(const std::string &filename) {
    std::ifstream file(filename, std::ios_base::in | std::ios_base::binary);
    if (!file) {
        LOGGER.e(0,"Cannot open file [" + filename + "].");
    }
    // Create a filtering stream for decompression
    boost::iostreams::filtering_stream<boost::iostreams::input> in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(file);

    const size_t bufferSize = 65536;  // 64 KB buffer size for faster reading
    char buffer[bufferSize];
    size_t totalLines = 0;

    // Reading and counting lines on-the-fly using OpenMP
    #pragma omp parallel
    {
        size_t localLines = 0;
        #pragma omp single nowait
        {
            while (in.read(buffer, bufferSize)) {
                std::streamsize bytesRead = in.gcount();

                // Process the buffer in parallel
                #pragma omp parallel for reduction(+:localLines)
                for (std::streamsize i = 0; i < bytesRead; ++i) {
                    if (buffer[i] == '\n') {
                        ++localLines;
                    }
                }
            }
        }

        // Combine the local line counts into the total count
        #pragma omp atomic
        totalLines += localLines;
    }
    return totalLines;
}


std::string Gadget::getSSEvar(){
        #ifdef __AVX512F__
    return "avx512";
    #elif __AVX2__
    return "avx2";
    #elif __AVX__
    return "avx";
    #elif __SSE4_1
    return "sse4";
    #elif __SSE2__
    return "sse2";
    #else
    return "unknown";
    #endif
}

std::string Gadget::getOSName(){
    #ifdef _WIN64
    return "Windows";
    #elif __linux__
    return "Linux";
    #elif __APPLE__ || __MACH__
    return "Mac";
    #else
    return "Other";
    #endif
}

std::string Gadget::getHostName(){
    char *temp = NULL;
    std::string computerName;

#if defined(WIN32) || defined(_WIN32) || defined(_WIN64)
    temp = getenv("COMPUTERNAME");
    if (temp != NULL) {
        computerName = temp;
    }
#else
    temp = getenv("HOSTNAME");
    if (temp != NULL) {
        computerName = temp;
    } else {
        temp = new char[512];
        if (gethostname(temp, 512) == 0) {
            computerName = temp;
        }
        delete []temp;
    }
#endif
    return computerName;
}


std::string Gadget::getGitCommitHash() {
    std::array<char, 128> buffer;
    std::string result;
    #ifdef _WIN32
    // On Windows, use _popen and _pclose
    std::unique_ptr<FILE, decltype(&_pclose)> pipe(_popen("git rev-parse HEAD", "r"), _pclose);
    #else
    // On Unix-like systems, use popen and pclose
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen("git rev-parse HEAD", "r"), pclose);
    #endif
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    // Remove the newline character from the end of the commit hash
    if (!result.empty() && result[result.length() - 1] == '\n') {
        result.erase(result.length() - 1);
    }
    return result;
}

std::string Gadget::getCurrentQuotes(){
    vector<string> quotes = {
        "From where we stand, the rain seems ran#pragma endregiondom\n. If we could stand somewhere else, we would see the order in it.\n --T. Hillerman (1990) ",
        "\"You can never know everything,\" Lan said quietly, \"and part of what you know is always wrong. Perhaps even the most important part."
        "A portion of wisdom lies in knowing that. A portion of courage lies in going on anyway.\"\n -- Robert Jordan",
        "Sometimes the Pattern has a randomness to it—to our eyes, at least—but what chance that you should meet a man who"
        "could guide you in this thing, and he one who could follow the guiding?\n -- Robert Jordan",
        "Everything should be made as simple as possible, but not simpler",
        "Models should be as simple as possible, but not more so.\n --Einstein",
        "In our lust for measurement, we frequently measure that which we can rather than that which we wish to measure... and forget that there is a difference.\n --George Udny Yule",
        "We must be careful not to confuse data with the abstractions we use to analyze them.\n --William James",
        "Natural selection is a mechanism for generating an exceedingly high degree of improbability.\n --R. A. Fisher",
        "Facts speak louder than statistics\n --Mr. Justice Streatfield",
        "An approximate answer to the right problem is worth a good deal more than an exact answer to an approximate problem.\n -- John Tukey",
        "Randomization is too important to be left to chance.\n --J. D. Petruccelli"
        "Probability theory is nothing but common sense reduced to calculation. — Pierre Laplace"
    };

    std::random_device rd;  
    std::mt19937 gen(rd());
    int min = 1; // Lower bound of the range (inclusive)
    int max = quotes.size(); // Upper bound of the range (inclusive)
    std::uniform_int_distribution<> distr(min, max);
    string currentQuote = quotes[distr(gen) -1];
    return currentQuote;
}




// // Function to convert Eigen::MatrixXd to Julia array
// jl_array_t* Gadget::eigen_matrix_to_julia(const Eigen::MatrixXd& matrix) {
//     // Get dimensions
//     size_t rows = matrix.rows();
//     size_t cols = matrix.cols();

//     // Create a new Julia array
//     jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_float64_type, 2);
//     jl_array_t* julia_array = jl_alloc_array_2d(array_type, rows, cols);

//     // Copy data from Eigen matrix to Julia array
//     double* julia_array_data = (double*)jl_array_data(julia_array);
//     std::memcpy(julia_array_data, matrix.data(), rows * cols * sizeof(double));

//     return julia_array;
// }

// // Assuming 'result' is the Julia array returned by your function
// Eigen::MatrixXd Gadget::julia_array_to_eigen_matrix(jl_array_t* result) {
//     size_t rows = jl_array_dim(result, 0);
//     size_t cols = jl_array_dim(result, 1);
//     double* data = reinterpret_cast<double*>(jl_array_data(result));
//     Eigen::Map<Eigen::MatrixXd> eigen_matrix(data, rows, cols);
//     return eigen_matrix;
// }
