#include <regex>
#include <unordered_map>
#include "MemoryBudget.hpp"
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BayesOmics, a statistical genetics software package
// developed by Shouye Liu.
//
// Copyright (C) 2025 Shouye Liu <syliu.xue@foxmail.com>
//
// BayesOmics is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// BayesOmics is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with BayesOmics. If not, see <https://www.gnu.org/licenses/>.


#include "Data.hpp"
#include <sys/stat.h>
#include <dirent.h>
// #include "H5Cpp.h"
#include "hdf5.h"
// using namespace H5;

///////////////////////////////////////////////////////////////////////////////////////
////////    Step 1. perform eigen-decomposition for ld blocks                   ///////
///////////////////////////////////////////////////////////////////////////////////////
void Data::readLDBlockInfoFile(const string &ldBlockInfoFile) {
    ifstream in(ldBlockInfoFile.c_str());
    if (!in) {
        LOGGER.e(0, "can not open the file [" + ldBlockInfoFile + "] to read.");
    }
    LOGGER << "Reading ld block info from file [" << ldBlockInfoFile << "]." << endl;

    ldBlockInfoVec.clear();
    ldBlockInfoMap.clear();

    string header;
    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t\n");

    getline(in, header);
    colData.getTokens(header, sep);

    const size_t ncol = colData.size();
    string firstHeader = ncol ? colData[0] : "";
    std::transform(
        firstHeader.begin(), firstHeader.end(), firstHeader.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); }
    );
    const bool isFourColumnBlock =
        ncol == 4 &&
        (firstHeader == "block" || firstHeader == "ldblock" ||
         firstHeader == "blockid");
    const bool mergeIntervals = ncol == 4 && !isFourColumnBlock;

    if (ncol == 3) {
        LOGGER << "LDdetect format is used." << endl;
    } else if (isFourColumnBlock) {
        LOGGER << "ref4cM LD block format is used." << endl;
    } else if (ncol == 4) {
        LOGGER << "Gene interval format is used." << endl;
    } else {
        LOGGER.e(0, "Wrong format.");
    }

    struct RawBlock {
        string id;
        int chr;
        int start;
        int stop;
    };

    vector<RawBlock> rawBlocks;
    int autoId = 1;
    while (getline(in, inputStr)) {
        if (inputStr.empty()) continue;
        colData.getTokens(inputStr, sep);
        if (colData.size() == 0) continue;
        RawBlock block;
        if (ncol == 3) {
            if (colData.size() != 3) {
                LOGGER.e(0, "Wrong number of columns in line: " + inputStr);
            }
            block.id = std::to_string(autoId++);
            block.chr = std::stoi(colData[0].substr(3));
            block.start = static_cast<int>(std::stod(colData[1]));
            block.stop = static_cast<int>(std::stod(colData[2]));
        } else {
            if (colData.size() != 4) {
                LOGGER.e(0, "Wrong number of columns in line: " + inputStr);
            }
            block.id = colData[0];
            block.chr = std::stoi(colData[1]);
            block.start = static_cast<int>(std::stod(colData[2]));
            block.stop = static_cast<int>(std::stod(colData[3]));
        }
        if (block.start > block.stop) {
            std::swap(block.start, block.stop);
        }
        rawBlocks.push_back(block);
    }
    in.close();
    // Sort by chromosome first, then by physical position within chromosome.
    std::sort(
        rawBlocks.begin(),
        rawBlocks.end(),
        [](const RawBlock &a, const RawBlock &b) {
            if (a.chr != b.chr) return a.chr < b.chr;
            if (a.start != b.start) return a.start < b.start;
            if (a.stop != b.stop) return a.stop < b.stop;
            return a.id < b.id;
        }
    );
    vector<RawBlock> outputBlocks;
    if (!mergeIntervals) {
        // Each row in an LD block file defines one block, even when two
        // blocks share a boundary or overlap.
        outputBlocks = rawBlocks;
    } else if (!rawBlocks.empty()) {
        // Gene intervals can overlap, so combine overlapping regions on the
        // same chromosome before constructing their LD matrices.
        RawBlock current = rawBlocks.front();
        for (size_t i = 1; i < rawBlocks.size(); ++i) {
            const RawBlock &next = rawBlocks[i];
            const bool sameChr = next.chr == current.chr;
            const bool overlap = next.start <= current.stop;
            if (sameChr && overlap) {
                current.start = std::min(current.start, next.start);
                current.stop = std::max(current.stop, next.stop);
            } else {
                outputBlocks.push_back(current);
                current = next;
            }
        }
        outputBlocks.push_back(current);
    }
    int idx = 0;
    for (const RawBlock &block : outputBlocks) {
        const string blockId = isFourColumnBlock
            ? block.id
            : std::to_string(idx + 1);

        LDBlockInfo *ld = new LDBlockInfo(
            idx,
            blockId,
            block.chr
        );
        ld->startPos = block.start;
        ld->endPos = block.stop;
        ldBlockInfoVec.push_back(ld);
        if (!ldBlockInfoMap.insert({blockId, ld}).second) {
            LOGGER.e(
                0,
                "Duplicate LD block ID found: \"" +
                blockId + "\"."
            );
        }

        idx++;
    }
    numLDBlocks = static_cast<unsigned>(ldBlockInfoVec.size());
    LOGGER << rawBlocks.size() << " input intervals read; "
           << numLDBlocks << " LD Blocks retained"
           << (mergeIntervals ? " after overlap merging" : "") << " from ["
           << ldBlockInfoFile
           << "]."
           << endl;
}

void Data::eigenDecomposition(const MatrixXf &X, const float &prop, VectorXf &values, MatrixXf &vectors, float &sumPositive) {
    if (X.rows() == 0 || X.rows() != X.cols() || !X.allFinite() || !(prop > 0 && prop <= 1))
        throw std::runtime_error("Invalid LD matrix or eigenvalue variance proportion (expected 0 < cutoff <= 1)");
    const MatrixXd symmetric = .5 * (X.cast<double>() + X.transpose().cast<double>());
    SelfAdjointEigenSolver<MatrixXd> solver(symmetric);
    if (solver.info() != Eigen::Success) throw std::runtime_error("LD eigen decomposition failed");
    const auto &lambda = solver.eigenvalues();
    if (lambda[lambda.size()-1] <= 0 || lambda[0] < -1e-4 * lambda[lambda.size()-1])
        throw std::runtime_error("LD matrix has no positive variance or is materially indefinite");
    int first = 0;
    while (first < lambda.size() && lambda[first] <= 1e-10) ++first;
    const double total = lambda.tail(lambda.size()-first).sum();
    int start = lambda.size(); double retained = 0;
    // Preserve the on-disk ascending order; include the minimal leading modes.
    while (start > first && (start == lambda.size() || prop == 1.f ||
        retained + 2*std::numeric_limits<float>::epsilon()*total < double(prop)*total)) retained += lambda[--start];
    values = lambda.tail(lambda.size()-start).cast<float>();
    vectors = solver.eigenvectors().rightCols(values.size()).cast<float>();
    sumPositive = static_cast<float>(total);
}

float computeAdaptiveLambda(const Eigen::MatrixXf& R, float eps = 1e-5f) {
    // Perform eigen decomposition for symmetric matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(R);
    // Fallback in case decomposition fails (rare but safe guard)
    if (es.info() != Eigen::Success) {
        return 1e-3f;
    }
    // Smallest eigenvalue (since eigenvalues are sorted in ascending order)
    float lmin = es.eigenvalues()(0);
    // If matrix is already sufficiently positive definite
    if (lmin >= eps) {
        // Return very small shrinkage (almost no change)
        return 1e-6f;
    }
    // Compute minimal lambda such that:
    // (1 - lambda) * lmin + lambda >= eps
    // => lambda >= (eps - lmin) / (1 - lmin)
    float denom = 1.0f - lmin;
    float lambda;
    if (denom <= 1e-8f) {
        // Numerical edge case (lmin ~ 1)
        lambda = 1e-3f;
    } else {
        lambda = (eps - lmin) / denom;
    }
    // Add small safety margin to avoid numerical instability
    lambda *= 1.05f;
    // Clamp lambda to a reasonable range to avoid over-shrinkage
    if (lambda < 1e-6f) lambda = 1e-6f;
    if (lambda > 0.1f)  lambda = 0.1f;
    return lambda;
}
////////////////////////////////////////////////////
// Shrink towards identity
////////////////////////////////////////////////////
void shrinkLD(MatrixXf& R, float lambda) {
    int n = R.rows();

    R *= (1.0f - lambda);

    for (int i = 0; i < n; ++i) {
        R(i, i) += lambda;
    }
}

MatrixXd Data::generateLDmatrixPerBlock(const std::string &bedFile, const std::vector<std::string> &ids, bool covariance) {
    if (ids.empty() || numKeptInds < 2) throw std::runtime_error("LD requires SNPs and at least two individuals");
    const size_t stride = (uint64_t(numInds)+3)/4;
    std::ifstream bed(bedFile, std::ios::binary);
    unsigned char header[3]; bed.read(reinterpret_cast<char*>(header),3);
    if (!bed || header[0]!=0x6c || header[1]!=0x1b || header[2]!=1)
        throw std::runtime_error("Invalid SNP-major BED header: "+bedFile);
    // BED positions are raw BIM indices; requested LD order may differ.
    std::unordered_map<std::string,size_t> rawIndex;
    for(size_t j=0;j<snpInfoVec.size();++j) rawIndex.emplace(snpInfoVec[j]->rsID,j);
    MatrixXd x(numKeptInds,ids.size());
    std::vector<unsigned char> packed(stride);
    const double decode[4]={2.,-9.,1.,0.};
    for(size_t j=0;j<ids.size();++j) {
        auto found=rawIndex.find(ids[j]);
        if(found==rawIndex.end()) throw std::runtime_error("LD SNP missing from BIM: "+ids[j]);
        bed.seekg(3+uint64_t(found->second)*stride);
        bed.read(reinterpret_cast<char*>(packed.data()),stride);
        if(!bed) throw std::runtime_error("Truncated BED while reading "+ids[j]);
        double sum=0; size_t observed=0;
        for(size_t i=0;i<indInfoVec.size();++i) if(indInfoVec[i]->kept) {
            const double value=decode[(packed[i/4] >> (2*(i%4)))&3];
            x(indInfoVec[i]->index,j)=value;
            if(value!=-9){sum+=value;++observed;}
        }
        if(observed<2) throw std::runtime_error("Insufficient observed genotypes for "+ids[j]);
        const double mean=sum/observed;
        for(size_t i=0;i<numKeptInds;++i) x(i,j)=x(i,j)==-9 ? 0 : x(i,j)-mean;
        const double ss=x.col(j).squaredNorm();
        if(!(ss>0)) throw std::runtime_error("Monomorphic SNP in LD reference: "+ids[j]);
        auto snp=snpInfoVec[found->second];
        snp->af=mean/2; snp->twopq=2*snp->af*(1-snp->af); snp->sampleSize=numKeptInds;
        if(!covariance) x.col(j)/=std::sqrt(ss);
    }
    // Exact sample correlation, without implicit ridge/shrinkage. Double
    // accumulation also avoids the N-dependent float summation error.
    MatrixXd result=x.transpose()*x;
    if(covariance) result/=numKeptInds; // established population covariance convention
    else result.diagonal().setOnes();
    return result;
}

// MatrixXd Data::generateLDmatrixPerBlock(const string &bedFile, const vector<string> &snplists,const bool isCovBool){
//     int numSnpInRange = snplists.size();
//     IndInfo *indi = NULL;
//     SnpInfo *snpj = NULL;
//     SnpInfo *snpk = NULL;
//     snpj = snpInfoMap.at(snplists[0]); // start
//     snpk = snpInfoMap.at(snplists[numSnpInRange - 1]); // end;
//     unsigned start = snpj->index;
//     unsigned end = snpk->index;
//     if (numIncdSnps == 0) LOGGER.e(0,"No SNP is retained for analysis.");
//     if (numKeptInds == 0) LOGGER.e(0,"No individual is retained for analysis.");
//     // if (start >= numIncdSnps) LOGGER.e(0,"Specified a SNP range of " + snpRange + " but " + to_string(static_cast<long long>(numIncdSnps)) + " SNPs are included.");
//     //////////////////////////////////////////////////////
//     // Step 1. read in the genotypes of SNPs in the given range
//     //////////////////////////////////////////////////////
//     const int bedToGeno[4] = {2, -9, 1, 0};
//     unsigned size = (numInds+3)>>2;
//     MatrixXd ZP(numSnpInRange, numKeptInds);  // SNP x Ind
//     VectorXd Dtmp;
//     Dtmp.setZero(numSnpInRange);
//     if (numKeptInds < 2) LOGGER.e(0, " Cannot calculate LD matrix with number of individuals < 2.");
//     FILE *in1 = fopen(bedFile.c_str(), "rb");
//     if (!in1) LOGGER.e(0, " can not open the file [" + bedFile + "] to read.");
//     // cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
//     char header[3];
//     fread(header, sizeof(header), 1, in1);
//     if (!in1 || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
//         cerr << "Error: Incorrect first three bytes of bed file: " << bedFile << endl;
//         exit(1);
//     }
//     int genoValue;
//     unsigned i, j, k;
//     unsigned incj, inck; // index of included SNP
//     unsigned long long skipj = 0;
//     unsigned nmiss;
//     float mean;
//     set<int> chromInRange;
//     for (j = 0, incj = 0; j < numSnps; j++) {
//         snpj = snpInfoVec[j];
//         if (snpj->index < start || !snpj->included) {
//             skipj += size;
//             continue;
//         }
//         // check if snp exist in snplist in case gap situation
//         auto it = std::find(snplists.begin(), snplists.end(), snpj->rsID);
//         if (it == snplists.end()) {
//             skipj += size;
//             continue;
//         } 
//         if (skipj) fseek(in1, skipj, SEEK_CUR);
//         skipj = 0;
//         char *bedLineIn = new char[size];
//         fread(bedLineIn, sizeof(char), size, in1);
//         chromInRange.insert(snpj->chrom);
//         mean = 0.0;
//         nmiss = 0;
//         for (i = 0; i < numInds; i++) {
//             indi = indInfoVec[i];
//             if (!indi->kept) continue;
//             genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
//             ZP(incj, indi->index) = genoValue;
//             if (genoValue == -9) ++nmiss;
//             else mean += genoValue;
//         }
//         delete[] bedLineIn;
//         // fill missing values with the mean
//         snpj->sampleSize = numKeptInds-nmiss;
//         mean /= float(snpj->sampleSize);
//         if (nmiss) {
//             for (i=0; i<numKeptInds; ++i) {
//                 if (ZP(incj, i) == -9) ZP(incj, i) = mean;
//             }
//         }
//         // compute allele frequency
//         snpj->af = 0.5f*mean;
//         snp2pq[incj] = snpj->twopq = 2.0f*snpj->af*(1.0f-snpj->af);
//         if (snp2pq[incj]==0) LOGGER.e(0, " " + snpj->rsID + " is a fixed SNP (MAF=0)!");
//         Dtmp[incj] = Gadget::calcVariance(ZP.row(incj)); // *numKeptInds;
//         // ZP.row(incj) = (ZP.row(incj).array() - ZP.row(incj).mean())/sqrt(Dtmp[incj]);

//         ZP.row(incj) = (ZP.row(incj).array() - ZP.row(incj).mean());
//         if(!isCovBool) ZP.row(incj) = ZP.row(incj)/sqrt(Dtmp[incj]);

//         if (++incj == numSnpInRange) break;
//     }
//     fclose(in1);
//     ZPZdiag = ZP.rowwise().squaredNorm();
//     //////////////////////////////////////////////////////
//     // Step 2. read in the bed file again to compute Z'Z
//     //////////////////////////////////////////////////////
//     MatrixXd denseZPZ;
//     denseZPZ.setZero(numSnpInRange, numSnpInRange);
//     VectorXd Zk(numKeptInds);
//     Dtmp.setZero(numSnpInRange);
//     FILE *in2 = fopen(bedFile.c_str(), "rb");
//     fseek(in2, 3, SEEK_SET);
//     unsigned long long skipk = 0;
//     set<int>::iterator setend = chromInRange.end();
//     if (numSkeletonSnps) {
//         for (k = 0, inck = 0; k < numSnps; k++) {
//             snpk = snpInfoVec[k];
//             // if (!snpk->included) {
//             if (snpk->index < start || !snpk->included) {
//                 skipk += size;
//                 continue;
//             }
//             // check if snp exist in snplist in case gap situation
//             auto it = std::find(snplists.begin(), snplists.end(), snpk->rsID);
//             if (it == snplists.end()) {
//                 skipk += size;
//                 continue;
//             } 
//             if (chromInRange.find(snpk->chrom) == setend && !snpk->skeleton) {
//                 skipk += size;
//                 ++inck;       // ensure the index is correct
//                 continue;
//             }
//             if (skipk) fseek(in2, skipk, SEEK_CUR);
//             skipk = 0;
//             char *bedLineIn = new char[size];
//             fread(bedLineIn, sizeof(char), size, in2);
//             mean = 0.0;
//             nmiss = 0;
//             for (i = 0; i < numInds; i++) {
//                 indi = indInfoVec[i];
//                 if (!indi->kept) continue;
//                 genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
//                 Zk[indi->index] = genoValue;
//                 if (genoValue == -9) ++nmiss;   // missing genotype
//                 else mean += genoValue;
//             }
//             delete[] bedLineIn;
//             // fill missing values with the mean
//             snpk->sampleSize = numKeptInds-nmiss;
//             mean /= float(snpk->sampleSize);
//             if (nmiss) {
//                 for (i=0; i<numKeptInds; ++i) {
//                     if (Zk[i] == -9) Zk[i] = mean;
//                 }
//             }
//             // compute allele frequency
//             snpk->af = 0.5f*mean;
//             snp2pq[inck] = snpk->twopq = 2.0f*snpk->af*(1.0f-snpk->af);
//             if (snp2pq[inck]==0) LOGGER.e(0, " " + snpk->rsID + " is a fixed SNP (MAF=0)!");
//             Dtmp[inck] = Gadget::calcVariance(Zk.row(inck)); //*numKeptInds;   // calculate  variances fro each snp;
//             // Zk = (Zk.array() - Zk.mean())/sqrt(Dtmp[inck]); // center and standardize genotype
//             Zk = (Zk.array() - Zk.mean());
//             if(!isCovBool) Zk = Zk/sqrt(Dtmp[inck]);

//             denseZPZ.col(inck) = ZP * Zk /numKeptInds;
//             if(isCovBool) denseZPZ.col(inck) = denseZPZ.col(inck)/numKeptInds;
//             //++inck;
//             if (++inck == numSnpInRange) break;
//         }
//     }
//     else {
//         for (k = 0, inck = 0; k < numSnps; k++) {
//             snpk = snpInfoVec[k];
//             // if (!snpk->included) {
//             if (snpk->index < start || !snpk->included) {
//                 skipk += size;
//                 continue;
//             }
//             // check if snp exist in snplist in case gap situation
//             auto it = std::find(snplists.begin(), snplists.end(), snpk->rsID);
//             if (it == snplists.end()) {
//                 skipk += size;
//                 continue;
//             } 
//             if (skipk) fseek(in2, skipk, SEEK_CUR);
//             skipk = 0;
//             char *bedLineIn = new char[size];
//             fread(bedLineIn, sizeof(char), size, in2);
//             mean = 0.0;
//             nmiss = 0;
//             for (i = 0; i < numInds; i++) {
//                 indi = indInfoVec[i];
//                 if (!indi->kept) continue;
//                 genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
//                 Zk[indi->index] = genoValue;
//                 if (genoValue == -9) ++nmiss;   // missing genotype
//                 else mean += genoValue;
//             }
//             delete[] bedLineIn;
//             // fill missing values with the mean
//             snpk->sampleSize = numKeptInds-nmiss;
//             mean /= float(snpk->sampleSize);
//             if (nmiss) {
//                 for (i=0; i<numKeptInds; ++i) {
//                     if (Zk[i] == -9) Zk[i] = mean;
//                 }
//             }
//             // compute allele frequency
//             snpk->af = 0.5f*mean;
//             snp2pq[inck] = snpk->twopq = 2.0f*snpk->af*(1.0f-snpk->af);
//             if (snp2pq[inck]==0) LOGGER.e(0, " " + snpk->rsID + " is a fixed SNP (MAF=0)!");
//             Dtmp[inck] = Gadget::calcVariance(Zk); //*numKeptInds;
//             // Zk = (Zk.array() - Zk.mean())/sqrt(Dtmp[inck]);
//             Zk = (Zk.array() - Zk.mean());
//             if(!isCovBool) Zk = Zk/sqrt(Dtmp[inck]);
            
//             denseZPZ.col(inck) = ZP * Zk /numKeptInds;
//             if(isCovBool) denseZPZ.col(inck) = denseZPZ.col(inck)/numKeptInds;
//             if (++inck == numSnpInRange) break;
//         }
//     }
//     fclose(in2);
//     for (k = 0, inck = 0; k < numSnps; k++) {
//         snpk = snpInfoVec[k];
//         // if (!snpk->included) {
//         if (snpk->index < start || !snpk->included) {
//             skipk += size;
//             continue;
//         }
//         if (skipk) fseek(in2, skipk, SEEK_CUR);
//        skipk = 0;
//        inck ++; 
//     }
//     cout << endl;
//     return denseZPZ;
// }


void Data::getEigenDataFromFullLDM(const string &filename, const float eigenCutoff){
    
    string outfilename = filename + ".eigen.bin";
    FILE *out3 = fopen(outfilename.c_str(), "wb");

    for (unsigned blk=0; blk < numKeptLDBlocks; blk++){
        
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal;
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(ZPZmat.cast<float>(), eigenCutoff, eigenVal, eigenVec, sumPosEigVal);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix
        
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = keptLdBlockInfoVec[blk]->numSnpInBlock;  // TMP

        cout << " Generate Eigen decomposition result for LD block " << blk << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << endl;
        
        // save summary
        // 1, the number of SNPs in the block
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out3);
        // 2, the number of eigenvalues at with the given cutoff
        fwrite(&numEigenValue, sizeof(int32_t), 1, out3);
        // 3. sum of all the positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, out3);
        // 4. eigenvalue cutoff based on the proportion of variance explained in LD
        fwrite(&eigenCutoff, sizeof(float), 1, out3);
        // 5. the selected eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out3);
        // 6. the selected eigenvector;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, out3);
        
    }
    
    fclose(out3);

    //outputEigenDataForLDM(filename);
    
}

void Data::mapSnpsToBlocks(int window) {
    std::map<int,std::vector<SnpInfo*>> chromosomes;
    for(auto snp:incdSnpInfoVec) {chromosomes[snp->chrom].push_back(snp);snp->isInBlock=false;}
    for(auto &chr:chromosomes) std::stable_sort(chr.second.begin(),chr.second.end(),
        [](SnpInfo* a,SnpInfo* b){return a->physPos<b->physPos;});
    for(auto block:ldBlockInfoVec) {
        block->snpNameVec.clear(); block->snpInfoVec.clear(); block->kept=false;
        const auto &snps=chromosomes[block->chrom];
        const int64_t start=int64_t(block->startPos)-window, end=int64_t(block->endPos)+window;
        auto it=std::lower_bound(snps.begin(),snps.end(),start,[](SnpInfo* s,int64_t bp){return s->physPos<bp;});
        // ref4cM blocks are half-open [start,end), so a boundary SNP has one owner.
        for(;it!=snps.end() && (*it)->physPos<end;++it) {
            auto snp=*it;
            if(snp->isInBlock) throw std::runtime_error("Overlapping LD intervals assign SNP twice: "+snp->rsID);
            block->snpNameVec.push_back(snp->rsID); block->snpInfoVec.push_back(snp);
            snp->isInBlock=true; snp->block=block->ID;
        }
        block->numSnpInBlock=block->snpInfoVec.size();
        if(block->numSnpInBlock) {
            block->kept=true;block->startSnpIdx=block->snpInfoVec.front()->index;
            block->endSnpIdx=block->snpInfoVec.back()->index;
        }
    }
    keptLdBlockInfoVec=makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks=keptLdBlockInfoVec.size();
    if(!numKeptLDBlocks) throw std::runtime_error("No SNP maps to the provided LD blocks");
    LOGGER << numKeptLDBlocks << " GWAS LD blocks retained." << endl;
}

void Data::makeBlockLDmatrix(const string &bedFile, const string &LDmatType, const unsigned block, const string &dirname, const bool writeLdmTxt, int ldBlockRegionWind,const bool isCovBool){
    string ldNameType = "LD",ldSuffixType = "ldm";
    if(isCovBool) {
        ldNameType = "genotype covariance";
        ldSuffixType = "covm";
    }
    LOGGER << "Making block "<< ldNameType << " matricies ..." << endl;
    
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        string create_cmd = "mkdir -p " + dirname;
        system(create_cmd.c_str());
        LOGGER << "Created folder [" << dirname << "] to store LD matrices." << endl;
    }
    
    mapSnpsToBlocks();
    
    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();

#pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < numKeptLDBlocks; i++) {
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        if (block && std::stoul(ldblock->ID) != block) {
            continue;
        }


        string outNameType;
        // if(isCovBool) {
            //  outNameType = "blockchr" + to_string(ldblock->chrom) + "_" + to_string(ldblock->startPos) + "_" + to_string(ldblock->endPos);
        // } else {
        outNameType = "block" + ldblock->ID;
        // }

        string outBinfile = dirname + "/" + outNameType + "." + ldSuffixType + ".bin";
        string outSnpfile = dirname + "/" + outNameType + ".snp.info";
        string outldmfile = dirname + "/" + outNameType + "." + ldSuffixType + ".info";

        FILE *outbin = fopen(outBinfile.c_str(), "wb");
        if (!outbin) {
            LOGGER.e(0, "Cannot open file [" + outBinfile + "] to write.");
        }
        ofstream outtxt;
        string outTxtfile;
        if (writeLdmTxt) {
            outTxtfile = dirname + "/" + outNameType + "." + ldSuffixType + ".txt";
            outtxt.open(outTxtfile.c_str());
        }

        if(!i) LOGGER << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
            
        MatrixXf rval = generateLDmatrixPerBlock(bedFile, ldblock->snpNameVec,isCovBool).cast<float>();
        
        unsigned numSnpInBlock = ldblock->numSnpInBlock;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numSnpInBlock;
        fwrite(rval.data(), sizeof(float), nElements, outbin);
        
        if (writeLdmTxt) {
            for (unsigned ii=0; ii<numSnpInBlock; ++ii){
                for (unsigned jj=0; jj<numSnpInBlock; ++jj) {
                    outtxt << ldblock->ID << "\t" << ldblock->snpNameVec[ii] << "\t" << ldblock->snpNameVec[jj] << "\t" << rval(ii,jj) << endl;
                }
            }
        }
        
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();
        
        outputBlockLDmatrixInfo(*ldblock, outSnpfile, outldmfile);
        
        if(!(i%1)) LOGGER << " computed block " << ldblock->ID << "\r" << flush;
        if (block) {
            LOGGER << "Written the " << ldNameType << " matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) LOGGER << "Written the " << ldNameType << " matrix into file [" << outTxtfile << "]." << endl;
            LOGGER << "Written the " << ldNameType << " matrix SNP info into file [" << outSnpfile << "]." << endl;
            LOGGER << "Written the " << ldNameType << " matrix ldm info into file [" << outldmfile << "]." << endl;
        }
    }
    
    if (!block) {
        LOGGER << "Written the " << ldNameType << " matrix into folder [" << dirname << "/block*." << ldSuffixType << ".bin]." << endl;
        if (writeLdmTxt) LOGGER << "Written the " << ldNameType << " matrix into text file [" << dirname << "/block*." << ldSuffixType << ".txt]." << endl;
        
        // if (chromInfoVec.size() >= 22) {  // genome-wide build of LD matrices
        //     mergeLdmInfo(LDmatType, dirname);
        // } else {
            LOGGER << "Written the " << ldNameType << " matrix into folder [" << dirname << "/block*.snp.info]." << endl;
            LOGGER << "Written the " << ldNameType << " matrix into folder [" << dirname << "/block*." << ldSuffixType << ".info]." << endl;
        // }
    }
}


void Data::outputBlockLDmatrixInfo(const LDBlockInfo &block, const string &outSnpfile, const string &outldmfile) const {
    // write snp info
    ofstream out1(outSnpfile.c_str());
    out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12s %10s %10s\n")
    % "Chrom"
    % "ID"
    % "Index"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "N"
    % "Block";
    SnpInfo *snp = NULL;
    for (unsigned i=0; i < block.numSnpInBlock; ++i) {
        snp = block.snpInfoVec[i];
        out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %22.17g %10s %10s\n")
        % snp->chrom
        % snp->rsID
        % snp->index
        % snp->genPos
        % snp->physPos
        % snp->a1
        % snp->a2
        % snp->af
        % numKeptInds
        % snp->block;
    }
    out1.close();

    // svd matrix for ld blocks here.
    ofstream out2(outldmfile.c_str());
    out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
    % "Block"
    % "Chrom"
    % "StartSnpIdx"
    % "StartSnpID"
    % "EndSnpIdx"
    % "EndSnpID"
    % "NumSnps";
    out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
    % block.ID
    % block.chrom
    % block.startSnpIdx
    % incdSnpInfoVec[block.startSnpIdx]->rsID
    % block.endSnpIdx
    % incdSnpInfoVec[block.endSnpIdx]->rsID
    % block.numSnpInBlock;
    out2.close();
}


void Data::impG(const unsigned block, double diag_mod){
    VectorXi numImpSnp;
    VectorXi numTypSnp;
    numImpSnp.setZero(numLDBlocks);
    numTypSnp.setZero(numLDBlocks);
    for (unsigned i = 0; i < numLDBlocks; i++ ){
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (!ldblock->kept) continue;
        for (unsigned j=0; j<ldblock->numSnpInBlock; ++j) {
            SnpInfo *snp = ldblock->snpInfoVec[j];
            if (snp->included) {
                ++numTypSnp[i];
            } else {
                ++numImpSnp[i];
            }
        }
    }
    unsigned totalNumImpSnp = numImpSnp.sum();
        
    LOGGER << "Imputing summary statistics for " << to_string(totalNumImpSnp) << " SNPs in the LD reference but not in the GWAS data file..." << endl;

    Gadget::Timer timer;
    timer.setTime();
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < numLDBlocks; i++ ){
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (!ldblock->kept) continue;
        
        if (numImpSnp[i]) {
            
            Stat::Normal normal;
            
            /// Step 1. construct LD
            MatrixXd LDPerBlock = eigenVecLdBlock[i] * eigenValLdBlock[i].asDiagonal() * eigenVecLdBlock[i].transpose();
            
            LDPerBlock.diagonal().array() += (double)diag_mod;
            /// Step 2. Construct the LD correlation matrix among the typed SNPs(LDtt) and the LD correlation matrix among the missing SNPs and typed SNPs (LDit).
            // Step 2.1 divide SNPs into typed and untyped SNPs
            VectorXi typedSnpIdx(numTypSnp[i]);
            VectorXi untypedSnpIdx(numImpSnp[i]);
            VectorXd zTypSnp(numTypSnp[i]);
            VectorXd nTypSnp(numTypSnp[i]);
            VectorXd varyTypSnp(numTypSnp[i]);
            for(unsigned j=0, idxTyp=0, idxImp=0; j < ldblock->numSnpInBlock; j++){
                SnpInfo *snp = ldblock->snpInfoVec[j];
                if(snp->included){
                    // typed snp
                    typedSnpIdx[idxTyp] = j;
                    zTypSnp[idxTyp] = snp->gwas_b / snp->gwas_se;
                    nTypSnp[idxTyp] = snp->gwas_n;
                    double hetj = 2.0 * snp->gwas_af * (1.0 - snp->gwas_af);
                    varyTypSnp[idxTyp] = hetj * (snp->gwas_n * snp->gwas_se * snp->gwas_se + snp->gwas_b * snp->gwas_b);
                    ++idxTyp;
                } else {
                    untypedSnpIdx[idxImp] = j;
                    ++idxImp;
                }
            }
            // Step 2.2 construct LDtt and LDit and Ztt.
            MatrixXd LDtt = LDPerBlock(typedSnpIdx,typedSnpIdx);
            MatrixXd LDit = LDPerBlock(untypedSnpIdx,typedSnpIdx);
            // Step 2.3 //  The Z score for the missing SNPs;
            VectorXd LDi_Z = LDtt.ldlt().solve(zTypSnp);
            VectorXd zImpSnp = LDit * LDi_Z;
            // Step 3. re-calcualte beta and se
            // if snp is missing use median to replace N
            std::sort(nTypSnp.data(), nTypSnp.data() + nTypSnp.size());
            double nMedian = nTypSnp[nTypSnp.size()/2];  // median
            // calcuate median of phenotypic variance
            std::sort(varyTypSnp.data(), varyTypSnp.data() + varyTypSnp.size());
            double varyMedian = varyTypSnp[varyTypSnp.size()/2];  // median
            // begin impute
            for(unsigned j = 0; j < numImpSnp[i]; j++){
                SnpInfo *snp = ldblock->snpInfoVec[untypedSnpIdx[j]];
                double base = sqrt(2.0 * snp->af *(1.0 - snp->af) * (nMedian + zImpSnp[j] * zImpSnp[j]));
                snp->gwas_b = zImpSnp[j] * sqrt(varyMedian)/base;
                snp->gwas_se = sqrt(varyMedian) / base;
                snp->gwas_n = nMedian;
                snp->gwas_af = snp->af;
                snp->gwas_pvalue = 2*(1.0-normal.cdf_01(abs(snp->gwas_b/snp->gwas_se)));
                snp->included = true;
                //LOGGER << "b " << snp->gwas_b << " se " << snp->gwas_se << " z " << snp->gwas_b/snp->gwas_se << " p " << snp->gwas_pvalue << endl;
            }
        }
        
        if(!(i%10)) LOGGER << " imputed block " << i << "\r" << flush;
        
        if (block) {
            string outfile = title + ".block" + ldblock->ID + ".imputed.ma";
            ofstream out(outfile.c_str());
            out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n") % "SNP" % "A1" % "A2" % "freq" % "b" % "se" % "p" % "N";
            for (unsigned i=0; i<ldblock->numSnpInBlock; ++i) {
                SnpInfo *snp = ldblock->snpInfoVec[i];
                out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n")
                % snp->rsID
                % snp->a1
                % snp->a2
                % snp->gwas_af
                % snp->gwas_b
                % snp->gwas_se
                % snp->gwas_pvalue
                % snp->gwas_n;
            }
            out.close();

            timer.getTime();
            LOGGER << "Imputation of summary statistics is completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
            LOGGER << "Summary statistics of all SNPs are save into file [" + outfile + "]." << endl;

        }

    }

    if (block) return;
    
    string outfile = title + ".imputed.ma";
    ofstream out(outfile.c_str());
    out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n") % "SNP" % "A1" % "A2" % "freq" % "b" % "se" % "p" % "N";
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n")
        % snp->rsID
        % snp->a1
        % snp->a2
        % snp->gwas_af
        % snp->gwas_b
        % snp->gwas_se
        % snp->gwas_pvalue
        % snp->gwas_n;
    }
    out.close();

    timer.getTime();
    LOGGER << "Imputation of summary statistics is completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
    LOGGER << "Summary statistics of all SNPs are save into file [" + outfile + "]." << endl;
}

void Data::getEigenDataForLDBlock(const string &bedFile, const string &ldBlockInfoFile, int ldBlockRegionWind, const string &filename, const float eigenCutoff){
    int i,j;
    vector<locus_bp> snpVec;
    SnpInfo *snp = NULL;

    map<int, string>  chrEndSnp;
    for (i = 1; i < numIncdSnps; i++) {
        snp = incdSnpInfoVec[i];
        if(incdSnpInfoVec[i]->chrom != incdSnpInfoVec[i-1]->chrom){
            chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[i - 1]->chrom,incdSnpInfoVec[i - 1]->rsID ));
        }
    }
    chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[numIncdSnps - 1]->chrom,incdSnpInfoVec[numIncdSnps - 1]->rsID ));
    /////////////////////////////////////////
    // Step 2. Map snps to blocks
    /////////////////////////////////////////
    vector<string> block2snp_1(numLDBlocks), block2snp_2(numLDBlocks);
    map<string,int> keptLdBlock2AllLdBlcokMap;
    vector<locus_bp>::iterator iter;
    map<int, string>::iterator chrIter;
    for (i = 0; i < numIncdSnps ; i++) {
        snp = incdSnpInfoVec[i];
        snpVec.push_back(locus_bp(snp->rsID, snp->chrom, snp->physPos ));
    }
#pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numLDBlocks; i++) {
        // find lowest snp_name in the block
        LDBlockInfo *ldblock = ldBlockInfoVec[i];

        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp( ldblock->ID ,ldblock->chrom, ldblock->startPos - ldBlockRegionWind));
        if (iter != snpVec.end()) block2snp_1[i] = iter->locusName;
        else block2snp_1[i] = "NA";
    }
#pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numLDBlocks; i++) {
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (block2snp_1[i] == "NA") {
            block2snp_2[i] = "NA";
            continue;
        }
        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp(ldblock->ID, ldblock->chrom, ldblock->endPos + ldBlockRegionWind));
        if (iter != snpVec.end()){
            if (iter->bp ==  ldblock->endPos + ldBlockRegionWind){
                block2snp_2[i] = iter->locusName;
            }else {
                if(iter!=snpVec.begin()){
                    iter--;
                    block2snp_2[i] = iter->locusName;
                }
                else block2snp_2[i] = "NA";
            }
        }
        else {
            chrIter = chrEndSnp.find(ldblock->chrom);
            if (chrIter == chrEndSnp.end()) block2snp_2[i] = "NA";
            else block2snp_2[i] = chrIter->second;
        }
    }
    int mapped = 0;
    for (i = 0; i < numLDBlocks; i++) {
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (block2snp_1[i] != "NA" && block2snp_2[i] != "NA")
        {
            mapped++;
            // ldblock->kept = true;
            keptLdBlock2AllLdBlcokMap.insert(pair<string, int>(ldblock->ID, i));
        } else {
            ldblock->kept = false;
        }
    }
    if (mapped < 1) LOGGER.e(0, "No SNP can be mapped to the provided LD block list. Please check the input data regarding chromosome and bp.");
    else cout << mapped << " LD blocks have at least one SNP." << endl;

    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();

    map<string, int>::iterator iter1, iter2;
    map<string, int> snpNameMap;
    VectorXf snpNumInldblock(numKeptLDBlocks);
    VectorXf eigenNumInldblock(numKeptLDBlocks);
    vector<VectorXf> cumsumNonNeg(numKeptLDBlocks);
    
    for (i = 0; i < numIncdSnps; i++) {
        SnpInfo *snp = incdSnpInfoVec[i];
        snpNameMap.insert(pair<string,int>(snp->rsID, i));
        
    }
      // eigenvector and eigenvalue
    //eigenValLdBlock.resize(numKeptLDBlocks);
    //eigenVecLdBlock.resize(numKeptLDBlocks);
    string lastSNPIDInPreLDblocks;
    string LDmatType = "block";
    string outfilename = filename + "." + LDmatType + ".eigen" + ".bin";
    FILE *out3 = fopen(outfilename.c_str(), "wb");

    bool readBedBool = true;
    for (i = 0; i < numKeptLDBlocks; i++) {
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        // cout << "ldblock id: " << ldblock->ID << endl;
        iter1 = snpNameMap.find(block2snp_1[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
        iter2 = snpNameMap.find(block2snp_2[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
        bool skip = false;
        if (iter1 == snpNameMap.end() || iter2 == snpNameMap.end() || iter1->second >= iter2->second) ldblock->kept = false;
        snpNumInldblock[i] = iter2->second - iter1->second + 1;
        // cout << "ldblock->kept: " << ldblock->kept << endl;
        if(!ldblock->kept) continue;
        vector<int> snp_indx;
        for (j = iter1->second; j <= iter2->second; j++) {
            if(snp->rsID == lastSNPIDInPreLDblocks){continue;} // if snp is matched to two ldblock
            snp_indx.push_back(j);
            ldblock->snpNameVec.push_back(incdSnpInfoVec[j]->rsID);
            if (j == iter2->second) { lastSNPIDInPreLDblocks = snp->rsID;}
        }
        if(readBedBool) {
            cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
            readBedBool = false;
        }
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal = 0;
        MatrixXf rval = generateLDmatrixPerBlock(bedFile, ldblock->snpNameVec).cast<float>();
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(rval, eigenCutoff,eigenVal, eigenVec, sumPosEigVal);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = ldblock->snpNameVec.size();
        eigenNumInldblock[i] = numEigenValue;
        // save summary
        // 1, nrow of eigenVecGene[i]
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out3);
        //        cout << "rval: " << rval << endl;
        // cout << "eigenVec: "  << eigenVec << endl;
        // cout << "";
        // 2, ncol of eigenVecGene[i]
        fwrite(&numEigenValue, sizeof(int32_t), 1, out3);
        // 3. sum of all positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, out3);
        // 4. eigenCutoff
        fwrite(&eigenCutoff, sizeof(float), 1, out3);
        // 5. eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out3);
        // 6. eigenvectors;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, out3);
        cout << " Generate and save Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << "\r" << flush;
    }
    fclose(out3);
    cout << "To explain " << eigenCutoff*100 << "% variance in LD, on average " << int(eigenNumInldblock.mean()) << " eigenvalues are selected across LD blocks (mean number of SNPs is " << int(snpNumInldblock.mean()) << ")." << endl;

}

void Data::readBlockLdmBinaryAndDoEigenDecomposition(const string &dirname, const unsigned block, const float eigenCutoff, const bool writeLdmTxt){
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        LOGGER.e(0,"Error: cannot find the folder [" + dirname + "]");
    }
    
    vector<int> numSnpInRegion;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        LDBlockInfo *block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
        
    keptLdBlockInfoVec = ldBlockInfoVec;
    
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        if (block && ldBlockInfoVec[i]->ID != std::to_string(block)) continue;
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[i];

        string infile = dirname + "/block" + blockInfo->ID + ".ldm.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){LOGGER.e(0,"Error: can not open the file [" + infile + "] to read.");}

        string outBinfile = dirname + "/block" + blockInfo->ID + ".eigen.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");
        if (!outbin) {
            LOGGER.e(0, "Cannot open file [" + outBinfile + "] to write.");
        }

        string outTxtfile;
        ofstream outtxt;
        if (writeLdmTxt) {
            outTxtfile = dirname + "/block" + blockInfo->ID + ".eigen.txt";
            outtxt.open(outTxtfile.c_str());
        }
        
        int32_t blockSize = blockInfo->numSnpInBlock;
        
        MatrixXf ldm(blockSize, blockSize);
        uint64_t nElements = (uint64_t)blockSize * (uint64_t)blockSize;
                
        if(fread(ldm.data(), sizeof(float), nElements, fp) != nElements){
            // cout << "fread(U.data(), sizeof(float), nElements, fp): " << fread(ldm.data(), sizeof(float), nElements, fp) << endl;
            cout << "nEle: " << nElements << " ldm.size: " << ldm.size() <<  " ldm.col: " << ldm.cols() << " row: " << ldm.rows() << endl;
            LOGGER.e(0,"In LD block " + blockInfo->ID + ",size error in " + outBinfile);
            // cout << "Read " << svdLDfile << " error (U)" << endl;
            // LOGGER.e(0,"read file error");
        }
        
        fclose(fp);
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal;  // sum of all positive eigenvalues
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(ldm, eigenCutoff, eigenVal, eigenVec, sumPosEigVal);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix

        blockInfo->sumPosEigVal = sumPosEigVal;
        
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = blockSize;

        //cout << " Generate Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << endl;
        
        // save summary
        // 1, the number of SNPs in the block
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, outbin);
        // 2, the number of eigenvalues at with the given cutoff
        fwrite(&numEigenValue, sizeof(int32_t), 1, outbin);
        // 3. sum of all the positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, outbin);
        // 4. eigenvalue cutoff based on the proportion of variance explained in LD
        fwrite(&eigenCutoff, sizeof(float), 1, outbin);
        // 5. the selected eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, outbin);
        // 6. the selected eigenvector;
        nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, outbin);
        
        if (writeLdmTxt) {
            outtxt << "Block " << blockInfo->ID << endl;
            outtxt << "numSnps " << numSnpInBlock << endl;
            outtxt << "numEigenvalues " << numEigenValue << endl;
            outtxt << "SumPositiveEigenvalues " << sumPosEigVal << endl;
            outtxt << "EigenCutoff " << eigenCutoff << endl;
            outtxt << "Eigenvalues\n" << eigenVal.transpose() << endl;
            outtxt << "Eigenvectors\n" << eigenVec << endl;
            outtxt << endl;
        }
        
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();

        if(!(i%10)) cout << " computed block " << blockInfo->ID << "\r" << flush;
        
        if (block) {
            cout << "Written the eigen data for block LD matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) cout << "Written the eigen data for block LD matrix into file [" << outTxtfile << "]." << endl;
        }
    }

    if (!block) {
        cout << "Written the eigen data for block LD matrix into file [" << dirname << "/block*.eigen.bin]." << endl;
        if (writeLdmTxt) cout << "Written the eigen data for block LD matrix into file [" << dirname << "/block*.eigen.txt]." << endl;
    }

}


void Data::readEigenMatrixBinaryFile(const string &dirname, const double eigenCutoff){
    if (!Gadget::directoryExist(dirname)) {
        LOGGER.e(0," cannot find the folder [" + dirname + "]");
    }
    
    vector<int>numSnpInRegion;
    LDBlockInfo * block;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
    eigenValLdBlock.resize(numLDBlocks);
    eigenVecLdBlock.resize(numLDBlocks);
        
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        LDBlockInfo * block;
        block = ldBlockInfoVec[i];
        
        if (!block->kept) continue;
        
        int32_t cur_m = 0;
        int32_t cur_k = 0;
        float sumPosEigVal = 0;
        float oldEigenCutoff =0;
        
        string infile = dirname + "/block" + block->ID + ".eigen.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){LOGGER.e(0,"can not open the file [" + infile + "] to read.");}

        // 1. marker number
        if(fread(&cur_m, sizeof(int32_t), 1, fp) != 1){
            LOGGER.e(0,"Read " + infile + " error (m)");
        }    
        if(cur_m != numSnpInRegion[i]){
            LOGGER.e(0,"In LD block " + block->ID + ", inconsistent marker number to marker information in " + infile);
        }
        // 2. ncol of eigenVec (number of eigenvalues)
        if(fread(&cur_k, sizeof(int32_t), 1, fp) != 1){
            LOGGER.e(0,"In LD block " + block->ID + ", error about number of eigenvalues in  " + infile);
        }
        // 3. sum of all positive eigenvalues
        if(fread(&sumPosEigVal, sizeof(float), 1, fp) != 1){
            LOGGER.e(0,"In LD block " + block->ID + ", error about the sum of positive eigenvalues in " + infile);
        }
        // 4. eigenCutoff
        if(fread(&oldEigenCutoff, sizeof(float), 1, fp) != 1){
            LOGGER.e(0,"In LD block " + block->ID + ", error about eigen cutoff used in " + infile);
        }
        // 5. eigenvalues
        if(cur_m<=0 || cur_k<=0 || cur_k>cur_m)throw std::runtime_error("Invalid LD eigen dimensions");
        MemoryBudget::require(MemoryBudget::bytes(cur_m,cur_k,64),"LD eigen block, transforms and copies");
        VectorXf lambda(cur_k);
        if(fread(lambda.data(), sizeof(float), cur_k, fp) != cur_k){
            LOGGER.e(0,"In LD block " + block->ID + ",size error about eigenvalues in " + infile);
        }
        // 6. eigenvector
        MatrixXf U(cur_m, cur_k);
        uint64_t nElements = (uint64_t)cur_m * (uint64_t)cur_k;
        if(fread(U.data(), sizeof(float), nElements, fp) != nElements){
            // LOGGER << "fread(U.data(), sizeof(double), nElements, fp): " << fread(U.data(), sizeof(float), nElements, fp) << endl;
            LOGGER << "nEle: " << nElements << " U.size: " << U.size() <<  " U.col: " << U.cols() << " row: " << U.rows() << endl;
            LOGGER.e(0,"In LD block " + block->ID + ",size error about eigenvectors in " + infile);
        }
        fclose(fp);
        bool haveValue = false;
        int revIdx = 0;
        if(oldEigenCutoff < eigenCutoff && i == 0){
            LOGGER << "Warning: current proportion of variance in LD block is set as " + to_string(eigenCutoff)+ ". But the proportion of variance is set as "<< to_string(oldEigenCutoff) + " in "  + infile + ".\n";
        }
        if (eigenCutoff < oldEigenCutoff) {
            truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda.cast<double>(), U.cast<double>(), eigenValLdBlock[i], eigenVecLdBlock[i]);
        } else {
            eigenValLdBlock[i] = lambda.cast<double>();
            eigenVecLdBlock[i] = U.cast<double>();
        }
        block->sumPosEigVal = sumPosEigVal;
        block->eigenvalues = lambda.cast<double>();
    }
    LOGGER << "GWAS LDM of" << numLDBlocks << " GWAS LD blocks to be include from [" << dirname << "]." << endl;
}

VectorXd Data::scoreOnLDScale(const VectorXd &score, const vector<string> &ids) const {
    if (!ldCorrelation || suppliedGenotypeVariance.empty()) return score;
    if (score.size()!=ids.size()) throw std::runtime_error("LD score/ID dimension mismatch");
    VectorXd result=score;
    for (size_t j=0;j<ids.size();++j) {
        const double variance=suppliedGenotypeVariance.at(ids[j]);
        if (!(variance>0) || !std::isfinite(variance)) throw std::runtime_error("Invalid LD genotype variance");
        result[j]/=std::sqrt(variance);
    }
    return result;
}

void Data::designOnGenotypeScale(MatrixXd &design, const vector<string> &ids) const {
    if (!ldCorrelation || suppliedGenotypeVariance.empty()) return;
    if (design.cols()!=ids.size()) throw std::runtime_error("LD design/ID dimension mismatch");
    for (size_t j=0;j<ids.size();++j) {
        const double variance=suppliedGenotypeVariance.at(ids[j]);
        if (!(variance>0) || !std::isfinite(variance)) throw std::runtime_error("Invalid LD genotype variance");
        design.col(j)*=std::sqrt(variance);
    }
}

void Data::readEigenMatrixBinaryFileAndMakeWandQ(const string &dirname, const double eigenCutoff, const vector<VectorXd> &GWASeffects, const double nGWAS, const bool noscale, const bool makePseudoSummary){
    if (!Gadget::directoryExist(dirname)) {
        LOGGER.e(0, " cannot find the folder [" + dirname + "]");
    }
    
    vector<int>numSnpInRegion(numKeptLDBlocks);
    
    for(int i = 0; i < numKeptLDBlocks;i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
    eigenValLdBlock.resize(numLDBlocks);
    eigenVecLdBlock.resize(numLDBlocks);
    wcorrBlocks.resize(numKeptLDBlocks);
    numSnpsBlock.resize(numKeptLDBlocks);
    numEigenvalBlock.resize(numKeptLDBlocks);
    Qblocks.resize(numKeptLDBlocks);

    
    //Constructing pseudo summary statistics for training and validation data sets, with 90% sample size for training and 10% for validation
    
    float n_trn, n_val;
    if (makePseudoSummary) {
        n_trn = 0.9*float(numKeptInds);
        n_val = numKeptInds - n_trn;
        pseudoGwasNtrn = n_trn;

        pseudoGwasEffectTrn.resize(numKeptLDBlocks);
        pseudoGwasEffectVal.resize(numKeptLDBlocks);
        b_val.setZero(numIncdSnps);
    }

// #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        int32_t cur_m = 0;
        int32_t cur_k = 0;
        float sumPosEigVal = 0;
        float oldEigenCutoff =0;
        
        string infile = dirname + "/block" + block->ID + ".eigen.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){LOGGER.e(0, " can not open the file [" + infile + "] to read.");}

        // 1. marker number
        if(fread(&cur_m, sizeof(int32_t), 1, fp) != 1){
            LOGGER.e(0, "Read " + infile + " error (m)");
        }
                
        if(cur_m != numSnpInRegion[i]){
            LOGGER.e(0, "In LD block " + block->ID + ", inconsistent marker number to marker information in " + infile);
        }
        // 2. ncol of eigenVec (number of eigenvalues)
        if(fread(&cur_k, sizeof(int32_t), 1, fp) != 1){
            LOGGER.e(0, "In LD block " + block->ID + ", error about number of eigenvalues in  " + infile);
        }
        // 3. sum of all positive eigenvalues
        if(fread(&sumPosEigVal, sizeof(float), 1, fp) != 1){
            LOGGER.e(0, "In LD block " + block->ID + ", error about the sum of positive eigenvalues in " + infile);
        }
        // 4. eigenCutoff
        if(fread(&oldEigenCutoff, sizeof(float), 1, fp) != 1){
            LOGGER.e(0, "In LD block " + block->ID + ", error about eigen cutoff used in " + infile);
        }
        // 5. eigenvalues
        if(cur_m<=0 || cur_k<=0 || cur_k>cur_m)throw std::runtime_error("Invalid LD eigen dimensions");
        MemoryBudget::require(MemoryBudget::bytes(cur_m,cur_k,64),"LD eigen block, transforms and copies");
        VectorXf lambda(cur_k);
        if(fread(lambda.data(), sizeof(float), cur_k, fp) != cur_k){
            LOGGER.e(0, "In LD block " + block->ID + ",size error about eigenvalues in " + infile);
        }
        // 6. eigenvector
        MatrixXf U(cur_m, cur_k);
        uint64_t nElements = (uint64_t)cur_m * (uint64_t)cur_k;
        if(fread(U.data(), sizeof(float), nElements, fp) != nElements){
            // cout << "fread(U.data(), sizeof(float), nElements, fp): " << fread(U.data(), sizeof(float), nElements, fp) << endl;
            cout << "nEle: " << nElements << " U.size: " << U.size() <<  " U.col: " << U.cols() << " row: " << U.rows() << endl;
            LOGGER.e(0, "In LD block " + block->ID + ",size error about eigenvectors in " + infile);
        }
        bool haveValue = false;
        int revIdx = 0;
        if(oldEigenCutoff < eigenCutoff & i == 0){
            cout << "Warning: current proportion of variance in LD block is set as " + to_string(eigenCutoff)+ ". But the proportion of variance is set as "<< to_string(oldEigenCutoff) + " in "  + infile + ".\n";
            // LOGGER.e(0, "");
        }
        if (eigenCutoff < oldEigenCutoff) {
            truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda.cast<double>(), U.cast<double>(), eigenValLdBlock[i], eigenVecLdBlock[i]);
        } else {
            eigenValLdBlock[i] = lambda.cast<double>();
            eigenVecLdBlock[i] = U.cast<double>();
        }
        block->sumPosEigVal = sumPosEigVal;
        block->eigenvalues = lambda.cast<double>();
        
        // make w and Q
        VectorXd sqrtLambda = eigenValLdBlock[i].array().sqrt();
        wcorrBlocks[i] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * scoreOnLDScale(GWASeffects[i],block->snpNameVec));
        //MatrixXf tmpQblocks = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        //MatrixDat matrixDat = MatrixDat(block->snpNameVec, tmpQblocks);
        Qblocks[i] = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        designOnGenotypeScale(Qblocks[i],block->snpNameVec);
        
        if (noscale) {
            VectorXd Dsqrt(block->numSnpInBlock);
            for (unsigned j=0; j<block->numSnpInBlock; ++j) {
                SnpInfo *snp = block->snpInfoVec[j];
                Dsqrt[j] = sqrt(snp->twopq);
            }
            Qblocks[i] = Qblocks[i] * Dsqrt.asDiagonal();
        }

        numSnpsBlock[i] = Qblocks[i].cols();
        numEigenvalBlock[i] = Qblocks[i].rows();
        
        // make pseudo summary data
        if (makePseudoSummary) {
            long size = eigenValLdBlock[i].size();
            VectorXd rnd(size);
            for (unsigned j=0; j<size; ++j) {
                rnd[j] = Stat::snorm();
            }
            
            pseudoGwasEffectTrn[i] = gwasEffectInBlock[i] + sqrt(1.0/n_trn - 1.0/nGWASblock[i]) * eigenVecLdBlock[i] * (eigenValLdBlock[i].array().sqrt().matrix().asDiagonal() * rnd);

            pseudoGwasEffectVal[i] = nGWASblock[i]/n_val * gwasEffectInBlock[i] - n_trn/n_val * pseudoGwasEffectTrn[i];
            b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
        }
        
        // eigenVecLdBlock[i].resize(0,0);
    }

    nGWASblock.resize(numKeptLDBlocks);
    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        nGWASblock[i] = nGWAS;
    }
}

void Data::UseLDBlockEigenMakeWAndQblocks(const vector<VectorXd> &GWASeffects, const double nGWAS, const bool noscale,const bool makePseudoSummary){
    ///// 
    float n_trn, n_val;
    if (makePseudoSummary) {
        n_trn = 0.9*float(numKeptInds);
        n_val = numKeptInds - n_trn;
        pseudoGwasNtrn = n_trn;
        pseudoGwasEffectTrn.resize(numKeptLDBlocks);
        pseudoGwasEffectVal.resize(numKeptLDBlocks);
        b_val.setZero(numIncdSnps);
    }
    wcorrBlocks.resize(numKeptLDBlocks);
    numSnpsBlock.resize(numKeptLDBlocks);
    numEigenvalBlock.resize(numKeptLDBlocks);
    Qblocks.resize(numKeptLDBlocks);
    VectorXd sqrtLambda;

    for(int i = 0,k = 0; i < numLDBlocks; i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        if(!block->kept) continue;
        sqrtLambda = eigenValLdBlock[i].array().sqrt();
        wcorrBlocks[k] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * scoreOnLDScale(GWASeffects[i],block->snpNameVec));
    
        Qblocks[k] = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        designOnGenotypeScale(Qblocks[k],block->snpNameVec);
        if (noscale) {
            VectorXd Dsqrt(block->numSnpInBlock);
            for (unsigned j=0; j<block->numSnpInBlock; ++j) {
                SnpInfo *snp = block->snpInfoVec[j];
                Dsqrt[j] = sqrt(snp->twopq);
            }
            Qblocks[k] = Qblocks[i] * Dsqrt.asDiagonal();
        }

        numSnpsBlock[k] = Qblocks[i].cols();
        numEigenvalBlock[k] = Qblocks[i].rows();

        // make pseudo summary data
        if (makePseudoSummary) {
            long size = eigenValLdBlock[i].size();
            VectorXd rnd(size);
            for (unsigned j=0; j<size; ++j) {
                rnd[j] = Stat::snorm();
            }
            pseudoGwasEffectTrn[k] = gwasEffectInBlock[i] + sqrt(1.0/n_trn - 1.0/nGWASblock[i]) * eigenVecLdBlock[i] * (sqrtLambda.matrix().asDiagonal() * rnd);
            pseudoGwasEffectVal[k] = nGWASblock[i]/n_val * gwasEffectInBlock[i] - n_trn/n_val * pseudoGwasEffectTrn[k];
            b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
        }
        k++;
        if(k == numKeptLDBlocks) break;
    }

    nGWASblock.resize(numKeptLDBlocks);
    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        nGWASblock[i] = nGWAS;
    }
}

void Data::truncateEigenMatrix(const float total, const double cutoff, const VectorXd &values, const MatrixXd &vectors, VectorXd &outValues, MatrixXd &outVectors) {
    if(values.size()==0 || vectors.cols()!=values.size() || !values.allFinite() || !vectors.allFinite() || !(total>0) || !(cutoff>0 && cutoff<=1))
        throw std::runtime_error("Invalid stored eigen decomposition or cutoff");
    std::vector<int> order(values.size());std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](int a,int b){return values[a]>values[b];});
    std::vector<int> kept;double sum=0;
    for(int j:order) {
        if(values[j]<=1e-10) continue;
        if(!kept.empty() && sum+2*std::numeric_limits<float>::epsilon()*total>=cutoff*total) break;
        kept.push_back(j);sum+=values[j];
    }
    if(kept.empty()) throw std::runtime_error("No positive LD eigenvalues");
    // Preserve the file's ordering after choosing the largest eigenvalues.
    std::sort(kept.begin(),kept.end());outValues.resize(kept.size());outVectors.resize(vectors.rows(),kept.size());
    for(size_t j=0;j<kept.size();++j){outValues[j]=values[kept[j]];outVectors.col(j)=vectors.col(kept[j]);}
}

void Data::readEigenMatrix(const string &dirname, const double eigenCutoff,const bool qcBool){
    LOGGER  << "............................" << endl;
    LOGGER << "Reading GWAS low-rank LD matrices..." << endl;
    LOGGER  << "............................" << endl;
    if(boost::filesystem::is_directory(dirname)) {
        LOGGER << "Reading GWAS low-rank LD matrices (SBayesRC format) ..." << endl;
        readSBayesRCBlockLdmInfoFile(dirname + "/ldm.info");
        readSBayesRCBlockLdmSnpInfoFile(dirname + "/snp.info");
        if(qcBool) readEigenMatrixBinaryFile(dirname,eigenCutoff);
    } else {
        string filename = dirname + ".eigen.ldblock";
        string LDMatType = "";
        readEigenMatLDBlockSnpInfoFile( filename + ".snp.info" ,LDMatType);
        readEigenMatLDBlockInfoFile(filename + ".info");
        readEigenMatBinFile(dirname,eigenCutoff,"ldblock");
    }
}

vector<LDBlockInfo*> Data::makeKeptLDBlockInfoVec(const vector<LDBlockInfo*> &ldBlockInfoVec){
    vector<LDBlockInfo*> keptLDBlock;
    ldblockNames.clear();
    LDBlockInfo * ldblock = NULL;
    for (unsigned i=0, j=0; i< numLDBlocks; ++i) {
        ldblock = ldBlockInfoVec[i];
        if(ldblock->kept) {
            ldblock->index = j++;  // reindex inds
            keptLDBlock.push_back(ldblock);
            ldblockNames.push_back(ldblock->ID);
        }
    }
    return keptLDBlock;
}

void Data::buildMMEigen(const string &dirname, const bool sampleOverlap, const double eigenCutoff, const bool noscale,const bool gwasInfoOnly){
    // check all available LD Blocks.
    includeMatchedBlocks();
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        if (snp->gwas_b == -999) {
            LOGGER.e(0," SNP " + snp->rsID + " in the LD reference has no summary data. Run --impute-summary first.");
        }
    }
    LOGGER  << "............................" << endl;
    LOGGER << "Building GWAS various input parameters based on " << numKeptLDBlocks << " LD Blocks with " << numIncdSnps << " SNPs ... " << endl;
    LOGGER  << "............................" << endl;
    if(gwasInfoOnly) {
        LOGGER << "Construct various maps to indicate the relationship between the complex trait and GWAS." << endl;
        ConstructGwasEqtlGeneMaps();
    }
    LOGGER << "Scale GWAS effect size assuming unit GWAS trait phenotypic variance (Jian Yang et.al. (2012))..." << endl;
    scaleGwasEffects();

    LOGGER << "Construct corrected GWAS trait phenotype vector based on GWAS low-rank LD matrices (Lloyd-Jones, Zeng et.al. (2019))..." << endl;

    // read eigen LD matrix here 
    if(boost::filesystem::is_directory(dirname)) {
        readEigenMatrixBinaryFileAndMakeWandQ(dirname, eigenCutoff, gwasEffectInBlock, numKeptInds, noscale, true);
     } else {
        UseLDBlockEigenMakeWAndQblocks(gwasEffectInBlock, numKeptInds, noscale);
     }

    if (matchedGWAS) {
        for (unsigned j=0; j<numIncdSnps; ++j)
            if (incdSnpInfoVec[j]->gwas_n != numKeptInds)
                throw string("Same-sample SVD GWAS requires common per-SNP sample counts.");
        nGWASblock.setConstant(numKeptLDBlocks, numKeptInds-1.0);
        ypy = (numKeptInds-1.0)*varPhenotypic;
    }

    QblocksDat.clear();
    SnpInfo *snp = NULL;
    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
        LDBlockInfo* ldblock = keptLdBlockInfoVec[i];
        MatrixDat matrixDat = MatrixDat(ldblock->snpNameVec, Qblocks[i]);
        QblocksDat.push_back(matrixDat);
        
    }

    if(gwasInfoOnly) {
        LOGGER  << "............................" << endl;
        LOGGER << "Summary of model input parameters:" << endl;
        LOGGER  << "............................" << endl;
        LOGGER << boost::format("%60s %8s %8s\n") %"" %"Mean" %"SD";
        LOGGER << boost::format("%60s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance" %Gadget::calcMean(varySnp) %sqrt(Gadget::calcVariance(varySnp));
        LOGGER << boost::format("%60s %8.3f %8.3f\n") %"GWAS SNP heterozygosity" %Gadget::calcMean(snp2pq) %sqrt(Gadget::calcVariance(snp2pq));
        LOGGER << boost::format("%60s %8.0f %8.0f\n") %"GWAS SNP sample size" %Gadget::calcMean(n) %sqrt(Gadget::calcVariance(n));
        LOGGER << boost::format("%60s %8.0f %8.0f\n") %"GWAS SNP allele flip number (ref: GWAS LD)" %numGWASFlip %0;
        LOGGER << boost::format("%60s %8.3f %8.3f\n") %"GWAS SNP effect (in genotype SD unit)" %Gadget::calcMean(b) %sqrt(Gadget::calcVariance(b));
        LOGGER << boost::format("%60s %8.3f %8.3f\n") %"GWAS SNP SE" %Gadget::calcMean(se) %sqrt(Gadget::calcVariance(se));
        LOGGER << boost::format("%60s %8.3f %8.3f\n") %"GWAS LD block size" %Gadget::calcMean(numSnpsBlock) %sqrt(Gadget::calcVariance(numSnpsBlock));
        LOGGER << boost::format("%60s %8.3f %8.3f\n") %"GWAS block rank" %Gadget::calcMean(numEigenvalBlock) %sqrt(Gadget::calcVariance(numEigenvalBlock));
        LOGGER << "............................" << endl;
        LOGGER << "Begin to model inference process..." << endl;
        LOGGER  << "............................" << endl;
    }
    lowRankModel = true;
}


void Data::includeMatchedBlocks(const bool haveEqtlInfo){
    // this step is to construct gwasSnp2geneVec
//    LOGGER << "Matching blocks..." << endl;
    SnpInfo * snp = NULL;
    LDBlockInfo * ldblock = NULL;
    
    for (unsigned i=0; i<numLDBlocks; ++i){
        ldblock = ldBlockInfoVec[i];
        ldblock->block2GwasSnpVec.clear();
        ldblock->snpInfoVec.clear();
    }
    for (unsigned j=0; j<numIncdSnps; ++j){
        snp = incdSnpInfoVec[j];
        ldblock = ldBlockInfoMap[snp->block];
        ldblock->block2GwasSnpVec.push_back(j);
        ldblock->snpInfoVec.push_back(snp);
        if(haveEqtlInfo){
            if(snp->iseQTL){
                ldblock->block2EqtlBoolVec.push_back(true);
            } else {
                ldblock->block2EqtlBoolVec.push_back(false);
            }
        }
    }
    for (unsigned i=0; i<numLDBlocks; ++i){
        ldblock = ldBlockInfoVec[i];
        if(ldblock->block2GwasSnpVec.size() == 0){
            ldblock->kept = false;
        } else {
            ldblock->startSnpIdx = ldblock->block2GwasSnpVec[0];
            ldblock->endSnpIdx = ldblock->block2GwasSnpVec[ldblock->numSnpInBlock-1];
            ldblock->kept = true;
        }
    }
        
    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();

    ldblock2gwasSnpMap.clear();
//    LOGGER << "Construct map from ld to snp" << endl;
    for(unsigned i = 0; i < numKeptLDBlocks; i++){
        ldblock = keptLdBlockInfoVec[i];
        ldblock2gwasSnpMap.insert(pair<int, vector<int> > (i,ldblock->block2GwasSnpVec));
    }
    
    // LOGGER << numKeptLDBlocks << " GWAS LD blocks are included." << endl;
}


//void Data::constructWandQ(const double eigenCutoff, const bool noscale){
//    VectorXd nMinusOne;
//    snp2pq.resize(numIncdSnps);
//    D.resize(numIncdSnps);
//   // ZPZdiag.resize(numIncdSnps);
//    ZPy.resize(numIncdSnps);
//    b.resize(numIncdSnps);
//    n.resize(numIncdSnps);
//    nMinusOne.resize(numIncdSnps);
//    se.resize(numIncdSnps);
//    tss.resize(numIncdSnps);
//    SnpInfo *snp;
//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        snp = incdSnpInfoVec[i];
//        snp->af = snp->gwas_af;
//        snp2pq[i] = snp->twopq = 2.0f*snp->gwas_af*(1.0-snp->gwas_af);
//        if(snp2pq[i]==0) LOGGER << "Error: SNP " << snp->rsID << " af " << snp->af << " has 2pq = 0." << endl;
//        D[i] = snp2pq[i]*snp->gwas_n;
//        b[i] = snp->gwas_b * sqrt(snp2pq[i]); // scale the marginal effect so that it's in per genotype SD unit
//        n[i] = snp->gwas_n;
//        nMinusOne[i] = snp->gwas_n - 1;
//        se[i]= snp->gwas_se * sqrt(snp2pq[i]);
//        tss[i] = D[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
//        ZPy[i] = n[i]*b[i];
//        //  D[i] = 1.0/(se[i]*se[i]+b[i]*b[i]/snp->gwas_n);  // NEW!
//        //  snp2pq[i] = snp->twopq = D[i]/snp->gwas_n;       // NEW!
//    }
//    LOGGER << endl;
//    LDBlockInfo * ldblock;
//    wcorrBlocks.resize(numKeptLDBlocks);
//    // Qblocks.resize(numKeptLDBlocks);
//    numSnpsBlock.resize(numKeptLDBlocks);
//    numEigenvalBlock.resize(numKeptLDBlocks);
//    Qblocks.clear();
//    VectorXd sqrtLambda;
//    // save gwas marginal effect into block
//    gwasEffectInBlock.resize(numKeptLDBlocks);
//    for (unsigned i = 0; i < numKeptLDBlocks; i++){
//        ldblock = keptLdBlockInfoVec[i];
//        gwasEffectInBlock[i] = b(ldblock->block2GwasSnpVec);
//        // calculate wbcorr and Qblocks
//        sqrtLambda = eigenValLdBlock[i].array().sqrt();
//        //LOGGER << eigenVecLdBlock[i].transpose().rows() << " " << eigenVecLdBlock[i].transpose().cols() << " " << gwasEffectInBlock[i].size() << endl;
//        wcorrBlocks[i] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * gwasEffectInBlock[i] );
//        // LOGGER << "eigenVecLdBlock[i]: " << eigenVecLdBlock[i] << endl;
//        // LOGGER << "wcorrBlocks[i]: " << wcorrBlocks[i] << endl;
//        // LOGGER << "sqrtLambda: " << sqrtLambda << endl;
//        // LOGGER << gwasEffectInBlock[i] << endl;
//        MatrixXd tmpQblocks = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
//        MatrixDat matrixDat = MatrixDat(ldblock->snpNameVec,tmpQblocks );
//        // LOGGER << "Qblock: " << endl;
//        // LOGGER << matrixDat.values << endl;
//        Qblocks.push_back(matrixDat);
//        numSnpsBlock[i] = Qblocks[i].ncol;
//        numEigenvalBlock[i] = Qblocks[i].nrow;
//    }
//    
//    //b.array() -= b.mean();  // DO NOT CENTER b
//    // estimate phenotypic variance based on the input allele frequencies in GWAS
//
//    //  Vp_buf = h_buf * N_buf * se_buf * se_buf + h_buf * b_buf * b_buf * N_buf / (N_buf - 1.0);
//    //VectorXd ypySrt = D.array()*n.array()*se.array().square() + D.array() * b.array().square() * n.array() / nMinusOne.array();
//    VectorXd ypySrt = D.array()*(n.array()*se.array().square()+b.array().square());
//    VectorXd varpSrt = ypySrt.array()/n.array();
//    std::sort(ypySrt.data(), ypySrt.data() + ypySrt.size());
//    std::sort(varpSrt.data(), varpSrt.data() + varpSrt.size());
//    ypy = ypySrt[ypySrt.size()/2];  // median
//    varPhenotypic = varpSrt[varpSrt.size()/2];
//    //LOGGER << "varPhenotypic: " << varPhenotypic << endl;
//    VectorXd nSrt = n;
//    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
//    numKeptInds = nSrt[nSrt.size()/2]; // median
//    
//    nGWASblock.resize(numKeptLDBlocks);
//    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
//        nGWASblock[i] = numKeptInds;
//    }
//
//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        snp = incdSnpInfoVec[i];
//        D[i] = varPhenotypic/(se[i]*se[i]+b[i]*b[i]/snp->gwas_n);  // NEW!
//        snp2pq[i] = snp->twopq = D[i]/snp->gwas_n;       // NEW!
//        tss[i] = D[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
//        // Need to adjust R and C models X'X matrix depending scale of genotypes or not
//        if (noscale == true) {
//            D[i] = snp2pq[i]*snp->gwas_n;
//        } else {
//            D[i] = snp->gwas_n;
//        }
//    }
//    // LOGGER << endl << "snp2pq: " << endl << snp2pq << endl;
//    //ypy = numKeptInds;
//    // NEW END
//    // data summary
//    LOGGER << "\nData summary:" << endl;
//    LOGGER << boost::format("%40s %8s %8s\n") %"" %"mean" %"sd";
//    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance" %Gadget::calcMean(varpSrt) %sqrt(Gadget::calcVariance(varpSrt));
//    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP heterozygosity" %Gadget::calcMean(snp2pq) %sqrt(Gadget::calcVariance(snp2pq));
//    LOGGER << boost::format("%40s %8.0f %8.0f\n") %"GWAS SNP sample size" %Gadget::calcMean(n) %sqrt(Gadget::calcVariance(n));
//    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP effect (in genotype SD unit)" %Gadget::calcMean(b) %sqrt(Gadget::calcVariance(b));
//    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP SE" %Gadget::calcMean(se) %sqrt(Gadget::calcVariance(se));
//    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"LD block size" %Gadget::calcMean(numSnpsBlock) %sqrt(Gadget::calcVariance(numSnpsBlock));
//    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"LD block rank" %Gadget::calcMean(numEigenvalBlock) %sqrt(Gadget::calcVariance(numEigenvalBlock));
//}


void Data::mergeLdmInfo(const string &outLDmatType, const string &dirname) {
    
    if (outLDmatType != "block") LOGGER.e(0, " --merge-ldm-info only works for block LD matrices at the moment!");
    
    string dir_path = dirname; // Replace with your folder path
    string search_str = outLDmatType;
    DIR* dirp = opendir(dir_path.c_str());
    
    if (dirp == NULL) {
        LOGGER.e(0, " opening directory [" + dirname + "]");
    }
    
    // find out all ldm in the folder
    vector<string> file_list;
    
    dirent* dp;
    while ((dp = readdir(dirp)) != NULL) {
        string file_name = dp->d_name;
        if (std::regex_match(file_name, std::regex("block[0-9]+\\.ldm\\.info"))) {
            file_list.push_back(file_name);
        }
    }
    
    closedir(dirp);
    
    set<unsigned> blockIdxSet;
    blockIdxSet.clear();
    
    for (vector<string>::iterator it = file_list.begin(); it != file_list.end(); ++it) {
        size_t block_pos = it->find(search_str);
        size_t dot_pos = it->find_first_of(".", block_pos);
        if (block_pos != string::npos && dot_pos != string::npos) {
            string block_num_str = it->substr(block_pos + search_str.size(), dot_pos - block_pos - search_str.size());
            int block_num = atoi(block_num_str.c_str());
            blockIdxSet.insert(block_num);
        }
    }
    
    if (blockIdxSet.size() == 0) {
        LOGGER.e(0, " there is no info file to merge in folder [" + dirname + "].");
    }
    
    unsigned nldm = blockIdxSet.size();
    
    set<unsigned>::iterator it = blockIdxSet.begin();
    
    string outSnpInfoFile = dirname + "/snp.info";
    string outldmInfoFile = dirname + "/ldm.info";
    
    ofstream out1(outSnpInfoFile.c_str());
    out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12s %10s %10s\n")
    % "Chrom"
    % "ID"
    % "Index"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "N"
    % "Block";
    
    ofstream out2(outldmInfoFile.c_str());
    out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
    % "Block"
    % "Chrom"
    % "StartSnpIdx"
    % "StartSnpID"
    % "EndSnpIdx"
    % "EndSnpID"
    % "NumSnps";
    
    
    
    unsigned snpIdx = 0;
    unsigned ldmIdx = 0;
    
    for (unsigned i=0; i<nldm; ++i) {
        
        string snpInfoFile = dirname + "/block" + to_string(*it) + ".snp.info";
        string ldmInfoFile = dirname + "/block" + to_string(*it) + ".ldm.info";
        
        // read snp info file
        ifstream in1(snpInfoFile.c_str());
        if (!in1) LOGGER.e(0, " can not open the file [" + snpInfoFile + "] to read.");
        cout << "Reading SNP info from file [" + snpInfoFile + "]." << endl;
        string header;
        string id, allele1, allele2;
        int chr, physPos,ld_n;
        float genPos;
        double allele1Freq;
        int idx = 0;
        int index;
        string blockID;
        map<string, int> snpID2index;
        getline(in1, header);
        while (in1 >> chr >> id >> index >> genPos >> physPos >> allele1 >> allele2 >> allele1Freq >> ld_n >> blockID) {
            out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %22.17g %10s %10s\n")
            % chr
            % id
            % snpIdx
            % genPos
            % physPos
            % allele1
            % allele2
            % allele1Freq
            % ld_n
            % blockID;
            snpID2index[id] = snpIdx;
            ++snpIdx;
        }
        in1.close();
        
        
        // read ldm info file
        ifstream in2(ldmInfoFile.c_str());
        if (!in2) LOGGER.e(0, " can not open the file [" + ldmInfoFile + "] to read.");
        cout << "Reading LDM info from file [" + ldmInfoFile + "]." << endl;
        
        int blockStart, blockEnd, snpNum;
        string startSnpID, endSnpID;
        getline(in2, header);
        while (in2 >> id >> chr >> blockStart >> startSnpID >> blockEnd >> endSnpID >> snpNum) {
            out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
            % id
            % chr
            % snpID2index[startSnpID]
            % startSnpID
            % snpID2index[endSnpID]
            % endSnpID
            % snpNum;
            
            ++ldmIdx;
        }
        in2.close();
        
        ++it;
        
        // Keep per-block metadata: merging is repeatable and non-destructive.

    }
    
    out1.close();
    out2.close();
    
    cout << "Written " << snpIdx << " SNPs info into file [" + outSnpInfoFile + "]." << endl;
    cout << "Written " << ldmIdx << " LDMs info into file [" + outldmInfoFile + "]." << endl;
    
}

 int Data::parseChrNumber(const std::string& s) {
    if (s == "X" || s == "x") return 23;
    if (s == "Y" || s == "y") return 24;
    if (s == "M" || s == "m" || s == "MT" || s == "mt") return 25;
    return std::stoi(s);
}

#include <dirent.h>
#include <sys/stat.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include "hdf5.h"

// ============================================================
// Parse block number from filename like:
//   block9.ldm.bin
// Return ID = "9"
// ============================================================
std::optional<LDBlockInfo>
Data::makeLDBlockFromFilename(const std::string& fname,
                              int index,
                              const std::string& marker) {
    // Require prefix, e.g. "block".
    if (fname.rfind(marker, 0) != 0) {
        return std::nullopt;
    }

    const std::string suffix = ".ldm.bin";

    // Require content between marker and suffix.
    if (fname.size() <= marker.size() + suffix.size()) {
        return std::nullopt;
    }

    // Require the expected suffix.
    if (fname.compare(fname.size() - suffix.size(),
                      suffix.size(),
                      suffix) != 0) {
        return std::nullopt;
    }
    // Extract the block ID between marker and ".ldm.bin".
    // Examples:
    //   block9.ldm.bin
    //       -> ID = "9"
    //   blockENSG00000283683.ldm.bin
    //       -> ID = "ENSG00000283683"
    //   blockENSG00000284922_ENSG00000149357.ldm.bin
    //       -> ID = "ENSG00000284922_ENSG00000149357"
    const std::string token = fname.substr(
        marker.size(),
        fname.size() - marker.size() - suffix.size()
    );

    if (token.empty()) {
        return std::nullopt;
    }

    // Store the ID as a string. It may be either a numeric block ID
    // or a single/merged gene ID.
    LDBlockInfo block(index, token, -1);
    block.startPos = -1;
    block.endPos = -1;

    return block;
}
// ============================================================
// Read metadata only from blockX.snp.info
// Fill:
//   chrom
//   startPos
//   endPos
//
// Note:
// Do NOT populate global snpInfoVec here.
// This function is only for block-level metadata.
// ============================================================
bool Data::loadBlockMetaFromSnpInfo(LDBlockInfo* ldblock,
                                    const std::string& snpInfoFile) {
    if (ldblock == nullptr) return false;

    std::ifstream in(snpInfoFile.c_str());
    if (!in) {
        LOGGER.e(0, "Error: cannot open SNP info file [" + snpInfoFile + "]");
        return false;
    }

    std::string header;
    std::getline(in, header);

    // block*.snp.info columns:
    // Chrom ID Index GenPos PhysPos A1 A2 A1Freq N Block
    int chr = -1;
    int physPos = -1;
    int index = 0;
    int ld_n = 0;
    float genPos = 0.0f;
    float af = 0.0f;
    std::string id, allele1, allele2, blockID;

    bool firstSnp = true;
    while (in >> chr >> id >> index >> genPos >> physPos
              >> allele1 >> allele2 >> af >> ld_n >> blockID) {
        if (firstSnp) {
            ldblock->chrom = chr;
            ldblock->startPos = physPos;
            firstSnp = false;
        }
        ldblock->endPos = physPos;
    }
    in.close();

    if (firstSnp) {
        LOGGER << "Warning: no SNPs found in [" << snpInfoFile << "]." << std::endl;
        return false;
    }

    return true;
}

// ============================================================
// Read full SNP list from blockX.snp.info
// Fill:
//   ldblock->snpNameVec
//   ldblock->block2GwasSnpVec
//
// Also append into global:
//   snpInfoVec
//   snpInfoMap
//
// Note:
// The SNP order in blockX.snp.info is assumed to follow physical position.
// This is the order kept in HDF5 snplist.
// ============================================================
bool Data::loadBlockSnpsFromSnpInfo(LDBlockInfo* ldblock,
                                    const std::string& snpInfoFile) {
    if (ldblock == nullptr) return false;

    std::ifstream in(snpInfoFile.c_str());
    if (!in) {
        LOGGER.e(0, "Error: cannot open SNP info file [" + snpInfoFile + "]");
        return false;
    }

    ldblock->snpNameVec.clear();
    ldblock->block2GwasSnpVec.clear();

    std::string header;
    std::getline(in, header);

    int chr = -1;
    int physPos = -1;
    int snpIndexInFile = 0;
    int ld_n = 0;
    float genPos = 0.0f;
    float af = 0.0f;
    std::string id, allele1, allele2, blockID;

    bool firstSnp = true;
    int localIdx = 0;

    while (in >> chr >> id >> snpIndexInFile >> genPos >> physPos
              >> allele1 >> allele2 >> af >> ld_n >> blockID) {
        SnpInfo* snp = new SnpInfo(localIdx++, id, allele1, allele2, chr, genPos, physPos);
        snp->af = af;
        snp->block = blockID;

        snpInfoVec.push_back(snp);
        ldblock->snpNameVec.push_back(id);
        ldblock->block2GwasSnpVec.push_back(physPos);
        chromosomes.insert(snp->chrom);

        if (!snpInfoMap.insert({id, snp}).second) {
            LOGGER.e(0, "Duplicate SNP ID: \"" + id + "\".");
        }

        if (firstSnp) {
            ldblock->chrom = chr;
            ldblock->startPos = physPos;
            firstSnp = false;
        }
        ldblock->endPos = physPos;
    }
    in.close();

    if (ldblock->snpNameVec.empty()) {
        LOGGER << "Warning: empty SNP list in [" << snpInfoFile << "]." << std::endl;
        return false;
    }

    return true;
}

// ============================================================
// Read LD matrix from blockX.ldm.bin
// ============================================================
bool Data::loadBlockLdmBinary(const std::string& ldmBinFile,
                              int32_t blockSize,
                              Eigen::MatrixXf& ldm) {
    FILE* fp = fopen(ldmBinFile.c_str(), "rb");
    if (!fp) {
        LOGGER.e(0, "Error: cannot open [" + ldmBinFile + "]");
        return false;
    }

    ldm.resize(blockSize, blockSize);
    const uint64_t nElem = static_cast<uint64_t>(blockSize) * static_cast<uint64_t>(blockSize);

    if (fread(ldm.data(), sizeof(float), nElem, fp) != nElem) {
        fclose(fp);
        LOGGER.e(0, "Size error in " + ldmBinFile);
        return false;
    }

    fclose(fp);
    return true;
}

// ============================================================
// Write one block into HDF5 group
//   /blk_k/snplist
//   /blk_k/ldblk
// ============================================================
bool Data::writeOneBlockToHdf5(hid_t file_id,
                               const std::string& groupName,
                               const std::vector<std::string>& snpNames,
                               const Eigen::MatrixXf& ldm) {
    hid_t group = H5Gcreate(file_id, groupName.c_str(),
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group < 0) {
        LOGGER.e(0, "Error: cannot create group [" + groupName + "]");
        return false;
    }

    const int32_t blockSize = static_cast<int32_t>(snpNames.size());

    // snplist
    hsize_t dims1[1] = { static_cast<hsize_t>(blockSize) };
    hid_t space1 = H5Screate_simple(1, dims1, NULL);

    hid_t strType = H5Tcopy(H5T_C_S1);
    H5Tset_size(strType, H5T_VARIABLE);

    std::vector<const char*> cstrs;
    cstrs.reserve(blockSize);
    for (const auto& s : snpNames) cstrs.push_back(s.c_str());

    hid_t dset1 = H5Dcreate(group, "snplist", strType, space1,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset1, strType, H5S_ALL, H5S_ALL, H5P_DEFAULT, cstrs.data());

    H5Dclose(dset1);
    H5Sclose(space1);
    H5Tclose(strType);

    // ldblk
    hsize_t dims2[2] = {
        static_cast<hsize_t>(blockSize),
        static_cast<hsize_t>(blockSize)
    };
    hid_t space2 = H5Screate_simple(2, dims2, NULL);

    hid_t dset2 = H5Dcreate(group, "ldblk", H5T_NATIVE_FLOAT, space2,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset2, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
             H5P_DEFAULT, ldm.data());

    H5Dclose(dset2);
    H5Sclose(space2);

    H5Gclose(group);
    return true;
}

// ============================================================
// Main function
//
// dirname : input folder containing
//           block*.ldm.bin / block*.snp.info
//
// outname : output prefix for per-chromosome HDF5 files,
//           e.g. "/tmp/ldblk_reference"
//
// Output:
//   outname_chr1.ld.hdf5
//   outname_chr2.ld.hdf5
//   ...
//   <outname parent directory>/ld.snpinfo
//
// Block IDs may be numeric or string-based gene IDs, for example:
//   block9.ldm.bin
//   blockENSG00000283683.ldm.bin
//   blockENSG00000284922_ENSG00000149357.ldm.bin
// ============================================================

void Data::mergeLdmInfoIntoPrsCSHdf5(
    const string &outLDmatType,
    const string &dirname,
    const string &outname
) {
    const std::string dirPath = dirname;
    const std::string blockPrefix = "block";
    const std::string ldmSuffix = ".ldm.bin";

    boost::filesystem::path outPrefix(outname);
    const std::string outPath = (outPrefix.parent_path() / "").string();

    DIR* dirp = opendir(dirPath.c_str());

    if (dirp == nullptr) {
        LOGGER.e(0, " opening directory [" + dirname + "]");
    }

    // --------------------------------------------------------
    // 1. Find all block*.ldm.bin files
    // --------------------------------------------------------
    std::vector<std::string> fileList;
    dirent* dp = nullptr;

    while ((dp = readdir(dirp)) != nullptr) {
        const std::string fileName = dp->d_name;

        const bool hasPrefix =
            fileName.rfind(blockPrefix, 0) == 0;

        const bool hasSuffix =
            fileName.size() > blockPrefix.size() + ldmSuffix.size() &&
            fileName.compare(
                fileName.size() - ldmSuffix.size(),
                ldmSuffix.size(),
                ldmSuffix
            ) == 0;

        if (hasPrefix && hasSuffix) {
            fileList.push_back(fileName);
        }
    }

    closedir(dirp);

    // Initial deterministic filename ordering.
    // Final HDF5 block ordering is determined later by chromosome
    // and physical position.
    std::sort(fileList.begin(), fileList.end());

    // --------------------------------------------------------
    // 2. Build block list from filenames and block*.snp.info
    // --------------------------------------------------------
    ldBlockInfoVec.clear();
    ldBlockInfoMap.clear();

    std::set<int> chromVec;
    std::multimap<int, LDBlockInfo*> chrToLDBlockMap;

    int idx = 0;

    for (const auto& fileName : fileList) {
        auto optBlock = Data::makeLDBlockFromFilename(
            fileName,
            idx,
            blockPrefix
        );

        if (!optBlock) {
            LOGGER << "Skip invalid LD filename ["
                   << fileName << "]." << std::endl;
            continue;
        }

        LDBlockInfo* ldblock =
            new LDBlockInfo(std::move(*optBlock));

        const std::string id = ldblock->ID;

        // Prevent duplicate block IDs.
        if (!ldBlockInfoMap.insert({id, ldblock}).second) {
            LOGGER << "Skip duplicated block ID ["
                   << id << "]." << std::endl;

            delete ldblock;
            continue;
        }

        const std::string snpInfoFile =
            dirname + "/" +
            blockPrefix + ldblock->ID +
            ".snp.info";

        if (!loadBlockMetaFromSnpInfo(ldblock, snpInfoFile)) {
            LOGGER << "Skip block [" << id
                   << "] because metadata cannot be loaded from ["
                   << snpInfoFile << "]." << std::endl;

            ldBlockInfoMap.erase(id);
            delete ldblock;
            continue;
        }

        ldBlockInfoVec.push_back(ldblock);
        chromVec.insert(ldblock->chrom);
        chrToLDBlockMap.insert({
            ldblock->chrom,
            ldblock
        });

        ++idx;
    }

    numLDBlocks =
        static_cast<unsigned>(ldBlockInfoVec.size());

    includeMatchedBlocks();

    LOGGER << numLDBlocks
           << " LD Blocks to be included from folder ["
           << dirname << "]."
           << std::endl;

    if (numLDBlocks == 0) {
        LOGGER.e(
            0,
            " there is no valid LD block to merge in folder [" +
            dirname + "]."
        );
    }

    // --------------------------------------------------------
    // 3. Reset SNP containers for merged SNP output
    // --------------------------------------------------------
    snpInfoVec.clear();
    snpInfoMap.clear();

    // --------------------------------------------------------
    // 4. Process each chromosome
    //
    // Chromosomes are traversed in ascending numeric order
    // because chromVec is std::set<int>.
    //
    // Within each chromosome, blocks are sorted by:
    //   1. startPos ascending
    //   2. endPos ascending
    //   3. string block ID ascending as tie-breaker
    //
    // Then valid blocks are written as:
    //   /blk_1
    //   /blk_2
    //   ...
    // --------------------------------------------------------
    for (const int chr : chromVec) {
        std::cout << "Chromosome " << chr << ":\n";

        const std::string hdf5File =
            outname +
            "_chr" +
            std::to_string(chr) +
            ".ld.hdf5";

        hid_t fileId = H5Fcreate(
            hdf5File.c_str(),
            H5F_ACC_TRUNC,
            H5P_DEFAULT,
            H5P_DEFAULT
        );

        if (fileId < 0) {
            LOGGER.e(
                0,
                "Error: cannot create HDF5 file [" +
                hdf5File + "]"
            );
            continue;
        }

        // Gather all blocks for the current chromosome.
        std::vector<LDBlockInfo*> chrBlocks;

        const auto range =
            chrToLDBlockMap.equal_range(chr);

        for (auto it = range.first;
             it != range.second;
             ++it) {
            chrBlocks.push_back(it->second);
        }

        // Sort blocks by physical position.
        // Gene IDs are strings and are used only as a final tie-breaker.
        std::sort(
            chrBlocks.begin(),
            chrBlocks.end(),
            [](const LDBlockInfo* a,
               const LDBlockInfo* b) {
                if (a->startPos != b->startPos) {
                    return a->startPos < b->startPos;
                }

                if (a->endPos != b->endPos) {
                    return a->endPos < b->endPos;
                }

                return a->ID < b->ID;
            }
        );

        int blkCount = 0;

        for (LDBlockInfo* ldblock : chrBlocks) {
            ldblock->snpNameVec.clear();
            ldblock->block2GwasSnpVec.clear();

            const std::string snpInfoFile =
                dirname + "/" +
                blockPrefix + ldblock->ID +
                ".snp.info";

            if (!loadBlockSnpsFromSnpInfo(
                    ldblock,
                    snpInfoFile
                )) {
                LOGGER << "Skip block ["
                       << ldblock->ID
                       << "] because SNP information cannot be loaded from ["
                       << snpInfoFile
                       << "]."
                       << std::endl;
                continue;
            }

            const int32_t blockSize =
                static_cast<int32_t>(
                    ldblock->snpNameVec.size()
                );

            if (blockSize <= 1) {
                LOGGER << "Skip block ["
                       << ldblock->ID
                       << "] because blockSize <= 1."
                       << std::endl;
                continue;
            }

            const std::string ldmBinFile =
                dirname + "/" +
                blockPrefix + ldblock->ID +
                ldmSuffix;

            Eigen::MatrixXf ldm;

            if (!loadBlockLdmBinary(
                    ldmBinFile,
                    blockSize,
                    ldm
                )) {
                LOGGER << "Skip block ["
                       << ldblock->ID
                       << "] because LD matrix cannot be loaded from ["
                       << ldmBinFile
                       << "]."
                       << std::endl;
                continue;
            }

            const std::string groupName =
                "/blk_" +
                std::to_string(blkCount + 1);

            if (!writeOneBlockToHdf5(
                    fileId,
                    groupName,
                    ldblock->snpNameVec,
                    ldm
                )) {
                LOGGER << "Skip block ["
                       << ldblock->ID
                       << "] because writing HDF5 group ["
                       << groupName
                       << "] failed."
                       << std::endl;
                continue;
            }

            // Only increment after successful HDF5 writing.
            ++blkCount;

            LOGGER << "Wrote block ["
                   << ldblock->ID
                   << "] as group ["
                   << groupName
                   << "] in chromosome "
                   << chr
                   << "."
                   << std::endl;
        }

        H5Fclose(fileId);

        LOGGER << "Chromosome "
               << chr
               << ": wrote "
               << blkCount
               << " blocks to file ["
               << hdf5File
               << "]."
               << std::endl;
    }

    // --------------------------------------------------------
    // 5. Write merged SNP information
    // --------------------------------------------------------
    numSnps =
        static_cast<unsigned>(snpInfoVec.size());

    includeMatchedSnp();

    const std::string outSnpInfoFile =
        outPath + "/ld.snpinfo";

    std::ofstream out(outSnpInfoFile.c_str());

    if (!out.is_open()) {
        LOGGER.e(
            0,
            "Error: cannot create SNP info file [" +
            outSnpInfoFile + "]"
        );
    }

    out << boost::format(
        "%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
    )
        % "CHR"
        % "SNP"
        % "BP"
        % "A1"
        % "A2"
        % "A1Freq"
        % "BLOCK";

    for (unsigned i = 0;
         i < numIncdSnps;
         ++i) {
        SnpInfo* snp = incdSnpInfoVec[i];

        out << boost::format(
            "%s\t%s\t%s\t%s\t%s\t%.6f\t%s\n"
        )
            % snp->chrom
            % snp->rsID
            % snp->physPos
            % snp->a1
            % snp->a2
            % snp->af
            % snp->block;
    }

    out.close();

    LOGGER << "Wrote total "
           << numIncdSnps
           << " SNP records into ["
           << outSnpInfoFile
           << "]."
           << std::endl;

    LOGGER
        << "\n## Standards for constructing LD blocks ##\n"
        << "1. SNPs must map to the same chromosome\n"
        << "2. Duplicate SNP entries are removed\n"
        << "3. SNPs within each block follow input physical-position order\n"
        << "4. Chromosomes are processed in ascending numeric order\n"
        << "5. Blocks within each chromosome are sorted by startPos and endPos\n"
        << "6. String block IDs are used only as a final ordering tie-breaker\n"
        << "7. Blocks containing one or zero SNPs are excluded\n";
}


void Data::mergeLdmInfoIntoCovmHdf5(const string &outLDmatType, const string &dirname, const string &outname) {
    
    string dir_path = dirname; // Replace with your folder path
    string search_str;
    if(outLDmatType == "covm") search_str = "block";
    DIR* dirp = opendir(dir_path.c_str());
    
    if (dirp == NULL) {LOGGER.e(0, " opening directory [" + dirname + "]");}
    // find out all ldm in the folder
    vector<string> file_list;
    
    dirent* dp;
    while ((dp = readdir(dirp)) != NULL) {
        string file_name = dp->d_name;
        if (std::regex_match(file_name, std::regex("block[0-9]+\\.ldm\\.info"))) {
            file_list.push_back(file_name);
        }
    }
    
    closedir(dirp);
    std::vector<LDBlockInfo> ldBlockInfo;
    // set<string> blockIdxSet;
    // blockIdxSet.clear();

    int idx = 0;
    LDBlockInfo *ldblock = NULL;
    ldBlockInfoVec.clear();
    ldBlockInfoMap.clear();
    for (const auto& fname : file_list) {
        if (auto optBlock = Data::makeLDBlockFromFilename(fname, idx)) {
            const std::string& id = optBlock->ID;
            // if (blockIdxSet.count(id)) continue;
            ldblock = new LDBlockInfo(std::move(*optBlock)); 
            if (ldBlockInfoMap.insert(pair<string, LDBlockInfo*>(id, ldblock)).second == false) continue;
            ldBlockInfoVec.push_back(ldblock);
            ++idx;
        }
    }

    numLDBlocks = ldBlockInfoVec.size();
    includeMatchedBlocks();
    LOGGER << numLDBlocks << " Genotype covariance Blocks to be included from folder [" + dirname + "]." << endl;


    // std::filesystem::path outSnpInfoFile = std::filesystem::path(outdirname) / "ldm.hdf5.snpinfo";

    if (numLDBlocks == 0) {LOGGER.e(0, " there is no info file to merge in folder [" + dirname + "].");}
    // std::filesystem::path hdf5File = std::filesystem::path(outdirname) / "ldm.hdf5";
    string outSnpInfoFile = outname + ".snp.info";
    
    // unsigned nldm = blockIdxSet.size();
    // set<string>::iterator it = blockIdxSet.begin();
    unsigned snpIdx = 0;
    unsigned ldmIdx = 0;
    snpInfoVec.clear();
    snpInfoMap.clear();
    SnpInfo *snp = NULL;
    // loop chroms

    std::string hdf5File = outname + ".covm.hdf5";
    /************************************************************
    * LD HDF5 Storage Specification (Column-major Upper-triangular)
    *
    * Each LD block is stored as a group `/blk_k` containing:
    *
    *   1. "snplist" : string array of SNP IDs, length = n.
    *   2. "ldblk_utri" : 1D float array, storing the upper-triangular
    *      (including diagonal) elements of the LD matrix.
    *      Flattening order is **column-major**:
    *         - Fix column j = 0..n-1
    *         - For each column j, store rows i = 0..j
    *      Example for n=3:
    *
    *           [ a11  a12  a13
    *                  a22  a23
    *                       a33 ]
    *
    *      Stored as: [a11, a12, a22, a13, a23, a33]
    *
    *   3. Attribute "nSNPs" : integer scalar attached to "ldblk_utri",
    *      indicating the number of SNPs (matrix dimension).
    *
    * Notes:
    *   - This saves ~50% memory compared to full n×n storage.
    *   - When reading in R/Python, reconstruct full symmetric LD:
    *       mat[i,j] = mat[j,i] = stored[k] (iterate in column-major order).
    *   - Keep SNP order consistent with "snplist".
    ************************************************************/


    // open new HDF5 file for this chromosome
    hid_t file_id = H5Fcreate(hdf5File.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) { LOGGER.e(0, "Error: cannot create HDF5 file [" + hdf5File + "]"); }

    // iterate over LD blocks for this chromosome
    for (unsigned lbs=0; lbs < numLDBlocks; ++lbs) {
        ldblock = ldBlockInfoVec[lbs];
        ldblock->snpNameVec.clear();
        // read SNP info
        std::string snpInfoFile = dirname + "/block" + ldblock->ID + ".snp.info";
        std::ifstream in1(snpInfoFile.c_str());
        if (!in1) { LOGGER.e(0, "Error: cannot open SNP info file [" + snpInfoFile + "]"); continue; }
        std::string header,id,allele1,allele2,blockID; int chr,physPos,ld_n,index,idx=0; float genPos,af;
        std::getline(in1, header);
        while (in1 >> chr >> id >> index >> genPos >> physPos >> allele1 >> allele2 >> af >> ld_n >> blockID) {
            snp = new SnpInfo(idx++, id, allele1, allele2, chr, genPos, physPos);
            snp->af=af; snp->block=blockID;
            snpInfoVec.push_back(snp);
            snp->windStart = ldblock->startPos;
            snp->windEnd = ldblock->endPos;
            snp->block = ldblock->ID;
            snp->ld_n = ld_n;
            snp->WindWidth = ldblock->endPos - ldblock->startPos;
            ldblock->snpNameVec.push_back(id);
            ldblock->block2GwasSnpVec.push_back(physPos);
            chromosomes.insert(snp->chrom);
            if (!snpInfoMap.insert({id, snp}).second) LOGGER.e(0,"Duplicate SNP ID: \""+id+"\".");
        }
        in1.close();

        // read LD matrix from binary file
        std::string ldmBinFile = dirname + "/block" + ldblock->ID + ".covm.bin";
        FILE* fp=fopen(ldmBinFile.c_str(),"rb");
        if(!fp){ LOGGER.e(0,"Error: cannot open ["+ldmBinFile+"]"); continue; }
        int32_t blockSize = ldblock->snpNameVec.size();
        Eigen::MatrixXf ldm(blockSize,blockSize);
        uint64_t nElem=(uint64_t)blockSize*blockSize;
        if(fread(ldm.data(),sizeof(float),nElem,fp)!=nElem) LOGGER.e(0,"Size error in "+ldmBinFile);
        fclose(fp);

        // create group for this block
        std::string groupName="/blk_"+std::to_string(ldblock->index+1);
        hid_t group=H5Gcreate(file_id,groupName.c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        if(group<0){ LOGGER.e(0,"Error: cannot create group "+groupName); continue; }

        // write SNP list as variable-length strings
        hsize_t dims1[1] = { static_cast<hsize_t>(blockSize) };
        hid_t space1=H5Screate_simple(1,dims1,NULL);
        hid_t strType=H5Tcopy(H5T_C_S1); H5Tset_size(strType,H5T_VARIABLE);
        std::vector<const char*> cstrs; cstrs.reserve(blockSize);
        for(auto& s:ldblock->snpNameVec) cstrs.push_back(s.c_str());
        hid_t dset1=H5Dcreate(group,"snplist",strType,space1,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Dwrite(dset1,strType,H5S_ALL,H5S_ALL,H5P_DEFAULT,cstrs.data());
        H5Dclose(dset1); H5Sclose(space1); H5Tclose(strType);

        std::vector<float> tri_data= Gadget::packUpperTriColMajor(ldm);
        // write LD as 1D
        hsize_t dims2[1] = { (hsize_t)tri_data.size() };
        hid_t space2 = H5Screate_simple(1,dims2,NULL);
        hid_t dset2 = H5Dcreate(group,"blk_covm_utri",H5T_NATIVE_FLOAT,space2,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Dwrite(dset2,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,H5P_DEFAULT,tri_data.data());

        // add attribute nSNPs
        hid_t attr_space=H5Screate(H5S_SCALAR);
        hid_t attr=H5Acreate(dset2,"nSNPs",H5T_NATIVE_INT,attr_space,H5P_DEFAULT,H5P_DEFAULT);
        int32_t n=ldm.rows(); H5Awrite(attr,H5T_NATIVE_INT,&n);
        H5Aclose(attr); 

        // add attribute chr
        hid_t attr_chr = H5Acreate(dset2,"chr",H5T_NATIVE_INT,attr_space,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(attr_chr,H5T_NATIVE_INT,&ldblock->chrom);
        H5Aclose(attr_chr);
        
        H5Sclose(attr_space);H5Dclose(dset2); H5Sclose(space2);


        H5Gclose(group);
        LOGGER << "Wrote SNP list and LD matrix for block ["+ldblock->ID+"]" << std::endl;
    }
    H5Fclose(file_id);
    LOGGER << "Written to file [" + hdf5File + "] (Column-major Upper-triangular for block matrices)." << std::endl;


    // update snp info
    numSnps = (unsigned) snpInfoVec.size();
    includeMatchedSnp();
    ofstream out1(outSnpInfoFile.c_str());
    // out1 << boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\n")
    //  % "CHR" % "SNP" % "BP" % "A1" % "A2" % "A1Freq" % "BLOCK";
    // for (unsigned i=0; i < numIncdSnps ; ++i) {
    //     snp = incdSnpInfoVec[i];
    //     out1 << boost::format("%s\t%s\t%s\t%s\t%s\t%.6f\t%s\n")
    //         % snp->chrom
    //         % snp->rsID
    //         % snp->physPos
    //         % snp->a1
    //         % snp->a2
    //         % snp->af
    //         % snp->block;
    // }
    // out1.close();
    out1 << boost::format("%6s %15s %10s %15s %6s %6s %12s %10s %10s %10s %10s %15s %10s %12s %12s\n")
    % "Chrom"
    % "ID"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "Index"
    % "WindStart"
    % "WindEnd"
    % "WindSize"
    % "WindWidth"
    % "N"
    % "SamplVar"
    % "LDsum";
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        out1 << boost::format("%6s %15s %10s %15s %6s %6s %12f %10s %10s %10s %10s %15s %10s %12.6f %12.6f\n")
        % snp->chrom
        % snp->rsID
        % snp->genPos
        % snp->physPos
        % snp->a1
        % snp->a2
        % snp->af
        % snp->index
        % snp->windStart
        % snp->windEnd
        % ldBlockInfoMap[snp->block]->snpNameVec.size()
        % snp->WindWidth
        % snp->ld_n
        % snp->ldSamplVar
        % snp->ldSum;
    }
    out1.close();

    LOGGER << "Write total " << numIncdSnps << " SNP info into the single file [" + outSnpInfoFile + "]." << endl;
    LOGGER << "\n##Standards for constructing LD blocks:##\n"
         << "1. SNPs must map within the same chromosome\n"
         << "2. Remove duplicate SNP entries\n"
         << "3. SNP must satisfy start <= pos <= stop\n"
         << "4. If SNP falls in multiple blocks, keep the later one\n"
         << "5. Drop blocks with only 1 SNP\n"
         << "6. Each SNP maps to only one block (unique mapping)\n"
         << "7. SNPs on boundaries are included, but drop block if only 1 boundary SNP\n"
         ;
}



void Data::mergeBlockGwasSummary(const string &gwasSummaryFile, const string &title) {
    string dir_path = "."; // Replace with your folder path
    string search_str1 = gwasSummaryFile + ".block";
    string search_str2 = ".ma";
    DIR* dirp = opendir(dir_path.c_str());
    
    if (dirp == NULL) {
        LOGGER.e(0,"Error opening directory [" + dir_path + "]");
    }
    
    // find out all .ma files in the folder
    vector<string> file_list;
    
    dirent* dp;
    while ((dp = readdir(dirp)) != NULL) {
        string file_name = dp->d_name;
        if (file_name.find(search_str1) != string::npos && file_name.find(search_str2) != string::npos) {
            file_list.push_back(file_name);
        }
    }
    
    closedir(dirp);
    
    map<unsigned, string> blockIdxMap;
    blockIdxMap.clear();
    
    for (vector<string>::iterator it = file_list.begin(); it != file_list.end(); ++it) {
        size_t block_pos = it->find(search_str1);
        size_t dot_pos = it->find_first_of(".", block_pos);
        if (block_pos != string::npos && dot_pos != string::npos) {
            string block_num_str = it->substr(block_pos + search_str1.size(), dot_pos - block_pos - search_str1.size());
            int block_num = atoi(block_num_str.c_str());
            blockIdxMap[block_num] = *it;
        }
    }
    
    if (blockIdxMap.size() == 0) {
        LOGGER.e(0,"there is no info file to merge in folder [" + dir_path + "].");
    }
    
    unsigned nBlk = blockIdxMap.size();
    
    LOGGER << "Merging GWAS summary statistics files across " + to_string(nBlk) + " blocks..." << endl;
    
    map<unsigned, string>::iterator it = blockIdxMap.begin();
    
    string outMaFile = dir_path + "/" + title + ".ma";
    
    ofstream out(outMaFile.c_str());
    out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n") % "SNP" % "A1" % "A2" % "freq" % "b" % "se" % "p" % "N";
    
    unsigned snpIdx = 0;
    
    for (unsigned i=0; i<nBlk; ++i) {
        
        string mafile = it->second;
        
        // read snp info file
        ifstream in(mafile.c_str());
        if (!in) LOGGER.e(0,"can not open the file [" + mafile + "] to read.");
        LOGGER << "Reading summary statistics from file [" + mafile + "]." << endl;

        string header;
        getline(in, header);

        string id, allele1, allele2, freq, b, se, pval, n;
        while (in >> id >> allele1 >> allele2 >> freq >> b >> se >> pval >> n) {
            out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n")
            % id
            % allele1
            % allele2
            % freq
            % b
            % se
            % pval
            % n;
            ++snpIdx;
        }
        in.close();
        
        ++it;
    }
    
    out.close();
    
    LOGGER << "Written " << snpIdx << " SNPs info into file [" + outMaFile + "]." << endl;
}

void Data::constructPseudoSummaryData(){
    LOGGER << "Constructing pseudo summary statistics for training and validation data sets, with 90% sample size for training and 10% for validation." << endl;
    
    pseudoGwasEffectTrn.resize(numKeptLDBlocks);
    pseudoGwasEffectVal.resize(numKeptLDBlocks);
    
    double n_trn = 0.9*double(numKeptInds);
    double n_val = numKeptInds - n_trn;
    pseudoGwasNtrn = n_trn;
    b_val.resize(numIncdSnps);

    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
        LDBlockInfo* block = keptLdBlockInfoVec[i];
        
        long size = eigenValLdBlock[i].size();
        VectorXd rnd(size);
        for (unsigned j=0; j<size; ++j) {
            rnd[j] = Stat::snorm();
        }
        
        pseudoGwasEffectTrn[i] = gwasEffectInBlock[i] + sqrt(1.0/n_trn - 1.0/nGWASblock[i]) * eigenVecLdBlock[i] * (eigenValLdBlock[i].array().sqrt().matrix().asDiagonal() * rnd);

        pseudoGwasEffectVal[i] = nGWASblock[i]/n_val * gwasEffectInBlock[i] - n_trn/n_val * pseudoGwasEffectTrn[i];
        b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
    }
    
    
}

void Data::constructWandQ(const vector<VectorXd> &GWASeffects, const double nGWAS, const bool noscale) {
    wcorrBlocks.resize(numKeptLDBlocks);
    numSnpsBlock.resize(numKeptLDBlocks);
    numEigenvalBlock.resize(numKeptLDBlocks);
    Qblocks.clear();

    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        // calculate wbcorr and Qblocks
        VectorXd sqrtLambda = eigenValLdBlock[i].array().sqrt();
        //LOGGER << eigenVecLdBlock[i].transpose().rows() << " " << eigenVecLdBlock[i].transpose().cols() << " " << gwasEffectInBlock[i].size() << endl;
        wcorrBlocks[i] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * GWASeffects[i] );
        // LOGGER << "eigenVecLdBlock[i]: " << eigenVecLdBlock[i] << endl;
        // LOGGER << "wcorrBlocks[i]: " << wcorrBlocks[i] << endl;
        // LOGGER << "sqrtLambda: " << sqrtLambda << endl;
        // LOGGER << gwasEffectInBlock[i] << endl;
//        MatrixXd tmpQblocks = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].trßanspose();
//        MatrixDat matrixDat = MatrixDat(ldblock->snpNameVec, tmpQblocks);
//        // LOGGER << "Qblock: " << endl;
//        // LOGGER << matrixDat.values << endl;
//        Qblocks.push_back(matrixDat);
        Qblocks[i] = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        
        if (noscale) {
            VectorXd Dsqrt(ldblock->numSnpInBlock);
            for (unsigned j=0; j<ldblock->numSnpInBlock; ++j) {
                SnpInfo *snp = ldblock->snpInfoVec[j];
                Dsqrt[j] = sqrt(snp->twopq);
            }
            Qblocks[i] = Qblocks[i] * Dsqrt.asDiagonal();
        }
        
        numSnpsBlock[i] = Qblocks[i].cols();
        numEigenvalBlock[i] = Qblocks[i].rows();
        
        eigenVecLdBlock[i].resize(0,0);
    }
    
    nGWASblock.resize(numKeptLDBlocks);
    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        nGWASblock[i] = nGWAS;
    }
}

void Data::scaleGwasEffects(){
    snp2pq.resize(numIncdSnps);
    scalingGWASFactorVec.resize(numIncdSnps);
    b.resize(numIncdSnps);
    n.resize(numIncdSnps);
    D.resize(numIncdSnps);
    se.resize(numIncdSnps);
    tss.resize(numIncdSnps); // only used in SBayesC
    ZPy.resize(numIncdSnps);
    varySnp.resize(numIncdSnps);

    SnpInfo *snp = NULL;
    int lineIdx = 0;
    auto startTime = std::chrono::steady_clock::now();
    for (unsigned i=0; i<numIncdSnps; ++i) {
        Gadget::showProgressBar(i, numIncdSnps, startTime,"Scale GWAS effects");
        double varpsPG=0,varpsDivid2pqPG = 0;
        snp = incdSnpInfoVec[i];
        snp->af = snp->gwas_af;
        snp2pq[i] = snp->twopq = 2.0f*snp->gwas_af*(1.0-snp->gwas_af);
        if(snp2pq[i]==0) LOGGER << "Error: SNP " << snp->rsID << " af " << snp->af << " has 2pq = 0." << endl;
        n[i] = snp->gwas_n;
        D[i] = snp2pq[i]*snp->gwas_n;
        b[i] = snp->gwas_b; // * sqrt(snp2pq[i]); // scale the marginal effect so that it's in per genotype SD unit
        se[i]= snp->gwas_se ; //* sqrt(snp2pq[i]);
        tss[i] = n[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
        ZPy[i] = n[i]*b[i];

        // start to scale
        if(lineIdx == 0) LOGGER.w(0,"Assuming non-scaled GWAS effect and se is used. Please double-check it.");
        // calculate scaling factor
        scalingGWASFactorVec[i] = sqrt(1/(n[i]* se[i] * se[i] + b[i]*b[i]));
        snp->scaleFactor = scalingGWASFactorVec[i];
        // check beta and se scaled or not
        // When generating beta effect under standardized genotype based on beta effect from GWAS using genotypes at 0/1/2 scale
        // one way is to time sqrt(snp2pq), in this way we also need to set phenotypic variance as the median of each SNP-based
        // phenotypic variance.
        // Another way is assuming phenotypic variance as 1, and then timeing a per-SNP constant scalingGWASFactorVec[i] 
        // if(false){
            // varpsPG =  D[i]*(n[i]*se[i]*se[i] + b[i]*b[i])/n[i];
            // varpsDivid2pqPG =  D[i]*(n[i]*se[i]*se[i] + b[i]*b[i])/(n[i] * snp2pq[i]);
            // b[i] = snp->gwas_b  * sqrt(snp2pq[i]); // scale the marginal effect so that it's in per genotype SD unit
            // se[i]= snp->gwas_se * sqrt(snp2pq[i]);
            // varySnp[i] = varpsPG;
        // }
        b[i] = snp->gwas_b  * (snp->scaleFactor); // scale the marginal effect so that it's in per genotype SD unit
        se[i]= snp->gwas_se * (snp->scaleFactor);
        varySnp[i] = 1.0;
        if (matchedGWAS) {
            if (!(suppliedPhenotypicVariance > 0.0) || n[i] <= 2)
                throw string("Same-sample SVD GWAS requires observed sample variance and n > 2.");
            const double sd = sqrt(suppliedPhenotypicVariance / (snp->gwas_b*snp->gwas_b + (n[i]-2.0)*snp->gwas_se*snp->gwas_se));
            scalingGWASFactorVec[i] = snp->scaleFactor = sd;
            b[i] = snp->gwas_b*sd;
            se[i] = snp->gwas_se*sd;
            snp2pq[i] = 1.0;
            varySnp[i] = suppliedPhenotypicVariance;
        }
        if (!suppliedGenotypeScale.empty()) {
            const double sd=suppliedGenotypeScale.at(snp->rsID);
            const double vx=suppliedGenotypeVariance.at(snp->rsID);
            scalingGWASFactorVec[i]=snp->scaleFactor=sd;
            b[i]=snp->gwas_b*sd*vx; // covariance score for Q'Q, not marginal slope
            se[i]=snp->gwas_se*sd;
            snp2pq[i]=vx;
            varySnp[i]=vx*sd*sd*(snp->gwas_b*snp->gwas_b+(n[i]-2)*snp->gwas_se*snp->gwas_se);
        }
        /////////////////////////////////////
        lineIdx ++;
    }
    // Estimate sample size
    VectorXd nSrt = n;
    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
    numKeptInds = Gadget::findMedian(nSrt);
    varPhenotypic = Gadget::findMedian(varySnp);
    // LOGGER << "varPhenotypic: " << varPhenotypic << endl;

    // divide gwas marginal effect into blocks
    if (numKeptLDBlocks) {
        gwasEffectInBlock.resize(numKeptLDBlocks);
        nGWASblock.resize(numKeptLDBlocks);
        for (unsigned i = 0; i < numKeptLDBlocks; i++){
            LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
            gwasEffectInBlock[i] = b(ldblock->block2GwasSnpVec);
            nGWASblock[i] = numKeptInds;
        }
    }
}
