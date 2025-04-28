// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BayesOmics, a statistical genetics software package
// developed by Shouye Liu.
//
// Copyright (C) 2025 Shouye Liu <syliu.xue@foxmail.com>
//
// This file is licensed under the GNU General Public License v3.0 or (at your option)
// any later version. You may redistribute and/or modify it under the terms of the GPL.
//
// BayesOmics is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the LICENSE file for details.



#ifndef macro_hpp
#define macro_hpp

// Define the macro only for Linux when compile static binary file
// macro ARMA_DONT_USE_WRAPPER should be placed before #include <armadillo>
#if defined(__linux__) || defined(__linux) || defined(linux) || defined(__unix__)
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>
////// general //////////
// GIT_COMMIT_HASH is passed during compilation with -D option
#ifndef GIT_COMMIT_HASH
#define GIT_COMMIT_HASH "unknown"
#endif
// oneline doc link
// #define DOC_ONLINE "https://shouyeliu.github.io/"
#define DOC_ONLINE ""
///////////////
#define MAXSNPNUMPERPROBEINSPARSE 0x300000
#define MAX_PROBE_NUM 0xF0000
#define MAX_SNP_NAME 64
#define DENSE_FILE_TYPE_1 0  // uint32 + float
#define SPARSE_FILE_TYPE_3F 0x40400000 // uint32 + uint64_t + uint64_ts + uint32_ts + float
#define SPARSE_FILE_TYPE_3 3 // 16*uint32s + uint64_t + uint64_ts + uint32_ts + float (indicator+samplesize+snpnumber+probenumber+ 6*-9s +valnumber+cols+rowids+betases) [default]
#define DENSE_FILE_TYPE_3 5  // 16*uint32s + doubles (indicator+samplesize+snpnumber+probenumber+ 6*-9s + values) [default]
#define RESERVEDUNITS 16
#define FNAMESIZE 4096
#define SNPMISSRATE 0.3
#define MIN_PVAL_ADJUSTED 1e-150
#define BayesOmics_VERSION "v0.0.1"
#define MISSING_XQTL_RATE_PER_GENE  0.3
// Define a macro for the time zone
#define TIME_ZONE "AEST" // Australian Eastern Standard Time
// Define color macros with "COLOR" suffix
#define RESET_COLOR   "\033[0m"
#define RED_COLOR     "\033[31m"
#define GREEN_COLOR   "\033[32m"
#define YELLOW_COLOR  "\033[33m"
#define BLUE_COLOR    "\033[34m"
#define MAGENTA_COLOR "\033[35m"
#define CYAN_COLOR    "\033[36m"
#define WHITE_COLOR   "\033[37m"


#endif /* macro_hpp */