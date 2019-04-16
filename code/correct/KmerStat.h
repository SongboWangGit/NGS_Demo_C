//
// Created by sbwang on 19-4-12.
//

#ifndef NGS_DEMO_KMERSTAT_H
#define NGS_DEMO_KMERSTAT_H


#include <string>
#include <iostream>
#include <vector>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include "util/nthash.hpp"
#include <stack>
#include "util/ntHashIterator.hpp"
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h>
#include <algorithm>
#include <cmath>
#include "util/Usage.hpp"
#include "util/ThreadPool.h"

using namespace std;

class KmerStat{

private:
    uint16_t ksize;
    samFile *bam_file;
    uint64_t table_size;
    uint32_t *count_table;
    uint16_t r;
    uint16_t s;

    uint64_t max_count;
public:
    KmerStat(string file, uint16_t k, uint16_t s_num, uint16_t r_num);

    uint64_t get_nthash(string read);
    void Update();
    void judge_and_set(uint64_t hash_num);
    void write_table_to_file();
    void write_f_to_file(double *f, uint64_t F0, uint64_t F1);
    uint64_t cal_F1(double *f);

    uint64_t ten_2_tow(uint64_t num);

    double* Estimate();
};
#endif //NGS_DEMO_KMERSTAT_H
