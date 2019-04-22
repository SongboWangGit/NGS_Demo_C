//
// Created by sbwang on 19-4-16.
//

#ifndef NGS_DEMO_KMERCHOOSE_H
#define NGS_DEMO_KMERCHOOSE_H

#include "KmerStat.h"

class KmerChoose{

public:

    void static kmer_stat_start(string in_file, int ksize, int s, int r);
    void kmer_stat_thread(string in_file, size_t thread_num, int s, int r);
    uint64_t choose_ksize();

};
#endif //NGS_DEMO_KMERCHOOSE_H
