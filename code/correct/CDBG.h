//
// Created by sbwang on 19-3-19.
//

#ifndef NGS_DEMO_CDBG_H
#define NGS_DEMO_CDBG_H

#include "util/Usage.hpp"
#include "util/nthash.hpp"
#include "util/ntHashIterator.hpp"
#include <iostream>
#include <cstdio>
#include <vector>
#include <unistd.h>
#include <algorithm>
#include <string>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <dirent.h>
#include <unistd.h>
#include <fstream>
#include <mutex>
#include <time.h>
#include <sys/stat.h>
#include <cstring>
#include "util/ThreadPool.h"
#include <unordered_map>
#include <cmath>
using namespace std;


class CDBG{

private:

    uint16_t lsize;
    uint16_t ksize;
    samFile *bam_file;

    uint16_t thread_num;
    uint8_t r;
    uint8_t threshold;
    unordered_map<string, string > buckets;
    struct dbg_node{
        string kmer;
        vector<uint64_t> pre_node;
        vector<uint64_t> suf_node;
        int flag;
        int left_lonely = 0;
        int right_lonely = 0;

    };


    struct lone_node{
        string left_lone = "";
        string right_lone = "";
        string kmer;
    };

//    static recursive_mutex reunited_lock;
//    static recursive_mutex all_unitig_lock;

    static unordered_map<string, lone_node> left_lone_unitigs;
    static unordered_map<string, lone_node> right_lone_unitigs;
    static unordered_map<string, lone_node> all_lone_unitigs;
    uint32_t *count_table;
//    static recursive_mutex left_lone_unitigs_lock;
//    static recursive_mutex right_lone_unitigs_lock;

public:
    CDBG();
    CDBG(std::string in_file, uint16_t ks, uint16_t ls, uint16_t thread_num, uint8_t thresh, uint8_t r);
    CDBG(CDBG &cdbg);

    string lmm(string kmer, int minl);
    string rmm(string kmer, int minl);
    string prefix(string str, uint32_t length);
    string sufix(string str, uint32_t length);

    void get_count_table();
    // 分割kmer到桶中所需要的功能函数
    void write_to_buckets(string out_file, string kmer);
    void write_buckets_to_file();
    void split_to_buckets();
    //

    void construct_dbg_with_thread();
    string compact_two_nodes(string u, string v);
    static void construct_dbg_in_bucket(CDBG* p, string file_name);
    void write_dbg_to_file(unordered_map<uint64_t, dbg_node>& bucket_dbg, string file_name);


    void compact_dbg_in_bucket(unordered_map<uint64_t, dbg_node>& bucket_dbg);
    void find_unitig_in_bucket(uint64_t last_num, dbg_node& cur_node, uint64_t cur_num, dbg_node& new_node, unordered_map<uint64_t, dbg_node> &bucket_dbg, uint64_t new_num);

    void erase_and_update_down(dbg_node cur_node, uint64_t cur_num, uint64_t new_num, unordered_map<uint64_t, dbg_node> &bucket_dbg);
    void erase_and_update_up_and_down(dbg_node cur_node, uint64_t cur_num, uint64_t new_num, unordered_map<uint64_t, dbg_node> &bucket_dbg);
    void set_lonely_in_bucket(unordered_map<uint64_t, dbg_node> &bucket_dbg, string file_name);

    void deal_left_lone();

    lone_node glue(lone_node u, lone_node v);

    void construct_final_dbg();

    void split_string(const std::string& s, std::vector<std::string>& v, const std::string& c);
    vector<string> get_files(string cate_dir);

};
#endif //NGS_DEMO_CDBG_H


//if (!left_k->second.right_lone.empty()){
//right_lone_unitigs[left_k->second.right_lone] = left_k->second;
//}
//// 虽然lonely找不到可以glue的node了
//else {
//all_file << left_k->second.kmer << endl;
//}
//left_lone_unitigs.erase(left_k);