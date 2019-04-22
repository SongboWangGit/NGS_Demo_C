//
// Created by sbwang on 19-4-9.
//

#ifndef NGS_DEMO_RMTODBG_H
#define NGS_DEMO_RMTODBG_H

#include <vector>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include "util/Usage.hpp"
#include <set>
#include <htslib/sam.h>
#include <htslib/hts.h>

using namespace std;

class RMtoDBG{

private:
    struct DBG_Node{
        string kmer;
        vector<uint64_t> pre_node;
        vector<uint64_t> suf_node;
        int flag;
        int left_lonely = 0;
        int right_lonely = 0;

    };
    uint16_t ksize;
    int max_step;

    unordered_map<uint64_t, DBG_Node> bucket_dbg;
    unordered_map<int, int> t;

    unordered_map<string, vector<uint64_t> > overlap_hash;
    string cdbg_file_path;

public:
    RMtoDBG(string file_name, uint16_t ksize, int m_step);
//    RMtoDBG(map<string , uint64_t> h);

    string prefix(string str, uint32_t length);
    string sufix(string str, uint32_t length);

    void get_from_file();
    void write_dbg_to_file(unordered_map<uint64_t, DBG_Node>& bucket_dbg, string file_name);

    void split_string(const std::string& s, std::vector<std::string>& v, const std::string& c);
    void map_read_muti_unitig(string read);
    void find_overlap(string overlap, vector<uint64_t > &overlap_pos);
    void find_all_paths(DBG_Node cur_node, string cur_path, unsigned long read_size, vector<string> &map_paths);

    void trim_paths(vector<string> all_map_paths, set<string> &all_trimed_paths, string start_str, uint64_t start_overlap_in_read_pos, int read_len);
    string find_best_path(string trimed_read, set<string> all_trimed_path);
};

#endif //NGS_DEMO_RMTODBG_H