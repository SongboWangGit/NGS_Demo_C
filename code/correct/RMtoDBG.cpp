//
// Created by sbwang on 19-4-9.
//

#include "RMtoDBG.h"
#include "util/nthash.hpp"
#include <unistd.h>

RMtoDBG::RMtoDBG(string file_name, uint16_t k) {
    cdbg_file_path = file_name;
    ksize = k;
}
//RMtoDBG::RMtoDBG(map<string, uint64_t> h) {
//    this->overlap_hash = h;
//}

/*
 * 函数功能：字符串分割函数
 */
void RMtoDBG::split_string(const std::string& s, std::vector<std::string>& v, const std::string& c)
{
    std::string::size_type pos1, pos2;
    pos2 = s.find(c);
    pos1 = 0;
    while(std::string::npos != pos2)
    {
        v.push_back(s.substr(pos1, pos2-pos1));

        pos1 = pos2 + c.size();
        pos2 = s.find(c, pos1);
    }
    if(pos1 != s.length())
        v.push_back(s.substr(pos1));
}

string RMtoDBG::prefix(string str, uint32_t length) {
    return str.substr(0, length);
}

string RMtoDBG::sufix(string str, uint32_t length) {
    return str.substr(str.size() - length, length);
}


/*
 * 把dbg写入文件中
 */
void RMtoDBG::write_dbg_to_file(unordered_map<uint64_t , DBG_Node> &bucket_dbg, string file_name) {

    auto dbg_iter = bucket_dbg.begin();
    ofstream fout(file_name);

    while(dbg_iter != bucket_dbg.end()) {

        string kmer = dbg_iter->second.kmer;

        uint64_t node_num = dbg_iter->first;
        vector<uint64_t> pre_node = dbg_iter->second.pre_node;
        vector<uint64_t> suf_node = dbg_iter->second.suf_node;

        fout << ">" << node_num << " pre:";
        for (int i = 0; i < pre_node.size(); i++){
            fout <<  pre_node[i] << ":";
        }
        fout << " suf:";
        for (int i = 0; i < suf_node.size(); i++){
            fout  << suf_node[i] << ":";
        }

        fout << " l_lone:" << dbg_iter->second.left_lonely << " r_lone:" << dbg_iter->second.right_lonely;
        fout << endl;
        fout << kmer << endl;

        dbg_iter ++;
    }

    fout.close();
}


void RMtoDBG::get_from_file() {

    ifstream cDBG_file(cdbg_file_path);
    string line;

    uint64_t node_num;

    // 从文件中获取cDBG
    while (getline(cDBG_file, line)){
        if (line[0] == '>'){

            line = line.substr(1, line.size() - 1);

            vector<string> split_line;
            split_string(line, split_line, " ");

            node_num = strtoull (split_line[0].c_str(), nullptr, 0);
            //获得kmer
            string kmer;
            getline(cDBG_file, kmer);

            // 获得pre nodes
            vector<string> pre_node_str;
            vector<uint64_t > pre_node;
            pre_node.reserve(pre_node_str.size());

            split_string(split_line[1], pre_node_str, ":");
            pre_node_str.erase(pre_node_str.begin());
            for (int i = 0; i < pre_node_str.size(); i ++){
                pre_node.emplace_back(strtoull (pre_node_str[i].c_str(), nullptr, 0));
            }

            if (!pre_node.empty()){
                overlap_hash[prefix(kmer, ksize - 1)].emplace_back(node_num);
            }


            // 获得suf nodes
            vector<string> suf_node_str;
            vector<uint64_t > suf_node;
            suf_node.reserve(suf_node_str.size());

            split_string(split_line[2], suf_node_str, ":");
            suf_node_str.erase(suf_node_str.begin());
            for (int i = 0; i < suf_node_str.size(); i++){
                suf_node.emplace_back(strtoull (suf_node_str[i].c_str(), nullptr, 0));
            }

            if (!suf_node.empty()){
                overlap_hash[sufix(kmer, ksize - 1)].emplace_back(node_num);
            }


            DBG_Node cur_node = {kmer, pre_node, suf_node, 0};

            bucket_dbg[node_num] = cur_node;

        }
    }

    write_dbg_to_file(bucket_dbg, "reunite/test.txt");

}


void RMtoDBG::map_read_muti_unitig(string read, int n) {
    cout << "all: " << overlap_hash.size() << endl;
    int start_flag = 0;
    int end_flag = 0;
    vector<uint64_t > start_overlap_pos;
    vector<uint64_t > end_overlap_pos;

    for (int i = 0; i < n; i++){
        if (start_flag == 0){

            string start_overlap = read.substr(i, ksize - 1);
            find_overlap(start_overlap, start_overlap_pos);

            // 如果begin 找到了对应的unitig
            if (!start_overlap_pos.empty()) {
                start_flag = 1;

                // 开始找end对应的unitig
                for (int j = 0; j < n; j++){
                    if (end_flag == 0){
                        string end_overlap = read.substr(read.size() - ksize + 1 - i, ksize - 1);

                        find_overlap(end_overlap, end_overlap_pos);

                        if (!end_overlap_pos.empty()){
                            end_flag = 1;
                            break;
                        }
                    }
                }
            }
        }
    }

    // start和end都找到才可以继续进行



}

void RMtoDBG::find_overlap(string overlap, vector<uint64_t > &overlap_pos) {

    overlap_pos.clear();
    auto overlap_iter = overlap_hash.find(overlap);

    if (overlap_iter != overlap_hash.end()){
        for (int i = 0; i < overlap_iter->second.size(); i++){
            overlap_pos.emplace_back(overlap_iter->second[i]);
            cout << "---" << overlap_iter->second[i] << endl;
        }
    }

}
//int main(){
//
//    time_t t = time(nullptr);
//
//
//    string read = "CTCCAATGCACTCCCTTCCATGAACTTCCATTTTTCCCCACTCCTTT";
//    int n = 3;
//    RMtoDBG rMtoDBG("reunite/final_dbg.txt", 16);
//    rMtoDBG.get_from_file();
//    rMtoDBG.map_read_muti_unitig(read, n);
//
//
//    showMemStat(getpid());
//    double cost_t = time(nullptr) - t;
//    cout << "all time : " << cost_t  << "s" << endl;
//}