//
// Created by sbwang on 19-4-9.
//

#include "RMtoDBG.h"
#include "util/nthash.hpp"
#include <unistd.h>

RMtoDBG::RMtoDBG(string cdbg_file, uint16_t k, int m_step, string err_file) {
    cdbg_file_path = cdbg_file;
    ksize = k;
    max_step = m_step;
    err_info_path = err_file;
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

//            if (!pre_node.empty()){
                overlap_hash[prefix(kmer, ksize - 1)].emplace_back(node_num);
//            }


            // 获得suf nodes
            vector<string> suf_node_str;
            vector<uint64_t > suf_node;
            suf_node.reserve(suf_node_str.size());

            split_string(split_line[2], suf_node_str, ":");
            suf_node_str.erase(suf_node_str.begin());
            for (int i = 0; i < suf_node_str.size(); i++){
                suf_node.emplace_back(strtoull (suf_node_str[i].c_str(), nullptr, 0));
            }

//            if (!suf_node.empty()){
                overlap_hash[sufix(kmer, ksize - 1)].emplace_back(node_num);
//            }


            DBG_Node cur_node = {kmer, pre_node, suf_node, 0};

            bucket_dbg[node_num] = cur_node;

        }
    }

    write_dbg_to_file(bucket_dbg, "reunite/test.txt");

}

void RMtoDBG::find_overlap(string overlap, vector<uint64_t > &overlap_pos) {

    overlap_pos.clear();
    auto overlap_iter = overlap_hash.find(overlap);

    if (overlap_iter != overlap_hash.end()){
        for (int i = 0; i < overlap_iter->second.size(); i++){
            overlap_pos.emplace_back(overlap_iter->second[i]);
//            cout << "---" << bucket_dbg[overlap_iter->second[i]].kmer << endl;
        }
    }
}

string RMtoDBG::map_read_muti_unitig(string read, uint64_t &start_overlap_in_read_pos) {
//    cout << "all: " << overlap_hash.size() << endl;
    int start_flag = 0;
    int end_flag = 0;

    vector<uint64_t > start_overlap_in_dbg_poses;
    vector<uint64_t > end_overlap_pos;


    int read_size = read.size();
    string start_overlap;

    for (int i = 0; i < max_step; i++){
        start_overlap = read.substr(i, ksize - 1);
        find_overlap(start_overlap, start_overlap_in_dbg_poses);
        if (not start_overlap_in_dbg_poses.empty()){
            start_overlap_in_read_pos = i;
            break;
        }
    }

    // 找到所有可能的路径
    if (not start_overlap_in_dbg_poses.empty()){
        vector<string> all_map_paths;

        for (int i = 0; i < start_overlap_in_dbg_poses.size(); i++){
            DBG_Node cur_node = bucket_dbg[start_overlap_in_dbg_poses[i]];
            vector<uint64_t > pre_nodes = cur_node.pre_node;
            string cur_path;

            find_all_paths(cur_node, cur_path, read.size(), all_map_paths);
//            if (pre_nodes.empty()){
//                find_all_paths(cur_node, cur_path, read.size(), all_map_paths);
//            } else {
//                // 加上前面一个节点来保证冗余，足够长
//                for (int j = 0; j < pre_nodes.size(); j++){
//                    cur_path = bucket_dbg[pre_nodes[j]].kmer;
//                    // 递归来寻找一条完整的路径
//                    find_all_paths(cur_node, cur_path, read.size(), all_map_paths);
//                }
//
//            }

        }

        string trimed_read = read.substr(start_overlap_in_read_pos);

        set<string> all_trimed_paths;
        trim_paths(all_map_paths, all_trimed_paths, start_overlap, start_overlap_in_read_pos, read_size);

        string best_path = find_best_path(trimed_read, all_trimed_paths);


        return best_path;
    }

    return "";

}


// 向后找节点加入路径
void RMtoDBG::find_all_paths(DBG_Node cur_node, string cur_path, unsigned long read_size, vector<string> &map_paths){
    vector<uint64_t> suf_nodes = cur_node.suf_node;

    cur_path += cur_node.kmer;

    // 总长度大于read的两倍的时候停止
    if (cur_path.size() > read_size * 2){
        map_paths.emplace_back(cur_path);
        return;
    }

    // 不再有后缀节点的时候停止
    if (suf_nodes.empty()){
        map_paths.emplace_back(cur_path);
        return;
    }
    // 对每一个后缀节点，使用递归来扩展路径
    for (int i = 0; i < suf_nodes.size(); i ++){
        DBG_Node next_node = bucket_dbg[suf_nodes[i]];
        // 递归来寻找一条完整的路径
        find_all_paths(next_node, cur_path, read_size, map_paths);
    }
}

/*
 * 函数作用：修剪找到的path，使得与read的长度相同
 */
void RMtoDBG::trim_paths(vector<string> all_map_paths, set<string> &all_trimed_paths, string start_overlap, uint64_t start_overlap_in_read_pos, int read_len) {

    for (int i = 0; i < all_map_paths.size(); i++){

        string cur_path = all_map_paths[i];
        if (cur_path.size() < read_len){
            all_trimed_paths.insert(cur_path);
            continue;
        }
        // 从overlap在path上的位置开始修剪
        unsigned long trim_pos = cur_path.find(start_overlap);
        // 加上前半段
//        string trimed_path = cur_path.substr(trim_pos - start_overlap.size(), start_overlap.size());
        // 加上后半段
        string trimed_path = cur_path.substr(trim_pos, read_len - start_overlap_in_read_pos);

        all_trimed_paths.insert(trimed_path);
//         cout << "trimed path: " << trimed_path << endl;
    }

}

/*
 * 函数功能：从所有修剪过得path中找到错配数最小的路径，并输出
 */
string RMtoDBG::find_best_path(string trimed_read, set<string> all_trimed_paths) {

    int min_mismatch = 200;
    string min_path;


    for(auto it = all_trimed_paths.begin (); it != all_trimed_paths.end ();it++){
        int cur_mismatch = 0;
        string cur_path = *it;

        for (int i = 0; i < trimed_read.size(); i++){
            if (trimed_read[i] != cur_path[i]){
                cur_mismatch ++;
            }

            if (cur_mismatch < min_mismatch){
                min_mismatch = cur_mismatch;
                min_path = *it;
            }
        }
    }

    return min_path;

}

void RMtoDBG::correct(string best_path, vector<string> err_segment, uint64_t err_pos, uint64_t bias) {

    ofstream err_correction("err_correction.txt", ios::app);
    err_correction << err_pos << " " << best_path;

    for (int i = 0; i < err_segment.size(); i ++){
        vector<Correct_Pair> correction;

        string cur_err_segment = err_segment[i];

        // 遍历每一条，每一条的纠正量大于阈值，则丢弃
        for (int j = 0; j < best_path.size() and correction.size() <= 5; j++){
            if (cur_err_segment[j] != best_path[j]){
                uint64_t err_bp_pos = err_pos + bias + j;
                Correct_Pair new_correct = {err_bp_pos, cur_err_segment[j], best_path[j]};
                correction.emplace_back(new_correct);
            }
        }

        //每一条的纠正量大于阈值，则丢弃, 否则将更改写入文件
        if (correction.size() >= 5){
            continue;
        } else if (correction.empty()) {
            err_correction << " None";
//            continue;
        }else {
                err_correction << " " << cur_err_segment << ":";
                for (int k = 0; k < correction.size(); k++){
                    err_correction << " " << correction[k].pos << "-" << correction[k].wrong_bp << "->" << correction[k].right_bp;
                }
        }

    }
    err_correction << endl;
    err_correction.close();
}
void RMtoDBG::map_start() {
    ifstream err_file(err_info_path);
    string line;

    // 创建输出文件
    ofstream err_correction("err_correction.txt");
    err_correction.close();

    while (getline(err_file, line)){
        vector<string> line_split;
        split_string(line, line_split, " ");

        vector<string> err_segment;
        err_segment.reserve(line_split.size() - 1);
        uint64_t err_pos = stoll(line_split[0]);

        for (int i = 1; i < line_split.size(); i ++){
            err_segment.emplace_back(line_split[i]);
        }

        uint64_t start_overlap_in_read_pos = 0;

        string best_path = map_read_muti_unitig(err_segment[0], start_overlap_in_read_pos);

        if (not best_path.empty()){
            for (int i = 0; i < err_segment.size(); i ++){
                err_segment[i] = err_segment[i].substr(start_overlap_in_read_pos);
            }

            correct(best_path, err_segment, err_pos, start_overlap_in_read_pos);
            cout << best_path << endl;

        }


//        break;
    }
}
int main(){

    time_t t = time(nullptr);

    int max_step = 20;
    RMtoDBG rMtoDBG("reunite/final_dbg.txt", 31, max_step, "err_info.txt");
    rMtoDBG.get_from_file();
    rMtoDBG.map_start();
//
//    uint64_t cnt = 0;
//    string in_file = "../../bams/sampleChr20_sorted.bam";
//    samFile *bam_file = hts_open(in_file.c_str(), "r");
//    bam1_t *aln = bam_init1(); //initialize an alignment
//    bam_hdr_t *bamHdr = sam_hdr_read(bam_file); //read header
//    while(sam_read1(bam_file, bamHdr, aln) > 0) {
//        //获得序列
//        int32_t len = aln->core.l_qseq;
//        uint8_t *q = bam_get_seq(aln);
//        string read;
//        for (int i = 0; i < len; i++) {
//            read += seq_nt16_str[bam_seqi(q, i)];
//        }
//
//        if (read.size() != 150){
//            continue;
//        }
//        read = read.substr(30, 40);
//        rMtoDBG.map_read_muti_unitig(read);
//
//        cnt ++;
//        if (cnt == 15){
//            cout << "tot: " << cnt << endl;
//            break;
//        }
//
//    }

    showMemStat(getpid());
    double cost_t = time(nullptr) - t;
    cout << "all time : " << cost_t  << "s" << endl;

}