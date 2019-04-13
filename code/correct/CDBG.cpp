// Created by sbwang on 19-3-19.
//

#include <sstream>
#include "CDBG.h"

//定义，分配内存，以后A每一个对象（实例）的创建都不再分配内存
unordered_map<string, CDBG::lone_node> CDBG::left_lone_unitigs;
unordered_map<string, CDBG::lone_node> CDBG::right_lone_unitigs;
recursive_mutex left_lone_unitigs_lock;
recursive_mutex right_lone_unitigs_lock;
recursive_mutex all_unitig_lock;
recursive_mutex reunited_lock;

CDBG::CDBG() {

}
CDBG::CDBG(string in_file, uint16_t ks, uint16_t ls, uint16_t t_num) {
    lsize = ls;
    ksize = ks;
    bam_file = hts_open(in_file.c_str(), "r");
    thread_num = t_num;
    semaphore = Semaphore(thread_num);



}

CDBG::CDBG(CDBG &cdbg) {

    this->lsize = cdbg.lsize;
    this->ksize = cdbg.ksize;

}

/*
 * 函数功能：得到文件夹中所有的文件名
 */
vector<string> CDBG::get_files(string cate_dir) {
    vector<string> files;//存放文件名

    DIR *dir;
    struct dirent *ptr;
    char base[1000];

    if ((dir=opendir(cate_dir.c_str())) == nullptr)
    {
        perror("Open dir error...");
        exit(1);
    }

    while ((ptr=readdir(dir)) != nullptr)
    {
        if(strcmp(ptr->d_name,".")==0 || strcmp(ptr->d_name,"..")==0)    ///current dir OR parrent dir
            continue;
        else if(ptr->d_type == 8)    ///file
            //printf("d_name:%s/%s\n",basePath,ptr->d_name);
            files.emplace_back(ptr->d_name);
        else if(ptr->d_type == 10)    ///link file
            //printf("d_name:%s/%s\n",basePath,ptr->d_name);
            continue;
        else if(ptr->d_type == 4)    ///dir
        {
            files.emplace_back(ptr->d_name);

        }
    }
    closedir(dir);
    return files;
}

/*
 * 函数功能：字符串分割函数
 */
void CDBG::split_string(const std::string& s, std::vector<std::string>& v, const std::string& c)
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

string CDBG::prefix(string str, uint32_t length) {
    return str.substr(0, length);
}

string CDBG::sufix(string str, uint32_t length) {
    return str.substr(str.size() - length, length);
}
string CDBG::lmm(string kmer, int minl) {

    string prefix_k_1 = prefix(kmer, ksize - 1);
    string min_lmer = prefix_k_1.substr(0, minl);

    for(int i = 1; i <= prefix_k_1.size() - minl; i++){
        string cur_lmer = prefix_k_1.substr(i, minl);

        if (cur_lmer < min_lmer){
            min_lmer = cur_lmer;
        }
    }
    return min_lmer;
}

string CDBG::rmm(string kmer, int minl) {
    string sufix_k_1 = sufix(kmer, ksize - 1);
    string min_rmer = sufix_k_1.substr(0, minl);

    for(int i = 1; i <= sufix_k_1.size() - minl; i++){
        string cur_lmer = sufix_k_1.substr(i, minl);

        if (cur_lmer < min_rmer){
            min_rmer = cur_lmer;
        }
    }
    return min_rmer;
}

/*
 * 函数功能：split_to_buckets的辅助函数，桶达到一定大小之后，将桶写入文件中
 */
void CDBG::write_buckets_to_file() {
    auto bucket_iter = buckets.begin();
    while(bucket_iter != buckets.end()) {
        string out_file = "buckets/" + bucket_iter->first + ".txt";
        string content = bucket_iter->second;
        ofstream fout;
        if (access(out_file.c_str(), 0) == 0){
            fout.open(out_file, ios::app);
        } else {
            fout.open(out_file);
        }

        if (fout) { // 如果创建成功
            int cnt = 0;
            for (int i = 0; i < content.size(); i++) {
                cnt ++;
                fout << content[i];
                if (cnt % ksize == 0){
                    fout << endl;
                    cnt = 0;
                }
            }
        }
        fout.close();
        bucket_iter++;
    }
    buckets.erase( buckets.begin(), buckets.end() );
}

/*
 * 函数功能：split_to_buckets的辅助函数，将kmer写到对应的桶中
 */
void CDBG::write_to_buckets(string out_file, string kmer){

    auto bucket_iter = buckets.find(out_file);
    if(bucket_iter != buckets.end()){
        buckets[out_file] += kmer;
    }
    else{
        buckets[out_file] = kmer;
    }

}

/*
 * 函数功能：读取bam文件，根据lmm和rmm将kmer分入不同的桶（文件中）
 */
void CDBG::split_to_buckets() {
    // 创建桶保存的文件夹
    string save_folder = "buckets";
    // 文件夹存在
    if (access(save_folder.c_str(), 0) == 0){
        system("rm -rf buckets");
    }
    mkdir(save_folder.c_str(), 0777);

    uint64_t cnt = 0;

    bam1_t *aln = bam_init1(); //initialize an alignment
    bam_hdr_t *bamHdr = sam_hdr_read(bam_file); //read header
    while(sam_read1(bam_file, bamHdr, aln) > 0){
        //获得序列
        int32_t len = aln->core.l_qseq;
        uint8_t *q = bam_get_seq(aln);
        string read;
        for(int i=0; i< len ; i++){
            read += seq_nt16_str[bam_seqi(q, i)];
        }


        cnt ++;
        cout << cnt << endl;

        if (cnt % 100 == 0){
            // 桶达到一定大小之后才写入文件存储，提高效率
            write_buckets_to_file();
        }
        if (cnt == 50002){
            break;
        }

        // 遍历每一个kmer
        for (uint32_t i = 0; i < read.size() - ksize + 1; i++){
            string kmer = read.substr(i, ksize);

            string kmer_lmm = lmm(kmer, lsize);
            string kmer_rmm = rmm(kmer, lsize);

            // 写入桶中
            write_to_buckets(kmer_lmm, kmer);

            if (kmer_lmm != kmer_rmm){
                // 写入桶中
                write_to_buckets(kmer_rmm, kmer);
            }

        }

    }
}


/*
 * 函数功能：对于每一个桶文件，开启多线程，开始连接DBG和压缩、标记lonely等后续功能
 */
void CDBG::construct_dbg_with_thread(){

    ThreadPool thread_pool(thread_num);

    std::vector< std::future<void > > results;

    string save_folder = "buckets_dbg";
    if (access(save_folder.c_str(), 0) == 0){
        system("rm -rf buckets_dbg");
    }
    mkdir(save_folder.c_str(), 0777);


    string reunite_folder = "reunite";
    // 文件夹存在
    if (access(reunite_folder.c_str(), 0) == 0){
        system("rm -rf reunite");
    }
    mkdir(reunite_folder.c_str(), 0777);



    vector<string> files = get_files("buckets");
    int c = 0;

    for (auto vec_iter = files.cbegin(); vec_iter != files.cend(); vec_iter++) {

        c ++;
        cout << c << endl;

        string file_name = *vec_iter;

        results.emplace_back(thread_pool.enqueue(construct_dbg_in_bucket, this, file_name));

//        construct_dbg_in_bucket(this, file_name);
    }

}


/*
 * 在每一个桶中针对每一个kmer开始生成dbg，并且启动压缩步骤、标记lonely步骤
 */
void CDBG::construct_dbg_in_bucket(CDBG *p, string file_name) {

    CDBG cdbg1(*p);

    string file_path = "buckets/" + file_name;
    //每一个unitig的编号，从0开始计算
    long int out_file_cnt = 0;

    string kmer;

    // 存储生成的dbg，键值对分别是编号和节点
    unordered_map<uint64_t, dbg_node> bucket_dbg;

    ifstream file(file_path);

    // 取出所有kmer
    vector<string> all_kmers;

    while(getline(file, kmer)) {

        // 去掉重复的kmer
        auto vec_iter = find(all_kmers.begin(), all_kmers.end(), kmer);
        if(vec_iter == all_kmers.end()){
            // 加入到kmer集中
            all_kmers.emplace_back(kmer);
            // 根据当前kmer创建dbg节点
            vector<uint64_t> pre_node;
            vector<uint64_t> suf_node;
            dbg_node cur_node = {kmer, pre_node, suf_node, 0};
            // 加入到dbg的节点中
            bucket_dbg[out_file_cnt] = cur_node;

            out_file_cnt ++;
        }
    }

    // 遍历kmer，得到前后节点关系
    for (int i = 0; i < all_kmers.size(); i++) {

        string kmer1 = all_kmers[i];

        dbg_node* kmer1_node = &bucket_dbg[i];

        // 再次遍历这个kmer集，找到他的前面的和后面的节点
        for (int j = i ; j < all_kmers.size(); j++){
            string kmer2 = all_kmers[j];

            dbg_node* kmer2_node = &bucket_dbg[j];
            // 检查是否为邻近节点
            // 1是2的前一个节点
            if (cdbg1.sufix(kmer1, cdbg1.ksize - 1) == cdbg1.prefix(kmer2, cdbg1.ksize - 1)){
                kmer1_node->suf_node.emplace_back(j);
                kmer2_node->pre_node.emplace_back(i);
            } else if (cdbg1.prefix(kmer1, cdbg1.ksize - 1) == cdbg1.sufix(kmer2, cdbg1.ksize - 1)){
                kmer2_node->suf_node.emplace_back(i);
                kmer1_node->pre_node.emplace_back(j);

            }
        }
    }
    // 启动压缩步骤
    cdbg1.compact_dbg_in_bucket(bucket_dbg);

    // 启动标记lonely步骤
    cdbg1.set_lonely_in_bucket(bucket_dbg, file_name);

    cdbg1.write_dbg_to_file(bucket_dbg, "buckets_dbg/" + file_name);

}

/*
 * 在每一个桶中，进行压缩每一个节点的步骤
 */
void CDBG::compact_dbg_in_bucket(unordered_map<uint64_t, dbg_node> &bucket_dbg) {


    while (true){
        int stop_flag = 1;
        // 对于每一个节点，都作为unitig的起始节点，开始向下寻找
        for (auto cur_node = bucket_dbg.begin(); cur_node != bucket_dbg.end(); cur_node ++){

            vector<uint64_t> pre_node = cur_node->second.pre_node;
            vector<uint64_t> suf_node = cur_node->second.suf_node;

            uint64_t cur_num = cur_node->first;

            uint32_t in_degree = pre_node.size();
            uint32_t out_degree = suf_node.size();

            if (cur_node->second.flag != 1){
                // 节点的入度=1
                if (in_degree <= 1) {
                    // 出度也等于1
                    if (out_degree == 1) {
                        stop_flag = 0;

                        dbg_node new_node{cur_node->second.kmer, cur_node->second.pre_node, cur_node->second.suf_node, 0};
                        uint64_t next_num = cur_node->second.suf_node[0];
                        dbg_node next_node = bucket_dbg[next_num];

                        find_unitig_in_bucket(cur_num, next_node, next_num, new_node, bucket_dbg, cur_num);
                        cur_node->second.flag ++;

                        if (cur_node->second.kmer != new_node.kmer){
                            bucket_dbg[cur_num] = new_node;
                        }

                    }
                        // 出度大于1,则遇到了分支节点， 出度=0，没有后续节点了
                    else {
                        cur_node->second.flag ++;
                        continue;
                    }
                }
                else if (in_degree > 1){
                    if (out_degree == 1){
                        stop_flag = 0;

                        dbg_node new_node{cur_node->second.kmer, cur_node->second.pre_node, cur_node->second.suf_node, 0};
                        uint64_t next_num = cur_node->second.suf_node[0];
                        dbg_node next_node = bucket_dbg[next_num];

                        // 访问次数加一
                        cur_node->second.flag ++;
                        find_unitig_in_bucket(cur_num, next_node, next_num, new_node, bucket_dbg, cur_num);

                        if (cur_node->second.kmer != new_node.kmer){
                            bucket_dbg[cur_num] = new_node;
                        }

                    }
                        // 入度出度都大于1,则遇到了分支节点， 出度=0，没有后续节点了
                    else {
                        cur_node->second.flag ++;
                        continue;
                    }
                }
            }
        }

        if (stop_flag == 1){
            break;
        }
    }
}

/*
 * 函数功能：压缩函数的辅助函数，一个递归函数，对于当前节点进行递归，直到不能再合并
 */
// 建立在当前节点可以向下合并的基础上，参数的cur_node 就是当前节点的下一个节点。
void CDBG::find_unitig_in_bucket(uint64_t last_num, CDBG::dbg_node &cur_node, uint64_t cur_num,
                                 CDBG::dbg_node &new_node, unordered_map<uint64_t, CDBG::dbg_node> &bucket_dbg, uint64_t new_num) {

    vector<uint64_t> cur_pre_node = cur_node.pre_node;
    vector<uint64_t> cur_suf_node = cur_node.suf_node;

    uint32_t in_degree = cur_pre_node.size();
    uint32_t out_degree = cur_suf_node.size();
//    cout << "----" << in_degree << " " << out_degree << endl;
    // 节点的入度=1
    if (in_degree == 1){
        // 出度也等于1
        if (out_degree == 1){
            // 合并并更新节点
            new_node.kmer = compact_two_nodes(new_node.kmer, cur_node.kmer);
            new_node.suf_node = cur_node.suf_node;

            // 找到下一个节点
            uint64_t next_num = cur_suf_node[0];
            dbg_node next_node = bucket_dbg[next_num];
            // 该节点访问次数加一
            cur_node.flag ++;

            // 合并节点后，更前前后缀列标的信息，因为是中间节点，只需要更新下一个节点的信息
            erase_and_update_down(cur_node, cur_num, new_num, bucket_dbg);

            // 递归
            find_unitig_in_bucket(cur_num, next_node, next_num, new_node, bucket_dbg, new_num);

        }
            // 出度大于1,则遇到了分支节点， 出度=0，没有后续节点了
        else{
            new_node.kmer = compact_two_nodes(new_node.kmer, cur_node.kmer);
            new_node.suf_node = cur_node.suf_node;
            // 当前节点访问次数加一
            cur_node.flag ++;
            erase_and_update_down(cur_node, cur_num, new_num, bucket_dbg);

            return;
        }

    }
        // 入度大于1
    else if (in_degree > 1){
        // 出度也大于1,遇到不可合并的节点
        if (out_degree > 1){
            return;
        }
        else if (out_degree == 1){
            return;
        }
    }


}


/*
 * 函数功能：合并两个节点
 */
string CDBG::compact_two_nodes(string u, string v) {

    string suf_u = sufix(u, ksize - 1);
    string pre_v = prefix(v, ksize - 1);

    if (suf_u == pre_v){
        string w = u + v.substr(ksize - 1, v.size() - ksize + 1);
        return w;
    }
    return "";
}

/*
 * 函数功能：针对DBG中当前节点合并后更新后一个节点的前缀信息并删除
 */
void CDBG::erase_and_update_down(CDBG::dbg_node cur_node, uint64_t cur_num, uint64_t new_num,
                                 unordered_map<uint64_t, CDBG::dbg_node> &bucket_dbg) {

    vector<uint64_t> suf_node = cur_node.suf_node;

    // 当前节点的所有后缀节点
    for (int i = 0; i < suf_node.size(); i++){
        // 找到每一个后缀节点的前缀节点，
        vector<uint64_t> suf_pre_node = bucket_dbg[suf_node[i]].pre_node;
        // 更前前缀节点的num变化
        auto vec_iter = find(suf_pre_node.begin(), suf_pre_node.end(), cur_num);
        bucket_dbg[suf_node[i]].pre_node[distance(suf_pre_node.begin(), vec_iter)] = new_num;

    }
    bucket_dbg.erase(cur_num);
}


void CDBG::erase_and_update_up_and_down(CDBG::dbg_node cur_node, uint64_t cur_num, uint64_t new_num,
                                        unordered_map<uint64_t, CDBG::dbg_node> &bucket_dbg) {

    vector<uint64_t> pre_node = cur_node.pre_node;
    vector<uint64_t> suf_node = cur_node.suf_node;
    vector<uint64_t>::iterator vec_iter;
    for (int i = 0; i < pre_node.size(); i++){
        vector<uint64_t> pre_suf_node = bucket_dbg[pre_node[i]].suf_node;

        vec_iter = find(pre_suf_node.begin(), pre_suf_node.end(), cur_num);
        bucket_dbg[pre_node[i]].suf_node[distance(pre_suf_node.begin(), vec_iter)] = new_num;
    }
    // 当前节点的所有后缀节点
    for (int i = 0; i < suf_node.size(); i++){
        // 找到每一个后缀节点的前缀节点，
        vector<uint64_t> suf_pre_node = bucket_dbg[suf_node[i]].pre_node;
        // 更前前缀节点的num变化
        vec_iter = find(suf_pre_node.begin(), suf_pre_node.end(), cur_num);
        bucket_dbg[suf_node[i]].pre_node[distance(suf_pre_node.begin(), vec_iter)] = new_num;

    }
}


void CDBG::set_lonely_in_bucket(unordered_map<uint64_t, CDBG::dbg_node> &bucket_dbg, string file_name) {
    string i = file_name.substr(0, lsize);

    auto dbg_iter = bucket_dbg.begin();

    while(dbg_iter != bucket_dbg.end()) {

        string kmer = dbg_iter->second.kmer;

        string kmer_lmm = lmm(kmer, lsize);
        string kmer_rmm = rmm(kmer, lsize);


        lone_node new_node{"", "", kmer};

        // lonely的话标记为1
        if (kmer_lmm != i){

            dbg_iter->second.left_lonely = 1;

            string lone_left = kmer.substr(0, ksize);

            // 查找右lone表
            right_lone_unitigs_lock.lock();
            // 没有找到
            auto unitig_iter = right_lone_unitigs.find(lone_left);
            if (unitig_iter == right_lone_unitigs.end()){
                new_node.left_lone = lone_left;
            } else {

                // 找到，则合并
                lone_node right_lone_unitig = unitig_iter->second;
                new_node = glue(right_lone_unitig, new_node);
                right_lone_unitigs.erase(lone_left);
            }
            right_lone_unitigs_lock.unlock();

        }

        if (kmer_rmm != i){

            dbg_iter->second.right_lonely = 1;

            string lone_right = kmer.substr(kmer.size() - ksize, ksize);

            left_lone_unitigs_lock.lock();
            auto unitig_iter = left_lone_unitigs.find(lone_right);
            if (unitig_iter == left_lone_unitigs.end()){
                new_node.right_lone = lone_right;
            } else {
                lone_node left_lone_node = unitig_iter->second;
                new_node = glue(new_node, left_lone_node);
                left_lone_unitigs.erase(lone_right);
            }
            left_lone_unitigs_lock.unlock();

        }

        // 检查新node是否依然有lone的端，如果有，就加入到对应的unitig中去，没有，就写入output
        if (!new_node.left_lone.empty()){
            left_lone_unitigs_lock.lock();
            if (left_lone_unitigs.find(new_node.left_lone) != left_lone_unitigs.end()){
                cout << "already_have...." << endl;
            }
            left_lone_unitigs[new_node.left_lone] = new_node;
            left_lone_unitigs_lock.unlock();

        } else if (!new_node.right_lone.empty()){

            right_lone_unitigs_lock.lock();
            if (right_lone_unitigs.find(new_node.right_lone) != right_lone_unitigs.end()){
                cout << "already_have...." << endl;
            }

            right_lone_unitigs[new_node.right_lone] = new_node;
            right_lone_unitigs_lock.unlock();

        } else {
            all_unitig_lock.lock();
            string out_file = "reunite/all_unitig.txt";
            ofstream fout(out_file, ios::app);
            fout << new_node.kmer << endl;

            fout.close();

            all_unitig_lock.unlock();
        }

        dbg_iter ++;
    }


}


void CDBG::deal_left_lone() {
    string out_file = "reunite/all_unitig.txt";
    ofstream all_file(out_file, ios::app);

    int cnt = 0;

//    cout << "-----------" << endl;
    // 将最后剩下的lone——node进一步处理写入文件
    while (true){

        if (left_lone_unitigs.empty()){
            break;
        }
        auto left_k = left_lone_unitigs.begin();
        auto right_k = right_lone_unitigs.find(left_k->first);
        if (right_k == right_lone_unitigs.end()){
            // 虽然lonely找不到可以glue的node了，此时写入输出文件
            all_file << left_k->second.kmer << endl;
            left_lone_unitigs.erase(left_k);
        } else {
            // 找到则合并
            lone_node new_node = glue(right_k->second, left_k->second);
            left_lone_unitigs.erase(left_k);
            right_lone_unitigs.erase(right_k);
            // 判断合并后的是否仍然有lonely
            if (!new_node.left_lone.empty()){
                left_lone_unitigs[new_node.left_lone] = new_node;
            } else if (!new_node.right_lone.empty()){
                right_lone_unitigs[new_node.right_lone] = new_node;
            } else {
                all_file << new_node.kmer << endl;
            }

        }

    }
    // 将剩余的右边lonely的写入输出文件
    for (auto right_k = right_lone_unitigs.begin(); right_k != right_lone_unitigs.end(); right_k ++){
        cnt ++;
        all_file << right_k->second.kmer << endl;
    }

    all_file.close();
//    cout << cnt << endl;
}
CDBG::lone_node CDBG::glue(CDBG::lone_node u, CDBG::lone_node v) {

    string kmer_u = u.kmer;
    string kmer_v = v.kmer;

    string new_kmer = kmer_u + kmer_v.substr(ksize, kmer_v.size() - ksize);

    lone_node new_node{u.left_lone, v.right_lone, new_kmer};

    return new_node;

}

void CDBG::construct_final_dbg() {
    ifstream file("reunite/all_unitig.txt");
    vector<string>::iterator vec_iter;
    vector<string> all_kmers;
    string kmer;
    unordered_map<uint64_t, dbg_node> bucket_dbg;
    long int out_file_cnt = 0;


    while(getline(file, kmer)) {

        // 去掉重复的kmer
        vec_iter = find(all_kmers.begin(), all_kmers.end(), kmer);
        if(vec_iter == all_kmers.end()){
            // 加入到kmer集中
            all_kmers.emplace_back(kmer);

            vector<uint64_t> pre_node;
            vector<uint64_t> suf_node;

            dbg_node cur_node = {kmer, pre_node, suf_node, 0};
            // 加入到dbg的节点中
            bucket_dbg[out_file_cnt] = cur_node;
            out_file_cnt ++;
        }
    }

    // 遍历kmer，得到前后节点关系
    for (int i = 0; i < all_kmers.size(); i++) {
        cout << i << endl;
        string kmer1 = all_kmers[i];


        dbg_node* kmer1_node = &bucket_dbg[i];

        // 再次遍历这个kmer集，找到他的前面的和后面的节点
        for (int j = i ; j < all_kmers.size(); j++){
            string kmer2 = all_kmers[j];

            dbg_node* kmer2_node = &bucket_dbg[j];
            // 检查是否为邻近节点
            // 1是2的前一个节点
            if (sufix(kmer1, ksize - 1) == prefix(kmer2, ksize - 1)){
                kmer1_node->suf_node.emplace_back(j);
                kmer2_node->pre_node.emplace_back(i);
            } else if (prefix(kmer1, ksize - 1) == sufix(kmer2, ksize - 1)){
                kmer2_node->suf_node.emplace_back(i);
                kmer1_node->pre_node.emplace_back(j);

            }
        }
    }
    cout << "writing" << endl;

    write_dbg_to_file(bucket_dbg, "reunite/final_dbg.txt");
}

/*
 * 把dbg写入文件中
 */
void CDBG::write_dbg_to_file(unordered_map<uint64_t , dbg_node> &bucket_dbg, string file_name) {

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

void start_method(){

    CDBG cdbg("../../bams/sampleChr20_sorted.bam", 16, 8, 6);

    cdbg.split_to_buckets();
    cdbg.construct_dbg_with_thread();
    cdbg.deal_left_lone();
    cdbg.construct_final_dbg();
}


int main(){
    //    g++ CDBG.cpp -lhts -pthread -o cdbg


    time_t t = time(nullptr);

//    start_method();

//    string seq = "ATCGT";
//    unsigned k = 5;
//
//    string kmer = seq.substr(0, k);
//    uint64_t hVal=0;
//    hVal = NTF64(kmer.c_str(), k); // initial hash value
//
//    for (size_t i = 0; i < seq.length() - k; i++)
//    {
//        hVal = NTF64(hVal, seq[i], seq[i+k], k); // consecutive hash values
//    }
//
//    cout << hVal << endl;



//    RMtoDBG rMtoDBG();
//    rMtoDBG().get_from_file();
    showMemStat(getpid());
    double cost_t = time(nullptr) - t;
    cout << "all time : " << cost_t  << "s" << endl;

}

