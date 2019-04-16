//
// Created by sbwang on 19-4-12.
//

#include "KmerStat.h"



KmerStat::KmerStat(string file, uint16_t k, uint16_t s_num, uint16_t r_num) {

    string save_folder = "kmer_count";
    // 文件夹不存在存在
    if (access(save_folder.c_str(), 0) != 0){
        mkdir(save_folder.c_str(), 0777);
    }

    bam_file = hts_open(file.c_str(), "r");
    ksize = k;

    r = r_num;
    s = s_num;

    table_size = uint64_t (pow(2, r));
    count_table = new uint32_t[table_size] ();

    max_count = 0;

}


uint64_t KmerStat::ten_2_tow(uint64_t n){

    uint64_t ans;
    stack<uint64_t> stk;
    while(n!=0){
        ans=n%2;
        n=n/2;
        stk.push(ans);
    }
    while(!stk.empty()){
        cout<<stk.top();
        stk.pop();
    }
    cout << endl;

}


void KmerStat::write_table_to_file() {
    ofstream fout("kmer_count/kmer_table.txt");

    for (int i = 0; i < table_size; i++){
        fout << count_table[i] << endl;
    }
    fout.close();

}


void KmerStat::write_f_to_file(double *f, uint64_t F0, uint64_t F1) {

    string out_file = "kmer_count/" + to_string(ksize) + "_f_table.txt";
    ofstream fout(out_file);
    fout << std::fixed;

    // 写入F0：总kmer种类;F1：总kmer数量
    fout << ksize << " " << F0 << " " << F1 << endl;

    uint64_t cout_size;
    if (max_count > 500){
        cout_size = 500;
    } else {
        cout_size = max_count;
    }

    for (int i = 1; i < cout_size; i++){
        if (uint64_t(f[i]) < 1844674407370954789){
            fout << i << " " << uint64_t(f[i]) << endl;
        }
    }
    fout.close();
}

uint64_t KmerStat::get_nthash(string read) {
//    string kmer = read.substr(0, ksize);
//    uint64_t hVal=0;
//    hVal = NTF64(kmer.c_str(), ksize); // initial hash value
//    judge_and_set(hVal);
//
//    for (size_t i = 0; i < read.length() - ksize; i++)
//    {
//        hVal = NTF64(hVal, read[i], read[i+ksize], ksize); // consecutive hash values
//        judge_and_set(hVal);
//    }

    /* k is the k-mer length */
    unsigned k = ksize;

    /* h is the number of hashes for each k-mer */
    unsigned h = 1;

    /* init ntHash state and compute hash values for first k-mer */
    ntHashIterator itr(read, h, k);
    while (itr != itr.end()) {
        judge_and_set((*itr)[0]);
        ++itr;
    }

}

void KmerStat::judge_and_set(uint64_t hash_num) {

    uint64_t end_r_bits = hash_num & uint64_t(pow(2, r) - 1);
    uint64_t start_s_bits = (hash_num >> (64 - r)) & uint64_t(pow(2, s) - 1);

    if (start_s_bits == 0){
        count_table[end_r_bits] ++;
        // 找最大值
        if (count_table[end_r_bits] > max_count){
            max_count = count_table[end_r_bits];
        }
    }


}
void KmerStat::Update() {

    uint64_t cnt = 0;

    bam1_t *aln = bam_init1(); //initialize an alignment
    bam_hdr_t *bamHdr = sam_hdr_read(bam_file); //read header
    while(sam_read1(bam_file, bamHdr, aln) > 0) {
        //获得序列
        int32_t len = aln->core.l_qseq;
        uint8_t *q = bam_get_seq(aln);
        string read;
        for (int i = 0; i < len; i++) {
            read += seq_nt16_str[bam_seqi(q, i)];
        }

        cnt++;
        cout << cnt << endl;

//        if (cnt == 1200000){
//            break;
//        }
        if (read.size() != 150){
            continue;
        }

        get_nthash(read);
    }

//    write_table_to_file();

}

double* KmerStat::Estimate() {

    auto *p = new double[max_count]();
    auto *f = new double[max_count];

    for (int i = 0; i < table_size; i++){
        p[count_table[i]] ++;
    }


    for (int i = 0; i < max_count; i++){
        p[i] = p[i] / pow(2, r);
    }

    // 估计F0
    double F0 = -(log(p[0]) * pow(2, s + r));
    cout << p[0] << endl;

    for (int i = 0; i < max_count; i++){
        f[i] = -(p[i] / (p[0] * log(p[0])));

        for (int j = 0; j < i - 1; j++){
            f[i] = f[i] - ((j * p[i - j] * f[j]) / p[0]) / i;
        }
    }

    for (int i = 0; i < max_count; i++){
        f[i] = f[i] * F0;
    }
    uint64_t F1 = cal_F1(f);

    write_f_to_file(f, uint64_t(F0), F1);
    return f;

}

uint64_t KmerStat::cal_F1(double *f) {

    double F1 = 0;

    for (int i = 1; i < max_count; i++){
        F1 += f[i];
    }
    return uint64_t(F1);
}

void kmer_stat_start(string in_file, int ksize, int s, int r){
    KmerStat kmerStat(in_file, ksize, s, r);
    kmerStat.Update();
    kmerStat.Estimate();
}
void kmer_stat_thread(string in_file, size_t thread_num, int s, int r){
    ThreadPool thread_pool(thread_num);
    for (int ksize = 21; ksize <=101; ksize += 10){
        thread_pool.enqueue(kmer_stat_start, in_file, ksize, s, r);
    }
}

uint64_t choose_ksize(){
    vector<vector<uint64_t>> all_hists;


    for (int ksize = 21; ksize <=81; ksize += 10){
        string table_file_name = "kmer_count/" + to_string(ksize) +"_f_table.txt";


        ifstream table_file(table_file_name);
        if (not table_file.is_open()){
            cout << "No such file" << endl;
            break;
        }
        string line01;
        getline(table_file, line01);


        string content_line;
        vector<uint64_t> kmer_hist;
        while (not table_file.eof()){
            string frequency;
            string num;

            getline(table_file, frequency, ' ');
            getline(table_file, num);
            if (num.size() > 0){
                kmer_hist.emplace_back(stoll(num));
            }

        }
        all_hists.emplace_back(kmer_hist);
    }
    double frac = 1.3;
    // 要求的最大threshold
    int max_threshold = 2;
    vector<uint64_t> minimumList;
    bool flag = true;
    for (int i = 0; i < all_hists.size() and flag; i++){
        for (int j = 0; j < all_hists[i].size() and flag; j++){
            // 找到第一个曲线回升的点
            if (all_hists[i][j] < all_hists[i][j+1] * frac){
                // 如果该点大于设置的最大阈值
//                cout << all_hists[i][j] << " " << all_hists[i][j+1] << endl;
//                cout << "------------ " << j << endl;

                if (j >= max_threshold){
                    minimumList.emplace_back(j);
                } else {
                    flag = false;
                }

                break;
            }

        }
    }
    uint64_t ksize_chosen = minimumList.size()*10 + 21;
    cout <<  ksize_chosen << endl;
    return ksize_chosen;
}
int main(){

//g++ KmerStat.cpp -o KmerStat -lhts -lpthread
    time_t t = time(nullptr);


//    kmer_stat_thread("../../bams/sampleChr20_sorted.bam", 8, 8, 20);
    choose_ksize();
    showMemStat(getpid());
    double cost_t = time(nullptr) - t;
    cout << "all time : " << cost_t  << "s" << endl;

    return 0;
}

