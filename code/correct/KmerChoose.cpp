//
// Created by sbwang on 19-4-16.
//

#include "KmerChoose.h"

void KmerChoose::kmer_stat_start(string in_file, int ksize, int s, int r){
    KmerStat kmerStat(in_file, ksize, s, r);
    kmerStat.Update();
    kmerStat.Estimate();
}
void KmerChoose::kmer_stat_thread(string in_file, size_t thread_num, int s, int r){
    ThreadPool thread_pool(thread_num);
    for (int ksize = 21; ksize <=101; ksize += 10){
        thread_pool.enqueue(kmer_stat_start, in_file, ksize, s, r);
    }
}

uint64_t KmerChoose::choose_ksize(){
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

    KmerChoose kmerChoose;

//    kmer_stat_thread("../../bams/sampleChr20_sorted.bam", 8, 8, 20);
    kmerChoose.choose_ksize();
    showMemStat(getpid());
    double cost_t = time(nullptr) - t;
    cout << "all time : " << cost_t  << "s" << endl;

    return 0;
}
