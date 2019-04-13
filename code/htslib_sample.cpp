#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "htslib/sam.h"
#include <iostream>
#include <string>
using namespace std;

//g++ main.cpp -L /home/sbwang/user/Micro/htslib-1.9 -lhts -o test

int main(){



    samFile *fp_in = hts_open("../bams/sampleChr20_sorted.bam", "r"); //open bam file
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment

    //头部信息
    char *tar = bamHdr->text ;
    uint32_t *tarlen = bamHdr->target_len ;

    int cnt = 0;
    while(sam_read1(fp_in,bamHdr,aln) > 0){

        cnt++;
        if (cnt == 10){
            break;
        }
        // read起始位置
        int32_t pos = aln->core.pos +1;
        // 染色体信息
        char *chr = bamHdr->target_name[aln->core.tid] ;
        // 序列的长度
        int32_t len = aln->core.l_qseq;

        // 序列名称
        char* qname  = bam_get_qname(aln);
        // CIGAR
        uint32_t *cigar = bam_get_cigar(aln);

        //比对质量
        uint32_t q2 = aln->core.qual ;

        //获得序列
        uint8_t *q = bam_get_seq(aln);

        char *read = (char *)malloc(len);
        for(int i=0; i< 150 ; i++){
            read[i] = seq_nt16_str[bam_seqi(q, i)];
        }

        uint8_t *qs = bam_get_qual(aln);

        cout << read << endl;


    }

    bam_destroy1(aln);
    sam_close(fp_in);

    return 0;

}