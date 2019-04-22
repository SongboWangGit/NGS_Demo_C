
import pysam

class FindCandidateErr():
    def __init__(self, file):
        self.in_file = file
        self.bam_file = pysam.AlignmentFile(self.in_file, 'rb')

        self.err_pos_dict = {}

    def judge_read1_or_2(self, flag):
        binNum = bin(flag).replace('0b', '')
        binNumStr = ''.join(binNum)
        length = len(binNumStr)
        # print(binNum)
        # print(length)
        # print(binNumStr[-length])
        if length == 7:
            return 1

        if binNumStr[-7] == '1':
            return 1
        elif binNumStr[-8] == '1':
            return 2

    def check_and_merge(self, err_pos):
        for i in range(err_pos - 40, err_pos + 40):
            if i in self.err_pos_dict.keys():
                pass

    def find_err_pos(self):
        for r in self.bam_file:
            rname = r.qname
            rflag = r.flag

            # 去掉read2
            # if self.judge_read1_or_2(rflag) == 1:
                # 获得bp串
            read = r.query_sequence
            ref_pos = r.pos


            # 获得数值质量
            qual = list(r.query_qualities)
            for i in range(len(qual)):
                if qual[i] < 27:
                    self.check_and_merge(ref_pos + i)

            print(ref_pos)
            break

if __name__ == '__main__':
    in_file = "../../bams/sampleChr20_sorted.bam"
    find_err_pos = FindCandidateErr(in_file)
    find_err_pos.find_err_pos()