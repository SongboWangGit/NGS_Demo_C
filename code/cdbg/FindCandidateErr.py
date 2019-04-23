
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

    def check_and_merge(self, err_pos, read, read_pos):
        flag = 1
        err_pos_on_ref = err_pos + read_pos
        # 找是否在已有范围内
        for key in self.err_pos_dict.keys():
            if err_pos_on_ref >= key - 20 and err_pos_on_ref <= key + 20:
                start_cut = err_pos - (20 - (key - err_pos_on_ref))
                end_cut = err_pos + (40 - (20 - (key - err_pos_on_ref)))
                if start_cut < 0:
                    flag = 0
                    break
                err_segment = read[start_cut: end_cut]
                if err_segment not in self.err_pos_dict[key]:
                    self.err_pos_dict[key].append(err_segment)
                flag = 0
                break
        # 不在已有范围内
        if flag == 1:
            self.err_pos_dict[err_pos_on_ref] = []
            err_segment = read[err_pos - 20: err_pos + 20]
            self.err_pos_dict[err_pos_on_ref].append(err_segment)



    def write_err_to_file(self):
        out_file = open("err_info.txt", 'w')
        for item in self.err_pos_dict.keys():
            out_file.write(str(item))
            out_file.write(" ")
            for s in self.err_pos_dict[item]:
                out_file.write(s)
                out_file.write(" ")
            out_file.write("\n")

        out_file.close()

    def trim_read_with_cigar(self, read, qual, cigar):
        for (type, num) in cigar:
            if type == 4:
                read = read[num:]
                qual = qual[num:]
                # read_pos += num

        return read, qual
    def find_err_pos(self):
        cnt = 0
        for r in self.bam_file:

            cnt += 1
            if cnt == 27000:
                break
            rname = r.qname
            rflag = r.flag
            # 去掉read2
            # if self.judge_read1_or_2(rflag) == 1:
            # 获得bp串
            read = r.query_sequence
            if len(read) < 150:
                continue
            read_pos = r.pos
            qual = list(r.query_qualities)

            # 4:S, 0:M
            cigar = r.cigartuples
            if cigar == None:
                continue
            # print(cigar)
            print(read_pos)
            read, qual = self.trim_read_with_cigar(read, qual, cigar)



            for i in range(20, len(qual)- 20):
                if qual[i] < 27:
                    self.check_and_merge(i, read, read_pos)

        self.write_err_to_file()

            # break

if __name__ == '__main__':
    in_file = "../../bams/sampleChr20_sorted.bam"
    find_err_pos = FindCandidateErr(in_file)
    find_err_pos.find_err_pos()