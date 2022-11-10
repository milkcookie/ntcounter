#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include <math.h>
#include "vafcaller.h"

void usage();
void ntcounts(const char *bam, const char *bedfile, uint32_t q, uint32_t F, const char *fafile, const char *op);
void printProgress(double percentage, int done, int tot);
int countline(const char *filename);

int main(int argc, char *argv[]) { 

    int opt;
    opterr = 0;
    char *bam, *bed, *fa, *op;
    while ((opt = getopt(argc, argv, "b:l:r:o:h")) != EOF)
        switch (opt) {
            case 'b': bam = optarg;
                break;
            case 'l': bed = optarg;
                break;
            case 'r': fa = optarg;
                break;
            case 'o': op = optarg;
                break;
            case 'h':
                usage();
                exit(1);
                break;
            default:
                usage();
                exit(1);
                break;
            case '?':
                usage();
                break;
        }
    if (argc == 0) {
        usage();
        exit(1);
    }
    uint32_t q = 10;
    uint32_t F = 1024;
    ntcounts(bam, bed, q, F, fa, op);
    return 0;
}

void usage() {
    printf("usage: vafcaller [options] -b bamfile -l locifile -r ref_fasta -o output file....\n");
    printf("-b the input file in bam format....\n");
    printf("-l bedfile, the input file , just with chrom and pos two columns....\n");
    printf("-r reference fasta file in fasta format....\n");
    printf("-o output file the result file....\n");
        // char *bam = "/PROJ/Proj_2022/20220906002MRD/SZL_WES_394434/Mapping/SZL-QP/SZL-QP.final.bam";
        // char *ref = "/Data/HUMAN_DB/HG19/ref_fa/human_hg19.fa";
        // char *maf = "/PROJ/Proj_2022/20220906002MRD/SZL_WES_394434/Somatic/SZL-QP/MAF/Somatic_mutation.maf";
        // char *result = "result";
        // uint32_t q = 10;
        // uint32_t F = 1024;
        // ntcounts(bam, maf, q, F, ref, result);
    }


int countline(const char *filename) {
    //计算bed文件行数
    FILE *fp = fopen(filename, "r");
    int ch = 0;
    int lines = 0;
    if (fp == NULL) {
        return 0;
    }

    while (!feof(fp)) {
        ch = fgetc(fp);
        if (ch == '\n') {
            lines ++;
        }
    }
    fclose(fp);
    return lines;
}

void printProgress(double percentage, int done, int tot) {
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  fprintf(stderr, "\r%3d%% [%.*s%*s] %d/%d", val, lpad, PBSTR, rpad, "", done, tot);
  fflush(stdout);
}

void ntcounts(const char *bam, const char *bedfile, uint32_t q, uint32_t F, const char *fafile, const char *op) {
    int vars_gt = 1;
    hts_verbose = 0;

    char tsv_file[1000];
    strcpy(tsv_file, op);
    strcat(tsv_file, ".tsv");
    char maf_file_out[1000];
    strcpy(maf_file_out, op);
    strcat(maf_file_out, ".maf");
    int nloci = countline(bedfile);

    FILE *bed_fp;
    bed_fp = fopen(bedfile, "r");
    char buff[1000];
    FILE *tsv_fp;
    tsv_fp = fopen(tsv_file, "w" );

    FILE *maf_fp;
    maf_fp = fopen(maf_file_out, "w");

    char *seq;
    faidx_t *fa = fai_load(fafile);  //fai_load 是faidx.h 里面导入的
    samFile *fp_in = hts_open(bam, "r");
    hts_idx_t *fp_idx = sam_index_load(fp_in, bam);
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    bam1_t *aln = bam_init1(); //返回一个空的 bam结构的对象 最后需要销毁（bam_destroy1) 主要用于传递给 sam_itr_next

    uint64_t n_mapped = 0;
    uint64_t n_unmapped = 0;
    uint64_t tot_mapped = 0;
    int res = 0;
    int32_t n_contigs = bamHdr->n_targets;

    for (int i=0; i < n_contigs; i++) {
        res = hts_idx_get_stat(fp_idx, i, &n_mapped, &n_unmapped); 
        if (res == 0) {
            //printf("%d\n", tot_mapped);
            tot_mapped = tot_mapped + n_mapped;
        }
    }
    fprintf(tsv_fp, "#idxstats_mapped_reads\t%d\n", tot_mapped);
    fprintf(tsv_fp, "loci\tfa_ref\tA\tT\tG\tC\tIns\tDel\n");
    char *snp_type = "SNP";
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    while (read = getline(&line, &len, bed_fp) != -1){ //按行读取文件
        int len = read;
        //printf("%s", buff);
        //printf("%s", buff_str);
        char *p_buff = line;
        int strlength = strlen(line);
        //printf("%d\n", strlength);
        char line_str[1000];
        sprintf(line_str, "%s", line);
        char *token, *chrom, ref_base, alt_base, *pos_char, *snvtype;
        double vaf;
        for (int i=0;i<37;i++) {
            if (i==4){
                chrom = strsep(&p_buff, "\t");
                if (strcmp(chrom, "Chromosome")==0) {
                    //printf("wtf!\n");
                    fputs(line_str, maf_fp);
                    break;
                }
                //printf("chrom=%s\n", chrom);
            } else if (i == 5) {
                pos_char = strsep(&p_buff, "\t");
                //printf("pos_char=%s\n", pos_char);
            } else if (i == 9) {
                snvtype = strsep(&p_buff, "\t");
            } else if (i == 10) {
                ref_base = strsep(&p_buff, "\t")[0];
                //printf("ref_base=%c\n", ref_base);
            } else if (i == 12) {
                alt_base = strsep(&p_buff, "\t")[0];
                //printf("alt_base=%c\n", alt_base);
            } else if (i == 36) {
                char *vaf_char = strsep(&p_buff, "\t");
                sscanf(vaf_char, "%.3lf", &vaf);
                break;
            } else {
            strsep(&p_buff, "\t");
            }
        }
        // strtok(buff, "\t");
        // strtok(NULL, "\t");
        // strtok(NULL, "\t");
        // strtok(NULL, "\t");
        // char *chrom = strtok(NULL, "\t");
        // char *start = strtok(NULL, "\t");
        char loci[250] = "";
        strcat(loci, chrom); strcat(loci, ":"); strcat(loci, pos_char);strcat(loci, "-");strcat(loci, pos_char);
        //printf("%s\t", loci);
        if (fa != NULL) {
            int templen = 100;
            seq = fai_fetch(fa, loci, &templen); //获取fasta文件目标片段信息
            //fprintf(tsv_fp, "%s:%s\t%s", chrom, pos_char, seq);
            free(seq);
        } else {
            //fprintf(tsv_fp, "%s:%s\tNA", chrom, pos_char);
        }
        int32_t target_pos = atoi(pos_char) - 1; //输入的位置是 1-based atoi函数是将string 转化成int

        hts_itr_t *samitr = sam_itr_querys(fp_idx, bamHdr, loci);
        //构建一个多重区域的迭代器 参数 idx是index， bamHdr是bam Header信息， loci是array of ref： interval region specifiers

        int32_t tot_reads = 0;
        float nt[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        int32_t alt_base_num = 0;
        vars_gt = vars_gt + 1; 
        while (sam_itr_next(fp_in, samitr, aln) > 0) {
            int32_t pos = aln->core.pos; //获取pos 即reads比对的起始位置
            uint32_t len = aln->core.l_qseq; //获取length
            uint32_t *cig = bam_get_cigar(aln);// 获取指向cigar array得一个指针
            uint8_t *qs = bam_get_seq(aln); //获取 seq
            

            if (aln->core.qual <= q) { //过滤低质量的reads
                continue;
            }

            if (aln->core.flag >= F) { //过滤不符合flag的reads
                continue;
            } 
            
            tot_reads = tot_reads + 1;
            char *qseq = (char *)malloc(len); //动态申请内存
            int i = 0;

            for (i=0;i<len;i++) {
                qseq[i] = seq_nt16_str[bam_seqi(qs, i)];
                //seq_nt16_str 这个定义了flow order 对应的IUPAC编码
                //seq_nt16_str = "=ACMGRSVTWYHKDBN"
                //bam_seqi返回4-bit的一个数字， 代表对应碱基的位置 4-bit二进制最大为15， 即最多只能有index到16.
            }

            int32_t pos_onread = 0; //定义当前cigar在read上的位置，初始值为0

            int k = 0;

            for (k=0;k< aln->core.n_cigar; ++k) { // n_cigar 类型是 uint32_t 指的是cigar个数
                int cop = cig[k] & BAM_CIGAR_MASK; // BAM_CIGAR_MASK 是 0xf 十进制15 二进制1111 按位与 & 获取cigar string 只保留前四位 因为前四位存储的是cigar信息的位置
                int cl = cig[k] >> BAM_CIGAR_SHIFT; //  BAM_CIGAR_SHIFT 这个值是4 二进制右移四位。只保留五位数以上的位 获取cigar length
                //这里是因为 cigar的顺序是数字加类型只要去除右边四个bit就能去掉字母代表的类型 例如 19S78M1D53M
                char cigar_type = BAM_CIGAR_STR[cop];
                if (cigar_type == 'M') {
                    pos_onread = pos_onread + cl;
                    pos = pos + cl;
                } else if (cigar_type == 'S') {
                    pos_onread = pos_onread + cl;
                } else if (cigar_type == 'I') {
                    pos_onread = pos_onread + cl;
                } else if (cigar_type == 'D') {
                    pos = pos + cl;
                }
                
                if (pos > target_pos) {
                    if (cigar_type == 'M') {
                        pos_onread = pos_onread - (pos - target_pos);
                        char base = qseq[pos_onread];
                        if (base == alt_base) {
                            alt_base_num += 1;
                        }
                        //printf("base=%c\n", base);
                        if (base == 'A') {
                            nt[0] += 1;
                        } else if (base == 'T') {
                            nt[1] += 1;
                        } else if ( base == 'G') {
                            nt[2] += 1;
                        } else if ( base == 'C') {
                            nt[3] += 1;
                        } else {
                            continue;
                        }
                    }
                } else if (pos == target_pos) {
                    if (cigar_type == 'I') {
                        nt[4] = nt[4] + 1;
                        break;
                    } else if (cigar_type == 'D') {
                        nt[5] = nt[5] + 1;
                        break;
                    }
                }
            }
            free(qseq);//释放动态申请的内存
        }
        hts_itr_destroy(samitr);
        double bam_vaf = (double)alt_base_num / tot_reads;
        //printf("alt_base_num: %d\n", alt_base_num);
        //printf("%.3lf\n", bam_vaf);
        // printf("snvtype=%d\n", sizeof(snvtype));
        // printf("snp_typr=%d\n", sizeof(snp_type));
        if (bam_vaf > 0.05 && strncmp(snvtype, "SNP", 3) == 0){
            fprintf(tsv_fp, "%s:%s\t%c\t%.f\t%.f\t%.f\t%.f\t%.f\t%.f\t%.3lf\n",  chrom, pos_char, ref_base, nt[0], nt[1], nt[2], nt[3], nt[4], nt[5], bam_vaf);
            fputs(line_str, maf_fp);
        }
    }
    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    // fai_destroy(fa);
    sam_close(fp_in);
    fclose(tsv_fp);
    fclose(maf_fp);
    //释放内存 关闭数据

}
