import os
import importlib
import sys
import subprocess
import sprint
from multiprocessing import Pool
from tqdm import tqdm
import re


# ===================== 全局并行函数（序列化兼容）=====================
def get_chrom_list(zz_file):
    """提取所有染色体列表（从合并的zz文件中）"""
    chrom_set = set()
    with open(zz_file, 'r') as f:
        for line in f:
            if line.strip():
                chrom = line.split('\t')[0]
                chrom_set.add(chrom)
    return sorted(list(chrom_set))


def process_snv_by_chrom(args):
    """单染色体SNV识别（并行任务）"""
    chrom, tmp_dir, refgenome, has_trans, baseq_cutoff = args
    chrom_prefix = chrom.replace('/', '_').replace(':', '_')  # 避免路径非法字符
    
    # 输出文件路径
    genome_mskAG_snv = f"{tmp_dir}/genome_mskAG_all.zz.dedup.{chrom_prefix}.snv"
    genome_mskTC_snv = f"{tmp_dir}/genome_mskTC_all.zz.dedup.{chrom_prefix}.snv"
    genome_all_snv = f"{tmp_dir}/genome_all.zz.dedup.{chrom_prefix}.snv"
    
    # 按染色体过滤zz文件并处理SNV
    for zz_type in ['genome_mskAG_all.zz.dedup', 'genome_mskTC_all.zz.dedup', 'genome_all.zz.dedup']:
        zz_file = f"{tmp_dir}/{zz_type}"
        filtered_zz = f"{zz_file}.{chrom_prefix}"
        snv_out = f"{zz_file}.{chrom_prefix}.snv"
        
        # 过滤当前染色体的zz记录
        with open(zz_file, 'r') as in_f, open(filtered_zz, 'w') as out_f:
            for line in in_f:
                if line.strip() and line.split('\t')[0] == chrom:
                    out_f.write(line)
        
        # 调用SNV识别
        sprint.mask_zz2snv(filtered_zz, snv_out, baseq_cutoff)
        
        # 清理临时过滤文件
        if os.path.exists(filtered_zz):
            os.remove(filtered_zz)
    
    # 转录组相关SNV处理（如有）
    trans_snv_results = []
    if has_trans:
        trans_files = [
            'transcript_mskAG_all.zz.dedup',
            'transcript_mskTC_all.zz.dedup',
            'transcript_all.zz.dedup'
        ]
        trans_loc = f"{refgenome}.trans.fa.loc"
        
        for trans_zz in trans_files:
            zz_file = f"{tmp_dir}/{trans_zz}"
            filtered_zz = f"{zz_file}.{chrom_prefix}"
            snv_out = f"{zz_file}.{chrom_prefix}.snv"
            genome_snv = f"{snv_out}.genome.snv"
            
            # 过滤染色体
            with open(zz_file, 'r') as in_f, open(filtered_zz, 'w') as out_f:
                for line in in_f:
                    if line.strip() and line.split('\t')[0] == chrom:
                        out_f.write(line)
            
            # 转录组SNV处理
            sprint.mask_zz2snv(filtered_zz, snv_out, baseq_cutoff)
            sprint.tzz2gzz(trans_loc, filtered_zz, f"{filtered_zz}.genome.zz")
            sprint.transcript_locator(snv_out, trans_loc, genome_snv)
            sprint.transcript_sort(genome_snv, f"{genome_snv}.sort")
            
            trans_snv_results.append((f"{genome_snv}.sort", chrom_prefix))
            
            # 清理临时文件
            for tmp_f in [filtered_zz, f"{filtered_zz}.genome.zz", snv_out, genome_snv]:
                if os.path.exists(tmp_f):
                    os.remove(tmp_f)
    
    return {
        'chrom': chrom,
        'genome_mskAG_snv': genome_mskAG_snv,
        'genome_mskTC_snv': genome_mskTC_snv,
        'genome_all_snv': genome_all_snv,
        'trans_snv': trans_snv_results
    }


def process_res_by_chrom(args):
    """单染色体RES识别（并行任务）"""
    chrom, tmp_dir, repeat_file, cluster_params, has_trans = args
    chrom_prefix = chrom.replace('/', '_').replace(':', '_')
    
    # 解析聚类参数
    (cluster_distance, cluster_size_alu_ad1, cluster_size_alu_ad2,
     cluster_size_nalurp, cluster_size_nrp, cluster_size_rg,
     cluster_size_hp, cluster_size_alu_hp, cluster_size_nalurp_hp,
     cluster_size_nrp_hp) = cluster_params
    
    # 输出路径
    regular_res = f"{tmp_dir}/regular.res.{chrom_prefix}"
    hyper_mskTC_res = f"{tmp_dir}/hyper_mskTC.res.{chrom_prefix}"
    hyper_mskAG_res = f"{tmp_dir}/hyper_mskAG.res.{chrom_prefix}"
    
    # 过滤当前染色体的SNV文件
    def filter_snv(snv_in, snv_out):
        with open(snv_in, 'r') as in_f, open(snv_out, 'w') as out_f:
            for line in in_f:
                if line.strip() and line.split('\t')[0] == chrom:
                    out_f.write(line)
    
    # 1. 常规RES处理
    regular_snv = f"{tmp_dir}/regular.snv"
    filtered_regular_snv = f"{regular_snv}.{chrom_prefix}"
    filter_snv(regular_snv, filtered_regular_snv)
    
    if repeat_file:
        # 带重复序列注释的RES识别
        anno_out = f"{filtered_regular_snv}.anno"
        sprint.annotate(filtered_regular_snv, repeat_file, anno_out)
        sprint.seperate(anno_out, f"{anno_out}.alu", f"{anno_out}.nalurp", f"{anno_out}.nrp", 'Alu')
        
        # Alu AD2过滤
        sprint.get_snv_with_ad(f"{anno_out}.alu", f"{anno_out}.alu.ad2", 2)
        
        # 聚类
        sprint.snv_cluster(f"{anno_out}.alu", f"{regular_res}.alu.ad1", cluster_distance, cluster_size_alu_ad1)
        sprint.snv_cluster(f"{anno_out}.alu.ad2", f"{regular_res}.alu.ad2", cluster_distance, cluster_size_alu_ad2)
        sprint.bed_or(f"{regular_res}.alu.ad1", f"{regular_res}.alu.ad2", f"{regular_res}.alu")
        sprint.snv_cluster(f"{anno_out}.nalurp", f"{regular_res}.nalurp", cluster_distance, cluster_size_nalurp)
        sprint.snv_cluster(f"{anno_out}.nrp", f"{regular_res}.nrp", cluster_distance, cluster_size_nrp)
        
        # 合并RES
        sprint.combine_res(f"{regular_res}.alu", f"{regular_res}.nalurp", f"{regular_res}.nrp", f"{regular_res}.split")
        cluster_size_regular_max = max([cluster_size_alu_ad1, cluster_size_alu_ad2, cluster_size_nalurp, cluster_size_nrp])
        sprint.combine_res(f"{anno_out}.alu", f"{anno_out}.nalurp", f"{anno_out}.nrp", f"{anno_out}.rmsrp")
        sprint.snv_cluster(f"{anno_out}.rmsrp", f"{regular_res}.overall", cluster_distance, cluster_size_regular_max)
        sprint.res_or(f"{regular_res}.split", f"{regular_res}.overall", regular_res)
        
        # 清理注释临时文件
        for f in [anno_out, f"{anno_out}.alu", f"{anno_out}.alu.ad2", f"{anno_out}.nalurp", f"{anno_out}.nrp", f"{anno_out}.rmsrp"]:
            if os.path.exists(f):
                os.remove(f)
    else:
        # 无重复序列注释的RES识别
        sprint.snv_cluster(filtered_regular_snv, f"{regular_res}.tmp", cluster_distance, cluster_size_rg)
        sprint.o2b(f"{regular_res}.tmp", regular_res)
        if os.path.exists(f"{regular_res}.tmp"):
            os.remove(f"{regular_res}.tmp")
    
    # 2. Hyper MSKTC RES处理
    hyper_mskTC_snv = f"{tmp_dir}/hyper_mskTC.snv"
    filtered_hyper_mskTC = f"{hyper_mskTC_snv}.{chrom_prefix}"
    filter_snv(hyper_mskTC_snv, filtered_hyper_mskTC)
    
    if repeat_file:
        anno_out = f"{filtered_hyper_mskTC}.anno"
        sprint.annotate(filtered_hyper_mskTC, repeat_file, anno_out)
        sprint.seperate(anno_out, f"{anno_out}.alu", f"{anno_out}.nalurp", f"{anno_out}.nrp", 'Alu')
        sprint.snv_cluster(f"{anno_out}.alu", f"{hyper_mskTC_res}.alu", cluster_distance, cluster_size_alu_hp)
        sprint.snv_cluster(f"{anno_out}.nalurp", f"{hyper_mskTC_res}.nalurp", cluster_distance, cluster_size_nalurp_hp)
        sprint.snv_cluster(f"{anno_out}.nrp", f"{hyper_mskTC_res}.nrp", cluster_distance, cluster_size_nrp_hp)
        sprint.combine_res(f"{hyper_mskTC_res}.alu", f"{hyper_mskTC_res}.nalurp", f"{hyper_mskTC_res}.nrp", f"{hyper_mskTC_res}.split")
        cluster_size_hyper_max = max([cluster_size_alu_hp, cluster_size_nalurp_hp, cluster_size_nrp_hp])
        sprint.combine_res(f"{anno_out}.alu", f"{anno_out}.nalurp", f"{anno_out}.nrp", f"{anno_out}.rmsrp")
        sprint.snv_cluster(f"{anno_out}.rmsrp", f"{hyper_mskTC_res}.overall", cluster_distance, cluster_size_hyper_max)
        sprint.res_or(f"{hyper_mskTC_res}.split", f"{hyper_mskTC_res}.overall", hyper_mskTC_res)
        
        # 清理临时文件
        for f in [anno_out, f"{anno_out}.alu", f"{anno_out}.nalurp", f"{anno_out}.nrp", f"{anno_out}.rmsrp"]:
            if os.path.exists(f):
                os.remove(f)
    else:
        sprint.snv_cluster(filtered_hyper_mskTC, hyper_mskTC_res, cluster_distance, cluster_size_hp)
    
    # 3. Hyper MSKAG RES处理
    hyper_mskAG_snv = f"{tmp_dir}/hyper_mskAG.snv"
    filtered_hyper_mskAG = f"{hyper_mskAG_snv}.{chrom_prefix}"
    filter_snv(hyper_mskAG_snv, filtered_hyper_mskAG)
    
    if repeat_file:
        anno_out = f"{filtered_hyper_mskAG}.anno"
        sprint.annotate(filtered_hyper_mskAG, repeat_file, anno_out)
        sprint.seperate(anno_out, f"{anno_out}.alu", f"{anno_out}.nalurp", f"{anno_out}.nrp", 'Alu')
        sprint.snv_cluster(f"{anno_out}.alu", f"{hyper_mskAG_res}.alu", cluster_distance, cluster_size_alu_hp)
        sprint.snv_cluster(f"{anno_out}.nalurp", f"{hyper_mskAG_res}.nalurp", cluster_distance, cluster_size_nalurp_hp)
        sprint.snv_cluster(f"{anno_out}.nrp", f"{hyper_mskAG_res}.nrp", cluster_distance, cluster_size_nrp_hp)
        sprint.combine_res(f"{hyper_mskAG_res}.alu", f"{hyper_mskAG_res}.nalurp", f"{hyper_mskAG_res}.nrp", f"{hyper_mskAG_res}.split")
        cluster_size_hyper_max = max([cluster_size_alu_hp, cluster_size_nalurp_hp, cluster_size_nrp_hp])
        sprint.combine_res(f"{anno_out}.alu", f"{anno_out}.nalurp", f"{anno_out}.nrp", f"{anno_out}.rmsrp")
        sprint.snv_cluster(f"{anno_out}.rmsrp", f"{hyper_mskAG_res}.overall", cluster_distance, cluster_size_hyper_max)
        sprint.res_or(f"{hyper_mskAG_res}.split", f"{hyper_mskAG_res}.overall", hyper_mskAG_res)
        
        # 清理临时文件
        for f in [anno_out, f"{anno_out}.alu", f"{anno_out}.nalurp", f"{anno_out}.nrp", f"{anno_out}.rmsrp"]:
            if os.path.exists(f):
                os.remove(f)
    else:
        sprint.snv_cluster(filtered_hyper_mskAG, hyper_mskAG_res, cluster_distance, cluster_size_hp)
    
    # 清理过滤的SNV临时文件
    for tmp_f in [filtered_regular_snv, filtered_hyper_mskTC, filtered_hyper_mskAG]:
        if os.path.exists(tmp_f):
            os.remove(tmp_f)
    
    return {
        'chrom': chrom,
        'regular_res': regular_res,
        'hyper_mskTC_res': hyper_mskTC_res,
        'hyper_mskAG_res': hyper_mskAG_res
    }


def merge_parallel_results(result_files, output_file):
    """合并并行任务结果文件"""
    with open(output_file, 'w') as out_f:
        for res_file in result_files:
            if os.path.exists(res_file):
                with open(res_file, 'r') as in_f:
                    out_f.write(in_f.read())
                # 清理单染色体结果文件
                os.remove(res_file)


# ===================== 主函数 =====================
def main():
    print('')
    print("##############################################################################################")
    print('')
    print("      SPRINT: SNP-free RNA editing Identification Toolkit")
    print("")
    print("      https://github.com/TianLab-Bioinfo/SPRINT/blob/master/sprint_py3/")
    print("")
    print("      Please contact limh25@m.fudan.edu.cn when questions arise.")
    print("")
    print("##############################################################################################")

    def help_doc():
        print("")
        print("   Attention:")
        print("")
        print("      Before using 'sprint main', please use 'sprint prepare' to build mapping index.")
        print("")
        print("   Usage:")
        print("")
        print("      sprint main   [options]   reference_genome(.fa)   output_path   bwa_path   samtools_path")
        print("")
        print("      options:")
        print("         -1       read1(.fq)       # Required !!!")
        print("         -2       read2(.fq)       # Optional")
        print(
            "         -rp      repeat_file      # Optional, you can download it from http://sprint.software/SPRINT/dbrep/")
        print(
            "         -ss      INT              # when input is strand-specific sequencing data, please clarify the direction of read1. [0 for antisense; 1 for sense] (default is 0)")
        print("         -c       INT              # Remove the fist INT bp of each read (default is 0)")
        print("         -p       INT              # Mapping CPU (default is 1)")
        print("         -cd      INT              # The distance cutoff of SNV duplets (default is 200)")
        print(
            "         -csad1   INT              # Regular - [-rp is required] cluster size - Alu - AD >=1 (default is 3)")
        print(
            "         -csad2   INT              # Regular - [-rp is required] cluster size - Alu - AD >=2 (default is 2)")
        print(
            "         -csnar   INT              # Regular - [-rp is required] cluster size - nonAlu Repeat - AD >=1 (default is 5)")
        print(
            "         -csnr    INT              # Regular - [-rp is required] cluster size - nonRepeat - AD >=1 (default is 7)")
        print("         -csrg    INT              # Regular - [without -rp] cluster size - AD >=1 (default is 5)")
        print(
            "         -csahp   INT              # Hyper - [-rp is required] cluster size - Alu - AD >=1 (default is 5)")
        print(
            "         -csnarhp INT              # Hyper - [-rp is required] cluster size - nonAlu Repeat - AD >=1 (default is 5)")
        print(
            "         -csnrhp  INT              # Hyper - [-rp is required] cluster size - nonRepeat - AD >=1 (default is 5)")
        print("         -cshp    INT              # Hyper - [without -rp] cluster size - AD >=1 (default is 5)")
        print("         -t      INT              # Parallel CPU for SNV/RES (default: same as -p)")
        print("")
        print("   Example:")
        print("")
        print(
            "       sprint main -rp hg38_repeat.txt -c 6 -p 6 -t 6 -1 read1.fq -2 read2.fq hg38.fa output ./bwa-0.7.19/bwa ./samtools-1.21/samtools")
        print("")
        print("       Notes: Default protocol of strand-specific RNA-seq is dUTP (read1: '-'; read2: '+')")
        print("")
        sys.exit(0)

    if len(sys.argv) < 2:
        help_doc()

    # 初始化参数
    read_format = 0
    cutbp = 0
    cluster_distance = 200
    cluster_size_alu_ad1 = 3
    cluster_size_alu_ad2 = 2
    cluster_size_nalurp = 5
    cluster_size_nrp = 7
    cluster_size_rg = 5
    cluster_size_hp = 5
    cluster_size_alu_hp = 5
    cluster_size_nalurp_hp = 5
    cluster_size_nrp_hp = 5
    strand_specify = 0
    mapcpu = 1
    parallel_cpu = 1  # 新增：SNV/RES并行CPU数
    var_limit = 20
    poly_limit = 10
    rm_multi = 0
    paired_end = False
    repeat = False
    options = []
    read2 = ''
    read1 = ''

    # 解析命令行参数（新增-t参数）
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '-1':
            try:
                read1 = sys.argv[i + 1]
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
        elif sys.argv[i] == '-2':
            paired_end = True
            try:
                read2 = sys.argv[i + 1]
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-rp':
            try:
                repeat = sys.argv[i + 1]
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-ss':
            try:
                strand_specify = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-c':
            try:
                cutbp = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-p':
            try:
                mapcpu = int(sys.argv[i + 1])
                parallel_cpu = mapcpu  # 默认并行CPU数等于mapping CPU
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-t':  # 新增：指定SNV/RES并行CPU数
            try:
                parallel_cpu = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-cd':
            try:
                cluster_distance = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-csad1':
            try:
                cluster_size_alu_ad1 = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-csad2':
            try:
                cluster_size_alu_ad2 = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-csnar':
            try:
                cluster_size_nalurp = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-csnr':
            try:
                cluster_size_nrp = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-csrg':
            try:
                cluster_size_rg = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-cshp':
            try:
                cluster_size_hp = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-csahp':
            try:
                cluster_size_alu_hp = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-csnarhp':
            try:
                cluster_size_nalurp_hp = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        elif sys.argv[i] == '-csnrhp':
            try:
                cluster_size_nrp_hp = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                sys.exit()
        i += 1

    # 提取非选项参数
    all_argv = []
    i = 1
    while i < len(sys.argv):
        if i not in options:
            all_argv.append(i)
        i += 1

    if len(all_argv) != 4 or read1 == '':
        help_doc()
        sys.exit()

    refgenome = sys.argv[all_argv[0]]
    output = sys.argv[all_argv[1]] + '/'
    tmp = output + '/tmp/'
    bwa = sys.argv[all_argv[2]]
    samtools = sys.argv[all_argv[3]]

    # 创建目录
    if not os.path.exists(output):
        os.mkdir(output)
    if not os.path.exists(tmp):
        os.mkdir(tmp)

    # 写入参数文件
    try:
        frc = open(tmp + 'PARAMETER.txt', 'w')
        frc.write(sys.argv[0])
        for one in sys.argv[1:]:
            frc.write('   ' + one)
        frc.write('\n')
        frc.close()
    except Exception as e:
        print(f"参数文件写入失败：{e}")
        sys.exit(1)

    # 核心修改：适配新版 Samtools 的 fq2sam 函数
    def fq2sam(TAG, paired_end, read1, read2, tmp, refgenome, bwa, samtools, mapcpu, read_format):
        if paired_end is True:
            mapcpu = max([int(int(mapcpu) / 2.0), 1])
        ori_tmp = tmp
        tmp = tmp + '/' + TAG + '/'
        if not os.path.exists(tmp):
            os.mkdir(tmp)

        # Step 1: BWA aln 比对
        step1_1 = subprocess.Popen(f"{bwa} aln -t {mapcpu} {refgenome} {read1} > {tmp}read1.sai", shell=True)
        if paired_end is True:
            step1_2 = subprocess.Popen(f"{bwa} aln -t {mapcpu} {refgenome} {read2} > {tmp}read2.sai", shell=True)

        step1_1.wait()
        if paired_end is True:
            step1_2.wait()

        # Step 2: BWA samse 生成SAM
        step1_3 = subprocess.Popen(f"{bwa} samse -n4 {refgenome} {tmp}read1.sai {read1} > {tmp}name_read1.sam",
                                   shell=True)
        if paired_end is True:
            step1_4 = subprocess.Popen(f"{bwa} samse -n4 {refgenome} {tmp}read2.sai {read2} > {tmp}name_read2.sam",
                                       shell=True)

        step1_3.wait()
        if paired_end is True:
            step1_4.wait()

        # 清理临时文件
        if os.path.exists(tmp + 'name_read1.sam'):
            if os.path.exists(tmp + 'read1.sai'):
                os.remove(tmp + 'read1.sai')
            if os.path.exists(ori_tmp + 'cut_read1.fastq'):
                os.remove(ori_tmp + 'cut_read1.fastq')
        if os.path.exists(tmp + 'name_read2.sam'):
            if os.path.exists(tmp + 'read2.sai'):
                os.remove(tmp + 'read2.sai')
            if os.path.exists(ori_tmp + 'cut_read2.fastq'):
                os.remove(ori_tmp + 'cut_read2.fastq')

        # Step 3: Samtools view -b 生成BAM（废弃 -S 参数）
        step1_7 = subprocess.Popen(f"{samtools} view -b {tmp}name_read1.sam > {tmp}name_read1.bam", shell=True)
        if paired_end is True:
            step1_8 = subprocess.Popen(f"{samtools} view -b {tmp}name_read2.sam > {tmp}name_read2.bam", shell=True)

        step1_7.wait()
        if paired_end is True:
            step1_8.wait()

        # Step 4: Samtools sort 排序（核心修改：使用 -o 指定输出文件）
        if paired_end is True:
            # 新版 sort 语法：samtools sort -o <output.bam> <input.bam>
            step1_9 = subprocess.Popen(f"{samtools} sort -@ {mapcpu} -o {tmp}name_read1_sorted.bam {tmp}name_read1.bam",
                                       shell=True)
            step1_10 = subprocess.Popen(
                f"{samtools} sort -@ {mapcpu} -o {tmp}name_read2_sorted.bam {tmp}name_read2.bam", shell=True)
            step1_9.wait()
            step1_10.wait()

            # Step 5: Samtools merge 合并（显式指定 -o 提升兼容性）
            step1_11 = subprocess.Popen(
                f"{samtools} merge -f -o {tmp}all.bam {tmp}name_read1_sorted.bam {tmp}name_read2_sorted.bam",
                shell=True)
            step1_11.wait()

            # 清理临时文件
            if os.path.exists(tmp + 'all.bam'):
                for f in ['name_read1.sam', 'name_read1.bam', 'name_read1_sorted.bam',
                          'name_read2.sam', 'name_read2.bam', 'name_read2_sorted.bam']:
                    if os.path.exists(tmp + f):
                        os.remove(tmp + f)
        else:
            # 单端数据排序（新版语法）
            step1_9 = subprocess.Popen(f"{samtools} sort -@ {mapcpu} -o {tmp}all.bam {tmp}name_read1.bam", shell=True)
            step1_9.wait()

            # 清理临时文件
            if os.path.exists(tmp + 'all.bam'):
                for f in ['name_read1.sam', 'name_read1.bam']:
                    if os.path.exists(tmp + f):
                        os.remove(tmp + f)

        # Step 6: 转换回SAM
        step2_2 = subprocess.Popen(f"{samtools} view -h -o {tmp}all.sam {tmp}all.bam", shell=True)
        step2_2.wait()
        subprocess.Popen(f"cp {tmp}/all.sam {ori_tmp}/{TAG}_all.sam", shell=True).wait()

        if os.path.exists(tmp + 'all.sam'):
            os.remove(tmp + 'all.sam')

    # 后续逻辑（预处理 + 比对）
    print('preprocessing...')
    if read_format != 0:
        subprocess.Popen(f"{samtools} view -h -o {tmp}read1.sam {read1}", shell=True).wait()
        sprint.sam2fq(tmp + 'read1.sam', tmp + 'read1.fq')
        read1 = tmp + 'read1.fq'
        sprint.cut(read1, tmp + 'cut_read1.fastq', cutbp, 'read1')
        if paired_end == True:
            subprocess.Popen(f"{samtools} view -h -o {tmp}read2.sam {read2}", shell=True).wait()
            sprint.sam2fq(tmp + 'read2.sam', tmp + 'read2.fq')
            read2 = tmp + 'read2.fq'
            sprint.cut(read2, tmp + 'cut_read2.fastq', cutbp, 'read2')
    else:
        if strand_specify == 0:
            sprint.cut(read1, tmp + 'cut_read1.fastq', cutbp, 'read1')
            if paired_end == True:
                sprint.cut(read2, tmp + 'cut_read2.fastq', cutbp, 'read2')
        else:
            sprint.cut(read1, tmp + 'cut_read1.fastq', cutbp, 'read2')
            if paired_end == True:
                sprint.cut(read2, tmp + 'cut_read2.fastq', cutbp, 'read1')

    sprint.get_baseq_cutoff(read1, tmp + 'baseq.cutoff')

    print('mapping...')
    TAG = 'genome'
    fq2sam(TAG, paired_end, tmp + 'cut_read1.fastq', tmp + 'cut_read2.fastq', tmp, refgenome, bwa, samtools, mapcpu,
           read_format)

    subprocess.Popen(f"{samtools} view -f4 {tmp}/{TAG}/all.bam > {tmp}/{TAG}_unmapped.sam", shell=True).wait()
    sprint.umsam2fq(tmp + '/' + TAG + '_unmapped.sam', tmp + '/' + TAG + '_unmapped.fq')

    has_trans = os.path.exists(refgenome + '.trans.fa')
    if has_trans:
        TAG = 'transcript'
        fq2sam(TAG, False, tmp + '/genome_unmapped.fq', read2, tmp, refgenome + '.trans.fa', bwa, samtools, mapcpu,
               read_format)

        subprocess.Popen(f"{samtools} view -f4 {tmp}/{TAG}/all.bam > {tmp}/{TAG}_unmapped.sam", shell=True).wait()
        sprint.umsam2fq(tmp + '/' + TAG + '_unmapped.sam', tmp + '/regular_unmapped.fq')
        sprint.maskfq(tmp + '/regular_unmapped.fq', 'A', 'G')
    else:
        sprint.umsam2fq(tmp + '/' + TAG + '_unmapped.sam', tmp + '/regular_unmapped.fq')
        sprint.maskfq(tmp + '/regular_unmapped.fq', 'A', 'G')

    TAG = 'genome_mskAG'
    fq2sam(TAG, False, tmp + '/regular_unmapped_A_to_G.fq', read2, tmp, refgenome + '.mskAG.fa', bwa, samtools, mapcpu,
           read_format)

    subprocess.Popen(f"{samtools} view -f4 {tmp}/{TAG}/all.bam > {tmp}/{TAG}_unmapped.sam", shell=True).wait()
    sprint.umsam2fq(tmp + '/' + TAG + '_unmapped.sam', tmp + '/' + TAG + '_unmapped.fq')

    TAG = 'genome_mskTC'
    fq2sam(TAG, False, tmp + '/regular_unmapped_A_to_G.fq', read2, tmp, refgenome + '.mskTC.fa', bwa, samtools, mapcpu,
           read_format)

    subprocess.Popen(f"{samtools} view -f4 {tmp}/{TAG}/all.bam > {tmp}/{TAG}_unmapped.sam", shell=True).wait()
    sprint.umsam2fq(tmp + '/' + TAG + '_unmapped.sam', tmp + '/' + TAG + '_unmapped.fq')

    if has_trans:
        TAG = 'transcript_mskAG'
        fq2sam(TAG, False, tmp + '/genome_mskAG_unmapped.fq', read2, tmp, refgenome + '.trans.fa.mskAG.fa', bwa,
               samtools, mapcpu, read_format)

        TAG = 'transcript_mskTC'
        fq2sam(TAG, False, tmp + '/genome_mskTC_unmapped.fq', read2, tmp, refgenome + '.trans.fa.mskTC.fa', bwa,
               samtools, mapcpu, read_format)

    # 清理临时文件
    if os.path.exists(tmp + 'genome_mskAG_unmapped.sam'):
        for f in ['cut_read1.fastq', 'cut_read2.fastq', 'genome_mskAG_unmapped.fq', 'genome_mskAG_unmapped.sam',
                  'genome_mskTC_unmapped.fq', 'genome_mskTC_unmapped.sam', 'genome_unmapped.fq', 'genome_unmapped.sam',
                  'transcript_unmapped_A_to_G.fq', 'transcript_unmapped.fq', 'transcript_unmapped.sam',
                  'regular_unmapped.fq', 'regular_unmapped_A_to_G.fq']:
            if os.path.exists(tmp + f):
                os.remove(tmp + f)

    # 处理SAM/ZZ文件（原有逻辑）
    if has_trans:
        sprint.recover_sam(tmp + 'transcript_mskAG_all.sam', tmp + 'transcript_mskAG_all.sam.rcv', var_limit,
                           poly_limit, rm_multi)
        sprint.sam2zz(tmp + 'transcript_mskAG_all.sam.rcv', refgenome + '.trans.fa', tmp + 'transcript_mskAG_all.zz')
        sprint.recover_sam(tmp + 'transcript_mskTC_all.sam', tmp + 'transcript_mskTC_all.sam.rcv', var_limit,
                           poly_limit, rm_multi)
        sprint.sam2zz(tmp + 'transcript_mskTC_all.sam.rcv', refgenome + '.trans.fa', tmp + 'transcript_mskTC_all.zz')
        sprint.sam2zz(tmp + 'transcript_all.sam', refgenome + '.trans.fa', tmp + 'transcript_all.zz')

        for f in ['transcript_mskAG_all.sam.rcv', 'transcript_mskAG_all.sam', 'transcript_mskTC_all.sam.rcv',
                  'transcript_mskTC_all.sam', 'transcript_all.sam']:
            if os.path.exists(tmp + f):
                os.remove(tmp + f)

    sprint.recover_sam(tmp + 'genome_mskAG_all.sam', tmp + 'genome_mskAG_all.sam.rcv', var_limit, poly_limit, rm_multi)
    sprint.sam2zz(tmp + 'genome_mskAG_all.sam.rcv', refgenome, tmp + 'genome_mskAG_all.zz')
    sprint.recover_sam(tmp + 'genome_mskTC_all.sam', tmp + 'genome_mskTC_all.sam.rcv', var_limit, poly_limit, rm_multi)
    sprint.sam2zz(tmp + 'genome_mskTC_all.sam.rcv', refgenome, tmp + 'genome_mskTC_all.zz')
    sprint.sam2zz(tmp + 'genome_all.sam', refgenome, tmp + 'genome_all.zz')

    for f in ['genome_mskAG_all.sam.rcv', 'genome_mskAG_all.sam', 'genome_mskTC_all.sam.rcv',
              'genome_mskTC_all.sam', 'genome_all.sam']:
        if os.path.exists(tmp + f):
            os.remove(tmp + f)

    if has_trans:
        sprint.dedup(tmp + 'transcript_mskAG_all.zz', tmp + 'transcript_mskAG_all.zz.dedup')
        sprint.dedup(tmp + 'transcript_mskTC_all.zz', tmp + 'transcript_mskTC_all.zz.dedup')
        sprint.dedup(tmp + 'transcript_all.zz', tmp + 'transcript_all.zz.dedup')

    sprint.dedup(tmp + 'genome_mskAG_all.zz', tmp + 'genome_mskAG_all.zz.dedup')
    sprint.dedup(tmp + 'genome_mskTC_all.zz', tmp + 'genome_mskTC_all.zz.dedup')
    sprint.dedup(tmp + 'genome_all.zz', tmp + 'genome_all.zz.dedup')

    for f in ['transcript_mskAG_all.zz', 'transcript_mskTC_all.zz', 'transcript_all.zz',
              'genome_mskAG_all.zz', 'genome_mskTC_all.zz', 'genome_all.zz']:
        if os.path.exists(tmp + f):
            os.remove(tmp + f)

    # 合并zz文件（原有逻辑）
    if has_trans:
        subprocess.Popen(
            f"cat {tmp}/genome_mskAG_all.zz.dedup {tmp}/genome_mskTC_all.zz.dedup {tmp}/genome_all.zz.dedup {tmp}/transcript_mskAG_all.zz.dedup.genome.zz {tmp}/transcript_mskTC_all.zz.dedup.genome.zz {tmp}/transcript_all.zz.dedup.genome.zz > {tmp}/all_combined.zz",
            shell=True).wait()
    else:
        subprocess.Popen(
            f"cat {tmp}/genome_mskAG_all.zz.dedup {tmp}/genome_mskTC_all.zz.dedup {tmp}/genome_all.zz.dedup > {tmp}/all_combined.zz",
            shell=True).wait()
    sprint.sort_zz(tmp + '/all_combined.zz', tmp + '/all_combined.zz.sorted')

    # ===================== 并行化SNV识别 =====================
    print('identifying SNVs (parallel)...')
    # 获取染色体列表
    chrom_list = get_chrom_list(tmp + '/all_combined.zz.sorted')
    baseq_cutoff = tmp + 'baseq.cutoff'
    
    # 构建并行任务参数
    snv_task_args = [
        (chrom, tmp, refgenome, has_trans, baseq_cutoff)
        for chrom in chrom_list
    ]
    
    # 启动并行池
    with Pool(processes=parallel_cpu) as pool:
        # 并行执行 + 进度条
        snv_results = list(tqdm(
            pool.imap_unordered(process_snv_by_chrom, snv_task_args),
            total=len(chrom_list),
            desc="SNV Processing"
        ))
    
    # 合并SNV结果
    # 1. 基因组SNV合并
    genome_mskAG_snv_files = [r['genome_mskAG_snv'] for r in snv_results]
    genome_mskTC_snv_files = [r['genome_mskTC_snv'] for r in snv_results]
    genome_all_snv_files = [r['genome_all_snv'] for r in snv_results]
    
    merge_parallel_results(genome_mskAG_snv_files, f"{tmp}/genome_mskAG_all.zz.dedup.snv")
    merge_parallel_results(genome_mskTC_snv_files, f"{tmp}/genome_mskTC_all.zz.dedup.snv")
    merge_parallel_results(genome_all_snv_files, f"{tmp}/genome_all.zz.dedup.snv")
    
    # 2. 转录组SNV合并（如有）
    if has_trans:
        trans_snv_types = ['transcript_mskAG', 'transcript_mskTC', 'transcript_all']
        for idx, trans_type in enumerate(trans_snv_types):
            trans_sort_files = [r['trans_snv'][idx][0] for r in snv_results if r['trans_snv']]
            merge_parallel_results(trans_sort_files, f"{tmp}/{trans_type}_all.zz.dedup.snv.genome.snv.sort")
        
        # 原有SNV_OR逻辑
        sprint.snv_or(f"{tmp}/transcript_all.zz.dedup.snv.genome.snv.sort", f"{tmp}/genome_all.zz.dedup.snv", f"{tmp}/regular.snv")
        sprint.snv_or(f"{tmp}/transcript_mskTC_all.zz.dedup.snv.genome.snv.sort", f"{tmp}/genome_mskTC_all.zz.dedup.snv", f"{tmp}/hyper_mskTC.snv")
        sprint.snv_or(f"{tmp}/transcript_mskAG_all.zz.dedup.snv.genome.snv.sort", f"{tmp}/genome_mskAG_all.zz.dedup.snv", f"{tmp}/hyper_mskAG.snv")
    else:
        subprocess.Popen(f"cp {tmp}/genome_all.zz.dedup.snv {tmp}/regular.snv", shell=True).wait()
        subprocess.Popen(f"cp {tmp}/genome_mskTC_all.zz.dedup.snv {tmp}/hyper_mskTC.snv", shell=True).wait()
        subprocess.Popen(f"cp {tmp}/genome_mskAG_all.zz.dedup.snv {tmp}/hyper_mskAG.snv", shell=True).wait()

    # ===================== 并行化RES识别 =====================
    print('identifying RESs (parallel)...')
    # 构建聚类参数元组
    cluster_params = (
        cluster_distance, cluster_size_alu_ad1, cluster_size_alu_ad2,
        cluster_size_nalurp, cluster_size_nrp, cluster_size_rg,
        cluster_size_hp, cluster_size_alu_hp, cluster_size_nalurp_hp,
        cluster_size_nrp_hp
    )
    
    # 构建RES并行任务参数
    res_task_args = [
        (chrom, tmp, repeat, cluster_params, has_trans)
        for chrom in chrom_list
    ]
    
    # 启动并行池
    with Pool(processes=parallel_cpu) as pool:
        # 并行执行 + 进度条
        res_results = list(tqdm(
            pool.imap_unordered(process_res_by_chrom, res_task_args),
            total=len(chrom_list),
            desc="RES Processing"
        ))
    
    # 合并RES结果
    regular_res_files = [r['regular_res'] for r in res_results]
    hyper_mskTC_res_files = [r['hyper_mskTC_res'] for r in res_results]
    hyper_mskAG_res_files = [r['hyper_mskAG_res'] for r in res_results]
    
    merge_parallel_results(regular_res_files, f"{tmp}/regular.res")
    merge_parallel_results(hyper_mskTC_res_files, f"{tmp}/hyper_mskTC.res")
    merge_parallel_results(hyper_mskAG_res_files, f"{tmp}/hyper_mskAG.res")
    
    # Hyper RES合并
    sprint.snv_or(tmp + '/hyper_mskTC.res', tmp + '/hyper_mskAG.res', tmp + '/hyper.res')

    # 清理临时注释文件
    try:
        subprocess.Popen(f"rm -rf {tmp}/*.anno.*", shell=True).wait()
    except Exception as e:
        pass

    # 输出结果文件（原有逻辑）
    sprint.get_depth(tmp + '/all_combined.zz.sorted', tmp + '/regular.res', tmp + '/regular.res.depth')
    subprocess.Popen(
        f"echo '#Chrom\tStart(0base)\tEnd(1base)\tType\tSupporting_reads\tStrand\tAD:DP' | cat - {tmp}/regular.res.depth > {output}/SPRINT_identified_regular.res",
        shell=True).wait()
    sprint.get_depth(tmp + '/all_combined.zz.sorted', tmp + '/hyper.res', tmp + '/hyper.res.depth')
    subprocess.Popen(
        f"echo '#Chrom\tStart(0base)\tEnd(1base)\tType\tSupporting_reads\tStrand\tAD:DP' | cat - {tmp}/hyper.res.depth > {output}/SPRINT_identified_hyper.res",
        shell=True).wait()
    subprocess.Popen(f"cp {tmp}/PARAMETER.txt {output}/PARAMETER.txt", shell=True).wait()

    sprint.snv_or(tmp + '/regular.res', tmp + '/hyper.res', tmp + '/all.res')
    sprint.get_depth(tmp + '/all_combined.zz.sorted', tmp + '/all.res', tmp + '/all.res.depth')
    subprocess.Popen(
        f"echo '#Chrom\tStart(0base)\tEnd(1base)\tType\tSupporting_reads\tStrand\tAD:DP' | cat - {tmp}/all.res.depth > {output}/SPRINT_identified_all.res",
        shell=True).wait()

    print('finished !')
    sys.exit(0)


if __name__ == '__main__':
    main()