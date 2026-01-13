import sys
import os
import subprocess
import multiprocessing
from tqdm import tqdm
import sprint  # 假设sprint模块已存在


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
        print("")
        print("   Example:")
        print("")
        print(
            "       sprint main -rp hg38_repeat.txt -c 6 -p 6 -1 read1.fq -2 read2.fq hg38.fa output ./bwa-0.7.19/bwa ./samtools-1.21/samtools")
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
    mapcpu = 1  # 映射CPU数，并行任务CPU数基于此调整
    var_limit = 20
    poly_limit = 10
    rm_multi = 0
    paired_end = False
    repeat = False
    options = []
    read2 = ''
    read1 = ''

    # 解析命令行参数
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

    # 适配新版 Samtools 的 fq2sam 函数（原逻辑不变）
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

    # -------------------------- 并行任务定义 --------------------------
    # 任务1：SNV识别子任务
    def process_snv_task(task_info):
        """
        并行处理SNV识别子任务
        task_info格式：(task_name, input_zz, output_snv, baseq_cutoff, refgenome, trans_loc=None)
        """
        task_name, input_zz, output_snv, baseq_cutoff, refgenome, trans_loc = task_info
        try:
            # 基础SNV识别
            sprint.mask_zz2snv(input_zz, output_snv, baseq_cutoff)

            # 如果是转录组任务，额外处理基因组坐标转换
            if trans_loc and os.path.exists(trans_loc):
                genome_snv = output_snv + '.genome.snv'
                sprint.transcript_locator(output_snv, trans_loc, genome_snv)
                sprint.transcript_sort(genome_snv, genome_snv + '.sort')
                return (task_name, "success", genome_snv + '.sort')
            return (task_name, "success", output_snv)
        except Exception as e:
            return (task_name, f"failed: {str(e)}", output_snv)

    # 任务2：RES识别子任务
    def process_res_task(task_info):
        """
        并行处理RES识别子任务
        task_info格式：(task_name, input_snv, repeat_file, output_res, cluster_distance, cluster_size, is_hyper=False)
        """
        task_name, input_snv, repeat_file, output_res, cluster_distance, cluster_size, is_hyper = task_info
        try:
            if repeat_file:
                # 带repeat注释的RES识别
                anno_snv = input_snv + '.anno'
                sprint.annotate(input_snv, repeat_file, anno_snv)

                # 拆分Alu/nonAlu Repeat/nonRepeat
                alu_snv = anno_snv + '.alu'
                nalurp_snv = anno_snv + '.nalurp'
                nrp_snv = anno_snv + '.nrp'
                sprint.seperate(anno_snv, alu_snv, nalurp_snv, nrp_snv, 'Alu')

                if not is_hyper:
                    # Regular模式：处理AD>=2的Alu位点
                    alu_ad2_snv = alu_snv + '.ad2'
                    sprint.get_snv_with_ad(alu_snv, alu_ad2_snv, 2)
                    alu_res_ad1 = output_res + '.alu.ad1'
                    alu_res_ad2 = output_res + '.alu.ad2'
                    sprint.snv_cluster(alu_snv, alu_res_ad1, cluster_distance, cluster_size['alu_ad1'])
                    sprint.snv_cluster(alu_ad2_snv, alu_res_ad2, cluster_distance, cluster_size['alu_ad2'])
                    sprint.bed_or(alu_res_ad1, alu_res_ad2, output_res + '.alu')

                    # 处理nonAlu Repeat/nonRepeat
                    nalurp_res = output_res + '.nalurp'
                    nrp_res = output_res + '.nrp'
                    sprint.snv_cluster(nalurp_snv, nalurp_res, cluster_distance, cluster_size['nalurp'])
                    sprint.snv_cluster(nrp_snv, nrp_res, cluster_distance, cluster_size['nrp'])

                    # 合并结果
                    sprint.combine_res(output_res + '.alu', nalurp_res, nrp_res, output_res + '.split')
                    cluster_size_max = max(cluster_size['alu_ad1'], cluster_size['alu_ad2'], cluster_size['nalurp'],
                                           cluster_size['nrp'])
                    sprint.combine_res(alu_snv, nalurp_snv, nrp_snv, anno_snv + '.rmsrp')
                    sprint.snv_cluster(anno_snv + '.rmsrp', output_res + '.overall', cluster_distance, cluster_size_max)
                    sprint.res_or(output_res + '.split', output_res + '.overall', output_res)
                else:
                    # Hyper模式：直接聚类
                    alu_res = output_res + '.alu'
                    nalurp_res = output_res + '.nalurp'
                    nrp_res = output_res + '.nrp'
                    sprint.snv_cluster(alu_snv, alu_res, cluster_distance, cluster_size['alu_hp'])
                    sprint.snv_cluster(nalurp_snv, nalurp_res, cluster_distance, cluster_size['nalurp_hp'])
                    sprint.snv_cluster(nrp_snv, nrp_res, cluster_distance, cluster_size['nrp_hp'])
                    sprint.combine_res(alu_res, nalurp_res, nrp_res, output_res + '.split')
                    cluster_size_max = max(cluster_size['alu_hp'], cluster_size['nalurp_hp'], cluster_size['nrp_hp'])
                    sprint.combine_res(alu_snv, nalurp_snv, nrp_snv, anno_snv + '.rmsrp')
                    sprint.snv_cluster(anno_snv + '.rmsrp', output_res + '.overall', cluster_distance, cluster_size_max)
                    sprint.res_or(output_res + '.split', output_res + '.overall', output_res)
            else:
                # 无repeat注释的RES识别
                sprint.snv_cluster(input_snv, output_res + '.tmp', cluster_distance, cluster_size)
                sprint.o2b(output_res + '.tmp', output_res)
            return (task_name, "success", output_res)
        except Exception as e:
            return (task_name, f"failed: {str(e)}", output_res)

    # -------------------------- 预处理（原逻辑不变） --------------------------
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

    if os.path.exists(refgenome + '.trans.fa'):
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

    if os.path.exists(refgenome + '.trans.fa'):
        TAG = 'transcript_mskAG'
        fq2sam(TAG, False, tmp + '/genome_mskAG_unmapped.fq', read2, tmp, refgenome + '.trans.fa.mskAG.fa', bwa,
               samtools, mapcpu, read_format)

        TAG = 'transcript_mskTC'
        fq2sam(TAG, False, tmp + '/genome_mskTC_unmapped.fq', read2, tmp, refgenome + '.trans.fa.mskTC.fa', bwa,
               samtools, mapcpu, read_format)

    # 清理临时文件（原逻辑不变）
    if os.path.exists(tmp + 'genome_mskAG_unmapped.sam'):
        for f in ['cut_read1.fastq', 'cut_read2.fastq', 'genome_mskAG_unmapped.fq', 'genome_mskAG_unmapped.sam',
                  'genome_mskTC_unmapped.fq', 'genome_mskTC_unmapped.sam', 'genome_unmapped.fq', 'genome_unmapped.sam',
                  'transcript_unmapped_A_to_G.fq', 'transcript_unmapped.fq', 'transcript_unmapped.sam',
                  'regular_unmapped.fq', 'regular_unmapped_A_to_G.fq']:
            if os.path.exists(tmp + f):
                os.remove(tmp + f)

    # 转录组SAM处理（原逻辑不变）
    if os.path.exists(refgenome + '.trans.fa'):
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

    if os.path.exists(refgenome + '.trans.fa'):
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

    # -------------------------- 并行执行SNV识别 --------------------------
    print('identifying SNVs (parallel)...')
    # 构建SNV任务列表
    snv_tasks = []
    baseq_cutoff = tmp + 'baseq.cutoff'
    trans_loc = refgenome + '.trans.fa.loc' if os.path.exists(refgenome + '.trans.fa') else None

    # 基因组SNV任务
    snv_tasks.append(('genome_mskAG', tmp + 'genome_mskAG_all.zz.dedup', tmp + 'genome_mskAG_all.zz.dedup.snv',
                      baseq_cutoff, refgenome, None))
    snv_tasks.append(('genome_mskTC', tmp + 'genome_mskTC_all.zz.dedup', tmp + 'genome_mskTC_all.zz.dedup.snv',
                      baseq_cutoff, refgenome, None))
    snv_tasks.append(
        ('genome', tmp + 'genome_all.zz.dedup', tmp + 'genome_all.zz.dedup.snv', baseq_cutoff, refgenome, None))

    # 转录组SNV任务（如果存在）
    if trans_loc:
        snv_tasks.append(('transcript_mskAG', tmp + 'transcript_mskAG_all.zz.dedup',
                          tmp + 'transcript_mskAG_all.zz.dedup.snv', baseq_cutoff, refgenome, trans_loc))
        snv_tasks.append(('transcript_mskTC', tmp + 'transcript_mskTC_all.zz.dedup',
                          tmp + 'transcript_mskTC_all.zz.dedup.snv', baseq_cutoff, refgenome, trans_loc))
        snv_tasks.append(('transcript', tmp + 'transcript_all.zz.dedup', tmp + 'transcript_all.zz.dedup.snv',
                          baseq_cutoff, refgenome, trans_loc))

    # 启动并行池（CPU数=mapcpu）
    pool = multiprocessing.Pool(processes=mapcpu)
    snv_results = []

    # 用tqdm展示并行进度
    for result in tqdm(pool.imap_unordered(process_snv_task, snv_tasks), total=len(snv_tasks), desc="SNV Processing"):
        snv_results.append(result)
        task_name, status, output_file = result
        if "failed" in status:
            print(f"Warning: Task {task_name} failed - {status}")

    # 关闭并行池
    pool.close()
    pool.join()

    # 检查SNV任务结果
    for res in snv_results:
        if "failed" in res[1]:
            print(f"Error: SNV task {res[0]} failed, exit!")
            sys.exit(1)

    # 转录组zz转基因组zz（原逻辑）
    if trans_loc:
        sprint.tzz2gzz(trans_loc, tmp + 'transcript_mskAG_all.zz.dedup',
                       tmp + 'transcript_mskAG_all.zz.dedup.genome.zz')
        sprint.tzz2gzz(trans_loc, tmp + 'transcript_mskTC_all.zz.dedup',
                       tmp + 'transcript_mskTC_all.zz.dedup.genome.zz')
        sprint.tzz2gzz(trans_loc, tmp + 'transcript_all.zz.dedup', tmp + 'transcript_all.zz.dedup.genome.zz')

    # 合并zz文件（原逻辑）
    if os.path.exists(refgenome + '.trans.fa'):
        subprocess.Popen(
            f"cat {tmp}/genome_mskAG_all.zz.dedup {tmp}/genome_mskTC_all.zz.dedup {tmp}/genome_all.zz.dedup {tmp}/transcript_mskAG_all.zz.dedup.genome.zz {tmp}/transcript_mskTC_all.zz.dedup.genome.zz {tmp}/transcript_all.zz.dedup.genome.zz > {tmp}/all_combined.zz",
            shell=True).wait()
    else:
        subprocess.Popen(
            f"cat {tmp}/genome_mskAG_all.zz.dedup {tmp}/genome_mskTC_all.zz.dedup {tmp}/genome_all.zz.dedup > {tmp}/all_combined.zz",
            shell=True).wait()
    sprint.sort_zz(tmp + '/all_combined.zz', tmp + '/all_combined.zz.sorted')

    # 生成SNV合并文件（原逻辑）
    if os.path.exists(refgenome + '.trans.fa'):
        sprint.snv_or(tmp + 'transcript_all.zz.dedup.snv.genome.snv.sort', tmp + 'genome_all.zz.dedup.snv',
                      tmp + 'regular.snv')
        sprint.snv_or(tmp + 'transcript_mskTC_all.zz.dedup.snv.genome.snv.sort', tmp + 'genome_mskTC_all.zz.dedup.snv',
                      tmp + 'hyper_mskTC.snv')
        sprint.snv_or(tmp + 'transcript_mskAG_all.zz.dedup.snv.genome.snv.sort', tmp + 'genome_mskAG_all.zz.dedup.snv',
                      tmp + 'hyper_mskAG.snv')
    else:
        subprocess.Popen(f"cp {tmp}/genome_all.zz.dedup.snv {tmp}/regular.snv", shell=True).wait()
        subprocess.Popen(f"cp {tmp}/genome_mskTC_all.zz.dedup.snv {tmp}/hyper_mskTC.snv", shell=True).wait()
        subprocess.Popen(f"cp {tmp}/genome_mskAG_all.zz.dedup.snv {tmp}/hyper_mskAG.snv", shell=True).wait()

    # -------------------------- 并行执行RES识别 --------------------------
    print('identifying RESs (parallel)...')
    # 构建RES任务列表
    res_tasks = []
    repeat_file = repeat if repeat else None

    # Regular RES任务
    if repeat_file:
        regular_cluster_size = {
            'alu_ad1': cluster_size_alu_ad1,
            'alu_ad2': cluster_size_alu_ad2,
            'nalurp': cluster_size_nalurp,
            'nrp': cluster_size_nrp
        }
        res_tasks.append(('regular', tmp + 'regular.snv', repeat_file, tmp + 'regular.res', cluster_distance,
                          regular_cluster_size, False))
    else:
        res_tasks.append(
            ('regular', tmp + 'regular.snv', None, tmp + 'regular.res', cluster_distance, cluster_size_rg, False))

    # Hyper RES任务
    hyper_cluster_size = {
        'alu_hp': cluster_size_alu_hp,
        'nalurp_hp': cluster_size_nalurp_hp,
        'nrp_hp': cluster_size_nrp_hp
    }
    res_tasks.append(('hyper_mskTC', tmp + 'hyper_mskTC.snv', repeat_file, tmp + 'hyper_mskTC.res', cluster_distance,
                      hyper_cluster_size, True))
    res_tasks.append(('hyper_mskAG', tmp + 'hyper_mskAG.snv', repeat_file, tmp + 'hyper_mskAG.res', cluster_distance,
                      hyper_cluster_size, True))

    # 启动并行池执行RES任务
    pool = multiprocessing.Pool(processes=mapcpu)
    res_results = []
    for result in tqdm(pool.imap_unordered(process_res_task, res_tasks), total=len(res_tasks), desc="RES Processing"):
        res_results.append(result)
        task_name, status, output_file = result
        if "failed" in status:
            print(f"Warning: Task {task_name} failed - {status}")

    # 关闭并行池
    pool.close()
    pool.join()

    # 检查RES任务结果
    for res in res_results:
        if "failed" in res[1]:
            print(f"Error: RES task {res[0]} failed, exit!")
            sys.exit(1)

    # 合并Hyper RES结果（原逻辑）
    sprint.snv_or(tmp + 'hyper_mskTC.res', tmp + 'hyper_mskAG.res', tmp + 'hyper.res')

    # 清理anno临时文件
    try:
        subprocess.Popen(f"rm -rf {tmp}/*.anno.*", shell=True).wait()
    except Exception as e:
        pass

    # -------------------------- 输出最终结果（原逻辑不变） --------------------------
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