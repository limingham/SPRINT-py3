import os
import importlib
import sys
import subprocess
import sprint


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

    # 初始化参数（略，与原代码一致）
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
    var_limit = 20
    poly_limit = 10
    rm_multi = 0
    paired_end = False
    repeat = False
    options = []
    read2 = ''
    read1 = ''

    # 解析命令行参数（略，与原代码一致）
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

    # 提取非选项参数（略，与原代码一致）
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

    # 创建目录（略，与原代码一致）
    if not os.path.exists(output):
        os.mkdir(output)
    if not os.path.exists(tmp):
        os.mkdir(tmp)

    # 写入参数文件（略，与原代码一致）
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

        # Step 1: BWA aln 比对（无变更）
        step1_1 = subprocess.Popen(f"{bwa} aln -t {mapcpu} {refgenome} {read1} > {tmp}read1.sai", shell=True)
        if paired_end is True:
            step1_2 = subprocess.Popen(f"{bwa} aln -t {mapcpu} {refgenome} {read2} > {tmp}read2.sai", shell=True)

        step1_1.wait()
        if paired_end is True:
            step1_2.wait()

        # Step 2: BWA samse 生成SAM（无变更）
        step1_3 = subprocess.Popen(f"{bwa} samse -n4 {refgenome} {tmp}read1.sai {read1} > {tmp}name_read1.sam",
                                   shell=True)
        if paired_end is True:
            step1_4 = subprocess.Popen(f"{bwa} samse -n4 {refgenome} {tmp}read2.sai {read2} > {tmp}name_read2.sam",
                                       shell=True)

        step1_3.wait()
        if paired_end is True:
            step1_4.wait()

        # 清理临时文件（无变更）
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

            # 清理临时文件（无变更）
            if os.path.exists(tmp + 'all.bam'):
                for f in ['name_read1.sam', 'name_read1.bam', 'name_read1_sorted.bam',
                          'name_read2.sam', 'name_read2.bam', 'name_read2_sorted.bam']:
                    if os.path.exists(tmp + f):
                        os.remove(tmp + f)
        else:
            # 单端数据排序（新版语法）
            step1_9 = subprocess.Popen(f"{samtools} sort -@ {mapcpu} -o {tmp}all.bam {tmp}name_read1.bam", shell=True)
            step1_9.wait()

            # 清理临时文件（无变更）
            if os.path.exists(tmp + 'all.bam'):
                for f in ['name_read1.sam', 'name_read1.bam']:
                    if os.path.exists(tmp + f):
                        os.remove(tmp + f)

        # Step 6: 转换回SAM（无变更）
        step2_2 = subprocess.Popen(f"{samtools} view -h -o {tmp}all.sam {tmp}all.bam", shell=True)
        step2_2.wait()
        subprocess.Popen(f"cp {tmp}/all.sam {ori_tmp}/{TAG}_all.sam", shell=True).wait()

        if os.path.exists(tmp + 'all.sam'):
            os.remove(tmp + 'all.sam')

    # 后续逻辑（与原代码一致，仅调整 print 为函数形式）
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

    # 后续文件清理、SNV识别、RES识别等逻辑（与原代码一致）
    if os.path.exists(tmp + 'genome_mskAG_unmapped.sam'):
        for f in ['cut_read1.fastq', 'cut_read2.fastq', 'genome_mskAG_unmapped.fq', 'genome_mskAG_unmapped.sam',
                  'genome_mskTC_unmapped.fq', 'genome_mskTC_unmapped.sam', 'genome_unmapped.fq', 'genome_unmapped.sam',
                  'transcript_unmapped_A_to_G.fq', 'transcript_unmapped.fq', 'transcript_unmapped.sam',
                  'regular_unmapped.fq', 'regular_unmapped_A_to_G.fq']:
            if os.path.exists(tmp + f):
                os.remove(tmp + f)

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

    print('identifying SNVs...')

    if os.path.exists(refgenome + '.trans.fa'):
        sprint.mask_zz2snv(tmp + 'transcript_mskAG_all.zz.dedup', tmp + 'transcript_mskAG_all.zz.dedup.snv',
                           tmp + 'baseq.cutoff')
        sprint.mask_zz2snv(tmp + 'transcript_mskTC_all.zz.dedup', tmp + 'transcript_mskTC_all.zz.dedup.snv',
                           tmp + 'baseq.cutoff')
        sprint.mask_zz2snv(tmp + 'transcript_all.zz.dedup', tmp + 'transcript_all.zz.dedup.snv', tmp + 'baseq.cutoff')

        sprint.tzz2gzz(refgenome + '.trans.fa.loc', tmp + 'transcript_mskAG_all.zz.dedup',
                       tmp + 'transcript_mskAG_all.zz.dedup.genome.zz')
        sprint.tzz2gzz(refgenome + '.trans.fa.loc', tmp + 'transcript_mskTC_all.zz.dedup',
                       tmp + 'transcript_mskTC_all.zz.dedup.genome.zz')
        sprint.tzz2gzz(refgenome + '.trans.fa.loc', tmp + 'transcript_all.zz.dedup',
                       tmp + 'transcript_all.zz.dedup.genome.zz')

    sprint.mask_zz2snv(tmp + 'genome_mskAG_all.zz.dedup', tmp + 'genome_mskAG_all.zz.dedup.snv', tmp + 'baseq.cutoff')
    sprint.mask_zz2snv(tmp + 'genome_mskTC_all.zz.dedup', tmp + 'genome_mskTC_all.zz.dedup.snv', tmp + 'baseq.cutoff')
    sprint.mask_zz2snv(tmp + 'genome_all.zz.dedup', tmp + 'genome_all.zz.dedup.snv', tmp + 'baseq.cutoff')

    if os.path.exists(refgenome + '.trans.fa'):
        subprocess.Popen(
            f"cat {tmp}/genome_mskAG_all.zz.dedup {tmp}/genome_mskTC_all.zz.dedup {tmp}/genome_all.zz.dedup {tmp}/transcript_mskAG_all.zz.dedup.genome.zz {tmp}/transcript_mskTC_all.zz.dedup.genome.zz {tmp}/transcript_all.zz.dedup.genome.zz > {tmp}/all_combined.zz",
            shell=True).wait()
    else:
        subprocess.Popen(
            f"cat {tmp}/genome_mskAG_all.zz.dedup {tmp}/genome_mskTC_all.zz.dedup {tmp}/genome_all.zz.dedup > {tmp}/all_combined.zz",
            shell=True).wait()
    sprint.sort_zz(tmp + '/all_combined.zz', tmp + '/all_combined.zz.sorted')

    if os.path.exists(refgenome + '.trans.fa'):
        sprint.transcript_locator(tmp + 'transcript_mskAG_all.zz.dedup.snv', refgenome + '.trans.fa.loc',
                                  tmp + 'transcript_mskAG_all.zz.dedup.snv.genome.snv')
        sprint.transcript_locator(tmp + 'transcript_mskTC_all.zz.dedup.snv', refgenome + '.trans.fa.loc',
                                  tmp + 'transcript_mskTC_all.zz.dedup.snv.genome.snv')
        sprint.transcript_locator(tmp + 'transcript_all.zz.dedup.snv', refgenome + '.trans.fa.loc',
                                  tmp + 'transcript_all.zz.dedup.snv.genome.snv')

        sprint.transcript_sort(tmp + 'transcript_all.zz.dedup.snv.genome.snv',
                               tmp + 'transcript_all.zz.dedup.snv.genome.snv.sort')
        sprint.transcript_sort(tmp + 'transcript_mskTC_all.zz.dedup.snv.genome.snv',
                               tmp + 'transcript_mskTC_all.zz.dedup.snv.genome.snv.sort')
        sprint.transcript_sort(tmp + 'transcript_mskAG_all.zz.dedup.snv.genome.snv',
                               tmp + 'transcript_mskAG_all.zz.dedup.snv.genome.snv.sort')

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

    print('identifying RESs...')

    if repeat is not False:
        sprint.annotate(tmp + 'regular.snv', repeat, tmp + 'regular.snv.anno')
        sprint.seperate(tmp + 'regular.snv.anno', tmp + 'regular.snv.anno.alu', tmp + 'regular.snv.anno.nalurp',
                        tmp + 'regular.snv.anno.nrp', 'Alu')
        sprint.get_snv_with_ad(tmp + 'regular.snv.anno.alu', tmp + 'regular.snv.anno.alu.ad2', 2)
        sprint.snv_cluster(tmp + 'regular.snv.anno.alu', tmp + 'regular_alu.res.ad1', cluster_distance,
                           cluster_size_alu_ad1)
        sprint.snv_cluster(tmp + 'regular.snv.anno.alu.ad2', tmp + 'regular_alu.res.ad2', cluster_distance,
                           cluster_size_alu_ad2)
        sprint.bed_or(tmp + 'regular_alu.res.ad1', tmp + 'regular_alu.res.ad2', tmp + 'regular_alu.res')
        sprint.snv_cluster(tmp + 'regular.snv.anno.nalurp', tmp + 'regular_nalurp.res', cluster_distance,
                           cluster_size_nalurp)
        sprint.snv_cluster(tmp + 'regular.snv.anno.nrp', tmp + 'regular_nrp.res', cluster_distance, cluster_size_nrp)
        sprint.combine_res(tmp + 'regular_alu.res', tmp + 'regular_nalurp.res', tmp + 'regular_nrp.res',
                           tmp + 'regular_split.res')
        cluster_size_regular_max = max(
            [cluster_size_alu_ad1, cluster_size_alu_ad2, cluster_size_nalurp, cluster_size_nrp])
        sprint.combine_res(tmp + 'regular.snv.anno.alu', tmp + 'regular.snv.anno.nalurp', tmp + 'regular.snv.anno.nrp',
                           tmp + 'regular.snv.anno.rmsrp')
        sprint.snv_cluster(tmp + 'regular.snv.anno.rmsrp', tmp + 'regular_overall.res', cluster_distance,
                           cluster_size_regular_max)
        sprint.res_or(tmp + 'regular_split.res', tmp + 'regular_overall.res', tmp + 'regular.res')

        sprint.annotate(tmp + 'hyper_mskTC.snv', repeat, tmp + 'hyper_mskTC.snv.anno')
        sprint.seperate(tmp + 'hyper_mskTC.snv.anno', tmp + 'hyper_mskTC.snv.anno.alu',
                        tmp + 'hyper_mskTC.snv.anno.nalurp', tmp + 'hyper_mskTC.snv.anno.nrp', 'Alu')
        sprint.snv_cluster(tmp + 'hyper_mskTC.snv.anno.alu', tmp + 'hyper_mskTC_alu.res', cluster_distance,
                           cluster_size_alu_hp)
        sprint.snv_cluster(tmp + 'hyper_mskTC.snv.anno.nalurp', tmp + 'hyper_mskTC_nalurp.res', cluster_distance,
                           cluster_size_nalurp_hp)
        sprint.snv_cluster(tmp + 'hyper_mskTC.snv.anno.nrp', tmp + 'hyper_mskTC_nrp.res', cluster_distance,
                           cluster_size_nrp_hp)
        sprint.combine_res(tmp + 'hyper_mskTC_alu.res', tmp + 'hyper_mskTC_nalurp.res', tmp + 'hyper_mskTC_nrp.res',
                           tmp + 'hyper_mskTC_split.res')
        cluster_size_hyper_max = max([cluster_size_alu_hp, cluster_size_nalurp_hp, cluster_size_nrp_hp])
        sprint.combine_res(tmp + 'hyper_mskTC.snv.anno.alu', tmp + 'hyper_mskTC.snv.anno.nalurp',
                           tmp + 'hyper_mskTC.snv.anno.nrp', tmp + 'hyper_mskTC.snv.anno.rmsrp')
        sprint.snv_cluster(tmp + 'hyper_mskTC.snv.anno.rmsrp', tmp + 'hyper_mskTC_overall.res', cluster_distance,
                           cluster_size_hyper_max)
        sprint.res_or(tmp + 'hyper_mskTC_split.res', tmp + 'hyper_mskTC_overall.res', tmp + 'hyper_mskTC.res')

        sprint.annotate(tmp + 'hyper_mskAG.snv', repeat, tmp + 'hyper_mskAG.snv.anno')
        sprint.seperate(tmp + 'hyper_mskAG.snv.anno', tmp + 'hyper_mskAG.snv.anno.alu',
                        tmp + 'hyper_mskAG.snv.anno.nalurp', tmp + 'hyper_mskAG.snv.anno.nrp', 'Alu')
        sprint.snv_cluster(tmp + 'hyper_mskAG.snv.anno.alu', tmp + 'hyper_mskAG_alu.res', cluster_distance,
                           cluster_size_alu_hp)
        sprint.snv_cluster(tmp + 'hyper_mskAG.snv.anno.nalurp', tmp + 'hyper_mskAG_nalurp.res', cluster_distance,
                           cluster_size_nalurp_hp)
        sprint.snv_cluster(tmp + 'hyper_mskAG.snv.anno.nrp', tmp + 'hyper_mskAG_nrp.res', cluster_distance,
                           cluster_size_nrp_hp)
        sprint.combine_res(tmp + 'hyper_mskAG_alu.res', tmp + 'hyper_mskAG_nalurp.res', tmp + 'hyper_mskAG_nrp.res',
                           tmp + 'hyper_mskAG_split.res')
        cluster_size_hyper_max = max([cluster_size_alu_hp, cluster_size_nalurp_hp, cluster_size_nrp_hp])
        sprint.combine_res(tmp + 'hyper_mskAG.snv.anno.alu', tmp + 'hyper_mskAG.snv.anno.nalurp',
                           tmp + 'hyper_mskAG.snv.anno.nrp', tmp + 'hyper_mskAG.snv.anno.rmsrp')
        sprint.snv_cluster(tmp + 'hyper_mskAG.snv.anno.rmsrp', tmp + 'hyper_mskAG_overall.res', cluster_distance,
                           cluster_size_hyper_max)
        sprint.res_or(tmp + 'hyper_mskAG_split.res', tmp + 'hyper_mskAG_overall.res', tmp + 'hyper_mskAG.res')

        sprint.snv_or(tmp + 'hyper_mskTC.res', tmp + 'hyper_mskAG.res', tmp + 'hyper.res')
    else:
        sprint.snv_cluster(tmp + 'regular.snv', tmp + 'regular.res_tmp', cluster_distance, cluster_size_rg)
        sprint.o2b(tmp + 'regular.res_tmp', tmp + 'regular.res')

        sprint.snv_cluster(tmp + 'hyper_mskTC.snv', tmp + 'hyper_mskTC.res', cluster_distance, cluster_size_hp)
        sprint.snv_cluster(tmp + 'hyper_mskAG.snv', tmp + 'hyper_mskAG.res', cluster_distance, cluster_size_hp)
        sprint.snv_or(tmp + 'hyper_mskTC.res', tmp + 'hyper_mskAG.res', tmp + 'hyper.res')

    try:
        subprocess.Popen(f"rm -rf {tmp}/*.anno.*", shell=True).wait()
    except Exception as e:
        pass

    # 输出结果文件
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