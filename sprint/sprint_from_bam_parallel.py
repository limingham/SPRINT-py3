import os
import importlib
import sys
import subprocess
import sprint
import multiprocessing
from multiprocessing import Pool, get_context
import shutil
from tqdm import tqdm
import traceback


# ==================== 移到全局作用域的核心函数 ====================
def process_chromosome(params):
    """
    处理单个染色体的BAM文件（全局函数，解决pickle问题）
    :param params: 包含所有必要参数的字典
    """
    chromosome = params['chromosome']
    bam = params['bam']
    samtools = params['samtools']
    tmp = params['tmp']
    refgenome = params['refgenome']
    repeat = params['repeat']
    cluster_distance = params['cluster_distance']
    cluster_size_alu_ad1 = params['cluster_size_alu_ad1']
    cluster_size_alu_ad2 = params['cluster_size_alu_ad2']
    cluster_size_nalurp = params['cluster_size_nalurp']
    cluster_size_nrp = params['cluster_size_nrp']
    cluster_size_rg = params['cluster_size_rg']
    mapcpu = params['mapcpu']
    chr_tmp = params['chr_tmp']

    error_msg = None
    try:
        # 创建该染色体的临时目录
        chr_process_tmp = tmp + f'chr_{chromosome}/'
        if not os.path.exists(chr_process_tmp):
            os.makedirs(chr_process_tmp)

        # 拆分出该染色体的BAM文件
        chr_bam = f"{chr_tmp}/{chromosome}.bam"
        if not os.path.exists(chr_bam):
            split_cmd = f"{samtools} view -@ {str(mapcpu)} -b {bam} {chromosome} > {chr_bam}"
            split_process = subprocess.Popen(split_cmd, shell=True, stderr=subprocess.PIPE)
            split_process.wait()
            if split_process.returncode != 0:
                raise Exception(f"Split BAM failed: {split_process.stderr.read().decode()}")

            # 建立索引
            index_cmd = f"{samtools} index {chr_bam}"
            index_process = subprocess.Popen(index_cmd, shell=True, stderr=subprocess.PIPE)
            index_process.wait()
            if index_process.returncode != 0:
                raise Exception(f"Index BAM failed: {index_process.stderr.read().decode()}")

        # BAM处理 - 单个染色体
        view_cmd = f"{samtools} view -@ {str(mapcpu)} {chr_bam} -o {chr_process_tmp}/aligned.sam"
        view_process = subprocess.Popen(view_cmd, shell=True, stderr=subprocess.PIPE)
        view_process.wait()
        if view_process.returncode != 0:
            raise Exception(f"Convert BAM to SAM failed: {view_process.stderr.read().decode()}")

        sprint.change_sam_read_name(chr_process_tmp + '/aligned.sam', chr_process_tmp + 'genome_all.sam', 'read1')
        sprint.sam2zz(chr_process_tmp + 'genome_all.sam', refgenome, chr_process_tmp + 'genome_all.zz')
        sprint.dedup(chr_process_tmp + 'genome_all.zz', chr_process_tmp + 'genome_all.zz.dedup')
        if os.path.exists(chr_process_tmp + 'genome_all.zz'):
            os.remove(chr_process_tmp + 'genome_all.zz')

        # 确定baseq cutoff
        fi = open(chr_process_tmp + '/genome_all.sam')
        fo = open(chr_process_tmp + '/baseq.cutoff', 'w')
        did = 0
        for line in fi:
            if did == 1:
                break
            if line.startswith('@'):
                continue
            qua = line.split('\t')[10]
            for i in qua:
                if ord(i) > 76:
                    fo.write('89')
                    did = 1
                    break
                if ord(i) < 60:
                    fo.write('58')
                    did = 1
                    break
        fi.close()
        fo.close()

        sprint.zz2snv(chr_process_tmp + 'genome_all.zz.dedup', chr_process_tmp + 'genome_all.zz.dedup.snv',
                      chr_process_tmp + 'baseq.cutoff')
        subprocess.Popen(f"cp {chr_process_tmp}/genome_all.zz.dedup.snv {chr_process_tmp}/regular.snv",
                         shell=True).wait()

        # 识别RESs
        if repeat is not False:
            sprint.annotate(chr_process_tmp + 'regular.snv', repeat, chr_process_tmp + 'regular.snv.anno')
            sprint.seperate(chr_process_tmp + 'regular.snv.anno', chr_process_tmp + 'regular.snv.anno.alu',
                            chr_process_tmp + 'regular.snv.anno.nalurp', chr_process_tmp + 'regular.snv.anno.nrp',
                            'Alu')
            sprint.get_snv_with_ad(chr_process_tmp + 'regular.snv.anno.alu',
                                   chr_process_tmp + 'regular.snv.anno.alu.ad2', 2)
            sprint.snv_cluster(chr_process_tmp + 'regular.snv.anno.alu', chr_process_tmp + 'regular_alu.res.ad1',
                               cluster_distance, cluster_size_alu_ad1)
            sprint.snv_cluster(chr_process_tmp + 'regular.snv.anno.alu.ad2',
                               chr_process_tmp + 'regular_alu.res.ad2',
                               cluster_distance, cluster_size_alu_ad2)
            sprint.bed_or(chr_process_tmp + 'regular_alu.res.ad1', chr_process_tmp + 'regular_alu.res.ad2',
                          chr_process_tmp + 'regular_alu.res')
            sprint.snv_cluster(chr_process_tmp + 'regular.snv.anno.nalurp', chr_process_tmp + 'regular_nalurp.res',
                               cluster_distance, cluster_size_nalurp)
            sprint.snv_cluster(chr_process_tmp + 'regular.snv.anno.nrp', chr_process_tmp + 'regular_nrp.res',
                               cluster_distance, cluster_size_nrp)
            sprint.combine_res(chr_process_tmp + 'regular_alu.res', chr_process_tmp + 'regular_nalurp.res',
                               chr_process_tmp + 'regular_nrp.res', chr_process_tmp + 'regular_split.res')
            cluster_size_regular_max = max(
                [cluster_size_alu_ad1, cluster_size_alu_ad2, cluster_size_nalurp, cluster_size_nrp])
            sprint.combine_res(chr_process_tmp + 'regular.snv.anno.alu',
                               chr_process_tmp + 'regular.snv.anno.nalurp',
                               chr_process_tmp + 'regular.snv.anno.nrp', chr_process_tmp + 'regular.snv.anno.rmsrp')
            sprint.snv_cluster(chr_process_tmp + 'regular.snv.anno.rmsrp', chr_process_tmp + 'regular_overall.res',
                               cluster_distance, cluster_size_regular_max)
            sprint.res_or(chr_process_tmp + 'regular_split.res', chr_process_tmp + 'regular_overall.res',
                          chr_process_tmp + 'regular.res')
        else:
            sprint.snv_cluster(chr_process_tmp + 'regular.snv', chr_process_tmp + 'regular.res_tmp',
                               cluster_distance, cluster_size_rg)
            sprint.o2b(chr_process_tmp + 'regular.res_tmp', chr_process_tmp + 'regular.res')

        # 清理注释临时文件
        try:
            subprocess.Popen(f"rm -rf {chr_process_tmp}/*.anno.*", shell=True).wait()
        except Exception as e:
            pass

        sprint.get_depth(chr_process_tmp + '/genome_all.zz.dedup', chr_process_tmp + '/regular.res',
                         chr_process_tmp + '/regular.res.depth')

        # 保存该染色体的结果
        chr_result = f"{tmp}/chr_{chromosome}_results.depth"
        cp_process = subprocess.Popen(f"cp {chr_process_tmp}/regular.res.depth {chr_result}", shell=True,
                                      stderr=subprocess.PIPE)
        cp_process.wait()
        if cp_process.returncode != 0:
            raise Exception(f"Copy result file failed: {cp_process.stderr.read().decode()}")

        return {
            'chromosome': chromosome,
            'result_file': chr_result,
            'error': None
        }

    except Exception as e:
        # 记录详细错误信息
        error_msg = f"Error processing {chromosome}: {str(e)}\n{traceback.format_exc()}"
        return {
            'chromosome': chromosome,
            'result_file': None,
            'error': error_msg
        }


# ==================== 原main函数（移除内部的process_chromosome） ====================
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
        print("")
        print("")
        print("   Usage:")
        print("")
        print(
            "      sprint_from_bam   [options]  alinged_reads(.bam)   reference_genome(.fa)   output_path   samtools_path")
        print("")
        print("      options:")
        print(
            "         -rp      repeat_file      # Optional, you can download it from http://sprint.software/SPRINT/dbrep/")
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
        print("         -t       INT              # Number of threads for parallel processing (default is 4)")
        print("")
        print("   Example:")
        print("")
        print(
            "       sprint_from_bam -rp hg38_repeat.txt -t 8 aligned_reads.bam hg38.fa output ./samtools-1.21/samtools")
        print("")
        print("")
        sys.exit(0)

    if len(sys.argv) < 2:
        help_doc()

    # 初始化参数
    n_threads = 4
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
                exit()
        elif sys.argv[i] == '-rp':
            try:
                repeat = sys.argv[i + 1]
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                exit()
        elif sys.argv[i] == '-cd':
            try:
                cluster_distance = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                exit()
        elif sys.argv[i] == '-csad1':
            try:
                cluster_size_alu_ad1 = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                exit()
        elif sys.argv[i] == '-csad2':
            try:
                cluster_size_alu_ad2 = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                exit()
        elif sys.argv[i] == '-csnar':
            try:
                cluster_size_nalurp = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                exit()
        elif sys.argv[i] == '-csnr':
            try:
                cluster_size_nrp = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                exit()
        elif sys.argv[i] == '-csrg':
            try:
                cluster_size_rg = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                exit()
        elif sys.argv[i] == '-t':
            try:
                n_threads = int(sys.argv[i + 1])
                options.append(i)
                options.append(i + 1)
            except Exception as e:
                print('options error!')
                help_doc()
                exit()

        i += 1

    all_argv = []
    i = 1
    while i < len(sys.argv):
        if i not in options:
            all_argv.append(i)
        i = i + 1

    if len(all_argv) != 4:
        help_doc()
        exit()

    refgenome = sys.argv[all_argv[1]]
    output = sys.argv[all_argv[2]] + '/'
    tmp = output + '/tmp/'
    chr_tmp = tmp + 'chr_split/'

    bam = sys.argv[all_argv[0]]
    samtools = sys.argv[all_argv[3]]

    # 创建目录
    if not os.path.exists(output):
        os.mkdir(output)
    if not os.path.exists(tmp):
        os.mkdir(tmp)
    if not os.path.exists(chr_tmp):
        os.makedirs(chr_tmp)

    # 保存参数文件
    frc = open(tmp + 'PARAMETER.txt', 'w')
    frc.write(sys.argv[0])
    for one in sys.argv[1:]:
        frc.write('   ' + one)
    frc.write('\n')
    frc.close()

    # fq2sam函数（保留在main内部，因为未被多进程调用）
    def fq2sam(TAG, paired_end, read1, read2, tmp, refgenome, bwa, samtools, mapcpu, read_format):
        if paired_end is True:
            mapcpu = max([int(int(mapcpu) / 2.0), 1])
        ori_tmp = tmp
        tmp = tmp + '/' + TAG + '/'
        if not os.path.exists(tmp):
            os.mkdir(tmp)

        step1_1_cmd = f"{bwa} aln -t {str(mapcpu)} {refgenome} {read1} > {tmp}read1.sai"
        step1_1 = subprocess.Popen(step1_1_cmd, shell=True)
        if paired_end is True:
            step1_2_cmd = f"{bwa} aln -t {str(mapcpu)} {refgenome} {read2} > {tmp}read2.sai"
            step1_2 = subprocess.Popen(step1_2_cmd, shell=True)

        step1_1.wait()
        if paired_end is True:
            step1_2.wait()

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

        # samtools view
        step1_7_cmd = f"{samtools} view -@ {str(mapcpu)} -bS {tmp}name_read1.sam -o {tmp}name_read1.bam"
        step1_7 = subprocess.Popen(step1_7_cmd, shell=True)
        if paired_end is True:
            step1_8_cmd = f"{samtools} view -@ {str(mapcpu)} -bS {tmp}name_read2.sam -o {tmp}name_read2.bam"
            step1_8 = subprocess.Popen(step1_8_cmd, shell=True)
        step1_7.wait()
        if paired_end is True:
            step1_8.wait()

        # samtools sort
        if paired_end is True:
            step1_9_cmd = f"{samtools} sort -@ {str(mapcpu)} {tmp}name_read1.bam -o {tmp}name_read1_sorted.bam"
            step1_9 = subprocess.Popen(step1_9_cmd, shell=True)
            step1_10_cmd = f"{samtools} sort -@ {str(mapcpu)} {tmp}name_read2.bam -o {tmp}name_read2_sorted.bam"
            step1_10 = subprocess.Popen(step1_10_cmd, shell=True)
            step1_9.wait()
            step1_10.wait()

            # samtools merge
            step1_11_cmd = f"{samtools} merge -f -o {tmp}all.bam {tmp}name_read1_sorted.bam {tmp}name_read2_sorted.bam"
            step1_11 = subprocess.Popen(step1_11_cmd, shell=True)
            step1_11.wait()

            # 清理临时文件
            if os.path.exists(tmp + 'all.bam'):
                for f in ['name_read1.sam', 'name_read1.bam', 'name_read1_sorted.bam',
                          'name_read2.sam', 'name_read2.bam', 'name_read2_sorted.bam']:
                    if os.path.exists(tmp + f):
                        os.remove(tmp + f)
        else:
            step1_9_cmd = f"{samtools} sort -@ {str(mapcpu)} {tmp}name_read1.bam -o {tmp}all.bam"
            step1_9 = subprocess.Popen(step1_9_cmd, shell=True)
            step1_9.wait()
            # 清理临时文件
            if os.path.exists(tmp + 'all.bam'):
                for f in ['name_read1.sam', 'name_read1.bam']:
                    if os.path.exists(tmp + f):
                        os.remove(tmp + f)

        # samtools view转sam
        step2_2_cmd = f"{samtools} view -@ {str(mapcpu)} -h {tmp}all.bam -o {tmp}all.sam"
        step2_2 = subprocess.Popen(step2_2_cmd, shell=True)
        step2_2.wait()

        # 复制sam文件
        cp_cmd = f"cp {tmp}/all.sam {ori_tmp}/{TAG}_all.sam"
        subprocess.Popen(cp_cmd, shell=True).wait()
        if os.path.exists(tmp + 'all.sam'):
            os.remove(tmp + 'all.sam')

    # 核心并行处理逻辑
    try:
        print("Getting chromosome list from BAM file...")
        # 获取染色体列表
        idxstats_cmd = f"{samtools} idxstats {bam} | cut -f1 | grep -v '^$' | grep -v '*'"
        result = subprocess.check_output(idxstats_cmd, shell=True, text=True)
        chromosomes = [chr.strip() for chr in result.split('\n') if chr.strip()]

        print(f"Found {len(chromosomes)} chromosomes: {chromosomes}")
        print(f"Starting parallel processing with {n_threads} threads...")

        # 准备每个染色体的参数
        process_params = []
        for chr_name in chromosomes:
            params = {
                'chromosome': chr_name,
                'bam': bam,
                'samtools': samtools,
                'tmp': tmp,
                'refgenome': refgenome,
                'repeat': repeat,
                'cluster_distance': cluster_distance,
                'cluster_size_alu_ad1': cluster_size_alu_ad1,
                'cluster_size_alu_ad2': cluster_size_alu_ad2,
                'cluster_size_nalurp': cluster_size_nalurp,
                'cluster_size_nrp': cluster_size_nrp,
                'cluster_size_rg': cluster_size_rg,
                'mapcpu': mapcpu,
                'chr_tmp': chr_tmp
            }
            process_params.append(params)

        # 创建进度条
        pbar = tqdm(
            total=len(chromosomes),
            desc="Processing chromosomes",
            unit="chr",
            ncols=80,
            dynamic_ncols=True,
            bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"
        )

        # 跨平台进程池
        if sys.platform == 'win32':
            ctx = get_context('spawn')
        else:
            ctx = get_context('fork')

        # 并行处理
        with ctx.Pool(processes=n_threads) as pool:
            results = []
            for result in pool.imap_unordered(process_chromosome, process_params):
                results.append(result)
                pbar.update(1)

        pbar.close()

        # 收集结果
        all_chr_results = []
        failed_chromosomes = []
        error_logs = []

        for res in results:
            chr_name = res['chromosome']
            if res['error'] is not None:
                failed_chromosomes.append(chr_name)
                error_logs.append(res['error'])
            elif res['result_file'] and os.path.exists(res['result_file']):
                all_chr_results.append(res['result_file'])

        # 输出统计信息
        if failed_chromosomes:
            print(f"\n⚠️  Warning: {len(failed_chromosomes)} chromosomes failed to process:")
            print(f"   Failed chromosomes: {failed_chromosomes}")
            error_log_file = f"{tmp}/processing_errors.log"
            with open(error_log_file, 'w') as f:
                f.write('\n'.join(error_logs))
            print(f"   Detailed error logs saved to: {error_log_file}")

        print(f"\n📊 Processing summary:")
        print(f"   Total chromosomes: {len(chromosomes)}")
        print(f"   Successfully processed: {len(all_chr_results)}")
        print(f"   Failed: {len(failed_chromosomes)}")

        # 合并结果
        print("\nMerging results from all chromosomes...")
        if all_chr_results:
            merge_cmd = f"cat {' '.join(all_chr_results)} > {tmp}/all_chromosomes.res.depth"
            merge_process = subprocess.Popen(merge_cmd, shell=True, stderr=subprocess.PIPE)
            merge_process.wait()
            if merge_process.returncode != 0:
                raise Exception(f"Merge results failed: {merge_process.stderr.read().decode()}")

            # 添加表头
            cat_cmd = f"echo \"#Chrom\tStart(0base)\tEnd(1base)\tType\tSupporting_reads\tStrand\tAD:DP\" | cat - {tmp}/all_chromosomes.res.depth > {output}/SPRINT_identified_regular.res"
            cat_process = subprocess.Popen(cat_cmd, shell=True, stderr=subprocess.PIPE)
            cat_process.wait()
            if cat_process.returncode != 0:
                raise Exception(f"Add header failed: {cat_process.stderr.read().decode()}")

            print(f"\n✅ Results merged successfully!")
            print(f"   Final output file: {output}/SPRINT_identified_regular.res")
        else:
            print("\n❌ Error: No results were generated from any chromosome!")
            sys.exit(1)

        # 复制参数文件
        subprocess.Popen(f"cp {tmp}/PARAMETER.txt {output}/PARAMETER.txt", shell=True).wait()

        print('\n🎉 All chromosomes processed and results merged!')
        print('finished !')
        sys.exit(0)

    except ImportError:
        print("\n❌ Error: tqdm library not found!")
        print("   Please install it first: pip install tqdm")
        sys.exit(1)
    except Exception as e:
        print('')
        print('❌ ERROR!')
        print('')
        print(e)
        print('')
        help_doc()


# ==================== 程序入口 ====================
if __name__ == "__main__":
    main()