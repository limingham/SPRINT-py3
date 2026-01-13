
############################
import subprocess
import os
import sys
import sprint


def prepare():
    # 移除冗余注释，保持代码整洁
    # Python3 核心修复：print 语句改为带括号的函数形式
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
        print("   Usage:")
        print("")
        print("      sprint prepare   [options]   reference_genome(.fa)   bwa_path")
        print("")
        print("      options:")
        print("         -t transcript_annotation(.gtf) #Optional")
        print("")
        print("   Example:")
        print("")
        print("      sprint prepare -t hg38.gtf hg38.fa ./bwa-0.7.19/bwa")
        print("")
        print("")
        sys.exit(0)

    # 修复参数判断逻辑，保持原功能不变
    if len(sys.argv) < 3:
        help_doc()

    gtf_file = False
    options = []
    read2 = ''
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == '-t':
            paired_end = True  # 保留原变量定义，不影响功能
            try:
                gtf_file = sys.argv[i+1]
                options.append(i)
                options.append(i+1)
            except Exception as e:  # Python2→Python3：异常捕获使用 as 关键字
                print('options error!')
                help_doc()
        i += 1  # 简化自增操作，符合 Python3 编码习惯

    all_argv = []
    i = 1
    while i < len(sys.argv):
        if i not in options:
            all_argv.append(i)
        i += 1

    try:
        refgenome = sys.argv[all_argv[0]]
        bwa = sys.argv[all_argv[1]]

        print('Masking A with G in reference genome...')  # Python3：print 改为函数
        sprint.maskAwithG(refgenome, refgenome + '.mskAG.fa')
        print('Masking T with C in reference genome...')
        sprint.maskTwithC(refgenome, refgenome + '.mskTC.fa')

        if gtf_file is not False:  # 优化布尔判断，更简洁易读
            print('Assembling transcripts...')
            sprint.transcript_assembler(refgenome, gtf_file, refgenome + '.trans.fa')
            print('Masking A with G in transcripts...')
            sprint.maskAwithG(refgenome + '.trans.fa', refgenome + '.trans.fa.mskAG.fa')
            print('Masking T with C in transcripts...')
            sprint.maskTwithC(refgenome + '.trans.fa', refgenome + '.trans.fa.mskTC.fa')

        print('Building BWA index...')
        # subprocess 调用逻辑完全兼容 Python3，无需修改
        step1 = subprocess.Popen(bwa + ' index -a bwtsw ' + refgenome, shell=True)
        step2 = subprocess.Popen(bwa + ' index -a bwtsw ' + refgenome + '.mskAG.fa', shell=True)
        step3 = subprocess.Popen(bwa + ' index -a bwtsw ' + refgenome + '.mskTC.fa', shell=True)
        if gtf_file is not False:
            step4 = subprocess.Popen(bwa + ' index -a bwtsw ' + refgenome + '.trans.fa', shell=True)
            step5 = subprocess.Popen(bwa + ' index -a bwtsw ' + refgenome + '.trans.fa.mskAG.fa', shell=True)
            step6 = subprocess.Popen(bwa + ' index -a bwtsw ' + refgenome + '.trans.fa.mskTC.fa', shell=True)
            step4.wait()
            step5.wait()
            step6.wait()
        step1.wait()
        step2.wait()
        step3.wait()
        sys.exit(0)
    except Exception as e:  # Python2→Python3：异常捕获使用 as 关键字
        print("ERROR!")
        print("")
        print(e)  # Python3：print 改为函数，输出异常信息
        help_doc()

# Python3 脚本规范：启用主程序入口，支持直接执行脚本
if __name__ == "__main__":
    prepare()