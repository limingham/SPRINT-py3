def mismatch_num(seqlen):
    if seqlen < 17:
        return 1
    elif seqlen < 38:
        return 2
    elif seqlen < 64:
        return 3
    elif seqlen < 93:
        return 4
    elif seqlen < 124:
        return 5
    elif seqlen < 157:
        return 6
    elif seqlen < 190:
        return 7
    elif seqlen < 225:
        return 8
    else:
        return 9


def mask_zz2snv(zz_in_dir=0, bed_out_dir=0, baseq_cutoff_dir=0):
    # 使用Python3推荐的with上下文管理器管理文件，自动关闭无需手动close()
    with open(zz_in_dir, 'r') as fi, \
         open(bed_out_dir, 'w') as fo, \
         open(baseq_cutoff_dir, 'r') as fqua:

        # 读取碱基质量阈值，兼容无换行符场景
        limitbasequa_line = fqua.readline().rstrip('\n')
        limitbasequa = int(limitbasequa_line)

        limitad = 1
        limitloc = 5
        limitmpqua = 0
        allsnv = {}

        for line in fi:
            truesnv = []
            line_stripped = line.rstrip()
            if not line_stripped:  # 跳过空行，提升健壮性
                continue
            seq = line_stripped.split('\t')
            mismatch = seq[4].split(';')
            basequa = seq[5].split(',')
            loc = seq[9].split(',')  # fragment-loc
            mpqua = int(seq[2])

            ####################################################################
            # change the sam flag 'seq[1]' when you didn't use "bwa -aln" as mapper
            seq[1] = int(seq[1])
            if len(bin(seq[1])) >= 7:
                if bin(seq[1])[-3] != '1':
                    if bin(seq[1])[-5] == '1':
                        seq[1] = '16'
                    else:
                        seq[1] = '0'
            elif len(bin(seq[1])) >= 5:
                if bin(seq[1])[-3] != '1':
                    seq[1] = '0'
            else:
                seq[1] = '0'
            #####################################################################

            if basequa[0] != '*' and mpqua >= limitmpqua and mpqua < 200:
                i = 0
                baseqlst = []
                mistype = {}
                while i < len(basequa):
                    # 增加异常捕获，避免非整数碱基质量值导致程序崩溃
                    try:
                        baseq_val = int(basequa[i])
                        mismatch_key = mismatch[i].split(':')[0]
                    except (ValueError, IndexError):
                        i += 1
                        continue

                    baseqlst.append(baseq_val)
                    # Python2→Python3 修复：异常捕获使用as关键字
                    try:
                        mistype[mismatch_key] += 1
                    except Exception as e:
                        mistype[mismatch_key] = 1

                    # 校验loc索引和转换，避免异常
                    try:
                        loc_val = int(loc[i])
                    except (ValueError, IndexError):
                        i += 1
                        continue

                    if baseq_val >= limitbasequa and loc_val > limitloc:
                        truesnv.append([mismatch[i], seq[1], seq[8]])

                    i += 1  # 简化自增，符合Python编码习惯

                miss = []
                for mis in mistype:
                    miss.append(mistype[mis])
                miss.sort()
                if len(miss) >= 2:
                    missnum = sum(miss[:-1])
                else:
                    missnum = 0

                # 计算平均碱基质量，避免除零异常
                if len(baseqlst) > 0:
                    avg_baseq = sum(baseqlst) / float(len(baseqlst))
                    # 注释部分保留原逻辑，如需启用可直接取消注释
                    # if avg_baseq >= limitbasequa and missnum <= mismatch_num(len(seq[7])):
                    if avg_baseq >= limitbasequa:
                        for snv in truesnv:
                            # 构造snv唯一键，使用f-string优化拼接（可选，保留原拼接方式兼容逻辑）
                            snv_key = f"{seq[0]}\t{snv[0].split(':')[0]}\t{snv[0].split(':')[1]}"
                            # 简化seq[1]整数转换的重复操作
                            seq1_int = int(seq[1])
                            # 条件判断保留原逻辑，优化格式提升可读性
                            condition1 = (len(bin(seq1_int)) > 5 and bin(seq1_int)[-5] == '1' and snv[2][-1] == '1')
                            condition2 = ((len(bin(seq1_int)) < 5 or bin(seq1_int)[-5] == '0') and snv[2][-1] == '2')
                            condition3 = (len(bin(seq1_int)) > 5 and bin(seq1_int)[-5] == '1' and snv[2][-1] == '2')
                            condition4 = ((len(bin(seq1_int)) < 5 or bin(seq1_int)[-5] == '0') and snv[2][-1] == '1')

                            try:
                                allsnv[snv_key][0] += 1
                                if condition1 or condition2:
                                    allsnv[snv_key][1] += 1
                                elif condition3 or condition4:
                                    allsnv[snv_key][2] += 1
                            except Exception as e:  # Python2→Python3 修复：异常捕获使用as关键字
                                if condition1 or condition2:
                                    allsnv[snv_key] = [1, 1, 0]
                                elif condition3 or condition4:
                                    allsnv[snv_key] = [1, 0, 1]

        # 构造SNV BED数据并排序
        snv_bed = []
        for snv in allsnv:
            seq_snv = snv.split('\t')
            # 索引校验，避免拆分后列数不足
            if len(seq_snv) < 3:
                continue
            snv_count = allsnv[snv]
            if snv_count[0] >= limitad:
                try:
                    loc_int = int(seq_snv[2])
                except ValueError:
                    continue
                if snv_count[1] > snv_count[2]:
                    snv_bed.append([seq_snv[0], loc_int, seq_snv[1], '+', snv_count[0]])
                elif snv_count[2] > snv_count[1]:
                    snv_bed.append([seq_snv[0], loc_int, seq_snv[1], '-', snv_count[0]])
                else:
                    snv_bed.append([seq_snv[0], loc_int, seq_snv[1], '.', snv_count[0]])

        # 排序并写入输出文件
        snv_bed.sort()
        for one in snv_bed:
            # 使用f-string优化字符串拼接，简洁高效
            output_line = f"{one[0]}\t{one[1]-1}\t{one[1]}\t{one[2]}\t{one[4]}\t{one[3]}\n"
            fo.write(output_line)
