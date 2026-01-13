def zz2snv(zz_in_dir=0, bed_out_dir=0, baseq_cutoff_dir=0):
    # 使用Python3推荐的with上下文管理器，同时管理3个文件，自动关闭无需手动close()
    with open(zz_in_dir, 'r') as fi, \
         open(bed_out_dir, 'w') as fo, \
         open(baseq_cutoff_dir, 'r') as fqua:

        # 读取碱基质量阈值，兼容最后一行无换行符的场景
        limitbasequa_line = fqua.readline().rstrip('\n')
        limitbasequa = int(limitbasequa_line)

        limitad = 1  # 2
        limitloc = 5
        limitmpqua = 20
        allsnv = {}

        for line in fi:
            truesnv = []
            # 替代line[0:-1]，精准去除换行符，兼容无换行符场景
            line_stripped = line.rstrip('\n')
            if not line_stripped:  # 跳过空行，提升健壮性
                continue
            seq = line_stripped.split('\t')
            # 列数校验，避免后续索引越界
            if len(seq) < 10:
                continue

            # 拆分各类信息，增加异常防护
            try:
                mismatch = seq[4].split(';')
                basequa = seq[5].split(',')
                loc = seq[9].split(',')  # fragment-loc
                mpqua = int(seq[2])
            except (IndexError, ValueError):
                continue

            ####################################################################
            # change the sam flag 'seq[1]' when you didn't use "bwa -aln" as mapper
            try:
                seq1_int = int(seq[1])
            except ValueError:
                seq[1] = '0'
                seq1_int = 0  # 初始化避免后续报错

            if len(bin(seq1_int)) >= 7:
                if bin(seq1_int)[-3] != '1':
                    if bin(seq1_int)[-5] == '1':
                        seq[1] = '16'
                    else:
                        seq[1] = '0'
            elif len(bin(seq1_int)) >= 5:
                if bin(seq1_int)[-3] != '1':
                    seq[1] = '0'
            else:
                seq[1] = '0'
            #####################################################################

            # 核心筛选条件，保留原逻辑
            if basequa[0] != '*' and mpqua >= limitmpqua and mpqua < 200:
                i = 0
                baseqlst = []
                while i < len(basequa):
                    # 捕获碱基质量/位置转换异常，避免程序崩溃
                    try:
                        baseq_val = int(basequa[i])
                        loc_val = int(loc[i])
                    except (ValueError, IndexError):
                        i += 1
                        continue

                    baseqlst.append(baseq_val)
                    # 筛选符合条件的真实SNV
                    if baseq_val >= limitbasequa and loc_val > limitloc:
                        # 校验mismatch索引，避免越界
                        if i < len(mismatch):
                            truesnv.append([mismatch[i], seq[1], seq[8]])

                    i += 1  # 简化自增，符合Python编码习惯

                # 计算平均碱基质量，避免除零异常
                if len(baseqlst) > 0:
                    avg_baseq = sum(baseqlst) / float(len(baseqlst))
                    if avg_baseq >= limitbasequa:
                        for snv in truesnv:
                            # 校验snv[0]拆分后的索引，避免异常
                            try:
                                snv_parts = snv[0].split(':')
                                if len(snv_parts) < 2:
                                    continue
                                # 使用f-string优化字符串拼接，替代原有的+拼接
                                snv_key = f"{seq[0]}\t{snv_parts[0]}\t{snv_parts[1]}"
                                seq1_int_current = int(seq[1])
                                # 提取条件判断，提升可读性
                                cond1 = (len(bin(seq1_int_current)) > 5 and bin(seq1_int_current)[-5] == '1' and snv[2][-1] == '1')
                                cond2 = ((len(bin(seq1_int_current)) < 5 or bin(seq1_int_current)[-5] == '0') and snv[2][-1] == '2')
                                cond3 = (len(bin(seq1_int_current)) > 5 and bin(seq1_int_current)[-5] == '1' and snv[2][-1] == '2')
                                cond4 = ((len(bin(seq1_int_current)) < 5 or bin(seq1_int_current)[-5] == '0') and snv[2][-1] == '1')

                                # 更新allsnv字典，Python3异常捕获语法修复
                                try:
                                    allsnv[snv_key][0] += 1
                                    if cond1 or cond2:
                                        allsnv[snv_key][1] += 1
                                    elif cond3 or cond4:
                                        allsnv[snv_key][2] += 1
                                except Exception as e:  # Python2→Python3 核心修复：except Exception as e
                                    if cond1 or cond2:
                                        allsnv[snv_key] = [1, 1, 0]
                                    elif cond3 or cond4:
                                        allsnv[snv_key] = [1, 0, 1]
                            except (IndexError, ValueError):
                                continue

        # 构造SNV BED数据
        snv_bed = []
        for snv in allsnv:
            seq_snv = snv.split('\t')
            # 列数校验，避免拆分后索引不足
            if len(seq_snv) < 3:
                continue
            snv_count = allsnv[snv]
            if snv_count[0] >= limitad:
                try:
                    loc_int = int(seq_snv[2])
                except ValueError:
                    continue
                # 按计数判断链方向，保留原逻辑
                if snv_count[1] > snv_count[2]:
                    snv_bed.append([seq_snv[0], loc_int, seq_snv[1], '+', snv_count[0]])
                elif snv_count[2] > snv_count[1]:
                    snv_bed.append([seq_snv[0], loc_int, seq_snv[1], '-', snv_count[0]])
                else:
                    snv_bed.append([seq_snv[0], loc_int, seq_snv[1], '.', snv_count[0]])

        # 排序并写入输出文件
        snv_bed.sort()
        for one in snv_bed:
            # 使用f-string优化输出拼接，简洁高效
            output_line = f"{one[0]}\t{one[1]-1}\t{one[1]}\t{one[2]}\t{one[4]}\t{one[3]}\n"
            fo.write(output_line)