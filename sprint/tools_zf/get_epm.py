def get_epm(bed_in_dir, epm_out_dir, flag, zz_in_dir):
    # 使用Python3推荐的with上下文管理器，同时管理3个文件，自动关闭无需手动close()
    with open(bed_in_dir, 'r') as fi, \
         open(epm_out_dir, 'w') as fo, \
         open(zz_in_dir, 'r') as fread:

        gene = {}
        for line in fi:
            # 去除换行符并拆分行数据，兼容最后一行无换行符场景
            line_stripped = line.replace('\n', '')
            seq = line_stripped.split('\t')
            i = 0
            while i < len(seq):
                if flag in seq[i]:
                    try:
                        # 拆分AD值并累加，保留原逻辑
                        ad_value = int(seq[4].split(':')[0])
                        gene[seq[i]][0] += ad_value
                        gene[seq[i]][1] += 1
                    except Exception as e:  # Python2→Python3 核心修复：异常捕获使用as关键字
                        ad_value = int(seq[4].split(':')[0])
                        gene[seq[i]] = [ad_value, 1]
                i += 1  # 简化自增，符合Python编码习惯

        # 统计zz文件的总行数（read_count）
        read_count = 0
        for _ in fread:  # 使用下划线替代未使用的变量j，更符合Python编码规范
            read_count += 1

        # 计算百万读段数
        read_million = float(read_count) / 1000000

        # 计算EPM并写入输出文件
        for one in gene:
            # 计算EPM值，保留原逻辑
            epm = float(gene[one][0]) / read_million
            # 使用Python3 f-string优化字符串拼接，简洁高效
            output_line = f"{one}\t{epm}\t{gene[one][0]}\t{gene[one][1]}\n"
            fo.write(output_line)
