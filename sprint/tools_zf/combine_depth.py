def combine_depth(bed_in_dir=0, bed_out_dir=0):
    # 使用Python3推荐的with上下文管理器，同时管理输入输出文件，自动关闭无需手动close()
    with open(bed_in_dir, 'r') as fi, open(bed_out_dir, 'w') as fo:
        for line in fi:
            # 去除行尾换行符，跳过空行，提升健壮性
            line_stripped = line.rstrip('\n')
            if not line_stripped:
                continue
            seq = line_stripped.split('\t')
            # 增加列数校验，避免seq[4]、seq[-1]索引越界
            if len(seq) < 6:
                continue
            try:
                # 拆分AD和DP值并转换为整数，捕获拆分/转换异常
                ad1 = int(seq[4].split(':')[0])
                ad2 = int(seq[-1].split(':')[0])
                dp1 = int(seq[4].split(':')[1])
                dp2 = int(seq[-1].split(':')[1])
            except (IndexError, ValueError):
                # 格式异常时跳过当前行，不中断程序
                continue
            # 计算累加后的AD和DP
            AD = ad1 + ad2
            DP = dp1 + dp2
            # 使用Python3 f-string优化字符串拼接，简洁高效
            output_line = f"{seq[0]}\t{seq[1]}\t{seq[2]}\t{seq[3]}\t{AD}:{DP}\t{seq[5]}\n"
            fo.write(output_line)