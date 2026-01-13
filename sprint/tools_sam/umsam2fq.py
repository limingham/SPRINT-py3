def umsam2fq(sam_in_dir, fq_out_dir):
    # 使用Python3推荐的with上下文管理器，同时管理输入输出文件，自动关闭无需手动close()
    with open(sam_in_dir, 'r') as fi, open(fq_out_dir, 'w') as fo:
        for line in fi:
            # 去除行尾换行符，跳过空行，提升健壮性
            line_rstripped = line.rstrip()
            if not line_rstripped:
                continue
            seq = line_rstripped.split('\t')
            # 增加列数校验，避免seq[1]索引越界
            if len(seq) < 11:
                continue
            try:
                # 尝试转换seq[1]为整数，避免非数字字符导致报错
                flag_int = int(seq[1])
            except ValueError:
                continue
            # 保留原逻辑的多条件判断：非头部行 + 标志位位运算判断
            if line[0] != '@' and len(bin(flag_int)) >= 5 and bin(flag_int)[-3] == '1':
                # 根据标志位进一步调整读取名称
                if len(bin(flag_int)) >= 9 and bin(flag_int)[-7] == '1':
                    seq[0] = seq[0][0:-2] + '_1'
                elif len(bin(flag_int)) >= 10 and bin(flag_int)[-8] == '1':
                    seq[0] = seq[0][0:-2] + '_2'
                elif line[-1] == '1':
                    seq[0] = seq[0][0:-2] + '_1'
                elif line[-1] == '2':
                    seq[0] = seq[0][0:-2] + '_2'
                # 使用Python3 f-string优化字符串拼接，简洁高效
                fastq_content = f"@{seq[0]}\n{seq[9]}\n+\n{seq[10]}\n"
                fo.write(fastq_content)