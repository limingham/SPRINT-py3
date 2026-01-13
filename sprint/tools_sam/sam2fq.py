def sam2fq(sam_in_dir, fq_out_dir):
    # 使用Python3推荐的with上下文管理器，同时管理输入输出文件，自动关闭无需手动close()
    with open(sam_in_dir, 'r') as fi, open(fq_out_dir, 'w') as fo:
        for line in fi:
            # 去除行尾换行符，跳过空行
            line_rstripped = line.rstrip()
            if not line_rstripped:
                continue
            seq = line_rstripped.split('\t')
            # 非头部行（不以@开头），转换为FastQ格式写入
            if line[0] != '@':
                # 增加列数校验，避免seq[9]、seq[10]索引越界
                if len(seq) >= 11:
                    # 使用f-string优化字符串拼接，简洁高效
                    fastq_content = f"@{seq[0]}\n{seq[9]}\n+\n{seq[10]}\n"
                    fo.write(fastq_content)