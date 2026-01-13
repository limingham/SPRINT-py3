def bed_sort(bed_in_dir, bed_out_dir):
    whole = []
    # 使用 Python3 推荐的 with 上下文管理器，自动关闭文件，避免资源泄露
    with open(bed_in_dir, 'r') as fi:
        for line in fi:
            # 去除首尾空白字符（含换行符），避免空行或多余空格干扰
            line_stripped = line.strip()
            if not line_stripped:  # 跳过空行，提升代码健壮性
                continue
            seq = line_stripped.split('\t')
            # 列数校验，避免因 BED 文件格式异常导致索引报错
            if len(seq) < 3:
                continue
            # 保留原有逻辑：存储染色体、起始位置、终止位置和原始行
            whole.append([seq[0], int(seq[1]), int(seq[2]), line])
    
    # 按染色体（字符串）、起始位置（整数）、终止位置（整数）自然排序
    whole.sort()
    
    # 写入文件，同样使用 with 上下文管理器
    with open(bed_out_dir, 'w') as fo:
        for one in whole:
            fo.write(one[3])