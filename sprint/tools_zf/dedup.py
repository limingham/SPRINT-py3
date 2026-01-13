def dedup(in_dir, out_dir):
    old = set()
    # 使用Python3推荐的with上下文管理器，同时管理输入输出文件，自动关闭无需手动close()
    with open(in_dir, 'r') as fi, open(out_dir, 'w') as fo:
        for line in fi:
            # 去除行尾空白，跳过空行，提升健壮性
            line_stripped = line.strip()
            if not line_stripped:
                continue
            seq = line.split('\t')
            # 增加列数校验，避免seq[0]、seq[3]、seq[7]索引越界
            if len(seq) < 8:
                continue
            # 构造去重唯一键，使用f-string优化拼接
            key = f"{seq[0]}:{seq[3]}:{seq[7]}"
            # 保留原逻辑：键不存在则写入并添加到集合
            if key not in old:
                fo.write(line)
                old.add(key)