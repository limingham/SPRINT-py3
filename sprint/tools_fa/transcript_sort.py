def transcript_sort(fbed_in_dir=0, fbed_out_dir=0):
    # 使用Python3推荐的with上下文管理器，同时管理输入输出文件，自动关闭无需手动close()
    with open(fbed_in_dir, 'r') as fi, open(fbed_out_dir, 'w') as fo:
        whole = {}
        for line in fi:
            # 去除行尾换行符，跳过空行，提升健壮性
            line_rstripped = line.rstrip()
            if not line_rstripped:
                continue
            seq = line_rstripped.split('\t')
            # 增加列数校验，避免索引越界
            if len(seq) < 6:
                continue
            # 构造字典唯一键
            key = f"{seq[0]}:{seq[2]}:{seq[3]}:{seq[5]}"
            try:
                # 键已存在则累加数值
                whole[key] += int(seq[4])
            except Exception as e:  # Python2→Python3 核心语法修复：异常捕获使用as关键字
                # 键不存在则初始化数值
                whole[key] = int(seq[4])

        tmp = []
        for one in whole:
            seq = one.split(':')
            # 增加拆分后列数校验，避免索引异常
            if len(seq) < 4:
                continue
            dep = str(whole[one])
            # 构造排序用列表，保留原始字段映射逻辑
            tmp.append([seq[0], int(seq[1]), seq[2], dep, seq[3]])
        
        # 对列表进行自然排序（与原逻辑一致）
        tmp.sort()
        
        for one in tmp:
            # 使用Python3 f-string优化字符串拼接，简洁高效
            output_line = f"{one[0]}\t{one[1]-1}\t{one[1]}\t{one[2]}\t{one[3]}\t{one[4]}\n"
            fo.write(output_line)