def res_or(bed_in_dir1, bed_in_dir2, bed_out_dir):
    whole = {}
    # 读取第一个文件，使用with上下文管理器自动关闭文件
    with open(bed_in_dir1, 'r') as f1:
        for line in f1:
            # 去除行尾换行符，跳过空行
            line_rstripped = line.rstrip()
            if not line_rstripped:
                continue
            seq = line_rstripped.split('\t')
            # 列数校验，确保满足索引要求（至少6列：seq[0]-seq[5]）
            if len(seq) < 6:
                continue
            # 构造字典键，保留原逻辑
            key = ':'.join(seq[0:4]) + ':' + seq[5]
            whole[key] = int(seq[4])
    
    # 读取第二个文件，同样使用with上下文管理器
    with open(bed_in_dir2, 'r') as f2:
        for line in f2:
            line_rstripped = line.rstrip()
            if not line_rstripped:
                continue
            seq = line_rstripped.split('\t')
            if len(seq) < 6:
                continue
            key = ':'.join(seq[0:4]) + ':' + seq[5]
            # 保留原逻辑：键已存在则跳过，不存在则添加
            if key not in whole:
                whole[key] = int(seq[4])
    
    lst = []
    # 遍历字典，构造排序用列表
    for one in whole:
        seq = one.split(':')
        # 列数校验，避免拆分后索引不足
        if len(seq) < 5:
            continue
        lst.append([seq[0], int(seq[1]), int(seq[2]), seq[3], str(whole[one]), seq[4]])
    
    # 对列表进行排序
    lst.sort()
    
    # 写入输出文件，自动关闭文件句柄
    with open(bed_out_dir, 'w') as fo:
        for one in lst:
            # 简化类型转换与拼接：直接生成字符串列表
            out = [str(i) for i in one]
            fo.write('\t'.join(out) + '\n')