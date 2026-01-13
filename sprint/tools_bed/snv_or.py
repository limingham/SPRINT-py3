def snv_or(bed_in_dir1, bed_in_dir2, bed_out_dir):
    whole = {}
    # 多文件嵌套 with 上下文管理器，自动关闭所有文件，无需手动调用 close()
    with open(bed_in_dir1, 'r') as f1, \
         open(bed_in_dir2, 'r') as f2, \
         open(bed_out_dir, 'w') as fo:

        # 读取第一个文件，构建字典并存储数据
        for line in f1:
            # 去除行尾换行符，跳过空行
            line_rstripped = line.rstrip()
            if not line_rstripped:
                continue
            seq = line_rstripped.split('\t')
            # 列数校验，确保满足 seq[0]-seq[5] 索引要求
            if len(seq) < 6:
                continue
            # 构造唯一键，保留原逻辑
            key = ':'.join(seq[0:4]) + ':' + seq[5]
            whole[key] = int(seq[4])

        # 读取第二个文件，存在相同键则累加值，不存在则新增
        for line in f2:
            line_rstripped = line.rstrip()
            if not line_rstripped:
                continue
            seq = line_rstripped.split('\t')
            if len(seq) < 6:
                continue
            key = ':'.join(seq[0:4]) + ':' + seq[5]
            if key in whole:
                whole[key] += int(seq[4])
            else:
                whole[key] = int(seq[4])

        # 构造排序用列表
        lst = []
        for one in whole:
            seq = one.split(':')
            # 列数校验，避免拆分后索引不足
            if len(seq) < 5:
                continue
            lst.append([seq[0], int(seq[1]), int(seq[2]), seq[3], str(whole[one]), seq[4]])

        # 对列表进行排序
        lst.sort()

        # 写入输出文件，使用列表推导式简化代码
        for one in lst:
            out = [str(i) for i in one]
            fo.write('\t'.join(out) + '\n')
