def hyper_or(bed_in_dir1, bed_in_dir2, bed_out_dir):
    whole = []
    # 读取第一个BED文件，使用with上下文管理器自动关闭文件，无需手动调用close()
    with open(bed_in_dir1, 'r') as f1:
        for line in f1:
            # 去除首尾空白字符（含换行符、空格），跳过无效空行
            line_clean = line.strip()
            if not line_clean:
                continue
            seq = line_clean.split('\t')
            # 列数校验，确保满足索引要求（至少需要6列：seq[0]、seq[2]、seq[3]、seq[5]、seq[-1]）
            if len(seq) < 6:
                continue
            whole.append([seq[0], int(seq[2]), seq[3], seq[-1], seq[5]])
    
    # 读取第二个BED文件，同样使用with上下文管理器
    with open(bed_in_dir2, 'r') as f2:
        for line in f2:
            line_clean = line.strip()
            if not line_clean:
                continue
            seq = line_clean.split('\t')
            if len(seq) < 6:
                continue
            whole.append([seq[0], int(seq[2]), seq[3], seq[-1], seq[5]])
    
    # 对合并后的所有数据进行排序
    whole.sort()
    old = set()  # 用于去重的集合，存储"染色体:位置"关键字
    
    # 写入输出文件，使用with上下文管理器自动关闭
    with open(bed_out_dir, 'w') as fo:
        for one in whole:
            # 使用Python3 f-string构造去重关键字，简洁高效
            key = f"{one[0]}:{one[1]}"
            if key not in old:
                # 同样使用f-string格式化输出行，替代原有字符串拼接
                output_line = f"{one[0]}\t{one[1]-1}\t{one[1]}\t{one[2]}\t{one[3]}\t{one[4]}\n"
                fo.write(output_line)
                old.add(key)  # 将已写入的关键字加入去重集合