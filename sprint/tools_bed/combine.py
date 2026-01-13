def combine(bed_in_dir1, bed_in_dir2, bed_in_dir3, bed_out_dir):
    whole = []
    # 读取第一个BED文件，使用with上下文管理器自动关闭文件
    with open(bed_in_dir1, 'r') as f1:
        for line in f1:
            # 去除首尾空白字符（含换行符），避免空行或多余空格干扰
            line_clean = line.strip()
            if not line_clean:  # 跳过空行，提升健壮性
                continue
            seq = line_clean.split('\t')
            # 列数校验，确保满足索引要求（至少需要6列：seq[0]、seq[2]、seq[3]、seq[5]、seq[-1]）
            if len(seq) < 6:
                continue
            whole.append([seq[0], int(seq[2]), seq[3], seq[5], seq[-1]])
    
    # 读取第二个BED文件
    with open(bed_in_dir2, 'r') as f2:
        for line in f2:
            line_clean = line.strip()
            if not line_clean:
                continue
            seq = line_clean.split('\t')
            if len(seq) < 6:
                continue
            whole.append([seq[0], int(seq[2]), seq[3], seq[5], seq[-1]])
    
    # 读取第三个BED文件
    with open(bed_in_dir3, 'r') as f3:
        for line in f3:
            line_clean = line.strip()
            if not line_clean:
                continue
            seq = line_clean.split('\t')
            if len(seq) < 6:
                continue
            whole.append([seq[0], int(seq[2]), seq[3], seq[5], seq[-1]])
    
    # 对合并后的数据按默认规则排序（染色体→位置→其他字段）
    whole.sort()
    
    # 写入输出文件，自动关闭文件句柄
    with open(bed_out_dir, 'w') as fo:
        for one in whole:
            # 使用Python3 f-string格式化字符串，简洁高效且可读性强
            output_line = f"{one[0]}\t{one[1]-1}\t{one[1]}\t{one[2]}\t{one[-1]}\t{one[3]}\n"
            fo.write(output_line)