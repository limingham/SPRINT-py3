def bed_or(bed_in_dir1, bed_in_dir2, bed_out_dir):
    whole = []
    # 自动关闭文件，无需手动调用close()
    with open(bed_in_dir1) as f1:
        for line in f1:
            seq = line.replace('\n', '').split('\t')
            whole.append([seq[0], int(seq[2]), seq[3], seq[4], seq[5]])
    
    with open(bed_in_dir2) as f2:
        for line in f2:
            seq = line.replace('\n', '').split('\t')
            whole.append([seq[0], int(seq[2]), seq[3], seq[4], seq[5]])
    
    whole.sort()
    old = set()
    # 写入文件时同样使用上下文管理器
    with open(bed_out_dir, 'w') as fo:
        for one in whole:
            key = one[0] + ':' + str(one[1])
            if key not in old:
                output_line = one[0] + '\t' + str(one[1]-1) + '\t' + str(one[1]) + '\t' + one[2] + '\t' + one[3] + '\t' + one[4] + '\n'
                fo.write(output_line)
                old.add(key)