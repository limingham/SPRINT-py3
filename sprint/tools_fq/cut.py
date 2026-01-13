def cut(fq_in_dir=0, fq_out_dir=0, cutnum=0, name='read1'):
    # 使用Python3推荐的with上下文管理器，自动关闭文件，无需手动调用close()
    with open(fq_in_dir, 'r') as fi, open(fq_out_dir, 'w') as fo:
        # 转换cutnum为整数，保留原逻辑
        cutnum = int(cutnum)
        # 读取FastQ文件的4行一组数据
        line1 = fi.readline()
        line2 = fi.readline()
        line3 = fi.readline()
        line4 = fi.readline()
        idd = 1
        # 循环处理，直到读取不到第一行（文件结束）
        while line1 != '':
            CELL_TAG = ''
            # 提取XC:Z:对应的CELL_TAG，保留原逻辑
            if "XC:Z:" in line1:
                seq = line1.split('_')
                for one in seq:
                    if one[:5] == 'XC:Z:':
                        CELL_TAG = one
            # 写入新的序列ID行
            if CELL_TAG != '':
                fo.write(f"@id_{idd}_{CELL_TAG}_{name}\n")
            else:
                fo.write(f"@id_{idd}_{name}\n")
            # 写入剪切后的序列和质量行，保留原逻辑
            fo.write(line2[cutnum:])
            fo.write(line3)
            fo.write(line4[cutnum:])
            # 读取下一组4行数据
            line1 = fi.readline()
            line2 = fi.readline()
            line3 = fi.readline()
            line4 = fi.readline()
            idd += 1