def antisense_reverse(read):
    read = read.upper()
    read_change_base = ""
    for one in read:
        if one == 'A':
            read_change_base += 'T'
        elif one == 'C':
            read_change_base += 'G'
        elif one == 'G':
            read_change_base += 'C'
        elif one == 'T':
            read_change_base += 'A'
        else:
            read_change_base += 'N'
    read_reverse = read_change_base[::-1]
    return read_reverse


def maskfq(fq_in_dir, mask_from, mask_to):
    mask_from = mask_from.upper()
    mask_to = mask_to.upper()
    # 构造输出文件名，保留原逻辑
    out_filename = fq_in_dir[0:-3] + '_' + mask_from + '_to_' + mask_to + '.fq'
    # 使用Python3推荐的with上下文管理器，自动关闭输入输出文件，无需手动close()
    with open(fq_in_dir, 'r') as fi, open(out_filename, 'w') as fo:
        # 读取FastQ第一组4行数据，去除换行符
        line1 = fi.readline().replace('\n', '')
        line2 = fi.readline().upper().replace('\n', '')
        line3 = fi.readline().replace('\n', '')
        line4 = fi.readline().replace('\n', '')
        
        # 循环处理，直到读取不到第一行（文件结束）
        while line1 != '':
            # 保留原逻辑：判断行尾是否为'1'，若是则进行反义链反转和质量行反转
            if line1[-1] == '1':
                line2 = antisense_reverse(line2)
                line4 = line4[::-1]
            
            record = "1"
            line2_new = ""
            # 遍历序列字符，构建记录字符串
            for one in line2.replace('\n', ''):
                if one == mask_from:
                    record = record + '1'
                else:
                    record = record + '0'
            
            # 拼接输出ID行，保留原逻辑
            output_line1 = line1 + '_|_' + mask_from + '_to_' + mask_to + '_|_' + str(int(record, 2)) + '_|_read2' + '\n'
            fo.write(output_line1)
            # 写入替换后的序列行、分隔行和质量行
            fo.write(line2.replace(mask_from, mask_to) + '\n')
            fo.write('+\n')
            fo.write(line4 + '\n')
            
            # 读取下一组4行数据
            line1 = fi.readline().replace('\n', '')
            line2 = fi.readline().upper().replace('\n', '')
            line3 = fi.readline().replace('\n', '')
            line4 = fi.readline().replace('\n', '')