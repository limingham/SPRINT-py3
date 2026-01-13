def poly_check(seq, poly_limit):
    # 保留原逻辑：检查序列中是否存在连续poly_limit个相同碱基
    if 'A' * poly_limit not in seq and 'T' * poly_limit not in seq and 'G' * poly_limit not in seq and 'C' * poly_limit not in seq:
        return True
    else:
        return False


def var_check(change_from, change_to, seq, var_limit):
    ALL = ['A', 'T', 'C', 'G']
    tmp = []
    for one in ALL:
        if one != change_from.upper() and one != change_to.upper():
            tmp.append(one)
    flag = 1
    for one in tmp:
        # Python3中除法默认返回浮点数，可移除多余的float()转换，保留原逻辑兼容
        if seq.count(one) < var_limit / (float(len(tmp)) + 2):
            flag = 0
    if flag == 1:
        return True
    else:
        return False


def reverse_base(base):
    base = base.upper()
    # 碱基互补配对逻辑，保留原功能
    if base == 'A':
        return 'T'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'
    elif base == 'T':
        return 'A'
    else:
        return 'N'


def recover_sam(sam_in_dir, sam_out_dir, var_limit=20, poly_limit=10, rm_multi=0):
    # 使用Python3推荐的with上下文管理器，自动关闭文件，无需手动close()
    with open(sam_in_dir, 'r') as fi, open(sam_out_dir, 'w') as fo:
        for line in fi:
            # 先拆分序列，避免重复split操作
            seq = line.split('\t')
            # 头部行（以@开头），原样写入
            if line[0] == '@':
                fo.write(line)
            # 未比对行，触发循环终止
            elif seq[1] == '4' and seq[2] == '*':
                break
            # 有效比对行，进行数据处理
            elif seq[1] != '4' and len(seq) >= 9:
                seq = line.split('\t')
                # 序列转为大写
                seq[9] = seq[9].upper()
                # 转换flag为整数进行位运算
                seq[1] = int(seq[1])
                # 位运算判断，修改flag值
                if len(bin(seq[1])) >= 7:
                    if bin(seq[1])[-3] != '1':
                        if bin(seq[1])[-5] == '1':
                            seq[1] = '16'
                        else:
                            seq[1] = '0'
                seq[1] = str(seq[1])
                
                # 提取记录、替换源碱基和目标碱基
                record = bin(int(seq[0].split('_|_')[2]))[3:]
                change_from = seq[0].split('_|_')[1].split('_')[0]
                change_to = seq[0].split('_|_')[1].split('_')[2]
                
                # 反向互补判断，调整碱基和记录
                if len(bin(int(seq[1]))) > 5 and bin(int(seq[1]))[-5] == '1':
                    change_from = reverse_base(change_from)
                    change_to = reverse_base(change_to)
                    record = record[::-1]
                else:
                    record = record
                
                # 序列恢复逻辑
                changed_read = seq[9]
                i = 0
                recovered_read = ''
                while i < len(seq[9]):
                    if record[i] == '1' and seq[9][i] == change_to:
                        recovered_read += change_from
                    elif record[i] == '1' and seq[9][i] != change_to:
                        recovered_read += 'N'
                    else:
                        recovered_read += seq[9][i]
                    i += 1  # 简化自增，符合Python习惯
                seq[9] = recovered_read
                
                # 过滤条件判断，满足条件则写入
                if len(record) == len(seq[9]) and len(record) - changed_read.count(change_to) > var_limit and poly_check(seq[9], poly_limit):
                    if (rm_multi == 1 and "XA:Z:" not in line) or rm_multi == 0:
                        fo.write(seq[0])
                        j = 1
                        while j < len(seq):
                            fo.write('\t' + seq[j])
                            j += 1  # 简化自增，符合Python习惯