
def get_baseq_cutoff(fq_in_dir=0, cutoff_out_dir=0):
    # 使用Python3推荐的with上下文管理器，同时管理输入输出文件，自动关闭无需手动close()
    with open(fq_in_dir, 'r') as fi, open(cutoff_out_dir, 'w') as fo:
        # 读取FastQ文件的第一组4行数据（ID行、序列行、分隔行、质量行）
        line1 = fi.readline()
        line2 = fi.readline()
        line3 = fi.readline()
        line4 = fi.readline()
        did = 0
        # 循环处理，直到读取不到第一行或触发did=1的终止条件
        while line1 != '':
            # 若did=1，直接终止循环（仅处理第一组有效数据）
            if did == 1:
                break
            # 去除质量行末尾的换行符，保留原逻辑
            qua = line4[0:-1]
            # 遍历质量行的每个字符，判断ASCII码值
            for i in qua:
                if ord(i) > 76:
                    fo.write('89')  # 64+25，对应质量值阈值判断
                    did = 1
                    break
                if ord(i) < 60:
                    fo.write('58')  # 33+25，对应质量值阈值判断
                    did = 1
                    break
            # 读取下一组4行数据（即使触发终止，仍执行一次读取，不影响功能）
            line1 = fi.readline()
            line2 = fi.readline()
            line3 = fi.readline()
            line4 = fi.readline()