def get_snv_with_ad(snv_in_dir=0, snv_out_dir=0, flag=0):
    # 使用Python3推荐的with上下文管理器，同时管理输入输出文件，自动关闭无需手动close()
    with open(snv_in_dir, 'r') as fi, open(snv_out_dir, 'w') as fo:
        for line in fi:
            # 去除行尾空白，跳过空行，提升健壮性
            line_stripped = line.strip()
            if not line_stripped:
                continue
            seq = line.split('\t')
            try:
                # 严格转换为整数后进行阈值判断，保留原逻辑
                if int(seq[4]) >= int(flag):
                    fo.write(line)
            except Exception as e:  # Python2→Python3 核心语法修复：异常捕获使用as关键字
                # Python2→Python3 核心修复：print改为函数，使用f-string优化字符串格式化
                error_info = f"{seq[0]}\t{seq[2]}\twithout\tAD flag\n"
                print(error_info)