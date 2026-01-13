def maskTwithC(fa_in_dir=0, fa_out_dir=0):
    if fa_in_dir == 0 or fa_out_dir == 0:
        # Python3 强制要求print为函数，添加圆括号（核心语法修复）
        print('fa_in_dir\tfa_our_dir')
        return 0
    # 使用Python3推荐的with上下文管理器，自动关闭文件，无需手动调用close()
    with open(fa_in_dir, 'r') as fi, open(fa_out_dir, 'w') as fo:
        for line in fi:
            # 保留原逻辑：识别FASTA文件的标识符行（以>开头）
            if line[0] == '>':
                fo.write(line)
            else:
                # 保留原逻辑：将所有T/t替换为C/c
                fo.write(line.replace('T', 'C').replace('t', 'c'))