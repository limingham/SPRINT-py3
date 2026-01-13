def maskAwithG(fa_in_dir=0, fa_out_dir=0):
    if fa_in_dir == 0 or fa_out_dir == 0:
        # Python3 中print是函数，必须添加圆括号
        print('fa_in_dir\tfa_our_dir')
        return 0
    # 使用with上下文管理器，自动关闭文件，无需手动调用close()
    with open(fa_in_dir, 'r') as fi, open(fa_out_dir, 'w') as fo:
        for line in fi:
            if line[0] == '>':
                fo.write(line)
            else:
                # 保留原逻辑：将所有A/a替换为G/g
                fo.write(line.replace('A', 'G').replace('a', 'g'))