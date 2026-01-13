
def change_sam_read_name(sam_in_dir=0, sam_out_dir=0, name='read1'):
    # 使用Python3推荐的with上下文管理器，同时管理输入输出文件，自动关闭无需手动close()
    with open(sam_in_dir, 'r') as fi, open(sam_out_dir, 'w') as fo:
        idd = 1
        # 遍历SAM文件每一行
        for line in fi:
            # 非头部行（不以@开头），修改读取名称
            if line[0] != '@':
                seq = line.split('\t')
                # 替换原始读取名称为新名称，保留后续内容
                new_line = f"id_{idd}_{name}{line.replace(seq[0], '')}"
                fo.write(new_line)
                idd += 1  # 简化变量自增，更符合Python编码习惯
            else:
                # 头部行（以@开头），原样写入
                fo.write(line)