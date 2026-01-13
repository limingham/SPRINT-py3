def remove(bed_in_dir, bed_out_dir, flag):
    # 使用 Python3 推荐的 with 上下文管理器，自动关闭文件，无需手动调用 close()
    with open(bed_in_dir, 'r') as fi:
        with open(bed_out_dir, 'w') as fo:
            for line in fi:
                # 保留原逻辑：判断 flag 不在当前行中，才将该行写入输出文件
                if flag not in line:
                    fo.write(line)
