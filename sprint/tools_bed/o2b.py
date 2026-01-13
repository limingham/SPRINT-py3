def o2b(bed_in, bed_out):
    # 使用 Python3 推荐的 with 上下文管理器，自动关闭文件，避免资源泄露
    with open(bed_in, 'r') as fi:
        # 写入文件同样使用 with 上下文管理器
        with open(bed_out, 'w') as fo:
            for line in fi:
                # 去除行尾换行符（与原代码 rstrip() 效果一致，保留行首/行中空白）
                line_rstripped = line.rstrip()
                if not line_rstripped:  # 跳过空行，提升代码健壮性
                    continue
                seq = line_rstripped.split('\t')
                # 列数校验，避免因 seq 不足6列导致 join 报错
                if len(seq) < 6:
                    continue
                # 保留原逻辑：拼接前6列并写入，使用 \t 分隔
                output_line = '\t'.join(seq[0:6]) + '\n'
                fo.write(output_line)