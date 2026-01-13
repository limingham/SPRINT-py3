def sort_zz(zz_in, zz_out):
    out = []
    # 使用Python3推荐的with上下文管理器，自动关闭文件无需手动调用close()
    with open(zz_in, 'r') as fi, open(zz_out, 'w') as fo:
        for line in fi:
            # 去除行尾空白，跳过空行，提升健壮性
            line_stripped = line.strip()
            if not line_stripped:
                continue
            seq = line.split('\t')
            # 增加列数和索引校验，避免后续拆分/转换异常
            try:
                # 校验seq[3]存在且能按:拆分
                inter_part = seq[3]
                inter_split = inter_part.split(':')
                if len(inter_split) < 2:
                    continue  # 格式异常则跳过当前行
                # 提取起始位置和终止位置并转换为整数
                start_pos = int(inter_split[0])
                end_pos = int(inter_split[-1])
                # 添加到排序列表中
                out.append([seq[0], start_pos, end_pos, line])
            except (IndexError, ValueError):
                # 捕获索引越界或整数转换异常，跳过无效行
                continue
        # 按默认规则排序（先染色体、再起始位置、再终止位置）
        out.sort()
        # 遍历写入排序后的内容
        for one in out:
            fo.write(one[3])