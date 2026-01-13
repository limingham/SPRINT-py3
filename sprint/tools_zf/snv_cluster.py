def snv_cluster(bed_in_dir=0, bed_out_dir=0, cluster_distance=-1, cluster_size=-1):
    # 使用Python3推荐的with上下文管理器，自动关闭文件无需手动调用close()
    with open(bed_in_dir, 'r') as fi, open(bed_out_dir, 'w') as fo:
        tmp = 'chr0:0:AA'
        limitdistance = int(cluster_distance)
        limitnum = int(cluster_size)
        lst = []
        for line in fi:
            # 去除行尾空白，跳过空行，提升健壮性
            line_stripped = line.strip()
            if not line_stripped:
                continue
            seq = line.split('\t')
            # 列数校验，避免后续索引越界
            if len(seq) < 4:
                continue
            tmpseq = tmp.split(':')
            # 索引安全校验，避免tmp拆分后列数不足
            if len(tmpseq) < 3:
                lst = [line]
                tmp = f"{seq[0]}:{seq[2]}:{seq[3]}"
                continue
            # 核心聚类判断逻辑，保留原业务逻辑
            try:
                chr_match = (seq[0] == tmpseq[0])
                distance_match = (int(seq[2]) - int(tmpseq[1])) <= limitdistance
                snv_match = (seq[3] == tmpseq[2])
                if chr_match and distance_match and snv_match:
                    lst.append(line)
                else:
                    # 满足聚类大小阈值则写入数据
                    if len(lst) >= limitnum:
                        # 避免lst为空或拆分后索引不足
                        if len(lst) == 0:
                            lst = [line]
                            tmp = f"{seq[0]}:{seq[2]}:{seq[3]}"
                            continue
                        first_line_seq = lst[0].split('\t')
                        last_line_seq = lst[-1].split('\t')
                        if len(first_line_seq) < 2 or len(last_line_seq) < 3:
                            lst = [line]
                            tmp = f"{seq[0]}:{seq[2]}:{seq[3]}"
                            continue
                        # 计算起始位置、终止位置和密度
                        begin = float(first_line_seq[1])
                        end = float(last_line_seq[2])
                        # 避免除零异常（起始位置等于终止位置）
                        density = len(lst) / (end - begin) if (end - begin) != 0 else 0.0
                        # 遍历写入聚类后的每条记录
                        for one in lst:
                            # 使用rstrip('\n')替代one[0:-1]，兼容最后一行无换行符场景
                            one_stripped = one.rstrip('\n')
                            output_line = f"{one_stripped}\t{len(lst)}\t{density}\n"
                            fo.write(output_line)
                    # 重置聚类列表
                    lst = []
                    lst.append(line)
                # 更新临时标记（使用f-string优化拼接）
                tmp = f"{seq[0]}:{seq[2]}:{seq[3]}"
            except (ValueError, IndexError):
                # 捕获整数转换/索引异常，重置聚类列表
                lst = [line]
                tmp = f"{seq[0]}:{seq[2]}:{seq[3]}"
                continue

        # 处理文件末尾剩余的聚类数据
        if len(lst) >= limitnum:
            if len(lst) == 0:
                return
            first_line_seq = lst[0].split('\t')
            last_line_seq = lst[-1].split('\t')
            if len(first_line_seq) >= 2 and len(last_line_seq) >= 3:
                begin = float(first_line_seq[1])
                end = float(last_line_seq[2])
                density = len(lst) / (end - begin) if (end - begin) != 0 else 0.0
                for one in lst:
                    one_stripped = one.rstrip('\n')
                    output_line = f"{one_stripped}\t{len(lst)}\t{density}\n"
                    fo.write(output_line)