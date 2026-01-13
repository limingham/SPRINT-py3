def tzz2gzz(trans_loc_file_in, transcript_zz_in, genome_zz_out):
    # 使用Python3推荐的with上下文管理器管理转录本位置文件，自动关闭
    with open(trans_loc_file_in, 'r') as fa:
        TRANS = {}
        l1 = fa.readline()
        l2 = fa.readline()
        while l1 != '':
            # 提取转录本名称（去除首字符和行尾空白）
            trans = l1[1:].rstrip()
            # 拆分序列并去除最后一个空元素（兼容原逻辑）
            seq = l2.split(';')[:-1]
            if trans:  # 避免空转录本名称
                TRANS[trans] = seq
            # 读取下一组两行数据
            l1 = fa.readline()
            l2 = fa.readline()

    def loc_t2g(tCHR, tLOC):
        # 提取基因组染色体名称
        chr_parts = tCHR.split('_|_')
        if len(chr_parts) < 2:
            return 0  # 格式异常返回0，避免后续报错
        gCHR = chr_parts[1]
        # 校验转录本是否存在于TRANS字典中
        if tCHR not in TRANS:
            return 0
        seq = TRANS[tCHR]
        try:
            tLOC = int(tLOC)
        except ValueError:
            return 0  # 转换失败返回0
        flag = 1
        tmp = 0
        i = 0
        while i < len(seq) and flag == 1:
            # 拆分起始和终止位置，增加异常捕获
            pos_parts = seq[i].split(',')
            if len(pos_parts) < 2:
                i += 1
                continue
            try:
                end = int(pos_parts[1])
                start = int(pos_parts[0])
            except ValueError:
                i += 1
                continue
            tmp += end - start + 1
            if tmp >= tLOC:
                flag = 0
            i += 1
        # 计算基因组位置，避免i=0导致j=-1
        if i == 0:
            return 0
        j = i - 1
        dis2end = tmp - tLOC
        # 再次校验j对应的位置格式
        j_pos_parts = seq[j].split(',')
        if len(j_pos_parts) < 2:
            return 0
        try:
            gLOC = int(j_pos_parts[1]) - dis2end
        except ValueError:
            return 0
        return gLOC

    def range_t2g(tCHR, tstart, tend):
        # 提取基因组染色体名称
        chr_parts = tCHR.split('_|_')
        if len(chr_parts) < 2:
            return []
        gCHR = chr_parts[1]
        # 转换起始和终止位置为整数，增加异常捕获
        try:
            tstart = int(tstart)
            tend = int(tend)
        except ValueError:
            return []
        # 生成转录本位置列表
        if tend < tstart:
            return []
        tloc = [tstart + i for i in range(tend - tstart + 1)]
        gloc = []
        for one in tloc:
            g_loc_val = loc_t2g(tCHR, one)
            if g_loc_val != 0:  # 过滤无效位置
                gloc.append(g_loc_val)
        if not gloc:
            return []
        # 合并连续的基因组位置区间
        i = 1
        record = [gloc[0]]
        out = []
        while i < len(gloc):
            if abs(gloc[i] - gloc[i-1]) > 1:
                out.append(f"{record[0]}:{record[-1]}")
                record = [gloc[i]]
            else:
                record.append(gloc[i])
            i += 1
        # 添加最后一个区间
        out.append(f"{record[0]}:{record[-1]}")
        return out

    # 使用with上下文管理器管理输入输出ZZ文件，自动关闭
    with open(transcript_zz_in, 'r') as fi, open(genome_zz_out, 'w') as fo:
        for line in fi:
            # 跳过空行，提升健壮性
            line_stripped = line.strip()
            if not line_stripped:
                continue
            seq = line.split('\t')
            # 列数校验，避免后续索引越界
            if len(seq) < 5:
                continue
            tCHR = seq[0]
            # 提取基因组染色体名称
            chr_parts = tCHR.split('_|_')
            if len(chr_parts) < 2:
                continue
            gCHR = chr_parts[1]

            # 处理SNV位点（转录本坐标转基因组坐标）
            if seq[4] != '*':
                snv = seq[4].split(';')
                gsnv = []
                for one in snv:
                    one_parts = one.split(':')
                    if len(one_parts) < 2:
                        continue
                    # 转换转录本SNV位置为基因组位置
                    g_loc_val = loc_t2g(tCHR, one_parts[1])
                    if g_loc_val != 0:
                        gsnv.append(f"{one_parts[0]}:{g_loc_val}")
                # 更新seq[4]为基因组SNV信息
                seq[4] = ';'.join(gsnv)

            # 处理区间信息（转录本坐标转基因组坐标）
            trange = seq[3].split(';')
            grange = []
            for one in trange:
                one_parts = one.split(':')
                if len(one_parts) < 2:
                    continue
                # 转换转录本区间为基因组区间
                new_one = range_t2g(tCHR, one_parts[0], one_parts[1])
                grange += new_one
            # 更新seq[3]为基因组区间信息
            seq[3] = ';'.join(grange)
            # 更新染色体名称为基因组染色体
            seq[0] = gCHR

            # 写入转换后的行
            fo.write('\t'.join(seq))