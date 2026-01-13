def transcript_locator(fbed_in_dir=0, ftransloc_in_dir=0, fbed_out_dir=0):
    if fbed_in_dir == 0 or ftransloc_in_dir == 0:
        # 关键修改1：Python3要求print为函数，必须加括号
        print('fbed_in_dir\tftransloc_in_dir\tfbed_out_dir')
        return 0

    # 关键修改2：使用with语句（上下文管理器）自动管理文件关闭，替代手动open/close
    # 同时兼容Python3的文件操作模式，明确指定文本模式't'（默认也是t，显式指定更规范）
    with open(ftransloc_in_dir, 'r', encoding='utf-8') as ftransloc, \
         open(fbed_in_dir, 'r', encoding='utf-8') as fi, \
         open(fbed_out_dir, 'w', encoding='utf-8') as fo:

        transcript = {}
        whole = ftransloc.read().split('>')[1:]
        for one in whole:
            line = one.split('\n')
            transcript[line[0]] = line[1]
            # 关键修改3：Python3 print函数需加括号
            # print(line[0])

        for line in fi:
            seq = line.rstrip().split('\t')
            # print(seq[0])  # 同样需加括号
            if seq[0] in transcript:
                transcript_id = seq[0]
                chrr = transcript_id.split('_|_')[1]
                trans_loc = int(seq[2])
                loc = []
                bande = transcript[seq[0]].split(';')[:-1]
                for one in bande:
                    be = one.split(',')
                    begin = int(be[0])
                    end = int(be[1])
                    tmp = [begin] * (end - begin + 1)
                    i = 0
                    while i < len(tmp):
                        tmp[i] = tmp[i] + i
                        i += 1
                    loc += tmp
                ref_loc = loc[trans_loc - 1]
                # 无需修改字符串拼接逻辑（Python3字符串拼接与Python2一致）
                fo.write(chrr + '\t' + str(ref_loc - 1) + '\t' + str(ref_loc) + '\t' + '\t'.join(seq[3:]) + '\n')
                # 't'+chrom[chrr][ref_loc-1]+':'+chrom_t[transcript_id][trans_loc-1]+'\n')

    # 原代码缺失函数返回值，补充None使函数行为更规范（可选，不影响功能）
    return None
    