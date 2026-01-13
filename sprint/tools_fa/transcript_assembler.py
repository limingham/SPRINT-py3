def transcript_assembler(fref_in_dir=0, fgtf_in_dir=0, ftrans_out_dir=0):
    if fref_in_dir == 0 or fgtf_in_dir == 0:
        # Python3 核心语法修复：print语句改为print函数，添加圆括号
        print('fref_in_dir\tfgtf_in_dir\tftrans_out_dir')
        return 0

    # 使用with上下文管理器自动关闭fref文件，避免句柄泄露
    with open(fref_in_dir, 'r') as fref:
        chrom = {}
        chrr = ''
        line = fref.read()
        line = line.split('>')
        for seq in line:
            if ' ' in seq:
                chrr = seq[0:seq.find(' ')]
            else:
                chrr = seq[0:seq.find('\n')]
            # 避免空序列导致的索引异常
            if seq.find('\n') != -1:
                chrom[chrr] = seq[seq.find('\n'):].replace('\n', '').upper()

    # 保留原始的反义链反转注释函数，无需修改（Python3 兼容）
    # def antisense_reverse(read):
    #        read=read.upper()
    #        read_change_base=""   
    #        for one in read:
    #                if one == 'A':
    #                        read_change_base += 'T'
    #                elif one == 'C':
    #                        read_change_base += 'G'
    #                elif one == 'G':
    #                        read_change_base += 'C'
    #                elif one == 'T':
    #                        read_change_base += 'A'
    #                else:  
    #                        read_change_base += 'N'
    #        read_reverse=read_change_base[::-1]
    #        return read_reverse

    # 使用with上下文管理器自动关闭fgtf文件
    with open(fgtf_in_dir, 'r') as fgtf:
        transcript = {}
        trc = []

        for line in fgtf:
            if '#' != line[0]:
                seq = line.rstrip().split('\t')
                # 增加seq索引校验，避免索引越界
                if len(seq) >= 9 and seq[0] in chrom and seq[2] == 'exon' and 'transcript_id ' in seq[8]:
                    try:
                        end_pos = int(seq[4])
                    except ValueError:
                        continue  # 跳过无效的位置数值
                    if end_pos <= len(chrom[seq[0]]):
                        chrr = seq[0]
                        begin = int(seq[3])
                        end = end_pos
                        # 提取transcript_id，保留原始逻辑
                        transcript_id = seq[8].split('transcript_id ')[1].split(';')[0].replace('"', '')
                        strand = seq[6]
                        if transcript_id not in transcript:
                            trc.append(transcript_id)
                            transcript[transcript_id] = [chrr, strand, [begin, end]]
                        else:
                            transcript[transcript_id].append([begin, end])

    # 对每个转录本的外显子位置进行排序
    for transcript_id in transcript:
        tmp = transcript[transcript_id][2:]
        tmp.sort()
        transcript[transcript_id][2:] = tmp
        # if transcript[transcript_id][1]=='-':
        #    transcript[transcript_id][2:]=transcript[transcript_id][2:][::-1]

    # 使用with上下文管理器同时管理两个输出文件，自动关闭
    with open(ftrans_out_dir, 'w') as ftrans, \
         open(ftrans_out_dir + '.loc', 'w') as ftransloc:

        for transcript_id in trc:
            # 构造转录本名称，保留原始拼接逻辑
            trans_name = '>' + transcript_id + '_|_' + transcript[transcript_id][0] + '_|_' + transcript[transcript_id][1] + '\n'
            chrr = transcript[transcript_id][0]
            strand = transcript[transcript_id][1]
            loc = transcript[transcript_id][2:]
            trans_seq = ''
            ftransloc.write(trans_name)
            for one in loc:
                ftransloc.write(str(one[0]) + ',' + str(one[1]) + ';')
                # 拼接外显子对应的序列，注意Python切片左闭右开特性（与原逻辑一致）
                trans_seq += chrom[chrr][one[0]-1:one[1]]
            ftransloc.write('\n')
            ftrans.write(trans_name)
            i = 0
            while i < len(trans_seq):
                ftrans.write(trans_seq[i])
                i += 1
                if i % 50 == 0:
                    ftrans.write('\n')
            if i % 50 != 0:
                ftrans.write('\n')