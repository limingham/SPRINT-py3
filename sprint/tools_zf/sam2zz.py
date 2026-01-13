import re

def sam2zz(sam_in_dir=0, fa_in_dir=0, zz_out_dir=0):
    # -----------------------------------------------------------
    # reading refgenome
    # 使用with上下文管理器管理参考基因组文件，自动关闭
    with open(fa_in_dir, 'r') as fref:
        chrom = {}
        chrr = ''
        line = fref.read()
        line = line.split('>')
        for seq in line:
            if ' ' in seq:
                chrr = seq[0:seq.find(' ')]
            else:
                # 兼容无空格场景，按换行符分割染色体名称
                chr_end_idx = seq.find('\n')
                if chr_end_idx != -1:
                    chrr = seq[0:chr_end_idx]
                else:
                    chrr = ''  # 空序列跳过后续处理
            # 提取染色体序列并去除换行符
            seq_start_idx = seq.find('\n')
            if seq_start_idx != -1 and chrr:
                chrom_seq = seq[seq_start_idx:].replace('\n', '')
                chrom[chrr] = chrom_seq

    # 0 base
    # ------------------------------------------------------------
    # donee:dolist;lst:donumlist
    def doCG(a):
        donee = []
        lst = re.findall(r'(\d+|\+|-|\*|/)', a)
        for i in a:
            if i in ['I', 'D', 'M', 'S', 'P', 'N']:
                donee.append(i)
        return donee, lst

    # donefunction
    def doneCG(CG, chrr, pos, seq, qseq):  # pos is 1 base
        donee, lst = doCG(CG)
        errorsite = ''
        intersite = ''
        quasite = ''
        locsite = ''
        pieceloc = ''
        refseq = ''
        seqseq = ''
        refpos = int(pos) - 1
        seqpos = 0
        step = 0
        while step < len(donee):
            if donee[step] == 'I':
                seqpos = seqpos + int(lst[step])
            elif donee[step] == 'D':
                refpos = refpos + int(lst[step])
            elif donee[step] == 'N':
                refpos = refpos + int(lst[step])
            elif donee[step] == 'S':
                seqpos = seqpos + int(lst[step])
            elif donee[step] == 'M':
                # 避免索引越界，先判断lst[step]是否为有效整数
                try:
                    step_len = int(lst[step])
                except (ValueError, IndexError):
                    step += 1
                    continue
                # 提取参考序列和读取序列
                if chrr in chrom:
                    ref_seq_slice = chrom[chrr][refpos:refpos + step_len]
                else:
                    ref_seq_slice = ''
                seq_seq_slice = seq[seqpos:seqpos + step_len]
                refseq += ref_seq_slice
                seqseq += seq_seq_slice

                j = refpos
                jj = seqpos
                while j < refpos + step_len:
                    try:
                        # 避免染色体不存在或索引越界
                        if chrr not in chrom:
                            break
                        ref_base = chrom[chrr][j].upper()
                        read_base = seq[jj].upper()
                        if ref_base != read_base and ref_base != 'N' and read_base != 'N':
                            errorsite += f"{ref_base}{read_base}:{j + 1};"
                            quasite += f",{ord(qseq[jj])}"
                            locsite += f",{jj + 1}"
                            min_loc = min(jj + 1 - seqpos, seqpos + step_len - jj)
                            pieceloc += f",{min_loc}"
                    except Exception as e:  # Python2→Python3 异常捕获修复：使用as关键字
                        pass
                        # 如需启用打印，使用Python3 f-string格式
                        # print(f"error with {chrr} {pos} {e}")
                    j += 1
                    jj += 1
                # 拼接区间信息
                intersite += f"{refpos + 1}:{refpos + step_len};"
                refpos += step_len
                seqpos += step_len
            step += 1  # 简化自增，符合Python编码习惯
        refseq = refseq.upper()
        seqseq = seqseq.upper()
        return refseq, seqseq, errorsite, intersite, quasite, locsite, pieceloc

    # ------------------------------------------------------------
    ###################additional
    '''
    whole={}
    fi=open(sam_in_dir)
    for line in fi:
        seq=line.split('\t')
        if line[0]!='@'and len(seq)>5:
            if seq[0][0]!='@' and seq[2]!='*' and seq[5]!='*':
                name=seq[0].split('_|_')[0]
                
                #tmp = [ seq[4].count(':') ,seq[0] ]
                if name in whole:
                #	if tmp[0] < whole[name][0]:
                #		whole[name]=tmp
                    
                    whole[name] +=1
                else:
                    whole[name]=1
                    #whole[name]=tmp
    fi.close()
    '''
    #######################
    # 使用with上下文管理器管理SAM输入文件和ZZ输出文件，自动关闭
    with open(sam_in_dir, 'r') as fi, open(zz_out_dir, 'w') as fo:
        for line in fi:
            seq = line.split('\t')
            if line[0] != '@' and len(seq) > 5:
                # 拆分读取名称，兼容无"_|_"的场景
                name_parts = seq[0].split('_|_')
                name = name_parts[0] if name_parts else seq[0]
                # 过滤有效比对记录
                if seq[0][0] != '@' and seq[2] != '*' and seq[5] != '*':  # and whole[name]<2:
                    # 调用doneCG处理CIGAR字符串
                    refseq, seqseq, errorsite, intersite, quasite, locsite, pieceloc = doneCG(
                        seq[5], seq[2], seq[3], seq[9], seq[10]
                    )
                    # 去除前缀逗号
                    quasite = quasite[1:] if len(quasite) > 0 else quasite
                    locsite = locsite[1:] if len(locsite) > 0 else locsite
                    pieceloc = pieceloc[1:] if len(pieceloc) > 0 else pieceloc

                    # 处理空值，赋值默认值
                    if len(intersite) > 0 and len(intersite[0:-1]) == 0:
                        intersite = '*;'
                    if len(errorsite) > 0 and len(errorsite[0:-1]) == 0:
                        errorsite = '*;'
                    if len(quasite) == 0:
                        quasite = '*'
                    if len(locsite) == 0:
                        locsite = '*'
                    if len(pieceloc) == 0:
                        pieceloc = '*'

                    # 处理SAM标志位，判断并拼接读取名称后缀
                    try:
                        seq1_int = int(seq[1])
                        bin_seq1 = bin(seq1_int)
                        # 判断标志位并写入输出
                        if len(bin_seq1) >= 9 and bin_seq1[-7] == '1':
                            output_line = (
                                f"{seq[2]}\t{seq[1]}\t{seq[4]}\t{intersite[0:-1]}\t{errorsite[0:-1]}\t"
                                f"{quasite}\t{locsite}\t{seq[9]}\t{seq[0]}_1\t{pieceloc}\n"
                            )
                            fo.write(output_line)
                        elif len(bin_seq1) >= 10 and bin_seq1[-8] == '1':
                            output_line = (
                                f"{seq[2]}\t{seq[1]}\t{seq[4]}\t{intersite[0:-1]}\t{errorsite[0:-1]}\t"
                                f"{quasite}\t{locsite}\t{seq[9]}\t{seq[0]}_2\t{pieceloc}\n"
                            )
                            fo.write(output_line)
                        else:
                            output_line = (
                                f"{seq[2]}\t{seq[1]}\t{seq[4]}\t{intersite[0:-1]}\t{errorsite[0:-1]}\t"
                                f"{quasite}\t{locsite}\t{seq[9]}\t{seq[0]}\t{pieceloc}\n"
                            )
                            fo.write(output_line)
                    except ValueError:
                        # 若seq[1]无法转换为整数，跳过当前行
                        continue