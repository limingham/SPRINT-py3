def get_depth(zz_in_dir=0, bed_in_dir=0, bed_out_dir=0):
    # 使用Python3推荐的with上下文管理器，同时管理3个文件，自动关闭无需手动close()
    with open(zz_in_dir, 'r') as fread, \
         open(bed_in_dir, 'r') as fsnv, \
         open(bed_out_dir, 'w') as fo:

        class Read:
            def __init__(self, read):
                # 初始化Read类属性，保留原逻辑
                self.snv = read.split('\t')[4].split(';')
                self.inter = read.split('\t')[3].split(';')
                self.direct = read.split('\t')[1]

            def locisin(self, loc):
                # 判断位置是否在区间内
                isin = 0
                for inter in self.inter:
                    inter = inter.split(':')
                    # 增加索引校验，避免拆分后索引不足
                    if len(inter) < 2:
                        continue
                    try:
                        loc_int = int(loc)
                        inter0_int = int(inter[0])
                        inter1_int = int(inter[1])
                    except ValueError:
                        continue
                    if loc_int <= inter1_int and loc_int >= inter0_int:
                        isin = 1
                        break
                return isin  # 简化返回逻辑，直接返回0/1

            def snvisin(self, snv):
                # 判断snv是否存在
                return 1 if snv in self.snv else 0

            def getmin(self):
                # 获取最小位置
                if not self.inter:
                    return 0
                inter0 = self.inter[0].split(':')
                if len(inter0) < 2:
                    return 0
                try:
                    return int(inter0[0])
                except ValueError:
                    return 0

            def getmax(self):
                # 获取最大位置
                if not self.inter:
                    return 0
                inter_last = self.inter[-1].split(':')
                if len(inter_last) < 2:
                    return 0
                try:
                    return int(inter_last[1])
                except ValueError:
                    return 0

        reads = {}
        for line in fread:
            seq = line.split('\t')
            # 增加列数校验，避免seq[0]索引越界
            if len(seq) < 1:
                continue
            read_key = seq[0]
            line_stripped = line.rstrip('\n')  # 替代line[0:-1]，兼容最后一行无换行符场景
            try:
                reads[read_key].append(Read(line_stripped))
            except Exception as e:  # Python2→Python3 核心修复：异常捕获使用as关键字
                print(f"{read_key} begin")  # Python3 f-string优化字符串格式化
                reads[read_key] = [Read(line_stripped)]

        top = 0
        chrr = ''
        for line in fsnv:
            line_stripped = line.rstrip('\n')  # 替代line[0:-1]，兼容性更强
            seq = line_stripped.split('\t')
            # 增加列数校验，避免后续索引越界
            if len(seq) < 4:
                continue
            deep = 0
            altdeep = 0
            snv = f"{seq[3]}:{seq[2]}"  # f-string优化拼接

            # 染色体切换时重置top
            if seq[0] != chrr:
                top = 0
                chrr = seq[0]

            # 深度统计逻辑，保留原业务逻辑
            if seq[0] in reads and top < len(reads[seq[0]]):
                # 跳过最大位置小于当前loc的read
                while (seq[0] == chrr and top < len(reads[seq[0]]) and
                       reads[seq[0]][top].getmax() < int(seq[2])):
                    top += 1  # 简化自增，符合Python习惯

                point = top
                # 统计符合条件的read深度
                while (seq[0] == chrr and point < len(reads[seq[0]]) and
                       reads[seq[0]][point].getmin() <= int(seq[2])):
                    if reads[seq[0]][point].locisin(seq[2]) == 1:
                        deep += 1
                    if reads[seq[0]][point].snvisin(snv) == 1:
                        altdeep += 1
                    point += 1  # 简化自增，符合Python习惯

            # 写入输出行，f-string优化拼接
            output_line = f"{line_stripped}\t{altdeep}:{deep}\n"
            fo.write(output_line)