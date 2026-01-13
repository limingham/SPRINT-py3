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
                # 判断位置是否在区间内，保留原逻辑并增加健壮性
                isin = 0
                for inter in self.inter:
                    inter_parts = inter.split(':')
                    # 增加索引校验，避免拆分后索引不足
                    if len(inter_parts) < 2:
                        continue
                    try:
                        loc_int = int(loc)
                        inter0_int = int(inter_parts[0])
                        inter1_int = int(inter_parts[1])
                    except (ValueError, TypeError):
                        continue
                    if loc_int <= inter1_int and loc_int >= inter0_int:
                        isin = 1
                        break
                return 1 if isin == 1 else 0  # 简化返回逻辑

            def snvisin(self, snv):
                # 判断snv是否存在，简化返回逻辑
                return 1 if snv in self.snv else 0

            def getmin(self):
                # 获取最小位置，增加异常防护
                if not self.inter:
                    return 0
                inter0_parts = self.inter[0].split(':')
                if len(inter0_parts) < 2:
                    return 0
                try:
                    return int(inter0_parts[0])
                except (ValueError, TypeError):
                    return 0

            def getmax(self):
                # 获取最大位置，增加异常防护
                if not self.inter:
                    return 0
                inter_last_parts = self.inter[-1].split(':')
                if len(inter_last_parts) < 2:
                    return 0
                try:
                    return int(inter_last_parts[1])
                except (ValueError, TypeError):
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
                reads[read_key].append(line_stripped)
            except Exception as e:  # Python2→Python3 核心修复：异常捕获使用as关键字
                # 注释保留，如需启用可修改为Python3 print格式
                # print(f"{read_key} begin")
                reads[read_key] = [line_stripped]

        top = 0
        chrr = ''
        for line in fsnv:
            line_stripped = line.rstrip()
            seq = line_stripped.split('\t')
            deep = 0
            altdeep = 0
            try:
                # 增加列数校验，避免seq[2]、seq[3]索引越界
                if len(seq) < 4:
                    raise IndexError("Line has insufficient columns")
                snv = f"{seq[3]}:{seq[2]}"  # Python3 f-string优化字符串拼接
                if seq[0] != chrr:
                    top = 0
                    chrr = seq[0]
                if seq[0] not in reads:
                    reads[seq[0]] = []
                if top < len(reads[seq[0]]):
                    # 跳过最大位置小于当前loc的read，简化自增
                    while (seq[0] == chrr and top < len(reads[seq[0]]) and
                           Read(reads[seq[0]][top]).getmax() < int(seq[2])):
                        top += 1

                    point = top
                    # 统计符合条件的read深度，简化自增
                    while (seq[0] == chrr and point < len(reads[seq[0]]) and
                           Read(reads[seq[0]][point]).getmin() <= int(seq[2])):
                        read_instance = Read(reads[seq[0]][point])  # 避免重复创建实例
                        if read_instance.locisin(seq[2]) == 1:
                            deep += 1
                        if read_instance.snvisin(snv) == 1:
                            altdeep += 1
                        point += 1
                # 写入输出行，使用f-string优化，兼容无换行符场景
                output_line = f"{line_stripped}\t{altdeep}:{deep}\n"
                fo.write(output_line)
            except Exception as e:
                # Python2→Python3 核心修复：print改为函数，使用f-string格式化
                print(f"Error processing line: {line}")