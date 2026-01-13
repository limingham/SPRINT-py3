def get_res(regular_alu, regular_nalurp, regular_nrp, hyper, output):
    # 使用Python3推荐的with上下文管理器，同时管理7个文件，自动关闭无需手动close()
    with open(regular_alu, 'r') as falu, \
         open(regular_nalurp, 'r') as fnalurp, \
         open(regular_nrp, 'r') as fnrp, \
         open(hyper, 'r') as fhyper, \
         open(output + '_A_to_I_regular.res', 'w') as fo_AI, \
         open(output + '_A_to_I_hyper.res', 'w') as fo_hAI, \
         open(output + '_C_to_U.res', 'w') as fo_CU:

        # 处理regular_alu文件，筛选A-to-I相关记录
        for line in falu:
            line_stripped = line.rstrip()
            if not line_stripped:  # 跳过空行，提升健壮性
                continue
            seq = line_stripped.split('\t')
            if len(seq) >= 4 and (seq[3] == 'AG' or seq[3] == 'TC'):  # 列数校验，避免索引越界
                fo_AI.write(line)

        # 处理regular_nalurp文件，分别筛选A-to-I和C-to-U相关记录
        for line in fnalurp:
            line_stripped = line.rstrip()
            if not line_stripped:
                continue
            seq = line_stripped.split('\t')
            if len(seq) >= 4:
                if seq[3] == 'AG' or seq[3] == 'TC':
                    fo_AI.write(line)
                if seq[3] == 'CT' or seq[3] == 'GA':
                    fo_CU.write(line)

        # 处理regular_nrp文件，分别筛选A-to-I和C-to-U相关记录
        for line in fnrp:
            line_stripped = line.rstrip()
            if not line_stripped:
                continue
            seq = line_stripped.split('\t')
            if len(seq) >= 4:
                if seq[3] == 'AG' or seq[3] == 'TC':
                    fo_AI.write(line)
                if seq[3] == 'CT' or seq[3] == 'GA':
                    fo_CU.write(line)

        # 处理hyper文件，筛选A-to-I相关记录
        for line in fhyper:
            line_stripped = line.rstrip()
            if not line_stripped:
                continue
            seq = line_stripped.split('\t')
            if len(seq) >= 4 and (seq[3] == 'AG' or seq[3] == 'TC'):
                fo_hAI.write(line)