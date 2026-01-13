import random

def random_select(file_in_dir=0, file_out_dir=0, ratio=1):
    # 使用Python3推荐的with上下文管理器，自动关闭文件无需手动调用close()
    with open(file_in_dir, 'r') as fi, open(file_out_dir, 'w') as fo:
        # 遍历输入文件每一行，按指定概率随机筛选写入
        for line in fi:
            if random.random() <= ratio:
                fo.write(line)
