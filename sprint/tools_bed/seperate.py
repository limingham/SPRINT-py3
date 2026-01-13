def seperate(bed_in_dir, flag_out_dir, rp_out_dir, nonrp_out_dir, flag):
    # 使用 Python3 推荐的 with 上下文管理器，嵌套管理多个文件，自动关闭无需手动调用 close()
    with open(bed_in_dir, 'r') as fi, \
         open(flag_out_dir, 'w') as fo_flag, \
         open(rp_out_dir, 'w') as fo_rp, \
         open(nonrp_out_dir, 'w') as fo_nonrp:
        
        for line in fi:
            # 去除首尾空白（含换行符），跳过空行，提升健壮性
            line_clean = line.strip()
            if not line_clean:
                continue
            
            # 保留原逻辑的分支判断，修复无效的 next 语句
            if 'Simple_repeat' in line or 'Low_complexity' in line:
                continue  # 替换原无效的 next 语句，实现跳过当前行的功能
            elif flag in line:
                fo_flag.write(line)
            elif 'Repeat_region' in line:
                fo_rp.write(line)
            else:
                fo_nonrp.write(line)