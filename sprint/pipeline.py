# 移除无效注释 #import sprint
def pipeline():
    import re
    import sys
    from sprint.sprint_prepare import prepare as PP
    from sprint.sprint_main import main as MA

    def help_doc():
        print('')
        print("##############################################################################################")
        print('')
        print("      SPRINT: SNP-free RNA editing Identification Toolkit")
        print("")
        print("      https://github.com/TianLab-Bioinfo/SPRINT/blob/master/sprint_py3/")
        print("")
        print("      Please contact limh25@m.fudan.edu.cn when questions arise.")
        print("")
        print("##############################################################################################")
        print("")
        print("")
        print('   Usage: ')
        print("")
        print('      "sprint_prepare"     #build mapping index')
        print("")
        print('      "sprint_main"        #identify RESs')
        print("")
        print('      "sprint_from_bam"    #identify RESs from aligned reads')
        print("")
        print('      "sprint_main_parallel"        #identify RESs in parallel')
        print("")
        print('      "sprint_from_bam_parallel"    #identify RESs from aligned reads in parallel')
        print("")
        print("")
        sys.exit(0)

    if len(sys.argv) < 2:
        help_doc()

#     if sys.argv[1] == 'sprint_prepare' or sys.argv[1] == 'prepare' or sys.argv[1] == 'pp':
#         sys.argv = sys.argv[1:]
#         PP()

#     elif sys.argv[1] == 'sprint_main' or sys.argv[1] == 'main' or sys.argv[1] == 'ma':
#         sys.argv = sys.argv[1:]
#         MA()
    else:
        help_doc()
