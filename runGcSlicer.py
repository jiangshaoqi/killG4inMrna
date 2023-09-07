import os

# rna序列文件名，无格式
seq_file_name = 'seqDemo2'
result_dir_name = seq_file_name + 'Rst'
# 第一次运行
dir_path = './' + result_dir_name + '/'
if not os.path.exists(dir_path):
    os.mkdir(dir_path)
seq_file_name_txt = seq_file_name + '.txt'
g4_command = 'python G4Hunter.py -i ' + seq_file_name_txt + ' -o ' + result_dir_name + ' -w 25 -s 1.4'
os.system(g4_command)
# 根据merge文件修改，并检测，直到对应的merge文件不包含不到序列
merge_file_name = dir_path + 'Results_' + seq_file_name + '/' + seq_file_name + '-Merged.txt'
isEmptyMerge = False
run_time = 1
while True:
    # 检测merged 文件行数，如merged无值，直接退出
    with open(merge_file_name, 'r') as merge_file:
        seq_name = merge_file.readline()
        key_line = merge_file.readline()
        value_line = merge_file.readline()
        if len(value_line) == 0:
            isEmptyMerge = True
    if isEmptyMerge:
        break
    # 如merged有值，开始运行
    gcSlicer_command = 'python gcSlicer.py ' + merge_file_name + ' ' + str(run_time) + ' ' + seq_file_name_txt
    os.system(gcSlicer_command)
    # 更新seq_file_name_txt，循环
    seq_file_name = seq_file_name + '_' + str(run_time)
    result_dir_name = seq_file_name + 'Rst'
    dir_path = './' + result_dir_name + '/'
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    seq_file_name_txt = seq_file_name + '.txt'
    g4_command = 'python G4Hunter.py -i ' + seq_file_name_txt + ' -o ' + result_dir_name + ' -w 25 -s 1.4'
    os.system(g4_command)
    merge_file_name = dir_path + 'Results_' + seq_file_name + '/' + seq_file_name + '-Merged.txt'

