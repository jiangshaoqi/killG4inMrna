from Bio import SeqIO
import sys


def find_longest_consecutive_string(s):
    longest_string = ""  # 最长连续字符串
    current_string = ""

    for char in s:
        if char == 'C' or char == 'G':
            if char == current_string[-1:]:
                current_string += char
            else:
                current_string = char
            if len(current_string) > len(longest_string):
                longest_string = current_string
        else:
            current_string = ''
    return longest_string


protein_rna_mapping = {"M": ["AUG"],  # Met
                       "T": ["ACA", "ACC", "ACG", "ACU"],  # Thr: AC4
                       "N": ["AAC", "AAU"],  # Asn: AA2
                       "K": ["AAA", "AAG"],  # Lys: AA2
                       # Ser: AG2, UC4
                       "S": ["AGC", "AGU", "UCA", "UCC", "UCG", "UCU"],
                       # Arg: AG2, CG4
                       "R": ["AGA", "AGG", "CGA", "CGC", "CGG", "CGU"],
                       # Val: GU4
                       "V": ["GUA", "GUC", "GUG", "GUU"],
                       # Ala: GC4
                       "A": ["GCA", "GCC", "GCG", "GCU"],
                       "D": ["GAC", "GAU"],  # Asp: GA2
                       "E": ["GAA", "GAG"],  # Glu: GA2
                       # Gly: GG4
                       "G": ["GGA", "GGC", "GGG", "GGU"],
                       "F": ["UUC", "UUU"],  # Phe: UU2
                       # Leu: UU2, CU4
                       "L": ["UUA", "UUG", "CUA", "CUC", "CUG", "CUU"],
                       "Y": ["UAC", "UAU"],  # Tyr: UA2
                       "C": ["UGC", "UGU"],  # Cys: UG2
                       "W": ["UGG"],  # Trp: UGG
                       # Pro: CC4
                       "P": ["CCA", "CCC", "CCG", "CCU"],
                       "H": ["CAC", "CAU"],  # His: CA2
                       "Q": ["CAA", "CAG"],  # Gln: CA2
                       "I": ["AUA", "AUC", "AUU"],  # Ile: AU3
                       # Sec: UA2, UGA
                       "*": ["UAA", "UAG", "UGA"]
                       }


def get_rna_candidate(protein, replace_idx, origin_rna):
    # replace_idx: 0, 1 or 2 最佳替换点在三个对应蛋白质的三个rna的位置
    # origin_r: 原本该位置的rna
    c1 = None
    c2 = None
    c3 = None
    c4 = None
    c5 = None
    c6 = None
    if origin_rna[replace_idx] == 'C':
        complement_r = 'G'
    else:
        complement_r = 'C'
    candidate_list = protein_rna_mapping[protein]
    if replace_idx == 0:
        for candidate in candidate_list:
            if candidate[0] == 'A' or candidate[0] == 'U':
                c1 = candidate
            if candidate[1] == 'A' or candidate[1] == 'U':
                c2 = candidate
            if candidate[2] == 'A' or candidate[2] == 'U':
                c3 = candidate
            if candidate[0] == complement_r:
                c4 = candidate
            if candidate[1] == complement_r:
                c5 = candidate
            if candidate[2] == complement_r:
                c6 = candidate
    if replace_idx == 1:
        for candidate in candidate_list:
            if candidate[1] == 'A' or candidate[1] == 'U':
                c1 = candidate
            if candidate[0] == 'A' or candidate[0] == 'U':
                c2 = candidate
            if candidate[2] == 'A' or candidate[2] == 'U':
                c3 = candidate
            if candidate[1] == complement_r:
                c4 = candidate
            if candidate[0] == complement_r:
                c5 = candidate
            if candidate[2] == complement_r:
                c6 = candidate
    if replace_idx == 2:
        for candidate in candidate_list:
            if candidate[2] == 'A' or candidate[2] == 'U':
                c1 = candidate
            if candidate[1] == 'A' or candidate[1] == 'U':
                c2 = candidate
            if candidate[0] == 'A' or candidate[0] == 'U':
                c3 = candidate
            if candidate[2] == complement_r:
                c4 = candidate
            if candidate[1] == complement_r:
                c5 = candidate
            if candidate[0] == complement_r:
                c6 = candidate
    if c1:
        return c1
    if c2:
        return c2
    if c3:
        return c3
    if c4:
        return c4
    if c5:
        return c5
    if c6:
        return c6
    return origin_rna


rna_protein_mapping = {}
for key in protein_rna_mapping:
    rna_list = protein_rna_mapping[key]
    for rna in rna_list:
        rna_protein_mapping[rna] = key

# 首先，读seq
# part 1: read g4Hunter result file
g4_merged_file_name = sys.argv[1]
run_time = int(sys.argv[2])
start_pos_list = []
new_rna_tri_list = []
with open(g4_merged_file_name, 'r') as g4_merged_file:
    # 第一行，名字
    seq_name = g4_merged_file.readline()
    seq_name = seq_name[1:]
    # 第二行，标题
    info_keys = g4_merged_file.readline().split()
    # 第三行，数值：读数值直到最后一行
    while True:
        info_vals = g4_merged_file.readline().split()
        if len(info_vals) == 0:
            break
        info_map = dict(zip(info_keys, info_vals))

        # 处理序列
        target_seq = info_map['Sequence']
        longest_seq = find_longest_consecutive_string(target_seq)
        idx_longest = target_seq.index(longest_seq)
        # 选中merged文件中，局部seq的index
        idx_offset = idx_longest + len(longest_seq) // 2
        # 计算该位置在总序列的位置
        idx_seq = idx_offset + int(info_map['Start'])
        # 计算要替换的序列位置
        idx_start_replace = idx_seq // 3 * 3
        # print(idx_start_replace)
        # 进行替换，首先找到局部序列的idx
        rna_tri_idx = idx_start_replace - int(info_map['Start'])
        rna_tri = target_seq[rna_tri_idx: rna_tri_idx + 3]
        target_protein = rna_protein_mapping[rna_tri]
        # print(idx_offset - rna_tri_idx)
        new_rna_tri = get_rna_candidate(target_protein, idx_offset - rna_tri_idx, rna_tri)
        # 记录 idx_start_replace和 new_rna_tri
        start_pos_list.append(idx_start_replace)
        new_rna_tri_list.append(new_rna_tri)

# part 2: 生成新的rna seq文件
replace_map = dict(zip(start_pos_list, new_rna_tri_list))
old_rna_file_name = sys.argv[3]
new_seq = ''
for seq_record in SeqIO.parse(old_rna_file_name, "fasta"):
    old_id = seq_record.id
    old_seq = seq_record.seq
    # print(old_seq)
    new_seq = str(seq_record.seq)
    for start_pos, rna_tri in replace_map.items():
        new_seq = new_seq[: start_pos] + rna_tri + new_seq[start_pos + 3:]
    # print(new_seq)

new_rna_file_name = old_rna_file_name.split('.')[0] + '_' + str(run_time) + '.txt'
with open(new_rna_file_name, 'w') as new_file:
    new_file.write('>' + seq_name)
    new_file.write(new_seq)

