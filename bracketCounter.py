def count_bracket(fold_seq):
    seq_list = []
    cur_count = 0
    for idx, s in enumerate(fold_seq):
        if idx == 0:
            cur_count = 1
            continue
        if fold_seq[idx-1] == s:
            cur_count += 1
        else:
            if s == '(':
                cur_count = cur_count * -1
            seq_list.append(cur_count)
            cur_count = 1
    seq_list.append(cur_count * -1)
    # 生成list，由连续同向括号组成
    print(seq_list)
    list_len = len(seq_list)
    for i in range(1, list_len):
        # 检查新值应大于0还是小于0
        p = 0
        if abs(seq_list[i-1]) > abs(seq_list[i]):
            p = seq_list[i-1]
        if abs(seq_list[i-1]) > abs(seq_list[i]):
            seq_list[i] = seq_list[i-1] + seq_list[i]



if __name__ == '__main__':
    test_seq = '((())()((())))(())'
    count_bracket(test_seq)
