n = 104
with open('combinetable.pdb70', 'r', encoding='utf-8') as gt_file:
    with open('combinetable.pdb70b', 'w', encoding='utf-8') as modilfied_gt_file:
        for line in gt_file:
            parts = line.split('\t')
            if(len(parts) > 2):
                parts[2] = parts[2].replace(' d.45.1.1', '')
                parts[1] = parts[1]
                left_part = parts[1][:n] 
                right_part = parts[1][n+1:]
                parts[1] = left_part + right_part
                modilfied_gt_file.write('\t'.join(parts))
