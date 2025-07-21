def calculate_utility(gene_sequence, target_sequence, weights):
    utility = 0
    maximum_length = max(len(gene_sequence), len(target_sequence))

    for i in range(maximum_length):
        if i < len(gene_sequence):
            gene_character = ord(gene_sequence[i])
        else:
            gene_character = 0
        
        if i < len(target_sequence):
            target_character = ord(target_sequence[i])
        else:
            target_character = 0
        
        if i < len(weights):
            weight = weights[i]
        else:
            weight = 1

        utility -= weight * abs(gene_character - target_character)

    return utility


gene = "ATCG"
target = "ATGC"
weights = [4, 0, 5, 2]

utility_score = calculate_utility(gene, target, weights)
print(f"Utility Score: {utility_score}")