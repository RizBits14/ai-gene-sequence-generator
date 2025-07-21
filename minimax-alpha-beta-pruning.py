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


class GeneSequenceGame:
    def __init__(self, pool, target, weights):
        self.pool = pool
        self.target = target
        self.weights = weights

    def alpha_beta(self, sequence, pool, maximizing, alpha, beta):
        if len(pool) == 0:
            return calculate_utility(sequence, self.target, self.weights), sequence
        
        if maximizing:
            max_eval = float('-inf')
            best_sequence = None
            for i in range(len(pool)):
                nucleotide = pool[i]
                new_sequence = sequence + nucleotide
                new_pool = pool[:i] + pool[i+1:]
                eval_score, _ = self.alpha_beta(new_sequence, new_pool, False, alpha, beta)

                if eval_score > max_eval:
                    max_eval = eval_score
                    best_sequence = new_sequence

                alpha = max(alpha, max_eval)
                if alpha >= beta:
                    break
            return max_eval, best_sequence
        
        else:
            min_eval = float('inf')
            best_sequence = None
            for i in range(len(pool)):
                nucleotide = pool[i]
                new_sequence = sequence + nucleotide
                new_pool = pool[:i] + pool[i+1:]
                eval_score, _ = self.alpha_beta(new_sequence, new_pool, True, alpha, beta)

                if eval_score < min_eval:
                    min_eval = eval_score
                    best_sequence = new_sequence

                beta = min(beta, min_eval)
                if beta <= alpha:
                    break
            return min_eval, best_sequence
