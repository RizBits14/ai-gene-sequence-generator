class UtilityCalculator:
    def __init__(self, target_sequence, weights):
        self.target_sequence = target_sequence
        self.weights = weights

    def calculate(self, gene_sequence):
        utility = 0
        maximum_length = max(len(gene_sequence), len(self.target_sequence))

        for i in range(maximum_length):
            if i < len(gene_sequence):
                gene_character = ord(gene_sequence[i])
            else:
                gene_character = 0
            
            if i < len(self.target_sequence):
                target_character = ord(self.target_sequence[i])
            else:
                target_character = 0
            
            if i < len(self.weights):
                weight = self.weights[i]
            else:
                weight = 1

            utility -= weight * abs(gene_character - target_character)

        return utility


class AlphaBetaPruning:
    def __init__(self, pool, utility_calculator):
        self.initial_pool = pool
        self.utility_calculator = utility_calculator

    def solve(self):
        return self._alpha_beta("", self.initial_pool, True, float('-inf'), float('inf'))

    def _alpha_beta(self, sequence, pool, maximizing, alpha, beta):
        if not pool:
            return self.utility_calculator.calculate(sequence), sequence

        if maximizing:
            max_eval = float('-inf')
            best_sequence = None
            for i in range(len(pool)):
                nucleotide = pool[i]
                new_sequence = sequence + nucleotide
                new_pool = pool[:i] + pool[i + 1:]
                eval_score, candidate_sequence = self._alpha_beta(new_sequence, new_pool, False, alpha, beta)

                if eval_score > max_eval:
                    max_eval = eval_score
                    best_sequence = candidate_sequence

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
                new_pool = pool[:i] + pool[i + 1:]
                eval_score, candidate_sequence = self._alpha_beta(new_sequence, new_pool, True, alpha, beta)

                if eval_score < min_eval:
                    min_eval = eval_score
                    best_sequence = candidate_sequence

                beta = min(beta, min_eval)
                if beta <= alpha:
                    break
            return min_eval, best_sequence


class GeneSequence:
    def __init__(self, pool, target, weights):
        self.utility_calculator = UtilityCalculator(target, weights)
        self.solver = AlphaBetaPruning(pool, self.utility_calculator)

    def run(self):
        utility_score, best_sequence = self.solver.solve()
        print(f"Best gene sequence generated: {best_sequence}")
        print(f"Utility score: {utility_score}")


# Sample Input
pool = ["A", "T", "C", "G"]
target = "GCAT"
weights = [8, 8, 1, 1]
game = GeneSequence(pool, target, weights)
game.run()