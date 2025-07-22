##################### Task - 1 #####################

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

# Sample Input for Task - 1
# pool = ["A", "T", "C", "G"]
# target = "GCAT"
# weights = [8, 8, 1, 1]
# game = GeneSequence(pool, target, weights)
# game.run()

##################### Task - 2 #####################


class UtilityBoosterCalc:
    def __init__(self, target_sequence, weights, booster_idx = None, mul = 1.0):
        self.target_sequence = target_sequence
        self.weights = weights
        self.booster_idx = booster_idx
        self.mul = mul

    def calculate(self, gene_sequence):
        utility = 0
        max_len = max(len(gene_sequence), len(self.target_sequence))

        for i in range(max_len):
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

            if self.booster_idx is not None:
                if i >= self.booster_idx:
                    weight = weight * self.mul
            
            utility -= weight * abs(gene_character - target_character)

        return round(utility, 2)


class AlphaBetaBooster:
    def __init__(self, pool, utility_calculator, enable_booster = False):
        self.initial_pool = pool
        self.utility_calculator = utility_calculator
        self.enable_booster = enable_booster

    def solve(self):
        return self._alpha_beta("", self.initial_pool, True, float('-inf'), float('inf'), False)

    def _alpha_beta(self, sequence, pool, maximizing, alpha, beta, booster_activated):
        if not pool:
            if booster_activated:
                booster_idx = sequence.index("S")
            else:
                booster_idx = None

            calculator = UtilityBoosterCalc(self.utility_calculator.target_sequence,self.utility_calculator.weights,booster_idx = booster_idx,mul = self.utility_calculator.mul)

            return calculator.calculate(sequence), sequence
        
        if maximizing:
            max_eval = float('-inf')
            best_sequence = None
            for i in range(len(pool)):
                nucleotide = pool[i]
                new_sequence = sequence + nucleotide
                new_pool = pool[:i] + pool[i + 1:]

                activated = booster_activated
                if self.enable_booster and nucleotide == "S":
                    activated = True

                eval_score, candidate_sequence = self._alpha_beta(new_sequence, new_pool, False, alpha, beta, activated)

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

                eval_score, candidate_sequence = self._alpha_beta(new_sequence, new_pool, True, alpha, beta, booster_activated)

                if eval_score < min_eval:
                    min_eval = eval_score
                    best_sequence = candidate_sequence
                
                beta = min(beta, min_eval)
                if beta <= alpha:
                    break
            return min_eval, best_sequence


class GeneSequenceWithBooster:
    def __init__(self, pool, target, student_id):
        self.pool = pool
        self.target = target
        self.sid_digits = student_id
        self.weights = student_id[-len(target):]  
        self.booster_multiplier = round((student_id[0] * 10 + student_id[1]) / 100, 2)

    def run(self):
        pool_without_s = [x for x in self.pool if x != "S"]
        calc_normal = UtilityBoosterCalc(self.target, self.weights)
        solver_normal = AlphaBetaBooster(pool_without_s, calc_normal)
        score_normal, sequence_normal = solver_normal.solve()

        print("Without special nucleotide:")
        print(f"Best gene sequence generated: {sequence_normal}")
        print(f"Utility score: {score_normal}\n")

        if "S" in self.pool:
            calc_booster = UtilityBoosterCalc(self.target, self.weights, mul=self.booster_multiplier)
            solver_booster = AlphaBetaBooster(self.pool, calc_booster, enable_booster=True)
            score_booster, sequence_booster = solver_booster.solve()

            print("With special nucleotide:")
            print(f"Best gene sequence generated: {sequence_booster}")
            print(f"Utility score: {score_booster}\n")

            if score_booster > score_normal:
                print("YES")
            else:
                print("NO")
        else:
            print("Special nucleotide 'S' not found in pool.")

# Sample Input for Task 2
nucleotide_pool = ["S", "A", "T", "G", "C"]
target_sequence = "GCAT"
student_id_digits = [2, 1, 2, 0, 1, 1, 2, 9]
game = GeneSequenceWithBooster(nucleotide_pool, target_sequence, student_id_digits)
game.run()