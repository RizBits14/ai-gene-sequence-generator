##################### Without Booster #####################

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

        if maximizing == True:
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
        self.weights = weights[-len(target):]
        self.utility_calculator = UtilityCalculator(target, self.weights)
        self.solver = AlphaBetaPruning(pool, self.utility_calculator)

    def run(self):
        utility_score, best_sequence = self.solver.solve()
        print(f"Best gene sequence generated: {best_sequence}")
        print(f"Utility score: {utility_score}")
        

# Sample Input W/O Booster
# pool = ["A", "T", "C", "G"]
# target = "GCAT"
# weights = [2, 3, 1, 8, 8, 8, 1, 1]
# game = GeneSequence(pool, target, weights)
# game.run()


##################### With Booster "S" #####################
class UtilityWithBooster:
    def __init__(self, target_sequence, weights, booster_index = None, multiplier = 1.0):
        self.target_sequence = target_sequence
        self.weights = weights
        self.booster_index = booster_index
        self.multiplier = multiplier

    def calculate(self, gene_sequence):
        total_utility = 0
        loop_range = max(len(gene_sequence), len(self.target_sequence))

        for i in range(loop_range):
            gene_char = ord(gene_sequence[i]) if i < len(gene_sequence) else 0
            target_char = ord(self.target_sequence[i]) if i < len(self.target_sequence) else 0
            if i < len(self.weights):
                weight = self.weights[i]
                if self.booster_index is not None and i >= self.booster_index:
                    weight *= self.multiplier
            else:
                weight = 1

            total_utility -= weight * abs(gene_char - target_char)

        return round(total_utility, 2)
        
class BoosterAlphaBeta:
    def __init__(self, pool, utility_object, booster_enabled = False):
        self.pool = pool
        self.utility_object = utility_object
        self.booster_enabled = booster_enabled

    def run_solver(self):  
        return self._minimax("", self.pool, True, float("-inf"), float("inf"), booster_triggered = False) 
    
    def _minimax(self, current_sequence, remaining_pool, maximizing_player, alpha, beta, booster_triggered):
        if not remaining_pool:
            if booster_triggered == True:
                booster_index = current_sequence.index("S")
            else:
                booster_index = None
            calculator = UtilityWithBooster(self.utility_object.target_sequence, self.utility_object.weights, booster_index, self.utility_object.multiplier)
            
            return calculator.calculate(current_sequence), current_sequence
        
        if maximizing_player:
            max_value = float('-inf')
            best_seq = ""
            for i in range(len(remaining_pool)):
                nucleotide = remaining_pool[i]
                new_sequence = current_sequence + nucleotide
                new_pool = remaining_pool[:i] + remaining_pool[i + 1:]
                is_boosted = booster_triggered or (self.booster_enabled and nucleotide == "S")
                eval_score, eval_sequence = self._minimax(new_sequence, new_pool, False, alpha, beta, is_boosted)

                if eval_score > max_value:
                    max_value = eval_score
                    best_seq = eval_sequence

                alpha = max(alpha, eval_score)

                if alpha >= beta:
                    break

            return max_value, best_seq
        
        else:
            min_value = float('inf')
            best_seq = ""
            for i in range(len(remaining_pool)):
                nucleotide = remaining_pool[i]
                new_sequence = current_sequence + nucleotide
                new_pool = remaining_pool[:i] + remaining_pool[i + 1:]
                eval_score, eval_sequence = self._minimax(new_sequence, new_pool, True, alpha, beta, booster_triggered)

                if eval_score < min_value:
                    min_value = eval_score
                    best_seq = eval_sequence

                beta = min(beta, eval_score)

                if beta <= alpha:
                    break

            return min_value, best_seq
        
class GeneBoosterRunner:
    def __init__(self, pool, target, student_id_digits):
        self.pool = pool
        self.target = target
        self.student_id_digits = student_id_digits
        self.weights = student_id_digits[-len(target):]
        self.boost_multiplier = round((student_id_digits[0] * 10 + student_id_digits[1]) / 100, 2)

    def execute(self):
        pool_without_s = [n for n in self.pool if n != "S"]
        no_booster_calc = UtilityWithBooster(self.target, self.weights)
        no_booster_solver = BoosterAlphaBeta(pool_without_s, no_booster_calc)
        score_without_s, best_seq_without_s = no_booster_solver.run_solver()

        print("Without special nucleotide:")
        print(f"Best gene sequence generated: {best_seq_without_s}")
        print(f"Utility score: {score_without_s}\n")

        if "S" in self.pool:
            booster_calc = UtilityWithBooster(self.target, self.weights, multiplier=self.boost_multiplier)
            booster_solver = BoosterAlphaBeta(self.pool, booster_calc, booster_enabled=True)
            score_with_s, best_seq_with_s = booster_solver.run_solver()
            if score_with_s > score_without_s:
                print("YES")
            else:
                print("NO")
            print("With special nucleotide:")
            print(f"Best gene sequence generated: {best_seq_with_s}")
            print(f"Utility score: {score_with_s}\n")
        else:
            print("Special nucleotide 'S' not found in pool.")

#Sample Input W Booster
nucleotide_pool = ["S", "A", "T", "G", "C"]
target_sequence = "GCAT"
student_id_digits = [2, 3, 1, 8, 8, 8, 1, 1]
gene_game = GeneBoosterRunner(nucleotide_pool, target_sequence, student_id_digits)
gene_game.execute()