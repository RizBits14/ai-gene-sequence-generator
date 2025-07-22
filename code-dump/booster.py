class UtilityCalculatorWithBooster:
    def __init__(self, target_sequence, weights, booster_index=None, multiplier=1.0):
        self.target_sequence = target_sequence
        self.weights = weights
        self.booster_index = booster_index
        self.multiplier = multiplier

    def calculate(self, gene_sequence):
        utility = 0
        max_len = max(len(gene_sequence), len(self.target_sequence))

        for i in range(max_len):
            gene_char = ord(gene_sequence[i]) if i < len(gene_sequence) else 0
            target_char = ord(self.target_sequence[i]) if i < len(self.target_sequence) else 0

            weight = self.weights[i] if i < len(self.weights) else 1

            if self.booster_index is not None and i >= self.booster_index:
                weight *= self.multiplier

            utility -= weight * abs(gene_char - target_char)

        return round(utility, 2)


class AlphaBetaSolverWithBooster:
    def __init__(self, pool, utility_calculator, enable_booster=False):
        self.initial_pool = pool
        self.utility_calculator = utility_calculator
        self.enable_booster = enable_booster

    def solve(self):
        return self._alpha_beta("", self.initial_pool, True, float('-inf'), float('inf'), 0, False)

    def _alpha_beta(self, sequence, pool, maximizing, alpha, beta, depth, booster_activated):
        if not pool:
            booster_index = sequence.index("S") if booster_activated else None
            calculator = UtilityCalculatorWithBooster(
                self.utility_calculator.target_sequence,
                self.utility_calculator.weights,
                booster_index=booster_index,
                multiplier=self.utility_calculator.multiplier
            )
            return calculator.calculate(sequence), sequence

        if maximizing:
            max_eval = float('-inf')
            best_sequence = None
            for i in range(len(pool)):
                nucleotide = pool[i]
                new_sequence = sequence + nucleotide
                new_pool = pool[:i] + pool[i + 1:]

                # Booster activates only when Agent 1 picks 'S'
                activated = booster_activated
                if self.enable_booster and nucleotide == "S":
                    activated = True

                eval_score, candidate_sequence = self._alpha_beta(
                    new_sequence, new_pool, False, alpha, beta, depth + 1, activated
                )

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

                # Minimizer does not activate booster
                eval_score, candidate_sequence = self._alpha_beta(
                    new_sequence, new_pool, True, alpha, beta, depth + 1, booster_activated
                )

                if eval_score < min_eval:
                    min_eval = eval_score
                    best_sequence = candidate_sequence

                beta = min(beta, min_eval)
                if beta <= alpha:
                    break
            return min_eval, best_sequence


class GeneSequenceWithBooster:
    def __init__(self, pool, target, student_id_digits):
        self.pool = pool
        self.target = target
        self.sid_digits = student_id_digits
        self.weights = student_id_digits[-len(target):]
        self.booster_multiplier = round((student_id_digits[0] * 10 + student_id_digits[1]) / 100, 2)

    def run(self):
        # Run without booster (remove S)
        pool_without_s = [x for x in self.pool if x != "S"]
        calc_normal = UtilityCalculatorWithBooster(self.target, self.weights)
        solver_normal = AlphaBetaSolverWithBooster(pool_without_s, calc_normal)
        score_normal, sequence_normal = solver_normal.solve()

        print("Without special nucleotide")
        print(f"Best gene sequence generated: {sequence_normal}")
        print(f"Utility score: {score_normal}")
        print()

        # Run with booster logic (if S exists)
        if "S" in self.pool:
            calc_booster = UtilityCalculatorWithBooster(self.target, self.weights, multiplier=self.booster_multiplier)
            solver_booster = AlphaBetaSolverWithBooster(self.pool, calc_booster, enable_booster=True)
            score_booster, sequence_booster = solver_booster.solve()

            print("With special nucleotide")
            print(f"Best gene sequence generated: {sequence_booster}")
            print(f"Utility score: {score_booster}")
            print()

            if score_booster > score_normal:
                print("YES")
            else:
                print("NO")
        else:
            print("Special nucleotide 'S' not found in pool.")


# ------------------------
# ✅ Sample Input for Testing
# ------------------------
nucleotide_pool = ["S", "A", "T", "G", "C"]
target_sequence = "GCAT"
student_id_digits = [2, 3, 1, 8, 8, 8, 1, 1]

game = GeneSequenceWithBooster(nucleotide_pool, target_sequence, student_id_digits)
game.run()
