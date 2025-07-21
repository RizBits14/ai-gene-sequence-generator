class GeneSequenceWithBooster:
    def __init__(self, pool, target, weights, student_id_first_two):
        self.original_pool = pool
        self.target = target
        self.weights = weights
        self.booster_multiplier = student_id_first_two / 100.0

    def run(self):
        # Normal run (without booster 'S')
        normal_game = GeneSequence(
            [p for p in self.original_pool if p != 'S'],
            self.target,
            self.weights
        )
        normal_utility_score, normal_best_sequence = normal_game.solver.solve()

        # Booster run (including 'S' with booster multiplier)
        booster_index = None
        for idx, nucleotide in enumerate(self.original_pool):
            if nucleotide == 'S':
                booster_index = idx
                break
        
        booster_game = GeneSequenceBoosterVersion(
            self.original_pool, self.target, self.weights, self.booster_multiplier, booster_index
        )
        booster_utility_score, booster_best_sequence = booster_game.run()

        if booster_utility_score > normal_utility_score:
            print("YES")
        else:
            print("NO")
        print("With special nucleotide")
        print(f"Best gene sequence generated: {booster_best_sequence}")
        print(f"Utility score: {round(booster_utility_score, 2)}")
