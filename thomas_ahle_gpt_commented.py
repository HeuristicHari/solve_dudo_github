
from collections import Counter 
from collections import defaultdict 
import itertools  
import sys 
from ortools.linear_solver import pywraplp 
from functools import lru_cache  
import time 


if len(sys.argv) < 4:
    print('Run {} [dice1] [dice2] [sides] mode'.format(sys.argv[0]))
    sys.exit()
else:
    DICE1 = int(sys.argv[1])
    DICE2 = int(sys.argv[2])
    SIDES = int(sys.argv[3])
    #FROMGPT: Parses command line arguments to set the number of dice and sides.


NORMAL, JOKER= range(2)
if len(sys.argv) >= 5:
    mode = {'normal': NORMAL, 'joker': JOKER}[sys.argv[4]]
else:
    mode = NORMAL
#FROMGPT: Determines the game mode based on a command line argument; defaults to NORMAL if not provided.


################################################################
# Game definition
################################################################

if mode == NORMAL:
    CALLS = [(count, side)
            for count in range(1, DICE1+DICE2+1)
            for side in range(1, SIDES+1)]
    #FROMGPT: For NORMAL mode, generates all possible calls (moves) using all dice counts and side values.
if mode == JOKER:
    # With jokers we can't call 1
    CALLS = [(count, side)
        for count in range(1, DICE1+DICE2+1)
        for side in range(2, SIDES+1)]
    #FROMGPT: In JOKER mode, the valid calls start at side 2 since side 1 is excluded.


ROLLS1 = list(itertools.product(range(1, SIDES+1), repeat=DICE1))
ROLLS2 = list(itertools.product(range(1, SIDES+1), repeat=DICE2))
#FROMGPT: Generates all possible outcomes for the two sets of dice (ROLLS1 for player 1 and ROLLS2 for player 2) using Cartesian product.

SNYD = None
#FROMGPT: Defines a special marker (SNYD) used to signal a particular game condition (e.g., ending the bidding).

def possible_calls(hist):
    if not hist:
        yield from CALLS
        return
    if hist[-1] is SNYD:
        return
    for call in CALLS:
        if call > hist[-1]:
            yield call
    yield SNYD
#FROMGPT: This function generates the valid calls (moves) after a given history `hist`. It yields all calls that are higher than the previous call and finally yields SNYD as a terminal move.

def is_correct_call(d1, d2, call):
    count, side = call
    if mode == JOKER:
        d1 = tuple(side if d == 1 else d for d in d1)
        d2 = tuple(side if d == 1 else d for d in d2)
    return not bool(Counter({side: count}) - Counter(d1 + d2))
#FROMGPT: Checks whether a call is valid based on the current dice outcomes (d1 and d2) and game mode adjustments. Uses Counter subtraction to verify that the count of a side meets the call's requirements.

def is_leaf(hist):
    assert not hist or hist[0] is not SNYD, "SNYD can't be first call"
    return hist and hist[-1] is SNYD
#FROMGPT: Determines if the history represents an end state (leaf) by checking if the last call is SNYD. Also asserts SNYD is not the first call.

def histories(hist=()):
    yield hist
    if not is_leaf(hist):
        for call in possible_calls(hist):
            yield from histories(hist + (call,)) #one element tuple: empty comma
#FROMGPT: Recursively generates all possible game histories starting from the current history. Continues until reaching a leaf state (where SNYD is the last call).

# xs (ys) are the states after a move by player + the root.
# Each of these are given a variable, since they are either leafs or parents to leafs.
# This is of course game specific, so maybe it's a bad way to do it...
xs = lambda: (h for h in histories() if not h or len(h) % 2 == 1)
ys = lambda: (h for h in histories() if not h or len(h) % 2 == 0)
#FROMGPT: Defines generators for game states: xs for states after player 1's move (or root) and ys for states after player 2's move (or root).

# The inners are those with children. Note the root is present in both.
inner_xs = lambda: (h for h in xs() if not is_leaf(h))
inner_ys = lambda: (h for h in ys() if not is_leaf(h))
#FROMGPT: Further filters xs and ys to include only non-leaf (inner) game states that still have continuation moves.

def score(d1, d2, hist):
    ''' Get the score in {-1,1} relative to player 1 '''
    assert hist and hist[-1] is SNYD
    # Player 2 called snyd
    if len(hist) % 2 == 0:
        res = is_correct_call(d1, d2, hist[-2])
    # Player 1 called snyd
    else:
        res = not is_correct_call(d1, d2, hist[-2])
    return int(res) * 2 - 1
#FROMGPT: Computes the game outcome score from player 1's perspective. The score is determined by the correctness of the call made just before SNYD, with a binary outcome (-1 or 1).

################################################################
# Write as matrix
################################################################

def initSolver():
    t = time.time()
    print('Creating variables...', file=sys.stderr)
    solver = pywraplp.Solver('', pywraplp.Solver.GLOP_LINEAR_PROGRAMMING)
    zvs = {(d2, x): solver.NumVar(-solver.infinity(), solver.infinity(),
        'z{}h{}'.format(d2, shist(x))) for x in inner_xs() for d2 in ROLLS2}
    xvs = {(d1, x): solver.NumVar(0, 1,
        'x{}h{}'.format(d1, shist(x))) for x in xs() for d1 in ROLLS1}
    print('Took {}s'.format(time.time() - t), file=sys.stderr)
    #FROMGPT: Initializes the linear programming solver and creates decision variables:
    # zvs for player 2’s value functions and xvs for player 1’s probability strategies, indexed by game state.

    t = time.time()
    print('Setting constraints 1', file=sys.stderr)

    # The thing to maximize: f.T@z
    objective = solver.Objective()
    for d2 in ROLLS2:
        objective.SetCoefficient(zvs[d2, ()], 1)
    objective.SetMaximization()
    #FROMGPT: Configures the objective function to maximize the sum of selected zvs variables at the root state for each dice outcome of player 2.

    print('Took {}s'.format(time.time() - t), file=sys.stderr)

    t = time.time()
    print('Setting constraints 2', file=sys.stderr)

    # Equalities: Ex = e
    for d1 in ROLLS1:
        # The root sums to 1
        constraint = solver.Constraint(1, 1)
        constraint.SetCoefficient(xvs[d1, ()], 1)
        #FROMGPT: Enforces that for each dice outcome in ROLLS1, the probability assigned at the root state sums to 1.
        # Make siblings sum to their parent
        for hist in inner_ys():
            constraint = solver.Constraint(0, 0)
            constraint.SetCoefficient(xvs[d1, hist[:-1]], 1)
            for call in possible_calls(hist):
                constraint.SetCoefficient(xvs[d1, hist + (call,)], -1)
            #FROMGPT: Ensures that for every inner node, the sum of the probabilities of its child nodes equals the probability of the parent node.

    print('Took {}s'.format(time.time() - t), file=sys.stderr)

  

    t = time.time()
    print('Setting constraints 3', file=sys.stderr)

    # Bound zT@F - xT@A <= 0
    # Bound F.T@z - A.Tx <= 0
    # z@F0 - x@A0 >= 0, ...
    for d2 in ROLLS2:
        for hist in ys():
            constraint = solver.Constraint(-solver.infinity(), 0)
            if hist == ():
                constraint.SetCoefficient(zvs[d2, ()], 1)
                for call in possible_calls(hist):
                    constraint.SetCoefficient(zvs[d2, hist + (call,)], 1)
                #FROMGPT: For the root state, sets constraints only on the zvs variables since there is no corresponding A matrix contribution.
                continue
            constraint.SetCoefficient(zvs[d2, hist[:-1]], -1)
            for call in possible_calls(hist):
                child = hist + (call,)
                if not is_leaf(child):
                    constraint.SetCoefficient(zvs[d2, child], 1)
            #FROMGPT: For non-root states, constructs constraints that link the zvs variables across parent and child histories.
            lhist = hist + (SNYD,) if hist[-1] is not SNYD else hist
            xhist = hist + (SNYD,) if hist[-1] is not SNYD else hist[:-1]
            for d1 in ROLLS1:
                constraint.SetCoefficient(xvs[d1, xhist], -score(d1, d2, lhist))
            #FROMGPT: Incorporates the xvs variables into the constraint with coefficients determined by the game score.
    print('Took {}s'.format(time.time() - t), file=sys.stderr)

    return solver, xvs, zvs
#FROMGPT: The initSolver() function sets up the complete linear programming model with decision variables and constraints, then returns the solver and variable dictionaries.


# Formatting of solution
def scall(call):
    if call is None:
        return "snyd"
    return '{}{}'.format(*call)
#FROMGPT: Formats a call into a string. If the call is None, returns "snyd"; otherwise, unpacks the tuple into a string.

def shist(hist):
    return ','.join(map(scall, hist))
#FROMGPT: Converts a history (sequence of calls) into a comma-separated string using scall().

def sfrac(val):
    return str(val)  #bro fuck fractions
#FROMGPT: Returns the string representation of a value; originally intended to use the fractions module but now just converts the value directly.


class CounterStrategy:
    def __init__(self, xvs):
        self.xvs = xvs
    #FROMGPT: The CounterStrategy class encapsulates methods for determining probabilities and optimal moves based on the LP solution.

    @lru_cache(maxsize=10**5)
    def findCallProb(self, d1, hist):
        ''' Return the probability that player 1 did the last move of hist '''
        assert len(hist) % 2 == 1
        xhis = self.xvs[d1, hist].solution_value()
        xpar = self.xvs[d1, hist[:-2]].solution_value()
        return xhis / xpar if xpar > 1e-10 else 0
    #FROMGPT: Calculates the probability that player 1 made the last move at a given history by using the cached LP solution values; employs division to compare the child and parent node probabilities.

    @lru_cache(maxsize=10**5)
    def findP2Call(self, d2, hist):
        ''' Find the best call for p2, choosing the optimal deterministic counter strategy '''
        assert len(hist) % 2 == 1
        if sum(self.findCallProb(d1, hist) for d1 in ROLLS1) < 1e-6:
            return next(possible_calls(hist))
        pd1s = self.estimateP1Rolls(hist)
        return min(possible_calls(hist), key=lambda call:
                sum(p * self.stateValue(d1, d2, hist + (call,)) for p, d1 in zip(pd1s, ROLLS1)))
    #FROMGPT: Determines the optimal call for player 2 by considering the expected state values for each possible move and choosing the one that minimizes player 1’s advantage.

    @lru_cache(maxsize=10**5)
    def stateValue(self, d1, d2, hist):
        ''' Return expected payoff for player 1 '''
        if hist and hist[-1] is SNYD:
            res = score(d1, d2, hist)
        # Player 1's turn to move
        elif len(hist) % 2 == 0:
            res = sum(self.stateValue(d1, d2, hist + (call,))
                      * self.findCallProb(d1, hist + (call,))
                      for call in possible_calls(hist))
        # Player 2's turn to move
        elif len(hist) % 2 == 1:
            p2call = self.findP2Call(d2, hist)
            res = self.stateValue(d1, d2, hist + (p2call,))
        return res
    #FROMGPT: Recursively computes the expected payoff (state value) for player 1 from the current history by considering the actions of both players and termination conditions.

    @lru_cache(maxsize=10**5)
    def estimateP1Rolls(self, hist):
        assert len(hist) % 2 == 1
        prob_hist_given_d = [self.findCallProb(d1, hist) for d1 in ROLLS1]
        if sum(prob_hist_given_d) < 1e-10:
            return [1 / len(ROLLS1) for _ in ROLLS1]
        return [p / sum(prob_hist_given_d) for p in prob_hist_given_d]
    #FROMGPT: Estimates the probability distribution over player 1's dice rolls given the history using a simple Bayesian normalization of the call probabilities.


def printTrees(cs):
    print('Trees:')
    for d1 in ROLLS1:
        for hist in histories():
            if not hist:
                avgValue = sfrac(sum(cs.stateValue(d1, d2, ()) for d2 in ROLLS2) / len(ROLLS2))
                values = ', '.join(sfrac(cs.stateValue(d1, d2, ())) for d2 in ROLLS2)
                print('Roll: {}, Expected: {}, Values: {}'.format(d1, avgValue, values))
                continue
            if any(cs.findCallProb(d1, hist[:j]) < 1e-8 for j in range(1, len(hist) + 1, 2)):
                continue
            s = '|  ' * len(hist) + (scall(hist[-1]) if hist else 'root')
            if len(hist) % 2 == 1:
                prob = sfrac(cs.findCallProb(d1, hist))
                print('{} p={}'.format(s, prob))
            else:
                tag = ''.join('_*'[hist[-1] == cs.findP2Call(d2, hist[:-1])]
                              for d2 in ROLLS2)
                print(s, tag)
#FROMGPT: This function prints a representation of the game tree, showing the expected values and decision probabilities for each history, to help with debugging and analysis.

def main():
    print('Setting up linear program', file=sys.stderr)
    solver, xvs, zvs = initSolver()

    t = time.time()
    print('Solving', file=sys.stderr)

    status = solver.Solve()
    if status != solver.OPTIMAL:
        print('Status:', status, file=sys.stderr)
        print(zvs[(0,), ()].solution_value())
        return

    print('Took {}s'.format(time.time() - t), file=sys.stderr)

    cs = CounterStrategy(xvs)
    printTrees(cs)

    res = sum(zv.solution_value() for (_, hist), zv in zvs.items() if hist == ())
    res /= len(ROLLS1) * len(ROLLS2)
    print('Value:', sfrac(res))

    res2 = sum(cs.stateValue(d1, d2, ()) for d1 in ROLLS1 for d2 in ROLLS2) / len(ROLLS1) / len(ROLLS2)
    print('Score:', sfrac(res2))
#FROMGPT: The main() function orchestrates the setup and solution of the linear programming model, then calculates and prints the overall game value and score based on the computed strategies.

if __name__ == '__main__':
    main()
#FROMGPT: Ensures main() is called only when this script is executed directly, not when imported as a module.
