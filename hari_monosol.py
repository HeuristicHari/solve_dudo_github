
from ortools.linear_solver import pywraplp


max_hist = 12 #1000000000000 is BS call #13 for 1v1
def pop_even(x):
    val=x
    tmp=True
    while (val>0):
        if val % 2 == 1:
            tmp = not tmp
        val = val // 2
    return tmp
def strip_one_off(x):
    if x == 0:
        return 0 #FML i can't do without this
    return x & ((1 << (x.bit_length() - 1)) - 1)

def eval(d1, d2, h): #I think should be minimum value p2 is guaranteed
    if h& (1<<max_hist) == 0:
        print(":(") 
        return
    hp = strip_one_off(h)
    max_call = hp.bit_length()
    die = max_call % 6
    quantity = 1 + max_call / 6 #TODO: NOT INT DIV
    #dith digit in base 5
    nonbs_call_correct = ((d1 / (5**die)) % 5) + ((d2 / (5**die)) % 5) >= quantity
    return 2 * (1 - (int) (nonbs_call_correct ^ pop_even(h)))-1 #val to p1


def main():
    priv_p1={}
    priv_p2={}

    player, d1, d2= map(int, input().split())

    #base 5 encoding of length 6
    if d1 == 1:
        for i in range (6):
            priv_p1[5**i]=1/6
    if d2 == 1:
        for i in range(6):
            priv_p2[5**i]=1/6 
    if d1 == 2:
        for i in range (6):
            for j in range(6):
                if (5**i + 5**j not in priv_p1):
                    priv_p1[5**i + 5**j] = 1/36
                else:
                    priv_p1[5**i + 5**j] += 1/36
    if d2 == 2:
        for i in range (6):
            for j in range(6):
                if (5**i + 5**j not in priv_p2):
                    priv_p2[5**i + 5**j] = 1/36
                else:
                    priv_p2[5**i + 5**j] += 1/36



    
    solver=pywraplp.Solver.CreateSolver("GLOP") 
    if not solver:
        return
    x = []


    primals = {}
    duals = {}
    even_histories = []
    odd_histories = [0] #idk why but apparently i gotta do this
    for i in range(1<<(max_hist+1)):
        if pop_even(i):
            even_histories.append(i)
        else:
            odd_histories.append(i)

     
    for history in odd_histories:
        for d1 in priv_p1:
            primals[(d1, history)] = solver.NumVar(0, 1, "")
        for d2 in priv_p2:
            if (history & (1 << max_hist)) == 0: #if not terminal
                duals[(d2, history)] = solver.NumVar(-solver.infinity(), solver.infinity(), "")
    objective = solver.Objective()
    for d2 in priv_p2:
        objective.SetCoefficient(duals[(d2, 0)], priv_p2[d2])
    objective.SetMaximization()
    #for something
        #solver.add(summation >= v)
    for d1 in priv_p1:
        solver.Add(primals[(d1, 0)] == priv_p1[d1] )
        for history in even_histories:
            p=strip_one_off(history)
            history_children = []

            for added_bit in range(max_hist):
                history_children.append(history | (1 << added_bit))
            if (history != 0):
                history_children.append(history | (1 << max_hist))
            solver.Add(primals[(d1, p)] == sum(primals[(d1, child)] for child in history_children))
    for d2 in priv_p2:
        for history in even_histories:
            history_children = []
            for added_bit in range(history.bit_length(), max_hist):#excludes BS
                history_children.append(history | (1 << added_bit))
            if history == 0:
                solver.Add(duals[(d2, history)] + sum(duals[(d2, child)] for child in history_children) <= 0)
                continue
            summation = sum(duals[d2,child] for child in history_children)
            p=strip_one_off(history)
            to_eval = history | (1 << max_hist) #BS call
            to_index = history ^ (1 << max_hist)#this gives an even history, toggles bs
            summation -= sum(primals[d1, to_index] * eval(d1, d2, to_eval) for d1 in priv_p1)
            solver.add(duals[(d2, p)] >= summation)

    status = solver.Solve()
    if status != pywraplp.Solver.OPTIMAL:
        print("yo wtf")
        return
    import csv

    # … after solver.Solve() and status check …

    with open('p1_strategy.csv','w', newline='') as fout:
        writer = csv.writer(fout)
        writer.writerow(['d1','history','entry_prob'])

        for (d1, h), var in primals.items():
            # only P1‐nodes (odd‐length histories)
            if not pop_even(h):
                # remove last two bits to get the P1‐parent
                p = strip_one_off(strip_one_off(h))
                pv = primals[(d1, p)].solution_value()
                if pv < 1e-8:
                    prob = 0.0
                else:
                    prob = var.solution_value() / pv
                writer.writerow([d1, h, f"{prob:.6f}"])
 

if __name__ == "__main__":
    main()