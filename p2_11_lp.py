
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
    if (x == 0):
        return 0
    return x - (1 << (x.bit_length() - 1))

def eval(d1, d2, h): #player one value
    if h& (1<<max_hist) == 0:
        print(":(") 
        return
    hp = strip_one_off(h)
    max_call = hp.bit_length()-1
    die = max_call % 6
    quantity = 1 + max_call // 6  
    #dith digit in base 7
    nonbs_call_correct = ( ((d1 // (7**die)) % 7) + ((d2 // (7**die)) % 7) ) >= quantity
    if nonbs_call_correct and pop_even(h) or not (nonbs_call_correct or pop_even(h)):
        return 1
    else:
        return -1


def main():
    global max_hist
    priv_p1={}
    priv_p2={}

    d1, d2= map(int, input().split())

    max_hist = 6 * (d1 + d2)

    #base 7 encoding of length 6
    if d1 == 1:
        for i in range (6):
            priv_p1[7**i]=1/6
    if d2 == 1:
        for i in range(6):
            priv_p2[7**i]=1/6 
    if d1 == 2:
        for i in range (6):
            for j in range(6):
                if (7**i + 7**j not in priv_p1):
                    priv_p1[7**i + 7**j] = 1/36
                else:
                    priv_p1[7**i + 7**j] += 1/36
    if d2 == 2:
        for i in range (6):
            for j in range(6):
                if (7**i + 7**j not in priv_p2):
                    priv_p2[7**i + 7**j] = 1/36
                else:
                    priv_p2[7**i + 7**j] += 1/36



    
    solver=pywraplp.Solver.CreateSolver("GLOP") 
    if not solver:
        return
    y = []
    p = []
    d1_hist_to_xindex = {}
    d2_hist_to_yindex = {} 
    yindex_to_d2_hist = {}


    y.append(solver.NumVar(1, 1, "yphi"))
    d1_hist_to_xindex[(0, 0)] = 0 #die 0: _. hist 0: nocalls
    d2_hist_to_yindex[(0, 0)] = 0 #
    yindex_to_d2_hist[0] = (0, 0)
    

    for history in range(1, 1<<(max_hist+1)):
        if (history == 1 << max_hist): #ban empty bs
            continue
        if pop_even(history):
            for d2 in priv_p2:
                y.append(solver.NumVar(0, 1, ""))
                curr_len = len(d2_hist_to_yindex)
                d2_hist_to_yindex[(d2, history)] = curr_len
                yindex_to_d2_hist[curr_len] = (d2, history)
                #it got appended, need to minus 1
        else:
            for d1 in priv_p1:
                curr_len = len(d1_hist_to_xindex)
                d1_hist_to_xindex[(d1, history)] = curr_len
    #implement flow constraints elsewhere -- lowkey chill


    #let's first do the flow constraints. Fy=f.
    solver.Add(y[0] == 1)
    for d2 in priv_p2:
        for move_1 in range(max_hist):
            solver.Add(-y[0] + sum(y[d2_hist_to_yindex[(d2, (1<<move_1) + (1<<move_2))]]
            for move_2 in range(move_1+1, max_hist+1))==0)
    for middle_history in range(1, 1<<max_hist):
        parent = strip_one_off(middle_history)
        if pop_even(middle_history) or (parent == 0):
            continue
        max_move_so_far = middle_history.bit_length()
        for d2 in priv_p2:
            solver.Add(-y[d2_hist_to_yindex[(d2, parent)]] + 
            sum(y[d2_hist_to_yindex[(d2, middle_history + (1<<move))]] 
            for move in range(max_move_so_far, max_hist+1)) == 0)
    
    #now, -Ay + E^Tp \geq 0

    #first, build p

    d1_yhistory_to_constraintindex = {}
    #put another way, we iterate through the variables
    p.append(solver.NumVar(-solver.infinity(), solver.infinity(), ""))#"root=1" constraint variable
    d1_yhistory_to_constraintindex[(0, 0)] = len(p)-1 #=0
    for d1 in priv_p1: #empty set splits into |d1| constraints, saying it splits
        p.append(solver.NumVar(-solver.infinity(), solver.infinity(), ""))
        d1_yhistory_to_constraintindex[(d1, 0)] = len(p) - 1
            #p variable corresponding to the constraint that the empty set splits into |d1| constraints
    for (d2, yhistory) in d2_hist_to_yindex.keys(): #inner_nodes
        if yhistory == 0 or yhistory & (1 << max_hist) != 0:
            continue
        if d2 != min(priv_p2): #just one for each history
            continue
        
        #xhistory = strip_one_off(yhistory) #just to clarify: this is what p really corresponds to
        for d1 in priv_p1:
            p.append(solver.NumVar(-solver.infinity(), solver.infinity(), ""))
            d1_yhistory_to_constraintindex[(d1, yhistory)] = len(p) - 1
    #now, p had an entry for every column of Et, so E^t p makes sense




    #now, we build E^t

    Et = [[] for _ in d1_hist_to_xindex] #variable to constraints where it exists
    
    #let's build by constraint.



    Et[0].append(('+', 0)) #root=1
  
    constraint_index = 0
    for d1 in priv_p1:
        constraint_index = d1_yhistory_to_constraintindex[(d1, 0)] #+1 because root is 0
        Et[0].append(('-', constraint_index)) #die splits into |d1| constraints
        for move_1 in range(max_hist):
            Et[d1_hist_to_xindex[(d1, (1<<move_1))]].append(('+', constraint_index))
    
    for (d2, yhistory) in d2_hist_to_yindex.keys():
        if yhistory == 0 or yhistory & (1 << max_hist) != 0: #needs to flow
            continue
        if d2!=min(priv_p2): #we dont care about d2, just want one for each history
            continue
        xhistory = strip_one_off(yhistory) #parent #since its a nonzero y history, has at least 2 
        #so xhistory is nonzero
        for d1 in priv_p1:
            constraint_index = d1_yhistory_to_constraintindex[(d1, yhistory)]
            Et[d1_hist_to_xindex[(d1, xhistory)]].append(('-', constraint_index))
            for move_1 in range(yhistory.bit_length(), max_hist+1):
                Et[d1_hist_to_xindex[(d1, yhistory + (1<<move_1))]].append(('+', constraint_index))





    for (d1, xhistory), x_index in d1_hist_to_xindex.items():
        summation = 0
        #enforce -Ay_{x_index} + E^Tp_{x_index} \geq 0
        yhistory = xhistory ^ (1 << max_hist) #the dude so that x,yhistories uniquely encode 
        #the relevant BS node
        if (xhistory != 0):
            for d2 in priv_p2:
                y_index = d2_hist_to_yindex[(d2, yhistory)]
                summation -= priv_p1[d1] * priv_p2[d2] * eval(d1, d2, max(xhistory, yhistory)) * y[y_index]
            #summation += E^tp_{x_index}'
        for (sign, index) in Et[x_index]:
            if sign == '+':
                summation += p[index]
            else:
                summation -= p[index]
        solver.Add(summation >= 0) #constraint (-Ay + E^Tp)_xindex \geq 0
    solver.Minimize(p[0]) #minimize e^tp
  



    status = solver.Solve()
    if status != pywraplp.Solver.OPTIMAL:
        print("yo wtf")
        return
    #dimensionality testing
   
    # print("numxhist: "+str(len(d1_hist_to_xindex)))
    # #seems legit
    # print ("numyhist: "+str(len(d2_hist_to_yindex)))
    # for ele in Et:
    #     for (sign, index) in ele:
    #         suit.add(index)
    # print("numconstraints: "+str(len(suit)) + "should be size p" + str(len(p)))
    #p size is correct
    # print(suit)
  
    
    print (f"{solver.Objective().Value():.10f}")
    for i in range(len(y)):
        num = y[i].solution_value()
        if num < 0.00000001:
            continue
        (d2, history) = yindex_to_d2_hist[i]
        grandparent = strip_one_off(strip_one_off(history))
        if grandparent == 0:
            den_index = d2_hist_to_yindex[(0, grandparent)]
        else:
            den_index = d2_hist_to_yindex[(d2, grandparent)]
        den = y[den_index].solution_value()
        val = 0
        if den > 0.00000001:
            val = num / den
        if (val > 0.00000001) and (d2 != 16807 and (history.bit_length() <=12) and (history.bit_length() >=7)):
            print (f"BELOW LINE TRASH occurs prob {num:.3}")
        (d2, history) = yindex_to_d2_hist[i]
        print (f"{d2} {history:b}: {val:.3}")


 
 

if __name__ == "__main__":
    main()