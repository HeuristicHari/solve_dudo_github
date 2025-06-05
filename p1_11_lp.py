
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
    x = []
    q = []
    d1_hist_to_xindex = {}
    d2_hist_to_yindex = {} 
    xindex_to_d1_hist = {}


    x.append(solver.NumVar(1, 1, "yphi"))
    d1_hist_to_xindex[(0, 0)] = 0 #die 0: _. hist 0: nocalls
    d2_hist_to_yindex[(0, 0)] = 0 #
    xindex_to_d1_hist[0] = (0, 0)
    

    for history in range(1, 1<<(max_hist+1)):
        if (history == 1 << max_hist): #ban empty bs
            continue
        if pop_even(history):
            for d2 in priv_p2:
                curr_len = len(d2_hist_to_yindex)
                d2_hist_to_yindex[(d2, history)] = curr_len
        else:
            for d1 in priv_p1:
                x.append(solver.NumVar(0, 1, ""))
                curr_len = len(d1_hist_to_xindex)
                d1_hist_to_xindex[(d1, history)] = curr_len
                xindex_to_d1_hist[curr_len] = (d1, history)


    #let's first do the flow constraints. Fy=f.
    solver.Add(x[0] == 1)
    for d1 in priv_p1:
        solver.Add(-x[0] + sum(x[d1_hist_to_xindex[(d1, 1<<move_1)]]
        for move_1 in range(max_hist))==0)
    for middle_history in range(1, 1<<max_hist):
        parent = strip_one_off(middle_history)
        if not pop_even(middle_history):
            continue
        max_move_so_far = middle_history.bit_length()
        for d1 in priv_p1:
            solver.Add(-x[d1_hist_to_xindex[(d1, parent)]] + 
            sum(x[d1_hist_to_xindex[(d1, middle_history + (1<<move))]] 
            for move in range(max_move_so_far, max_hist+1)) == 0)
    
    #now, A^Tx + F^Tq \geq 0

    #first, build q

    d2_xhistory_to_constraintindex = {}
    #put another way, we iterate through the variables
    q.append(solver.NumVar(-solver.infinity(), solver.infinity(), ""))#"root=1" constraint variable
    d2_xhistory_to_constraintindex[(0, 0)] = len(q)-1 #=0
    for (d1, xhistory) in d1_hist_to_xindex.keys(): #inner_nodes
        if xhistory == 0 or xhistory & (1 << max_hist) != 0:
            continue
        if d1 != min(priv_p1): #just one for each history
            continue
        for d2 in priv_p2:
            q.append(solver.NumVar(-solver.infinity(), solver.infinity(), ""))
            d2_xhistory_to_constraintindex[(d2, xhistory)] = len(q) - 1




    #now, we build F^t

    Ft = [[] for _ in d2_hist_to_yindex] #variable to constraints where it exists
    
    #let's build by constraint.



    Ft[0].append(('+', 0)) #root=1
  
    
    for (d1, xhistory) in d1_hist_to_xindex.keys():
        if xhistory == 0 or xhistory & (1 << max_hist) != 0: #needs to flow
            continue
        if d1!=min(priv_p1): #we dont care about d1, just want one for each history
            continue
        yhistory = strip_one_off(xhistory) #parent #since its a nonzero y history, has at least 2 
        #so xhistory is nonzero
        for d2 in priv_p2:
            constraint_index = d2_xhistory_to_constraintindex[(d2, xhistory)]
            if yhistory == 0:
                Ft[d2_hist_to_yindex[(0, yhistory)]].append(('-', constraint_index)) #0,0
            else:
                Ft[d2_hist_to_yindex[(d2, yhistory)]].append(('-', constraint_index))
            for move_1 in range(xhistory.bit_length(), max_hist+1):
                Ft[d2_hist_to_yindex[(d2, xhistory + (1<<move_1))]].append(('+', constraint_index))





    for (d2, yhistory), y_index in d2_hist_to_yindex.items():
        summation = 0
        #enforce A^Tx_{y_index} + F^Tq_{y_index} \geq 0
        xhistory = yhistory ^ (1 << max_hist) #the dude so that x,yhistories uniquely encode 
        #the relevant BS node
        if (yhistory != 0):
            for d1 in priv_p1:
                x_index = d1_hist_to_xindex[(d1, xhistory)]
                summation += priv_p1[d1] * priv_p2[d2] * eval(d1, d2, max(xhistory, yhistory)) * x[x_index]
                #+ as its A^Tx + F^T q \geq 0

        #summation += F^tp_{x_index}'
        for (sign, index) in Ft[y_index]:
            if sign == '+':
                summation += q[index]
            else:
                summation -= q[index]
        solver.Add(summation >= 0) #constraint (-Ay + E^Tp)_xindex \geq 0

   


    solver.Maximize(-q[0]) #maximize -q^Tf
  



    status = solver.Solve()
    if status != pywraplp.Solver.OPTIMAL:
        print("yo wtf")
        print (f"Status: {status}")
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
    for i in range(len(x)):
        num = x[i].solution_value()
        if num < 0.00000001:
            continue
        (d1, history) = xindex_to_d1_hist[i]
        grandparent = strip_one_off(strip_one_off(history))
        if grandparent == 0:
            den_index = d1_hist_to_xindex[(0, grandparent)]
        else:
            den_index = d1_hist_to_xindex[(d1, grandparent)]
        den = x[den_index].solution_value()
        val = 0
        if den > 0.00000001:
            val = num / den
        if (val > 0.00000001) and ((history.bit_length() <=12) and (history.bit_length() >=7)):
            print (f"BELOW LINE TRASH occurs prob {num:.3}")
        (d1, history) = xindex_to_d1_hist[i]
        print (f"{d1} {history:b}: {val:.3}")


 
 

if __name__ == "__main__":
    import time
    #t = time.time()
    main()
    #print ('Total time:', time.time() - t, 'seconds')