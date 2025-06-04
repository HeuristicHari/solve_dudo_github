
from ortools.linear_solver import pywraplp
import pygambit as gbt


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

def eval(d1, d2, h): #p1 guaranteed value
    if h& (1<<max_hist) == 0:
        print(":(") 
        return
    hp = strip_one_off(h)
    max_call = hp.bit_length()
    die = max_call % 6
    quantity = 1 + max_call // 6 #TODO: NOT INT DIV
    #dith digit in base 5. TODO if you can have more than 5 of a kind, then this is cooked. 
    #but I assume I will never have more than 5 of a kind
    nonbs_call_correct = ((d1 // (7**die)) % 7) + ((d2 // (7**die)) % 7) >= quantity
    if nonbs_call_correct and pop_even(h):
        return 1
    else:
        return -1


def main():
    priv_p1={}
    priv_p2={}
    d1, d2= map(int, input().split())
    max_hist = 6 * (d1 + d2)

    #base 7 one-hot encoding of length 6. This facilitates up to 6 of a kind on 6 different dice.
    if d1 == 1:
        for i in range (6):
            priv_p1[7**i]=gbt.Rational(1,6)
    if d2 == 1:
        for i in range(6):
            priv_p2[7**i]=gbt.Rational(1,6)
    if d1 == 2:
        for i in range (6):
            for j in range(6):
                if (7**i + 7**j not in priv_p1):
                    priv_p1[7**i + 7**j] = gbt.Rational(1,36)
                else:
                    priv_p1[7**i + 7**j] += gbt.Rational(1,36)
    if d2 == 2:
        for i in range (6):
            for j in range(6):
                if (7**i + 7**j not in priv_p2):
                    priv_p2[7**i + 7**j] = gbt.Rational(1,36)
                else:
                    priv_p2[7**i + 7**j] += gbt.Rational(1,36)
    
    g = gbt.Game.new_tree(players = ["P1","P2"], 
                          title = "Perudo")
    g.append_move(g.root, g.players.chance,[str(die) for die in priv_p1] )

    g.set_chance_probs(g.root.infoset, [priv_p1[die] for die in priv_p1])

    for node in g.root.children:
        g.append_move(node, g.players.chance, [str(die) for die in priv_p2])    
        g.set_chance_probs(node.infoset, [priv_p2[die] for die in priv_p2])
    node_container = {}

    tmpFirst=0
    for roll1 in priv_p1:
        tmpSecond=0
        for roll2 in priv_p2:
            nodes = []
            for _ in range(1 << (max_hist + 1)):
                nodes.append(None)
            node_container[(roll1, roll2)]=nodes
            node_container[(roll1,roll2)][0]=g.root.children[tmpFirst].children[tmpSecond]
            tmpSecond+=1
        tmpFirst+=1
            #put in the root node for each roll pair, as the first node after indexing twice

    empty_string = ""
    for hist in range(0, 1<< max_hist): #add children for all nonterminal nodes
        if (hist % 100 == 0):
            print(hist)
        if node_container[(1,1)][hist] is None:
            continue
        if (hist % 100 == 0):
            print(hist)
        p1_to_play = pop_even(hist)
        player = "P1" if p1_to_play else "P2"
        curr_state = hist.bit_length() #1000 -> 4
        #first, in the 0,0 case, make a move
        if (player == "P1"):
            for roll1 in priv_p1:
                if curr_state > 6:
                    g.append_move(node_container[(roll1,1)][hist], player, [empty_string])
                    node_container[(roll1,1)][hist + (1<<max_hist)]=node_container[(roll1,1)][hist].children[0]
                    replicas = [node_container[(roll1, roll2)][hist] for roll2 in priv_p2 if roll2 != 1]
                    g.append_infoset(replicas, node_container[(roll1, 1)][hist].infoset)   
                    for roll2 in priv_p2:
                        if (roll2 == 1):
                            continue
                        node_container[(roll1, roll2)][hist + (1<<max_hist)] = node_container[(roll1, 1)][hist].children[0]
                else:
                    g.append_move(node_container[(roll1,1)][hist], player, [empty_string for j in range(curr_state, max_hist+1)])
                    for new_move in range(curr_state, max_hist + 1):

                        node_container[(roll1, 1)][hist + (1<<new_move)]=node_container[(roll1, 1)][hist].children[new_move - curr_state]
                    replicas = [node_container[(roll1, roll2)][hist] for roll2 in priv_p2 if roll2 != 1]
                    g.append_infoset(replicas, node_container[(roll1, 1)][hist].infoset)
                    for roll2 in priv_p2:
                        if (roll2 == 1):
                            continue
                        for new_move in range(curr_state, max_hist + 1):
                            node_container[(roll1, roll2)][hist + (1<<new_move)] = node_container[(roll1, roll2)][hist].children[new_move - curr_state]

        else:#player == "P2"
            for roll2 in priv_p2:
                if curr_state > 6:
                    g.append_move(node_container[(1,roll2)][hist], player, [empty_string])
                    node_container[(1,roll2)][hist + (1<<max_hist)]=node_container[(1,roll2)][hist].children[0]
                    replicas = [node_container[(roll1, roll2)][hist] for roll1 in priv_p1 if roll1 != 1]
                    g.append_infoset(replicas, node_container[(1, roll2)][hist].infoset)   
                    for roll1 in priv_p1:
                        if (roll1 == 1):
                            continue
                        node_container[(roll1, roll2)][hist + (1<<max_hist)] = node_container[(1, roll2)][hist].children[0]
                else:
                    g.append_move(node_container[(1,roll2)][hist], player, [empty_string for j in range(curr_state, max_hist+1)])
                    for new_move in range(curr_state, max_hist + 1):

                        node_container[(1, roll2)][hist + (1<<new_move)]=node_container[(1, roll2)][hist].children[new_move - curr_state]
                    replicas = [node_container[(roll1, roll2)][hist] for roll1 in priv_p1 if roll1 != 1]
                    g.append_infoset(replicas, node_container[(1, roll2)][hist].infoset)
                    for roll1 in priv_p1:
                        if (roll1 == 1):
                            continue
                        for new_move in range(curr_state, max_hist + 1):
                            node_container[(roll1, roll2)][hist + (1<<new_move)] = node_container[(roll1, roll2)][hist].children[new_move - curr_state]
   #set outcome for each terminal node
    p1_wins = g.add_outcome([1,-1], label="p1_wins")
    p2_wins = g.add_outcome([-1,1], label="p2_wins")
    for roll1 in priv_p1:
        for roll2 in priv_p2:
            for hist in range(1<<max_hist, 1<<(max_hist+1)): #+1: we mess with EMPTY BS
                if node_container[(roll1, roll2)][hist] is None:
                    continue
                node = node_container[(roll1, roll2)][hist]
                if (hist == 1<<max_hist):
                    g.set_outcome(node, p2_wins)    
                    continue
                if eval(roll1, roll2, hist)==1:
                    g.set_outcome(node, p1_wins)
                else:
                    g.set_outcome(node, p2_wins)
    result = gbt.nash.lp_solve(g)
    assert (result is not None)

    eq=result.equilibria[0]

    import csv
    with open("equilibrium.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
    # header
        writer.writerow(["player", "strategy_label", "probability"])

    # for each player in the original game
        for pl in g.players:
            mix_probs = eq[pl]         # e.g. a list of Rational or float
            for strat_obj, prob in zip(pl.strategies, mix_probs):
                # strat_obj.label is the string you passed when building the game
                writer.writerow([pl.label, strat_obj.label, float(prob)])

    print("Wrote equilibrium.csv")



    
    




    
    
 

if __name__ == "__main__":
    main()