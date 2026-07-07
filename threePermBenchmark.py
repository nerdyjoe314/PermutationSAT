#https://helda.helsinki.fi/server/api/core/bitstreams/3f1f286b-3def-49e9-98ba-f887b1bc250e/content
# based on the two-trees paper by Scott Pakin (page 51/52)

from pysat.solvers import Solver
import time
import itertools
import matplotlib.pyplot as plt

import numpy
rng = numpy.random.default_rng(seed = 22)
 
n=177
names = ['cd19',]
repeats = 40


base_scores = {}
for i in range(repeats):
    base_scores[i] = {}
    for name in names:
       base_scores[i][name] = 0
for repeat_id in range(repeats):
    local_n = n
    tree_one = list(range(1,local_n+1))
    tree_two = list(range(1,local_n+1))
    rng.shuffle(tree_one)
    rng.shuffle(tree_two)
    #really get 3 permutations, but due to symmetry, first one can be identity.
    #and second one only cares about cycle length, through re-labeling
    #could even just specify cycle lengths (all >=2, sum to n)
    perm1 = list(range(local_n))
    perm2 = list(range(1,local_n))+[0] #cycle structure of one block "n"
#fully random
    perm3 = list(rng.permutation(range(local_n)))

#ALTERNATE perm3 choice: Go forward a bit every time
    perm4 = list(rng.permutation(range(local_n)))
    to_choose = list(range(local_n))
    logn = int(numpy.log(local_n))
    for i in range(local_n):
        tripled_choice = list(to_choose) +[i+local_n for i in to_choose] +[i+2*local_n for i in to_choose]
        start_point = 0
        while tripled_choice[start_point]<i+logn:
            start_point += 1
        tail = [i+2*local_n for i in to_choose]
        scale = 2
        while len(tail)<4*logn:
            scale += 1
            tail = tail+[i+scale*local_n for i in to_choose]
        tail = list(to_choose) +[i+local_n for i in to_choose] + tail
        perm4[i] = int(rng.choice(tail[start_point:start_point+4*logn]))%local_n
        to_choose.remove(perm4[i])


#ALTERNATE perm3 choices: window blocking
    #width = int(numpy.sqrt(n/2)) # first approximation. overlaps, blah blah blah
    width = 7# for this testing, going up to 9 introduces high failure rate. There's not currently good protections.
    while True:
        perm5= [-1]*n
        valid_positions = numpy.ones((n,n), dtype = int)
        for i in range(n):
            for m in range(-width-1,width+1):
                valid_positions[i,(i+m+n)%n] = 0
        while numpy.sum(valid_positions)>0:
            next_pos_index = rng.integers(numpy.where(valid_positions>0)[0].shape)
            selection = int(numpy.where(valid_positions>0)[0][next_pos_index][0]),int(numpy.where(valid_positions>0)[1][next_pos_index][0])
            perm5[selection[0]] = selection[1]
            # block out these from being chosen again
            valid_positions[selection[0],:] = 0
            valid_positions[:,selection[1]] = 0
            for i in range(1,width):
                for j in range(1,width-i):
                    valid_positions[(selection[0]+i+n)%n,(selection[1]+j+n)%n] = 0
                    valid_positions[(selection[0]+i+n)%n,(selection[1]-j+n)%n] = 0
                    valid_positions[(selection[0]-i+n)%n,(selection[1]-j+n)%n] = 0
                    valid_positions[(selection[0]-i+n)%n,(selection[1]+j+n)%n] = 0
        if min(perm5) == 0:
            break
    
    

    binary_assignment = rng.choice([0,1],local_n+1)

    clausesTwoTrees = []
    clauses3 = []
    clauses4 = []
    clauses5 = []
    for i in range(local_n):
        sign_perms3 = binary_assignment[perm1[i]+1]+binary_assignment[perm2[i]+1]+binary_assignment[perm3[i]+1]
        sign_perms4 = binary_assignment[perm1[i]+1]+binary_assignment[perm2[i]+1]+binary_assignment[perm4[i]+1]
        sign_perms5 = binary_assignment[perm1[i]+1]+binary_assignment[perm2[i]+1]+binary_assignment[perm5[i]+1]
        if sign_perms3 % 2 == 1:
            clauses3.append([-int(perm1[i]+1),-int(perm2[i]+1),-int(perm3[i]+1)])
            clauses3.append([-int(perm1[i]+1),int(perm2[i]+1),int(perm3[i]+1)])
            clauses3.append([int(perm1[i]+1),-int(perm2[i]+1),int(perm3[i]+1)])
            clauses3.append([int(perm1[i]+1),int(perm2[i]+1),-int(perm3[i]+1)])
        else:
            clauses3.append([int(perm1[i]+1),int(perm2[i]+1),int(perm3[i]+1)])
            clauses3.append([-int(perm1[i]+1),-int(perm2[i]+1),int(perm3[i]+1)])
            clauses3.append([int(perm1[i]+1),-int(perm2[i]+1),-int(perm3[i]+1)])
            clauses3.append([-int(perm1[i]+1),int(perm2[i]+1),-int(perm3[i]+1)])
        if sign_perms4 % 2 == 1:
            clauses4.append([-int(perm1[i]+1),-int(perm2[i]+1),-int(perm4[i]+1)])
            clauses4.append([-int(perm1[i]+1),int(perm2[i]+1),int(perm4[i]+1)])
            clauses4.append([int(perm1[i]+1),-int(perm2[i]+1),int(perm4[i]+1)])
            clauses4.append([int(perm1[i]+1),int(perm2[i]+1),-int(perm4[i]+1)])
        else:
            clauses4.append([int(perm1[i]+1),int(perm2[i]+1),int(perm4[i]+1)])
            clauses4.append([-int(perm1[i]+1),-int(perm2[i]+1),int(perm4[i]+1)])
            clauses4.append([int(perm1[i]+1),-int(perm2[i]+1),-int(perm4[i]+1)])
            clauses4.append([-int(perm1[i]+1),int(perm2[i]+1),-int(perm4[i]+1)])
        if sign_perms5 % 2 == 1:
            clauses5.append([-int(perm1[i]+1),-int(perm2[i]+1),-int(perm5[i]+1)])
            clauses5.append([-int(perm1[i]+1),int(perm2[i]+1),int(perm5[i]+1)])
            clauses5.append([int(perm1[i]+1),-int(perm2[i]+1),int(perm5[i]+1)])
            clauses5.append([int(perm1[i]+1),int(perm2[i]+1),-int(perm5[i]+1)])
        else:
            clauses5.append([int(perm1[i]+1),int(perm2[i]+1),int(perm5[i]+1)])
            clauses5.append([-int(perm1[i]+1),-int(perm2[i]+1),int(perm5[i]+1)])
            clauses5.append([int(perm1[i]+1),-int(perm2[i]+1),-int(perm5[i]+1)])
            clauses5.append([-int(perm1[i]+1),int(perm2[i]+1),-int(perm5[i]+1)])
    for i in range(int(local_n/2)):
        sign_tree_one = binary_assignment[tree_one[i]]+binary_assignment[tree_one[2*i+1]]+binary_assignment[tree_one[2*i+2]]
        sign_tree_two = binary_assignment[tree_two[i]]+binary_assignment[tree_two[2*i+1]]+binary_assignment[tree_two[2*i+2]]
        if sign_tree_one % 2 == 1:
            clausesTwoTrees.append([tree_one[i],tree_one[2*i+1],tree_one[2*i+2]])
            clausesTwoTrees.append([-tree_one[i],-tree_one[2*i+1],tree_one[2*i+2]])
            clausesTwoTrees.append([tree_one[i],-tree_one[2*i+1],-tree_one[2*i+2]])
            clausesTwoTrees.append([-tree_one[i],tree_one[2*i+1],-tree_one[2*i+2]])
        else:
            clausesTwoTrees.append([-tree_one[i],-tree_one[2*i+1],-tree_one[2*i+2]])
            clausesTwoTrees.append([-tree_one[i],tree_one[2*i+1],tree_one[2*i+2]])
            clausesTwoTrees.append([tree_one[i],-tree_one[2*i+1],tree_one[2*i+2]])
            clausesTwoTrees.append([tree_one[i],tree_one[2*i+1],-tree_one[2*i+2]])
        if sign_tree_two % 2 == 1:
            clausesTwoTrees.append([tree_two[i],tree_two[2*i+1],tree_two[2*i+2]])
            clausesTwoTrees.append([-tree_two[i],-tree_two[2*i+1],tree_two[2*i+2]])
            clausesTwoTrees.append([tree_two[i],-tree_two[2*i+1],-tree_two[2*i+2]])
            clausesTwoTrees.append([-tree_two[i],tree_two[2*i+1],-tree_two[2*i+2]])
        else:
            clausesTwoTrees.append([-tree_two[i],-tree_two[2*i+1],-tree_two[2*i+2]])
            clausesTwoTrees.append([-tree_two[i],tree_two[2*i+1],tree_two[2*i+2]])
            clausesTwoTrees.append([tree_two[i],-tree_two[2*i+1],tree_two[2*i+2]])
            clausesTwoTrees.append([tree_two[i],tree_two[2*i+1],-tree_two[2*i+2]])

#write the clauses to file
    #with open("ThreeRegPerm_v"+str(local_n)+"_c"+str(len(clauses3))+"_i"+str(repeat_id)+".cnf",'w') as out:
    #    out.write("c ThreeRegPerm_v"+str(local_n)+"_c"+str(len(clauses3))+"_i"+str(repeat_id)+".cnf\n")
    #    out.write("c "+str(binary_assignment)+"\n")
    #    out.write("p cnf "+str(local_n)+" "+str(len(clauses3))+"\n")
    #    for clause in clauses3:
    #         out.write(str(clause[0])+" "+str(clause[1])+" "+str(clause[2])+" 0\n")

    for name in names:
        results_dict = {}
        s = Solver(name = name, bootstrap_with = clausesTwoTrees, use_timer=True)
        s.solve()
        temp_results = dict(s.accum_stats())
        for key in temp_results.keys():
            results_dict[key+str(2)] = temp_results[key]
        results_dict['time2'] = s.time_accum()

        s = Solver(name = name, bootstrap_with = clauses3, use_timer=True)
        s.solve()
        temp_results = dict(s.accum_stats())
        for key in temp_results.keys():
            results_dict[key+str(3)] = temp_results[key]
        results_dict['time3'] = s.time_accum()

        s = Solver(name = name, bootstrap_with = clauses4, use_timer=True)
        s.solve()
        temp_results = dict(s.accum_stats())
        for key in temp_results.keys():
            results_dict[key+str(4)] = temp_results[key]
        results_dict['time4'] = s.time_accum()

        s = Solver(name = name, bootstrap_with = clauses5, use_timer=True)
        s.solve()
        temp_results = dict(s.accum_stats())
        for key in temp_results.keys():
            results_dict[key+str(5)] = temp_results[key]
        results_dict['time5'] = s.time_accum()

        results_dict['perm3'] = [int(j) for j in perm3]
        results_dict['perm4'] = [int(j) for j in perm4]
        results_dict['perm5'] = [int(j) for j in perm5]
        base_scores[repeat_id][name] = results_dict
window = []
lognth = []
randot = []
twotre = []
for repeat_id in range(repeats):
    local_dict = base_scores[repeat_id][names[0]]
    window.append(local_dict['time5'])
    lognth.append(local_dict['time4'])
    randot.append(local_dict['time3'])
    twotre.append(local_dict['time2'])
fig, axs = plt.subplots(1,4,sharey=True)
fig.suptitle('Time to solve n='+str(n)+" in seconds, sample size 40")
axs[0].boxplot(window)
axs[0].set_title('Blocking')
axs[1].boxplot(lognth)
axs[1].set_title('Log Ahead')
axs[2].boxplot(randot)
axs[2].set_title('Random')
axs[3].boxplot(twotre)
axs[3].set_title('Two Trees')
plt.show()