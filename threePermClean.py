#https://helda.helsinki.fi/server/api/core/bitstreams/3f1f286b-3def-49e9-98ba-f887b1bc250e/content
# based on the two-trees paper by Scott Pakin (page 51/52)

from pysat.solvers import Solver
import numpy

rng = numpy.random.default_rng(seed = 22)
 
n=157
names = ['cd19',]
repeats = 40


base_scores = {}
for i in range(repeats):
    base_scores[i] = {}
    for name in names:
       base_scores[i][name] = 0
for repeat_id in range(repeats):
    local_n = n
    binary_assignment = rng.choice([0,1],local_n+1)
    #really get 3 permutations, but due to symmetry, first one can be identity.
    #and second one only cares about cycle length, through re-labeling
    #could even just specify cycle lengths (all >=2, sum to n)
    perm1 = list(range(local_n))
    perm2 = list(range(1,local_n))+[0] #cycle structure of one block "n"

#ALTERNATE perm3 choice: Go forward a bit every time
    perm3 = list(rng.permutation(range(local_n)))
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
        perm3[i] = int(rng.choice(tail[start_point:start_point+4*logn]))%local_n
        to_choose.remove(perm3[i])

    clauses3 = []
    for i in range(local_n):
        sign_perms3 = binary_assignment[perm1[i]+1]+binary_assignment[perm2[i]+1]+binary_assignment[perm3[i]+1]
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

#write the clauses to file
    #with open("ThreeRegPerm_v"+str(local_n)+"_c"+str(len(clauses3))+"_i"+str(repeat_id)+".cnf",'w') as out:
    #    out.write("c ThreeRegPerm_v"+str(local_n)+"_c"+str(len(clauses3))+"_i"+str(repeat_id)+".cnf\n")
    #    out.write("c "+str(binary_assignment)+"\n")
    #    out.write("p cnf "+str(local_n)+" "+str(len(clauses3))+"\n")
    #    for clause in clauses3:
    #         out.write(str(clause[0])+" "+str(clause[1])+" "+str(clause[2])+" 0\n")

    for name in names:
        results_dict = {}

        s = Solver(name = name, bootstrap_with = clauses3, use_timer=True)
        s.solve()
        temp_results = dict(s.accum_stats())
        for key in temp_results.keys():
            results_dict[key+str(3)] = temp_results[key]
        results_dict['time3'] = s.time_accum()

        results_dict['perm3'] = [int(j) for j in perm3]
        base_scores[repeat_id][name] = results_dict
times = []
for repeat_id in range(repeats):
    local_dict = base_scores[repeat_id][names[0]]
    times.append(local_dict['time3'])
print(times)