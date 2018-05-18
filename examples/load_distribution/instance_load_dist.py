import numpy as np
import itertools

def instantiate_load_dist(n_gen, n_load, connection, gen_cap, load_demand, max_fault,
                          classes, class_weights, init_load, init_len, crit_load, flag_switch):
    """Create an instance of load distribution problem with given paremeters

    Arguments:
    n_gen -- number of power supplies
    n_load -- number of loads
    connection -- connection of a power supply to loads
    gen_cap -- capacity of power supplies
    load_demand -- power demand of each load
    max_fault -- bound on number of faulty power supplies
    classes -- correlation amongst the loads
    class_weights -- weights associated with classes
    init_load -- list of loads that should be powered initially
    init_len -- list of length of initial period for initial loads
    crit_load -- a boolean list determining the critical loads
    flag_switch -- whether to include soft constraints from switching limit
    """

    # availability lists using positive integers
    ap = [[i+1 for i, ii in enumerate(connection[j]) if ii!=0] for j in range(n_gen)]
    al = [[j+1 for j, jj in enumerate([connection[jjj][i] for jjj in range(n_gen)])
           if jj!=0] for i in range(n_load)]
    e_bits = int(np.ceil(np.log2(n_gen+1)))

    # hard constraints: "critical loads should be powered"
    flag_crit_load_power = (sum(crit_load) > 0)
    crit_load_power = 'G ( ' +\
                      ' & '.join(['( '+' | '.join(['s_'+str(i+1)+'_'+str(j) for j in al[i]])+' )'
                      for i in range(n_load) if crit_load[i]]) +\
                      ' )'

    # soft constraints: "non-critical loads should be powered"
    flag_noncrit_load_power = (sum(crit_load) < n_load)
    noncrit_load_power = [' | '.join(['s_'+str(i+1)+'_'+str(j) for j in al[i]])
                          for i in range(n_load) if not crit_load[i] if (i+1) not in init_load]

    # hard constraints: "an intitial load must be active during initial time steps"
    flag_init_power = (len(init_load) > 0)
    init_power = '( ' +\
                 ' & '.join([' & '.join(['X'*t+'('+' | '.join([
                 's_'+str(l_init)+'_'+str(i)
                 for i in al[l_init-1]])+')'
                 for t in range(t_init)])
                 for l_init, t_init in zip(init_load, init_len)]) +\
                 ' )'

    # hard constraints: "a load should only be assigned to one power supply"
    flag_unique_supp = True
    unique_supp = 'G ( ' +\
                  ' & '.join([' & '.join(['(' +
                  's_'+str(i+1)+'_'+str(j1) +
                  ' -> ' +
                  ' & '.join(['~'+'s_'+str(i+1)+'_'+str(j2)
                  for j2 in al[i] if j2!=j1]+['true']) + ')'
                  for j1 in al[i]]) for i in range(n_load) if al[i]]) +\
                  ' )'

    # hard constraints: "do not exceed capacity of power supplies"
    flag_limit_cap = True
    limit_cap = 'G ( ' +\
                ' & '.join([' & '.join(['(' +
                ' & '.join(['s_'+str(i1)+'_'+str(j+1) for i1 in l_sub]) +
                ' -> ' +
                ' & '.join(['~'+'s_'+str(i2)+'_'+str(j+1) for i2 in list(set(ap[j])-set(l_sub))]+['true']) + ')'
                for l_sub in itertools.combinations(ap[j], gen_cap[j])])
                for j in range(n_gen)]) +\
                ' )'

    # hard constraints: "disconnect loads from faulty power supplies"
    flag_fault_disc = (max_fault > 0)
    fault_disc = 'G ( ' +\
                 ' & '.join([' & '.join(['(' +
                 ' & '.join(['e_'+str(k)+'_'+str(b) if format(j+1, "0"+str(e_bits)+"b")[b]=='1'
                             else '~'+'e_'+str(k)+'_'+str(b) for b in range(e_bits)]) +
                 ' -> ' +
                 ' & '.join(['~'+'s_'+str(i)+'_'+str(j+1) for i in ap[j]]+['true']) + ')'
                 for j in range(n_gen)])
                 for k in range(1,max_fault+1)]) +\
                 ' )'

    # soft constraints: "lower-bound the intervals between switching a load by 2"
    flag_switch_bound = flag_switch
    switch_bound = ['s_'+str(i+1)+'_'+str(j)+' & '+\
                    '~'+'('+' | '.join(['('+' & '.join(['e_'+str(k)+'_'+str(b)
                                                        if format(j, "0"+str(e_bits)+"b")[b]=='1'
                                                        else '~'+'e_'+str(k)+'_'+str(b) for b in range(e_bits)])+')'
                                        for k in range(1,max_fault+1)])+')'+\
                    ' -> '+'X'+'s_'+str(i+1)+'_'+str(j)
                    for i in range(n_load) for j in al[i]]

    # select desired specifications
    f_hard_list = [(crit_load_power, flag_crit_load_power),
                   (init_power, flag_init_power),
                   (unique_supp, flag_unique_supp),
                   (limit_cap, flag_limit_cap),
                   (fault_disc, flag_fault_disc)]
    f_hard = [' & '.join([spec[0] for spec in f_hard_list if spec[1]])]

    f_soft_list = [(noncrit_load_power, flag_noncrit_load_power),
                   (switch_bound, flag_switch_bound)]
    f_soft = [subspec for spec in f_soft_list for subspec in spec[0] if spec[1]]

    return (f_hard, f_soft)

def generate_var_name(n_gen, n_load, connection, max_fault):
    """Generate the names of input and output variables for an instance of load distribution"""

    al = [[j+1 for j, jj in enumerate([connection[jjj][i] for jjj in range(n_gen)])
           if jj!=0] for i in range(n_load)]
    e_bits = int(np.ceil(np.log2(n_gen+1)))

    i_vars = ['e_'+str(f)+'_'+str(b) for f in range(1, max_fault+1) for b in range(e_bits)]
    o_vars = ['s_'+str(i+1)+'_'+str(j) for i in range(n_load) for j in al[i]]

    return (i_vars, o_vars)