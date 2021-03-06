import sys
sys.path.insert(0,'/home/ubuntu/src') # Use appropriate address for source codes
sys.path.insert(1,'/home/ubuntu/examples/load_distribution') # Use appropriate address for instance codes
import time
from random import randint
from collections import OrderedDict
import json
from max_synthesis import *
from aut_translator import *
from instance_load_dist import *

def dict_to_str(in_dict):
    """Convert keys and values of a dictionary into string"""

    out_dict = {str(k): str(v) for k, v in in_dict.iteritems()}

    return out_dict

if __name__ == '__main__':

    t_s_total = time.time()
    n_instance = 5
    n_gen_all = [3, 3, 3, 3, 3]
    n_load_all = [3, 4, 5, 6, 7]
    gen_cap_all = [1, 1, 2, 2, 3]
    crit_load_all = [[0,0,0],
                     [0,0,0,0],
                     [0,1,0,0,0],
                     [0,1,0,0,0,0],
                     [0,1,1,0,0,0,0]]
    tr_bound_all = [8, 8, 8, 8, 8]

    print "--------------------------------------------------" + "\n"
    for inst in range(n_instance):

        print "Current instance: " + str(inst+1) + "\n"
        t_s_g = time.time()

        n_gen = n_gen_all[inst]
        n_load = n_load_all[inst]
        connection = [[1]*n_load]*n_gen
        gen_cap = [gen_cap_all[inst]]*n_gen
        load_demand = [1]*n_load
        max_fault = 1
        classes = None
        class_weights = None
        init_load = [1]
        init_len = [2]
        crit_load = crit_load_all[inst]
        flag_switch = False

        instance_name = 'ld'+'_'+str(n_gen)+'_'+str(n_load)+'_'+\
                        ''.join(map(str, init_load))+'_'+\
                        ''.join(map(str, crit_load))+'_'+\
                        'r'+str(randint(1, 1000))

        # hard and soft specifications
        (f_hard, f_soft) = instantiate_load_dist(n_gen, n_load, connection, 
                                                 gen_cap, load_demand, max_fault,
                                                 classes, class_weights, 
                                                 init_load, init_len, crit_load,
                                                 flag_switch)
        n_soft_const = len(f_soft)
        print "Specifications generated." + "\n"

        # definition of input and output variables
        (i_vars, o_vars) = generate_var_name(n_gen, n_load, connection, max_fault)
        (n_i, n_o) = (len(i_vars), len(o_vars))
        (I, I_val, io_vars) = gen_ex_io(n_i, i_vars, n_o, o_vars)

        # generate raw automata for secifications
        (B_h_UCBA, B_s_UBA, B_s_UCBA, (dt_aut_gen_hard, dt_aut_gen_soft)) = instance_gen(f_hard, f_soft)
        print "Automata generated." + "\n"

        # define implementation parameters
        tr_bound = tr_bound_all[inst] # bound on the size of the transition systems
        n_t = 1 # size of transition system; starts from 2
        f_real = False # flag for "sufficient" realizability
        f_max = False # flag for passing the bound on implementation

        # define statistics variables
        n_var_all = {} # number of variables
        n_clause_all = {} # number of clauses
        sol_value_all = {} # sum of weights of satisfied soft clauses
        dt_enc_all = {} # times for encoding under each bound on implementation size
        dt_maxsat_all = {} # times for solving maxsat instance under each bound on implementation size
        dt_onebound_all = {} # times for synthesis under each bound on implementation size

        while (f_real == False) and (f_max == False):
            
            t_s_onebound = time.time()
            # n_t = n*2 # exponential search
            n_t += 1 # linear search
            print "Current implementation bound: " + str(n_t) + "\n"
            if n_t == tr_bound:
                f_max = True
                print "Maximum implementation bound reached!" + "\n"

            T = ['t'+str(j) for j in range(n_t)]
            T0 = 't0'

            (A_h_UCBA, A_s_UBA, A_s_UCBA) = automata_for_encoding(B_h_UCBA, B_s_UBA, B_s_UCBA,
                                                                  I, I_val, io_vars, T)
            (f_real, (n_var, n_clause), sol_value, dt_enc, dt_maxsat) = synthesize(i_vars, I, I_val, o_vars, 
                                                                        A_h_UCBA, A_s_UBA, A_s_UCBA, T, n_t, n_o)
            f_real = False
            t_f_onebound = time.time()
            n_var_all[n_t] = n_var
            n_clause_all[n_t] = n_clause
            sol_value_all[n_t] = sol_value
            dt_enc_all[n_t] = dt_enc
            dt_maxsat_all[n_t] = dt_maxsat
            dt_onebound_all[n_t] = t_f_onebound - t_s_onebound

        t_f_g = time.time()
        dt_g = t_f_g - t_s_g # time for solving the whole instance

        # write results and statistics into a log file
        data_log = OrderedDict()

        data_log['Parameters'] = OrderedDict([
                                    ('n_gen', str(n_gen)),
                                    ('n_load', str(n_load)),
                                    ('connection', 'full'),
                                    ('gen_cap', str(gen_cap[0])),
                                    ('load_demand', '1'),
                                    ('max_fault', str(max_fault)),
                                    ('classes', 'none'),
                                    ('class_weights', 'none'),
                                    ('init_load', '-'.join(map(str, init_load))),
                                    ('init_len', str(init_len)),
                                    ('crit_load', '-'.join(map(str, crit_load))),
                                    ('flag_switch', str(flag_switch)),
                                    ('n_soft_const', str(n_soft_const))
                                    ])

        data_log['Time Statistics'] = OrderedDict([
                                          ('ltl to automaton', 
                                               OrderedDict([
                                                   ('hard constraints', str(dt_aut_gen_hard)),
                                                   ('soft constraints', str(dt_aut_gen_soft)),
                                                   ])),
                                          ('individual bound',
                                               OrderedDict([
                                                   ('encoding', OrderedDict(sorted(dict_to_str(dt_enc_all).iteritems()))),
                                                   ('Max-SAT solver', OrderedDict(sorted(dict_to_str(dt_maxsat_all).iteritems()))),
                                                   ('total', OrderedDict(sorted(dict_to_str(dt_onebound_all).iteritems())))
                                                   ])),
                                          ('total time', str(dt_g))
                                          ])

        data_log['Encoding Statistics'] = OrderedDict([
                                              ('implementation bound', str(tr_bound)),
                                              ('n_var', OrderedDict(sorted(dict_to_str(n_var_all).iteritems()))),
                                              ('n_clause', OrderedDict(sorted(dict_to_str(n_clause_all).iteritems())))
                                              ])

        data_log['Solution'] = OrderedDict([
                                   ('upperbound on value', str(n_soft_const * int(math.pow(n_soft_const,2)+n_soft_const+1))),
                                   ('solution value', OrderedDict(sorted(dict_to_str(sol_value_all).iteritems())))
                                   ])

        log_path = "/home/ubuntu/examples/load_distribution/log_files/"+instance_name
        with open(log_path, 'w') as f_log:
            json.dump(data_log, f_log, indent=4)

        print "Instance " + str(inst+1) + " completed." + "\n"
        print "--------------------------------------------------" + "\n"

    t_f_total = time.time()
    dt_total = t_f_total - t_s_total # time for solving all instances
    print "Total time: " + str(dt_total) + "\n"
