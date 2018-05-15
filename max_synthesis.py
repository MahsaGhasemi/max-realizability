import numpy as np
import sympy as sp
import math
import subprocess
import time

windows = False # OS determiner
maxsat_solver = "/home/ubuntu/open-wbo2.0/open-wbo_static" # Use appropriate address for MaxSAT solver

def negate(clauses):
    """Negate a (set) of clauses without double negation"""

    n_c = len(clauses)
    new_clauses = []
    for i in range(n_c):
        if clauses[i][0] == 'T':
            new_clauses.append('F')
        elif clauses[i][0] == 'F':
            new_clauses.append('T')
        elif clauses[i][0] == '~':
            new_clauses.append(clauses[i][1:])
        else:
            new_clauses.append('~' + clauses[i])
    return new_clauses

def bitvector_ineq(ineq_type, n_bit, l_type, t_names, q_names, index):
    """Generate bit vector inequalities in CNF"""

    clause_ineq = [] # clauses associated with inequality encoding

    # ">" inequality
    if ineq_type == 'g':
        temp_g = ['g'+l_type+str(b)+'_'+t_names[1]+'_'+q_names[1]+'_'+\
                  t_names[0]+'_'+q_names[0]+'_'+str(index) \
                  for b in range(n_bit)] # vector of g's
        clause_ineq.append(' '.join(temp_g))
        for b in range(n_bit):
            temp_l1 = 'l'+l_type+str(b)+'_'+t_names[1]+'_'+q_names[1]+\
                      '_'+str(index)
            temp_l2 = 'l'+l_type+str(b)+'_'+t_names[0]+'_'+q_names[0]+\
                      '_'+str(index)
            clause_ineq.append(temp_l1+' ~'+temp_l2+' '+' '.join(temp_g[:b]))

    # ">=" inequality
    elif ineq_type == 'ge':
        temp_g = ['g'+l_type+str(b)+'_'+t_names[1]+'_'+q_names[1]+'_'+\
                  t_names[0]+'_'+q_names[0]+'_'+str(index) \
                  for b in range(n_bit)] # vector of g's
        for b in range(n_bit):
            temp_l1 = 'l'+l_type+str(b)+'_'+t_names[1]+'_'+q_names[1]+\
                      '_'+str(index)
            temp_l2 = 'l'+l_type+str(b)+'_'+t_names[0]+'_'+q_names[0]+\
                      '_'+str(index)
            clause_ineq.append(temp_l1+' ~'+temp_l2+' '+' '.join(temp_g[:b]))

    # "=" equality to zero --> resetting the counter at next state
    elif ineq_type == 'e':
        clause_ineq.extend('~l'+l_type+str(b)+'_'+t_names[1]+
                           '_'+q_names[1]+'_'+str(index)
                           for b in range(n_bit))

    # unrecognized inequality
    else:
        raise Exception("Unexpected inequality type")

    return clause_ineq

def inc_annot(q, t, i, qq, tt, l_type, n, n_bit, Q_rej, Q_rej_UCBA):
    """Generate the clauses related to increment of annotations"""

    if l_type == 's':
        cnf_temp = ['lb_'+tt+'_'+qq+'_'+str(n)]
        if qq in Q_rej:
            cnf_temp.extend(bitvector_ineq('g', n_bit, 's', [t,tt], [q,qq], n))
        else:
            cnf_temp.extend(bitvector_ineq('ge', n_bit, 's', [t,tt], [q,qq], n))

    elif l_type == 'i':
        cnf_temp = []
        cnf_temp.append('~lei_'+t+'_'+q+'_'+str(n)+' '+
                        'lei_'+tt+'_'+qq+'_'+str(n))

        cnf_temp.extend(
                        [
                         '~lei_'+t+'_'+q+'_'+str(n)+' '+
                         'bad_'+t+'_'+q+'_'+i+'_'+qq+'_'+str(n)+' '+
                         cl_temp
                         for cl_temp in bitvector_ineq('ge', n_bit, 'i',
                                                       [t,tt], [q,qq], n)])

        cnf_temp.extend(
                        [
                         '~lei_'+t+'_'+q+'_'+str(n)+' '+
                         '~bad_'+t+'_'+q+'_'+i+'_'+qq+'_'+str(n)+' '+
                         cl_temp
                         for cl_temp in bitvector_ineq('g', n_bit, 'i',
                                                       [t,tt], [q,qq], n)])

    elif l_type == 'f':
        cnf_temp = []
        cnf_temp.append('~lef_'+t+'_'+q+'_'+str(n)+' '+
                        'lef_'+tt+'_'+qq+'_'+str(n))

        if qq in Q_rej_UCBA[n-1]:
            cnf_temp.extend('~lef_'+t+'_'+q+'_'+str(n)+' '+
                            cl_temp
                            for cl_temp in bitvector_ineq('g', n_bit, 'f',
                                                          [t,tt], [q,qq], n))
        else:
            cnf_temp.extend('~lef_'+t+'_'+q+'_'+str(n)+' '+
                            cl_temp
                            for cl_temp in bitvector_ineq('ge', n_bit, 'f',
                                                          [t,tt], [q,qq], n))

    else:
        raise Exception("Unexpected annotation label")

    return cnf_temp

def soft_const(n_soft_const, n_bit_s, t0, Q0_UBA, Q0_UCBA):
    """Generate soft clauses for the soft constraints"""
    # t0 -> one string
    # Q0 -> list of list strings corresponding to initial states of soft automata
    # assumes only one initial state per automata

    hard_clauses = []
    soft_clauses = []
    soft_weights = []
    for n in range(1,n_soft_const+1):
        q0_UBA = Q0_UBA[n-1][0]
        q0_UCBA = Q0_UCBA[n-1][0]

        n_bit = n_bit_s[n-1][0]
        # encoding of s1 <-> li = m
        hard_clauses.append('s1_'+str(n)+' '+\
                            '~lei'+'_'+t0+'_'+q0_UBA+'_'+str(n)+' '+\
                            ' '.join(['~li'+str(b)+'_'+t0+'_'+q0_UBA+'_'+str(n)
                            for b in range(n_bit)]))
        hard_clauses.extend(['~s1_'+str(n)+' li'+str(b)+'_'+t0+'_'+q0_UBA+'_'+str(n)
                             for b in range(n_bit)])
        hard_clauses.append('~s1_'+str(n)+' '+\
                            'lei'+'_'+t0+'_'+q0_UBA+'_'+str(n))

        # encoding of s2 <-> li = 0
        hard_clauses.append('s2_'+str(n)+' '+\
                            '~lei'+'_'+t0+'_'+q0_UBA+'_'+str(n)+' '+\
                            ' '.join(['li'+str(b)+'_'+t0+'_'+q0_UBA+'_'+str(n)
                            for b in range(n_bit)]))
        hard_clauses.extend(['~s2_'+str(n)+' ~li'+str(b)+'_'+t0+'_'+q0_UBA+'_'+str(n)
                             for b in range(n_bit)])
        hard_clauses.append('~s2_'+str(n)+' '+\
                            'lei'+'_'+t0+'_'+q0_UBA+'_'+str(n))

        n_bit = n_bit_s[n-1][1]
        # encoding of s3 <-> lef_t0_q0_str(n)
        hard_clauses.append('s3_'+str(n)+' '+\
                            '~lef'+'_'+t0+'_'+q0_UCBA+'_'+str(n))
        hard_clauses.append('~s3_'+str(n)+' '+\
                            'lef'+'_'+t0+'_'+q0_UCBA+'_'+str(n))

        # prioritize level of satisfaction
        soft_clauses.extend([
                             's1_'+str(n),
                             's1_'+str(n)+' '+'s2_'+str(n),
                             's1_'+str(n)+' '+'s2_'+str(n)+' '+'s3_'+str(n)])

        # assign weights to soft clauses
        soft_weights.extend([1, n_soft_const, int(math.pow(n_soft_const,2))])

    return (hard_clauses, soft_clauses, soft_weights)

def gen_mapping(I, n_o, Q, Q_UBA, bad_tra_UBA, Q_UCBA, T, n_bit_h, n_bit_s):
    """Generates the mapping from original variables to integers"""

    n_soft_const = len(Q_UBA)

    # mapping binary variables to integers
    count_int = 1 # counter for integers
    mapping = {} # mapping dictionary

    # variables for transition relation
    for t in T:
        for i in I:
            for tt in T:

                mapping['tau_'+t+'_'+i+'_'+tt] = str(count_int)
                count_int += 1

    # variables for labeling function
    for t in T:
        for i in I:
            for o_ind in range(n_o):

                mapping['o'+str(o_ind)+'_'+t+'_'+i] = str(count_int)
                count_int += 1

    # variables for annotation
    n_bit = n_bit_h
    for t in T:
        for q in Q:

            # reachability variables
            mapping['lb_'+t+'_'+q+'_'+'1'] = str(count_int)
            count_int += 1

            # counting variables
            for b in range(n_bit):

                mapping['ls'+str(b)+'_'+t+'_'+q+'_'+'1'] = str(count_int)
                count_int += 1

                # representative variables
                for tt in T:
                    for qq in Q:
                        for l_type in ['s']:

                            mapping['g'+l_type+str(b)+'_'+t+'_'+q+'_'+\
                                    tt+'_'+qq+'_'+'1'] = str(count_int)
                            count_int += 1

    for n in range(1,n_soft_const+1):

        for t in T:
            # for UBA
            n_bit = n_bit_s[n-1][0]

            for q in Q_UBA[n-1]:

                # reachability variables
                mapping['lei_'+t+'_'+q+'_'+str(n)] = str(count_int)
                count_int += 1

                # counting variables
                for b in range(n_bit):
                    mapping['li'+str(b)+'_'+t+'_'+q+'_'+str(n)] = str(count_int)
                    count_int += 1

                    # representative variables
                    for tt in T:
                        for qq in Q_UBA[n-1]:
                            for l_type in ['i']:
                                mapping['g'+l_type+str(b)+'_'+t+'_'+q+'_'+\
                                        tt+'_'+qq+'_'+str(n)] = str(count_int)
                                count_int += 1

            # for UCBA
            n_bit = n_bit_s[n-1][1]
            for q in Q_UCBA[n-1]:

                # reachability variables
                mapping['lef_'+t+'_'+q+'_'+str(n)] = str(count_int)
                count_int += 1

                # counting variables
                for b in range(n_bit):
                    mapping['lf'+str(b)+'_'+t+'_'+q+'_'+str(n)] = str(count_int)
                    count_int += 1

                    # representative variables
                    for tt in T:
                        for qq in Q_UCBA[n-1]:
                            for l_type in ['f']:
                                mapping['g'+l_type+str(b)+'_'+t+'_'+q+'_'+\
                                        tt+'_'+qq+'_'+str(n)] = str(count_int)
                                count_int += 1

    for n in range(1,n_soft_const+1):
        mapping['s1_'+str(n)] = str(count_int)
        count_int += 1
        mapping['s2_'+str(n)] = str(count_int)
        count_int += 1
        mapping['s3_'+str(n)] = str(count_int)
        count_int += 1

    ## should generalize
    for n in range(n_soft_const):
        for k, v in bad_tra_UBA[n].iteritems():
            mapping[k] = str(count_int)
            count_int += 1

    count_int -= 1 # final value

    return (count_int, mapping)

def clause_check(cl_list):
    """Prepares a list of clause to bo included in the DIMACS file"""

    # format of cl_list = ['x y ...', 'x z w ...', ...]
    cl_list_mod = []
    for cl in cl_list:
        if 'T' in cl:
            continue
        else:
            cl_mod = cl
            cl_mod = cl_mod.replace(' F', '')
            cl_mod = cl_mod.replace('F ', '')
            cl_mod = cl_mod.replace('~', '-')
            cl_list_mod.append(cl_mod)

    return cl_list_mod

def clause_map(cl_list, mapping):
    """Prepares a list of clause to bo included in the DIMACS file"""

    # format of cl_list = ['x y ...', 'x z w ...', ...]
    cl_list_rep = []
    for cl in cl_list:
        lits = [l for l in cl.split(' ') if l] # literals, remove empty ones
        if lits:
            cl_list_rep.append(' '.join(['-'+mapping[l[1:]] if l[0] in ['~','-']
                                         else mapping[l] for l in lits]))

    return cl_list_rep

def encode_synthesis(i_vars, I, I_val, o_vars,
                     A_h_UCBA, A_s_UBA, A_s_UCBA, T, n_t, count_int, mapping,n_bit_h, n_bit_s):
    """Generates Max-SAT encoding of synthesis problem"""

    # extract automaton definition
    (Q,Q0,delta,Q_rej) = A_h_UCBA[0]
    #
    n_soft_const = len(A_s_UBA)
    #
    Q_UBA = [A_s_UBA[i][0] for i in range(n_soft_const)]
    Q0_UBA = [A_s_UBA[i][1] for i in range(n_soft_const)]
    delta_UBA = [A_s_UBA[i][2] for i in range(n_soft_const)]
    Q_rej_UBA = [A_s_UBA[i][3] for i in range(n_soft_const)]
    bad_tra_UBA = [A_s_UBA[i][4] for i in range(n_soft_const)]
    #
    Q_UCBA = [A_s_UCBA[i][0] for i in range(n_soft_const)]
    Q0_UCBA = [A_s_UCBA[i][1] for i in range(n_soft_const)]
    delta_UCBA = [A_s_UCBA[i][2] for i in range(n_soft_const)]
    Q_rej_UCBA = [A_s_UCBA[i][3] for i in range(n_soft_const)]

    # a negated version of delta in CNF
    delta_neg = {}
    for k, v in delta.items():
        delta_neg[k] = [' '.join(negate(vv.split(' '))) for vv in v]

    delta_UBA_neg = []
    for n in range(n_soft_const):
        temp_dict = {}
        for k, v in delta_UBA[n].items():
            temp_dict[k] = [' '.join(negate(vv.split(' '))) for vv in v]
        delta_UBA_neg.append(temp_dict)

    delta_UCBA_neg = []
    for n in range(n_soft_const):
        temp_dict = {}
        for k, v in delta_UCBA[n].items():
            temp_dict[k] = [' '.join(negate(vv.split(' '))) for vv in v]
        delta_UCBA_neg.append(temp_dict)

    ##### CNF-FORMAT FILE GENERATION #####
    with open('init_encoding.txt', 'w') as f, open('dimacs_encoding_nohead.txt', 'w') as f_dimacs:
        n_clause = 0 # counter for number of clauses
        max_weight = n_soft_const * int(math.pow(n_soft_const,2)+n_soft_const+1) + 1 # weight of hard clauses

        # initial state must be reachable
        for q0 in Q0:
            cl = 'lb_t0_'+q0+'_1'
            cl_dimacs = str(max_weight) + ' ' +\
                        mapping['lb_t0_'+q0+'_1'] + ' 0\n'

            n_clause += 1
            f.write(cl + '\n')
            f_dimacs.write(cl_dimacs)

        # there should be a successor per state and input in transition system
        for t in T:
            for i in I:
                cl = ' '.join(['tau_'+t+'_'+i+'_'+tt for tt in T])
                cl_dimacs = str(max_weight) + ' ' +\
                            ' '.join([mapping['tau_'+t+'_'+i+'_'+tt] for tt in T]) +\
                            ' 0\n'

                n_clause += 1
                f.write(cl + '\n')
                f_dimacs.write(cl_dimacs)

        # equivalency of representative variables for bit comparison
        # for hard specification
        l_type = 's'
        n_bit = n_bit_h
        for b in range(n_bit):
            for t in T:
                for q in Q:
                    for tt in T:
                        for qq in Q:
                            temp_l1 = 'l'+l_type+str(b)+'_'+t+'_'+q+'_'+'1'
                            temp_l2 = 'l'+l_type+str(b)+'_'+tt+'_'+qq+'_'+'1'
                            temp_g = 'g'+l_type+str(b)+'_'+t+'_'+q+'_'+\
                                     tt+'_'+qq+'_'+'1'
                            cl = '~'+temp_g+' '+temp_l1+'\n'+\
                                 '~'+temp_g+' '+'~'+temp_l2+'\n'+\
                                 '~'+temp_l1+' '+temp_l2+' '+temp_g+'\n'
                            cl_dimacs = str(max_weight) + ' ' +\
                                        '-'+mapping[temp_g]+' '+mapping[temp_l1]+' 0\n'+\
                                        str(max_weight) + ' ' +\
                                        '-'+mapping[temp_g]+' '+'-'+mapping[temp_l2]+' 0\n'+\
                                        str(max_weight) + ' ' +\
                                        '-'+mapping[temp_l1]+' '+mapping[temp_l2]+' '+mapping[temp_g]+' 0\n'

                            n_clause += 3
                            f.write(cl)
                            f_dimacs.write(cl_dimacs)

        # for soft specification and UBA
        l_type = 'i'
        for n in range(1,n_soft_const+1):
            n_bit = n_bit_s[n-1][0]

            for b in range(n_bit):
                for t in T:
                    for q in Q_UBA[n-1]:
                        for tt in T:
                            for qq in Q_UBA[n-1]:
                                temp_l1 = 'l'+l_type+str(b)+'_'+t+'_'+q+'_'+str(n)
                                temp_l2 = 'l'+l_type+str(b)+'_'+tt+'_'+qq+'_'+str(n)
                                temp_g = 'g'+l_type+str(b)+'_'+t+'_'+q+'_'+\
                                         tt+'_'+qq+'_'+str(n)
                                cl = '~'+temp_g+' '+temp_l1+'\n'+\
                                     '~'+temp_g+' '+'~'+temp_l2+'\n'+\
                                     '~'+temp_l1+' '+temp_l2+' '+temp_g+'\n'
                                cl_dimacs = str(max_weight) + ' ' +\
                                            '-'+mapping[temp_g]+' '+mapping[temp_l1]+' 0\n'+\
                                            str(max_weight) + ' ' +\
                                            '-'+mapping[temp_g]+' '+'-'+mapping[temp_l2]+' 0\n'+\
                                            str(max_weight) + ' ' +\
                                            '-'+mapping[temp_l1]+' '+mapping[temp_l2]+' '+mapping[temp_g]+' 0\n'

                                n_clause += 3
                                f.write(cl)
                                f_dimacs.write(cl_dimacs)

        # for soft specification and UCBA
        l_type = 'f'
        for n in range(1,n_soft_const+1):
            n_bit = n_bit_s[n-1][1]

            for b in range(n_bit):
                for t in T:
                    for q in Q_UCBA[n-1]:
                        for tt in T:
                            for qq in Q_UCBA[n-1]:
                                temp_l1 = 'l'+l_type+str(b)+'_'+t+'_'+q+'_'+str(n)
                                temp_l2 = 'l'+l_type+str(b)+'_'+tt+'_'+qq+'_'+str(n)
                                temp_g = 'g'+l_type+str(b)+'_'+t+'_'+q+'_'+\
                                         tt+'_'+qq+'_'+str(n)
                                cl = '~'+temp_g+' '+temp_l1+'\n'+\
                                     '~'+temp_g+' '+'~'+temp_l2+'\n'+\
                                     '~'+temp_l1+' '+temp_l2+' '+temp_g+'\n'
                                cl_dimacs = str(max_weight) + ' ' +\
                                            '-'+mapping[temp_g]+' '+mapping[temp_l1]+' 0\n'+\
                                            str(max_weight) + ' ' +\
                                            '-'+mapping[temp_g]+' '+'-'+mapping[temp_l2]+' 0\n'+\
                                            str(max_weight) + ' ' +\
                                            '-'+mapping[temp_l1]+' '+mapping[temp_l2]+' '+mapping[temp_g]+' 0\n'

                                n_clause += 3
                                f.write(cl)
                                f_dimacs.write(cl_dimacs)

        n_clause_mod = n_clause # #(clause) after removing redundancies

        # equivalency of representative variables for bad transitions
        for n in range(n_soft_const):
            for var, lab in bad_tra_UBA[n].iteritems():
                cl = '\n'.join(lab) + '\n'
                cl_dimacs = ' 0\n'.join([str(max_weight)+' '+c for c in
                                         clause_map(clause_check(lab), mapping)]) + ' 0\n'

                n_clause += len(lab)
                f.write(cl)
                n_clause_mod += len(clause_map(clause_check(lab), mapping))
                f_dimacs.write(cl_dimacs)

        # propagation of reachability and counting by annotations
        # for hard constraints
        n_bit = n_bit_h
        for q in Q:
            for t in T:
                for i in I:
                    for qq in Q:
                        # corresponding delta formula
                        del_temp = delta_neg['del_'+t+'_'+q+'_'+i+'_'+qq]

                        for tt in T:
                            # clauses resulting from annotation validity
                            cnf_temp = inc_annot(q,t,i,qq,tt,'s',1,n_bit,Q_rej,Q_rej_UCBA)

                            for d in del_temp:
                                for c in cnf_temp:
                                    cl = '~lb_' + t + '_' + q + '_1'\
                                         ' ' +\
                                         d +\
                                         ' ' +\
                                         '~' + 'tau_' + t + '_' + i + '_' + tt +\
                                         ' ' +\
                                         c
                                    cl_dimacs = clause_map(clause_check([cl]), mapping)

                                    n_clause += 1
                                    f.write(cl + '\n')
                                    if cl_dimacs:
                                        n_clause_mod += 1
                                        f_dimacs.write(str(max_weight) + ' ' + cl_dimacs[0] + ' 0\n')

        # for soft constraints and UBA
        for n in range(1,n_soft_const+1):
            n_bit = n_bit_s[n-1][0]

            for q in Q_UBA[n-1]:
                for t in T:
                    for i in I:
                        for qq in Q_UBA[n-1]:
                            # corresponding delta formula
                            del_temp = delta_UBA_neg[n-1]['del_'+t+'_'+q+'_'+\
                                                          i+'_'+qq+'_'+str(n)]

                            for tt in T:
                                # clauses resulting from annotation validity
                                cnf_temp = inc_annot(q,t,i,qq,tt,'i',n,n_bit,Q_rej,Q_rej_UCBA)

                                for d in del_temp:
                                    for c in cnf_temp:
                                        cl = d +\
                                             ' ' +\
                                             '~' + 'tau_' + t + '_' + i + '_' + tt +\
                                             ' ' +\
                                             c # disjunctive clause resulting from annotations
                                        cl_dimacs = clause_map(clause_check([cl]), mapping)

                                        n_clause += 1
                                        f.write(cl + '\n')
                                        if cl_dimacs:
                                            n_clause_mod += 1
                                            f_dimacs.write(str(max_weight) + ' ' + cl_dimacs[0] + ' 0\n')

        # for soft constraints and UCBA
        for n in range(1,n_soft_const+1):
            n_bit = n_bit_s[n-1][1]

            for q in Q_UCBA[n-1]:
                for t in T:
                    for i in I:
                        for qq in Q_UCBA[n-1]:
                            # corresponding delta formula
                            del_temp = delta_UCBA_neg[n-1]['del_'+t+'_'+q+'_'+\
                                                           i+'_'+qq+'_'+str(n)]

                            for tt in T:
                                # clauses resulting from annotation validity
                                cnf_temp = inc_annot(q,t,i,qq,tt,'f',n,n_bit,Q_rej,Q_rej_UCBA)

                                for d in del_temp:
                                    for c in cnf_temp:
                                        cl = d +\
                                             ' ' +\
                                             '~' + 'tau_' + t + '_' + i + '_' + tt +\
                                             ' ' +\
                                             c # disjunctive clause resulting from annotations
                                        cl_dimacs = clause_map(clause_check([cl]), mapping)

                                        n_clause += 1
                                        f.write(cl + '\n')
                                        if cl_dimacs:
                                            n_clause_mod += 1
                                            f_dimacs.write(str(max_weight) + ' ' + cl_dimacs[0] + ' 0\n')

        # clauses from soft constraints
        (hard_cl_temp, soft_clauses, soft_weights) =\
            soft_const(n_soft_const, n_bit_s, 't0', Q0_UBA, Q0_UCBA)
        for c in hard_cl_temp:
            c_dimacs = clause_map([c], mapping)[0]
            f.write(c + '\n')
            f_dimacs.write(str(max_weight) + ' ' + c_dimacs + ' 0\n')
        for ind_c, c in enumerate(soft_clauses):
            c_dimacs = clause_map([c], mapping)[0]
            f.write(c + '\n')
            f_dimacs.write(str(soft_weights[ind_c]) + ' ' + c_dimacs + ' 0\n')

        n_clause_mod += (len(hard_cl_temp) + len(soft_clauses))
        n_hard = n_clause + len(hard_cl_temp) # number of hard clauses
        n_soft = len(soft_clauses) # number of soft clauses
        n_clause = n_hard + n_soft # total number of clauses

    with open('dimacs_encoding_nohead.txt', 'r') as f_in, open('dimacs_encoding.txt', 'w') as f_out:
        f_out.write('p wcnf ' + str(count_int) + ' ' + str(n_clause_mod) +\
                    ' ' + str(max_weight) + '\n')
        for line in f_in:
            f_out.write(line)

    return (n_hard, n_soft, soft_weights)

def call_solver(count_int, mapping, n_soft_const):
    """Call Max-SAT solver with and interpret its result"""

    # call the Max-SAT solver and output in a text file
    if windows:
        t_s_maxsat = time.time()
        subprocess.call("bash ./maxsat_solver.sh")
        t_f_maxsat = time.time()
        dt_maxsat = t_f_maxsat - t_s_maxsat # time for solving Max-SAT instance by Open-WBO
    else:
        call_maxsat = maxsat_solver + " ./dimacs_encoding.txt > ./maxsat_output.txt"
        t_s_maxsat = time.time()
        subprocess.call(call_maxsat, shell=True)
        t_f_maxsat = time.time()
        dt_maxsat = t_f_maxsat - t_s_maxsat # time for solving Max-SAT instance by Open-WBO
    
    # get variable interpretation from "maxsat_output.txt"
    with open('maxsat_output.txt', 'r') as f_int:
        l =  ' '
        error_flag = False
        while l[0] != 's':
            l = f_int.readline()
            if l == '':
                error_flag = True
                break
        if error_flag != True:
            if l[2] == 'O': # optimum found --> satisfiable
                f_real = True

                l = f_int.readline()
                l = l[2:-2]
                interpretation = [int(j) for j in
                                  np.sign([int(val) for val in
                                           l.split(' ')])/2.0 + 0.5]

                # Max-SAT solver might output 1 redundant variable --> remove it
                if len(interpretation) != count_int:
                    interpretation = interpretation[:-1]

                # keep variables assignment with their actual name
                int_map = dict(zip(mapping.keys(),
                                   [interpretation[int(j)-1] for j in mapping.values()]))

                with open('var_assign.txt', 'w') as f_assign:
                    for k, v in sorted(int_map.items()):
                        f_assign.write(k + ' ---> '+ str(v) + '\n')

                # compute optimal value of objective function in Max-SAT instance
                sol_value = 0
                for n in range(1, n_soft_const+1):
                    if int_map['s1_'+str(n)] == 1:
                        sol_value += int(math.pow(n_soft_const,2)+n_soft_const+1)
                    elif int_map['s2_'+str(n)] == 1:
                        sol_value += int(math.pow(n_soft_const,2)+n_soft_const)
                    elif int_map['s3_'+str(n)] == 1:
                        sol_value += int(math.pow(n_soft_const,2))
                    else:
                        pass

            elif l[2] == 'U': # unsatisfiable
                f_real = False
                sol_value = 0
                print("The hard constraints are unsatisfiable!")

            else: # unknown solver output
                f_real = False
                sol_value = -1
                raise Exception("Unexpected output from the solver!")

        else:
            f_real = False
            sol_value = -1
            raise Exception("Parsing error!")

    return (f_real, sol_value, dt_maxsat)

def gen_ex_io(n_i, i_vars, n_o, o_vars):
    """Generate the input-output variables specific to the problem"""

    I = ['i'+str(j) for j in range(int(math.pow(2,n_i)))]
    I_val = {}
    for j in range(int(math.pow(2,n_i))):
        I_val[I[j]] = tuple(int(jj) for jj in tuple(format(j, "0"+str(n_i)+"b")))

    io_vars = (n_i, i_vars, I_val, n_o, o_vars)

    return (I, I_val, io_vars)

def synthesize(i_vars, I, I_val, o_vars,
               A_h_UCBA, A_s_UBA, A_s_UCBA, T, n_t,n_o):
  
    t_s_enc = time.time()
    # extract states, rejecting states and bad transitions
    Q = A_h_UCBA[0][0]
    Q_rej = A_h_UCBA[0][3]
    #
    n_soft_const = len(A_s_UBA)
    #
    Q_UBA = [A_s_UBA[i][0] for i in range(n_soft_const)]
    Q_rej_UBA = [A_s_UBA[i][3] for i in range(n_soft_const)]
    bad_tra_UBA = [A_s_UBA[i][4] for i in range(n_soft_const)]
    #
    Q_UCBA = [A_s_UCBA[i][0] for i in range(n_soft_const)]
    Q_rej_UCBA = [A_s_UCBA[i][3] for i in range(n_soft_const)]

    # bounds on annotation variables
    annot_bound_h = n_t*len(Q_rej)
    annot_bound_s = [(n_t, n_t*max(1,len(Q_rej_UCBA[n]))) for n in range(n_soft_const)] # list of (bound(li), bound(lf))
    n_bit_h = int(np.ceil(np.log2(annot_bound_h))) # size of bit vector for hard annotation
    n_bit_s = map(lambda x: ( int(np.ceil(np.log2(x[0]))), int(np.ceil(np.log2(x[1]))) ),annot_bound_s) # size of bit vector for soft annotation
                     
    # run the synthesis procedure
    (count_int, mapping) = gen_mapping(I, n_o, Q, Q_UBA, bad_tra_UBA, Q_UCBA, T, n_bit_h, n_bit_s)
    (n_hard, n_soft, soft_weights) = encode_synthesis(i_vars, I, I_val, o_vars,
                                                      A_h_UCBA, A_s_UBA, A_s_UCBA, T, n_t, count_int, mapping, n_bit_h, n_bit_s)
    n_clause = n_hard + n_soft # total number of clauses
    t_f_enc = time.time()
    dt_enc = t_f_enc - t_s_enc # time for encoding
    print count_int; print n_clause
    (f_real, sol_value, dt_maxsat) = call_solver(count_int, mapping, n_soft_const)

    return (f_real, (count_int, n_clause), sol_value, dt_enc, dt_maxsat)