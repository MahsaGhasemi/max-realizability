import re
import sympy as sp
import subprocess
import time

windows = False # OS determiner
spot = "/home/ubuntu/spot-2.4.4/bin/ltl2tgba" # Use appropriate address for automaton generator

def find_all(string, sub):
    """Find all indices of a (non-overlapping) substring occurance
    in a string"""
    s_ind = 0
    while True:
        s_ind = string.find(sub, s_ind)
        if s_ind == -1:
            return
        yield s_ind
        s_ind += len(sub)


def automaton_gen(ltl_f):
    """Generate and format an automaton for the specifications"""

    # Call SPOT
    if windows:
        with open('aut_gen_caller.sh', 'w') as f_in:
            f_in.write('#! /bin/sh' + '\n')
            f_in.write('ltl2tgba -B -f')
            f_in.write(' "' + ltl_f + '" ')
            f_in.write('-H --output=aut_output.txt')

        t_s_aut_gen = time.time()
        subprocess.call("bash ./aut_gen_caller.sh")
        t_f_aut_gen = time.time()
        dt_aut_gen = t_f_aut_gen - t_s_aut_gen # time for generating automaton by SPOT
    else:
        call_spot =  spot+" -B -f" +" '" + ltl_f +"' "+ "-H --output=./aut_output.txt"
        t_s_aut_gen = time.time()
        subprocess.call(call_spot, shell=True)
        t_f_aut_gen = time.time()
        dt_aut_gen = t_f_aut_gen - t_s_aut_gen # time for generating automaton by SPOT

    # Parse output
    #   transitions: list of 3-tuples (state,label,state)
    error_flag = False
    with open('aut_output.txt', 'r') as f_hoa:
        l = f_hoa.readline()

        if l == '':
            error_flag = True
        else:
            while l.find('States') < 0:
                l = f_hoa.readline()

            n_states = l[8:-1]

            while l.find('Start') < 0:
                l = f_hoa.readline()

            Q0 = list(l[7:-1]) # for 1 initial state

            while l.find('AP') < 0:
                l = f_hoa.readline()

            n_var = int(l[4:4+l[4:].find(' ')])
            var_all = [] # the name of all variables
            map_var = {} # mapping of variables to integers from 0 to n-1
            ind_b = [i+2 for i in list(find_all(l, ' "'))]
            ind_f = list(find_all(l, '" '))
            ind_f.append(l.rfind('"'))
            if len(ind_b) != len(ind_f):
                error_flag = True
            else:
                var_all = [l[ind_b[i]:ind_f[i]] for i in range(len(ind_b))]
                map_var = {str(i):var_all[i] for i in range(n_var)}

            while l.find('Acceptance') < 0:
                l = f_hoa.readline()

            if l[11:14] != ' 1 ':
                error_flag = True
            Q_acc = []
            Q_acc_lab = l[l.find('Inf')+4:l.find(')')] # for 1 accepting label

            l = f_hoa.readline()
            while l.find('BODY') < 0:
                l = f_hoa.readline()
            l = f_hoa.readline()

            delta = []
            n_tr = 0 # number of transitions
            while l.find('END') < 0:
                if l[7:].find(' ') < 0:
                    s_end = 7 + l[7:].find('\n')
                else:
                    s_end = 7 + l[7:].find(' ')
                s_temp = l[7:s_end]
                if l.find('{'+Q_acc_lab+'}') >= 0:
                    Q_acc.append(s_temp)
                l = f_hoa.readline()
                while (l.find('State') < 0) & (l.find('END') < 0):
                    tr_temp = l[1:l.find(']')]

                    # replace with variable names
                    tr_temp = tr_temp.replace('t','T') # variable names should not contain t
                    tr_temp = tr_temp.replace('!','~')
                    pattern = re.compile("|".join(map_var.keys()))
                    tr_temp = pattern.sub(lambda m: map_var[m.group(0)], tr_temp)

                    tr_temp =  tr_temp.split(' | ')
                    tr_temp = [' '.join(cl.split('&')) for cl in tr_temp]

                    next_s_temp = l[l.find(']')+2:-1]
                    delta.append((s_temp, tr_temp, next_s_temp))
                    n_tr += 1
                    l = f_hoa.readline()

            Q = list(set([tr[0] for tr in delta]))
            Q.sort()

            # initiate state names by 'q'
            Q = ['q'+q for q in Q]
            Q0 = ['q'+q for q in Q0]
            Q_acc = ['q'+q for q in Q_acc]
            delta = [('q'+tr[0], tr[1], 'q'+tr[2]) for tr in delta]

    return (Q, Q0, delta, Q_acc, error_flag, dt_aut_gen)

def instance_gen(f_hard, f_soft):
    """Generate an instance set of parameters for a problem"""

    if len(f_hard) > 1:
        f_hard = ' & '.join(f_hard)
    elif len(f_hard) == 1:
        f_hard = f_hard[0]
    else:
        raise Exception("No hard constraints")

    f_hard = '~( ' + f_hard + ' )'
    aut_temp = automaton_gen(f_hard)
    A_h_UCBA = [aut_temp[:-1]]
    dt_aut_gen_hard = aut_temp[-1] # time for generating automata for hard constraints

    dt_aut_gen_soft = 0 # time for generating automata for soft constraints
    A_s_UBA = []
    A_s_UCBA = []
    for fs in f_soft:
        aut_temp = automaton_gen(fs)
        A_s_UBA.append(aut_temp[:-1])
        dt_aut_gen_soft += aut_temp[-1]
        aut_temp = automaton_gen('~( ' + 'GF ( ' + fs + ' )' + ' )')
        A_s_UCBA.append(aut_temp[:-1])
        dt_aut_gen_soft += aut_temp[-1]

    return (A_h_UCBA, A_s_UBA, A_s_UCBA, (dt_aut_gen_hard, dt_aut_gen_soft))

def aut_modification(Q, Q0, delta):
    """Modifies the UBA automaton and outputs modified transition and bad transitions"""

    delta_mod = delta[:]
    delta_add = [] # transitions to be added
    delta_bad_all = [] # all bad transitions
    if len(Q0) != 1:
        error_flag = True
    else:
        error_flag = False
        q0 = Q0[0]

    for q in Q:

        trs_exist = [tr[1] for tr in delta_mod if tr[0]==q]
        trs_exist_toq0 = [tr[1] for tr in delta_mod if tr[0]==q if tr[2]==q0]

        trs_exist = ' | '.join([' | '.join(['('+' & '.join(tr_cl.split(' '))+')'
                                            for tr_cl in tr])
                                for tr in trs_exist])
        trs_exist = trs_exist.replace('T','True')
        trs_exist_toq0 = ' | '.join([' | '.join(['('+' & '.join(tr_cl.split(' '))+')'
                                                 for tr_cl in tr])
                                     for tr in trs_exist_toq0])
        trs_exist_toq0 = trs_exist_toq0.replace('T','True')

        trs_complement = str(sp.simplify('~( ' + trs_exist +' )'))
        # account for strange outputs of 'simplify'
        # can avoid this by using 'true' and false' instead
        if trs_complement == '-1':
            trs_complement = 'True'
        elif trs_complement == '-2':
            trs_complement = 'False'

        if trs_complement == 'True':
            delta_bad_all.append((q,['T'],q0))
        elif trs_complement == 'False':
            pass
        else:
            delta_bad_all.append((q,trs_complement,q0))

        if q == q0:
            trs_add = 'True'
        elif trs_exist_toq0:
            trs_add = trs_exist_toq0 + ' | ' + trs_complement
        else:
            trs_add = trs_complement

        ind_rep = [ind for ind, tr in enumerate(delta_mod) if tr[0]==q if tr[2]==q0] # existing transition to q0 that should be replaced
        trs_add = str(sp.to_dnf(trs_add, True))
        if trs_add == 'True':
            add_flag = True
            trs_add_dnf = ['T']
        elif trs_add == 'False':
            add_flag = False
        else:
            add_flag = True
            trs_add = trs_add.split(' | ')
            trs_add_dnf = []
            for cl in trs_add:
                if cl[0] == '(':
                    trs_add_dnf.append(' '.join(cl[1:-1].split(' & ')))
                else:
                    trs_add_dnf.append(' '.join(cl.split(' & ')))
        # add complementary transitions
        if add_flag:
            if ind_rep: # remove existing transition
                delta_add.append((q,trs_add_dnf,q0))
                del delta_mod[ind_rep[0]]
            else:
                delta_add.append((q,trs_add_dnf,q0))

    delta_mod.extend(delta_add)

    return (delta_mod, delta_bad_all, error_flag)

def rep_label(tra, t, i, io_vars):
    """Instances the label with input and output variables"""

    (n_i, i_vars, I_val, n_o, o_vars) = io_vars

    for i_ind in range(n_i):
        del_temp = [] # labels with valued environment variables
        for d in tra:
            d_temp = d
            if d_temp.find('~'+i_vars[i_ind]) >= 0:
                if I_val[i][i_ind] == 1:
                    d_temp = d_temp.replace('~'+i_vars[i_ind], 'F')
                else:
                    d_temp = d_temp.replace('~'+i_vars[i_ind], 'T')
            if d_temp.find(i_vars[i_ind]) >= 0:
                if I_val[i][i_ind] == 1:
                    d_temp = d_temp.replace(i_vars[i_ind], 'T')
                else:
                    d_temp = d_temp.replace(i_vars[i_ind], 'F')

            del_temp.append(d_temp)

        tra = del_temp

    for o_ind in range(n_o):
        del_temp = []
        o_rep = 'o'+str(o_ind)+'_'+t+'_'+i
        for d in tra:
            del_temp.append(d.replace(o_vars[o_ind],o_rep))

        tra = del_temp

    return del_temp

def gen_ex_del_h(Q, I, io_vars, delta_template, T):
    """Generate the extended delta for hard constraints specific to the problem"""

    # delta is initially in DNF
    #     ' ' --> and
    #      ,  --> or
    delta = {}
    for t in T:
        for q in Q:
            for i in I:
                for qq in Q:
                    k = 'del_'+t+'_'+q+'_'+i+'_'+qq
                    trs_temp = [tr for tr in delta_template
                                if tr[0]==q and tr[2]==qq]
                    if not trs_temp:
                        v = ['F']
                    elif len(trs_temp) == 1:
                        v = rep_label(trs_temp[0][1], t, i, io_vars)
                    else:
                        raise Exception("More than 1 transition between states")
                    delta[k] = v

    return delta

def gen_ex_del_bad_s(n_soft_const, Q_UBA, delta_template_UBA, delta_bad_all_UBA,
                    Q_UCBA, delta_template_UCBA, I, io_vars, T):
    """Generate the extended delta and bad transitions
       for soft constraints specific to the problem"""

    # delta is initially in DNF
    #     ' ' --> and
    #      ,  --> or
    delta_UBA = [{} for n in range(n_soft_const)]
    for n in range(n_soft_const):
        for t in T:
            for q in Q_UBA[n]:
                for i in I:
                    for qq in Q_UBA[n]:
                        k = 'del_'+t+'_'+q+'_'+i+'_'+qq+'_'+str(n+1)
                        trs_temp = [tr for tr in delta_template_UBA[n]
                                    if tr[0]==q and tr[2]==qq]
                        if not trs_temp:
                            v = ['F']
                        elif len(trs_temp) == 1:
                            v = rep_label(trs_temp[0][1], t, i, io_vars)
                        else:
                            raise Exception("More than 1 transition between states")
                        delta_UBA[n][k] = v

    # bad_tra is in CNF and represents the whole equivalency formula
    #     ' ' --> or
    #      ,  --> and
    bad_tra_template_UBA = [[] for n in range(n_soft_const)]
    for n in range(n_soft_const):
        for tr in delta_bad_all_UBA[n]:
            eq_formula = '( bad_t_q_i_qq_n >> ( ' + tr[1] + ' ) )' + ' & '+\
                         '( bad_t_q_i_qq_n << ( ' + tr[1] + ' ) )'
            eq_formula = str(sp.to_cnf(eq_formula, True))
            bad_tra_template_UBA[n].append((tr[0], eq_formula, tr[2]))
    
    bad_tra_UBA = [{} for n in range(n_soft_const)]
    for n in range(n_soft_const):
        for t in T:
            for q in Q_UBA[n]:
                for i in I:
                    for qq in Q_UBA[n]:
                        k = 'bad_'+t+'_'+q+'_'+i+'_'+qq+'_'+str(n+1)
                        trs_temp = [tr for tr in bad_tra_template_UBA[n]
                                    if tr[0]==q and tr[2]==qq]
                        if not trs_temp:
                            v = ['~'+k]
                        elif len(trs_temp) == 1:
                            v_str = rep_label([trs_temp[0][1]], t, i, io_vars)[0]
                            v_str = v_str.replace('bad_t_q_i_qq_n',k)
                            v_str = v_str.split(' & ')
                            v = []
                            for cl in v_str:
                                if cl[0] == '(':
                                    v.append(' '.join(cl[1:-1].split(' | ')))
                                else:
                                    v.append(' '.join(cl.split(' | ')))
                        else:
                            raise Exception("More than 1 bad transition between states")
                        bad_tra_UBA[n][k] = v

    # delta is initially in DNF
    #     ' ' --> and
    #      ,  --> or
    delta_UCBA = [{} for n in range(n_soft_const)]
    for n in range(n_soft_const):
        for t in T:
            for q in Q_UCBA[n]:
                for i in I:
                    for qq in Q_UCBA[n]:
                        k = 'del_'+t+'_'+q+'_'+i+'_'+qq+'_'+str(n+1)
                        trs_temp = [tr for tr in delta_template_UCBA[n]
                                    if tr[0]==q and tr[2]==qq]
                        if not trs_temp:
                            v = ['F']
                        elif len(trs_temp) == 1:
                            v = rep_label(trs_temp[0][1], t, i, io_vars)
                        else:
                            raise Exception("More than 1 transition between states")
                        delta_UCBA[n][k] = v

    return (delta_UBA, bad_tra_UBA, delta_UCBA)

def automata_for_encoding(B_h_UCBA, B_s_UBA, B_s_UCBA, I, I_val, io_vars, T):
    """Prepares and processes the automata for encoding of max-realizability for an implementation size"""
    
    # automaton for hard specification
    (Q,Q0,delta_template,Q_rej,err_flag) = B_h_UCBA[0]
    delta = gen_ex_del_h(Q, I, io_vars, delta_template, T)
    A_h_UCBA = [(Q,Q0,delta,Q_rej)]

    n_soft_const =  len(B_s_UBA)

    # automaton for soft specification
    Q_UBA =  [B_s_UBA[j][0] for j in range(n_soft_const)]
    Q0_UBA = [B_s_UBA[j][1] for j in range(n_soft_const)]
    Q_rej_UBA = [[] for j in range(n_soft_const)]
    delta_template_UBA_unmodified = [B_s_UBA[j][2] for j in range(n_soft_const)]

    # modify automaton for soft specification
    delta_template_UBA = []
    delta_bad_all_UBA = []
    for j in range(n_soft_const):
        (delta_mod, delta_bad_all, error_flag) = aut_modification(Q_UBA[j], Q0_UBA[j], delta_template_UBA_unmodified[j])
        delta_template_UBA.append(delta_mod)
        delta_bad_all_UBA.append(delta_bad_all)

    # universal co-Buchi automaton for recurrence
    Q_UCBA =  [B_s_UCBA[j][0] for j in range(n_soft_const)]
    Q0_UCBA = [B_s_UCBA[j][1] for j in range(n_soft_const)]
    Q_rej_UCBA = [B_s_UCBA[j][3] for j in range(n_soft_const)]
    delta_template_UCBA = [B_s_UBA[j][2] for j in range(n_soft_const)]

    (delta_UBA, bad_tra_UBA, delta_UCBA) = gen_ex_del_bad_s(n_soft_const, Q_UBA, delta_template_UBA, delta_bad_all_UBA,
                                                            Q_UCBA, delta_template_UCBA, I, io_vars, T)

    A_s_UBA = [(Q_UBA[n], Q0_UBA[n], delta_UBA[n], Q_rej_UBA[n], bad_tra_UBA[n])
                    for n in range(n_soft_const)]
    A_s_UCBA = [(Q_UCBA[n], Q0_UCBA[n], delta_UCBA[n], Q_rej_UCBA[n])
                    for n in range(n_soft_const)]
         
    return (A_h_UCBA, A_s_UBA, A_s_UCBA)