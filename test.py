from max_synthesis import *
from aut_translator import *

def test_inst_gen():
    f_hard = ['G aa']
    f_soft = ['!aa']
    (A_h_UCBA,A_s_UBA,A_s_UCBA)=instance_gen(f_hard, f_soft)
    print (A_h_UCBA[0])
    print (A_s_UBA[0])
    print (A_s_UCBA[0])

def test_aut(ltl_f):
    if not windows:
        call_spot =  spot+" -B -H" +" '" + ltl_f +"' "
        subprocess.call(call_spot,shell=True)
    (Q,Q0,delta,Q_acc,error_flag) = automaton_gen(ltl_f)
    print ("-----------------------------------")
    print ("Q: "+str(Q))
    print ("Q0: "+str(Q0))
    print ("delta: "+str(delta))
    print ("Q_acc: "+str(Q_acc))
    print ("err_flag: "+str(error_flag))
    return (Q,Q0,delta,Q_acc,error_flag)

def test_modification(Q, Q0, delta):
    (delta_mod, delta_bad_all, error_flag) = aut_modification(Q, Q0, delta)
    print ("-----------------------------------")
    print ("delta_mod: "+str(delta_mod))
    print ("delta_bad_all: "+str(delta_bad_all))
    print ("err_flag: "+str(error_flag))
    return (delta_mod, delta_bad_all, error_flag)

def test_gen_delta_h(Q, I, io_vars, delta_template, T):
    delta_h = gen_ex_del_h(Q, I, io_vars, delta_template, T)
    print ("-----------------------------------")
    print ("delta_h:")
    for (k,v) in delta_h.iteritems():
        print str(k)+"    "+str(v)

def test_gen_delta_s(n_soft_const, Q_UBA, delta_template_UBA, delta_bad_all_UBA,
                    Q_UCBA, delta_template_UCBA, I, io_vars, T):
    (delta_UBA, bad_tra_UBA, delta_UCBA) = gen_ex_del_bad_s(n_soft_const, Q_UBA, delta_template_UBA, delta_bad_all_UBA,
                    Q_UCBA, delta_template_UCBA, I, io_vars, T)
    for i in range(n_soft_const):
        print ("-----------------------------------")
        print ("delta template UBA: ")
        print (delta_template_UBA[i])
        print ("delta UBA: ")
        for (k,v) in delta_UBA[i].iteritems():
            print str(k)+"    "+str(v)
        print ("-----------------------------------")
        print ("delta bad all UBA: ")
        print (delta_bad_all_UBA[i])
        print ("bad tra UBA: ")
        for (k,v) in bad_tra_UBA[i].iteritems():
            print str(k)+"    "+str(v)
        print ("-----------------------------------")
        print ("delta template UCBA: ")
        print (delta_template_UCBA[i])
        print ("delta UCBA: ")
        for (k,v) in delta_UCBA[i].iteritems():
            print str(k)+"    "+str(v)


if __name__ == '__main__':
    f_g = 'G(aa->Xbb)'
    f_ngf = '~GF(aa->Xbb)'
    (Q, Q0, delta_template, Q_acc, error_flag) = test_aut(f_g)

    i_vars = ['aa'] # name of input variables
    n_i = len(i_vars) # number of binary varibles of the input
    o_vars = ['bb'] # name of output variables
    n_o = len(o_vars) # number of binary varibles of the input

    (I, I_val, O, O_val, io_vars) = gen_ex_io(n_i, i_vars, n_o, o_vars)
    T = ['t']

    test_gen_delta_h(Q, I, io_vars, delta_template, T)

    (delta_mod, delta_bad_all, error_flag) = test_modification(Q, Q0, delta_template)

    (Q_UCBA, Q0_UCBA, delta_template_UCBA, Q_acc_UCBA, error_flag) = test_aut(f_ngf)

    test_gen_delta_s(1, [Q], [delta_mod], [delta_bad_all], [Q_UCBA], [delta_template_UCBA], I, io_vars, T)