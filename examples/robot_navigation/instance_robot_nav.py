def instantiate_robot_nav(scenario):
    """Create an instance of robotic navigation problem with given paremeters"""

    # hard constraint: 
    start_point = 'point_0'

    # hard constraint:
    valid_space = 'G ('+'('+'point_0'+' -> '+'('+' & '.join(['~corr_1', '~corr_2', '~ex_1', '~ex_2', '~passa', '~off', '~lib'])+')'+')'+' & '+\
                        '('+'corr_1'+' -> '+'('+' & '.join(['~point_0', '~corr_2', '~ex_1', '~ex_2', '~passa', '~off', '~lib'])+')'+')'+' & '+\
                        '('+'corr_2'+' -> '+'('+' & '.join(['~point_0', '~corr_1', '~ex_1', '~ex_2', '~passa', '~off', '~lib'])+')'+')'+' & '+\
                        '('+'ex_1'+' -> '+'('+' & '.join(['~point_0', '~corr_1', '~corr_2', '~ex_2', '~passa', '~off', '~lib'])+')'+')'+' & '+\
                        '('+'ex_2'+' -> '+'('+' & '.join(['~point_0', '~corr_1', '~corr_2', '~ex_1', '~passa', '~off', '~lib'])+')'+')'+' & '+\
                        '('+'passa'+' -> '+'('+' & '.join(['~point_0', '~corr_1', '~corr_2', '~ex_1', '~ex_2', '~off', '~lib'])+')'+')'+' & '+\
                        '('+'off'+' -> '+'('+' & '.join(['~point_0', '~corr_1', '~corr_2', '~ex_1', '~ex_2', '~passa', '~lib'])+')'+')'+' & '+\
                        '('+'lib'+' -> '+'('+' & '.join(['~point_0', '~corr_1', '~corr_2', '~ex_1', '~ex_2', '~passa', '~off'])+')'+')'+')'

    # hard constraint: 
    allowed_tra = 'G ( ' +\
                  '('+'point_0'+' -> '+' X '+\
                      '('+' | '.join(['point_0', 'corr_1'])+')'+')'+' & '+\
                  '('+'ex_1'+' -> '+' X '+\
                      '('+' | '.join(['ex_1', 'corr_1', 
                                      'passa', 'lib'])+')'+')'+' & '+\
                  '('+'ex_2'+' -> '+' X '+\
                      '('+' | '.join(['ex_2', 'passa',
                                      'lib', 'corr_2'])+')'+')'+' & '+\
                  '('+'corr_1'+' -> '+' X '+\
                      '('+' | '.join(['corr_1', 'ex_1', 
                                      'off'])+')'+')'+' & '+\
                  '('+'corr_2'+' -> '+' X '+\
                      '('+' | '.join(['corr_2', 'ex_2', 
                                      'point_0'])+')'+')'+' & '+\
                  '('+'passa'+' -> '+' X '+\
                      '('+' | '.join(['passa', 'ex_1', 
                                      'ex_2'])+')'+')'+' & '+\
                  '('+'off'+' -> '+' X '+\
                      '('+' | '.join(['off', 'corr_1'])+')'+')'+' & '+\
                  '('+'lib'+' -> '+' X '+\
                      '('+' | '.join(['lib', 'ex_1', 
                                      'ex_2'])+')'+')'+\
                  ' )'

    # hard constraint: 
    tour = '('+' & '.join(['GF '+'ex_1', 'GF '+'ex_2'])+')'

    # hard constraint: 
    ordering = 'G ( '+'ex_1'+' -> '+'X '+'('+'~'+'point_0'+' U '+'ex_2'+')'+' )'+' & '+\
               'G ( '+'ex_2'+' -> '+'X '+'('+'~'+'ex_1'+' U '+'point_0'+')'+' )'+' & '+\
               'G ( '+'point_0'+' -> '+'X '+'('+'~'+'ex_2'+' U '+'ex_1'+')'+' )'

    # hard constraint:
    access = '('+'~'+'ex_2'+' U '+'off'+')'

    # hard constraint:
    nouse_ret = 'G (lib | passa -> X ~ ex_1) & G (ex_2 -> X ~(lib | passa))'

    # soft constraint: 
    off_no_pass = '('+'corr_1'+' -> '+'X '+'~'+'off'+')'

    # soft constraint: 
    lib_no_pass = 'ex_1'+' | '+'ex_2'+' -> '+'X '+'~'+'lib'

    # soft constraint:
    pass_no_pass = '('+'ex_1'+' | '+'ex_2'+')'+' & '+'X '+'pass_occ'+' -> '+'X '+'~'+'passa'

    f_hard = [' & '.join([start_point, valid_space, allowed_tra, tour, access, nouse_ret])]
    f_soft = [off_no_pass, lib_no_pass, pass_no_pass]

    return (f_hard, f_soft)

def generate_var_name():
    """Generate the names of input and output variables for an instance of robotic navigation"""

    i_vars = ['pass_occ']
    o_vars = ['point_0', 'corr_1', 'corr_2', 'ex_1', 'ex_2', 'passa', 'off', 'lib']

    return (i_vars, o_vars)