""" ===============================================
                  INITIALIZATION
=============================================== """

import init
import simulation.direct_noCoast_injection as direct_noCoast_injection
import simulation.single_run as single_run
import simulation.single_run_full as single_run_full

if __name__ == '__main__':
    if init.SYM_TYPE == 1:
        single_run.execute()
    elif init.SYM_TYPE == 2:
        single_run_full.execute()
    elif init.SYM_TYPE == 3:
        direct_noCoast_injection.execute()