""" ===============================================
                  INITIALIZATION
=============================================== """

import init
import simulation.direct_noCoast_injection as direct_noCoast_injection
import simulation.single_run as single_run

if __name__ == '__main__':
    if init.SYM_TYPE == 1:
        single_run.execute()
    elif init.SYM_TYPE == 2:
        direct_noCoast_injection.execute()