""" ===============================================
                  INITIALIZATION
=============================================== """

import init
import simulation.direct_noCoast_injection as direct_noCoast_injection
import simulation.single_run as single_run
import simulation.coasting_single_burn as coasting_single_burn

if __name__ == '__main__':
    if init.SYM_TYPE == 1 or init.SYM_TYPE == 2:
        single_run.execute()
    elif init.SYM_TYPE == 3:
        direct_noCoast_injection.execute()
    elif init.SYM_TYPE == 4:
        coasting_single_burn.execute()