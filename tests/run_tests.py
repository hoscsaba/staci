import os
import sys
from printf import printf


def getres():
    with open('tmp.dat', 'r') as f:
        for line in f:
            num = float(line)
    return num


def checksol(sollwert, iswert):
    if (abs(sollwert) < 1.e-3):
        err = abs(sollwert-iswert)
        printf('\n ==> expected: %g, got: %g', sollwert, iswert)
    else:
        err = abs((sollwert-iswert)/sollwert)*100.
        printf('\n ==> expected: %g, got: %g, error: %g %%', sollwert, iswert, err)
    
    if (err < 0.1):
        printf('\tOK.\n')
    else:
        printf('\n\n')
        printf('Error too large!')
        sys.exit('Error too large!')


log_file = open("run_tests.log", "w")
sys.stdout = log_file


# printf('\n************** TEST 1: channel.spr (teltszelveny) ************** \n')
# os.system('cp channel.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL16 -p bottom_level -n 2.0 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL16 -p water_level  -n 1.0 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL19 -p bottom_level -n 0.0 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL19 -p water_level  -n 2.0 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('rm tmp1.spr')
# os.system('new_staci -s tmp.spr')
# os.system('new_staci -g tmp.spr -e CHANNEL10 -p mass_flow_rate')
# checksol(11473, getres())

# printf('\n********* TEST 2: channel.spr, backwards (case teltszelveny) ****** \n')
# os.system('cp channel.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL16 -p bottom_level -n 1.0 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL16 -p water_level  -n 1.0 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL19 -p bottom_level -n 0.0 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL19 -p water_level  -n 3.0 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('rm tmp1.spr')
# os.system('new_staci -s tmp.spr')
# os.system('new_staci -g tmp.spr -e CHANNEL10 -p mass_flow_rate')
# checksol(-9484.48, getres())

# printf('\n************ TEST 3: channel.spr, empty (case 0.a.) ************ \n')
# os.system('cp channel.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL16 -p bottom_level -n 0.0 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL16 -p water_level  -n 0.5 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL19 -p bottom_level -n 0.0 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL19 -p water_level  -n 0.5 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('rm tmp1.spr')
# os.system('new_staci -s tmp.spr')
# os.system('new_staci -g tmp.spr -e CHANNEL10 -p mass_flow_rate')
# checksol(0.0, getres())

# printf('\n********* TEST 4: channel.spr, backflow (case 0.b.i & ii) ********** \n')
# os.system('cp channel.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL16 -p bottom_level -n 0.0 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL16 -p water_level  -n 0.5 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL19 -p bottom_level -n 0.0 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL19 -p water_level  -n 1.1 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('rm tmp1.spr')
# os.system('new_staci -s tmp.spr')
# os.system('new_staci -g tmp.spr -e CHANNEL10 -p mass_flow_rate')
# checksol(-22.3555, getres())

# printf('\n*********** TEST 5: channel.spr, backflow (case 1.) ************ \n')
# os.system('cp channel.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL16 -p bottom_level -n 0.0 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL16 -p water_level  -n 1.5 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL19 -p bottom_level -n 0.0 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('new_staci -m tmp.spr -e POOL19 -p water_level  -n 1.9 -o tmp1.spr')
# os.system('cp tmp1.spr tmp.spr')
# os.system('rm tmp1.spr')
# os.system('new_staci -s tmp.spr')
# os.system('new_staci -g tmp.spr -e CHANNEL10 -p mass_flow_rate')
# checksol(-1711.25, getres())

printf('\n*********** TEST 6: channel.spr, yn>yc (case 2.a.i) ************ \n')
os.system('cp channel.spr tmp.spr')
os.system('new_staci -m tmp.spr -e POOL16 -p bottom_level -n 0.0 -o tmp1.spr')
os.system('cp tmp1.spr tmp.spr')
os.system('new_staci -m tmp.spr -e POOL16 -p water_level  -n 1.5 -o tmp1.spr')
os.system('cp tmp1.spr tmp.spr')
os.system('new_staci -m tmp.spr -e POOL19 -p bottom_level -n 0.0 -o tmp1.spr')
os.system('cp tmp1.spr tmp.spr')
os.system('new_staci -m tmp.spr -e POOL19 -p water_level  -n 1.3 -o tmp1.spr')
os.system('cp tmp1.spr tmp.spr')
os.system('rm tmp1.spr')
os.system('new_staci -s tmp.spr')
os.system('new_staci -g tmp.spr -e CHANNEL10 -p mass_flow_rate')
checksol(3155.29, getres())

log_file.close()
