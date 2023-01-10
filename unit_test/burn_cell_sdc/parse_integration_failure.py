"""Parse the error message from the code and output the uhit_test
parameters needed to run that state.

We expect a VODE failure of the form:

DVODE: error test failed repeatedly or with abs(H) = HMIN
ERROR: integration failed in net
istate = -4
time = 0.0008849480914
dens start = 10525608.07
temp start = 3347664045
xn start = 1.000000073e-30 1.53526048e-15 1.560960121e-05 6.387658931e-07 4.571561053e-21 3.666142474e-05 1.20404776e-08 2.397887145e-06 0.07070344325 0.04565960122 0.009637796102 0.008085542001 9.801361948e-06 0.001118282854 0.005675213992 0.5829008684 0.2761473469 1.246029397e-13 6.784182992e-06
aux start = 0.4892089316 49.27736109 8.674515486
dens current = 10510108.57
temp current = 3379430345
xn current = -2.67730649e-16 7.796607819e-16 -0.0001043861767 6.816071749e-07 -2.402840773e-17 3.889369267e-05 2.908401171e-08 2.565283008e-06 0.0724569667 0.04737593398 0.01026178688 0.008415833899 1.025251957e-05 0.001100710668 -0.009999999999 0.5972518763 0.2737852418 -0.0006008627266 5.016066797e-06
aux current = 0.4892429892 50.69851419 8.681072691
A(rho) = -17514590.12
A(rho e) = -5.796702982e+24
A(rho X_k) = 2.628156108e-25 6.732326552e-07 -2008.845301 312.8627233 8.338951261e-13 26612.63061 8.64245201 891.6504097 25737257.52 16443175.59 3439400.749 3146205.09 2294.978782 -192362.4729 -824393.0185 -31906149.68 -33385044.41 0.0001508623748 -791.4056806
A(rho aux_k) = -8114515.432 -9919671.364 -16129133.32

or alternately, the output from a burn_t << to stdout of the form:

rho = 10525608.0691
T =   3347664045.4
xn = 1.0000000735e-30 1.53526047976e-15 1.56096012093e-05 6.38765893051e-07 4.57156105345e-21 3.66614247379e-05 1.20404775959e-08 2.39788714461e-06 0.0707034432503 0.0456
596012213 0.00963779610159 0.00808554200089 9.80136194845e-06 0.00111828285449 0.00567521399192 0.582900868414 0.276147346901 1.24602939655e-13 6.78418299203e-06
aux = 0.489208931565 49.2773610892 8.67451548604
y[SRHO] = 10525608.0691
y[SEINT] = 4.67355786331e+24
y[SFS:] = 1.05256088427e-23 1.61595500939e-08 164.300544444 6.72339943814 4.81184599124e-14 385.883788045 0.126733348139 25.239220278 744196.732786 480595.067045 101443.6
64415 85105.2461273 103.165294412 11770.6070367 59735.078187 6135386.08404 2906618.74279 1.31152170706e-06 71.4076512429
ydot_a[SRHO] = -17514590.1181
ydot_a[SEINT] = -5.79670298203e+24
ydot_a[SFS:] = 2.62815610774e-25 6.73232655207e-07 -2008.84530074 312.862723324 8.33895126133e-13 26612.6306148 8.64245201022 891.650409654 25737257.5194 16443175.5881 34
39400.74876 3146205.08998 2294.97878169 -192362.472871 -824393.018453 -31906149.6784 -33385044.4088 0.000150862374768 -791.405680585

"""

import sys

def doit(string):
    """break down the SDC VODE integration failure message"""

    # figure out if it is the VODE failre or burn_t that was provided

    is_vode = False
    rhoe = None
    tmax = None

    aux = None
    A_aux_k = None

    for line in string:
        if line.startswith("dens start"):
            is_vode = True
            break

    if is_vode:
        for line in string:
            if line.startswith("dens start"):
                density = float(line.split("=")[-1])

            elif line.startswith("temp start"):
                temperature = float(line.split("=")[-1])

            elif line.startswith("rhoe start"):
                rhoe = float(line.split("=")[-1])

            elif line.startswith("dt"):
                tmax = float(line.split("=")[-1])

            elif line.startswith("xn start"):
                _tmp = line.split("=")[-1]
                xn = [float(q) for q in _tmp.split()]

            elif line.startswith("aux start"):
                _tmp = line.split("=")[-1]
                aux = [float(q) for q in _tmp.split()]

            elif line.startswith("A(rho)"):
                A_rho = float(line.split("=")[-1])

            elif line.startswith("A(rho e)"):
                A_rhoe = float(line.split("=")[-1])

            elif line.startswith("A(rho X_k)"):
                _tmp = line.split("=")[-1]
                A_X_k = [float(q) for q in _tmp.split()]

            elif line.startswith("A(rho aux_k)"):
                _tmp = line.split("=")[-1]
                A_aux_k = [float(q) for q in _tmp.split()]

    else:

        for line in string:
            if line.startswith("rho"):
                density = float(line.split("=")[-1])

            elif line.startswith("T"):
                temperature = float(line.split("=")[-1])

            elif line.startswith("xn"):
                _tmp = line.split("=")[-1]
                xn = [float(q) for q in _tmp.split()]

            elif line.startswith("aux"):
                _tmp = line.split("=")[-1]
                aux = [float(q) for q in _tmp.split()]

            elif line.startswith("ydot_a[SRHO]"):
                A_rho = float(line.split("=")[-1])

            elif line.startswith("ydot_a[SEINT]"):
                A_rhoe = float(line.split("=")[-1])

            elif line.startswith("ydot_a[SFS:]"):
                _tmp = line.split("=")[-1]
                A_X_k = [float(q) for q in _tmp.split()]

            elif line.startswith("ydot_a[SFX:]"):
                _tmp = line.split("=")[-1]
                A_aux_k = [float(q) for q in _tmp.split()]

    print(f"unit_test.density = {density}")
    print(f"unit_test.temperature = {temperature}")
    if rhoe:
        print(f"unit_test.rhoe = {rhoe}")

    if tmax:
        print(f"unit_test.tmax = {tmax}")

    print("")

    for n, X in enumerate(xn):
        print(f"unit_test.X{n+1} = {X}")

    print("")

    if aux:
        for n, a in enumerate(aux):
            print(f"unit_test.Aux{n+1} = {a}")

    print("")

    print(f"unit_test.Adv_rho = {A_rho}")
    print(f"unit_test.Adv_rhoe = {A_rhoe}")

    print("")

    for n, A in enumerate(A_X_k):
        print(f"unit_test.Adv_X{n+1} = {A}")

    if A_aux_k:
        print("")

        for n, A in enumerate(A_aux_k):
            print(f"unit_test.Adv_Aux{n+1} = {A}")


if __name__ == "__main__":
    err_file = sys.argv[1]
    with open(err_file, "r") as f:
        lines = f.readlines()

    doit(lines)
