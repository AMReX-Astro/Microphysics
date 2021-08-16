import sys

def doit(string):
    """break down the SDC VODE integration failure message"""

    for line in string:
        if line.startswith("dens start"):
            density = float(line.split("=")[-1])

        elif line.startswith("temp start"):
            temperature = float(line.split("=")[-1])

        elif line.startswith("xn start"):
            _tmp = line.split("=")[-1]
            xn = [float(q) for q in _tmp.split()]

        elif line.startswith("A(rho)"):
            A_rho = float(line.split("=")[-1])

        elif line.startswith("A(rho e)"):
            A_rhoe = float(line.split("=")[-1])

        elif line.startswith("A(rho X_k)"):
            _tmp = line.split("=")[-1]
            A_X_k = [float(q) for q in _tmp.split()]

    print(f"unit_test.density = {density}")
    print(f"unit_test.temperature = {temperature}")

    print("")

    for n, X in enumerate(xn):
        print(f"unit_test.X{n+1} = {X}")

    print("")

    print(f"unit_test.Adv_rho = {A_rho}")
    print(f"unit_test.Adv_rhoe = {A_rhoe}")

    print("")

    for n, A in enumerate(A_X_k):
        print(f"unit_test.Adv_X{n+1} = {A}")


if __name__ == "__main__":
    err_file = sys.argv[1]
    with open(err_file, "r") as f:
        lines = f.readlines()

    doit(lines)
