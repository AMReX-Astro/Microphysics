# Read in the Helmholtz EOS table and create a
# Fortran module containing the data.

# Number of density rows
imax = 541

# Number of temperature columns
jmax = 201

table_name = 'helm_table.dat'

table = open(table_name, 'r')

# Free energy and its derivatives

f      = [[None] * jmax for _ in range(imax)]
fd     = [[None] * jmax for _ in range(imax)]
ft     = [[None] * jmax for _ in range(imax)]
fdd    = [[None] * jmax for _ in range(imax)]
ftt    = [[None] * jmax for _ in range(imax)]
fdt    = [[None] * jmax for _ in range(imax)]
fddt   = [[None] * jmax for _ in range(imax)]
fdtt   = [[None] * jmax for _ in range(imax)]
fddtt  = [[None] * jmax for _ in range(imax)]

# Pressure derivatives

dpdf   = [[None] * jmax for _ in range(imax)]
dpdfd  = [[None] * jmax for _ in range(imax)]
dpdft  = [[None] * jmax for _ in range(imax)]
dpdfdt = [[None] * jmax for _ in range(imax)]

# Electron chemical potential

ef     = [[None] * jmax for _ in range(imax)]
efd    = [[None] * jmax for _ in range(imax)]
eft    = [[None] * jmax for _ in range(imax)]
efdt   = [[None] * jmax for _ in range(imax)]

# Number density

xf     = [[None] * jmax for _ in range(imax)]
xfd    = [[None] * jmax for _ in range(imax)]
xft    = [[None] * jmax for _ in range(imax)]
xfdt   = [[None] * jmax for _ in range(imax)]

for j in range(jmax):
    for i in range(imax):

        line = table.readline().split()

        f[i][j]      = line[0]
        fd[i][j]     = line[1]
        ft[i][j]     = line[2]
        fdd[i][j]    = line[3]
        fdt[i][j]    = line[4]
        ftt[i][j]    = line[5]
        fddt[i][j]   = line[6]
        fdtt[i][j]   = line[7]
        fddtt[i][j]  = line[8]

for j in range(jmax):
    for i in range(imax):

        line = table.readline().split()

        dpdf[i][j]   = line[0]
        dpdfd[i][j]  = line[1]
        dpdft[i][j]  = line[2]
        dpdfdt[i][j] = line[3]

for j in range(jmax):
    for i in range(imax):

        line = table.readline().split()

        ef[i][j]     = line[0]
        efd[i][j]    = line[1]
        eft[i][j]    = line[2]
        efdt[i][j]   = line[3]

for j in range(jmax):
    for i in range(imax):

        line = table.readline().split()

        xf[i][j]     = line[0]
        xfd[i][j]    = line[1]
        xft[i][j]    = line[2]
        xfdt[i][j]   = line[3]

table.close()


# Now write out the module

module_file = 'helm_table.F90'
module = open(module_file, 'w')

module.write('module helm_table_module\n\n')

module.write('  implicit none\n\n')

module.write('  integer, parameter :: imax = {}\n'.format(imax))
module.write('  integer, parameter :: jmax = {}\n\n'.format(jmax))

module.write('  double precision :: f(imax,jmax) = reshape([ &\n    ')

n_per_line = 4

k = 0

for j in range(jmax):
    for i in range(imax):
        module.write('{}'.format(f[i][j].replace('E','d')))

        if not (i == imax-1 and j == jmax-1):
            module.write(', ')
        else:
            module.write('], &\n')

        k += 1

        if k == n_per_line:
            module.write(' &\n    ')
            k = 0

module.write('    [imax, jmax])\n\n')

module.write('end module helm_table_module\n')

module.close()

