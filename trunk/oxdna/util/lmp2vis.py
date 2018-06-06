import sys, math, subprocess

# converts quaternion DOF into local body reference frame
def q_to_exyz(q1,q2,q3,q4):

    q=[q1, q2, q3, q4]
    ex= [0, 0, 0]
    ey = [0, 0,0]
    ez = [0, 0, 0]

    ex[0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3]
    ex[1]=2*(q[1]*q[2]+q[0]*q[3])
    ex[2]=2*(q[1]*q[3]-q[0]*q[2])

    ey[0]=2*(q[1]*q[2]-q[0]*q[3])
    ey[1]=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3]
    ey[2]=2*(q[2]*q[3]+q[0]*q[1])

    ez[0]=2*(q[1]*q[3]+q[0]*q[2])
    ez[1]=2*(q[2]*q[3]-q[0]*q[1])
    ez[2]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3]

    return ex,ey,ez

# processes line by line of LAMMPS output
def transform(line): 

    list1 = line.split()
    ident, mol, typ = int(list1[0]), int(list1[1]), int(list1[2])
    typb = typ + 10  #defines new backbone types
    x, y, z = float(list1[3]), float(list1[4]), float(list1[5])
    c_quat1, c_quat2, c_quat3, c_quat4 = \
            float(list1[6]), float(list1[7]), float(list1[8]), float(list1[9])
    ex, ey, ez = q_to_exyz(c_quat1, c_quat2, c_quat3, c_quat4)

    # position of sugar-phosphate backbone interaction site in oxDNA2
    x1, y1, z1 = x -0.34*ex[0]+0.3408*ey[0],y -0.34*ex[1]+0.3408*ey[1], z-0.34*ex[2]+0.3408*ey[2]

    # position of base interaction site in oxDNA2
    x2, y2, z2 = x +0.4*ex[0], y + 0.4*ex[1], z+0.4*ex[2]

    # compose basic output data: id, molecule id, type, position, quaternion
    line1 = '%d'%(2*ident-1) +' '+ '%d'%mol +' '+ '%d'%typb +' '+\
        '%1.6e'%(x1) +' '+ '%1.6e'%(y1) +' '+ '%1.6e'%(z1) +' '+\
        '%1.6e'%(c_quat1) +' '+ '%1.6e'%(c_quat2) +' '+ '%1.6e'%(c_quat3) +' '+ '%1.6e'%(c_quat4)
    line2 = '%d'%(2*ident) +' '+ '%d'%mol +' '+ '%d'%typ +' '+\
        '%1.6e'%(x2) +' '+ '%1.6e'%(y2) +' '+ '%1.6e'%(z2) +' '+\
        '%1.6e'%(c_quat1) +' '+ '%1.6e'%(c_quat2) +' '+' %1.6e'%(c_quat3) +' '+ '%1.6e'%(c_quat4)

    # for ovito we use oblate particles for the bases
    shape_sphere = ' 0.4 0.4 0.4'
    shape_ellipsoid = ' 0.5 0.2 0.1'

    if vismethod == 'ovito':
        line1 += shape_sphere +' '
        line2 += shape_ellipsoid +' '

    # append remaining output data
    for i in range(10, len(list1)):
        line1 += ' '+ '%1.6e'%float(list1[i])
        line2 += ' '+ '%1.6e'%float(list1[i])
    line1 += '\n'
    line2 += '\n'

    line= line1 + line2
    return line

# main part

if len(sys.argv)!=4:
    print("Syntax: $> python lmp2vis.py visualisation_method(vmd or ovito) input_filename output_filename")
    sys.exit(1)

vismethod = sys.argv[1]
r=open(sys.argv[2],'r')
w=open(sys.argv[3],'w+')

if (sys.argv[1]!='vmd' and sys.argv[1]!='ovito'):
    print("Please select visualisation method: vmd or ovito")
    print("Syntax: $> python lmp2vis.py visualisation_method input_filename output_filename")
    sys.exit(1)

print('# Converting LAMMPS output for visualisation with %s' % sys.argv[1])     

# count lines in output file for progress report
n = 0

try:
    result = subprocess.run(['wc', '-l', '%s'%sys.argv[2]], stdout=subprocess.PIPE)
    reslist=str(result).split()
    nlines=float(reslist[5])
except:
    nlines = 100

line=r.readline()

while line != '':

    sys.stdout.write('# Processed %3d %%\r' % (100*n/nlines))     

    if line.find('NUMBER OF ATOMS') != -1: 
        w.write(line)
        N=int(r.readline())
        w.write('%d'%int(2*N)+'\n')
        line=r.readline()
    if line.find('id mol type x y z') != -1:
        if vismethod == 'ovito':
            linestring=line.split()
            line = linestring[0]
            for i in range(1, 12):
                line += ' '+ linestring[i]
            line += ' shape[0] shape[1] shape[2]'
            for i in range(12, len(linestring)):
                line += ' '+ linestring[i]
            line += '\n'
        i=0
        w.write(line)
        while i<N:
            line=r.readline()
            w.write(transform(line))
            i +=1
            n += 1
        else:
            line=r.readline()
    else:
        w.write(line)
        line=r.readline()

print('# Done                                ')     

r.close()
w.close()
