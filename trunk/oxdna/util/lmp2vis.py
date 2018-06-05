import sys, math

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
    ident, mol, typ = int(list1[0]), list1[1],list1[2]
    typb = typ + '1'  #defines new backbone types
    x, y, z = float(list1[3]), float(list1[4]), float(list1[5])
    c_quat1, c_quat2, c_quat3, c_quat4 = \
            float(list1[6]), float(list1[7]), float(list1[8]), float(list1[9])
    etc = list1[10]
    ex, ey, ez = q_to_exyz(c_quat1, c_quat2, c_quat3, c_quat4)

    # position of sugar-phosphate backbone interaction site in oxDNA2
    x1, y1, z1 = x -0.34*ex[0]+0.3408*ey[0],y -0.34*ex[1]+0.3408*ey[1], z-0.34*ex[2]+0.3408*ey[2]

    # position of base interaction site in oxDNA2
    x2, y2, z2 = x +0.4*ex[0], y + 0.4*ex[1], z+0.4*ex[2]

    line1 = str(2*ident-1)+' '+ mol+' '+typb+' '+\
        '%1.6e'%(x1)+' '+'%1.6e'%(y1)+' '+'%1.6e'%(z1)+' '+\
        '%1.6e'%(c_quat1)+' '+'%1.6e'%(c_quat2)+' '+'%1.6e'%(c_quat3)+' '+'%1.6e'%(c_quat4)
    line2 = str(2*ident)+' '+mol+' '+typ+' '+\
        '%1.6e'%(x2)+' '+'%1.6e'%(y2)+' '+'%1.6e'%(z2)+' '+\
        '%1.6e'%(c_quat1)+' '+'%1.6e'%(c_quat2)+' '+'%1.6e'%(c_quat3)+' '+'%1.6e'%(c_quat4)

    for i in range(10, len(list1)):
        line1 += ' '+ '%1.6e'%float(list1[i])
        line2 += ' '+ '%1.6e'%float(list1[i])
    line1 += '\n'
    line2 += '\n'

    line= line1 + line2
    return line

# main part

if len(sys.argv)!=3:
    print("Syntax: $ python lmp2vis.py input_filename output_filename")

r=open(sys.argv[1],'r')
w=open(sys.argv[2],'w+')

line=r.readline()

while line != '':
    if line.find('NUMBER OF ATOMS') != -1: 
        w.write(line)
        N=int(r.readline())
        w.write(str(2*N)+'\n')
        line=r.readline()
    if line.find('id mol type x y z') != -1:
       i=0
       w.write(line)
       while i<N:
           line=r.readline()
           w.write(transform(line))
           i +=1
       else:
           line=r.readline()
    else:
        w.write(line)
        line=r.readline()

r.close()
w.close()
