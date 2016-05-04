# Setup tool for oxDNA input in LAMMPS format.

import math,numpy as np,sys,os

r0 = 0.7

# definition of single untwisted strand
def single():

  strand = inp[1].split(':')

  com_start=strand[0].split(',')

  posx=float(com_start[0])
  posy=float(com_start[1])
  posz=float(com_start[2])
  risex=0
  risey=0
  risez=r0

  strandstart=len(nucleotide)+1

  for letter in strand[1]:
    temp=[]

    temp.append(nt2num[letter])
    temp.append([posx,posy,posz])
    vel=[0,0,0,0,0,0]
    temp.append(vel)
    temp.append(shape)

    quat=[1,0,0,0]
    temp.append(quat)

    posx=posx+risex
    posy=posy+risey
    posz=posz+risez

    if (len(nucleotide)+1 > strandstart):
      topology.append([1,len(nucleotide),len(nucleotide)+1])

    nucleotide.append(temp)

  return

# definition of single twisted strand
def single_helix():

  strand = inp[1].split(':')

  com_start=strand[0].split(',')
  twist=float(strand[1])

  posx = float(com_start[0])
  posy = float(com_start[1])
  posz = float(com_start[2])
  risex=0
  risey=0
  risez=math.sqrt(r0**2-4.0*math.sin(0.5*twist)**2) 

  dcomh=0.76
  axisx=dcomh + posx
  axisy=posy

  strandstart=len(nucleotide)+1
  quat=[1,0,0,0]

  qrot0=math.cos(0.5*twist)
  qrot1=0
  qrot2=0
  qrot3=math.sin(0.5*twist)

  for letter in strand[2]:
    temp=[]

    temp.append(nt2num[letter])
    temp.append([posx,posy,posz])
    vel=[0,0,0,0,0,0]
    temp.append(vel)
    temp.append(shape)

    temp.append(quat)

    quat0 = quat[0]*qrot0 - quat[1]*qrot1 - quat[2]*qrot2 - quat[3]*qrot3 
    quat1 = quat[0]*qrot1 + quat[1]*qrot0 + quat[2]*qrot3 - quat[3]*qrot2 
    quat2 = quat[0]*qrot2 + quat[2]*qrot0 + quat[3]*qrot1 - quat[1]*qrot3 
    quat3 = quat[0]*qrot3 + quat[3]*qrot0 + quat[1]*qrot2 + quat[2]*qrot1 

    quat = [quat0,quat1,quat2,quat3]

    posx=axisx - dcomh*(quat[0]**2+quat[1]**2-quat[2]**2-quat[3]**2)
    posy=axisy - dcomh*(2*(quat[1]*quat[2]+quat[0]*quat[3]))
    posz=posz+risez

    if (len(nucleotide)+1 > strandstart):
      topology.append([1,len(nucleotide),len(nucleotide)+1])

    nucleotide.append(temp)

  return

nt2num = {'A':1, 'C':2, 'G':3, 'T':4}
shape = [1.1739845031423408,1.1739845031423408,1.1739845031423408]

nucleotide=[]
topology=[]

seqfile = open(sys.argv[1],'r')

# process sequence file line by line
for line in seqfile:

  inp = line.split()
  if inp[0] == 'single':
    single()
  if inp[0] == 'single_helix':
    single_helix()

# output atom data in LAMMPS format
out = open(sys.argv[2],'w')

out.write('%d atoms\n' % len(nucleotide))
out.write('%d ellipsoids\n' % len(nucleotide))
out.write('%d bonds\n' % len(topology))
out.write('\n')
out.write('4 atom types\n')
out.write('1 bond types\n')
out.write('\n')
out.write('-20.0 20.0 xlo xhi\n')
out.write('-20.0 20.0 ylo yhi\n')
out.write('-20.0 20.0 zlo zhi\n')
out.write('\n')
out.write('Masses\n')
out.write('\n')
out.write('1 3.1575\n')
out.write('2 3.1575\n')
out.write('3 3.1575\n')
out.write('4 3.1575\n')

out.write('\n')
out.write('Atoms\n')
out.write('\n')
for ib in range(len(nucleotide)):
  out.write(str(ib+1) + " " + str(nucleotide[ib][0]) + " " + str(nucleotide[ib][1][0]) + " " + str(nucleotide[ib][1][1]) + " " + str(nucleotide[ib][1][2]) + " 1 1 1\n")

out.write('\n')
out.write('Velocities\n')
out.write('\n')
for ib in range(len(nucleotide)):
  out.write(str(ib+1) + " " + str(nucleotide[ib][2][0]) + " " + str(nucleotide[ib][2][1]) + " " + str(nucleotide[ib][2][2]) +  " " + str(nucleotide[ib][2][3]) + " " + str(nucleotide[ib][2][4]) + " " + str(nucleotide[ib][2][5]) + "\n")

out.write('\n')
out.write('Ellipsoids\n')
out.write('\n')
for ib in range(len(nucleotide)):
  out.write(str(ib+1) + " " + str(nucleotide[ib][3][0]) + " " + str(nucleotide[ib][3][1]) + " " + str(nucleotide[ib][3][2]) + " " + str(nucleotide[ib][4][0]) + " " + str(nucleotide[ib][4][1]) + " " + str(nucleotide[ib][4][2]) + " " + str(nucleotide[ib][4][3])  + "\n")

out.write('\n')
out.write('Bonds\n')
out.write('\n')
for ib in range(len(topology)):
  out.write(str(ib+1) + " " + str(topology[ib][0]) + " " + str(topology[ib][1]) + " " + str(topology[ib][2])+ "\n")

out.close() 

seqfile.close()
sys.exit(0)


