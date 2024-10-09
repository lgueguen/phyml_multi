#! /usr/bin/env python
#the input of the program is the output lk file from phyml_multi. 


import sys

import os
import partition
import sequence
import matrice
import lexique
import math

from scipy.optimize import minimize

argv=sys.argv[1:]
file=argv[0]
if len(argv)>=2:
  lamb = float(argv[1])
else:
  lamb=-1
  
out=argv[0].split('.')[0]+'.mat'
outvit=argv[0].split('.')[0]+'_vit.ps'
outfb=argv[0].split('.')[0]+'_fb.ps'
partvit=argv[0].split('.')[0]+'_vit.part'
partfb=argv[0].split('.')[0]+'_fb.part'


#Formatting the output file from phyml_multi for Sarment.


try:
    f=open(file, 'r')
except IOError:
    print("Unknown file: ",file)
    sys.exit()

i=0
for l in f:
    i=i+1


try:
    fout=open(out, 'w')
except IOError:
    print("Unknown file: ",out)
    sys.exit()

length=i-1

fout.write(str(i-1)+"\n")
f.close()

try:
    f=open(file, 'r')
except IOError:
    print("Unknown file: ",file)
    sys.exit()


i=0
numMod=0
for l in f:
    if i==0:
        liste=l.split()
        numMod=round((len(liste)-1)/2)
        print(numMod)
        mod=numMod * [0]
        vitmod=numMod * [0]
        fbmod=numMod * [0]
        for k in range(numMod):
            mod[k]=length*[0]
            vitmod[k]=length*[0]
            fbmod[k]=length*[0]
        lexiqueString=""
        for j in range(numMod):
            lexiqueString+=str(j)+":#"+str(j)+" "
            fout.write("#"+str(j)+"\t")
        fout.write("\n")
        i=i+1
    else:
        liste=l.split()
        for j in range(numMod):
            fout.write(str(liste[1+numMod+j])+"\t")
            mod[j][i-1]=float(liste[1+numMod+j])
        fout.write("\n")
        i=i+1

f.close()
fout.close()





lx=lexique.Lexique(str=lexiqueString)#"0:#0 1:#1 ")
#print lx

import matrice

m=matrice.Matrice(fic=out)

def forward(x):
  m_f=matrice.Matrice()
  for i in range(numMod):
    for j in range(numMod):
      if (i==j) :
        lx.g_inter(i,j,x+(1-x)/numMod)
      else :
        lx.g_inter(i,j,(1-x)/(numMod)) 
  m_f.forward(m,lx)
  return(-sum(m_f.line(len(m_f)-1).values()))

cons = ({'type': 'ineq',
         'fun' : lambda x: x - 0.001,
         'jac' : lambda x: 1},
        {'type': 'ineq',
         'fun' : lambda x: 0.9999-x,
         'jac' : lambda x: -1})

  
if lamb<=0 or lamb>=1:
  print("optimization of lambda:")
  res = minimize(forward, x0=0.5, constraints=cons, method="SLSQP")
  lamb=res.x[0]
  print("\t lambda=%f"%(lamb))

for i in range(numMod):
   for j in range(numMod):
       if (i==j) :
           lx.g_inter(i,j,lamb+(1-lamb)/numMod)
       else :
           lx.g_inter(i,j,(1-lamb)/(numMod)) 



p_vit=partition.Partition()
p_vit.viterbi(m,lx)

print("\n\nViterbi : ")
#print len(p_vit)
#print p_vit.len_don()
print(p_vit)



p_vit.draw_nf(outvit,num=1)
# FOR PDF OUTPUT, UNCOMMENT THE NEXT LINE
#os.system("ps2pdf "+outvit)
m_fb=matrice.Matrice()
m_fb.fb(m,lx)
p_fb=partition.Partition()
p_fb.read_Matrice(m_fb)
p_fb.draw_nf(outfb,num=1)
print("\n\nForward-Backward : ")
print(p_fb)

outfb=argv[0].split('.')[0]+'_fb.mat'
foutfb=open(outfb,"w")
foutfb.write(str(m_fb))
foutfb.close()

outfb=argv[0].split('.')[0]+'_fb.mat'
foutfb=open(outfb,"w")
foutfb.write(str(m_fb))
foutfb.close()

del(lx)
del(m)
del(m_fb)

# FOR PDF OUTPUT, UNCOMMENT THE NEXT LINE
#os.system("ps2pdf "+outfb)

