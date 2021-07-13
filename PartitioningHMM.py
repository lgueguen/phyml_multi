#! /usr/bin/env python
#the input of the program is the output lk file from phyml_multi. 
import sys
import os
import partition
import sequence
import matrice
import lexique
import math


argv=sys.argv[1:]
file=argv[0]
lamb = float(argv[1])
out=argv[0].split('.')[0]+'.mat'
outvit=argv[0].split('.')[0]+'_vit.ps'
outfb=argv[0].split('.')[0]+'_fb.ps'
partvit=argv[0].split('.')[0]+'_vit.part'
partfb=argv[0].split('.')[0]+'_fb.part'


#Formatting the output file from phyml_multi for Sarment.

try:
    f=open(file, 'r')
except IOError, e:
    print "Unknown file: ",file
    sys.exit()

i=0
for l in f:
    i=i+1


try:
    fout=open(out, 'w')
except IOError, e:
    print "Unknown file: ",out
    sys.exit()

length=i-1

fout.write(str(i-1)+"\n")
f.close()



try:
    f=open(file, 'r')
except IOError, e:
    print "Unknown file: ",file
    sys.exit()


i=0
numMod=0
for l in f:
    if i==0:
        liste=l.split()
        numMod=(len(liste)-1)/2
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


# lx=lexique.Lexique(fprop="prob")

# m=matrice.Matrice()
# s=sequence.Sequence(fic="lambda.seq")
# m.prediction(s,lx)

# lx2=lexique.Lexique(str="1:#1 2:#2 3:#3 4:#4")
# lx2.init_trans()

# lx2.g_trans(1,1,math.log(0.9985))
# lx2.g_trans(1,2,math.log(0.0005))
# lx2.g_trans(1,3,math.log(0.0005))
# lx2.g_trans(1,4,math.log(0.0005))
# lx2.g_trans(2,2,math.log(0.9985))
# lx2.g_trans(2,1,math.log(0.0005))
# lx2.g_trans(2,3,math.log(0.0005))
# lx2.g_trans(2,4,math.log(0.0005))
# lx2.g_trans(3,3,math.log(0.9985))
# lx2.g_trans(3,2,math.log(0.0005))
# lx2.g_trans(3,1,math.log(0.0005))
# lx2.g_trans(3,4,math.log(0.0005))
# lx2.g_trans(4,4,math.log(0.9985))
# lx2.g_trans(4,2,math.log(0.0005))
# lx2.g_trans(4,3,math.log(0.0005))
# lx2.g_trans(4,1,math.log(0.0005))

#print str(lamb)
#print type(lamb)


#Here we need to specify the probabilities of transition between the states.

for i in range(numMod):
   for j in range(numMod):
       if (i==j) :
           lx.g_inter(i,j,lamb+(1-lamb)/numMod)
       else :
           lx.g_inter(i,j,(1-lamb)/(numMod)) 





#lx.g_inter(0,0,lamb)
#lx.g_inter(0,1,(1-lamb)/2)
#lx.g_inter(0,2,(1-lamb)/2)
#lx.g_inter(1,1,lamb)
#lx.g_inter(1,0,(1-lamb)/2)
#lx.g_inter(1,2,(1-lamb)/2)
#lx.g_inter(2,2,lamb)
#lx.g_inter(2,0,(1-lamb)/2)
#lx.g_inter(2,1,(1-lamb)/2)


#print m[:10]
#print lx

p_vit=partition.Partition()
p_vit.viterbi(m,lx)

print "\n\nViterbi : "
#print len(p_vit)
#print p_vit.len_don()
print p_vit



p_vit.draw_nf(outvit,num=1)
# FOR PDF OUTPUT, UNCOMMENT THE NEXT LINE
#os.system("ps2pdf "+outvit)
m_fb=matrice.Matrice()
m_fb.fb(m,lx)
p_fb=partition.Partition()
p_fb.read_Matrice(m_fb)
p_fb.draw_nf(outfb,num=1)
print "\n\nForward-Backward : "
print p_fb

# FOR PDF OUTPUT, UNCOMMENT THE NEXT LINE
#os.system("ps2pdf "+outfb)

