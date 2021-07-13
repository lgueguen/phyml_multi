#the input of the program is the output lk file from phyml_multi. The output is a well formatted file that is input into Sarment. We also output the number of partitions.
#! /usr/bin/env python

import sys
import math
import os


#THREE USEFUL FUNCTIONS
def logsumexp(l): #computes the logarithm of the sum of exponentials
    m=max(l)
    return(math.log(reduce(lambda x,y:x+math.exp(y-m),l,0))+m)

def calc_pk(l,ls,k):
    " Retourne P(P_k|S) normalise"
    sl=logsumexp(map(lambda x: x[0]-x[1], zip(l,ls)))
    return(l[k-1]-ls[k-1]-sl)

def max_list(l):
	"retourne le max d'une liste"
	max = -100000000000
	index=0
	for i in range(len(l)) :
		#print str(i)+" "+str(l[i])+" "+str(max)
		if (l[i]>max):
			max = l[i]
			index=i
		#print max
	return index+1

#AND NOW THE REAL SCRIPT

argv=sys.argv[1:]
file=argv[0]
out=argv[0].split('.')[0]+'.mat'
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
lineMod=""
summ=0
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
        for j in range(numMod):
            fout.write("#"+str(j)+"\t")
        fout.write("\n")
        i=i+1
    else:
        liste=l.split()
        for j in range(numMod):
            fout.write(str(liste[1+numMod+j])+"\t")
            mod[j][i-1]=float(liste[1+numMod+j])
            summ=summ+float(liste[1+numMod+j])
        fout.write("\n")
        i=i+1

f.close()
fout.close()

if (numMod==1) :
    print summ
    

else :
    for j in range(numMod):
        lineMod+="#"+str(j)+"\t"

    print lineMod


 




    file=out

    nb_cl=100
    nshuffle=50

    import lexique
    #lx=lexique.Lexique(str="#0 #1 #2 #3 #4")
    #lx=lexique.Lexique(str="#0 #1 #2")
    lx=lexique.Lexique(str=lineMod)
    print lx
    import matrice
    m=matrice.Matrice(fic=file)

    #res=lx.probability(m,nb_cl)

    res=lx.log_likelihood(m,nb_cl)


    f=open(file+"PartitionProbabilities", 'w')
    for x in res:
	f.write("%f\t"%(x))

    f.write("\n")
    f.flush()

    totals=nb_cl * [0]

    for n in range(nshuffle):
	print "replicate "+str(n)+"\r"
	#for i in range(1,nb_cl+1):
	m.shuffle()
	#ress=lx.probability(m,nb_cl)
        ress=lx.log_likelihood(m,nb_cl)
	for i in range(1,nb_cl+1):
		temp = calc_pk(res,ress,i)
		totals[i-1]=totals[i-1]+temp
		f.write("%f\t"%(temp))
        
        f.write("\n")
        f.flush()

    f.close()

    for i in range(1,nb_cl+1):
	totals[i-1] = totals[i-1]/nshuffle

    f=open(file+"NormalizedPartitionProbabilities", 'w')
    for i in range(1,nb_cl+1):
	f.write("%f\t"%(totals[i-1]))


    f.close()
    NumPart = max_list(totals)

    m=matrice.Matrice(fic=file)
    res=lx.log_likelihood(m,nb_cl)


    numModels=100


    import parti_simp
    ps=parti_simp.Parti_simp()
    ps.mpp(m,lx,numModels)
    #print ps
    ps.draw_nf(file.split('.')[0]+"_Partitioned.ps")
    # FOR PDF OUTPUT, UNCOMMENT THE NEXT LINE
    #os.system("ps2pdf "+file.split('.')[0]+"_Partitioned.ps")
    if NumPart >1:
        print NumPart
        print str(ps[NumPart-1])
        liste=str(ps[NumPart-1]).split(" XXX ")

    f=open(file.split('.')[0]+"_PartitionNumbers", 'w')
    for i in range(len(ps)):
        f.write(str(ps[i])+"\n")
    f.close()


    f=open(file.split('.')[0]+"_PartitionBestNumber", 'w')
    f.write("%d partitions\t\n"%(NumPart))
    for i in range(NumPart):
	f.write(liste[i]+"\n")

    f.close()



# import lexique
# l=lexique.Lexique(str="#0 #1 #2 #3 #4")
# print l
# import matrice
# m=matrice.Matrice(fic=file)
# lp=l.probability(m,30)
# print lp
# import parti_simp
# ps=parti_simp.Parti_simp()
# ps.mpp(m,l,30)
# ps.draw_nf(file+"Partitioned.ps")
# fout=open(file+"PartitionProbabilities", 'w')
# fout.write("numPartitions\tAverageLogLk\tMaxlogLk\n")
# for i in range(len(lp)):
# 	fout.write(str(i)+"\t"+str(lp[i])+"\t"+str(ps.ls_val()[i])+"\n")

# fout.close()

# fout=open(file+"PartitionGeography", 'w')
# liste=str(ps[4]).split(" XXX ")
# for i in range(5):
# 	fout.write(liste[i]+"\n")

# fout.close()
