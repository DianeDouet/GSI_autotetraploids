# -*- coding: utf-8 -*-
"""
Created on Sun Jul  6 15:23:35 2025

@author: douet

"""

"""
Correction apportée: seules les fitness des tétraploïdes =! 0 (càd qu'il y a des tétraploïdes dans la population) sont comptées dans la moyenne
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 10:18:20 2025

@author: diane_d
"""

"""
Ce code a les changements suivants: 
    - Limitation pollinique: une fois la mère choisie, le père peut être changé un certain nombre de fois correspondant à maxAttempt
    - Lorsque le père est reselectionné, le fait que la gamète de la mère soit réduit ou non réduit et également tiré à nouveau
    - On met en sortie la fitness des diploïdes (à l'équilibre après les 4000 première générations), et celles de tétraploïdes (moyenne sur 10000 générations à la fin de la simu)
    - On se place ici dans le cas où le taux d'autofécondation est fixé et ne peut pas muter
    - On veut tester différentes scénarios de limitations polliniques: 1, 5 et 20 essais de pères avant de recommencer et de sélectionner une nouvelle mère                                             
"""

import numpy as np
import random as rd


#performs mutations on the L fitness loci in the genome
def mutation(haplotype, U, L, varAddEff):
    nbMut = np.random.poisson(U) #number of mutations
    posMut =  rd.sample([i for i in range(L)],nbMut) #positions of the mutations

    if nbMut !=0:
        for p in posMut:
            haplotype[p] = haplotype[p] + np.random.normal(0,varAddEff**0.5)

    return np.array(haplotype)


#peforms mutations on the S-genotypes: each S-allele can mutate to another S-allele
def mutation_S(genotypes_S, U_SI, k): #genotype_S here is the list of S alleles for each of the N individuals of the population
    flat_genoS = [] #flatten the list genotypes_S
    for i in range(Npop):
        flat_genoS += genotypes_S[i]

    list_S = [i for i in range(k)]
    nbMut = np.random.binomial(len(flat_genoS),U_SI) #number of mutations
    for _ in range(nbMut):
        indMut = rd.sample([i for i in range(Npop)],1)[0]#choosing the individual
        allMut = rd.sample([i for i in range(len(genotypes_S[indMut]))],1)[0] #chossing the allele of that individual that will mutate
        mutant_S = rd.sample(list_S,1)[0] #Choosing a random S allele to replace the existing one
        while mutant_S == genotypes_S[indMut][allMut]: #Making sure that the new S allele is different from the old one
            mutant_S = rd.sample(list_S,1)[0]
        genotypes_S[indMut][allMut] = mutant_S #Replacing with mutant S allele

    return genotypes_S



#returns the fitness of an individual
def fitness(phenotype,opt,om_2):
    return np.exp(-((phenotype-opt)**2)/(2*om_2))



#produces the offspring of the next generation
def reproduction(Npop, phenotypes, om_2, opt, genotypes_S, ploidy, genome,phase, prob_unreduced, list_selfRates):
    listOff = [] #the list containing the genomes of the next generation
    Snext = [] #The list containing the S alleles of the next generation
    fitnessMax = np.max(fitness(phenotypes,opt,om_2))


    Npop_off = 0 #number of offspring
    Npop_self = 0 

    while Npop_off != Npop:
        par1 = np.random.randint(0,Npop) #Selecting a parent randomly
        fitnessPar1 = fitness(phenotypes[par1],opt,om_2)

        while fitnessPar1/fitnessMax < np.random.rand(): #A new maternal parent is selected if the fitness of the first selected parent is too low 
            par1 = np.random.randint(0,Npop)
            fitnessPar1 = fitness(phenotypes[par1],opt,om_2)

        ploidy_par1 = len(genotypes_S[par1]) #ploidy level of the mother
        listAlleles_par1 = genotypes_S[par1] #S alleles in the mother

        if np.random.rand() <= prob_unreduced and phase == 2: #Producing an unreduced gamete
            Sallele_par1 = listAlleles_par1 #creating the S_alleles transmitted by the mother to the next generation
            gamete_par1 = ploidy_par1 #number of chromosomes in the gamete of par1
        else: #Producing a reduced gamete
            Sallele_par1 = rd.sample(listAlleles_par1, ploidy_par1//2) 
            gamete_par1 = ploidy_par1//2



        selfRate = list_selfRates[par1]
        isSC =  selfing(par1, genotypes_S) #Checking if the parent is self-compatible or not
        maxAttempt = 1 #Determines the degree of pollen limitation


        #OUTCROSSING
        if selfRate <= np.random.rand() or isSC == False:
            par2 = np.random.randint(0,Npop) #A second parent is selected
            ploidy_par2 = len(genotypes_S[par2])
            fitnessPar2 = fitness(phenotypes[par2],opt,om_2)
            listAlleles_par2 = genotypes_S[par2]
            if np.random.rand() <= prob_unreduced and phase == 2: #Producing an unreduced gamete
                Sallele_par2 = rd.sample(listAlleles_par2, ploidy_par2) #creating the S_alleles transmitted by the pollen to the next generation
                gamete_par2 = ploidy_par2
            else: #producing a reduced gamete
                Sallele_par2 = rd.sample(listAlleles_par2, ploidy_par2//2)
                gamete_par2 = ploidy_par2//2

            nbAllele_diff = len(set(Sallele_par2))

            OK_parents = False

            if par2 != par1: #Parent1 different from parent2
                if fitnessPar2/fitnessMax >= np.random.rand(): #Fitness of parent 2 high enough
                    if gamete_par1 == gamete_par2 and gamete_par1 <=2: #Two parents have the same number of chromosomes in their gametes (at most diploid)
                        if nbAllele_diff >= 2 or Sallele_par2[0] not in listAlleles_par1: #Parent 2 compatible with parent 1 (at least two different S alleles or if the same, it needs to be different than all of the S alelles of the mother)
                            OK_parents = True

            loop = 0
            #Selecting a new paternal parent maxAttempt times (depending on the degree of pollen limitation) or until the production of a new individual
            while OK_parents == False and loop < maxAttempt:
                loop+=1
                #Choosing the size of the first parent's gamete
                if np.random.rand() <= prob_unreduced and phase == 2: #unreduced gamete
                    Sallele_par1 = listAlleles_par1 #creating the S_alleles transmitted by the mother to the next generation
                    gamete_par1 = ploidy_par1 #number of chromosomes in the gamete of par1
                else: #reduced gamete
                    Sallele_par1 = rd.sample(listAlleles_par1, ploidy_par1//2) 
                    gamete_par1 = ploidy_par1//2
                par2 = np.random.randint(0,Npop) #A second parent is selected
                ploidy_par2 = len(genotypes_S[par2])
                fitnessPar2 = fitness(phenotypes[par2],opt,om_2)
                listAlleles_par2 = genotypes_S[par2]
                #Choosing the size of the second parent's gamete
                if np.random.rand() <= prob_unreduced and phase == 2: #unreduced gamete
                    Sallele_par2 = rd.sample(listAlleles_par2, ploidy_par2) #creating the S_alleles transmitted by the pollen to the next generation
                    gamete_par2 = ploidy_par2
                else: #reduced gamete
                    Sallele_par2 = rd.sample(listAlleles_par2, ploidy_par2//2)
                    gamete_par2 = ploidy_par2//2
                nbAllele_diff = len(set(Sallele_par2))
                if par2 != par1: #Parent1 different from parent2 (outcrossing)
                    if fitnessPar2/fitnessMax >= np.random.rand(): #Fitness of parent 2 high enough
                        if gamete_par1 == gamete_par2 and gamete_par1 <= 2: #Both gametes of the same size and at most diploid
                            if nbAllele_diff >= 2 or Sallele_par2[0] not in listAlleles_par1: #Compatibility: diploid gamete with two different S alleles, or gamete with one or two same S alleles not present in parent 1
                                OK_parents = True


            if OK_parents == True: #creating a new individual for the next generation
              Npop_off += 1
              #Creating and adding the S genotype of the offspring
              Snext.append(Sallele_par1+Sallele_par2)
              #Creating the genome of the offspring
              listChromosomes1 = genome[par1]
              off_geno = crossover(listChromosomes1, L, ploidy_par1)[:gamete_par1]
              listChromosomes2 = genome[par2]
              off_geno += crossover(listChromosomes2, L, ploidy_par2)[:gamete_par2]
              new_ind = []
              #Mutation on the offspring's genome
              for i in range(len(off_geno)):
                  if phase >= 1:
                      off_geno[i] = mutation(off_geno[i], U, L, varAddEff)
                  new_ind.append(off_geno[i])
              listOff.append(new_ind) #adding the list of chromosomes of the new individual



        #SELFING
        else:
            loop = 0
            repro = False
            #maxAttempt = 50
            while loop < maxAttempt and repro == False:
                loop+=1
                #Choosing the first gamete
                if np.random.rand() <= prob_unreduced and phase == 2:
                    Sallele_par1 = listAlleles_par1 #creating the S_alleles of the next generation
                    gamete_par1 = ploidy_par1 #number of chromosomes in the gamete of par1
                else:
                    Sallele_par1 = rd.sample(listAlleles_par1, ploidy_par1//2) #creating the S_alleles of the next generation
                    gamete_par1 = ploidy_par1//2
                #Choosing the second gamete
                if np.random.rand() <= prob_unreduced and phase == 2:
                    gamete2_par1 = ploidy_par1 #number of chromosomes in the gamete of par1
                else:
                    gamete2_par1 = ploidy_par1//2
                if gamete_par1 == 2 and gamete_par1 == gamete2_par1: #The offspring can be tetraploids at most
                    Spollen_par1 = rd.sample([genotypes_S[par1][c] for c in range(ploidy_par1)], gamete_par1) #Choosing the S-alleles of the pollen
                    while Spollen_par1[0] == Spollen_par1[1]: #Self-incompatible pollen
                        Spollen_par1 = rd.sample([genotypes_S[par1][c] for c in range(ploidy_par1)], gamete_par1)

                    repro = True
                    #Adding a new individual to the next generation
                    Npop_off +=1
                    Npop_self+=1
                    Snext.append(Sallele_par1 + Spollen_par1)
                    listChromosomes1 = genome[par1] 
                    off_geno = crossover(listChromosomes1, L, ploidy_par1)[:gamete_par1]
                    off_geno += crossover(listChromosomes1, L, ploidy_par1)[:gamete_par1]
                    new_ind = []
                    #Mutation on the offspring's genome
                    for i in range(len(off_geno)):
                        if phase >= 1:
                            off_geno[i] = mutation(off_geno[i], U, L, varAddEff)
                        new_ind.append(off_geno[i])
                    listOff.append(new_ind)
                        

    k = 100 #number of possible S-alleles
    Snext = mutation_S(Snext, U_SI, k) #Mutations at the S locus
    freq_self = Npop_self/Npop_off

    return listOff, Snext, freq_self


#creates the new chromosomes obtained after random permutation
def crossover(listChromosomes, L, ploidy):
    all_haplo = [[] for _ in range(ploidy)]
    chosen_allels = [i for i in range(ploidy)]
    for i in range(L):
        np.random.shuffle(chosen_allels)
        for h in range(ploidy):
            all_haplo[h].append(listChromosomes[chosen_allels[h]][i])
    return all_haplo


#returns the genotyped of the offspring
def offspringGenotype(genome, Npop):
    off_geno = []

    for i in range(Npop):
        sum_geno = 0
        ploidy = len(genome[i])
        if ploidy == 2:
            dosage =1
        elif ploidy == 4:
            dosage = 0.67
        for c in range(ploidy):
            sum_geno += np.sum(genome[i][c])
        off_geno.append(dosage*sum_geno)

    return np.array(off_geno)


#returns the phenotyped of the offspring
def offspringPhenotype(genotype, Npop):
    off_pheno = [0 for _ in range(Npop)]

    for i in range(Npop):
        off_pheno[i] = genotype[i]+np.random.normal(0,1)

    return np.array(off_pheno)


#returns True if the individual i is self-compatible
def selfing(i, genotypes_S): #i is the individuals number in the list (similar to par in selection)
    ploidy = len(genotypes_S[i])
    if ploidy == 2:
        return False
    Salleles = [genotypes_S[i][c] for c in range(ploidy)]
    return len(set(Salleles)) != 1


#returns the frequencies of tetraploid individuals, of SC individuals in the population, and of SC pollen in SC individuals 
def frequencies(genotypes_S, Npop):
    nb_tetra = 0
    nbSC_ind = 0
    list_SCpollen = []
    for n in range(Npop):
        #frequency of tetraploids in the population
        if len(genotypes_S[n]) == 4:
            nb_tetra +=1
            ploidy = 4
            isSC = selfing(n, genotypes_S)
            #if the tetraploids is SC, compute nb of SC individualas and proportion of SC pollen
            if isSC is True:
                nbSC_ind +=1
                nbSC_pollen = 0
                Salleles = [genotypes_S[n][c] for c in range(ploidy)] #list of the S alleles of an individual
                maleGenotypes = [] #Listing all the possible pollen genotypes for that individual
                for j in range(ploidy-1):
                    for k in range(j+1,ploidy):
                        pollen = [Salleles[j], Salleles[k]]
                        maleGenotypes.append(pollen)
                for i in range(len(maleGenotypes)):
                    if maleGenotypes[i][0]!= maleGenotypes[i][1]: #If the two S alleles in the pollen are different, then the individual is SC
                        nbSC_pollen+=1
                list_SCpollen.append(nbSC_pollen/len(maleGenotypes)) 
        freqSC_pollen = np.sum(list_SCpollen)/nbSC_ind #Frequency of SC pollen in an SC individual
    else:
        freqSC_pollen = 0
    freq_tetra =  nb_tetra/Npop
    freqSC_ind = nbSC_ind/Npop
    return freq_tetra, freqSC_ind, freqSC_pollen


#computes inbreeding depression using a sample of the population
def inbreedingDepression(genome, om_2, Npop,L, genotypes_S):
    nbSample = 100
    if nbSample> Npop:
        nbSample = Npop

    fitness_allof = 0
    fitness_autof = 0

    #outcrossing
    for _ in range(nbSample):
        par1 = np.random.randint(0,Npop)
        par2 = np.random.randint(0,Npop)
        while par1 == par2 or len(genome[par1]) != len(genome[par2]):
            par1 = np.random.randint(0,Npop)
            par2 = np.random.randint(0,Npop)

        ploidy = len(genome[par1])
        if ploidy == 2:
            dosage = 1
        else: #tetraploids
            dosage = 0.67

        listChromosomes1 = genome[par1]
        haplo = crossover(listChromosomes1, L, ploidy)[:ploidy//2]

        listChromosomes2 = genome[par2]
        haplo += crossover(listChromosomes2, L, ploidy)[:ploidy//2]

        genotype_allof = 0

        for i in range(ploidy):
            genotype_allof += np.sum(haplo[i])
        genotype_allof *= dosage
        fitness_allof += np.exp(-(genotype_allof**2)/(2*om_2))


    #selfing
    for _ in range(nbSample):
        par1 = np.random.randint(0,Npop)
        par2 = par1

        ploidy = len(genome[par1])
        if ploidy == 2:
            dosage = 1
        else: #tetraploids
            dosage = 0.67

        listChromosomes1 = genome[par1]
        haplo = crossover(listChromosomes1, L, ploidy)[:ploidy//2]

        listChromosomes2 = genome[par2]
        haplo += crossover(listChromosomes2, L, ploidy)[:ploidy//2]

        genotype_autof = 0

        for i in range(ploidy):
            genotype_autof += np.sum(haplo[i])
        genotype_autof *= dosage
        fitness_autof += np.exp(-(genotype_autof**2)/(2*om_2))

    inbreeding_depression = 1 - fitness_autof/fitness_allof
    return inbreeding_depression



#Computes the mean fitness of diploid and tetraploid individuals
def get_meanFitnes(genotypes_S, phenotypes,opt,om_2):
    sumFitDiplo = 0
    sumFitTetra = 0
    nbDiplo, nbTetra = 0,0
    for i in range(Npop):
        if len(genotypes_S[i]) == 2:
            sumFitDiplo += fitness(phenotypes[i],opt,om_2)
            nbDiplo += 1
        else:
            sumFitTetra += fitness(phenotypes[i],opt,om_2)
            nbTetra += 1
    if nbDiplo == 0:
        fitDiplo = 0
    else:
        fitDiplo = sumFitDiplo/nbDiplo
        fitDiplo = fitDiplo.tolist()
    if nbTetra == 0:
        fitTetra = 0
    else:
        fitTetra = sumFitTetra/nbTetra
        fitTetra = fitTetra.tolist()
    return fitDiplo, fitTetra





def simulation(nbSim, L, varAddEff, dosage, Npop, U, om_2 , ploidy, U_SI, prob_unreduced, alpha):
    list_nbSalleles = [] #contains the number of S alleles at equilibrium for each simulation
    list_freqSC = [] #contains the final frequency of SC individuals in the pop after each simulation
    list_freqPollenSC =[] #contains the final proportion of SC pollen on average in an SC individual after each simulation
    inDep = [] #contains inbreeding depression for each simulation
    list_freqTetra = [] #contains the frequency of tetraploids in the population for each simulation
    list_fitDiplo, list_fitTetra = [],[] #contains the fitness of diploids and tetraploids for each simulation



    for k in range(nbSim):
        it = 0
        freqSC = []
        freqPollenSC = []
        freqTetra = []
        list_ID = []
        meanFit = []
        list_nb_Salleles = []
        print('Sim nb=', k)
        meanFitness = []
        opt = 0

        ploidy = 2

        #Initialisation (OK)
        genoInit = [[[0 for _ in range(L)],[0 for _ in range(L)]] for _ in range(Npop)] #initial diploid genomes for all individuals
        genome = np.array(genoInit)
        phenotype = np.random.normal(0,1,size = Npop) #initial phenotypes of all individuals
        #Initialisation at the S locus
        nbS_init = 100
        genotypes_S = []
        for _ in range(Npop):
            list_Sind = []
            for _ in range(ploidy):
                Sallele = np.random.choice(nbS_init)
                list_Sind.append(Sallele)
            genotypes_S.append(list_Sind)
        list_selfRates = [alpha for _ in range(Npop)]

        #Checking that the number of S-alleles in the population is sufficient
        flat_genoS=[]
        for i in range(len(genotypes_S)):
            flat_genoS += genotypes_S[i]
        if len(set(flat_genoS))<=2:
            print('No fertilization possible')
            return 'e','e','e'


        generation = 0
        phase = 0 #starting with the diploid phase

        while phase==0 or phase == 1:
            generation += 1
            if generation%100==0:
                print('generation=', generation)

            fit = fitness(phenotype,opt,om_2) #list of fitness of all individuals
            meanFitness.append(np.mean(fit)) #mean fitness of the population

            flat_genoS=[]
            for i in range(len(genotypes_S)):
                flat_genoS += genotypes_S[i]
            list_nb_Salleles.append(len(set(flat_genoS)))


            #creating the next generation
            genome, genotypes_S, freq_self = reproduction(Npop, phenotype, om_2, opt, genotypes_S, ploidy, genome,phase,prob_unreduced, list_selfRates)
            #computing the genotypes of the offspring
            genotype = offspringGenotype(genome, Npop)
            #computing their phenotypes
            phenotype = offspringPhenotype(genotype, Npop)

            flat_genoS=[]
            for i in range(len(genotypes_S)):
                flat_genoS += genotypes_S[i]
            if len(set(flat_genoS))<=2:
                print('No fertilization possible')
                return 'e','e','e'

            if generation > 2000:
                if phase == 0:
                    print('Going to phase 1')
                    phase+= 1
                    generation = 0
                if phase == 1 and generation != 0: 
                    print('Going to phase 2')
                    generation = 0 
                    #Counting the number of S alleles in the population
                    flat_genoS=[]
                    for i in range(len(genotypes_S)):
                        flat_genoS += genotypes_S[i]
                    nbSalleles = len(set(flat_genoS))
                    phase += 1 
                    freq_tetra, freqSC_ind, freqSC_pollen = frequencies(genotypes_S, Npop)
                    it = 0
                    meanFitDiplo, notused = get_meanFitnes(genotypes_S, phenotype,opt,om_2) #Computing the mean fitness of diploids after phase 1

        while phase == 2: #starting production of unreduced gametes, and introduction of tetraploids into the population
            max_gen = 200000
            generation += 1
            if generation%100==0:
                print(generation)
            fit = fitness(phenotype,opt,om_2) #list of fitness of all individuals
            meanFitness.append(np.mean(fit)) #mean fitness of the population

            #creating the next generation
            genome, genotypes_S, freq_self = reproduction(Npop, phenotype, om_2, opt, genotypes_S, ploidy,genome,phase, prob_unreduced, list_selfRates)
            #computing the genotypes of the offspring
            genotype = offspringGenotype(genome, Npop)
            #computing their phenotypes
            phenotype = offspringPhenotype(genotype, Npop)


            #Computing the frequencies of SC individuals in the population and of the SC pollen in an SC individual on average
            freq_tetra, freqSC_ind, freqSC_pollen = frequencies(genotypes_S, Npop)
            if generation == 1:
              print("freqTetra_gen1=", freq_tetra)
            ID = inbreedingDepression(genome, om_2, Npop,L, genotypes_S)



            if (generation%100==0 and generation >= max_gen-10000) or (generation%100==0 and freq_tetra>=1):
                it+=1
                freqSC.append(freqSC_ind)
                freqPollenSC.append(freqSC_pollen)
                list_ID.append(ID)
                freqTetra.append(freq_tetra)
                notused2, meanFitTetra = get_meanFitnes(genotypes_S, phenotype,opt,om_2)
                if meanFitTetra != 0:
                    meanFit.append(meanFitTetra)
                if it >=100: #Stop if tetraploids have invaded for more than 10000 generations
                  flat_genoS=[]
                  for i in range(len(genotypes_S)):
                      flat_genoS += genotypes_S[i]
                  nbSalleles_end = len(set(flat_genoS))
                  phase += 1

            if generation%1000 == 0 and generation >= max_gen - 1000:
                print("generation =",generation)
                if generation>= max_gen:
                    notused2, meanFitTetra = get_meanFitnes(genotypes_S, phenotype,opt,om_2)
                    flat_genoS=[]
                    for i in range(len(genotypes_S)):
                        flat_genoS += genotypes_S[i]
                    nbSalleles_end = len(set(flat_genoS))
                    phase += 1 #Stops the simulation
                    
                    
                    
                    
        list_nbSalleles.append([nbSalleles, nbSalleles_end])
        list_freqSC.append(sum(freqSC)/len(freqSC)) #average over 10000 generations
        list_freqPollenSC.append(sum(freqPollenSC)/len(freqPollenSC)) #average over 10000 generations
        list_freqTetra.append(freqTetra)
        list_fitDiplo.append(meanFitDiplo)
        list_fitTetra.append(sum(meanFit)/len(meanFit))
        inDep.append(sum(list_ID)/len(list_ID))
    print("Done !")

    return list_nbSalleles, list_freqSC, list_freqPollenSC, list_freqTetra, inDep, list_fitTetra, list_fitDiplo





#Parameters
L = 50 #number of loci
varAddEff = 0.05 #extent of the mutation
dosage = 0.67
Npop = 1000 #population size
U = 0.005 #mutation rate
U_SI = 1e-5 #mutation rate at the S-locus
om_2 = 1 #width of the fitness function, i.e. strength of selection
ploidy = 2

alpha= 0.1 #selfing rate of the entire population
prob_unreduced = 0.05 #probability of producing unreduced gametes

nbSim = 5 #number of simulations





nbSalleles, freqSC, freqPollenSC, freqTetra, inDep, list_fitTetra, list_fitDiplo = simulation(nbSim, L, varAddEff, dosage, Npop, U, om_2 , ploidy, U_SI, prob_unreduced, alpha)

print(freqTetra)





