# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 13:13:05 2025

@author: douet
"""

import numpy as np
import random as rd



#Creates the offspring, returns the list of ploidy for the next generation
def reproduction(Npop, prob_unreduced, selfRate, SCrate, list_ploidy):
    listOff = [] #contains the ploidy levels of the next generation

    
    #Creating the list of SC individuals
    list_tetra = []
    for i in range(Npop):
        if list_ploidy[i] == 4:
            list_tetra.append(i)
    list_SC = list_tetra
        
    
    Npop_off = 0
    
    while Npop_off != Npop:
        par1 = np.random.randint(0,Npop) #Selecting a maternal parent randomly
        ploidy_par1 = list_ploidy[par1]
        
        if np.random.rand() <= prob_unreduced : #unreduced gamete
            gamete_par1 = ploidy_par1 #number of chromosomes in the gamete of par1
        else: #reduces gamete
            gamete_par1 = ploidy_par1//2

        
        #OUTCROSSING
        if selfRate <= np.random.rand() or par1 not in list_SC: #The parent is either SC or its selfing rate is too low
            par2 = np.random.randint(0,Npop) #A second paternal parent is selected
            while par2==par1: 
                par2 = np.random.randint(0,Npop)
            ploidy_par2 = list_ploidy[par2]

            if np.random.rand() <= prob_unreduced : #unreduced gamete
                gamete_par2 = ploidy_par2      
            else: #reduced gamete
                gamete_par2 = ploidy_par2//2


            #Check if both gametes are compatible, else no new individual is added to the list of offspring
            if gamete_par1 == gamete_par2 and gamete_par1 <= 2: #Gametes of same size and at most diploid
                Npop_off += 1
                listOff.append(gamete_par1 + gamete_par2) #adding the ploidy of the new individual
                            
        #SELFING
        else:
            if np.random.rand() <= prob_unreduced : #unreduced gamete
                gamete2_par1 = ploidy_par1      
            else: #reduced gamete
                gamete2_par1 = ploidy_par1//2
                
            if gamete_par1 <= 2 and gamete_par1==gamete2_par1: #The offspring can be tetraploids at most
                Npop_off += 1
                listOff.append(gamete_par1 + gamete2_par1)       

    return listOff



#returns the frequency of tetraploids in the population
def frequencies(list_ploidy, Npop):
    nb_tetra = 0
    for n in range(Npop):
        #frequency of tetraploids in the population
        if list_ploidy[n] == 4:
            nb_tetra +=1
    freq_tetra =  nb_tetra/Npop
    return freq_tetra




def simulation(nbSim, Npop, prob_unreduced, selfRate, SCrate):
    Freq_tetra, Eq_tetra = [],[] #contain the frequency of tetraploids throughout and at the end of each simulation
    
    for k in range(nbSim):
        print(k)
        
        #initialization
        list_ploidy = [2 for i in range(Npop)] 
        
        nb_gen = 200000
        list_ft = []
        
        for g in range(nb_gen):
            list_ploidy = reproduction(Npop, prob_unreduced, selfRate, SCrate, list_ploidy) #Computing the next generation
            
            if g%100 == 0:
                print(g)
                freq_tetra = frequencies(list_ploidy, Npop)
                list_ft.append(freq_tetra)
        feq_tetra = sum(list_ft[-10:])/len(list_ft[-10:]) #average of the frequency of tetraploids over the last 1000 generations
        
        Freq_tetra.append(list_ft)
        Eq_tetra.append(feq_tetra)
            
    
    print("Done !")
        
    return Freq_tetra, Eq_tetra





#Parameters
Npop = 1000 #population size 
prob_unreduced = 0.05 #probability of producing unreduced gametes
SCrate = 1
selfRate = 0.3


nbSim = 30 #number of simulations



list_freq_tetra, feq_tetra = simulation(nbSim, Npop, prob_unreduced, selfRate, SCrate)






