from numpy import *
temp = 37
kT = (temp + 273.15) * 8.314462175# kT in Joules/mol
kcal = 4.184E3 # kcal in Joules
Na = 6.0221412927E23 # Avogadro number
V = 4./3 * pi * 7.5E-6**3
# concentration to copy number
c0 = Na*V*1000

#print "cell volume of 15um diameter %.2f um^3" % (V*1e6**3)

def cell_conc(n):
    return n/V/Na/1000 # 1000L = 1 m^3

def kcal_to_Kd(E,temp=37):
    kT = (temp + 273.15) * 8.314462175# kT in Joules/mol
    
    # kcal/mol to Kd in nM
    return exp(E*kcal/kT)*1e9

def Kd_to_kcal(K,temp=37):
    kT = (temp + 273.15) * 8.314462175# kT in Joules/mol
    
    # Kd in nM to E in kcal/mol
    return log(K/1E9)*kT/kcal
  
def zamore_energy(struct):
    """
    accepts RNAhybrid 2ndary structure of target&mirna .
    returns binding energy in kcal/mol based on base-pairing
    in the different segments of the guide.
    """
    contrib = array(list(struct.split('&')[1])) != "."
    #for s,e in zip(struct.split('&')[1],mir_pairing_energy):
        #print s,e
    return (E_AGO_mouse[:len(contrib)]*contrib).sum()

mir_pairing_energy =  zeros(26,dtype=double)    

def setup_mir_pairing_energies():
    # energy model inspired by Wee et al. 2012 (Cell, Zamore lab)
    global mir_pairing_energy

    # seed match (2:9) only
    s_seed7 = ".))))))).............."
    E_seed7 = Kd_to_kcal(.026,temp=25)
    # contribution of a match inside the seed
    seed_match = E_seed7/7.
    mir_pairing_energy[1:8] = seed_match
    
    # mouse Ago2 with let-7 pairing to 7-mer seed 
    # w/ additional 3'pairing
    s_supp = ".))))))).....))))....."
    E_supp = Kd_to_kcal(.013,temp=25)
    # contribution of a match in g13-g16
    supp_pairing = (E_supp-E_seed7)/4.
    mir_pairing_energy[12:17] = supp_pairing
   
    # assign to pos. g9 quarter the binding energy of seed
    # and pos g1,g10 the supp. 3' pairing quantity
    # this is ad-hoc but probably OK
    mir_pairing_energy[8] = .25*seed_match
    mir_pairing_energy[0] = supp_pairing
    mir_pairing_energy[9] = supp_pairing

    # perfect 21 nt is actually bound *less* stable.
    s_perfect = "))))))))))))))))))))))"
    E_perfect = Kd_to_kcal(.020,temp=25)

    # assign the apparent energy penalty to g17-g21 (and to the rest)
    extensive_pairing_penalty = (E_perfect-zamore_energy('&'+s_perfect))/5
    #mir_pairing_energy[8:13] = penalty
    mir_pairing_energy[17:] = extensive_pairing_penalty

E_AGO_mouse = []
E_AGO_fly = []

def setup_mir_pairing_energies2():
    """ constructed from the supplement of Wee et al. 2012 (Zamore lab) """
    global E_AGO_mouse,E_AGO_fly
    pos = [
        'g1','g1-g2','g2-g3','g3-g4','g4-g5','g5-g6','g6-g7','g7-g8',
        'g8','g9','g8-g9','g9-g10','g10-g11','g12-g13','g13-g14','g15','g16_1','g16_2','g16_3','g14-g15','g15-g16_1','g15-g16_2','g17','g16-g17','g17-g18','g18-g19','g19-g20','g20-g21','g17-g21']

    Kmm = array([17,42,250,320,2400,1100,390,88,74,30,6.0,17,8.3,34,68,12,13,14,7.0,36,20,20,23,22,19,33,19,29,15])
    Kperfect = array([17,26,23,31,32,51,15,24,49,16,15,20,8,20,19,7,4.3,4.3,3.4,10,6.7,5.4,12,6.9,12,22,13,16,11])

    # log-ratios of mismatched and matched Km should correspond to the energy differences
    dE = -log(Kmm/Kperfect)

    Ep = {}
    for p,e,mm in zip(pos,dE,Kmm):
        Ep[p] = e
        #print p,mm,e

    # build a pseudo-energy model per position by averaging the available experimental data
    E = zeros(25)
    E[0] = (Ep['g1'])# + Ep['g1-g2']/2.)/2.
    E[1] = (Ep['g1-g2'] + Ep['g2-g3'])/4.
    E[2] = (Ep['g2-g3'] + Ep['g3-g4'])/4.
    E[3] = (Ep['g3-g4'] + Ep['g4-g5'])/4.
    E[4] = (Ep['g4-g5'] + Ep['g5-g6'])/4.
    E[5] = (Ep['g5-g6'] + Ep['g6-g7'])/4.
    E[6] = (Ep['g6-g7'] + Ep['g7-g8'])/4.
    E[7] = (Ep['g7-g8']/2. + Ep['g8-g9']/2. + Ep['g8'])/3.
    E[8] = (Ep['g8-g9']/2. + Ep['g9-g10']/2. + Ep['g9'])/3.
    E[9] = (Ep['g9-g10'] + Ep['g10-g11'])/4.
    E[10] = Ep['g10-g11'] -E[9]
    E[11] = Ep['g12-g13']/2.
    E[12] = (Ep['g12-g13'] + Ep['g13-g14'])/4.
    E[13] = (Ep['g13-g14'] + Ep['g14-g15'])/4.
    E[14] = (Ep['g14-g15']/2. + Ep['g15-g16_1']/2. + Ep['g15-g16_2']/2. + Ep['g15'])/4.
    E[15] = (Ep['g15-g16_1']/2. +Ep['g15-g16_2']/2. + Ep['g16-g17']/2. + Ep['g16_1'] + Ep['g16_2'] + Ep['g16_3'])/6.
    E[16] = (Ep['g17'] + Ep['g16-g17']/2. + Ep['g17-g18']/2. + Ep['g17-g21']/5.)/4.
    E[17] = (Ep['g17-g18']/2. + Ep['g18-g19']/2. + Ep['g17-g21']/5.)/3.
    E[18] = (Ep['g18-g19']/2. + Ep['g19-g20']/2. + Ep['g17-g21']/5.)/3.
    E[19] = (Ep['g19-g20']/2. + Ep['g20-g21']/2. + Ep['g17-g21']/5.)/3.
    E[20] = (Ep['g20-g21']/2. + Ep['g17-g21']/5.)/2.
    #print E.sum()
    #print E

    # the best fly-Ago Kd is 3.7 pM for a 21nt perfect match.
    # Use that energy to scale the pseudo-energies
    E_fly = Kd_to_kcal(0.0037,temp=25)
    #print "E_fly",E_fly
    E_AGO_fly = E*(E_fly/E.sum())
    
    # now try to estimate a model for mouse (mammals)
    # using the three measured duplexes
    #print "ENERGY MODEL"
    #".))))))).............."
    E_seed = Kd_to_kcal(0.026,temp=25)
    #".))))))).....))))....."
    E_best = Kd_to_kcal(0.013,temp=25)
    #"))))))))))))))))))))))"
    E_perfect = Kd_to_kcal(0.020,temp=25)

    deltaE = E_best - E_seed
    #print "deltaE",deltaE

    # the three guide segments: Seed, supplementary 3' pairing, 
    # and everything else, which in mouse was observed to *reduce* affinity!
    C1 = E_AGO_fly[1:8].sum()
    C2 = E_AGO_fly[13:17].sum()
    C3 = E_AGO_fly.sum()-C1-C2-E_AGO_fly[8]

    #print "delta_perf",delta_perf
    s1 = E_seed/C1
    s2 = deltaE/C2
    delta_perf = E_perfect - E_best-E_AGO_fly[8]*s1
    s3 = delta_perf/C3

    #print "scaling factors",s1,s2,s3
    # rescale the energy model in each segment to get 
    # an approximate model for mouse Ago
    E_AGO_mouse = array(E_AGO_fly)
    E_AGO_mouse[1:9] = s1*E_AGO_fly[1:9]
    E_AGO_mouse[13:17] = s2*E_AGO_fly[13:17]
    E_AGO_mouse[9:13] = s3*E_AGO_fly[9:13]
    E_AGO_mouse[17:] = s3*E_AGO_fly[17:]
    
setup_mir_pairing_energies2()

if __name__=="__main__":
    #print cell_conc(1)

    
    print "1000 molecules are %.3f nM" % (cell_conc(1000)/1E-9)
    print "1nM is %.1f molecules per cell" % (1e-9*Na*V*1000)
    print "concentration of all binding sites in 3'UTRs in uM", cell_conc(1E8)*1e6
    print "Kd of -15 kcal/mol in nM", kcal_to_Kd(-15)
    print "Kd of -17 kcal/mol in nM", kcal_to_Kd(-17)
    print "Kd of -17 kcal/mol in nM", exp(-17*kcal/kT)*1e9
    print "binding energy corr. to 10nM Kd", log(20E-9)/kcal*kT
    print "1 kcal/mol in kT @37 C", kcal,kT,kcal/kT*11
    #print Kd_to_kcal(kcal_to_Kd(-15))

    # units in nM
    zamore_mouse_21 = .020 #
    zamore_mouse_7 = .026 # seed7 only
    zamore_mouse_11 = .013 # seed7 + g13-16

    #print Kd_to_kcal(zamore_mouse_21)
    #print Kd_to_kcal(zamore_mouse_7)
    #print Kd_to_kcal(zamore_mouse_11)
    #print Kd_to_kcal(0.004)
        
    print "testing energy model"
    for e,s in [(Kd_to_kcal(.013,temp=25),"&.))))))).....))))....."),(Kd_to_kcal(.026,temp=25),"&.))))))).............."),(Kd_to_kcal(.020,temp=25),"&))))))))))))))))))))))")]:
        print s, "should have ",e
        print zamore_energy(s)

    #print cell_conc(1E6)
