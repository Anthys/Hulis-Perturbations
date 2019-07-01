import numpy as np
import matplotlib as plt
import socket, select, string, sys, time, math

A = "Matrice des hamiltoniens (liaisons ou non)"

B = "Matrice des orbitales"

B_E = "Liste des energies des orbitales"

P = "Matrice des perturbations"

FRAGMENTS = {"L1":[], "L2":[]}

LA = "Liste des atomes"

O1 = "Orbi1, max E"

O2 = "Orbi2, min E"


#ALPHABETICAL ORDER FOR ATOM PAIRS
TABLE = {"CO2": 0.66, "CC":1, "CN2":0.89, "CN1":1.02, "CO1":1.06}

atoms_perturbed = []

visited = []

s = s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

# main function
def telnetmain():

    global s

    host = "127.0.0.1"
    port = 2345

    s.settimeout(2)

    # connect to remote host
    try:
        s.connect((host, port))
    except:
        print('Unable to connect')
        sys.exit()

    print('Connected to remote host')

    DICO = ["HAMILT", "ORBI_E", "ORBI_N", "ORBI_C", "ATOMS"]
    FINAL= {}
    indx = 0
    for i in ["HAMD", "ENED", "OCCD", "ORBD", "ATOM"]:
        msg = i
        s.send(msg.encode())
        time.sleep(0.2)
        data = ""
        socket_list = [sys.stdin, s]
        # Get the list sockets which are readable
        read_sockets, write_sockets, error_sockets = select.select(socket_list, [], [])
        for sock in read_sockets:
            data += s.recv(4096).decode()
            FINAL[DICO[indx]] = data.replace("EOF", "")
        indx+=1

    for i in FINAL:
        print(i, FINAL[i])

    return FINAL

def algorithme_parcours_de_graphs():
    global visited
    explor(0)
    visit1 = visited[:]
    visited=[]
    for i in range(len(ATOMS)):
        if ATOMS[i] not in visit1:
            explor(i)
            break
    if not visited:
        print("ONLY ONE FRAGMENT")
        clean_exit()
    for i in ATOMS:
        if i not in visited and i not in visit1:
            print("MORE THAN TWO FRAGMENTS")
            clean_exit()
    # visited = visit1[:]

def clean_exit():
    s.send("CLOSE".encode())
    sys.exit(0)

def explor(x):
    global visited
    visited.append(ATOMS[x])
    for j in range(len(HAMILT[x])):
        if HAMILT[x][j]>0.02 and not ATOMS[j] in visited:
            explor(j)

def get_perturb(i, j):
    i, j = "".join([a for a in i if not a.isdigit()]), "".join([a for a in j if not a.isdigit()])
    lis = sorted([i, j])
    str_fin = ''.join([i for i in lis])
    return TABLE[str_fin]

def find_orbitals():
    global FRAGMENTS

    first_find = ""
    sec = False
    for j in range(len(ORBI_E)-1, -1, -1):
        # Go through all orbitals, from the most instable to the most stable, to bet both Highest Occupied
        if ORBI_E[j] <= 0 and (not first_find) and ORBI_N[j]!=0:
            for i in range(len(ORBI_C)):
                if abs(ORBI_C[i][j])>0.002:
                    first_find ="L1" if ATOMS[i] in FRAGMENTS["L1"]["A"] else "L2"
                    # abs_max_col(j) is just a guess of the atom that will participate.
                    # If it is right, it saves a lot of time, if not, it will be corrected later
                    FRAGMENTS[first_find]["HO"] = [j, abs_max_col(j)]
                    break
        elif first_find:
            for i in range(len(ORBI_C)):
                if abs(ORBI_C[i][j]) > 0.002 and ATOMS[i] not in FRAGMENTS[first_find]["A"]:
                    # Verify that this is an orbital from a different fragment than the first
                    FRAGMENTS["L1" if ATOMS[i] in FRAGMENTS["L1"]["A"] else "L2"]["HO"] = [j, abs_max_col(j)]
                    sec = True
                    break
        if sec:
            break

    first_find = ""
    sec = False
    for j in range(len(ORBI_E)):
        # Go from most stable to less stable to get both Lowest Vacant
        if ORBI_E[j] >= 0 and (not first_find) and ORBI_N[j]==0:
            for i in range(len(ORBI_C)):
                if abs(ORBI_C[i][j]) > 0.002:
                    first_find = "L1" if ATOMS[i] in FRAGMENTS["L1"]["A"] else "L2"
                    FRAGMENTS[first_find]["BV"] = [j, abs_max_col(j)]
                    print(first_find, j)
                    break
        elif first_find:
            for i in range(len(ORBI_C)):
                if abs(ORBI_C[i][j]) > 0.002 and ATOMS[i] not in FRAGMENTS[first_find]["A"]:
                    FRAGMENTS["L1" if ATOMS[i] in FRAGMENTS["L1"]["A"] else "L2"]["BV"] = [j, abs_max_col(j)]
                    sec=True
                    break
        if sec:
            break


TAB = [["L1", "BV", "L2", "HO"], ["L1", "HO", "L2", "BV"]]

import sys

def get_atoms_participating():

    # Correct and get every atoms participating in the perturbations

    atoms_perturbed = []

    # if system is 2 atoms encountering 2 atoms (which may be the only case the program can't solve every time)
    if len(FRAGMENTS["L1"]["A"])==len(FRAGMENTS["L2"]["A"])==2:
        temp_i_l1 = 0
        temp_i_l2 = 0
        for yes in range(4):
            if ORBI_C[yes][FRAGMENTS["L1"]["BV"][0]]<0:
                temp_i_l1 = yes
                break
        for yes2 in range(4):
            if ORBI_C[yes2][FRAGMENTS["L2"]["BV"][0]]<0:
                temp_i_l2 = yes2
                break
        print(temp_i_l1, temp_i_l2)
        atoms_perturbed.append(ATOMS.index(FRAGMENTS["L1"]["A"][0] if FRAGMENTS["L1"]["A"][0]!=ATOMS[temp_i_l1] else FRAGMENTS["L1"]["A"][1]))
        atoms_perturbed.append(ATOMS.index(FRAGMENTS["L2"]["A"][0] if FRAGMENTS["L2"]["A"][0]==ATOMS[temp_i_l2] else FRAGMENTS["L2"]["A"][1]))
        atoms_perturbed.append(ATOMS.index(FRAGMENTS["L1"]["A"][0] if FRAGMENTS["L1"]["A"][0]==ATOMS[temp_i_l1] else FRAGMENTS["L1"]["A"][1]))
        atoms_perturbed.append(ATOMS.index(FRAGMENTS["L2"]["A"][0] if FRAGMENTS["L2"]["A"][0]!=ATOMS[temp_i_l2] else FRAGMENTS["L2"]["A"][1]))
        return atoms_perturbed




    for i in range(len(TAB)):
        temp_tab = FRAGMENTS[TAB[i][0]][TAB[i][1]]
        temp_tab2 = FRAGMENTS[TAB[i][2]][TAB[i][3]]
        c_a1, c_a2 = ORBI_C[temp_tab[1]][temp_tab[0]], ORBI_C[temp_tab2[1]][temp_tab2[0]]
        if sign(c_a1) == sign(c_a2) and not temp_tab2[1] in atoms_perturbed and not temp_tab[1] in atoms_perturbed:
            # if atoms already are correct
            atoms_perturbed.append(temp_tab[1])
            atoms_perturbed.append(temp_tab2[1])
        else:
            # if atoms aren't correct
            print("CORRECTING")
            done = False
            tempdic = np.copy(np_orbi[:,temp_tab[0]])
            tempdic2 = np.copy(np_orbi[:,temp_tab2[0]])
            while not done:
                # Awful loop correcting the atoms
                done2 = False
                while not done2:
                    print("YES", tempdic)

                    # Choose the maximum absolute value of coefficients in the orbitals of the two fragments
                    # and suppose it can pair it with another atom from the other fragment.
                    max_c = tempdic[np.argmax([abs(val) for val in tempdic])] if abs(tempdic[np.argmax([abs(val) for val in tempdic])]) > abs(tempdic[np.argmax([abs(val) for val in tempdic])]) else tempdic[np.argmax([abs(val) for val in tempdic])]

                    # Get corresponding informations
                    if max_c in tempdic:
                        max_i, max_L, float_nivel = list(tempdic).index(max_c), "L1", TAB[i][3]
                    else:
                        max_i, max_L, float_nivel = list(tempdic2).index(max_c), "L2", TAB[i][1]

                    # Check for problems, and if both maximum atoms are incorrect, set their coefficients to 0.
                    if max_i in atoms_perturbed:
                        #max_c = tempdic[np.argmax([abs(val) for val in tempdic])] if max_c == tempdic2[np.argmax([abs(val) for val in tempdic2])] else tempdic2[np.argmax([abs(val) for val in tempdic2])]
                        if max_c in tempdic:
                            tempdic[max_i] = 0
                            #max_i, max_L, float_nivel = list(tempdic).index(max_c), "L1", TAB[i][3]
                        else:
                            tempdic2[max_i] = 0
                            #max_i, max_L, float_nivel = list(tempdic2).index(max_c), "L2", TAB[i][1]
                    if max_i not in atoms_perturbed:
                        done2 = True
                    elif max_c==0:
                        # In case all atoms were tested, sys_exit
                        print("ERROR-1 FIRST IMPOSSIBLE")
                        sys.exit(0)
                    else:
                        # In case the top atom isn't right, restart the loop
                        print("ERROR-1 FIRST INADEQUATE")

                # One atom is now correct, research of its pair
                (float_frag, float_nivel) = "L2" if max_L=="L1" else "L2", float_nivel
                t_c=0
                t_ci=0
                print(float_frag, float_nivel)
                print(max_c)
                print("YESYES", np_orbi[:,FRAGMENTS[float_frag][float_nivel][0]])
                for c_i in range(len(np_orbi[:,FRAGMENTS[float_frag][float_nivel][0]])):
                    # Check for > or < based on sign of max_c
                    # Get the best pair available, that is in phase with the first atom
                    if (np_orbi[:,FRAGMENTS[float_frag][float_nivel][0]][c_i] > t_c if sign(max_c)==1 else np_orbi[:,FRAGMENTS[float_frag][float_nivel][0]][c_i] < t_c) and c_i not in atoms_perturbed:
                        t_c = np_orbi[:,FRAGMENTS[float_frag][float_nivel][0]][c_i]
                        t_ci=c_i

                # Fill the atoms_perturbed list in a certain order: L1 atoms then L2 ones
                if t_c == 0:
                    # In case every second atom has been tested
                    print("ERROR-2 SECOND IMPOSSIBLE")
                    #sys.exit(0)
                    print(max_i)
                    tempdic[max_i]=0
                elif max_L == "L1":
                    done = True
                    atoms_perturbed.append(max_i)
                    atoms_perturbed.append(t_ci)
                else:
                    done = True
                    atoms_perturbed.append(t_ci)
                    atoms_perturbed.append(max_i)
            print("On " + float_nivel + " " + float_frag + ", " + str(t_ci) + " -> " + str(max_i))
        print(atoms_perturbed, TAB[i])
    return atoms_perturbed




def sign(numb):
    return -1 if numb<0 else 1

def abs_max_col(col):
    return np.argmax(np.array([abs(i) for i in np_orbi[:, col]]))

def normalise(psi):
    squared=sum([i**2 for i in psi])
    squared = math.sqrt(1/squared)
    return [i*squared for i in psi]

def case_degenerate(Prtb, instr, energ):
    # From Orbitales frontières, Nguyên Trong Anh, InterEditions
    energ = -energ
    E1, E2 = energ + Prtb, energ - Prtb
    psi1 = [0.707*i for i in (sum_orb(FRAGMENTS[instr[0]][instr[1]][0], FRAGMENTS[instr[2]][instr[3]][0]))]
    psi2 = [0.707*i for i in (sum_orb(FRAGMENTS[instr[0]][instr[1]][0], FRAGMENTS[instr[2]][instr[3]][0]), -1)]
    return [psi1, psi2], [E1,E2]

def sum_orb(i1, i2, optional=1):
    new_orb = [0 for i in range(len(np_orbi[:,i1]))]
    for i in range(len(new_orb)):
        new_orb[i] = np_orbi[:,i1][i] + np_orbi[:,i2][i]*optional
    return new_orb

def case_non_degenerate(Prtb, instr):
    # From Orbitales frontières, Nguyên Trong Anh, InterEditions
    N=1

    Emin = min(ORBI_E[FRAGMENTS[instr[0]][instr[1]][0]],ORBI_E[FRAGMENTS[instr[2]][instr[3]][0]])
    Emax = max(ORBI_E[FRAGMENTS[instr[0]][instr[1]][0]],ORBI_E[FRAGMENTS[instr[2]][instr[3]][0]])

    Emin, Emax = -Emin, -Emax

    E1, E2 = Emin + Prtb**2/(Emin-Emax), Emax + Prtb**2/(Emax-Emin)

    psi1 = [N*i for i in (sum_orb(FRAGMENTS[instr[0]][instr[1]][0], FRAGMENTS[instr[2]][instr[3]][0], Prtb/(Emin-Emax)))]
    psi2 = [N*i for i in (sum_orb(FRAGMENTS[instr[2]][instr[3]][0], FRAGMENTS[instr[0]][instr[1]][0], Prtb/(Emax-Emin)))]

    return [normalise(psi1), normalise(psi2)], [E1,E2]


def r_atom(atom):
    while atom[0] in "0123456789":
        atom=atom[1:]
    return atom


P = [0, 0]
PSI = [[],[]]
E = [0, 0]

TAB = [["L1", "BV", "L2", "HO"], ["L1", "HO", "L2", "BV"]]

def find_perturbations():
    for i in range(len(TAB)):
        t_str = sorted([r_atom(ATOMS[FRAGMENTS[TAB[i][0]][TAB[i][1]][1]]), r_atom(ATOMS[FRAGMENTS[TAB[i][2]][TAB[i][3]][1]])])
        t_str = "".join([i for i in t_str])
        try:
            P[i] = TABLE[t_str]
        except Exception as e:
            print(e)
            print("ERROR, BINDING OF ATOMS MISSING IN TABLE")
            sys.exit(0)

        if abs(ORBI_E[FRAGMENTS[TAB[i][0]][TAB[i][1]][0]]-ORBI_E[FRAGMENTS[TAB[i][2]][TAB[i][3]][0]])<0.002:
            PSI[i], E[i] = case_degenerate(P[i], TAB[i], ORBI_E[FRAGMENTS[TAB[i][0]][TAB[i][1]][0]])
        else:
            PSI[i], E[i] = case_non_degenerate(P[i], TAB[i])

    return PSI, E

ORBI_C = [[0.84, 0, 0.54, 0],
        [0.71, 0, 0.84, 0],
        [0, 0.71, 0, 0.71],
        [0, 0.71, 0, 0.71]]
np_orbi = np.array(ORBI_C)

je_suis_sûr_de_moi = True

import json

if __name__ == "__main__":

    try:
        a = telnetmain()
    except Exception as e:
        print(e)
    HAMILT, ORBI_E, ORBI_N, ORBI_C = json.loads(a["HAMILT"]), json.loads(a["ORBI_E"]), json.loads(a["ORBI_N"]), json.loads(a["ORBI_C"])

    for i in range(len(ORBI_E)):
        if abs(ORBI_E[i])<0.02:
            ORBI_E[i]=0
        ORBI_E[i]=-ORBI_E[i]

    for i in range(len(HAMILT)):
        HAMILT[i][i]=0

    for i in range(len(ORBI_C)):
        for j in range(len(ORBI_C[i])):
            if abs(ORBI_C[i][j])<0.05:
                ORBI_C[i][j]=0

    if "ATOMS" in a:
        ATOMS = json.loads(a["ATOMS"])
        for i in range(len(ATOMS)):
            ATOMS[i] = str(i+1) + ATOMS[i]

    print("--")
    print(len(HAMILT), len(ORBI_C), len(ORBI_N), len(ORBI_E))
    print("ATOMS", ATOMS)
    for i in ORBI_C:
        print(i)

    FRAGMENTS = {"L1": {"A": [], "BV": 0, "HO": 0}, "L2": {"A": [], "BV": 0, "HO": 0}}
    #ATOMS = ["1C", "2C", "3C", "4C", "5O1"]
    algorithme_parcours_de_graphs()
    FRAGMENTS["L2"]["A"] = visited
    FRAGMENTS["L1"]["A"] = [i for i in ATOMS if i not in visited]
    np_orbi = np.array(ORBI_C)

    find_orbitals()
    print(FRAGMENTS)
    if je_suis_sûr_de_moi:
        # Correcting atoms
        tempor = get_atoms_participating()
        FRAGMENTS["L1"]["BV"][1]=tempor[0]
        FRAGMENTS["L2"]["HO"][1]=tempor[1]
        FRAGMENTS["L1"]["HO"][1]=tempor[2]
        FRAGMENTS["L2"]["BV"][1]=tempor[3]
    print(FRAGMENTS)
    # Calculations
    PSI, E = find_perturbations()
    print('----')
    print("PSI", PSI)
    print("E", E)
    #PSI, E = json.dumps(PSI), json.dumps(E)
    print('----')
    print("RSORBITALS", PSI)
    print("RSENERGIES", E)

    # PSI = [{"ATOMS":[{"ORB":2, "ATOM":"C2"}], "PRTRB_ORBS":[[[0,0,0,0],[0,0,0,0]], []], "PRTB_E":[[],[]]}]

    fintab = ["BV", "HO"]
    output = [0, 0]
    for i in range(2):
        output[i]= {"ATOMS": [{"ORB_INDX": FRAGMENTS["L1"][fintab[i]][0], "ATOM":ATOMS[FRAGMENTS["L1"][fintab[i]][1]]},
                              {"ORB_INDX": FRAGMENTS["L2"][fintab[(i+1)%2]][0], "ATOM":ATOMS[FRAGMENTS["L2"][fintab[(i+1)%2]][1]]}],
                    "PRTRB_ORBS": PSI[i],
                    "PRTRB_E": E[i]}
    print(output)

    output = json.dumps(output)
    print(len(output))

    msg = "RSPERTURBATION"
    s.send(msg.encode())
    s.send(output.encode())
    time.sleep(0.2)
    msg = "HOUR"
    s.send(msg.encode())
    time.sleep(0.2)
    s.send("CLOSE".encode())
    #FRAGMENTS={"L1":{"ATOMS INDEX COMPOSING FRAGMENT", "BV":[INDEX OF BV IN ORBI_E, INDEX OF ATOM PARTICIPATING IN ATOMS}}

sys.exit(0)


TAB = [["L1", "BV", "L2", "HO"], ["L1", "HO", "L2", "BV"]]

FRAGMENTS = {"L1": {"A":["4C", "3O2", "6C"], "BV":0, "HO":0}, "L2":{"A":["1C","2C", "5C"], "BV":0, "HO":0}}
ATOMS = ["1C","2C", "3O1", "4C", "5C", "6C"]
HAMILT = [[0, 1, 0, 0, 0, 0],
          [1, 0, 0, 0, 1, 0],
          [0, 0, 0, 1.06, 0, 0],
          [0, 0, 1.06, 0, 0, 1],
          [0, 1, 0, 0, 0, 0],
          [0, 0, 0, 1, 0, 0]]
algorithme_parcours_de_graphs()
FRAGMENTS["L2"]["A"]=visited
FRAGMENTS["L1"]["A"]= [i for i in ATOMS if i not in visited]
ORBI_E = [-1.84,  -1.41,   -0.41,  0,     1.28,   1.41]
ORBI_N = [2,       2,      2,     0,     0,       0]

ORBI_C = [[0,      0.5,    0,     -0.71,  0,     -0.50],
          [0,      0.71,   0,     0,     0,      0.71],
          [0.73,   0,      0.59,  0,    -0.35,   0],
          [0.60,   0,      0.31,  0,     0.74,   0],
          [0,      0.5,   -0.0,   0.71,  0,      0.5],
          [0.33,   0,     -0.75,  0,    -0.58,   0]]



FRAGMENTS = {"L1": {"A":[], "BV":0, "HO":0}, "L2":{"A":[], "BV":0, "HO":0}}
ATOMS = ["1C","2C", "3C", "4C"]
HAMILT = [[0, 1, 0, 0],
          [1, 0, 0, 0],
          [0, 0, 0, 1],
          [0, 0, 1, 0]]
algorithme_parcours_de_graphs()
FRAGMENTS["L2"]["A"]=visited
FRAGMENTS["L1"]["A"]= [i for i in ATOMS if i not in visited]
ORBI_E = [-1, -1, 1, 1]
ORBI_N = [2,       2,      0,       0]

ORBI_C = [[0.71,   0,      0,     0.71],
          [0.71,   0,      0,    -0.71],
          [0,      0.71,   0.71,  0],
          [0,      0.71,  -0.71,  0]]