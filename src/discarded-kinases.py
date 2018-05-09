import pssms

AMINO_ACIDS = {"G": 0, "A": 1, "V": 2, "L": 3, "I": 4, "P": 5, "F": 6, "Y": 7, "W": 8, "S": 9, "T": 10, "C": 11, "M": 12, "N": 13, "Q": 14,"K": 15, "R": 16, "H": 17, "D": 18, "E": 19} 

if __name__ == "__main__":
    disc = "out/discarded-kinases.tsv"
    kept = "out/kept-kinases.tsv"


    d =pssms.read_kin_sub_data(disc)
    k=pssms.read_kin_sub_data(kept)

    daa = {}
    kaa = {}

    for kinase in d:
        noas = [0] *15 
        for i in range(len(d[kinase][0][3])):
            colaa = set()
            if i == 7:
                continue 	
            for j in range(len(d[kinase])):
             
                if d[kinase][j][3][i] not in colaa and d[kinase][j][3][i] in AMINO_ACIDS:
                    noas[i] = noas[i] +1 
                    colaa.add(d[kinase][j][3][i])
        if noas[i] != len(colaa):
            break
        daa[kinase] = (sum(noas)/14)/len(d[kinase])

    with open("out/discarded-kinases-colno.tsv","w") as f:
        for kinase in daa:
            f.write(kinase + "\t" + str(daa[kinase])+ "\n")
         

    for kinase in k:
        noas = [0] *15
        for i in range(len(k[kinase][0][3])):
            colaa = set()
            if i == 7:
                continue
            for j in range(len(k[kinase])):

                if k[kinase][j][3][i] not in colaa and k[kinase][j][3][i] in AMINO_ACIDS:
                    noas[i] = noas[i] +1
                    colaa.add(k[kinase][j][3][i])
        if noas[i] != len(colaa):
            break
        kaa[kinase] = (sum(noas)/14)/len(k[kinase])

    with open("out/kept-kinases-colno.tsv","w") as f:
        for kinase in kaa:
            f.write(kinase + "\t" + str(kaa[kinase])+ "\n")
 

   
           
