#!/bin/env python3

#################################################################################################################################
#Script Name    : EvoSim.py
#Description    : simulate evolution along a tree using a model dependent on protein stability and population genetics
#Author         : Michael A Sennett
#################################################################################################################################

import gzip
import numpy as np
import random as rd
import pypdb as pp
import dendropy as dp
import os
import argparse
import matplotlib.pyplot as plt

        
#################################################################################################################################

def parse_pdb_file(pdb):
    
    #grab some pdb file
    pdb_file=pp.get_pdb_file(str(pdb))
    pdb_lines=pdb_file.split('\n')
        
    #parse through the pdb file
    tmp_pdb_lines=[]
    for line in pdb_lines:
        if line.startswith('ATOM'):
            tmp_pdb_lines.append(line[0:5]) #'ATOM '
            tmp_pdb_lines.append(line[6:11]) #atom serial num
            tmp_pdb_lines.append(line[13:15]) #element on residue
            tmp_pdb_lines.append(line[17:20]) #residue
            tmp_pdb_lines.append(line[21]) #chain
            tmp_pdb_lines.append(line[22:26]) #residue integer
            tmp_pdb_lines.append(line[30:37]) #x coord
            tmp_pdb_lines.append(line[38:45]) #y coord
            tmp_pdb_lines.append(line[46:53]) #z coord
        else:
            pass
    
    pdb_2d_arr=np.reshape(tmp_pdb_lines, (-1,9))
    pdb_2d_arr=np.char.strip(pdb_2d_arr)
    
    
    #find the number of chains in the pdb file
    chain = []
    chain_len = []
    for i in pdb_2d_arr[:,4]:
        if i not in chain:
            chain.append(i)
    pdb_2d_arr_A=[]        
    for line in pdb_2d_arr:
        if line[4] != chain[0] or line[2] != 'CA':
            pass
        else:
            pdb_2d_arr_A.append(list(line))
    
    pdb_2d_arr_A=np.reshape(pdb_2d_arr_A, (-1,9))
    pdb_2d_arr_A=np.char.strip(pdb_2d_arr_A)
    
    #determine the length of each chain
    for j in range(len(chain)):
        res_num=[]
        for k in range(len(pdb_2d_arr)):
            if pdb_2d_arr[k][4] == chain[j]:
                res_num.append(int(pdb_2d_arr[k][5]))
        chain_len.append(max(res_num))    
        
    #determine the sequence and corresponding residue number of the first chain in each pdb file
    wt_chain=chain[0]
    res_nums=[]
    wt_seq=[]
    chain_A=[]
    for l in range(len(pdb_2d_arr)):
        if pdb_2d_arr[l,4] != wt_chain:
            pass
        else:
            if pdb_2d_arr[l,5] not in res_nums:
                res_nums.append(pdb_2d_arr[l,5])
                wt_seq.append(pdb_2d_arr[l,3])

    
    AA_3=['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    AA_1=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    
    AA_dict=dict(zip(AA_3, AA_1))
    wt_seq_1=[]
    for aa in wt_seq:
        wt_seq_1.append(AA_dict.get(aa))
    
    return wt_seq_1, res_nums, pdb_2d_arr_A

#################################################################################################################################

#create a random sequence
def gen_rd_seq(seq_len):
    
    rd_seq=list()
    for res in range(len(seq_len)):
        rd_seq.append(rd.sample(AA_miya, 1)[0])

    return rd_seq

#################################################################################################################################

#grab list of sites that are near residue of interest in a single chain
def nn_sites(pdb_2d_arr):

    pdb_2d_arr=pdb_2d_arr
    
    #create a look-up table for each site in a sequence
    site=[]
    nn_site=[]
    for line in pdb_2d_arr:
        site_coord=np.zeros((1,3))
    #get the ref residue coord
        if line[5] not in site:
            site.append(line[5])
            site_coord[0,0]=float(line[6])
            site_coord[0,1]=float(line[7])
            site_coord[0,2]=float(line[8])

    #check to see if next residue is close
            tmp_nn_sites=[]    
            for new_line in pdb_2d_arr:
                nn_coord=np.zeros((1,3))
                if new_line[5] != line[5]:

                    nn_coord[0,0]=float(new_line[6])
                    nn_coord[0,1]=float(new_line[7])
                    nn_coord[0,2]=float(new_line[8])
    #calc dist bw sites
                    dist=np.sum((site_coord-nn_coord)**2,1)**0.5
                    if dist < 7:
                        tmp_nn_sites.append(new_line[5])
                    else:
                        pass
                else:
                    pass
            nn_site.append(tmp_nn_sites)
        else:
            pass

    lookup_table=dict(zip(site,nn_site))
    
    return lookup_table            

#################################################################################################################################

def trav_tree(sim_tree, all_tax):
    
    N=0
    tree1=dp.Tree(sim_tree)
    taxa=[]
    for node in sim_tree.preorder_node_iter():
        
        if node.is_internal():
            
            if node.taxon is None:
                x=dp.datamodel.taxonmodel.Taxon(label="N"+str(N))
                node.taxon=x
                all_tax.add_taxon(node.taxon)
                N+=1
            
            taxa.append(node.taxon.label)
        
        else:
            taxa.append(node.taxon.label)
    
    edges=[]
    for edge in sim_tree.preorder_edge_iter():
        
        if edge.length == None:
            edges.append(0)
        
        else:
            edges.append(edge.length)
    
    return taxa, edges
    
#################################################################################################################################

def get_DNA_CDS(prot_seq):
    
    cds_seq=['X']*len(prot_seq)
    i=0
    for aa in prot_seq:
        cdn=CODON[aa]
        tmp_cdn=rd.sample(cdn, 1)
        cds_seq[i]=tmp_cdn[0]
        i+=1
    
    return cds_seq

def get_cdn_key(val):
    for keys, values in CODON.items():
        for value in values:
            if val == value:
                return keys
        

#################################################################################################################################

def calc_ddG(tmp_seq, lookup_table, res_nums, AA_miya, seq):

    #calculate the ddG of a mutation
    
    res_dct=dict(zip(res_nums, range(len(res_nums))))
    site_mut=rd.sample(range(1,len(tmp_seq)),1)
    
    if seq == 2:
        cdn_seq=tmp_seq
        #wt codon identity, 3L
        wt_cd=cdn_seq[site_mut[0]]
        #pos in wt codon to be mutated
        pos=rd.sample(range(3), 1)
        cd_pos=pos[0]
        #sample nuc to mutate
        mut_num_lst=rd.sample(range(len(DNA)), 4)
        #make sure nuc actually would change
        for mut_num in mut_num_lst:
            #change nuc in codon
            mut_cd=list(wt_cd)
            mut_cd[cd_pos]=DNA[mut_num]
            mut_cd=''.join(mut_cd)
            #if no mut or mut to stop then pass, else loop is over
            if DNA[mut_num] == wt_cd[cd_pos]:
                pass
            elif get_cdn_key(mut_cd) == 'stop':
                pass
            else:
                break
        
        var_res=get_cdn_key(mut_cd)
        wt_res=get_cdn_key(wt_cd)
        
        #create an aa seq for later
        tmp_seq=[]
        for cdn in cdn_seq:
            aa = get_cdn_key(cdn)
            tmp_seq.append(aa)
        
                
    if seq < 2:    
        wt_res=tmp_seq[site_mut[0]]
        mut_num_lst=rd.sample(range(len(AA_miya)), 20)

        #intro mutation according to poisson process for aa
        for mut_num in mut_num_lst:
            if AA_miya[mut_num] == wt_res:
                pass
            else:
                break

        var_res=AA_miya[mut_num]
        #dummy codon
        mut_cd='NNN'
    
    struc_site_mut=res_nums[site_mut[0]]
    nn_res_3d=lookup_table.get(struc_site_mut)

    
    cp_wt=[]
    cp_var=[]
    
    for res in nn_res_3d:
        prim_site=res_dct.get(res)
        prim_aa=tmp_seq[prim_site]
        nn_id=AA_miya_dct.get(prim_aa)
        wt_id=AA_miya_dct.get(wt_res)
        var_id=AA_miya_dct.get(var_res)
        cp_wt.append(miyazawa[wt_id][nn_id])
        cp_var.append(miyazawa[var_id][nn_id])

    ddG=sum(cp_var)-sum(cp_wt)

    return ddG, wt_id, var_id, mut_num, site_mut[0], mut_cd
        
#################################################################################################################################

def calc_Pfix(dG_i, ddG, N_eff, fit):
    
    #calculate the probability of fixation  
    if fit == 2:
        mu=1
        sigma=1.4
        ddG=rd.gauss(mu, sigma)
        
        dG_V=dG_i+ddG
        wv=1/(1+np.exp(dG_V))
        wi=1/(1+np.exp(dG_i))

        s_mut=np.log(wv)-np.log(wi)

        if N_eff > 10:
            scale=N_eff/30
            fv=np.expm1(-s_mut)/np.expm1(-1*N_eff*s_mut)*scale

        else:
            fv=np.expm1(-s_mut)/np.expm1(-1*N_eff*s_mut)

        return fv, dG_V, s_mut, ddG

    
    if fit == 1:
        dG_V=dG_i+ddG
        wv=1/(1+np.exp(dG_V))
        wi=1/(1+np.exp(dG_i))

        s_mut=np.log(wv)-np.log(wi)

        if N_eff > 10:
            scale=N_eff/30
            fv=np.expm1(-s_mut)/np.expm1(-1*N_eff*s_mut)*scale

        else:
            fv=np.expm1(-s_mut)/np.expm1(-1*N_eff*s_mut)

        return fv, dG_V, s_mut, ddG
    
    else:
        
        wi=1/(1+np.exp(dG_i))
        c=10/N_eff
        dw=rd.uniform(-0.5,0.5)*c
        wv=wi+dw
        if wv >= 1.0:
            wv=2.0-wv
            
        s_mut=np.log(wv)-np.log(wi)
        if N_eff > 10:
            scale=N_eff/100
            fv=np.expm1(-s_mut)/np.expm1(-1*N_eff*s_mut)*scale
        else:
            fv=np.expm1(-s_mut)/np.expm1(-1*N_eff*s_mut)

        dG_V=np.log(1-wv)-np.log(wv)
        ddG=dG_V-dG_i

        return fv, dG_V, s_mut, ddG

#################################################################################################################################

def EvoSim(wt_seq, res_nums, pdb_2d_arr, tree, Neff, seq, fit):

    #user defined initial variables from input
    N_eff=Neff
    fit=fit
    seq=seq
    
    
    #non-user defined initial variables from input
    res_dct=dict(zip(res_nums, range(len(res_nums))))
    seq_len=len(res_nums)
    dG_WT=(-1*np.log(N_eff)-0.57)/1
    lookup_table=nn_sites(pdb_2d_arr)

    #sample a random number and then throw into seed so we can track
    seed=rd.randint(0,100000000)
    print("random seed number:", seed)
    rd.seed(seed)
    
    #create a random sequence if you want to
    if seq == 0:
        rd_seq=gen_rd_seq(seq_len)
        wt_seq=rd_seq

    #create a coding sequence if you want to
    elif seq == 2:
        cd_seq=get_DNA_CDS(wt_seq)
        wt_seq=cd_seq

    #get the tree to sim evolution along
    sim_tree=dp.Tree.get(path=tree,
                        schema="newick",
                        preserve_underscores=True)

    sim_tree=sim_tree.extract_tree()
    all_tax=dp.TaxonNamespace(label="taxa")
    all_tax.add_taxa(sim_tree.taxon_namespace.labels())    

    #traverse the tree and name internal nodes that are unnamed
    taxa, edges=trav_tree(sim_tree, all_tax)

    #initialize storage lists to track seqs and values
    aln={}
    seq_dGs={}
    af_probs=[]
    af_prob_of_fix=[]
    af_slxns=[]
    af_dGs=[]
    af_ddGs=[]
    mut_mat=np.zeros((20,20))
    bf_dGs=[]
    bf_ddG=[]

    #traverse the tree and sim evo
    dG_i=dG_WT

    for node in sim_tree.preorder_node_iter():

        if node.taxon.label not in aln:

            aln[node.taxon.label]="".join(wt_seq)
            seq_dGs[node.taxon.label]=dG_i
            
        x=node.child_nodes()

        #for each of the child nodes determine the numb of mutations expected
        for i in range(len(x)):
            
            exp_mut=int(x[i].edge_length*len(wt_seq))
            tmp_dGs=[]

            if seq == 2:
                tmp_seq=aln.get(node.taxon.label)
                new_tmp=[]
                for j in range(0, len(tmp_seq), 3):
                    new_tmp.append(tmp_seq[j:j+3])
                tmp_seq=new_tmp
                
            else:
                tmp_seq=list(aln.get(node.taxon.label))

            dG_i=seq_dGs.get(node.taxon.label)

            if x[i].taxon.label not in aln:

                print(x[i].taxon.label, ", add as child and mutate")

                m=0

                if exp_mut == 0:
                    exp_mut=1
                
                while m < exp_mut:
                    
                    ddG, wt_id, var_id, aa_mut_num, site_mut, mut_cd=calc_ddG(tmp_seq, lookup_table, res_nums, AA_miya, seq)
                    fv, dG_V, s_mut, ddG=calc_Pfix(dG_i, ddG, N_eff, fit)

                    bf_dGs.append(dG_V)
                    bf_ddG.append(ddG)
                    dice=rd.uniform(0,1)
                    if dice < fv:
                        mut_mat[wt_id][var_id]+=1
                        dG_i=dG_V
                        #if codon seq, sub codon, else sub in amino acid
                        if seq == 2:
                            tmp_seq[site_mut]=mut_cd
                        else:
                            tmp_seq[site_mut]=AA_miya[aa_mut_num] 
                        af_probs.append(dice)
                        af_prob_of_fix.append(fv)
                        af_slxns.append(s_mut)
                        tmp_dGs.append(dG_i)
                        af_ddGs.append(ddG)
                        m+=1
                    else:
                        dG_i=dG_i
                        tmp_dGs.append(dG_i)
                aln[x[i].taxon.label]="".join(tmp_seq)

            else:
                pass

            seq_dGs[x[i].taxon.label]=tmp_dGs[-1]
            af_dGs.append(tmp_dGs)
           
    out_tree=sim_tree.as_string('newick', unquoted_underscores=True)

    return m, aln, af_slxns, af_probs, af_prob_of_fix, af_dGs, af_ddGs, mut_mat, bf_dGs, bf_ddG, out_tree

#################################################################################################################################
#write the tree as a text file
def write_tree(tree, Neff, pdb, direct, add):
    
    title=pdb+'_Neff_'+str(Neff)+'_outtree_'+add+'.txt'
    with open(os.path.join(direct,title), 'w+') as outfile:
        outfile.write(tree)
        outfile.close()

#write two alignments, one for the ancestors and one for the extant sequences

def write_aln(aln, Neff, pdb, direct, add):

    title=pdb+'_Neff_'+str(Neff)+'_extants_'+add+'.txt'
    with open(os.path.join(direct,title), 'w+') as outfile:
        for item in aln:
            if not item.startswith("N"):
                outfile.write(">"+item+"\n")
                outfile.write(aln.get(item)+"\n")
            else:
                pass
        outfile.close()

    title=pdb+'_Neff_'+str(Neff)+'_ancs_'+add+'.txt'
    with open(os.path.join(direct,title), 'w+') as outfile:
        for item in aln:
            if item.startswith("N"):
                outfile.write(">"+item+"\n")
                outfile.write(aln.get(item)+"\n")
            else:
                pass
        outfile.close()

#################################################################################################################################
#plt the slxn coefficients, sampled probs, and scaled probabilities of fixation
        
def plt_slxn_pfix(s, p, pfix, direct, Neff, pdb, add):

    #plt.rc('font', family='Arial')
    plt.rcParams.update({'font.size': 16})

    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    plt.subplots_adjust(wspace=0.4)

    #plt the data
    count_s=ax[0].hist(s, int(np.sqrt(len(s))))
    count_p=ax[1].hist(p, int(np.sqrt(len(p))))
    count_f=ax[2].hist(pfix, int(np.sqrt(len(pfix))))

    #set the axis labels
    ax[0].set_xlabel("selection coefficients")
    ax[1].set_xlabel("sampled probabilities")
    ax[2].set_xlabel("fixation probabilities")
    ax[0].set_ylabel("counts")

    #include relevant stats
    ax[0].text(max(count_s[1])*0.25, max(count_s[0])*.8, 
               "Avg S: "+str(round(np.average(s),8))+"\n"
               +"Std S: "+str(round(np.std(s), 8))+"\n"
               +"Substitutions: "+str(len(s)),
               fontsize=14
              )
    ax[1].text(max(count_p[1])*0.25, max(count_p[0])*.8, 
               "Avg S: "+str(round(np.average(p),8))+"\n"
               +"Std S: "+str(round(np.std(p), 8))+"\n",
               fontsize=14
              )
    ax[2].text(max(count_f[1])*0.25, max(count_f[0])*.8, 
               "Avg S: "+str(round(np.average(pfix),8))+"\n"
               +"Std S: "+str(round(np.std(pfix), 8))+"\n"
               +"Max Prob: "+str(round(max(pfix),3)),
               fontsize=14
              )

    title=pdb+'_Neff_'+str(Neff)+'_slxn_probs_fix_'+add+'.png'
    plt.savefig(os.path.join(direct,title), bbox_inches='tight', pad_inches=0.3, dpi=300) 
    
    title2=pdb+'_Neff_'+str(Neff)+'_slxn_probs_fix_'+add+'.txt'
    with open(os.path.join(direct,title2), 'w+') as outfile:
        outfile.write('{0}            {1}           {2}\n'.format("slxn","probs","pfix"))
        for i in range(len(s)):
            outfile.write("%4.4E\t%4.4E\t%4.4E\n"%(s[i],p[i],pfix[i]))
        outfile.close()

#################################################################################################################################
#plt the dG of protein mutants before selection occurs
    
def plt_dG_b4_slxn(pdb, Neff, b4_dG, direct, add):
    
    #plt.rc('font', family='Arial')
    plt.rcParams.update({'font.size': 14})

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    plt.subplots_adjust(wspace=0.3)

    count_dG=plt.hist(b4_dG, int(np.sqrt(len(b4_dG))))
    plt.xlabel('dG mut')
    plt.ylabel('count')

    plt.text(max(count_dG[1])*0.25, max(count_dG[0])*.8, 
               "Avg dG: "+str(round(np.average(b4_dG),4))+"\n"
               +"Std dG: "+str(round(np.std(b4_dG), 4))+"\n"
               +"Tried dGs: "+str(len(b4_dG)),
               fontsize=14
              )

    title=pdb+'_Neff_'+str(Neff)+'_b4_dGs_'+add+'.png'
    plt.savefig(os.path.join(direct,title), bbox_inches='tight', pad_inches=0.3, dpi=300)


#################################################################################################################################
#plt the dG of protein mutants after selection occurs

def plt_dG_af_slxn(pdb, Neff, dGs, direct, add):
    
    #plt.rc('font', family='Arial')
    plt.rcParams.update({'font.size': 14})

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    plt.subplots_adjust(wspace=0.3)

    all_dGs=[]
    for i in range(len(dGs)):
        for j in range(len(dGs[i])):
            all_dGs.append(dGs[i][j])

    count_dG=plt.hist(all_dGs, int(np.sqrt(len(all_dGs))/10))
    plt.xlabel('dG mut')
    plt.ylabel('count')

    plt.text(max(count_dG[1])*0.25, max(count_dG[0])*.8, 
               "Avg dG: "+str(round(np.average(all_dGs),4))+"\n"
               +"Std dG: "+str(round(np.std(all_dGs), 4))+"\n"
               +"Accepted dGs: "+str(len(all_dGs)),
               fontsize=14
              )

    title=pdb+'_Neff_'+str(Neff)+'_af_dGs_'+add+'.png'
    plt.savefig(os.path.join(direct,title), bbox_inches='tight', pad_inches=0.3, dpi=300)

################################################################################################################################# 

def plt_ddG_b4_slxn(pdb, Neff, b4_ddG, direct, add):

    #plt.rc('font', family='Arial')
    plt.rcParams.update({'font.size': 14})

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    plt.subplots_adjust(wspace=0.3)

    count_ddG=plt.hist(b4_ddG, int(np.sqrt(len(b4_ddG))))
    plt.xlabel('ddG mut')
    plt.ylabel('count')

    plt.text(max(count_ddG[1])*0.25, max(count_ddG[0])*.8, 
               "Avg ddG: "+str(round(np.average(b4_ddG),4))+"\n"
               +"Std ddG: "+str(round(np.std(b4_ddG), 4))+"\n"
               +"Tried ddGs: "+str(len(b4_ddG)),
               fontsize=14
              )

    title=pdb+'_Neff_'+str(Neff)+'_b4_ddGs_'+add+'.png'
    plt.savefig(os.path.join(direct,title), bbox_inches='tight', pad_inches=0.3, dpi=300)

#################################################################################################################################

def plt_ddG_af_slxn(pdb, Neff, af_ddG, direct, add):

    #plt.rc('font', family='Arial')
    plt.rcParams.update({'font.size': 14})

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    plt.subplots_adjust(wspace=0.3)

    count_ddG=plt.hist(af_ddG, int(np.sqrt(len(af_ddG))))
    plt.xlabel('ddG mut')
    plt.ylabel('count')

    plt.text(max(count_ddG[1])*0.25, max(count_ddG[0])*.8, 
               "Avg ddG: "+str(round(np.average(af_ddG),4))+"\n"
               +"Std ddG: "+str(round(np.std(af_ddG), 4))+"\n"
               +"Accepted ddGs: "+str(len(af_ddG)),
               fontsize=14
              )

    title=pdb+'_Neff_'+str(Neff)+'_af_ddGs_'+add+'.png'
    plt.savefig(os.path.join(direct,title), bbox_inches='tight', pad_inches=0.3, dpi=300)

#################################################################################################################################

#contact energies

#################################################################################################################################

AA_miya = ['C', 'M', 'F', 'I', 'L', 'V', 'W', 'Y', 'A', 'G', 'T', 'S', 'Q', 'N', 'E', 'D', 'H', 'R', 'K', 'P']

AA_miya_dct=dict(zip(AA_miya, range(len(AA_miya))))

miyazawa = [[-1.06,0.19,-0.23,0.16,-0.08,0.06,0.08,0.04,0.00,-0.08,0.19,-0.02,0.05,0.13,0.69,0.03,-0.19,0.24,0.71,0.00],
            [0.19,0.04,-0.42,-0.28,-0.20,-0.14,-0.67,-0.13,0.25,0.19,0.19,0.14,0.46,0.08,0.44,0.65,0.99,0.31,0.00,-0.34],
            [-0.23,-0.42,-0.44,-0.19,-0.30,-0.2,-0.16,0.00,0.03,0.38,0.31,0.29,0.49,0.18,0.27,0.39,-0.16,0.41,0.44,0.20],
            [0.16,-0.28,-0.19,-0.22,-0.41,-0.25,0.02,0.11,-0.22,0.25,0.14,0.21,0.36,0.53,0.35,0.59,0.49,0.42,0.36,0.25],
            [-0.08,-0.20,-0.30,-0.41,-0.27,-0.29,-0.09,0.24,-0.01,0.23,0.20,0.25,0.26,0.30,0.43,0.67,0.16,0.35,0.19,0.42],
            [0.06,-0.14,-0.22,-0.25,-0.29,-0.29,-0.07,0.02,-0.10,0.16,0.25,0.18,0.24,0.50,0.34,0.58,0.19,0.30,0.44,0.09],
            [0.08,-0.67,-0.16,0.02,-0.09,-0.07,-0.12,-0.04,-0.09,-0.18,0.22,0.34,0.08,0.06,0.29,0.24,-0.12,-0.16,0.22,-0.28],
            [0.04,-0.13,0.00,0.11,0.24,0.02,-0.04,-0.06,0.09,0.14,0.13,0.09,-0.20,-0.20,-0.10,0.00,-0.34,-0.25,-0.21,-0.33],
            [0.00,0.25,0.03,-0.22,-0.01,-0.10,-0.09,0.09,-0.13,-0.07,-0.09,-0.06,0.08,0.28,0.26,0.12,0.34,0.43,0.14,0.10],
            [-0.08,0.19,0.38,0.25,0.23,0.16,0.18,0.14,-0.07,-0.38,-0.26,-0.16,-0.06,-0.14,0.25,-0.22,0.20,-0.04,0.11,-0.11],
            [0.19,0.19,0.31,0.14,0.20,0.25,0.22,0.13,-0.09,-0.26,0.03,-0.08,-0.14,-0.11,0.00,-0.29,-0.19,-0.35,-0.09,-0.07],
            [-0.02,0.14,0.29,0.21,0.25,0.18,0.34,0.09,-0.06,-0.16,-0.08,-0.20,-0.14,-0.14,0.26,-0.31,-0.05,0.17,-0.13,0.01],
            [0.05,0.46,0.49,0.36,0.26,0.24,0.08,-0.20,0.08,-0.06,-0.14,-0.14,0.29,-0.25,-0.17,-0.17,-0.02,-0.52,-0.38,-0.42],
            [0.13,0.08,0.18,0.53,0.30,0.50,0.06,-0.20,0.28,-0.14,-0.11,-0.14,-0.25,-0.53,-0.32,-0.30,-0.24,-0.14,-0.33,-0.18],
            [0.69,0.44,0.27,0.35,0.43,0.34,0.29,-0.10,0.26,0.25,0.00,-0.26,-0.17,-0.32,-0.03,-0.15,-0.45,-0.74,-0.97,-0.10],
            [0.03,0.65,0.39,0.59,0.67,0.58,0.24,0.00,0.12,-0.22,-0.29,-0.31,-0.17,-0.30,-0.15,0.04,-0.39,-0.72,-0.76,0.04],
            [-0.19,0.99,-0.16,0.49,0.16,0.19,-0.12,-0.34,0.34,0.20,-0.19,-0.05,-0.02,-0.24,-0.45,-0.39,-0.29,-0.12,0.22,-0.21],
            [0.24,0.31,0.41,0.42,0.35,0.30,-0.16,-0.25,0.43,-0.04,-0.35,0.17,-0.52,-0.14,-0.74,-0.72,-0.12,0.11,0.75,-0.38],
            [0.71,0.00,0.44,0.36,0.19,0.44,0.22,-0.21,0.14,0.11,-0.09,-0.13,-0.38,-0.33,-0.97,-0.76,0.22,0.75,0.25,0.11],
            [0.00,-0.34,0.20,0.25,0.42,0.09,-0.28,-0.33,0.10,-0.11,-0.07,0.01,-0.42,-0.18,-0.10,0.04,-0.21,-0.38,0.11,0.26]]

DNA = ['A', 'C', 'T', 'G']

CODON = {'F':['TTT','TTC'], 'L':['TTA','TTG','CTT','CTC','CTA','CTG'],
         'I':['ATT','ATC','ATA'], 'M':['ATG'], 'V':['GTT','GTC','GTA','GTG'],
         'S':['TCT','TCC','TCA','TCG','AGT','AGC'], 'P':['CCT','CCC','CCA','CCG'],
         'T':['ACT','ACC','ACA','ACG'], 'A':['GCT','GCC','GCA','GCG'],'Y':['TAT','TAC'], 
         'stop':['TAA','TAG','TGA'], 'H':['CAT','CAC'], 'Q':['CAA','CAG'], 'N':['AAT','AAC'], 
         'K':['AAA','AAG'], 'D':['GAT','GAC'], 'E':['GAA','GAG'], 'C':['TGT','TGC'], 'W':['TGG'], 
         'R':['CGT','CGC','CGA','CGG','AGA','AGG'],'G':['GGT','GGC','GGA','GGG']
         }

#################################################################################################################################            

#Argument Parser

#################################################################################################################################

parser = argparse.ArgumentParser()

parser.add_argument("-D", "--direct", default=os.getcwd(), 
                    help="where to direct files, default cwd")
parser.add_argument("-N", "--Neff", default=1000000,
                    help="effective population size, default=1000000")
parser.add_argument("-F", "--fitness", default=1,
                    help="if 0 then fitness determined by res contact potentials,\
if 1 then fitness determined from  contact potentials,\
and if 2 then fitness determined by sampling ddG, default=1")
parser.add_argument("-S", "--sequence", default=1,
                    help="if 0 then random seq is seeded at root of tree,\
if 1 then the pdb seq is seeded at root of tree,\
and if 2 then corresponding codon sequence is seeded at root of tree, default=1")
parser.add_argument("-A", "--addendum", type=str, default='',
                    help="extra identifier, like gene, prot fam, tree, etc...")


requiredNamed = parser.add_argument_group('required arguments')

requiredNamed.add_argument("-T", "--tree",
                           help="treefile, newick format, to perform evolution across")
requiredNamed.add_argument("-P", "--pdb-file",
                           help="name of structure model in pdb that evolution will be performed on")


args = parser.parse_args()

#################################################################################################################################            

#Main

#################################################################################################################################

def main():
    #perform evolution
    seq, res, pdb_coord = parse_pdb_file(args.pdb_file)
    tree_path=os.path.join(args.direct,args.tree)
    mut, aln, slxns, probs, prob_of_fix, dGs, ddGs, mut_mat, b4_dG, b4_ddG, tree = EvoSim(seq, res, pdb_coord, tree_path, args.Neff, args.sequence, args.fitness)

    #make some alnmts & plts
    write_tree(tree, args.Neff, args.pdb_file, args.direct, args.addendum) 
    write_aln(aln, args.Neff, args.pdb_file, args.direct, args.addendum)
    plt_slxn_pfix(slxns, probs, prob_of_fix, args.direct, args.Neff, args.pdb_file, args.addendum)
    plt_dG_b4_slxn(args.pdb_file, args.Neff, b4_dG, args.direct, args.addendum)
    plt_dG_af_slxn(args.pdb_file, args.Neff, dGs, args.direct, args.addendum)
    plt_ddG_b4_slxn(args.pdb_file, args.Neff, b4_ddG, args.direct, args.addendum)
    plt_ddG_af_slxn(args.pdb_file, args.Neff, ddGs, args.direct, args.addendum)

if __name__ == '__main__':
    main()
