#!/usr/bin/env python3

#################################################################################################################################
#Script Name    : EvoSim.py
#Description    : simulate evolution along a tree using a model dependent on protein stability and population genetics
#Author         : Michael A Sennett & Douglas Theobald
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
        if line[3] != 'GLY':
            if line[4] != chain[0] or line[2] != 'CB':
                pass
            else:
                pdb_2d_arr_A.append(list(line))
        else:
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

def get_rdn_coil_pdbs():
    rc_pdb=['6BQV', '6R8B', '6WBG', '7M3F', '6X2K', '6VV5', '7VHM', '6S6L', '6YJ6', '8DG4',
        '7P5C', '7DFV', '5IS0', '7W01', '6UX4', '7OQH', '7JZ3', '6U88', '7R91', '7V3U',
        '7QKW', '7PBD', '7WJO', '7ZZ2', '7WIU', '6VZ1', '7KF6', '6UW4', '7XNR', '7RHL',
        '6ZP2', '7UNM', '7UNQ', '6PM2', '7JQ9', '6WXI', '7N6B', '7WZ2', '7WUB', '6OG1',
        '6U5N', '6M15', '7RJT', '7S89', '7WS8', '6VXN', '7LH2', '5IRZ', '7UV0', '7T4N',
        '7EZ0', '7B5P', '6NPL']

    #seq, res, pdb_coord = parse_pdb_file(pdb)

    rc_pdb_coord=[]
    rc_res=[]
    for rc in rc_pdb:
        try:
            tmp_seq, tmp_res, tmp_pdb_coord = parse_pdb_file(rc)
            rc_pdb_coord.append(tmp_pdb_coord)
            rc_res.append(tmp_res)
        except:
            pass

    return rc_res, rc_pdb_coord

#################################################################################################################################

#create a random sequence
def gen_rd_seq(seq_len):
    
    rd_seq=list()
    for res in range(len(seq_len)):
        rd_seq.append(rd.sample(AA_miya, 1)[0])

    return rd_seq

#################################################################################################################################

#grab list of sites that are near residue of interest in a single chain
def nn_sites(pdb_2d_arr, seq_len):

    pdb_2d_arr=pdb_2d_arr
    
    #create a look-up table for each site in a sequence
    site=[]
    nn_site=[]
    for line in pdb_2d_arr:
        site_coord=np.zeros((1,3))
    #get the ref residue coord
        
        if line[5] not in site and int(line[5]) < seq_len:
            site_coord[0,0]=float(line[6])
            site_coord[0,1]=float(line[7])
            site_coord[0,2]=float(line[8])

    #check to see if next residue is close
            tmp_nn_sites=[]    
            for new_line in pdb_2d_arr:
                nn_coord=np.zeros((1,3))
                if new_line[5] != line[5] and int(new_line[5]) < seq_len:
                    
                    nn_coord[0,0]=float(new_line[6])
                    nn_coord[0,1]=float(new_line[7])
                    nn_coord[0,2]=float(new_line[8])
    #calc dist bw sites
                    dist=np.sum((site_coord-nn_coord)**2,1)**0.5
                    if dist < 7:
                        if int(line[5]) < int(new_line[5]):
                            site.append(line[5])
                            nn_site.append(new_line[5])
                        else:
                            site.append(new_line[5])
                            nn_site.append(line[5])
                    else:
                        pass
                else:
                    pass
        else:
            pass

    lookup_table=list(zip(site,nn_site))
    nr_lookup_table=[]
    for item in lookup_table:
        if item not in nr_lookup_table:
            nr_lookup_table.append(item)
        
        else:
            pass
        
    
    return nr_lookup_table     

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
        
def cdn_2_aa(cdn_seq):
    #create an aa seq for later
    tmp_seq=[]
    for cdn in cdn_seq:
        aa = get_cdn_key(cdn)
        tmp_seq.append(aa)

    return tmp_seq

#################################################################################################################################

def calc_G_F(tmp_seq, lookup_table, res_nums, AA_miya):


    res_dct=dict(zip(res_nums, range(len(res_nums))))

    #calculate GNS, native state free energy, for some struct/seq
    fold_GNS=[]
    for res_pair in lookup_table:
        #convert to a int
        one=res_pair[0]
        two=res_pair[1]
        #find int in dict to get linear seq number
        aa_one=res_dct.get(one)
        aa_two=res_dct.get(two)
        #identify aa from linear seq
        seq_res_one=tmp_seq[aa_one]
        seq_res_two=tmp_seq[aa_two]
        #calc CP between two residues
        miya_id_one=AA_miya_dct.get(seq_res_one)
        miya_id_two=AA_miya_dct.get(seq_res_two)
        fold_GNS.append(miyazawa[miya_id_one][miya_id_two])

    tot_fold_GNS=sum(fold_GNS)

    return tot_fold_GNS

#################################################################################################################################

def intro_mut(cdn_seq):

    #position in seq to mutate

    site_mut=rd.sample(range(1,len(cdn_seq)),1)
    wt_cd=cdn_seq[site_mut[0]]

    roll = rd.uniform(0,1)
    if roll < 0.9:
        #pos in wt codon to be mutated
        pos=rd.sample(range(3), 1)
        cd_pos=pos[0]
        #ACTG
        if wt_cd[cd_pos] == 'A':
            mut_num_lst=rd.sample([1,2,3,3], 4)
        elif wt_cd[cd_pos] == 'C':
            mut_num_lst=rd.sample([0,2,2,3], 4)
        elif wt_cd[cd_pos] == 'T':
            mut_num_lst=rd.sample([0,1,1,3], 4)
        else:
            mut_num_lst=rd.sample([0,0,1,2], 4)

        #sampling from nuc dist with equal probs

        #propose mutation
        for mut_num in mut_num_lst:
            mut_cd=list(wt_cd)
            mut_cd[cd_pos]=DNA[mut_num]
            mut_cd=''.join(mut_cd)
            if DNA[mut_num] == wt_cd[cd_pos]:
                pass
            elif get_cdn_key(mut_cd) == 'stop':
                pass
            else:
                break

        var_cdn_seq=list('X')*len(cdn_seq)

        for cdn in range(len(cdn_seq)):
            if cdn != site_mut[0]:
                var_cdn_seq[cdn]=cdn_seq[cdn]
            else:
                var_cdn_seq[cdn]=mut_cd
    else:
        cdn_lst=[]
        for val in CODON.values():
            for v in val:
                cdn_lst.append(v)
        cdn_smp=rd.sample(cdn_lst, len(cdn_lst))

        #propose mutation
        for c in cdn_smp:
            if c == wt_cd:
                pass
            elif get_cdn_key(c) == 'stop':
                pass
            else:
                break
        mut_cd = c

        var_cdn_seq=list('X')*len(cdn_seq)

        for cdn in range(len(cdn_seq)):
            if cdn != site_mut[0]:
                var_cdn_seq[cdn]=cdn_seq[cdn]
            else:
                var_cdn_seq[cdn]=mut_cd


    return wt_cd, mut_cd, site_mut[0], var_cdn_seq

#################################################################################################################################

def calc_Pfix(wt_aa_seq, var_aa_seq, lookup_table, res_nums, rc_lookup_tables, rc_res_dct, N_eff, m, fit, scale):

    if fit == 0:

        seq_len=len(wt_aa_seq)
        #calc the native state G for the initial seq/struct
        G_i=calc_G_F(wt_aa_seq, lookup_table, res_nums, AA_miya)
        #calc the average random coil G for the initial seq/struct
        G_i_rc=[]
        for j in range(len(rc_lookup_tables)):
            G_i_tmp=calc_G_F(wt_aa_seq, rc_lookup_tables[j], rc_res_dct[j], AA_miya)
            G_i_rc.append(G_i_tmp)
        avg_G_i_rc=np.mean(G_i_rc)
        std_G_i_rc=np.std(G_i_rc)
        RT=0.593
        dG_i=G_i+RT*np.log(10.**160)+((std_G_i_rc**2)-(2.*RT*avg_G_i_rc))/(2.*RT)
        #calc the native state G for the final seq/struct
        G_v=calc_G_F(var_aa_seq, lookup_table, res_nums, AA_miya)
        #calc the average random coil G for the final seq/struct
        G_v_rc=[]
        for j in range(len(rc_lookup_tables)):
            G_v_tmp=calc_G_F(var_aa_seq, rc_lookup_tables[j], rc_res_dct[j], AA_miya)
            G_v_rc.append(G_v_tmp)
        avg_G_v_rc=np.mean(G_v_rc)
        std_G_v_rc=np.std(G_v_rc)
        dG_v=G_v+RT*np.log(10.**160)+((std_G_v_rc**2)-(2.*RT*avg_G_v_rc))/(2.*RT)
    
        ddG=dG_v-dG_i
        wv=np.longdouble(1)/(np.longdouble(1)+np.exp(dG_v/RT))
        wi=np.longdouble(1)/(np.longdouble(1)+np.exp(dG_i/RT))
    
        s_mut=(wv-wi)/wi
    
        if s_mut == 0:
            fv=1/N_eff*scale
        else:
            fv=np.expm1(-2*s_mut)/np.expm1(-4*N_eff*s_mut)*scale

        return fv, dG_i, G_i, avg_G_i_rc, std_G_i_rc, dG_v, G_v, avg_G_v_rc, std_G_v_rc, ddG, s_mut

#################################################################################################################################

#################################################################################################################################

def EvoSim(wt_seq, res_nums, pdb_2d_arr, rc_res_lst, rc_pdb_2d_arr, tree, Neff, seq, fit, burn):

    #user defined initial variables from input
    N_eff=Neff
    fit=fit
    seq=seq
    burn_mut=burn
    
    
    #non-user defined initial variables from input
    res_dct=dict(zip(res_nums, range(len(res_nums))))
    seq_len=len(res_nums)
    dG_WT=(-1*np.log(N_eff)-0.57)/1
    lookup_table=nn_sites(pdb_2d_arr, seq_len)
    rc_res_dct=[]
    for rc_res in rc_res_lst:
        tmp_res_dct=dict(zip(rc_res[0:seq_len], range(0,seq_len)))
        rc_res_dct.append(tmp_res_dct)

    rc_lookup_tables=[]
    for rc in rc_pdb_2d_arr:
        tmp_LT=nn_sites(rc, seq_len)
        rc_lookup_tables.append(tmp_LT)

    #sample a random number and then throw into seed so we can track
    seed=rd.randint(0,100000000)
    print("random seed number:", seed)
    rd.seed(seed)
    
    #create a random sequence if you want to
    if seq == 0:
        rd_seq=gen_rd_seq(seq_len)
        wt_seq=rd_seq
        cd_seq=get_DNA_CDS(wt_seq)

    #create a coding sequence if you want to
    else:
        cd_seq=get_DNA_CDS(wt_seq)

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
    af_G=[]
    af_aG_rc=[]
    af_sG_rc=[]
    eq_res=[]

    #traverse the tree and sim evo
    dG_i=dG_WT


    #### burn-in here so that when traverse tree it will start from eq ####
    b_bf_dGs=[]
    burn_mut=100000
    tmp_b_dGs=[]
    b_tmp_cd_seq="".join(cd_seq)
    
    b_new_tmp=[]
    for j in range(0, len(b_tmp_cd_seq), 3):
        b_new_tmp.append(b_tmp_cd_seq[j:j+3])
    b_tmp_seq=b_new_tmp
    b=0
    fit=0
    scale=1
    while b < burn_mut:
        wt_b_cdn_seq=b_tmp_seq
        wt_b_aa_seq=cdn_2_aa(wt_b_cdn_seq)
        wt_b_cdn, var_b_cdn, site_b, var_b_cdn_seq=intro_mut(wt_b_cdn_seq)
        var_b_aa_seq=cdn_2_aa(var_b_cdn_seq)
        wt_b_aa=get_cdn_key(wt_b_cdn)
        var_b_aa=get_cdn_key(var_b_cdn)
        if var_b_aa == wt_b_aa and b > 0:
            fv=1/N_eff*scale
            b_ddG=np.longdouble(0.0)
            b_s_mut=np.longdouble(0.0)
            b_dG_i=tmp_b_dGs[-1]
            b_dG_v=tmp_b_dGs[-1]
        else:
            fv, dG_i, G_i, avg_G_i_rc, std_G_i_rc, dG_v, G_v, avg_G_v_rc, std_G_v_rc, ddG, s_mut = calc_Pfix(wt_b_aa_seq, var_b_aa_seq, lookup_table, res_nums, rc_lookup_tables, rc_res_dct, N_eff, b, fit, scale)
        b_bf_dGs.append(np.longdouble(dG_v))
        
        dice=rd.uniform(0,1)
        if dice < fv:
            dG_i=dG_v
            b_tmp_seq=var_b_cdn_seq
            tmp_b_dGs.append(np.longdouble(dG_i))
        else:
            dG_i=dG_i
            tmp_b_dGs.append(np.longdouble(dG_i))
            
        if b%(10000) == 0:
            #scale=scale*10
            print(int(b/burn_mut*100), '% done')
        else:
            pass
        
        b+=1
        
    cd_seq=b_tmp_seq

    ######  end burn in #####

    for node in sim_tree.preorder_node_iter():

        if node.taxon.label not in aln:

            aln[node.taxon.label]="".join(cd_seq)
            seq_dGs[node.taxon.label]=dG_i
            
        x=node.child_nodes()

        #for each of the child nodes determine the numb of mutations expected
        for i in range(len(x)):
            
            exp_mut=int(x[i].edge_length*len(wt_seq))
            tmp_dGs=[]

            tmp_seq=aln.get(node.taxon.label)
            new_tmp=[]
            for j in range(0, len(tmp_seq), 3):
                new_tmp.append(tmp_seq[j:j+3])
            tmp_seq=new_tmp

            dG_i=seq_dGs.get(node.taxon.label)

            if x[i].taxon.label not in aln:

                print(x[i].taxon.label, ", add as child and mutate")

                m=0

                if exp_mut == 0:
                    exp_mut=1
                
                while m < exp_mut:
                    scale=N_eff/10
                    #intro mut and convert codon seq to aa seq
                    wt_cdn_seq=tmp_seq
                    wt_aa_seq=cdn_2_aa(tmp_seq)
                    wt_cdn, var_cdn, site, var_cdn_seq=intro_mut(wt_cdn_seq)
                    var_aa_seq=cdn_2_aa(var_cdn_seq)
                    wt_aa=get_cdn_key(wt_cdn)
                    var_aa=get_cdn_key(var_cdn)

                    if var_aa == wt_aa and m > 0:

                        fv=1/N_eff*scale
                        ddG=np.longdouble(0.0)
                        s_mut=np.longdouble(0.0)
                        dG_i=tmp_dGs[-1]
                        dG_v=tmp_dGs[-1]

                    else:
                        
                        fv, dG_i, G_i, avg_G_i_rc, std_G_i_rc, dG_v, G_v, avg_G_v_rc, std_G_v_rc, ddG, s_mut = calc_Pfix(wt_aa_seq, var_aa_seq, lookup_table, res_nums, rc_lookup_tables, rc_res_dct, N_eff, m, fit, scale)
                    bf_dGs.append(np.longdouble(dG_v))
                    bf_ddG.append(np.longdouble(ddG))
                    dice=rd.uniform(0,1)

                    if dice < fv:

                        G_i=G_v
                        avg_G_i_rc=avg_G_v_rc
                        std_G_i_rc=std_G_v_rc
                        af_G.append(G_i)
                        af_aG_rc.append(avg_G_i_rc)
                        af_sG_rc.append(std_G_i_rc)
                        wt_id=AA_miya_dct.get(wt_aa)
                        var_id=AA_miya_dct.get(var_aa)
                        mut_mat[wt_id][var_id]+=1
                        eq_res.append(var_aa)
                        dG_i=dG_v
                        tmp_seq=var_cdn_seq
                        af_probs.append(dice)
                        af_prob_of_fix.append(np.longdouble(fv))
                        af_slxns.append(np.longdouble(s_mut))
                        tmp_dGs.append(np.longdouble(dG_i))
                        af_ddGs.append(np.longdouble(ddG))

                    else:
                        dG_i=dG_i
                        tmp_dGs.append(np.longdouble(dG_i))
                        af_G.append(G_i)
                        af_aG_rc.append(avg_G_i_rc)
                        af_sG_rc.append(std_G_i_rc)
                        eq_res.append(wt_aa)

                    m+=1


                aln[x[i].taxon.label]="".join(tmp_seq)

            else:
                pass

            seq_dGs[x[i].taxon.label]=tmp_dGs[-1]
            af_dGs.append(tmp_dGs)

    out_tree=sim_tree.as_string('newick')

    return m, aln, af_slxns, af_probs, af_prob_of_fix, af_dGs, af_ddGs, mut_mat, eq_res, bf_dGs, bf_ddG, out_tree, af_G, af_aG_rc, af_sG_rc

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

    title1=pdb+'_Neff_'+str(Neff)+'_b4_dGs_'+add+'.png'
    plt.savefig(os.path.join(direct,title1), bbox_inches='tight', pad_inches=0.3, dpi=300)

    title2=pdb+'_Neff_'+str(Neff)+'_b4_dGs_'+add+'.txt'
    with open(os.path.join(direct, title2), 'w+') as outfile:
        for line in b4_dG:
            outfile.write(str(line)+'\n')
        outfile.close()


#################################################################################################################################
#plt the dG of protein mutants after selection occurs

def plt_dG_af_slxn(pdb, Neff, dGs, direct, add):
    
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

    title1=pdb+'_Neff_'+str(Neff)+'_af_dGs_'+add+'.png'
    plt.savefig(os.path.join(direct,title1), bbox_inches='tight', pad_inches=0.3, dpi=300)

    title2=pdb+'_Neff_'+str(Neff)+'_af_dGs_'+add+'.txt'
    with open(os.path.join(direct, title2), 'w+') as outfile:
        for line in all_dGs:
            outfile.write(str(line)+'\n')
        outfile.close()

################################################################################################################################# 

def plt_ddG_b4_slxn(pdb, Neff, b4_ddG, direct, add):

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

    title1=pdb+'_Neff_'+str(Neff)+'_b4_ddGs_'+add+'.png'
    plt.savefig(os.path.join(direct,title1), bbox_inches='tight', pad_inches=0.3, dpi=300)

    title2=pdb+'_Neff_'+str(Neff)+'_b4_ddGs_'+add+'.txt'
    with open(os.path.join(direct, title2), 'w+') as outfile:
        for line in b4_ddG:
            outfile.write(str(line)+'\n')
        outfile.close()

#################################################################################################################################

def plt_ddG_af_slxn(pdb, Neff, af_ddG, direct, add):

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

    title1=pdb+'_Neff_'+str(Neff)+'_af_ddGs_'+add+'.png'
    plt.savefig(os.path.join(direct,title1), bbox_inches='tight', pad_inches=0.3, dpi=300)

    title2=pdb+'_Neff_'+str(Neff)+'_af_ddGs_'+add+'.txt'
    with open(os.path.join(direct, title2), 'w+') as outfile:
        for line in af_ddG:
            outfile.write(str(line)+'\n')
        outfile.close()

#################################################################################################################################

def plt_sub_mat(pdb, Neff, mut_mat, direct, add):
    
    new_mat=np.zeros((20,20))
    off_by=[]
    for i in range(len(mut_mat)):
        for j in range(len(mut_mat)):
            #eq, lt, gt
            if j == i:
                pass
            elif j < i:
                new_mat[i,j]=mut_mat[i,j]+mut_mat[j,i]
            else:
                pass

    plt.rc('font', family='courier')
    plt.rcParams.update({'font.size': 12})
    plt.figure(figsize=[5,5])

    x_=[]
    y_=[]
    dots=[]
    for i in range(len(AA_miya)):
        for j in range(len(AA_miya)):
            y_.append(i)
            x_.append(j)
            dots.append(new_mat[i,j])

    plt.scatter(x_,y_,color='black',s=dots)
    plt.yticks(ticks=range(len(AA_miya)),labels=AA_miya)
    plt.xticks(ticks=range(len(AA_miya)),labels=AA_miya)
    plt.title('substitutions')
    plt.grid(which='major', color='black', linewidth=0.5)
    plt.ylim(top=-1, bottom=20)

    title1=pdb+'_Neff_'+str(Neff)+'_submat_'+add+'.png'
    plt.savefig(os.path.join(direct,title1), bbox_inches='tight', pad_inches=0.3, dpi=300)

    title2=pdb+'_Neff_'+str(Neff)+'_submat_'+add+'.txt'
    with open(os.path.join(direct, title2), 'w+') as outfile:
        for line in new_mat:
            outfile.write(str(line)+'\n')
        outfile.close()

#################################################################################################################################

def plt_eq_freq(pdb, Neff, eq_res, direct, add):
    freqs=[]
    for a in AA_miya:
        c=0
        for res in eq_res:
            if a == res:
                c+=1
            else:
                pass
        freqs.append(c/len(eq_res))


    plt.rc('font', family='courier')
    plt.rcParams.update({'font.size': 12})
    plt.figure(figsize=[8,4.5])

    aa_range=np.arange(len(AA_miya))

    plt.bar(aa_range, freqs, width=0.4, label='Sim eq freqs')
    plt.ylabel('frequency')
    plt.xlabel('amino acid')
    plt.xticks(aa_range, AA_miya)
    plt.legend(loc='upper left')


    title1=pdb+'_Neff_'+str(Neff)+'_eqfreq_'+add+'.png'
    plt.savefig(os.path.join(direct,title1), bbox_inches='tight', pad_inches=0.3, dpi=300)

    title2=pdb+'_Neff_'+str(Neff)+'_eqfreq_'+add+'.txt'
    with open(os.path.join(direct, title2), 'w+') as outfile:
        outfile.write(str(AA_miya)+'\n')
        for line in freqs:
            outfile.write(str(line)+'\n')
        outfile.close()

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

RC_pdb = ['6BQV', '6R8B', '6WBG', '7M3F', '6X2K', '6VV5', '7VHM', '6S6L', '6YJ6', '8DG4', 
          '7P5C', '7DFV', '5IS0', '7W01', '6UX4', '7OQH', '7JZ3', '6U88', '7R91', '7V3U', 
          '7QKW', '7PBD', '7WJO', '7ZZ2', '7WIU', '6VZ1', '7KF6', '6UW4', '7XNR', '7RHL', 
          '6ZP2', '7UNM', '7UNQ', '6PM2', '7JQ9', '6WXI', '7N6B', '7WZ2', '7WUB', '6OG1', 
          '6U5N', '6M15', '7RJT', '7S89', '7WS8', '6VXN', '7LH2', '5IRZ', '7UV0', '7T4N', 
          '7EZ0', '7B5P', '6NPL']

#################################################################################################################################            

#Argument Parser

#################################################################################################################################

parser = argparse.ArgumentParser()

parser.add_argument("-D", "--direct", default=os.getcwd(), 
                    help="where to direct files, default cwd")
parser.add_argument("-N", "--Neff", default=1000000,
                    help="effective population size, default=1000000")
parser.add_argument("-F", "--fitness", default=0,
                    help="if 0 then fitness determined by res contact potentials,\
default=0")
parser.add_argument("-S", "--sequence", default=1,
                    help="if 0 then random AA seq is seeded at root of tree,\
if 1 then a random codon sequence corresponding to the pdb AA seq is seeded at root of tree,\
default=1")
parser.add_argument("-A", "--addendum", type=str, default='',
                    help="extra identifier, like gene, prot fam, tree, etc...")
parser.add_argument("-B", "--burn-in", type=int, default=100000,
                    help="number of proposed mutations to reach equilibrium,\
the default=100000")


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
    rc_res, rc_pdb_coord = get_rdn_coil_pdbs()    
    tree_path=os.path.join(args.direct,args.tree)
    mut, aln, slxns, probs, prob_of_fix, dGs, ddGs, mut_mat, eq_res, b4_dG, b4_ddG, tree, af_G, af_aG_rc, af_sG_rc=EvoSim(seq, res, pdb_coord, rc_res, rc_pdb_coord, tree_path, args.Neff, args.sequence, args.fitness, args.burn_in)
    
    #make some alnmts & plts
    write_tree(tree, args.Neff, args.pdb_file, args.direct, args.addendum) 
    write_aln(aln, args.Neff, args.pdb_file, args.direct, args.addendum)
    plt_slxn_pfix(slxns, probs, prob_of_fix, args.direct, args.Neff, args.pdb_file, args.addendum)
    plt_dG_b4_slxn(args.pdb_file, args.Neff, b4_dG, args.direct, args.addendum)
    plt_dG_af_slxn(args.pdb_file, args.Neff, dGs, args.direct, args.addendum)
    plt_ddG_b4_slxn(args.pdb_file, args.Neff, b4_ddG, args.direct, args.addendum)
    plt_ddG_af_slxn(args.pdb_file, args.Neff, ddGs, args.direct, args.addendum)
    plt_sub_mat(args.pdb_file, args.Neff, mut_mat, args.direct, args.addendum)
    plt_eq_freq(args.pdb_file, args.Neff, eq_res, args.direct, args.addendum)

if __name__ == '__main__':
    main()
    

