import os, sys, glob, re
import numpy as np
import pandas as pd
import gzip
import subprocess

# target tissues
tissues = ["Adrenal Gland", "Artery - Tibial", "Kidney - Cortex", "Heart - Left Ventricle"]
tissues_abbr = ["Adrenal", "ArteryTi", "KidneyCo", "HeartLV"]

##################
# Input files
##################

# GENCODE annotations (gene body, exon)
f_gencode = "gencode.v19.annotation.gene.bed"
f_exon_mapping = "exon_mapping_gencode_v19.bed"

# ensemble transcript/gene id mapping
f_gene_to_ens = "gene_id_ens_maping.txt"

# Gene-wise processed data of GTEx RNA-seq (median TPM)
tpm_fn = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"

# top-100k CRE annotations
files_cres = list(map(lambda x: "%s_CREs_100k.bed" % x, tissues_abbr))

##################
# list up gencode genes
##################

f = open(f_gencode)
gencode_genes = []
for line in f.readlines():
    gencode_genes.append(line.split()[-1])
f.close()
gencode_genes = list(set(gencode_genes))


##################
# read TPM data
##################

f = gzip.open(tpm_fn)
f.readline()
f.readline()

columns = f.readline().decode("ascii").strip().split('\t')
tissues_i = []
for tissue in tissues:
    tissues_i.append(columns.index(tissue))
    
tissue_genes = {}
for tissue in tissues:
    tissue_genes[tissue] = []

gene_exp_pair = {}
gene_exp_dist = {}
gene_exp_dic = {}
for tissue in tissues:
    gene_exp_dist[tissue] = []
    gene_exp_pair[tissue] = []
    gene_exp_dic[tissue] = {}
    
for line in f.readlines():
    line_tab = line.decode("ascii").split()
    gene = line_tab[1]
    if not gene in gencode_genes:
        continue
    for i, tissue in enumerate(tissues):
        exp_val = float(line_tab[tissues_i[i]])
        addFlag = True
        if gene in gene_exp_dic[tissue]:
            if exp_val < gene_exp_dic[tissue][gene]:
                addFlag = False
        if addFlag:
            gene_exp_dic[tissue][gene] = exp_val
f.close()

gene_exp_pair = {}
for tissue in gene_exp_dic:
    dic = gene_exp_dic[tissue]
    gene_exp_pair[tissue] = list(dic.items())


##################
# Get Top-N expressed genes for tissues
##################

#rank_tissue_genes = {}
for rank in [5000, 6000, 7000, 8000, 9000, 10000]:

    print("ranked cutoff: %d" % rank)

    tissue_genes = {}
    for tissue in gene_exp_pair:
        l_gene_exp_pair = list(filter(lambda x: gene_annot_dic[x[0]][0][-1].isdigit(), gene_exp_pair[tissue]))
        genes = list(map(lambda x: x[0], sorted(l_gene_exp_pair, key= lambda x: x[1], reverse=True)))[:rank]
        tissue_genes[tissue] = genes
        print("%s: TPM: %.2f" % (tissue, sorted(gene_exp_pair[tissue], key= lambda x: x[1], reverse=True)[rank][1]))
    print("")

    tissue_genes_n = {}
    for at, bt in zip(tissues_abbr, tissues):
        tissue_genes_n[at] = tissue_genes[bt]

        fo = open("GTEx_genes_%s_top%dk.txt" % (at, rank / 1000), "w")
        print("Actual counts: " + str(len(tissue_genes[bt])))
        fo.write('\n'.join(tissue_genes[bt]))
        fo.close()
    

##################
# Get gene body of Top-N expressed genes for tissues
##################

# create bed files
for padding in [10000, 50000]:
    
    # Get annotation profile
    gene_annot_dic = {}
    f = open(f_gencode)
    for line in f.readlines():
        line_tab = line.split()
        gene_annot_dic[line_tab[3]] = [
            line_tab[0],\
            int(line_tab[1]) - padding if int(line_tab[1]) > padding else 0,\
            int(line_tab[2]) + padding
        ]
    f.close()
    
    for rank in [5000, 6000, 7000, 8000, 9000, 10000]:
        for tissue, abbr in zip(tissues, tissues_abbr):
            l_gene_exp_pair = list(filter(lambda x: gene_annot_dic[x[0]][0][-1].isdigit(), gene_exp_pair[tissue]))
            genes = list(map(lambda x: x[0], sorted(l_gene_exp_pair, key= lambda x: x[1], reverse=True)))[:rank]
            
            fo = open("%s_genes_top%dk_%dkb.bed" % (abbr, rank / 1000, padding / 1000), "w")
            for gene in genes[:rank]:
                fo.write("%s\t%d\t%d\t%s\n" % tuple(gene_annot_dic[gene] + [gene]))
            fo.close()

# multiIntersect Bed
for padding in [10000, 50000]:
    for rank in [5000, 6000, 7000, 8000, 9000, 10000]:
        files = list(map(lambda x: "%s_genes_top%dk_%dkb.bed" % (x, rank / 1000, padding / 1000), tissues_abbr))
        files_nsrt = list(map(lambda x: "%s_genes_top%dk_%dkb_nsrt.bed" % (x, rank / 1000, padding / 1000), tissues_abbr))
        files_nsrt_m = list(map(lambda x: "%s_genes_top%dk_%dkb_nsrt.m.bed" % (x, rank / 1000, padding / 1000), tissues_abbr))
        for fi, fi_nsrt in zip(files, files_nsrt):
            os.system("sortBed -i %s > %s" % (fi, fi_nsrt))
        
        for fi_nsrt, fi_nsrt_m in zip(files_nsrt, files_nsrt_m):
            os.system("mergeBed -i %s > %s" % (fi_nsrt, fi_nsrt_m))

        annot_chr = " ".join(tissues_abbr)
        files_chr = " ".join(files_nsrt_m)
        group_fi = "tissues_mintersect_top%dk_%dkb.bed" % (rank / 1000, padding / 1000)
        r = subprocess.getoutput("multiIntersectBed -i %s -names %s > %s" % (files_chr, annot_chr, group_fi))
        
        
##################
# 5-way analysis of
# open-chromatin/intergenic, open-chromatin/intragenic,
# closed-exon, closed-intron, flanking regions
##################

# transcript/gene dic
enst2id_dic = {}
f = open(f_gene_to_ens)
f.readline()
for line in f.readlines():
    symbol, ensg, enst = line.split("\t")[:3]
    enst2id_dic[enst.strip()] = symbol
f.close()

# get the exon dic
gene_exon_dic = {}
f = open(f_exon_mapping)
for line in f.readlines():
    chrn, start, end, exon_annot, _, _ = line.split()
    enst = exon_annot.split('.')[0]
    gene = enst2id_dic[enst]
    if not gene in gene_exon_dic:
        gene_exon_dic[gene] = []
    gene_exon_dic[gene].append((chrn, int(start), int(end), exon_annot))
f.close()

# no padding

# Get annotation profile
gene_annot_dic = {}
f = open(f_gencode)
for line in f.readlines():
    line_tab = line.split()
    gene_annot_dic[line_tab[3]] = [line_tab[0], int(line_tab[1]), int(line_tab[2])]
f.close()

for rank in [5000, 6000, 7000, 8000, 9000, 10000]:
    
    for tissue, abbr in zip(tissues, tissues_abbr):
        l_gene_exp_pair = list(filter(lambda x: gene_annot_dic[x[0]][0][-1].isdigit(), gene_exp_pair[tissue]))
        genes = list(map(lambda x: x[0], sorted(l_gene_exp_pair, key= lambda x: x[1], reverse=True)))[:rank]

        # stricted gene body
        fo = open("%s_genes_top%dk.bed" % (abbr, rank / 1000), "w")
        for gene in genes[:rank]:
            fo.write("%s\t%d\t%d\t%s\n" % tuple(gene_annot_dic[gene] + [gene]))
        fo.close()
        
        # exon
        fo = open("%s_exons_top%dk.bed" % (abbr, rank / 1000), "w")
        for gene in genes[:rank]:
            for exon in gene_exon_dic.get(gene, []):
                fo.write("%s\t%d\t%d\t%s\n" % exon)
        fo.close()
        
# multiIntersect Bed
import subprocess
for rank in [5000, 6000, 7000, 8000, 9000, 10000]:
    print(rank)
    
    # genes
    files_genes = list(map(lambda x: "%s_genes_top%dk.bed" % (x, rank / 1000), tissues_abbr))
    files_genes_nsrt = list(map(lambda x: "%s_genes_top%dk_nsrt.bed" % (x, rank / 1000), tissues_abbr))
    files_genes_nsrt_m = list(map(lambda x: "%s_genes_top%dk_nsrt.m.bed" % (x, rank / 1000), tissues_abbr))
    
    # exons
    files_exons = list(map(lambda x: "%s_exons_top%dk.bed" % (x, rank / 1000), tissues_abbr))
    files_exons_nsrt = list(map(lambda x: "%s_exons_top%dk_nsrt.bed" % (x, rank / 1000), tissues_abbr))
    files_exons_nsrt_m = list(map(lambda x: "%s_exons_top%dk_nsrt.m.bed" % (x, rank / 1000), tissues_abbr))
    
    print("sort genes and exons bed")
    for fi, fi_nsrt in zip(files_genes + files_exons, files_genes_nsrt + files_exons_nsrt):
        os.system("sortBed -i %s > %s" % (fi, fi_nsrt))

    print("merge genes and exons bed")
    for fi_nsrt, fi_nsrt_m in zip(files_genes_nsrt + files_exons_nsrt, files_genes_nsrt_m + files_exons_nsrt_m):
        os.system("mergeBed -i %s > %s" % (fi_nsrt, fi_nsrt_m))
    
    # introns
    print("create introns bed")
    files_introns = list(map(lambda x: "%s_introns_top%dk.bed" % (x, rank / 1000), tissues_abbr))
    files_introns_nsrt = list(map(lambda x: "%s_introns_top%dk_nsrt.bed" % (x, rank / 1000), tissues_abbr))
    files_introns_nsrt_m = list(map(lambda x: "%s_introns_top%dk_nsrt.m.bed" % (x, rank / 1000), tissues_abbr))
    
    for fi_genes, fi_exons, fi_introns in zip(files_genes_nsrt_m, files_exons_nsrt_m, files_introns):
        os.system("subtractBed -a %s -b %s > %s" % (fi_genes, fi_exons, fi_introns))
        
    for fi, fi_nsrt in zip(files_introns, files_introns_nsrt):
        os.system("sortBed -i %s > %s" % (fi, fi_nsrt))

    for fi_nsrt, fi_nsrt_m in zip(files_introns_nsrt, files_introns_nsrt_m):
        os.system("mergeBed -i %s > %s" % (fi_nsrt, fi_nsrt_m))
    
    # flanking regions
    print("flanks bed")
    for padding in [10000, 50000]:
        print(padding)
        files_padding_nsrt_m = list(map(lambda x: "%s_genes_top%dk_%dkb_nsrt.m.bed" % (x, rank / 1000, padding / 1000), tissues_abbr))
        files_flanks = list(map(lambda x: "%s_flanks_top%dk_%dkb.bed" % (x, rank / 1000, padding / 1000), tissues_abbr))
        files_flanks_nsrt = list(map(lambda x: "%s_flanks_top%dk_%dkb_nsrt.bed" % (x, rank / 1000, padding / 1000), tissues_abbr))
        files_flanks_nsrt_m = list(map(lambda x: "%s_flanks_top%dk_%dkb_nsrt.m.bed" % (x, rank / 1000, padding / 1000), tissues_abbr))
        
        print("create")
        for fi_padding, fi_genes, fi_flanks in zip(files_padding_nsrt_m, files_genes_nsrt_m, files_flanks):
            os.system("subtractBed -a %s -b %s > %s" % (fi_padding, fi_genes, fi_flanks))
            
        print("sort")
        for fi, fi_nsrt in zip(files_flanks, files_flanks_nsrt):
            os.system("sortBed -i %s > %s" % (fi, fi_nsrt))

        print("merge")
        for fi_nsrt, fi_nsrt_m in zip(files_flanks_nsrt, files_flanks_nsrt_m):
            os.system("mergeBed -i %s > %s" % (fi_nsrt, fi_nsrt_m))

        annot_chr = " ".join(tissues_abbr)
        files_chr = " ".join(files_nsrt_m)
        group_fi = "../data/tissues_mintersect_%s_top%dk.bed" % (elements, padding / 1000)
        r = subprocess.getoutput("multiIntersectBed -i %s -names %s > %s" % (files_chr, annot_chr, group_fi))

# CREs_intergenic
files_cres_inter = list(map(lambda x: "%s_CREs_inter_100k.bed" % x, tissues_abbr))
files_cres_intra = list(map(lambda x: "%s_CREs_intra_100k.bed" % x, tissues_abbr))

for fi_cres, fi_cres_inter, fi_cres_intra in zip(files_cres, files_cres_inter, files_cres_intra):
    os.system("subtractBed -a %s -b %s > %s" % (fi_cres, f_gencode, fi_cres_inter))
    os.system("intersectBed -a %s -b %s > %s" % (fi_cres, f_gencode, fi_cres_intra))
    
for rank in [5000, 6000, 7000, 8000, 9000, 10000]:
    print(rank)
    
    # genes - flanks 50kb
    files_genes_nsrt_m = list(map(lambda x: "%s_genes_top%dk_50kb_nsrt.m.bed" % (x, rank / 1000), tissues_abbr))
    
    # CREs intragenic/intergenic
    files_cres_inter_nsrt_m = list(map(lambda x: "%s_CREs_inter_top%dk_nsrt.m.bed" % (x, rank / 1000), tissues_abbr))
    files_cres_intra_nsrt_m = list(map(lambda x: "%s_CREs_intra_top%dk_nsrt.m.bed" % (x, rank / 1000), tissues_abbr))

    print("Get intra/inter CREs with respect to the gene bodies of tissue-expressed genes")
    for fi_cres, fi_genes, fi_cres_inter, fi_cres_intra in zip(files_cres, files_genes_nsrt_m, files_cres_inter_nsrt_m, files_cres_intra_nsrt_m):
        os.system("subtractBed -a %s -b %s > %s" % (fi_cres, fi_genes, fi_cres_inter))
        os.system("intersectBed -a %s -b %s > %s" % (fi_cres, fi_genes, fi_cres_intra))
    
    # exons/introns/flanks
    files_exons_nsrt_m = list(map(lambda x: "%s_exons_top%dk_nsrt.m.bed" % (x, rank / 1000), tissues_abbr))
    files_introns_nsrt_m = list(map(lambda x: "%s_introns_top%dk_nsrt.m.bed" % (x, rank / 1000), tissues_abbr))
    files_flanks_nsrt_m = list(map(lambda x: "%s_flanks_top%dk_50kb_nsrt.m.bed" % (x, rank / 1000), tissues_abbr))
    
    # exons/introns/flanks-closed
    files_exons_closed_nsrt_m = list(map(lambda x: "%s_exons_closed_top%dk_nsrt.m.bed" % (x, rank / 1000), tissues_abbr))
    files_introns_closed_nsrt_m = list(map(lambda x: "%s_introns_closed_top%dk_nsrt.m.bed" % (x, rank / 1000), tissues_abbr))
    files_flanks_closed_nsrt_m = list(map(lambda x: "%s_flanks_closed_top%dk_nsrt.m.bed" % (x, rank / 1000), tissues_abbr))
    
    print("Get exon/intron/flanks-closed regions")
    for fi, fi_cres_intra, fi_closed in zip(files_exons_nsrt_m + files_introns_nsrt_m + files_flanks_nsrt_m,\
                                            files_cres_intra_nsrt_m + files_cres_intra_nsrt_m + files_cres_intra_nsrt_m,\
                                            files_exons_closed_nsrt_m + files_introns_closed_nsrt_m + files_flanks_closed_nsrt_m):
        os.system("subtractBed -a %s -b %s > %s" % (fi, fi_cres_intra, fi_closed))