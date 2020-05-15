import pandas as pd
import numpy as np
from math import log, floor, ceil
from scipy.stats import pearsonr, gaussian_kde, ks_2samp
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from sklearn.metrics import auc, precision_recall_curve, roc_curve, roc_auc_score, average_precision_score
import os, sys, svgwrite

# This function returns the reverse complement of a given sequence
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N','K':'M','M':'K','R':'Y','Y':'R','S':'S','W':'W','B':'V','V':'B','H':'D','D':'H','-':'-'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

def GetMostSevereClinSig(string):
    if type(string) == float:
        return 'None'
    type_list = string.split(';')
    if 'Pathogenic' in type_list:
        return 'Pathogenic'
    elif 'Likely pathogenic' in type_list:
        return 'Likely pathogenic'
    elif 'Pathogenic/Likely pathogenic' in type_list:
        return 'Likely pathogenic'
    elif 'Uncertain significance' in type_list:
        return 'Uncertain significance'
    elif 'Conflicting interpretations of pathogenicity' in type_list:
        return 'Conflicting interpretations of pathogenicity'
    elif 'None' in type_list:
        return 'Not in ClinVar'
    elif 'Likely benign' in type_list:
        return 'Likely benign'
    elif 'Benign/Likely benign' in type_list:
        return 'Likely benign'
    elif 'Benign' in type_list:
        return 'Benign'

def GetGeneCategory(construct_id):
    if 'NO_SITE' in construct_id:
        return 'NO_SITE'
    elif 'ONE_NON-GENE_SITE' in construct_id:
        return 'ONE_NON-GENE_SITE'
    elif 'Panlethal splice donor' in construct_id:
        return 'Panlethal splice donor'
    elif 'Essential splice site' in construct_id:
        return 'Essential splice site'
    elif 'Nonessential splice site' in construct_id:
        return 'Nonessential splice site'    
    elif 'Mutant' in construct_id:
        return 'Mutant drug target'    
    else:
        return ''

def GetMostSevereMutationType(string):
    if type(string) == float:
        return 'No edits'
    elif 'Nonsense' in string:
        return 'Nonsense'
    elif 'Splice-acceptor' in string or 'Splice-donor' in string:
        return 'Splice site'
    elif 'Missense' in string:
        return 'Missense'
    elif 'Intron' in string:
        return 'Intron'        
    elif 'Silent' in string:
        return 'Silent'
    elif 'UTR' in string:
        return 'UTR'

def make_ridgeplots(data, col, cutoff, lims, ticks, geneset, control_geneset,
                    xlabel='z-score',
                    type_list=['Targeting controls','No edits','Silent','Missense','Splice site','Nonsense']):

    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    control_data = data[data['Gene symbol'].isin(control_geneset)][col]
    if cutoff < 0:
        ha = 'left'
        fraction_ha ='left'
        text_start = 0
        fraction_start = 0
    elif cutoff > 0:
        ha = 'right'
        fraction_ha = 'right'
        text_start = 0.98
        fraction_start = 0.98

    fig,ax = plt.subplots(6,figsize=(3,2),sharex=True)
    for i,mut_type in enumerate(type_list):
        # Plot controls first
        if i == 0:
            curr_data = control_data
        else:
            curr_data = data.loc[(data['Mutation bin'] == mut_type) & (data['Gene symbol'].isin(geneset)),col]
        sns.kdeplot(curr_data,color=sns.color_palette("Set2")[i],ax=ax[i],shade=True,legend=False,clip_on=False,alpha=0.7,lw=0.5,zorder=1)
        sns.kdeplot(curr_data,clip_on=False,color='black',lw=0.25,ax=ax[i],legend=False,zorder=1)
        ax[i].set_title('')
        ax[i].set_xticks([])
        ax[i].set_yticks([])
        # displays the fraction of guides to the left or right of the cutoff (depending on the percentile cutoff)
        if cutoff < 0:
            fraction = float(len([datum for datum in curr_data if datum < cutoff]))/float(len(curr_data)) 
        if cutoff > 0:
            fraction = float(len([datum for datum in curr_data if datum > cutoff]))/float(len(curr_data))                     
        ax[i].text(fraction_start,0.05,"{:.1%}".format(fraction),fontsize=6,ha=fraction_ha,transform=ax[i].transAxes) 

        # tidy up plot               
        ax[i].set_xlim(lims)
        ax[i].text(text_start,0.4,type_list[i]+', n = %i' % (len(curr_data)),ha=ha,fontsize=6,transform=ax[i].transAxes)
        ax[i].spines['bottom'].set_linewidth(0.5)
        
    plt.gca().axvline(x=cutoff,ymin=0,ymax=5,c="black",linestyle='dashed',linewidth=0.5,zorder=3, clip_on=False)
    plt.gca().set_xticks(ticks)
    plt.xticks(fontsize=6)
    plt.gca().set_xlabel(xlabel,fontsize=6)
    sns.despine(top=True,bottom=False,left=True,right=True)

    # Set the subplots to overlap
    plt.subplots_adjust(hspace=-.25)

    plt.show()
    plt.tight_layout()
    return fig,ax

def GetResidues(string):
    new_string = ''
    if type(string) != float:
        edits = string.split(';')
        for edit in edits:
            if edit.startswith('Exon'):
                new_string += 'intron;'
            elif edit == 'utr':
                new_string += 'utr;'
            elif edit == '':
                continue
            else:
                # Strip off all non-digit characters
                for character in edit:
                    if character.isdigit():        
                        new_string += character
                new_string += ';'
    return new_string

'''
GetMedianResidues_v2 is updated so that sgRNAs that are binned as "Missense" but contain
intronic (not splice site!) edits still get a "median residue" and can appear on a protein plot. This
is consistent with the ordering of GetMostSevereMutationType, in which Missense > Intron.

Note that sgRNAs with 'Mutation bin' == 'Splice site' still do NOT receive a "median residue," consistent
with Splice site > Missense in mutation bin ordering.
'''
    
def GetMedianResidues_v2(row):
    residues = row['Residues'].split(';')
    if (row['Mutation bin'] == 'Splice site') or (row['Mutation bin'] == 'UTR'):
        return np.nan
    residues = [int(res) for res in residues if res not in ['','intron','utr']]
    if len(residues) != 0:
        return np.median(residues)
    else:
        return np.nan    

def get_z_score_v3(data,col,ctrl_col,controls):
    mean = data.loc[data[ctrl_col].isin(controls),col].mean()
    std = data.loc[data[ctrl_col].isin(controls),col].std()
    data.loc[:,str(col + ';z-score')] = data.loc[:,col].apply(lambda x: (x-mean)/std)
    return data

def get_most_severe_consequence(string):
    consequence_list = ['splice_acceptor_variant',
                        'splice_donor_variant',
                        'stop_gained',
                        'frameshift_variant',
                        'stop_lost',
                        'start_lost',
                        'transcript_amplification',
                        'inframe_insertion',
                        'inframe_deletion',
                        'missense_variant',
                        'protein_altering_variant',
                        'splice_region_variant',
                        'incomplete_terminal_codon_variant',
                        'start_retained_variant',
                        'stop_retained_variant',
                        'synonymous_variant',
                        'coding_sequence_variant',
                        'mature_miRNA_variant',
                        '5_prime_UTR_variant',
                        '3_prime_UTR_variant',
                        'non_coding_transcript_exon_variant',
                        'intron_variant',
                        'NMD_transcript_variant',
                        'non_coding_transcript_variant',
                        'upstream_gene_variant',
                        'downstream_gene_variant',
                        'TFBS_ablation',
                        'TFBS_amplification',
                        'TF_binding_site_variant',
                        'regulatory_region_ablation',
                        'regulatory_region_amplification',
                        'feature_elongation',
                        'regulatory_region_variant',
                        'feature_truncation',
                        'intergenic_variant']
    if type(string) == float:
        return 'None'
    type_list = string.split(',')
    for consequence in consequence_list:
        if consequence in type_list:
            return consequence

def make_stacked_bar(cat1,cat2,cat3,cat_labels,xlabels,colors=['r','b','g'],figsize=(2.5,2.5),width=0.35,linewidth=0):
    ###### USAGE EXAMPLE ######
    # cat1 = [0.2,0.7,0.1] # the first list you want to stack
    # cat2 = [0.3,0.1,0.5] # the second
    # cat3 = [0.5,0.2,0.4] # the third
    # cat_labels = ['cat1','cat2','cat3']
    # xlabels = ['A','B','C']
    # make_stacked_bar(cat1,cat2,cat3,cat_labels,xlabels)
    ###########################
    # This is written for 3 categories to stack (on y axis) and you have to modify if you want more or less
    # However, it can accept any number of items on x axis
    
    fig,ax = plt.subplots(figsize=figsize)
    xvals = np.arange(len(xlabels))
    p1 = plt.bar(x=xvals, height=cat1, width=width, color=colors[0],linewidth=linewidth,alpha=0.7)
    p2 = plt.bar(x=xvals, height=cat2, width=width, bottom=cat1,color=colors[1],linewidth=linewidth,alpha=0.7)
    p3 = plt.bar(x=xvals, height=cat3, width=width, bottom=[a+b for a,b in zip(cat1, cat2)],color=colors[2],linewidth=linewidth,alpha=0.7)
    plt.xticks(xvals, xlabels)
    plt.legend((p1[0], p2[0], p3[0]), (cat_labels[0],cat_labels[1],cat_labels[2]))
    sns.despine()
    plt.show()
    return fig,ax


def filter_pDNA_v2(scores):
    mean = scores['pDNA;lognorm'].mean()
    std = scores['pDNA;lognorm'].std()
    scores.loc[:,'pDNA_filter'] = scores.loc[:,'pDNA;lognorm'].apply(lambda x: (x > (mean - 3*std)) and (x < (mean + 3*std)))
    return scores

def filter_off_targets_v2(row):
    if row > 5:
        return False
    else:
        return True
    
def get_lfc_dropout(cell_line,condition,reps,df):
    for rep in reps:
        col = ';'.join([cell_line,condition,'Rep'+rep])
        df.loc[:,col+';LFC_dropout'] = df.loc[:,col+';lognorm'] - df.loc[:,cell_line+';Dropout;Rep'+rep+';lognorm']
    return df

def calc_lfc_from_pdna(df):
    cols = list(df.columns)
    cols = [x for x in cols if ('lognorm' in x)&('pDNA' not in x)]
    for c in cols:
        df.loc[:,c[:-8]+';LFC_pdna'] = df.loc[:,c] - df.loc[:,'pDNA;lognorm']
    return df

def calc_lognorm(df):
    cols = list(df.columns)
    cols = cols[1:]
    for c in cols:
        col_sum = np.sum(df[c])
        df.loc[:,c+';lognorm'] = [np.log2((x*1000000/float(col_sum))+1) for x in df[c]]
    return df

def get_colnames_from_supp_table(filename, sheetname):
    df = pd.read_excel(filename, nrows=2, sheet_name=sheetname)
    cols = ['Construct Barcode', 'pDNA'] 
    for i, r in df.T.iterrows():
        if (i != 'Construct') and (i != 'Cell line') and ('pDNA' not in i):
            colname = i.split('.')[0]+';'+ ';'.join(list(r))
            cols.append(colname)
    df = pd.read_excel(filename, header=2, sheet_name=sheetname)
    df.columns = cols
    return df