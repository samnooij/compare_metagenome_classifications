#! /usr/bin/env python

"""
    classification_comparison
    Copyright (C) 2019  Sam Nooij

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""

#Merge contig annotations from CAT classifications ("add_names.txt" file) and Jovian ("all_taxClassified.tsv" + "all_taxUnclassified.tsv").
# Primary goal: match reads numbers to contigs to quantify CAT profiles and generate Krona plots.

import pandas as pd

if __name__ == "__main__":
    #Open the files as Pandas dataframes:
    jovian_clas_df = pd.read_csv(snakemake.input["jovian_clas"], delimiter='\t')
    jovian_unclas_df = pd.read_csv(snakemake.input["jovian_unclas"], delimiter='\t')
    cat_clas_df = pd.read_csv(snakemake.input["names"], delimiter='\t')
    
    #Make sure the CAT dataframe has the same column name as the Jovian dataframe and that it is sorted correctly!
    cat_clas_df["# contig"] = cat_clas_df["# contig"].astype(str)
    cat_clas_df.rename(columns={ "# contig" : "#ID" }, inplace=True)
    cat_clas_df['sort'] = cat_clas_df['#ID'].str.extract('(\d+)', expand=False).astype(int)
    cat_clas_df.sort_values('sort',inplace=True, ascending=True)
    cat_clas_df = cat_clas_df.drop('sort', axis=1)
    cat_clas_df.reset_index(inplace=True, drop=True)
    
    #Concatenate Jovian dataframes into one long dataframe:
    jovian_contigs_df = pd.concat([jovian_clas_df, jovian_unclas_df], sort=True)
    
    #Calculate the number of reads mapped to each contig:
    jovian_contigs_df["reads"] = jovian_contigs_df.Plus_reads + jovian_contigs_df.Minus_reads
    minimal_jovian_df = jovian_contigs_df[["#ID", "reads"]]
    minimal_jovian_df.columns = [ "#ID", "reads" ]
    
    #Merge Jovian annotations to CAT classification:
    merged_df = minimal_jovian_df.merge(cat_clas_df, left_on="#ID", right_on="#ID", sort=True)
    
    #Remove unnecessary columns:
    merged_df.drop(columns=["classification", 
                            "number of ORFs on contig", 
                            "number of ORFs classification is based on", 
                            "lineage scores", "lineage"], inplace=True)
    
    #Remove the fractions listed with each taxonomic rank:
    for column in [ "superkingdom", "phylum", "class", "order", "family", "genus", "species" ]:
        #Remove the trailing ": [number]" from CAT's output, e.g. from "Bacteria: 1.00" to "Bacteria":
        merged_df[column].replace(to_replace=':\s[\d].+', value='', regex=True, inplace=True)
    
    #Write the resulting dataframe to a file:
    merged_df.to_csv(snakemake.output[0], sep='\t', index=False)
