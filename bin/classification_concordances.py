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


Author: Sam Nooij
Date: 2019-01-22

Script goal:
Count the number of concordantly classified contigs (and later also reads?)
between classification algorithms. Visualise the results as (1) tables, 
(2) bar graphs, (3) Venn diagrammes, and (4) lists of taxa by each/either 
method.

Ideally this script should be able to:
 - detect if the input is from 1 or more samples to generate sample-specific 
   statistics or "overall" stats
 - count concordances between the classification methods
 - make tables and figures with those concordance-numbers
 - draw Venn diagrammes

INPUT:
 - tab-separated tables of metagenomics classifications divided in taxonomic
   lineages, e.g. a table with columns "superkingdom", "phylum", ..., "genus",
   "species" AND a column with contig names (to compare contig-to-contig 
   classifications by e.g. PZN and CAT) 
   
   ## OR -pehaps a later addition- a column with number of reads 
   classified as the given taxon (to compare read-based classifications to 
   e.g. MetaMeta)

OUTPUT:
 - tables of classification concordance (absolute and relative=percentages)
 - bargraphs of classification concordance (absolute and relative)
 - Venn diagram of classification concordance (absolute, per taxon)
 - lists of taxa found by either method, or both(, or with 3 methods?)
"""

### STEP 1: Load libraries -------------------------------------------

import pandas as pd
from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
from bokeh.core.properties import value
from bokeh.io import save, output_file
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Tabs, Panel
from bokeh.plotting import figure
from itertools import zip_longest as zl #for writing tables of unequal length

### STEP 2: Define global variables ----------------------------------

RANKS = [ "superkingdom", "phylum", "class", 
         "order", "family", "genus", "species" ]

### STEP 3: Define functions -----------------------------------------

def concat_files(file_list):
    """
    input: list of (tab-separated) files
    output: concatenated pandas dataframe
    """
    df_list = []
    for file in file_list:
        df = pd.read_csv(file, delimiter='\t')
        df_list.append(df)
        
    return(pd.concat(df_list, sort=True))

def find_consistencies(df):
    """
    input: merged dataframe with classifications per method per taxonomic rank
    output: dataframe with numbers of (per rank!):
        1. number of contigs classified by each method
        2. number of contigs consistently classified by both methods
        3. number of inconsistently classified contigs (i.e. both methods classified the contig, but to different taxa - one says 'Viruses', the other 'Bacteria')
        4. number of contigs consistently UNclassified
    """
        
    concordances = { "Rank" : [],
                   "classified_contigs" : [],
                   "consistently_classified" : [],
                   "inconsistently_classified" : [],
                   "consistently_unclassified" : [],
                   "only_pzn" : [],
                   "only_cat" : []
                   }
    
    for rank in RANKS:
        concordances["Rank"].append(rank)
        #If both methods did not report the contig as "not classified", it was classified by both methods:
        classified_contigs = sum((df["%s_pzn" % rank] != "not classified") & (df["%s_cat" % rank] != "not classified"))
        concordances["classified_contigs"].append(classified_contigs)
        #If one of the methods (PZN) did not report "not classified" and 
        # the classifications of PZN and CAT are identical, the contig has been consistently classified:
        consistently_classified = sum((df["%s_pzn" % rank] != "not classified") & (df["%s_pzn" % rank] == df["%s_cat" % rank]))
        concordances["consistently_classified"].append(consistently_classified)
        #The inconsistently classified contigs are those that were classified by both methods, 
        # but with unequal taxa (total classified - consistent):
        inconsistently_classified = classified_contigs - consistently_classified
        concordances["inconsistently_classified"].append(inconsistently_classified)
        #Then there are contigs that PZN could classify, but CAT not:
        only_pzn = sum((df["%s_pzn" % rank] != "not classified") & (df["%s_cat" % rank] == "not classified"))
        concordances["only_pzn"].append(only_pzn)
        #And contigs that only CAT could classify, and PZN not:
        only_cat = sum((df["%s_pzn" % rank] == "not classified") & (df["%s_cat" % rank] != "not classified"))
        concordances["only_cat"].append(only_cat)
        #And the contigs for which both methods report "not classified" are consistently unclassified:
        consistently_unclassified = sum((df["%s_pzn" % rank] == "not classified") & (df["%s_pzn" % rank] == df["%s_cat" % rank]))
        concordances["consistently_unclassified"].append(consistently_unclassified)
        
    ## THANKS askewchan (https://stackoverflow.com/a/16343791) for telling about the use of '&' with conditionals in dataframes!
        
    return(pd.DataFrame.from_dict(concordances))

def plot_bars(df, perc, outfile, total):
    """
    input: dataframe with ranks and number of concordantly classified contigs,
        name for an output file, total number of contigs (as input for analysis),
        number of contigs classified by both/all methods, and whether or not
        to create a plot of percentages
    output: bar graph as html file (saved as 'outfile' -> see input)
    """
    output_file(outfile)

    def find_suitable_y(n):
        """
        Find a proper upper limit for the y-axis, based on the maximum value to plot.
        """
        if n > 20000:
            return(n) #just use the maximum value in case the numbers get really big
        elif n > 10000:
            return(20000)
        elif n > 5000:
            return(10000)
        elif n > 2000:
            return(5000)
        elif n > 1000:
            return(2000)
        elif n > 500:
            return(1000)
        elif n > 200:
            return(500)
        elif n > 100:
            return(200)
        else:
            return(100) #100 will probably do as lower limit
    
    nr_source = ColumnDataSource(data=dict(ranks = df["Rank"],
                                           classified_contigs = df["classified_contigs"],
                                           consistently_classified = df["consistently_classified"],
                                           inconsistently_classified = df["inconsistently_classified"],
                                           consistently_unclassified = df["consistently_unclassified"],
                                           only_pzn = df["only_pzn"],
                                           only_cat = df["only_cat"]
                                          )
                                )
        
    title = "Classified contigs comparison (PZN-CAT)"
    ymax = find_suitable_y(total)

    colors = [ "#0000A6","#63FFAC","#B79762","#004D43", "#FFDBE5" ]
    parts = [ "consistently_classified", "inconsistently_classified",
             "only_pzn", "only_cat", "consistently_unclassified" ]

    nr_fig = figure(x_range=df["Rank"], plot_height = 600, plot_width = 1000,
                   title = title, toolbar_location=None, tools = "hover, pan",
                   tooltips="@ranks $name: @$name")
    
    nr_fig.vbar_stack(parts, x='ranks', width=0.9, color=colors, source=nr_source,
                 legend=[value(x) for x in parts])

    nr_fig.xgrid.grid_line_color = None
    nr_fig.legend.orientation = "horizontal"
    nr_fig.legend.location = "top_center"
    
    nr_panel = Panel(child=nr_fig, title="Absolute numbers of contigs")
    
    perc_source = ColumnDataSource(data=dict(ranks = perc["Rank"],
                                           classified_contigs = perc["classified_contigs"],
                                           consistently_classified = perc["consistently_classified"],
                                           inconsistently_classified = perc["inconsistently_classified"],
                                           consistently_unclassified = perc["consistently_unclassified"],
                                           only_pzn = perc["only_pzn"],
                                           only_cat = perc["only_cat"]
                                          )
                                )
    
    perc_fig = figure(x_range=perc["Rank"], plot_height = 600, plot_width = 1000,
                   title = title, toolbar_location=None, tools = "hover, pan",
                   tooltips="@ranks $name: @$name %")
    
    perc_fig.vbar_stack(parts, x='ranks', width=0.9, color=colors, source=perc_source,
                 legend=[value(x) for x in parts])

    perc_fig.xgrid.grid_line_color = None
    perc_fig.legend.orientation = "horizontal"
    perc_fig.legend.location = "top_center"
    
    perc_panel = Panel(child=perc_fig, title="Percentages of contigs")
    
    tabs = Tabs(tabs=[nr_panel, perc_panel])
    
    save(tabs)
    
    return(None)

def find_taxa(df):
    """
    List all taxa identified by the different methods, per taxonomic rank.
    input: dataframe with all classifications
    output: dictionary with deduplicated lists of taxa per rank
    """
    set_dict = {}
    for column in df.columns:
        if column != "#ID" and column.split('_')[0] in RANKS and "concordance" not in column:
            #exclude "not classified" from this list:
            taxa_list = set(df[column])
            taxa_list.discard("not classified")
            set_dict[column] = taxa_list
        else:
            pass
        
    return(set_dict)

def draw_venn(taxa_dict):
    """
    Make Venn diagrammes for the taxa identified by different methods (i.e. PZN and CAT).
    input: dictionary with list of taxa classified per taxonomic rank, per method
    output: Venn diagrammes as png file
    """
    plt.rcParams['figure.figsize'] = [48, 24]
    #Setting plot size tip by tacaswell, see https://stackoverflow.com/a/36368418
    plt.rcParams.update({'font.size': 24})
    #Increase font size with Herman Schaaf's help (https://stackoverflow.com/a/3900167)

    #Make subplots as shown by Jarad (https://stackoverflow.com/a/39133654)
    fig = plt.figure()
    coordinates = [ 241, 242, 243, 244, 245, 246, 247 ]
    rank_plurals = [ "Superkingdoms", "Phyla", "Classes",
                    "Orders", "Families", "Genera", "Species" ]
    ### DEBUG ###                
    #print(taxa_dict)
    
    for i in range(len(RANKS)):
    #Make Venn diagrammes for each rank as subplot (7 figures in 1 file)
        ax = fig.add_subplot(coordinates[i])
        diag = venn2([taxa_dict["%s_pzn" % RANKS[i]], taxa_dict["%s_cat" % RANKS[i]]], set_labels = ('PZN', 'CAT'))
        diag.get_patch_by_id('010').set_color('skyblue')
        diag.get_patch_by_id('110').set_color('palegoldenrod')
        diag.get_patch_by_id('100').set_color('sienna')
        ax.title.set_text("%s classified" % rank_plurals[i])
        
    plt.savefig(snakemake.output["venn"], dpi=200)
    plt.show()
    
    return(None)

def write_tables(taxa_dict):
    """
    Write the taxa that have been classified by different methods, and the overlap between the two
    to tab-separated text files.
    input: dictionary with taxa per rank per method
    output: tab-separated files (e.g. {sample}.PZN-CAT.comparison.species.tsv)
    """
    for rank in RANKS:
        overlap = sorted(list(set(taxa_dict["%s_pzn" % rank]) & set(taxa_dict["%s_cat" % rank])))
        #overlap thanks to Mark Byers: https://stackoverflow.com/a/3697438
        unique_pzn = sorted(list(set(taxa_dict["%s_pzn" % rank]) - set(taxa_dict["%s_cat" % rank])))
        unique_cat = sorted(list(set(taxa_dict["%s_cat" % rank]) - set(taxa_dict["%s_pzn" % rank])))
        #uniques thanks to Javed: https://stackoverflow.com/a/47264446
    
        with open("%s%s.tsv" % (snakemake.params["list_prefix"], rank), 'w') as outfile:
            outfile.write("Overlap\tPZN\tCAT\n")
            for a, b, c in zl(overlap, unique_pzn, unique_cat, fillvalue=''):
                outfile.write("{0:20s}\t{1:20s}\t{2:20s}\n".format(a, b, c))
                
    return(None)

def main(cat_df):
    """
    Main execution of the script: does all the collection of data from files,
    counting of contigs and taxonomic classifications and concordances
    between PZN and CAT.
    Input:
        Pandas DataFrame with CAT classifications per sample of concatenated
        samples for an "OVERALL" comparison.
    Output:
        1. Table with classifications comparisons (counts of similarities/
        differences between the classification of contigs on different
        taxonomic ranks of PZN (Megablast) and CAT (Diamond).
        2. Table with classification comparisons as percentages.
        3. Stacked bar graph visualising the comparisons from the table.
        (both in absolute numbers and in percentages)
        4. Venn diagrammes that show the overlap in assigned taxa
        between PZN and CAT at different taxonomic ranks.        
    """
    pzn_df = concat_files(snakemake.input["pzn"])
    
    merged_df = pzn_df.merge(cat_df, on='#ID', how="inner", suffixes=('_pzn', '_cat'), sort=True)
    merged_df.fillna("not classified", inplace=True)

    total_contigs = len(merged_df)
    ## DEBUG:
    #print(merged_df.info(), merged_df.head())

    #Now find the list of taxa that have been classified by each method:
    taxa_per_method_dict = find_taxa(merged_df)
    #Draw Venn diagrammes using these lists of taxa:
    draw_venn(taxa_per_method_dict)
    #And finally write these to tab-separated files:
    write_tables(taxa_per_method_dict)

    #Count the number of consistently (un)classified contigs per rank:
    concordance_df = find_consistencies(merged_df)

    #Export the number of concordances to a file:
    concordance_df.to_csv(snakemake.output["table_nr"], sep='\t', index = False)

    #Calculate percentages:
    concordance_percentages = pd.DataFrame()
    concordance_percentages["Rank"] = concordance_df["Rank"]
    for column in concordance_df.columns:
        if column != "Rank":
            concordance_percentages[column] = concordance_df[column] / total_contigs * 100
        else:
            pass
    #Also export percentages to a file:
    concordance_percentages.to_csv(snakemake.output["table_pc"], sep='\t', index = False)

    #And plot the concordances in a bar graph:
    plot_bars(df = concordance_df,
              perc = concordance_percentages,
              outfile = snakemake.output["figure_nr"], 
              total = total_contigs) #provide the total number of contigs used in the analysis = length of the dataframe

    return(None)
                               
### STEP 4: Execute script -------------------------------------------

if __name__ == "__main__":
    if isinstance(snakemake.input["cat"], str):
#OPTION 1: single-sample comparison
        cat_df = pd.read_csv(snakemake.input["cat"], delimiter='\t')
        main(cat_df)
        
    elif isinstance(snakemake.input["cat"], list):
#OPTION 2: multi-sample "overall" comparison
        cat_df = concat_files(snakemake.input["cat"])
        main(cat_df)
            
    else:
#OPTION 3: invalid input (need one file or a list of files)
        raise InputError("Need a file (as string) or multiple files (as list) as input from Snakemake.")
