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
Date: 2019-01-31

Snakemake-linked script to draw stacked bar graphs for the 
taxonomic classifications of PZN and CAT, to compare the
output of the two methods.

INPUT:
  -  table of PZN classified superkingdoms (tab-separated)
  -  list of samples analysed by CAT
  -  list of tables of CAT classified contigs for each
     analysed sample
OUTPUT:
  -  grouped stacked bar graph of superkingdom classifications
     by PZN and CAT, quantified as number of reads (to
     facilitate comparison between classification possibilities)
"""

### STEP 1: Load libraries -------------------------------------------

import pandas as pd
from bokeh.core.properties import value
from bokeh.io import save, output_file
from bokeh.models import ColumnDataSource, FactorRange
from bokeh.plotting import figure

### STEP 2: Define global variables ----------------------------------

SAMPLES = snakemake.params["samples"]
COLOURS = [ "#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF" ]
TAXA = [ "Archaea", "Bacteria", "Eukaryota", "Viruses", "unclassified" ]

### STEP 3: Define functions -----------------------------------------

def sum_superkingdoms(classified_file):
    """
    Input:
        taxonomic classifications and quantifications for
         CAT (i.e. {sample}.CAT.Krona.info-monimal.tsv)
    Output:
        dataframe with the number of reads assigned to each
        superkindgom per sample
    """
    clas_df = pd.read_csv(classified_file, delimiter = '\t')
    clas_df["Sample_name"] = classified_file.split('.')[0]
    clas_df.fillna("not classified", inplace=True)
    #Count the reads assigned to Archaea, Bacteria, Eukaryota and Viruses per sample:
    superkingdom_sums = pd.DataFrame(clas_df.groupby(
    [ "Sample_name", "superkingdom" ]).sum()
                         [[ "reads" ]])

    #And make these into a proper dataframe:
    superkingdom_sums.reset_index(inplace=True)
    superkingdom_sums = pd.DataFrame(
        superkingdom_sums.pivot(index="Sample_name",
                                columns="superkingdom",
                                values="reads"))
    superkingdom_sums.reset_index(drop=True, inplace=True)
    
    for superkingdom in TAXA[:-1]:
        if superkingdom not in superkingdom_sums.columns:
            superkingdom_sums[superkingdom] = 0
        else:
            pass
    superkingdom_sums["Sample_name"] = classified_file.split('.')[0]
        
    return(superkingdom_sums)

def draw_stacked_bars(pzn_df, cat_df):
    """
    Draw grouped stacked bar charts
    Input:
        dataframe with pzn quantified classifications,
        dataframe with cat quantified classifications
    Output:
        interactive html figure with stacked bars
    """
    factors = []
    arch = []
    bact = []
    euk = []
    vir = []
    uncl = []

    for sample in SAMPLES:
        for tool in [ "PZN", "CAT" ]:
            factors.append((sample, tool))
            if tool == "PZN":
                arch.append(int(pzn_df[pzn_df["Sample_name"] == sample]["Archaea"]))
                bact.append(int(pzn_df[pzn_df["Sample_name"] == sample]["Bacteria"]))
                euk.append(int(pzn_df[pzn_df["Sample_name"] == sample]["Eukaryota"]))
                vir.append(int(pzn_df[pzn_df["Sample_name"] == sample]["Viruses"]))
                uncl.append(int(pzn_df[pzn_df["Sample_name"] == sample]["not classified"]))
            elif tool == "CAT":
                arch.append(int(cat_df[cat_df["Sample_name"] == sample]["Archaea"]))
                bact.append(int(cat_df[cat_df["Sample_name"] == sample]["Bacteria"]))
                euk.append(int(cat_df[cat_df["Sample_name"] == sample]["Eukaryota"]))
                vir.append(int(cat_df[cat_df["Sample_name"] == sample]["Viruses"]))
                uncl.append(int(cat_df[cat_df["Sample_name"] == sample]["not classified"]))
                
    source = ColumnDataSource(data=dict(
        x=factors,
        Archaea=arch,
        Bacteria=bact,
        Eukaryota=euk,
        Viruses=vir,
        unclassified=uncl
    ))
    
    fig = figure(x_range = FactorRange(*factors), plot_height = 600, plot_width = 1200,
                 title = "Profile of classified contigs per method, quantified as reads",
                 toolbar_location=None, tools = "hover, pan", tooltips = "@x $name: @$name reads")
    
    fig.vbar_stack(TAXA, x='x', width = 0.9, color = COLOURS, source = source,
                   legend = [value(x) for x in TAXA])
    
    fig.y_range.start = 0
    fig.x_range.range_padding = 0.1
    fig.xaxis.major_label_orientation = 1
    fig.xgrid.grid_line_color = None
    fig.legend.location = "top_center"
    fig.legend.orientation = "horizontal"

    output_file(snakemake.output[0])

    save(fig)
    
    return(None)

def main():
    """
    Main script function: does all the work from importing data to
    visualising the graph with Bokeh.
    """
    ## 1. import and select PZN input data:
    pzn_file = snakemake.input["pzn"]
    pzn_df = pd.read_csv(pzn_file)
    pzn_df.rename(columns={ "Unclassified" : "not classified",
                           "Sample" : "Sample_name" }, inplace=True)
    pzn_df = pzn_df.loc[ pzn_df["Sample_name"].isin(SAMPLES) ]
    
    ## 2. import CAT input data:
    cat_files = snakemake.input["cat"]
    cat_df = pd.DataFrame()
    for f in cat_files:
        superkingdom_sums = sum_superkingdoms(f)
        cat_df = pd.concat( [cat_df, superkingdom_sums], axis = 0,
                           join = "outer", sort = True )
    
    draw_stacked_bars(pzn_df, cat_df)
    
    return(None)

### STEP 4: Execute script -------------------------------------------

if __name__ == "__main__":
    main()
