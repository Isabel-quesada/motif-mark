#!/usr/bin/env python

import cairo
import argparse
import re
import random

'''Python script using object-oriented code to visualize motifs on sequences.'''

def get_args(): 
    parser = argparse.ArgumentParser(description="Tool to visualize motifs on sequences.")
    parser.add_argument("-f", "--fasta", help="designates absolute file path to fasta file", type=str,required=True)
    parser.add_argument("-m", "--motifs", help="designates absolute file path to motifs file", type=str, required=True)
    return parser.parse_args()

args = get_args()
input_fasta = args.fasta
input_motifs = args.motifs


# Getting the prefix of the input fasta file to use when generating the ouput png. This accounts for absolute paths as input and 
# also for when the input files are in the working directory.
filepath = vars(args)['fasta'] 
file_name = "" # variable to hold the prefix of the input fasta file

if filepath.startswith('/'): # absolute path was given for input fasta file
    prefix = filepath.split('/')
    prefix.reverse()
    file_name = prefix[0].split('.')[0]
else: 
    # input files are in the same directory as motif-mark-oop.py
    file_name = filepath.split('.')[0]


class Locate: 
    '''Locate the features: motifs, exon, intron, gene name'''
    def __init__(self):
        self.feature = None

    ## Methods ##
    def find_motif(self,data,regex_motifs,file_name): 
        '''Find all instances of input motifs in each sequence. Assign each unique motif a random (r,g,b) color.'''
        pos_motif = 200
        for seq in data.values(): 
            sequence = seq.upper()
            for i in regex_motifs: 
                # i = motif
                pattern = regex_motifs[i] #regex pattern for the motif
                matches = [(m.start(0), m.end(0)) for m in re.finditer(pattern, sequence)]
                if matches == []: 
                    # motif is not present in the sequence
                    pass
                else: 
                    for pos in matches: 
                        # pos = (start,end)
                        start = pos[0] + 50
                        end = pos[1] + 50
                        if i not in colors.keys(): 
                            random_color = figure.generate_random_color(colors)
                            colors[i]=random_color
                            figure.motif(start,end,pos_motif,colors[i][0],colors[i][1],colors[i][2],file_name)
                        else: 
                            figure.motif(start,end,pos_motif,colors[i][0],colors[i][1],colors[i][2],file_name)
            pos_motif += 150

    def find_intron(self,data,file_name):
        '''Determine the length of the sequence'''
        pos = 200
        for key in data: 
            length = len(data[key]) + 50
            figure.intron(length, pos,file_name)
            pos += 150
    
    def find_exon(self,data,file_name): 
        '''Find the exon/exons present in the sequence'''
        pos_exon = 200
        for seq in data.values(): 
            matches = [(m.start(0), m.end(0)) for m in re.finditer('[A-Z+]+', seq)]
            start = matches[0][0] + 50
            end = matches[0][1] + 50
            figure.exon(start, end, pos_exon,file_name)
            pos_exon += 150

    def gene_name(self,data,file_name):
        '''Get the gene name for each sequence'''
        pos_gene = 150
        for key in data: 
            figure.gene_name(key,pos_gene,file_name)
            pos_gene += 150

class Draw: 
    '''Draw the features: gene name, legend, motifs, intron, exon'''
    surface = cairo.ImageSurface(cairo.FORMAT_RGB24, 1500, 3000)
    context = cairo.Context(surface)
    context.set_source_rgb(1,1,1)
    context.paint()

    def __init__(self): 
        self.feature = None

    ## METHODS ##
    def generate_random_color(self,colors): 
        '''Generates random (r,g,b) to be used as the motif colors. Ensures the same color is not used for different motifs.'''
        while True: 
            r = round(random.random(), 2) 
            g = round(random.random(), 2)
            b = round(random.random(), 2)
            if (r,g,b) not in colors.values():
                return (r,g,b)
    
    def gene_name(self,gene,pos_gene,file_name):
        '''Add gene name for each sequence to the figure'''
        self.context.set_source_rgb(0,0,0)
        self.context.set_font_size(20)
        self.context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        self.context.move_to(50,pos_gene)
        self.context.show_text(f'Gene {gene}')
        self.context.stroke()
        self.surface.write_to_png(f'{file_name}.png')

    def legend(self,colors,file_name,x=10,y=10,width=20,height=10,spacing=10,font_size=12): 
        '''Add a legend to show the color of each motif'''
        self.context.set_font_size(font_size)
        for motif, color in colors.items():
            r = color[0]
            g = color[1]
            b = color[2]
            self.context.set_source_rgb(r,g,b)
            self.context.rectangle(x,y,width,height)
            self.context.fill()
            self.context.move_to(x+width+spacing,y+height/2 + font_size/2)
            self.context.set_source_rgb(0,0,0)
            self.context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            self.context.show_text(motif)
            y+= height+spacing
            self.surface.write_to_png(f'{file_name}.png')
        
    def motif(self,start,end,pos_motif,r,g,b,file_name): 
        '''Draws motif to scale'''
        self.context.set_source_rgb(r,g,b)
        self.context.set_line_width(45)
        self.context.move_to(start,pos_motif)
        self.context.line_to(end,pos_motif)
        self.context.stroke()
        self.surface.write_to_png (f'{file_name}.png')

    def intron(self,length,pos,file_name): 
        '''Draws intron to scale'''
        self.context.set_source_rgb(0,0,0)
        self.context.set_line_width(3)
        self.context.move_to(50,pos)
        self.context.line_to(length,pos)
        self.context.stroke()
        self.surface.write_to_png (f'{file_name}.png')

    def exon(self,start,end,pos_exon,file_name): 
        '''Draws exon to scale'''
        self.context.set_source_rgb(0,0,0)
        self.context.set_line_width(45)
        self.context.move_to(start,pos_exon)
        self.context.line_to(end,pos_exon)
        self.context.stroke()
        self.surface.write_to_png (f'{file_name}.png')


# one line fasta file
seq = ""
with open(input_fasta, "r") as fh, open("single_line.fasta", "w") as out:
    '''Put all the sequence lines onto one line. Now the fasta will have 2 line records (one header, one sequence line)''' 
    first_line = True
    for line in fh:
        if first_line == True: 
            out.write(line)
            first_line = False
        elif first_line == False and not line.startswith('>'): 
            seq+=line.strip('\n')
        else: 
            out.write(seq+"\n")
            out.write(line)
            seq = ""
    out.write(seq)


# initializing empty dictionary
data = {}
# loop through single_line.fasta and fill in key:value for data: data = {'gene':'sequence'}
with open("single_line.fasta","r") as fh: 
    while True: 
        header = fh.readline().strip()
        if header == "":
            #EOF
            break
        x = re.split(r'([^>]\S+)', header)
        x = list(filter(None, x))   # example output for x: ['>', 'INSR', ' chr19:7150261-7150808', ' (reverse', ' complement)']
        gene = x[1]
        sequence = fh.readline().strip()
        data[gene]=sequence


# initializing empty list 
motifs = []
# loop through input motifs and fill in values for motifs
with open(input_motifs, "r") as fh: 
    while True: 
        line = fh.readline().strip().upper()
        if line == "": 
            # EOF
            break
        motifs.append(line)


# IUPAC degenerate base symbols stored in dictionary iupac = {'Symbol': 'regex pattern for bases represented'}
iupac = {'A':'[A]','C':'[C]','G':'[G]', 'T':'[T]', 'U':'[U]', 'W':'[AT]', 'S':'[CG]', 'M':'[AC]', 'K':'[GT]', 
                'R':'[AG]', 'Y':'[CT]', 'B':'[CGT]', 'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]', 'N':'[ACGT]'}


# input motifs and the regex pattern for them are stored in the dictionary regex_motifs = {'motif': 'regex pattern'}
regex_motifs = {}
for i in motifs: 
    pattern = ''
    for base in i: 
        value = iupac[base]
        pattern+=value
    regex_motifs[i] = pattern
    pattern = ''

colors = {} # initializing an empty dictionary, color, to store the colors of each motif. colors = {'motif':(r,g,b)}


# Creating objects of both classes
figure = Draw()
search = Locate()

search.find_intron(data,file_name)
search.find_exon(data,file_name)
search.find_motif(data,regex_motifs,file_name)
search.gene_name(data,file_name)
figure.legend(colors,file_name)