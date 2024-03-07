#!/usr/bin/env python

import cairo
import argparse
import re
import random

'''Write a python script using object-oriented code to visualize motifs on sequences.'''

def get_args(): 
    parser = argparse.ArgumentParser(description="Tool to visualize motifs on sequences.")
    parser.add_argument("-f", "--fasta", help="designates absolute file path to fasta file", type=str, required=False)
    parser.add_argument("-m", "--motifs", help="designates absolute file path to motifs file", type=str, required=False)
    return parser.parse_args()
    
args = get_args()
input_fasta = args.fasta
input_motifs = args.motifs


class Locate: 
    def __init__(self):
        '''Locate the features: motifs, exon, intron'''
        self.feature = None

    ## Methods ##
        '''Include the most basic qualities of each feature (DONT WORRY ABOUT THE DETAILS (EX: MOTIF COLOR))'''
    def motif(self): 
        pass

    def exon(self): 
        pass

    def intron(self):
        pass

class Draw: 
    surface = cairo.ImageSurface(cairo.FORMAT_RGB24, 1000, 600)
    context = cairo.Context(surface)
    context.set_source_rgb(1,1,1)
    context.paint()

    def __init__(self): 
        '''Draw the features: motifs, intron, exon'''
        self.feature = None

    ## METHODS ##
    def generate_random_color(self,colors): 
        while True: 
            r = round(random.random(), 2) 
            g = round(random.random(), 2)
            b = round(random.random(), 2)
            if (r,g,b) not in colors.values():
                return (r,g,b)
    
    def motif(self, start, end, pos_motif, r, g, b): 
        self.context.set_source_rgb(r,g,b)
        self.context.set_line_width(45)
        self.context.move_to(start,pos_motif)
        self.context.line_to(end,pos_motif)
        self.context.stroke()
        self.surface.write_to_png ("motif_mark.png")

    def intron(self,length, pos): 
        self.context.set_source_rgb(0, 0, 0)
        self.context.set_line_width(3)
        self.context.move_to(50,pos)
        self.context.line_to(length,pos)
        self.context.stroke()
        self.surface.write_to_png ("motif_mark.png")

    def exon(self,start,end,pos_exon): 
        self.context.set_source_rgb(0, 0, 0)
        self.context.set_line_width(45)
        self.context.move_to(start,pos_exon)
        self.context.line_to(end,pos_exon)
        self.context.stroke()
        self.surface.write_to_png ("motif_mark.png")



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
# loop through single_line.fasta and fill in key:value for data: data = {gene:sequence}
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
# print(motifs)


# Draws a horizontal line the length of the sequence
pos = 100
test = Draw()
for key in data: 
    length = len(data[key])
    # print(f'The sequence length is {length}')
    test.intron(length,pos)
    pos += 100


# Draw the exons
pos_exon = 100
for seq in data.values(): 
    matches = [(m.start(0), m.end(0)) for m in re.finditer('[A-Z+]+', seq)]
    start = matches[0][0] + 50
    end = matches[0][1] + 50
    test.exon(start,end, pos_exon)
    pos_exon += 100
    # print(start)
    # print(end)
    # print('\n')


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


# Draw the motifs 
pos_motif = 100
colors = {} # initializing an empty dictionary, color, to store the colors of each motif. colors = {'motif':(r,g,b)}
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
                    random_color = test.generate_random_color(colors)
                    colors[i]=random_color
                    test.motif(start,end,pos_motif,colors[i][0],colors[i][1],colors[i][2])
                else: 
                    test.motif(start,end,pos_motif,colors[i][0],colors[i][1],colors[i][2])
    pos_motif += 100
print(colors)
    