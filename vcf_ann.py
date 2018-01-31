import argparse
import urllib2

parser = argparse.ArgumentParser(prog='VCF Annotation Tool', description='A program that parses a given VCF file and outputs a variant annotations table.')
parser.add_argument('input_file', help = 'Name of the vcf input file to be parsed')
parser.add_argument('output_file', help = 'Name of the output file created. VCF file suggested')
args = parser.parse_args()

f = open(args.output_file, "w")

var_type = "TYPE"
depth_of_seq = "DP"
reads_supp_var = "AO"
reads_supp_ref = "RO"
allele_freq = 'allele_freq"'
add_info = '"Consequence"'
gene = '"Gene"'

def find_info(field_name):
    indices = [i for i, j in enumerate(info_subfields) if field_name in j]
    #in case where RO is found in PRO use second index
    if field_name.startswith('RO'):
        i = indices[1]
    else:
        i = indices[0]
    x = info_subfields[i][len(field_name)+1:]
    #if multiple counts of variant reads are present, use first count
    #unless field is variant type
    if x.count(',') >= 1 and field_name != 'TYPE':
        reads = x.split(',')
        x = reads[0]
    return x
    
def percent_calc(x, y):
    var_reads = (int(x)/(int(y)+float(x)))*100
    ref_reads = (int(y)/(int(y)+float(x)))*100
    return str(round(var_reads, 2)) + ':' + str(round(ref_reads, 2))
    
def exac_file_search(keyword):
    u = urllib2.urlopen("http://exac.hms.harvard.edu/rest/variant/variant/" + chr + "-" + pos + "-" + ref + "-" + alt + "").read() 
    #searches through website to find first instance of field keyword,
    #pulls keyword between comma indeces
    if u.find(keyword) > 0:
        index = u.find(keyword)
        end = u.find(",", index, index+100)
        x = u[index+14:end]
        x = x.replace('_', ' ').replace('"', '')
    #if no keyword is determined by ExAC . placeholder is written
    else:
        x = "."
    return x
    
with open(args.input_file, "r") as file:
    print "Annotating files..."
    for line in file:
        if line[1] == "#": 
            #parses through meta-information of the input file
            continue
        if line[1] == "C" and line[0] == "#": 
            #once header line of input file is reached, print new header to output file
            f.write('{:10s} {:15s} {:10s} {:20s} {:25s} {:20s} {:20s}'.format('VAR TYPE', 'DEPTH SEQ COV', '# READS', '% READS (VAR VS ALT)', 'ALLELE FREQ', 'GENE', 'ADDITIONAL INFO'))
            continue
        else: 
            linelist = line.strip().split() 
            #defines each variable split in data line fields from line of input file
            (chr, pos, id, ref, alt, qual, filter, info) = linelist[0:8] 
            info_subfields = info.split(';')
            
            #1.writes type of variation to table in output file
            #type can be complex, del, ins, mnp, or snp ranked most to least serious
            vt = find_info(var_type)
            if vt.count(',') >= 1:
                types = vt.split(',')
                vt = types[0]
        
            #2.writes depth of sequence coverage to output file
            ds = find_info(depth_of_seq)
            
            #3.writes number of reads supporting variant to output file
            rsv = find_info(reads_supp_var)
            
            #4.writes % of supporting reads, variant vs reference to output file
            rsr = find_info(reads_supp_ref)
            pr = percent_calc(rsv, rsr)            
    
            #5.writes ExAC browser allele frequency to output file
            af = exac_file_search(allele_freq)
            
            #6.additional info, pull out 'Consequence' and 'Gene' of variant from exac file
            ai = exac_file_search(add_info)
            ge = exac_file_search(gene)
            
            f.write('{:10s} {:15s} {:10s} {:20s} {:25s} {:20s} {:20s}'.format(vt, ds, rsv, pr, af, ge, ai))
            
f.close()
