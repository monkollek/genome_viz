#!/usr/bin/env python

import sys
import os.path
import argparse
import re
import gzip
import pysam

def get_junction_reads(file,chr,start,end):

	junction_reads = {}
	junctions1 = re.compile(r'^([0-9]+)M([0-9]+)N([0-9]+)M$')
	#Need to do two junctions!

	samfile = pysam.AlignmentFile(file, "rb")

	for read in samfile.fetch(chr, start, end):

	    if read.cigarstring:
	 		#match = re.search(r'^([0-9]+)M([0-9]+)N([0-9]+)M$', read.cigarstring)
	 		match = junctions1.search(read.cigarstring)

	 		if match:
	 			junction_start = read.reference_start + eval(match.group(1))
	 			junction_end = junction_start + eval(match.group(2))
	 			junction_key = "%d_%d" % (junction_start, junction_end)

	 			if junction_key in junction_reads:
	 				junction_reads[junction_key] += 1
	 			else:
	 				junction_reads[junction_key] = 1



	samfile.close()

	return junction_reads

def get_gene_positions(file):
	fin = gzip.open(file,'rt')
	gene_positions = {}

	for line in fin:
		line = line.strip()
		fields = line.split('\t')
		gene_positions[fields[3]] = fields


	fin.close()

	return gene_positions



def exon_track(file, chr, start, end, gene_name, scale, x_offset, y_offset, fout):
	
	gencode_file = pysam.TabixFile(file)

	for tabix_line in gencode_file.fetch(chr, start, end):
		fields = tabix_line.split('\t')
		pos = (eval(fields[1]) - x_offset)/scale
		width = (eval(fields[2]) - eval(fields[1]))/scale

		fout.write("<rect class=\"exon_rect\" x=\"%.1f\" y=\"%d\" width=\"%.1f\" height=\"20\" />\n" % (pos, y_offset, width))

	fout.write("<line class=\"exon_line\" x1=\"{0:.1f}\" y1=\"{1}\" x2=\"{2:.1f}\" y2=\"{1}\" />\n"
		.format((start-x_offset)/scale, y_offset+10, (end-x_offset)/scale))
	fout.write("<text class=\"gene_text\" x=\"{0:.1f}\" y=\"{1}\">{2}</text>\n"
		.format((start-x_offset)/scale, y_offset+40, gene_name))


def ruler_track(chr, start, end, scale, x_offset, y_offset, fout):

	region_width = end - start

	if region_width < 1000:
		tick_interval = 100
		tick_text_interval = 100
	elif region_width < 10000:
		tick_interval = 1000
		tick_text_interval = 1000
	elif region_width < 100000:
		tick_interval = 10000
		tick_text_interval = 10000		
	else:
		tick_interval = 10000
		tick_text_interval = 10000

	#Ruler line
	fout.write("<line class=\"ruler_line\" x1=\"{0:.1f}\" y1=\"{1}\" x2=\"{2:.1f}\" y2=\"{1}\" />\n"
		.format((start-x_offset)/scale, y_offset, (end-x_offset)/scale))

	tick_at = start - start%tick_interval + tick_interval 
	
	#First tick
	tick_at_pos = (tick_at-x_offset)/scale
	fout.write("<text class=\"tick_text\" x=\"{0:.1f}\" y=\"{1}\">{2:,}</text>\n".format(tick_at_pos,y_offset-10,tick_at))

	while tick_at < end:
		tick_at_pos = (tick_at-x_offset)/scale

		fout.write("<line class=\"ruler_line\" x1=\"{0:.1f}\" y1=\"{1}\" x2=\"{0:.1f}\" y2=\"{2}\" />\n"
			.format(tick_at_pos,y_offset+5,y_offset-5))

		if tick_at%tick_text_interval == 0:
			fout.write("<text class=\"tick_text\" x=\"{0:.1f}\" y=\"{1}\">{2:,}</text>\n"
				.format(tick_at_pos,y_offset-10,tick_at))

		tick_at = tick_at + tick_interval 


	#Last tick	
	tick_at = tick_at - tick_interval
	if tick_at%tick_text_interval != 0:
		tick_at_pos = (tick_at-x_offset)/scale
		fout.write("<text class=\"tick_text\" x=\"{0:.1f}\" y=\"{1}\">{2:,}</text>\n"
			.format(tick_at_pos,y_offset-10,tick_at))


def splicing_track(chr, start, end, junction_reads, scale, x_offset, y_offset, fout):
	read_threshold = 10
	arc_length = 25

	for junc_key,reads in junction_reads.iteritems():
		junc_fields = junc_key.split('_')
		#print junc_fields[0], " ", junc_fields[1], " ", reads

		if reads > read_threshold:
			arc_start = (eval(junc_fields[0]) - x_offset)/scale
			arc_end = (eval(junc_fields[1])-x_offset)/scale - arc_start
			fout.write("<path class=\"junc_arc1\" d=\"M {0:.1f} {1} c 0 -{2} {3:.1f} -{2} {3:.1f} 0\" />\n"
				.format(arc_start, y_offset, arc_length, arc_end))
			fout.write("<text class=\"junc_arc_text1\" x=\"{0:.1f}\" y=\"{1}\">{2}</text>\n"
				.format(arc_start + (arc_end)/2, y_offset - 25, reads))

		#print "%s %d" % (junc_key, reads)
		#print "Junction: %s Reads: %d" % (k, v)

def calculate_offset(region_chr, region_start, region_end, scale):

	region_width = region_end - region_start

	if region_width < 10000:
		interval_length = 100
	elif region_width < 500000:
		interval_length = 1000		
	else:
		interval_length = 100000

	if region_start%interval_length/scale < 10:
		offset = region_start - region_start%interval_length - 5*scale
	else:
		offset = region_start - region_start%interval_length

	#print "Width: %d offset: %d" % ((region_end - region_start),offset)

	return offset


def find_gene(file, chr, start, end):

	hgnc_file = pysam.TabixFile(file)
	gene = "unknown"

	for tabix_line in hgnc_file.fetch(chr, start, end):
		fields = tabix_line.split('\t')
		gene = fields[3]

	return gene

def check_files(args):
	files_not_found = 0;

	if not os.path.isfile(args.hgnc_file):
		print "Unable to find file: %s" % args.hgnc_file
		files_not_found = 1

	if not os.path.isfile(args.gencode_bed_file):
		print "Unable to find file: %s" % args.gencode_bed_file
		files_not_found = 1

	if not os.path.isfile(args.bam):
		print "Unable to find file: %s" % args.bam
		files_not_found = 1


	if files_not_found == 1:
		sys.exit()



def print_svg_header(width, height, fout):
	fout.write("<?xml version=\"1.0\"?>\n")
	#fout.write("<?xml-stylesheet type=\"text/css\" href=\"svg.css\" ?>\n")
	fout.write("<svg width=\"%d\" height=\"%d\" xmlns=\"http://www.w3.org/2000/svg\">\n" % (width,height))

	fout.write("<style type=\"text/css\">\n")
	css_file = open(args.css_file)

	for line in css_file:
		fout.write(line)

	fout.write("</style>\n")	



	#print css
	#print "in here"

def print_svg_footer(fout):	
	fout.write("</svg>\n")


def main(args):

	all_gene_positions = get_gene_positions(args.hgnc_file)
	gene = args.gene

	if args.region:

		if args.gene:
			print "Region specified. Ignoring --gene %s option" % args.gene

		match = re.search(r'^(?P<region_chr>[0-9XY]+):(?P<region_start>[0-9]+)-(?P<region_end>[0-9]+)$',args.region)	

		if match:
			#print "%s %s %s" % (match.group('region_chr'),match.group('region_start'),match.group('region_end'))
			region_chr = match.group('region_chr')
			region_start = eval(match.group('region_start'))
			region_end = eval(match.group('region_end'))
			gene = find_gene(args.hgnc_file,region_chr,region_start,region_end)

		else:
			print "Genomic region in wrong format"
			sys.exit()

	elif gene in all_gene_positions:	
		gene_position = all_gene_positions[gene]
		region_chr = gene_position[0]
		region_start = eval(gene_position[1])
		region_end = eval(gene_position[2])
	
	else:
		print "Gene: %s not found!" % (gene)
		sys.exit()


	fout = open(args.out,'wt') if args.out != "" else sys.stdout

	svg_scale = (region_end - region_start) / 1200.0
	x_offset = calculate_offset(region_chr, region_start, region_end, svg_scale)

	print_svg_header(args.width,args.height,fout)

	exon_track(args.gencode_bed_file,region_chr, region_start, region_end, gene, svg_scale, x_offset, 150, fout)
	ruler_track(region_chr,region_start,region_end,svg_scale,x_offset,100,fout)

	junc_reads = get_junction_reads(args.bam,region_chr, region_start, region_end)
	splicing_track(region_chr, region_start, region_end, junc_reads, svg_scale, x_offset, 150,fout)

	print_svg_footer(fout)


#/Applications/Inkscape.app/Contents/Resources/bin/inkscape -b "white" -z -e  /Users/monkol/dev/visualizations/test.png /Users/monkol/dev/visualizations/test4.svg
#./ngs_viz.py -r 19:39062658-39078204 -o test4.svg -b GTEX-N7MS_muscle_RYR1.bam

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument('--bam', '-b', help='BAM File', required=True)
	parser.add_argument('--gene', '-g', help='HGNC Gene Name', required=False)
	parser.add_argument('--region', '-r', help='chr:start-end format, eg. 1:100-200', default="", required=False)
	parser.add_argument('--hgnc_file', '-hf', help='hgnc file', default="./resources/hgnc_anno.txt.gz", required=False)
	parser.add_argument('--gencode_bed_file', '-gbf', help='Gencode BED file', default="./resources/gencode.v19.exon.bed.gz", required=False)
	parser.add_argument('--css_file', '-cf', help='CSS file', default="./resources/svg.css", required=False)

	parser.add_argument('--width', '-iw', help='SVG image height', default=1280)
	parser.add_argument('--height', '-ih', help='SVG image height', default=200)


	parser.add_argument('--out', '-o', help='Output SVG file', default="")

	args = parser.parse_args()

	check_files(args)
	main(args)


	
