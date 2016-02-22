#!/usr/bin/env python

import sys
import os.path
import argparse
import re


def get_cytobands(file):
	fin = open(file,'rt')
	chr_cytobands = {}

	for line in fin:
		line = line.strip()
		fields = line.split('\t')

		if fields[0] in chr_cytobands:
			chr_cytobands[fields[0]].append(fields)
		else:
			chr_cytobands[fields[0]] = []
			chr_cytobands[fields[0]].append(fields)

	fin.close()

	return chr_cytobands

def chrom_diagram(cytobands,fout):
	
	chr_width = 40
	centromere = 0
	scale_factor = 250000

	for band in cytobands:

		gpos_match = re.search(r'gpos',band[4])	
		acen_match = re.search(r'acen',band[4])


		if gpos_match:
			#print band
			start = eval(band[1])
			end = eval(band[2])
			length = end - start
			fout.write("<path class=\"locus\" d=\"M0 {0:.2f} h{1} v{2:.2f} h-{1} z\"\/>\n"
				.format(start/scale_factor, chr_width, length/scale_factor))

		elif acen_match and centromere == 0:
			centromere = band[2]




def check_files(args):
	files_not_found = 0;

	if not os.path.isfile(args.cytoband_file):
		print "Unable to find file: %s" % args.hgnc_file
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

def print_svg_footer(fout):	
	fout.write("</svg>\n")



def main(args):

	fout = open(args.out,'wt') if args.out != "" else sys.stdout

	chr_cytobands = get_cytobands(args.cytoband_file)
	cytobands = chr_cytobands["chr20"]

	chrom_diagram(cytobands,fout)

	#for bands in cytobands:
	#	print bands


	#print_svg_header(args.width,args.height,fout)

	#print_svg_footer(fout)	



if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument('--region', '-r', help='chr:start-end format, eg. 1:100-200', default="", required=False)
	parser.add_argument('--cytoband_file', '-cf', help='cytoband file', default="cytoBand.txt", required=False)

	parser.add_argument('--width', '-iw', help='SVG image height', default=1280)
	parser.add_argument('--height', '-ih', help='SVG image height', default=200)

	parser.add_argument('--out', '-o', help='Output SVG file', default="")

	args = parser.parse_args()

	check_files(args)
	main(args)
