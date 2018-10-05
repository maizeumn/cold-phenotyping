#! /usr/bin/perl -w

# generate_qrcode_R_script.pl

# This script will take a tab delimited file and print the corresponding R code
# to generate the necessary qr code page to use with phenotyping pictures for sample tracking.

# Cory Hirsch
# Feb. 17 2016

use strict;
use Getopt::Long();

my $usage = "\nUsage: $0 --in <tab delimented file with information for qrcode page> --out <output file containing R code to generate sample tracking page> --help <help on usage of script>\n\n";

my ($in, $out, $help);

Getopt::Long::GetOptions('in=s' => \$in,
			       'out=s' => \$out
				 );

if (defined($help)) {
	die $usage;
}

# Open in file and output file and print R Code to generate qrcode sample tracking page
open (my $in_fh, '<', $in) or die "Can't open file $in\n\n";
open (my $out_fh, '>', $out) or die "Can't open file $out\n\n";

# Skip header line
my $header = <$in_fh>;

# Print initial R environment Setup
# Load necessary R libraries
print $out_fh "# Load necessary R libraries\n";
print $out_fh "library(fields)\n";
print $out_fh "library(qrcode)\n\n";

# Change to output directory for qrcodes
print $out_fh "# Change to output directory for qrcodes\n";
print $out_fh "setwd(\"C:/Users/taraa/Documents/QR_code/\")\n\n";

while (my $line = <$in_fh>) {
	chomp $line;

	# Split line and store sample data in strings
	my ($plot, $experiment, $dateplanted, $seedsource, $seedyear, $genotype, $treatment) = split ("\t", $line);

	# store information for qrcode
	my $qrcode = "Plot:$plot;Experiment:$experiment;Planted:$dateplanted;SeedSource:$seedsource;SeedYear:$seedyear;Genotype:$genotype;Treatment:$treatment";

###### Printing R script
	# Store the qrcode as a matrix
	print $out_fh "# Store the qrcode as a matrix\n";
	print $out_fh "qr<-qrcode_gen('$qrcode', mask=1, dataOutput=TRUE, ErrorCorrectionLevel=\"L\")\n\n";

	# Setup pdf output -> file name: experiment#_qrcode_plot# (ie: 1_qrcode_1.pdf)
	print $out_fh "# Setup pdf output\n";
	print $out_fh "pdf(\"$experiment\_qrcode_$plot.pdf\", height=8.5, width=11, pointsize=18)\n\n";

	# Setup blank plot: remove as much border on page as possible
	print $out_fh "# Setup blank plot: remove as much border on page as possible\n";
	print $out_fh "par(mai=c(0, 0, .1, 0))\n";
	print $out_fh "plot(c(0,10), c(0,10), type=\"n\", axes=F, xlab=\"\", ylab=\"\")\n\n";

	# Add border around plot -> probably remove in final version
	print $out_fh "# Add border around plot -> probably remove in final version\n";
	print $out_fh "rect(xleft=0, ybottom=0, xright=10, ytop=10, border=\"red\", lwd=15)\n\n";

	# Add the qr image to plot in lower left hand
	print $out_fh "# Add the qr image to plot in lower left hand\n";
	print $out_fh "add.image(1.5, 8.25, qr, col=gray((32:0)/32), image.width=.25, image.height=.25)\n\n";

	# Add red lines to mark regions of page
	print $out_fh "# Add red lines to mark regions of page\n";
	print $out_fh "lines(c(3, 3), c(6.5,10), type=\"l\", lwd=15, col=\"red\", lend=1)\n";
	print $out_fh "lines(c(0,10), c(6.5,6.5), type=\"l\", lwd=15, col=\"red\", lend=1)\n\n";

	# Add boxes to mark day of picture, with text inside
	print $out_fh "# Add boxes to mark day of picture, with text inside\n";
	print $out_fh "rect(xleft=0.5, ybottom=5, xright=1.5, ytop=6, border=\"blue\", lwd=10)\n";
	print $out_fh "text(1,5.5, labels=\"Day4\")\n\n";

	print $out_fh "rect(xleft=0.5, ybottom=3.5, xright=1.5, ytop=4.5, border=\"blue\", lwd=10)\n";
	print $out_fh "text(1,4, labels=\"Day5\")\n\n";

	print $out_fh "rect(xleft=0.5, ybottom=2, xright=1.5, ytop=3, border=\"blue\", lwd=10)\n";
	print $out_fh "text(1,2.5, labels=\"Day6\")\n\n";

	print $out_fh "rect(xleft=0.5, ybottom=.5, xright=1.5, ytop=1.5, border=\"blue\", lwd=10)\n";
	print $out_fh "text(1,1, labels=\"Day7\")\n\n";

	print $out_fh "rect(xleft=2, ybottom=5, xright=3, ytop=6, border=\"blue\", lwd=10)\n";
	print $out_fh "text(2.5,5.5, labels=\"Day8\")\n\n";

	print $out_fh "rect(xleft=2, ybottom=3.5, xright=3, ytop=4.5, border=\"blue\", lwd=10)\n";
	print $out_fh "text(2.5,4, labels=\"Day9\")\n\n";

	print $out_fh "rect(xleft=2, ybottom=2, xright=3, ytop=3, border=\"blue\", lwd=10)\n";
	print $out_fh "text(2.5,2.5, labels=\"Day10\")\n\n";

	print $out_fh "rect(xleft=2, ybottom=.5, xright=3, ytop=1.5, border=\"blue\", lwd=10)\n";
	print $out_fh "text(2.5,1, labels=\"Day11\")\n\n";

	print $out_fh "rect(xleft=3.5, ybottom=5, xright=4.5, ytop=6, border=\"blue\", lwd=10)\n";
	print $out_fh "text(4,5.5, labels=\"Day12\")\n\n";

	print $out_fh "rect(xleft=3.5, ybottom=3.5, xright=4.5, ytop=4.5, border=\"blue\", lwd=10)\n";
	print $out_fh "text(4,4, labels=\"Day13\")\n\n";

	print $out_fh "rect(xleft=3.5, ybottom=2, xright=4.5, ytop=3, border=\"blue\", lwd=10)\n";
	print $out_fh "text(4,2.5, labels=\"Day14\")\n\n";

	print $out_fh "rect(xleft=3.5, ybottom=.5, xright=4.5, ytop=1.5, border=\"blue\", lwd=10)\n";
	print $out_fh "text(4,1, labels=\"Day15\")\n\n";

	print $out_fh "rect(xleft=5, ybottom=5, xright=6, ytop=6, border=\"blue\", lwd=10)\n";
	print $out_fh "text(5.5,5.5, labels=\"Day16\")\n\n";

	print $out_fh "rect(xleft=5, ybottom=3.5, xright=6, ytop=4.5, border=\"blue\", lwd=10)\n";
	print $out_fh "text(5.5,4, labels=\"Day17\")\n\n";

	print $out_fh "rect(xleft=5, ybottom=2, xright=6, ytop=3, border=\"blue\", lwd=10)\n";
	print $out_fh "text(5.5,2.5, labels=\"Day18\")\n\n";

	print $out_fh "rect(xleft=5, ybottom=.5, xright=6, ytop=1.5, border=\"blue\", lwd=10)\n";
	print $out_fh "text(5.5,1, labels=\"Day19\")\n\n";

	print $out_fh "rect(xleft=6.5, ybottom=5, xright=7.5, ytop=6, border=\"blue\", lwd=10)\n";
	print $out_fh "text(7,5.5, labels=\"Day20\")\n\n";

	print $out_fh "rect(xleft=6.5, ybottom=3.5, xright=7.5, ytop=4.5, border=\"blue\", lwd=10)\n";
	print $out_fh "text(7,4, labels=\"Day21\")\n\n";

	print $out_fh "rect(xleft=6.5, ybottom=2, xright=7.5, ytop=3, border=\"blue\", lwd=10)\n";
	print $out_fh "text(7,2.5, labels=\"Day22\")\n\n";

	print $out_fh "rect(xleft=6.5, ybottom=.5, xright=7.5, ytop=1.5, border=\"blue\", lwd=10)\n";
	print $out_fh "text(7,1, labels=\"Day23\")\n\n";

	print $out_fh "rect(xleft=8, ybottom=5, xright=9, ytop=6, border=\"blue\", lwd=10)\n";
	print $out_fh "text(8.5,5.5, labels=\"Day24\")\n\n";

	print $out_fh "rect(xleft=8, ybottom=3.5, xright=9, ytop=4.5, border=\"blue\", lwd=10)\n";
	print $out_fh "text(8.5,4, labels=\"Day25\")\n\n";

	print $out_fh "rect(xleft=8, ybottom=2, xright=9, ytop=3, border=\"blue\", lwd=10)\n";
	print $out_fh "text(8.5,2.5, labels=\"Day26\")\n\n";

	print $out_fh "rect(xleft=8, ybottom=.5, xright=9, ytop=1.5, border=\"blue\", lwd=10)\n";
	print $out_fh "text(8.5,1, labels=\"Day27\")\n\n";

	# Write text in top panel for manual backup checking of picture if need to
	print $out_fh "# Write text in top panel for manual backup checking of picture if need to\n";
	print $out_fh "text(6.5,8.25, labels=\"Plot: $plot; Experiment: $experiment;\\nPlanted: $dateplanted;\\nSeedSource: $seedsource;\\nSeedYear: $seedyear; Genotype: $genotype;\\nTreatment: $treatment\", col=\"black\", cex=1.5)\n\n";

	# Turn off pdf
	print $out_fh "# Turn off pdf\n";
	print $out_fh "dev.off()\n\n";
	print $out_fh "########\n\n";
}
close $in_fh;
close $out_fh;

exit;
