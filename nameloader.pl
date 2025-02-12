use strict;
use warnings;

die "usage: $0 <download dir> <gene names>\n" unless @ARGV == 2;
my ($DIR, $IN) = @ARGV;

my $dir = "$DIR/gene_names";
`mkdir -p $dir`;
open(my $fh, $IN) or die;
my $base = "https://rest.ensembl.org/xrefs/id";
my $header = <$fh>;
while (<$fh>) {
	my ($gene, $txlen) = split;
	next if -s "$dir/$gene.json";
	my $url = "\"$base/$gene?content-type=application/json\"";
	print $url, "\n";
	`wget -O $dir/$gene.json $url`;
	sleep(0.5)
}