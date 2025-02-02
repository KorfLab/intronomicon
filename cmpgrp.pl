use strict;
use warnings;

die "usage: $0 meta-training/*/*edited*\n" unless @ARGV;

# read all of the html files
my %data;
foreach my $file (@ARGV) {
	chomp $file;
	my ($name, $gse) = $file =~ /meta-training\/(\w+)\/(\w+)/;
	my $html = `cat $file`;
	while ($html =~ /(<dl.+?<\/dl>)/gs) {
		my $grptxt = $1;
		my %group;
		while ($grptxt =~ /(<dd.+?<\/dd>)/gs) {
			my $gsmtxt = $1;
			my ($gsm) = $gsmtxt =~ /(GSM\d+)/;
			$group{$gsm} = 1;
		}
		my ($tagtxt) = $grptxt =~ /Tags: (.+?)</;
		my @tag;
		while ($tagtxt =~ /(\w+)/gs) {push @tag, $1}
		push @{$data{$gse}{$name}}, {group => \%group, tags => \@tag};
	}
}

# report as yaml
foreach my $gse (sort keys %data) {
	next unless scalar keys %{$data{$gse}} > 1;
	print "- $gse:\n";
	foreach my $name (sort keys %{$data{$gse}}) {
		print "  - $name:\n";
		foreach my $grp (@{$data{$gse}{$name}}) {
			my @gsm = sort keys %{$grp->{group}};
			my @tag = @{$grp->{tags}};
			print "    - { ", "group: [", join(", ", @gsm), "], ",
				"tags: [", join(", ", @tag), "] }\n";
		}
	}
}

