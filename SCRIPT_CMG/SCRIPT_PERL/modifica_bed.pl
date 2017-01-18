use warnings;
use strict;

my $one = 'target_modificato.txt';
my @vcf = &leggi_correggi($one);

#estraggo i dati dal file vcf e li salvo.
sub leggi_correggi {
                open(FILE,"<$one"); #apro in sola lettura
                my @leggiStr = <FILE>; #salvo il contenuto del file in @leggiStr
		close(FILE);
		return @leggiStr;
}

foreach (@vcf) {
$_ =~ s/\.chr[0-9]\.[0-9]+\.[0-9]+//;
print $_;
}
