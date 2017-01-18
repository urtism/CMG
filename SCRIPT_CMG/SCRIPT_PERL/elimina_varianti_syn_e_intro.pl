use warnings;
use strict;

my $one = 'ANN_exome.vcf';
my @vcf = &leggi_correggi($one);
my @header = &estrai_intestazione($one);
my $intestazione = join ('',@header);

#estraggo i dati dal file vcf e li salvo.
sub leggi_correggi {
                open(FILE,"<$one"); #apro in sola lettura
                my @leggiStr = <FILE>; #salvo il contenuto del file in @leggiStr
		close(FILE);
		return @leggiStr;
}

#estraggo l'intestazione
sub estrai_intestazione {
		my @try;
                open(FILE,"<$one"); #apro in sola lettura
                my @leggiStr = <FILE>; #salvo il contenuto del file in @leggiStr
		my $a = 0;
		my $b = 203;
		while ( $a<=$b ) { #con questo ciclo tolgo le prime righe dal file txt fino a raggiungere quella con #CHROM POS ID etc...
                my $file = shift @leggiStr; #tolgo la prima riga
		push @try, $file;
               	$a = $a+1;
		}
		close(FILE);
		return @try;
}


#-----------------con questo metto in ogni elemento dell'array le righe che contengono ENSG---------------------#
my @array;
foreach (@vcf) {
my $pp = $_;
#print $pp;
if ($pp =~ /ENSG[0-9]/) {
push @array,$pp;
	}
}

#---------------------------------elimino gli elementi contenenti varianti introniche ed esoniche------------------
my @path;
foreach (@array) {
	my $stringa = $_;
	if (($stringa !~ /\|synonymous_variant\|/g) and ($stringa !~ /\|intron_variant\|/g)) {
	push @path, $stringa;
	}
}

my $aa = join('',@path); #unisco tutto il testo in un'unica stringa

#--------------------------------stampo il nuovo file vcf--------------------------

my @save;
push @save, $intestazione, $aa; #Inserisco le righe con intestazione ## del vcf nell'array
my $file = join ('',@save);

my $out = 'D2_no_intron_synonymous_variants.vcf';
&salva_file($out,$file);


#-----------------------------------le mie sub----------------------------

#salvo le stringhe in un nuovo file.
sub salva_file {
        my ($name_file,$aaa) =@_;
        open FILE,">$name_file";
        print FILE $aaa;
        close (FILE);
}
