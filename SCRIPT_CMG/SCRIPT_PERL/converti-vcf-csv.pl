use warnings;
use strict;

my $one = 'prova.vcf';
my @vcf = &leggi_correggi($one);
#print @vcf; #Ogni elemento del vettore contiene una riga del mio vcf.


#tolgo le prime righe dal file vcf
sub leggi_correggi {
                open(FILE,"<$one"); #apro in sola lettura
                my @leggiStr = <FILE>; #salvo il contenuto del file in @leggiStr
		my $a = 0;
		my $b = 203;
		while ( $a<=$b ) { #con questo ciclo tolgo le prime righe dal file txt fino a raggiungere quella con #CHROM POS ID etc...
                my $file = shift @leggiStr; #tolgo la prima riga
		#print $file; contiene le righe eliminate
                $a = $a+1;
		}
		close(FILE);
		return @leggiStr;
}

my $aa = join('',@vcf);
#print $aa; #metto tutto il testo in un'unica variabile
#my @bb = split ('\t|;|:',$aa);
#print @bb[50]; #in ogni elemento dell'array abbaimo un solo elemento
#my $num = @bb; #contiene il numero totale di elementi nell'array
$aa =~ s/AC=.+SOR=\d+\.+[0-9]+;//g; #levo ogni carattere a partire da AC= fino SOR=...; così da avere solo l'annotazione nella colonna info
$aa =~ s/INFO/ANNOTATION/; #Sostituisco ad INFO la sigla ANNOTATION
$aa =~ s/FORMAT/\tGT\tAD\DP/; #Sostituisco ad INFO la sigla ANNOTATION
$aa =~ s/:GQ:PL//g; #Elimino i caratteri :GQ:PL

#$aa =~ s/\:[0-9\,]+\:[0-9\,]+\:[0-9\,]+\t/\t/g;

$aa =~ s/\:[0-9\,]+\:[0-9\,]+\t/\t/g; #Elimino i numeri corrispondenti alle tag sopra elencate
print $aa;

# $aa =~ /\t[0-9]\:[0-9]\:/
# $aa =~ s/\t([0-9])\/([0-9])\:/\t$1\/$2\t/g

# $aa =~ s/\t(\d+)\/(\d+)\:/\t$1\/$2\t/g
# $aa =~ s/\t[(\d+)\,]+\:(\d+)/\t$1\t$2\t/g
# $aa =~ s/\t\d+\,//



#Salvo i risultati su un file .txt
my $newfile = 'Results.csv';
&salva_file($newfile,$aa);

sub salva_file {
        my ($name_file) =@_;
        open FILE,">$name_file";
        print FILE $aa;
        close (FILE);
}

