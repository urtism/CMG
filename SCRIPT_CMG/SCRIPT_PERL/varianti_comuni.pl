use warnings;
use strict;

#--------con questo script estraggo i geni che presentano le stesse varianti per due campioni ma il terzo ha almeno una variante su quel gene-----------

my $one = 'common_D1_01_02.vcf';
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
		my $b = 209;
		while ( $a<=$b ) { #con questo ciclo tolgo le prime righe dal file txt fino a raggiungere quella con #CHROM POS ID etc...
                my $file = shift @leggiStr; #tolgo la prima riga
		push @try, $file;
               	$a = $a+1;
		}
		close(FILE);
		return @try;
}

#-----------------con questo metto in ogni elemento dell'array le righe che contengono ENSG---------------------#
#my @array;
#foreach (@vcf) {
#my $pp = $_;
#print $pp;
#if ($pp =~ /ENSG[0-9]/) {
#push @array,$pp;
#	}
#}

#--------------metto in un array solo le parole ENSG (geni)------------
my $aa = join('',@vcf); #unisco tutto il testo in un'unica stringa
my @array = $aa =~ /ENSG[0-9]+/g; #cerco solo ENSG
my %hash = map { $_ , 1 } @array; #inserisco in un hash gli elementi dell'array che sono unici

#----------inserisco nelle values le sole righe appartenenti al gene presente nella keys---------
foreach my $k(keys %hash) {
my @new;
	foreach (@vcf) {
		my $bb = $_;
			if ($bb =~ /$k/) {
			push @new,$bb;
			}
		my $str = join ('',@new);	
		$hash{$k} = $str;	
	}
}
#print %hash;
#--------------------------controllo che su ogni gene tutti i pazienti abbiano almeno una variante-----------------------------------
#ricordarsi di implementare tutto in un unico ciclo for!
#in questo caso il vcf contiene le varianti che sono comuni ai primi due cammpioni. Quindi controllo che ci sia la variante sul terzo campione
#sul gene. se il campione 3 presenta la variante sul gene, allora mantengo il gene altrimenti lo elimino

#PAZIENTE SU TERZA COLONNA
foreach my $k(keys %hash) {
my @split;
@split = split /\t|\n/ , $hash{$k};
my $num = scalar(@split)-1;
my @GT;
my $m=11;
	while ($m<=$num) {
	if ( $split[$m] =~ /0\/0|\.\/\.|[0-9]\/[0-9]/) {
		my $stringa = $hash{$k};
		push @GT,$split[$m];
		}
	$m = $m+15;
	}
my $stringa = join ('',@GT);
if ($stringa !~ /[0-9]+\/[1-9]\:/) {
		delete $hash{$k};
		}
}


#------------------------------Estraggo una lista di geni comuni dal vcf------------------------------
my @genes = keys (%hash);

#------------------------------estrai nome del gene-----------------------------

my @nome_gene = $aa =~ /\|[A-Z0-9\.\-]+\|ENSG[0-9]+/g; #cerco solo ENSG
my %dash = map { $_ , 1 } @nome_gene; #inserisco in un hash gli elementi dell'array che sono unici
my @prova = keys (%dash);

my @push;
foreach (@prova) {
	my $string = $_;
	$string =~ /ENSG[0-9]+/;
	my $boh = $&; #prendo la sola parte che matcha allo step precedente.
	#print "\n$boh\n";
	foreach (@genes) { #faccio matchare la stringa trovata nel precedente ciclo for con le stringhe presenti nelle keys dellìhash precedente.
		my $tr= $_;
		if ($tr eq $boh) { #se le stringhe sono uguali, salvo il codice ENSG e il nome del gene
		push @push, $string;
		}
	}
}

my @free;
foreach (@push) {
	my $true = $_;
	$true =~ s/(ENSG[0-9]{11})/\t$1\n/;
	$true =~ s/\|//g;
	push @free, $true;
}

my @ord = sort @free;
my $ENSG = join ('', @ord);

#--------------------------------stampo il nuovo file vcf--------------------------

my @save;
push @save, $intestazione; #Inserisco le righe con intestazione ## del vcf nell'array

foreach my $k(keys %hash) {
my $test = $hash{$k};
push @save, $test; #inserisco l'hash nell'array
}

my $def= join ('',@save);
my $newfile = 'D1_geni_con_varianti_comuni.vcf';
&salva_file($newfile,$def);

#----------------------------------stampo la lista dei geni comuni------------------------------

my $list = 'lista_geni_common_D1_01_02.txt';
&salva_file($list,$ENSG);


#-----------------------------------le mie sub----------------------------

#salvo le stringhe in un nuovo file.
sub salva_file {
        my ($name_file,$aaa) =@_;
        open FILE,">$name_file";
        print FILE $aaa;
        close (FILE);
}
