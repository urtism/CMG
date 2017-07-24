VCF ANNOTATOR

Author: Matteo Di Giovannantonio
Contact: matteodeg@gmail.com

DESCRIPTION:
VCF annotator lavora su un file vcf in input, annotato con VEP, per aggiungere dati estratti da database estratti. Con VCF ANNOTATOR
puoi inserire le features dei seguenti database nel formato di sotto richiesto nel vcf di output:

-CLINVAR (formato vcf)
-HUMSAVAR (tab delimited)
-ESP Exome Variant Server (formato vcf)
-GERP score (tab delimited)
-PHASTCONS (tab delimited)
-PHYLOP (tab delimited)

Con vcf annotator puoi aggiungere i campi che vuoi dai database. Ad esempio, in command line scrivere:

python vcf_Annotator -i input.vcf -o output.vcf -clnv path/to/clinvar.vcf,CLNSIG,CLNDB

In questo caso, annoto il vcf con gli attributi richiesti in ingresso a clinvar, cioè CLNSIG e CLNDB