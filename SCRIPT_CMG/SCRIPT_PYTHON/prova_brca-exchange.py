from ga4gh.client import client
c = client.HttpClient("http://brcaexchange.org/backend/data/ga4gh/v0.6.0a7/")
dataset = c.search_datasets().next()
#print dataset



individual_dataset = c.get_dataset(dataset_id="brca")
print individual_dataset

# variant_sets = [i for i in c.search_variant_sets(dataset_id=dataset.id)]
# Sets = {}
# for variantSets in variant_sets:
#     Sets[variantSets.id] = {"Name" : variantSets.name, "Reference Set Id" : variantSets.reference_set_id,
#                             "Data Set Id" : variantSets.dataset_id}
#     print"Variant Set Id: {}\n\tName: {}\n\tReference Set Id: {}\n\tData Set Id: {}\n".format(variantSets.id,
#         variantSets.name, variantSets.reference_set_id, variantSets.dataset_id)


# Varset = c.get_variant_set(variant_set_id="brca-hg37")
# print "Variant Id: {}\nName: {}\nDataset Id: {}\nReference Set Id: {}\n".format(Varset.id, Varset.name, Varset.dataset_id, Varset.reference_set_id)
# for i in Varset.metadata:
#     print "Metadata Field: {} ;  Value: {} ;  Type: {}".format(i.key, i.value, i.type)




# counter = 0
# for Vars in c.search_variants(reference_name="chr17", variant_set_id="brca-hg37", start=41246794, end=41296814):
#     print "Variant Id: {},\tVariant Set Id: {},\tchr: {}\n\tVariant Start: {},\tVariant End: {}\n\tReference Bases: {},\tAlternate Bases: {},\t".format(
#         Vars.id, Vars.variant_set_id,Vars.reference_name,Vars.start, Vars.end,Vars.reference_bases,Vars.alternate_bases)
#     if counter >= 5:
#         break
#     counter += 1





#SingleVar = c.get_variant(variant_id="hg37-32793")
#print "Variant Id: {},\tVariant Set Id: {},\tReference Name: {}\n\tVariant Start: {},\tVariant End: {}\n\tReference Bases: {},\tAlternate Bases: {},\t".format(SingleVar.id, #SingleVar.variant_set_id,SingleVar.reference_name,SingleVar.start, SingleVar.end,SingleVar.reference_bases,SingleVar.alternate_bases)
#for i in SingleVar.info:
    #print "{}: \t{}".format(i, SingleVar.info[str(i)].values[0].string_value or SingleVar.info[str(i)].values[0].number_value)
    #print i

#print SingleVar.info['EAS_Allele_frequency_1000_Genomes'].values[0].string_value
#print SingleVar




