from repDNA.nac import Kmer
kmer = Kmer(k=2,normalize=True)
vec=kmer.make_kmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
print '\t'.join(['AA', 'AC', 'AG', 'AT', 
'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'])
#print kmer[0]
print '\t'.join(str(i) for i in vec[0])

from repDNA.nac import IDkmer

idkmer = IDkmer()
vec=idkmer.make_idkmer_vec('AAAA', 'AAAT','AAAA')

idkmer = IDkmer(k=2)
idkmer.make_idkmer_vec('AAAA', 'AAAT','AAAA')

idkmer = IDkmer(k=2, upto=False)

idkmer.make_idkmer_vec('AAAA', 'AAAT','AAAA')



















# from repDNA.psenac import SCPseTNC 
# >>> sc_psetnc = SCPseTNC() 
# >>> vec = 
# sc_psetnc.make_scpsetnc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], 
# phyche_index=['Dnase I', 'Nucleosome']) 
# >>> len(vec[0]) 
# 66 
# >>> sc_psetnc = SCPseTNC(lamada=2, w=0.05) 
# >>> vec = 
# sc_psetnc.make_scpsetnc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], 
# all_property=True) 
# >>> len(vec[0]) 


# import repDNA
# from repDNA.psenac import PseKNC

# #pseknc = PseKNC(k=2, lamada=1, w=0.05) 
# pseknc = PseKNC(k=4)
# phyche_index=['Twist']
# # vec = pseknc.make_pseknc_vec(['GACTG'])
# # vec1 =  pseknc.make_pseknc_vec(['GAATG'])
# # vec2 = pseknc.make_pseknc_vec(['GAGTG'])
# # vec3 = pseknc.make_pseknc_vec(['GATTG'])

# #vec = pseknc.make_pseknc_vec(['GAG'])
# vec = pseknc.make_pseknc_vec(['AAAAAC'])
# vec1 = pseknc.make_pseknc_vec(['AACAAC'])
# vec2 =  pseknc.make_pseknc_vec(['AAGAAC'])
# vec3 = pseknc.make_pseknc_vec(['AATAAC'])
# print '\t'.join(['AA', 'AC', 'AG', 'AT', 
# 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'])
# print '\t'.join(str(i) for i in vec[0])
# print '\t'.join(str(i) for i in vec1[0])
# print '\t'.join(str(i) for i in vec2[0])
# print '\t'.join(str(i) for i in vec3[0])
# #print vec1
# #print vec2
# #print vec3
# #print vec4

# print len(vec[0])



# >>> from repDNA.util import normalize_index 
# >>> 
# vec = 
# pseknc.make_pseknc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], 
# extra_phyche_index=normalize_index(phyche_index,is_convert_dict=True)) 