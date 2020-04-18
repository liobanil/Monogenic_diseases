####################################### README de comandos de bash

######################### para filtrar OMIM
tail -n +4 OMIM-Entry-Retrieval_singlegene.txt | cut -f-1,4,5 > OMIM-Entry-Retrieval_singlegene2.txt
grep -E "[:digit:]{1,2}\.[:digit:]{8}-[0-9]{8}" OMIM-Entry-Retrieval_monogen2.txt 

sed 's/-/\t/g' OMIM-Entry-Retrieval_singlegene2.txt | sed 's/:/\t/g' > OMIM-Entry-Retrieval_singlegene3.txt


########################## para filtrar KEGG
grep -E "[HSA:[0-9]{4,}]" DISEASE_\(11\)_kegg.txt | sed '/^hsa/d' | sed '/^H[0-9]/d' | sed -e 's/\[HSA:.*//g'
