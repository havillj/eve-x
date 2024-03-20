## Sample BLAST database

This BLAST viral database (named virusdb) contains complete sequences for all known RNA virus families that contain arboviruses.  Specifically, it contains sequences from the families/orders in the following table:

| Family/Order | NCBI taxid |
|-|-|
| Bromoviridae | 39740 |
| Chrysoviridae | 249310 |
| Flaviviridae | 11050 |
| Nodaviridae | 12283 |
| Orthomyxoviridae | 11308 |
| Partitiviridae | 11012 |
| Picornaviridae | 12058 |
| Sedoreovirinae | 689832 |
| Togaviridae | 11018 |
| Totiviridae | 11006 |
| Tymoviridae | 249184 |
| Mesoniviridae | 1312872 |
| Bunyavirales | 1980410 |
| Mononegavirales | 11157 |

It was created using the following command:

```
makeblastdb -in virusdb_no_retro_or_unverified.fasta -dbtype nucl -parse_seqids -title "virusdb" -out virusdb
```
