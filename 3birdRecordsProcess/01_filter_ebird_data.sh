# how to execute this script in a linux or mac terminal:
# option 1: run it directly using bash
# bash 01_filter_ebird_data.sh
# 
# option 2: give it execution permissions and run it
# chmod +x 01_filter_ebird_data.sh
# ./01_filter_ebird_data.sh

# why this bash method is used instead of r:
# the global ebird file is too massive to extract or load into ram.
# awk processes the data line by line keeping memory usage near zero,
# filtering and recompressing on the fly without creating huge temp files.

# filter 1: european countries
LC_ALL=C awk -F'\t' -v codes='RU,AL,AT,BA,BE,BG,BY,CH,DK,ES,EE,GR,CZ,DE,FR,FI,HR,HU,IE,IM,IT,LT,LU,LV,MK,MD,ME,NL,NO,RO,RS,PL,PT,SE,SI,UA,SK,GB,AD,LI,SM' '
    BEGIN { n = split(codes, a, ","); for(i=1; i<=n; i++) keep[a[i]] = 1 }
    NR == 1 || (toupper($17) in keep)
' ebd_relJun-2025.txt \
| awk -F'\t' '
    # filter 2: year >= 2000
    NR == 1 || substr($31, 1, 4) >= 2000
' \
| awk -F'\t' '
    # filter 3: non-incidental protocol
    NR == 1 || $37 != "Incidental"
' \
| awk -F'\t' '
    # filter 4: clean problematic species names (scientific name is column 7)
    NR == 1 || ($7 !~ / sp\./ && $7 !~ /\// && $7 !~ / x / && $7 !~ /\(Domestic type\)/)
' \
| awk -F'\t' -v OFS='\t' '
    # final selection of 12 columns
    { print $7, $11, $17, $29, $30, $31, $37, $38, $41, $42, $45, $14 }
' \
| gzip -c > ebd_europe_prefiltered_clean.tsv.gz