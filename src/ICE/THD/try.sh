#!/bin/bash


file=$1


filer=`echo ${file} | cut -d'.' -f1`



# 1: remove all blanks:
sed -e s/' '/''/g ${file} > file1.tmp

# 2: remove everything commented and other useless lines:
cat file1.tmp | grep -Ev '^!|^USE|^REAL|^INTEGER|^DO|^END|^CASE' | cut -d'!' -f1 > file2.tmp

# 3: remove all blanks:
sed -e s/'ELSE'/''/g file2.tmp > file3.tmp


#######################
#### Updated ARRAYS ###
#######################
# 4: here keep only lines in which something is given a value:
cat file3.tmp | grep '\=' | grep -Ev '==|<=|>=|/=' > file_equal.tmp

# 5: keep only what seems to be an array given a value, excluding local arrays (begining with 'z'):
cat file_equal.tmp | cut -d'=' -f1 | grep '(' | grep -v '^z' > file_array_equal.tmp

# 6: remove the content between parenthesis:
cat file_array_equal.tmp | cut -d '(' -f1 > file_array_equal2.tmp

# 7: remove duplicate words:
sort file_array_equal2.tmp | uniq > list_updated_arrays_${filer}.tmp


############################
#### All ARRAYS present  ###
############################

# 4b: replace operators by a new line:
sed  -e s/'*'/'\n'/g  -e s/'+'/'\n'/g  -e s/'-'/'\n'/g  -e s/'\/'/'\n'/g  -e s/'<'/'\n'/g  file3.tmp  > file3b.tmp

# 4b: keep only what seems to be an array, excluding local arrays (begining with 'z'):
cat file3b.tmp | grep -v '^1' | cut -d'=' -f1 | grep '(' | grep -v '^z'     | grep '(j' > file_array.tmp

# 5b: this needs some cleaning, and replace space by a newline:
sed -e s/'IF('/''/g -e s/')THEN'/''/g -e s/').AND.('/' '/g -e s/'.AND.'/' '/g -e s/'&'/''/g -e s/'\/'/''/g  -e s/'^('/''/g \
    -e s/'MAX(0._wp,'/''/g  -e s/'SIGN(1._wp,'/''/g   -e s/'MAX('/''/g  -e s/'3,'/''/g \
    file_array.tmp | tr ' ' '\n' > file_array2.tmp

# 6b: replace operators by a new line, and remove '(' at the begining:
sed -e s/'>'/'\n'/g file_array2.tmp -e s/'<'/'\n'/g  file_array2.tmp  > file_array3.tmp

# 7b: again, keep lines that contain array stuff:
cat file_array3.tmp | grep '(j' > file_array4.tmp

# 8b: now there should be one array per line and we can dump everything after '(', and ignore local arrays:
cat file_array4.tmp | cut -d'(' -f1 | grep -v '^z' >file_array5.tmp

# 9b: remove duplicate words:
sort file_array5.tmp | uniq > list_present_arrays_${filer}.tmp


###################################################
#### Now list of array present but not updated ####
###################################################
comm -23 list_present_arrays_${filer}.tmp list_updated_arrays_${filer}.tmp > list_NotUpdated_arrays_${filer}.tmp


# Finally making appropriate for pasting:
for ff in "list_present_arrays" "list_updated_arrays" "list_NotUpdated_arrays"; do
    tr '\n' ',' < ${ff}_${filer}.tmp > ${ff}_${filer}.txt
done
