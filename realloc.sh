#!/bin/bash

for file in `find . -name *.c`
do
	echo "$file"

	# var = (type*(*...*)) realloc(var , (...) * sizeof(type(*...*)))
	cat $file | sed -e 's/\(\s*\)\([[:alnum:][:punct:][:space:]]\+\)\s*=.*\s*realloc\s*([[:alnum:][:punct:][:space:]]\+,\([[:alnum:][:punct:][:space:]]\+\)\*\s*sizeof\s*(\([[:alnum:][:punct:][:space:]]\+\))\s*)/\1PDM_realloc(\2,\2,\3,\4)/g' > $file.1
  mv $file.1 $file

	# var = (type*(*...*)) realloc(var, sizeof(type(*...*)) * (...))
	cat $file | sed -e 's/\(\s*\)\([[:alnum:][:punct:][:space:]]\+\)=.*realloc\s*([[:alnum:][:punct:][:space:]]\+,\s*sizeof\s*(\([[:alnum:][:punct:][:space:]]\+\))\s*\*\([[:alnum:][:punct:][:space:]]\+\))/\1PDM_realloc(\2,\2,\4,\3)/g'  > $file.1
  mv  $file.1 $file


done
