#!/bin/sh

filename=test.txt

content=$(< $filename)

while read -r line
do
 echo "$line"
done < $filename

echo $content

