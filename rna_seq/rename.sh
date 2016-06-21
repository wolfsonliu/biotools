#! /bin/bash

# for cell in Ieko Rajo; do
#     filename=$(ls | grep ${cell})
#     if [ ${cell} == "Ieko" ]; then
#         changed="jeko"
#     else 
#         changed="raji"
#     fi
#     for name in ${filename}; do
#         mv ${name} ./${name/${cell}/${changed}}
#     done
# done

for filename in $(ls); do
    mv ${filename} ./$(echo $filename | sed -e 's/_HJ2WHCCXX_L[[:digit:]]//')
done

