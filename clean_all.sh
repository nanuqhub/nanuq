#!/bin/bash

list_clean="`\find ./cfgs -name BLD` "

list_clean+="`\find ./cfgs -name WORK` "

list_clean+="`\find ./cfgs -name EXP00`"

echo "  ===> delete ${list_clean} !"
sleep 2

rm -rf ${list_clean}

rm -f ./mk/full_key_list.txt

rm -f ./mk/arch_nemo.fcm
