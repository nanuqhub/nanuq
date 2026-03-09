#!/bin/bash

list_clean="`\find ./cfgs -name BLD` "

list_clean+="`\find ./cfgs -name WORK` "

list_clean+="`\find ./cfgs -name EXP00`"

list_clean+="`\find ./cfgs -name MY_SRC`"

echo "  ===> delete ${list_clean} !"
sleep 2

rm -rf ${list_clean}

rm -f ./mk/full_key_list.txt

rm -f mk/arch.history mk/arch_nemo.fcm mk/cpp.fcm mk/cpp.history

