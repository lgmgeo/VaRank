#!/bin/bash

dir=$1

MOUNTPOINT=$(df $dir | tail -n1 | awk '{print $1}')

FSTYPE=$(mount -l|grep $MOUNTPOINT | sed 's/.* type \([^ ]*\) .*/\1/')

echo $FSTYPE
