#!/bin/bash
pkill gvfs
D=$(date +"%m-%d-%y")
folder=$HOME'/maizeData/seedlingData/'$D'/'
echo $folder
mkdir -p $folder
tmpLocation=$HOME'/tmpData/'
mkdir -p $tmpLocation
#gphoto2 --set-config /main/capturesettings/f-number=16
#gphoto2 --set-config /main/capturesettings/shutterspeed=31
gphoto2 --capture-image-and-download --force-overwrite --filename=$tmpLocation/tmpImage.nef

