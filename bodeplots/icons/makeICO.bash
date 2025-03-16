#!/bin/bash

# Just pass the original .png image as the only parameter to this script.
SOURCE="$1"
BASE=`basename "${SOURCE}" .png`

convert "${SOURCE}" -thumbnail 16x16 "${BASE}-16.png"
convert "${SOURCE}" -thumbnail 32x32 "${BASE}-32.png"
convert "${SOURCE}" -thumbnail 64x64 "${BASE}-64.png"
convert "${SOURCE}" -thumbnail 128x128 "${BASE}-128.png"
convert "${SOURCE}" -thumbnail 256x256 "${BASE}-256.png"
convert "${SOURCE}" -thumbnail 512x512 "${BASE}-512.png"
icotool -c -o "${BASE}.ico" "${BASE}"-{16,32,64,128,256,512}.png

