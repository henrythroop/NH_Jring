# Shell script from Simon to retrieve latest OD kernels from ixion, and save locally.
# Also converts from .xsp to bsp
# 1-Dec-2018

cspice=/usr/local/naif/latest/icy

idir="$(ssh ixion.boulder.swri.edu ls -1d /home/nhsciops/kernels/nav_doppler/solutionsets/OD???|tail -1)/"
rsync --exclude "*.txt" -av ixion.boulder.swri.edu:$idir kernels/

for ifn in kernels/*.xsp; do
    ofn="$(echo "$ifn" | cut -f 1 -d '.').bsp"
    if [ -f $ofn ]; then
        continue
    fi
    $cspice/exe/tobin $ifn
done
