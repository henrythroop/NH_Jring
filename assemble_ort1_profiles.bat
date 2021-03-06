# Create the directory locally, and clear it out

cd ~/git/nh_rings
rm -rf out/ort1_profiles

# Create the subdirs

mkdir out/ort1_profiles
mkdir out/ort1_profiles/HAZ01
mkdir out/ort1_profiles/HAZ02
mkdir out/ort1_profiles/HAZ03
mkdir out/ort1_profiles/HAZ04

# Copy local files from various places, into the out/ dir locally

cd /Users/throop/Data/ORT1/throop/backplaned/K1LR_HAZ01/
\cp stack*z4*profile.txt stack*z4*png stack*z4*profile.txt stack*z4*fits ~/git/nh_rings/out/ort1_profiles/HAZ01/

cd /Users/throop/Data/ORT1/throop/backplaned/K1LR_HAZ02/
\cp stack*z4*profile.txt stack*z4*png stackz4**profile.txt stack*z4*fits ~/git/nh_rings/out/ort1_profiles/HAZ02/

cd /Users/throop/Data/ORT1/throop/backplaned/K1LR_HAZ03/
\cp stack*profile.txt stack*z4*png stack*z4*profile.txt stack*z4*fits ~/git/nh_rings/out/ort1_profiles/HAZ03/

cd /Users/throop/Data/ORT1/throop/backplaned/K1LR_HAZ04/
\cp stack*profile.txt stack*z4*png stack*z4*profile.txt stack*z4*fits ~/git/nh_rings/out/ort1_profiles/HAZ04/

# Copy everything from local machine to ixion

scp -r ~/git/nh_rings/out/ort1_profiles ixion:ORT1/astrom/ort1_track1

scp -r ~/git/nh_rings/out/ort1_throop_track4.tar.gz ixion:ORT1/astrom/track4
