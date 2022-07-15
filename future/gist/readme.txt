#move file "setup" into a new empty directory



# Run new fresh Ubuntu 16.04 docker image
docker run -it /home/username/gen/SRR/SingleCellRNASeq:/home ubuntu:16.04

cd /home
chmod +x setup

#For mouse genome
./setup m

#For human genome
./setup h

# Both mouse and human
./setup hm

# Wait...

# When setup is done run an example
cd /home
./calculate m SRR1784315 799 /home/SingleCellRNASeq/dropEst/configs/indrop_v1_2.xml clean
# Which mean to preprocess raw read counts from SRR1784315, 
# where we have 799 cells
# use barcoding scheme from indrop_v1_2.xml
# and perform cleaning of all downloaded files after obtaining read counts matrix

# Save container for a future use
docker commit
