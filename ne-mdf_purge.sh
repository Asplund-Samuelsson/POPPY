#!/bin/bash

# Infiles and outfiles
mdf_script=`echo $0 | rev | cut -f 2- -d \/ | rev | tr "\n" "/"; echo "mdf.py"`
drg_script=`echo $0 | rev | cut -f 2- -d \/ | rev | tr "\n" "/"; echo "poppy_rank.py"`
network_in=$1
bounds=$2
ratios=$3
gibbs=$4
ph=$5
network_out=$6
bounds_out=$7

# Additional parameters
cmax=0.1
cmin=0.0000001
n_samples=10000
cpu=32

# Create a copy of the network file
cp $network_in /tmp/network
network="/tmp/network"

# Export variables
export mdf_script drg_script network bounds ratios gibbs ph cmax cmin

# Clear out /tmp
rm /tmp/mdf*

### Section 1: Purge network

# 1.1) Create a concentration bounds file that contains only the H2O bounds
grep C00001 $bounds > /tmp/bounds

round=0

while true; do

  let round=round+1
  echo ""
  echo "Reaction purge round $round: Performing $n_samples random samples..."

  # 1.2) Perform n random samples
  # 1.2.1) Randomly pick 50 reactions to represent the network.
  # 1.2.2) Create drG file
  # 1.2.3) Randomly pick a “pathway”: Set of 5 reactions from the network to do MDF optimization on.

  seq 0 $n_samples | parallel --no-notice --bar --jobs $cpu 'shuf -n 70 $network > /tmp/net{};
  shuf -n 5 /tmp/net{} | cut -f 1 > /tmp/pw{};
  $drg_script --write_gibbs --gibbs $gibbs --pH $ph /tmp/net{} /tmp/drg{};
  $mdf_script --constraints /tmp/bounds --ratios $ratios --min_conc $cmin --max_conc $cmax --pathway /tmp/pw{} /tmp/net{} /tmp/drg{} /tmp/mdf{} > /dev/null;'

  # 1.3) Count failures
  n_fail=`cat /tmp/mdf* | rev | cut -f 1 -d , | rev | grep -v MDF | grep NA | wc -l`

  # 1.3.X) Stop if there are no more failures
  if [[ $n_fail == 0 ]]; then
    break
  fi

  # 1.4) Otherwise remove the most failing reaction and continue
  rxn_fail_champion=`ls /tmp/mdf* | while read File; do cat $File | tr "\n" "\t" | sed -e "s/\t$/\n/" | grep -P "NA$" | tr "," "\n" | grep drGopt | cut -f 2 -d \_; done | sort | uniq -c | sort -n | tail -1 | sed -e 's/^ \+//g' | cut -f 2 -d \ `
  grep -vP "^${rxn_fail_champion}\t" $network > /tmp/network2
  mv /tmp/network2 $network

  echo "Removed $rxn_fail_champion."

done

# 1.5) Save the purged network
cp $network $network_out


### Section 2: Purge concentration bounds
