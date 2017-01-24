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

# Set the failure limit
limit=`echo $n_samples/20 | bc` # <5%

# Create a random hash so that multiple scripts may run simultaneously
h=`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 16 | sed -e "s/^/_/"`

# Create a copy of the network file
network="/tmp/network${h}"
cp $network_in $network

# Export variables
export mdf_script drg_script network bounds ratios gibbs ph cmax cmin h

echo ""
echo "NE-MDF purge ${h}"
echo ""

### Section 1: Purge network

# Create a concentration bounds file that contains only the H2O bounds
grep C00001 $bounds > /tmp/bounds${h}

round=0

while true; do

  let round=round+1
  echo ""

  # Perform control NE-MDF for randomly sampled pathways and full network

  # Create drG file
  $drg_script --write_gibbs --gibbs $gibbs --pH $ph $network /tmp/drg${h}_control

  # Perform sampling
  echo "Reaction purge round $round: Control sampling with full network..."

  seq -f "%05g" 1 $n_samples | parallel --no-notice --bar --jobs $cpu 'shuf -n 5 $network | cut -f 1 > /tmp/pw${h}_{}_control;
  $mdf_script --constraints /tmp/bounds${h} --ratios $ratios --min_conc $cmin --max_conc $cmax --pathway /tmp/pw${h}_{}_control $network /tmp/drg${h}_control /tmp/mdf${h}_{}_control > /dev/null;'

  # Count failures
  n_fail=`cat /tmp/mdf${h}*_control | grep -c "NA$"`

  # Stop if there are less than the set number of failures
  if [[ $n_fail -lt $limit ]]; then
    echo "$n_fail / $n_samples failed. The network passed the limit of $limit failures."
    break
  fi

  # Otherwise continue with purge
  echo "$n_fail / $n_samples failed. The network did not pass the limit of $limit failures."
  echo "Reaction purge round $round: Sampling for reactions causing failures..."

  # Perform n random samples
  # Randomly pick 70 reactions to represent the network.
  # Create drG file
  # Randomly pick a “pathway”: Set of 5 reactions from the network to do MDF optimization on.

  seq -f "%05g" 1 $n_samples | parallel --no-notice --bar --jobs $cpu 'shuf -n 70 $network > /tmp/net${h}_{};
  shuf -n 5 /tmp/net${h}_{} | cut -f 1 > /tmp/pw${h}_{};
  $drg_script --write_gibbs --gibbs $gibbs --pH $ph /tmp/net${h}_{} /tmp/drg${h}_{};
  $mdf_script --constraints /tmp/bounds${h} --ratios $ratios --min_conc $cmin --max_conc $cmax --pathway /tmp/pw${h}_{} /tmp/net${h}_{} /tmp/drg${h}_{} /tmp/mdf${h}_{}_sample > /dev/null;'

  # Count failures
  n_fail=`cat /tmp/mdf${h}*_sample | grep -c "NA$"`

  # Stop if there are no more failures
  if [[ $n_fail -eq 0 ]]; then
    break
  fi

  # Otherwise remove the most failing reaction
  rxn_fail_champion=`ls /tmp/mdf${h}*_sample | while read File; do cat $File | tr "\n" "\t" | sed -e "s/\t$/\n/" | grep -P "NA$" | tr "," "\n" | grep drGopt | cut -f 2 -d \_; done | sort | uniq -c | sort -n | tail -1 | sed -e 's/^ \+//g' | cut -f 2 -d \ `
  grep -vP "^${rxn_fail_champion}\t" $network > /tmp/network2${h}
  mv /tmp/network2${h} $network

  echo "$n_fail / $n_samples failed. Removed $rxn_fail_champion."

done

# Save the purged network
echo ""
echo "Reaction purge finished. Saving output."
cp $network $network_out


### Section 2: Purge concentration bounds

# Create a bounds file without C00001
grep -v C00001 $bounds > /tmp/bounds2${h}

round=0

while true; do

  let round=round+1
  echo ""

  # Perform control NE-MDF for randomly sampled pathways and full network

  # Create drG file
  $drg_script --write_gibbs --gibbs $gibbs --pH $ph $network /tmp/drg${h}_control

  # Create full bounds file
  cat /tmp/bounds${h} /tmp/bounds2${h} > /tmp/bounds${h}_control

  # Perform sampling
  echo "Concentration purge round $round: Control sampling with full network..."

  seq -f "%05g" 1 $n_samples | parallel --no-notice --bar --jobs $cpu 'shuf -n 5 $network | cut -f 1 > /tmp/pw${h}_{}_control;
  $mdf_script --constraints /tmp/bounds${h}_control --ratios $ratios --min_conc $cmin --max_conc $cmax --pathway /tmp/pw${h}_{}_control $network /tmp/drg${h}_control /tmp/mdf${h}_{}_control > /dev/null;'

  # Count failures
  n_fail=`cat /tmp/mdf${h}*_control | grep -c "NA$"`

  # Stop if there are less than the set number of failures
  if [[ $n_fail -lt $limit ]]; then
    echo "$n_fail / $n_samples failed. The concentrations passed the limit of $limit failures."
    break
  fi

  # Otherwise continue with purge
  echo "$n_fail / $n_samples failed. The concentrations did not pass the limit of $limit failures."
  echo "Concentration purge round $round: Performing $n_samples random samples..."

  # Create drG file
  $drg_script --write_gibbs --gibbs $gibbs --pH $ph $network /tmp/drg${h}

  # Perform n random samples
  # Randomly pick 50 concentrations
  # Randomly pick a “pathway”: Set of 5 reactions from the network to do MDF optimization on.

  seq -f "%05g" 1 $n_samples | parallel --no-notice --bar --jobs $cpu '(cat /tmp/bounds${h}; shuf -n 50 /tmp/bounds2${h}) > /tmp/bnd${h}_{};
  shuf -n 5 $network | cut -f 1 > /tmp/pw${h}_{};
  $mdf_script --constraints /tmp/bnd${h}_{} --ratios $ratios --min_conc $cmin --max_conc $cmax --pathway /tmp/pw${h}_{} $network /tmp/drg${h} /tmp/mdf${h}_{}_sample > /dev/null;'

  # Count failures
  n_fail=`cat /tmp/mdf${h}*_sample | grep -c "NA$"`

  # Stop if there are no more failures
  if [[ $n_fail -eq 0 ]]; then
    break
  fi

  # Otherwise remove the most failing concentration and continue
  cpd_fail_champion=`grep "NA$" /tmp/mdf${h}*_sample | cut -f 1 -d \: | sed -e 's/mdf/bnd/' | while read File; do cat $File; done | cut -f 1 | grep -v C00001 | sort | uniq -c | sort -n | tail -1 | sed -e 's/^ \+//g' | cut -f 2 -d \ `
  grep -vP "^${cpd_fail_champion}\t" /tmp/bounds2${h} > /tmp/bounds3${h}
  mv /tmp/bounds3${h} /tmp/bounds2${h}

  echo "$n_fail / $n_samples failed. Removed $cpd_fail_champion bounds."

done

# Save the purged bounds
echo ""
echo "Concentration purge finished. Saving output."
cat /tmp/bounds${h} /tmp/bounds2${h} > $bounds_out
