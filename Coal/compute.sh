#!/bin/bash
#cd deal1
for((ii=0;ii<20;ii++))
do
var=$(printf "%02d" "$ii")
mkdir coaled_results/ampt_zpc_00010$var
cp -r ../JFF_AMPT/amptdata21_nojet_pPb5p02TeV_0mb_zpc/ampt_zpc_00010$var.dat ./
./main  ampt_zpc_00010$var.dat

mv remnant_jet_parton.dat coaled_results/ampt_zpc_00010$var/
mv thermal_jet.dat coaled_results/ampt_zpc_00010$var/
rm -r ampt_zpc_00010$var.dat
done
