#!/bin/sh
if [ $# -ne 1 ]; then
	echo "$0 file"
	exit
fi
file=$1
#destdir=$2
#echo `pwd`
#echo $file
rsync -avzP -e "ssh -p 20184" $file wbzhao@localhost:/lustre/work2/xnwang/wbzhao/CoLBT_running/sampler_UrQMD/Zyvisc_sampler_UrQMD/coal_frag_with_thth/

