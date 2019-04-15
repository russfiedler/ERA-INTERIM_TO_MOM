#!/bin/bash

# vars are the ERA -Interim names
#varf are the new

# acc=n this is an instantaneous/mean value in the original
# acc=y this is an accumulated value in the original which need to be converted to a flux
# assumes values are stored according to year/variable

#hack this to your heart's delight!

vars=(u10 v10 d2m t2m ci strd msl tp ssrd)
varf=(10u 10v 2d  2t  ci strd msl tp ssrd)
acc=(n   n   n   n   n     y    n   y  y   )
fill=(y   y   y   y   n    n    n   y  n   )

for y in 2019
do
for i in {0..8}
#for i in 9 11 12
#for i in {0..9}
do
vs=${vars[${i}]}
vf=${varf[${i}]}
ac=${acc[${i}]}
fi=${fill[${i}]}
  echo $vs $vf $ac $fi
  fin=/home/datalib/reanalyses/ERA-interim/daily_raw/lores/ei_${vf}_${y}.nc
  fout=era_${vs}_${y}.nc
  /bin/cp era_${vs}_template.nc $fout
  ./make_era_interim  -v $vs -x $vs -a $ac -f $fi -i $fin -o $fout
done
done
