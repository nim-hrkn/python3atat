export PATH=/home/kino/kino/kit/atat/bin:$PATH
#extract_vasp
#cp str_relax.out str.out
python /home/kino/kino/kit/python-atat/atat/poscar2str.py opt.vasp
cp str.out str_relax.out
fitfc -er=12 -ns=1 -dr=0.1 -nrr

for name in vol_0/p*; 
do (echo $name; cd $name; python /home/kino/kino/kit/python-atat/atat/str2poscar.py str.out );
done

